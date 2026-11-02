"""Empirical 95% CI coverage: structeval (BCa) vs the original (percentile).
Known true recall p; many datasets; fraction of intervals containing p."""
import numpy as np, pandas as pd, time
import structeval as se
from structeval.uncertainty import bootstrap_ci, wilson, z_for
from structeval.aggregate import set_metrics
import llmvalidate.validation as ov

def recall_fn(): return lambda t: set_metrics(t[0],t[1],t[2],t[3],t[4])["recall"]

def run(p_true=0.7, n=50, trials=250, n_boot=400, seed=0):
    rng = np.random.default_rng(seed); z = z_for(0.95); fn = recall_fn()
    hit = {"structeval BCa":0, "original percentile":0, "Wilson (ref)":0}
    width = {k:0.0 for k in hit}
    for _ in range(trials):
        cor = rng.binomial(1, p_true, n).astype(float)
        counts = np.column_stack([cor, np.zeros(n), 1-cor, np.zeros(n), np.zeros(n)])
        # structeval BCa
        r = bootstrap_ci(counts, fn, 0.95, n_boot, rng)
        if not np.isnan(r["ci_low"]):
            hit["structeval BCa"] += r["ci_low"] <= p_true <= r["ci_high"]
            width["structeval BCa"] += r["ci_high"] - r["ci_low"]
        # Wilson reference (clean binomial: k=cor, n=n)
        lo, hi = wilson(int(cor.sum()), n, z)
        hit["Wilson (ref)"] += lo <= p_true <= hi; width["Wilson (ref)"] += hi-lo
        # original percentile bootstrap
        res = pd.DataFrame({'f':['a']*n,'Cor: f':cor,'Inc: f':0,'Mis: f':1-cor,'Spu: f':0,'Par: f':0})
        ci = ov.bootstrap_CI(res, ['f'], n_bootstrap=n_boot, random_state=int(rng.integers(1e9)))
        row = ci[ci.field=='f'].iloc[0]
        lo2, hi2 = row['recall (micro): lower'], row['recall (micro): upper']
        if not np.isnan(lo2):
            hit["original percentile"] += lo2 <= p_true <= hi2
            width["original percentile"] += hi2-lo2
    print(f"\n--- p_true={p_true}, n={n}, trials={trials}, n_boot={n_boot} ---")
    print(f"{'method':22} {'coverage':>9} {'mean width':>11}")
    for k in hit:
        print(f"{k:22} {hit[k]/trials:9.3f} {width[k]/trials:11.3f}")

if __name__ == "__main__":
    t=time.time()
    run(p_true=0.7, n=50)
    run(p_true=0.9, n=40)
    print(f"\ntotal {time.time()-t:.1f}s")
