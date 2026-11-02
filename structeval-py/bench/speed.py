"""Speed: structeval's vectorised BCa vs the original's per-draw groupby bootstrap."""
import time, numpy as np, pandas as pd
import structeval as se
from structeval.uncertainty import bootstrap_ci
from structeval.aggregate import set_metrics
import llmvalidate.validation as ov

def recall_fn(): return lambda t: set_metrics(t[0],t[1],t[2],t[3],t[4])["recall"]

N, B = 2000, 2000
rng = np.random.default_rng(0)
cor = rng.binomial(1, 0.7, N).astype(float)
counts = np.column_stack([cor, np.zeros(N), 1-cor, np.zeros(N), np.zeros(N)])

t = time.time()
for _ in range(3):
    bootstrap_ci(counts, recall_fn(), 0.95, B, np.random.default_rng(1))
se_t = (time.time()-t)/3

res = pd.DataFrame({'f':['a']*N,'Cor: f':cor,'Inc: f':0,'Mis: f':1-cor,'Spu: f':0,'Par: f':0})
t = time.time()
ov.bootstrap_CI(res, ['f'], n_bootstrap=B, random_state=1)
ov_t = time.time()-t

print(f"N={N} rows, B={B} draws, single field")
print(f"  original bootstrap_CI : {ov_t:7.3f}s")
print(f"  structeval BCa        : {se_t:7.3f}s")
print(f"  speedup               : {ov_t/se_t:7.1f}x")
