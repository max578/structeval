# =============================================================================
# Figure 2: size and power of the threshold and asymmetry tests
#
# Author      : Max Moldovan
# Version date: 29 October 2026
# Usage       : Rscript sim/figures/fig2_sizepower.R   (from the repository root)
#
# Data        : sim/results/{grid_a,grid_b,planner_validation}.csv
# Requirements: structeval, ggplot2, patchwork, ragg
# Outputs     : sim/figures/fig2_sizepower.{pdf,png,tiff}, 180 mm wide
#
# Structure   : 1. Read the simulation results
#               2. Build the four panels
#               3. Write the figure
# =============================================================================

require(structeval)
require(ggplot2)
require(patchwork)
require(ragg)

rm(list = ls())

ptm1 <- proc.time()   # measure execution time

dir_res <- "sim/results/"
dir_fig <- "sim/figures/"

# categorical slots 1-4 in fixed order (validated: adjacent CVD dE >= 9.1)
col_vec <- c("#2a78d6", "#eb6834", "#1baf7a", "#eda100")
ink_muted <- "#52514e"
base_size <- 8

theme_fig <- theme_minimal(base_size = base_size) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(colour = "#e4e3df", linewidth = 0.3),
    axis.line = element_line(colour = ink_muted, linewidth = 0.3),
    axis.ticks = element_line(colour = ink_muted, linewidth = 0.3),
    axis.text = element_text(colour = ink_muted),
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.key.width = unit(5, "mm"),
    legend.margin = margin(0, 0, 0, 0),
    strip.text = element_text(size = base_size, hjust = 0),
    plot.tag = element_text(face = "bold", size = base_size + 1),
    plot.background = element_rect(fill = "white", colour = NA)
  )

# =============================================================================
# 1. Read the simulation results
# =============================================================================
print("==== 1. Read the simulation results ====")
ptm_local <- proc.time()

tab_a <- read.csv(paste0(dir_res, "grid_a.csv"))
tab_b <- read.csv(paste0(dir_res, "grid_b.csv"))
tab_p <- read.csv(paste0(dir_res, "planner_validation.csv"))


time_section <- proc.time() - ptm_local
print("Section: 1. Read the simulation results")
print(time_section)

# =============================================================================
# 2. Build the four panels
# =============================================================================
print("==== 2. Build the four panels ====")
ptm_local <- proc.time()

# (a) exact power of the Wilson decision, simulated points at true 0.95
########################################################
n_vec <- 20:400
theta1_vec <- c(0.93, 0.95, 0.97)
tab_exact <- do.call(rbind, lapply(theta1_vec, function(th) {
  pw <- vapply(n_vec, function(nn) {
    k <- 0:nn
    lower <- structeval:::.wilson_ci(k, nn, 0.95)[, "lower"]
    sum(stats::dbinom(k[lower > 0.90], nn, th))
  }, numeric(1L))
  data.frame(n = n_vec, power = pw, truth = sprintf("true %.2f", th))
}))
tab_sim_a <- tab_a[tab_a$type == "proportion" & tab_a$theta == 0.95 &
  tab_a$method == "wilson", ]
tab_sim_a$truth <- "true 0.95"
# direct labels placed in clear space beside each curve
tab_lab_a <- data.frame(n = c(395, 300, 95), power = c(0.32, 0.62, 0.93),
  hj = c(1, 0, 1), truth = sprintf("true %.2f", theta1_vec))

fig_a <- ggplot(tab_exact, aes(n, power, colour = truth)) +
  geom_hline(yintercept = 0.8, linetype = "dashed", colour = ink_muted,
    linewidth = 0.3) +
  geom_line(linewidth = 0.5) +
  geom_errorbar(data = tab_sim_a, aes(ymin = power - 1.96 * power_se,
    ymax = power + 1.96 * power_se), width = 8, linewidth = 0.4,
    show.legend = FALSE) +
  geom_point(data = tab_sim_a, shape = 21, fill = "white", size = 1.6,
    stroke = 0.6, show.legend = FALSE) +
  geom_text(data = tab_lab_a, aes(label = truth, hjust = hj),
    colour = "#0b0b0b", size = 2.6, show.legend = FALSE) +
  scale_colour_manual(values = col_vec[1:3]) +
  scale_x_continuous(breaks = c(20, 100, 200, 300, 400),
    expand = expansion(mult = c(0.02, 0.04))) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  labs(x = "Notes (n)", y = "Power to certify 0.90") +
  theme_fig

# (b) notes required for a list F1 of 0.90 against 0.95, with validation
########################################################
tab_plan <- expand.grid(m = 1:4, rho = c(0, 0.2, 0.5))
tab_plan$n_notes <- n_threshold(0.90, 0.95, m = tab_plan$m,
  rho = tab_plan$rho)$n_notes
tab_plan$rho_lab <- factor(sprintf("rho = %.1f", tab_plan$rho),
  levels = sprintf("rho = %.1f", c(0, 0.2, 0.5)))
tab_val <- tab_p[tab_p$theta0 == 0.90 & tab_p$theta1 == 0.95 &
  tab_p$n_run == tab_p$n_planned, ]
tab_val$rho_lab <- factor(sprintf("rho = %.1f", tab_val$rho),
  levels = levels(tab_plan$rho_lab))
tab_val$lab <- sprintf("%.2f", tab_val$power)
# the m = 1 point is shared by every rho; draw it once, in the rho = 0 colour
tab_val$rho_lab[tab_val$m == 1] <- "rho = 0.0"

fig_b <- ggplot(tab_plan, aes(m, n_notes, colour = rho_lab)) +
  geom_line(linewidth = 0.5) +
  geom_point(shape = 16, size = 1.2) +
  geom_point(data = tab_val, aes(m, n_planned), shape = 23, fill = "white",
    size = 2.4, stroke = 0.7) +
  geom_text(data = tab_val, aes(m, n_planned, label = lab), colour = "#0b0b0b",
    size = 2.4, nudge_x = 0.2, hjust = 0, show.legend = FALSE) +
  scale_colour_manual(values = col_vec[1:3],
    labels = c(expression(rho == 0), expression(rho == 0.2),
      expression(rho == 0.5))) +
  scale_x_continuous(breaks = 1:4, expand = expansion(add = c(0.15, 0.45))) +
  scale_y_continuous(limits = c(0, 260), breaks = seq(0, 250, 50)) +
  labs(x = "Items per note (m)", y = "Notes required") +
  theme_fig

# (c) size of T1 in the deployment-like cell
########################################################
method_lab <- c(naive_wilson = "Naive Wilson", percentile = "Percentile",
  cluster_bca = "Cluster BCa", cluster_delta = "Cluster delta")
tab_c <- tab_a[tab_a$type == "list" & tab_a$n == 100 & tab_a$theta == 0.90 &
  tab_a$m == 4 & tab_a$prevalence == 1, ]
tab_c$method_f <- factor(method_lab[tab_c$method], levels = method_lab)
dodge <- position_dodge(width = 0.07)

fig_c <- ggplot(tab_c, aes(rho, size, colour = method_f, shape = method_f)) +
  geom_hline(yintercept = 0.025, linetype = "dashed", colour = ink_muted,
    linewidth = 0.3) +
  geom_line(position = dodge, linewidth = 0.5) +
  geom_errorbar(aes(ymin = pmax(0, size - 1.96 * size_se),
    ymax = size + 1.96 * size_se), width = 0.03, linewidth = 0.4,
    position = dodge) +
  geom_point(position = dodge, size = 1.8, fill = "white", stroke = 0.6) +
  scale_colour_manual(values = col_vec) +
  scale_shape_manual(values = c(16, 17, 21, 24)) +
  scale_x_continuous(breaks = c(0, 0.2, 0.5)) +
  scale_y_continuous(limits = c(0, NA)) +
  labs(x = expression("Within-note correlation " * rho),
    y = "Empirical size of T1") +
  guides(colour = guide_legend(nrow = 2), shape = guide_legend(nrow = 2)) +
  theme_fig

# (d) T3 rejection in the same cell: size at r = 1, power at r = 0.6
########################################################
t3_lab <- c(binomial = "Binomial", cluster_t = "Cluster t",
  signflip = "Sign-flip")
tab_d <- tab_b[tab_b$n == 100 & tab_b$theta == 0.90 & tab_b$m == 4 &
  tab_b$spu_ratio %in% c(1, 0.6), ]
tab_d$method_f <- factor(t3_lab[tab_d$method], levels = t3_lab)
tab_d$facet <- factor(ifelse(tab_d$spu_ratio == 1, "Size (r = 1)",
  "Power (r = 0.6)"), levels = c("Size (r = 1)", "Power (r = 0.6)"))
tab_nom <- data.frame(facet = factor("Size (r = 1)",
  levels = levels(tab_d$facet)), y = 0.05)

fig_d <- ggplot(tab_d, aes(rho, reject, colour = method_f,
  shape = method_f)) +
  geom_hline(data = tab_nom, aes(yintercept = y), linetype = "dashed",
    colour = ink_muted, linewidth = 0.3) +
  geom_line(position = dodge, linewidth = 0.5) +
  geom_errorbar(aes(ymin = pmax(0, reject - 1.96 * reject_se),
    ymax = pmin(1, reject + 1.96 * reject_se)), width = 0.03,
    linewidth = 0.4, position = dodge) +
  geom_point(position = dodge, size = 1.8, fill = "white", stroke = 0.6) +
  facet_wrap(~facet, scales = "free_y") +
  scale_colour_manual(values = col_vec[1:3]) +
  scale_shape_manual(values = c(16, 17, 21)) +
  scale_x_continuous(breaks = c(0, 0.2, 0.5)) +
  scale_y_continuous(limits = c(0, NA)) +
  labs(x = expression("Within-note correlation " * rho),
    y = "Rejection rate of T3") +
  theme_fig


time_section <- proc.time() - ptm_local
print("Section: 2. Build the four panels")
print(time_section)

# =============================================================================
# 3. Write the figure
# =============================================================================
print("==== 3. Write the figure ====")
ptm_local <- proc.time()

fig_all <- (fig_a + fig_b) / (fig_c + fig_d) +
  plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")")
width_mm <- 180
height_mm <- 150

ggsave(paste0(dir_fig, "fig2_sizepower.pdf"), fig_all, width = width_mm,
  height = height_mm, units = "mm", device = cairo_pdf)
ggsave(paste0(dir_fig, "fig2_sizepower.png"), fig_all, width = width_mm,
  height = height_mm, units = "mm", dpi = 300, device = ragg::agg_png)
ggsave(paste0(dir_fig, "fig2_sizepower.tiff"), fig_all, width = width_mm,
  height = height_mm, units = "mm", dpi = 600, device = ragg::agg_tiff,
  compression = "lzw")


time_section <- proc.time() - ptm_local
print("Section: 3. Write the figure")
print(time_section)

time_total <- proc.time() - ptm1
print("Total time")
print(time_total)
