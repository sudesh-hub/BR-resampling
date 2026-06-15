## highdim_illustration.R
## High-dimensional illustration consistent with manuscript Section 6 and Appendix C

suppressPackageStartupMessages({
  library(data.table)
  library(future.apply)
  library(ggplot2)
  library(patchwork)
})

if (!file.exists("R/br_functions.R")) stop("Input file 'R/br_functions.R' not found.")
source("R/br_functions.R")

if (exists("computepvalue")) {
  compute_pvalue_hd <- computepvalue
} else if (exists("compute_pvalue")) {
  compute_pvalue_hd <- compute_pvalue
} else {
  stop("Neither computepvalue() nor compute_pvalue() was found after sourcing R/br_functions.R.")
}

if (exists("bhrejections")) {
  bh_rejections_hd <- bhrejections
} else if (exists("bh_rejections")) {
  bh_rejections_hd <- bh_rejections
} else {
  stop("Neither bhrejections() nor bh_rejections() was found after sourcing R/br_functions.R.")
}

set.seed(123)
N <- 40
n <- 20
p <- 5000
alpha <- 0.05
b_prime <- 300
S_ba <- 10000
S_mc_benchmark <- 1000000
max_swaps <- 10000
n_cores <- max(1, parallel::detectCores(logical = FALSE) - 1)
future::plan(future::multisession, workers = n_cores)

true_effects <- c(rep(0.8, 10), rep(0, p - 10))
group <- c(rep(0, n), rep(1, n))

X <- matrix(rnorm(p * N), nrow = p, ncol = N)
X[, group == 1] <- X[, group == 1] + true_effects

des <- generate_near_balanced_design(
  N = N,
  n = n,
  b_prime = b_prime,
  max_swaps = max_swaps,
  seed = 123
)

br_time <- system.time({
  gene_pvals_br <- future_apply(X, 1, function(row) {
    compute_pvalue_hd(row, n, design = des$support, method = "BR")
  }, future.seed = TRUE)
})[3]

ba_time <- system.time({
  gene_pvals_ba <- future_apply(X, 1, function(row) {
    compute_pvalue_hd(row, n, method = "BA", S = S_ba)
  }, future.seed = TRUE)
})[3]

bm_time <- system.time({
  gene_pvals_bm <- future_apply(X, 1, function(row) {
    compute_pvalue_hd(row, n, method = "MC_BENCH", S = S_mc_benchmark)
  }, future.seed = TRUE)
})[3]

raw_br <- gene_pvals_br <= alpha
raw_ba <- gene_pvals_ba <= alpha
raw_bm <- gene_pvals_bm <= alpha
bh_br <- bh_rejections_hd(gene_pvals_br, alpha = alpha)
bh_ba <- bh_rejections_hd(gene_pvals_ba, alpha = alpha)
bh_bm <- bh_rejections_hd(gene_pvals_bm, alpha = alpha)

summary_hd <- data.table(
  Method = c("BR", "BA", "Benchmark"),
  Raw_agreement = c(mean(raw_br == raw_bm), mean(raw_ba == raw_bm), 1),
  FDR_agreement = c(mean(bh_br == bh_bm), mean(bh_ba == bh_bm), 1),
  Runtime_seconds = as.numeric(c(br_time, ba_time, bm_time))
)
summary_hd[, Relative_runtime_BA_100 := (Runtime_seconds / Runtime_seconds[Method == "BA"]) * 100]

dir.create("results/highdim", recursive = TRUE, showWarnings = FALSE)
dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)

fwrite(summary_hd, "results/highdim/highdim_summary.csv")
saveRDS(list(
  gene_pvals_br = gene_pvals_br,
  gene_pvals_ba = gene_pvals_ba,
  gene_pvals_benchmark = gene_pvals_bm,
  design = des,
  summary = summary_hd,
  params = list(
    N = N,
    n = n,
    p = p,
    alpha = alpha,
    b_prime = b_prime,
    S_ba = S_ba,
    S_mc_benchmark = S_mc_benchmark,
    max_swaps = max_swaps,
    n_cores = n_cores
  )
), "results/highdim/highdim_results.rds")

viz_dt <- data.table(
  p_br = gene_pvals_br,
  p_bm = gene_pvals_bm,
  Fold_Change = true_effects
)

viz_dt[, `:=`(
  FDR_BR = p.adjust(p_br, method = "BH"),
  FDR_BM = p.adjust(p_bm, method = "BH"),
  reject_BR = FDR_BR < alpha,
  reject_BM = FDR_BM < alpha
)]

viz_dt[, Decision_Agreement := fifelse(
  reject_BR & reject_BM, "Both Reject",
  fifelse(reject_BR & !reject_BM, "BR Reject Only",
          fifelse(!reject_BR & reject_BM, "BM Reject Only", "Neither Reject"))
)]

viz_dt[, Decision_Agreement := factor(
  Decision_Agreement,
  levels = c("Both Reject", "Neither Reject", "BR Reject Only", "BM Reject Only")
)]

viz_dt[, `:=`(
  log10p_BR = -log10(pmax(p_br, .Machine$double.xmin)),
  log10p_BM = -log10(pmax(p_bm, .Machine$double.xmin))
)]

cols_agreement <- c(
  "Both Reject" = "#0072B2",
  "Neither Reject" = "#999999",
  "BR Reject Only" = "#D55E00",
  "BM Reject Only" = "#009E73"
)

volcano_plot <- ggplot(viz_dt, aes(x = Fold_Change, y = log10p_BR, color = Decision_Agreement)) +
  geom_point(size = 1.8, alpha = 0.75, stroke = 0, shape = 16) +
  geom_hline(yintercept = -log10(alpha / p), linetype = "dotted", color = "gray50") +
  geom_vline(xintercept = 0, color = "black", linewidth = 0.5) +
  scale_color_manual(values = cols_agreement) +
  theme_bw(base_size = 12) +
  labs(
    title = "Volcano Agreement: Balanced Resampling vs. Benchmark FDR Control",
    subtitle = "Color indicates agreement in FDR rejection decision (p = 5000 genes)",
    x = "True Mean Difference (Fold Change)",
    y = "-log10(BR p-value)",
    color = "Statistical Agreement"
  ) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold"),
    strip.background = element_rect(fill = "gray95")
  )

ggsave("results/figures/volcano_agreement.png", volcano_plot, width = 11, height = 7, dpi = 300)

p_top <- ggplot(viz_dt, aes(x = Fold_Change, fill = Decision_Agreement, color = Decision_Agreement)) +
  geom_density(alpha = 0.3, linewidth = 0.6) +
  scale_fill_manual(values = cols_agreement) +
  scale_color_manual(values = cols_agreement) +
  theme_void() +
  theme(legend.position = "none")

p_side <- ggplot(viz_dt, aes(y = log10p_BR, fill = Decision_Agreement, color = Decision_Agreement)) +
  geom_density(alpha = 0.3, linewidth = 0.6) +
  scale_fill_manual(values = cols_agreement) +
  scale_color_manual(values = cols_agreement) +
  theme_void() +
  theme(legend.position = "none")

final_plot <- (p_top + plot_spacer()) / (volcano_plot + p_side) +
  plot_layout(widths = c(4, 1), heights = c(1, 4), guides = "collect") &
  theme(legend.position = "bottom")

ggsave("results/figures/volcano_marginal_density.png", final_plot, width = 12, height = 8, dpi = 300)

print(volcano_plot)
print(final_plot)