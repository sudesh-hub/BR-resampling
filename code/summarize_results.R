## scripts/summarize_results.R
## Post-simulation summaries for the main two-sample t-statistic study

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

input_file <- "results/main/simulation_results.rds"
if (!file.exists(input_file)) stop(sprintf("Input file '%s' not found.", input_file))

results <- as.data.table(readRDS(input_file))

required_cols <- c(
  "dist", "regime", "hyp", "is_null", "N", "n", "replicate",
  "p_br", "p_mc", "p_ba", "p_ref", "ref_method",
  "time_br", "time_mc", "time_ba", "time_ref",
  "discrepancy", "max_first_dev", "max_second_dev",
  "br_support_size", "total_allocations"
)

missing_cols <- setdiff(required_cols, names(results))
if (length(missing_cols) > 0) {
  stop(sprintf("Missing required columns: %s", paste(missing_cols, collapse = ", ")))
}

alpha <- 0.05
results[, `:=`(
  reject_br = p_br < alpha,
  reject_mc = p_mc < alpha,
  reject_ba = p_ba < alpha,
  reject_ref = p_ref < alpha,
  agree_br = (p_br < alpha) == (p_ref < alpha),
  agree_mc = (p_mc < alpha) == (p_ref < alpha),
  agree_ba = (p_ba < alpha) == (p_ref < alpha)
)]

rmse_tbl <- results[is_null == TRUE, .(
  BR = sqrt(mean((p_br - p_ref)^2, na.rm = TRUE)),
  MC = sqrt(mean((p_mc - p_ref)^2, na.rm = TRUE)),
  BA = sqrt(mean((p_ba - p_ref)^2, na.rm = TRUE)),
  N_replicates = .N
), by = .(dist, regime)]

power_tbl <- results[is_null == FALSE, .(
  BR = mean(reject_br, na.rm = TRUE),
  MC = mean(reject_mc, na.rm = TRUE),
  BA = mean(reject_ba, na.rm = TRUE),
  N_replicates = .N
), by = .(dist, regime)]

agreement_tbl <- results[, .(
  BR = mean(agree_br, na.rm = TRUE),
  MC = mean(agree_mc, na.rm = TRUE),
  BA = mean(agree_ba, na.rm = TRUE),
  N_replicates = .N
), by = .(dist, regime, hyp)]

runtime_tbl <- results[, .(
  BR = mean(time_br, na.rm = TRUE),
  MC = mean(time_mc, na.rm = TRUE),
  BA = mean(time_ba, na.rm = TRUE),
  REF = mean(time_ref, na.rm = TRUE),
  N_replicates = .N
), by = .(dist, regime)]

design_tbl <- results[, .(
  mean_discrepancy = mean(discrepancy, na.rm = TRUE),
  mean_max_first_dev = mean(max_first_dev, na.rm = TRUE),
  mean_max_second_dev = mean(max_second_dev, na.rm = TRUE),
  mean_br_support_size = mean(br_support_size, na.rm = TRUE),
  min_br_support_size = min(br_support_size, na.rm = TRUE),
  max_br_support_size = max(br_support_size, na.rm = TRUE),
  N_replicates = .N
), by = .(dist, regime)]

reference_tbl <- results[, .(N_replicates = .N), by = .(dist, regime, ref_method)]

allocation_tbl <- results[, .(
  mean_total_allocations = mean(total_allocations, na.rm = TRUE),
  min_total_allocations = min(total_allocations, na.rm = TRUE),
  max_total_allocations = max(total_allocations, na.rm = TRUE),
  N_replicates = .N
), by = .(dist, regime)]

dir.create("results/main", recursive = TRUE, showWarnings = FALSE)

fwrite(rmse_tbl, "results/main/table_rmse.csv")
fwrite(power_tbl, "results/main/table_power.csv")
fwrite(agreement_tbl, "results/main/table_agreement.csv")
fwrite(runtime_tbl, "results/main/table_runtime.csv")
fwrite(design_tbl, "results/main/table_design_diagnostics.csv")
fwrite(reference_tbl, "results/main/table_reference_method_counts.csv")
fwrite(allocation_tbl, "results/main/table_allocation_space_summary.csv")

rmse_long <- melt(
  rmse_tbl,
  id.vars = c("dist", "regime", "N_replicates"),
  measure.vars = c("BR", "MC", "BA"),
  variable.name = "method",
  value.name = "rmse"
)

power_long <- melt(
  power_tbl,
  id.vars = c("dist", "regime", "N_replicates"),
  measure.vars = c("BR", "MC", "BA"),
  variable.name = "method",
  value.name = "power"
)

viz_df <- merge(rmse_long, power_long, by = c("dist", "regime", "method"), all = TRUE)

viz_df[, method := factor(
  method,
  levels = c("BR", "MC", "BA"),
  labels = c("Balanced BR", "Monte Carlo MC", "Bootstrap BA")
)]
viz_df[, regime := factor(regime, levels = c("small", "moderate", "large"))]

p <- ggplot(viz_df, aes(x = rmse, y = power, color = method, shape = regime)) +
  geom_point(size = 3.2, alpha = 0.9) +
  geom_path(aes(group = interaction(dist, method)), alpha = 0.25, linewidth = 0.5) +
  facet_wrap(~ dist, scales = "free") +
  scale_color_manual(values = c(
    "Balanced BR" = "#0072B2",
    "Monte Carlo MC" = "#D55E00",
    "Bootstrap BA" = "#009E73"
  )) +
  labs(
    title = "Accuracy vs power across resampling methods",
    x = "RMSE relative to reference p-value",
    y = "Power at alpha = 0.05",
    color = "Method",
    shape = "Regime"
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "bottom",
    strip.background = element_rect(fill = "gray95"),
    plot.title = element_text(face = "bold")
  )

dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)
ggsave("results/figures/fig_power_vs_rmse.png", p, width = 11, height = 7, dpi = 300)