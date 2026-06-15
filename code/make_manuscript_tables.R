## make_manuscript_tables.R
## Manuscript-ready table generation for the main simulation study

suppressPackageStartupMessages({
  library(data.table)
})

fmt3 <- function(x) sprintf("%.3f", x)
fmt2 <- function(x) sprintf("%.2f", x)
pct1 <- function(x) sprintf("%.1f\\%%", 100 * x)

# Write a LaTeX-style table fragment: columns separated by " & ", rows end with "\\"
save_latex_fragment <- function(dt, filename) {
  dt_copy <- copy(dt)
  last_col <- names(dt_copy)[ncol(dt_copy)]
  # Append '\\' to the last column for LaTeX
  dt_copy[, (last_col) := paste0(get(last_col), " \\\\")]
  # Build lines: header + one line per row
  header_line <- paste(names(dt_copy), collapse = " & ")
  body_lines <- apply(dt_copy, 1L, function(row) paste(row, collapse = " & "))
  lines <- c(header_line, body_lines)
  writeLines(lines, filename)
}

sim_file <- "results/main/simulation_results.rds"
if (!file.exists(sim_file)) {
  stop(sprintf("Input file '%s' not found.", sim_file))
}

x <- as.data.table(readRDS(sim_file))
required_cols <- c(
  "dist", "regime", "hyp", "is_null",
  "p_br", "p_mc", "p_ba", "p_ref", "ref_method",
  "time_br", "time_mc", "time_ba", "time_ref"
)
missing_cols <- setdiff(required_cols, names(x))
if (length(missing_cols) > 0) {
  stop(sprintf(
    "Missing required columns in simulation results: %s",
    paste(missing_cols, collapse = ", ")
  ))
}

alpha <- 0.05
x[, `:=`(
  reject_br = p_br < alpha,
  reject_mc = p_mc < alpha,
  reject_ba = p_ba < alpha,
  reject_ref = p_ref < alpha
)]

x[, `:=`(
  agree_br = reject_br == reject_ref,
  agree_mc = reject_mc == reject_ref,
  agree_ba = reject_ba == reject_ref
)]

rmse <- x[is_null == TRUE, .(
  BR = sqrt(mean((p_br - p_ref)^2, na.rm = TRUE)),
  MC = sqrt(mean((p_mc - p_ref)^2, na.rm = TRUE)),
  BA = sqrt(mean((p_ba - p_ref)^2, na.rm = TRUE)),
  N_replicates = .N
), by = .(dist, regime)]

power <- x[is_null == FALSE, .(
  BR = mean(reject_br, na.rm = TRUE),
  MC = mean(reject_mc, na.rm = TRUE),
  BA = mean(reject_ba, na.rm = TRUE),
  N_replicates = .N
), by = .(dist, regime)]

runtime <- x[, .(
  BR = mean(time_br, na.rm = TRUE),
  MC = mean(time_mc, na.rm = TRUE),
  BA = mean(time_ba, na.rm = TRUE),
  REF = mean(time_ref, na.rm = TRUE),
  N_replicates = .N
), by = .(dist, regime)]

reference_mix <- x[, .(
  N_replicates = .N
), by = .(dist, regime, ref_method)]

agreement <- x[, .(
  BR = mean(agree_br, na.rm = TRUE),
  MC = mean(agree_mc, na.rm = TRUE),
  BA = mean(agree_ba, na.rm = TRUE),
  N_replicates = .N
), by = .(dist, regime, hyp)]

dir.create("results/main", recursive = TRUE, showWarnings = FALSE)

fwrite(rmse,           "results/main/manuscript_rmse.csv")
fwrite(power,          "results/main/manuscript_power.csv")
fwrite(runtime,        "results/main/manuscript_runtime.csv")
fwrite(reference_mix,  "results/main/manuscript_reference_mix.csv")
fwrite(agreement,      "results/main/manuscript_agreement.csv")

rmse_tex <- copy(rmse)
rmse_tex[, c("BR", "MC", "BA") := lapply(.SD, fmt3), .SDcols = c("BR", "MC", "BA")]

power_tex <- copy(power)
power_tex[, c("BR", "MC", "BA") := lapply(.SD, pct1), .SDcols = c("BR", "MC", "BA")]

runtime_tex <- copy(runtime)
runtime_tex[, c("BR", "MC", "BA", "REF") := lapply(.SD, fmt2),
            .SDcols = c("BR", "MC", "BA", "REF")]

agreement_tex <- copy(agreement)
agreement_tex[, c("BR", "MC", "BA") := lapply(.SD, pct1),
              .SDcols = c("BR", "MC", "BA")]

save_latex_fragment(rmse_tex,      "results/main/manuscript_rmse.tex")
save_latex_fragment(power_tex,     "results/main/manuscript_power.tex")
save_latex_fragment(runtime_tex,   "results/main/manuscript_runtime.tex")
save_latex_fragment(agreement_tex, "results/main/manuscript_agreement.tex")

highdim_candidates <- c(
  "results/highdim/highdim_results.rds",
  "results/highdim_results.rds",
  "results/main/highdim_results.rds"
)
highdim_file <- highdim_candidates[file.exists(highdim_candidates)][1]

if (!is.na(highdim_file)) {
  hd_raw <- readRDS(highdim_file)
  hd <- if (is.list(hd_raw) && "summary" %in% names(hd_raw)) {
    as.data.table(hd_raw$summary)
  } else {
    as.data.table(hd_raw)
  }

  fwrite(hd, "results/main/manuscript_highdim.csv")

  hd_tex <- copy(hd)
  cols_pct <- intersect(c("Raw_agreement", "FDR_agreement"), names(hd_tex))
  if (length(cols_pct) > 0) {
    hd_tex[, (cols_pct) := lapply(.SD, function(z)
      ifelse(is.na(z), "--", pct1(z))), .SDcols = cols_pct]
  }
  if ("Runtime_seconds" %in% names(hd_tex)) {
    hd_tex[, Runtime_seconds := ifelse(is.na(Runtime_seconds),
                                       "--", fmt2(Runtime_seconds))]
  }
  if ("Relative_runtime_BA_100" %in% names(hd_tex)) {
    hd_tex[, Relative_runtime_BA_100 := ifelse(
      is.na(Relative_runtime_BA_100),
      "--",
      sprintf("%.1f\\%%", Relative_runtime_BA_100)
    )]
  }

  save_latex_fragment(hd_tex, "results/main/manuscript_highdim.tex")
}

cat("Manuscript-ready table files created in results/main/.\n")