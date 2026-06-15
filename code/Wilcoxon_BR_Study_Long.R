## =============================================================================
## scripts/Wilcoxon_BR_Study_Long.R
##
## Full Wilcoxon rank-sum robustness study for RStudio.
## Corrected version:
##   1) fixes scalar coercion for parallel workers
##   2) fixes grid column extraction by quoted names
## =============================================================================

rm(list = ls())
gc()

suppressPackageStartupMessages({
  library(data.table)
})

have_future <- requireNamespace("future", quietly = TRUE) &&
  requireNamespace("future.apply", quietly = TRUE) &&
  requireNamespace("progressr", quietly = TRUE)

if (have_future) {
  suppressPackageStartupMessages({
    library(future)
    library(future.apply)
    library(progressr)
  })
}

## -----------------------------------------------------------------------------
## 1. Parallel setup
## -----------------------------------------------------------------------------

n_cores <- max(1L, parallel::detectCores(logical = FALSE) - 1L)
use_parallel <- have_future

if (use_parallel) {
  future::plan(future::multisession, workers = n_cores)
  progressr::handlers(global = TRUE)
} else {
  message("Parallel packages not available; running sequentially.")
}

## -----------------------------------------------------------------------------
## 2. Global settings
## -----------------------------------------------------------------------------

set.seed(20260423)

quick_test <- FALSE

if (quick_test) {
  S_mc <- 1000L
  S_benchmark <- 5000L
  S_ba <- 1000L
  M_reps <- 5L
  max_exact_alloc <- 5000
} else {
  S_mc <- 5000L
  S_benchmark <- 50000L
  S_ba <- 5000L
  M_reps <- 200L
  max_exact_alloc <- 50000
}

alpha <- 0.05

regimes <- list(
  Small    = list(N = 16L, n = 8L),
  Moderate = list(N = 26L, n = 13L),
  Large    = list(N = 36L, n = 18L)
)

distribution_names <- c("Normal", "Uniform", "Exponential", "Gamma", "Lognormal")

## -----------------------------------------------------------------------------
## 3. Distribution dictionary
## -----------------------------------------------------------------------------

distributions <- list(
  Normal = list(
    r = function(n) {
      n <- as.integer(n)[1]
      if (is.na(n) || n < 0L) stop("Invalid n in Normal generator: ", n)
      rnorm(n, mean = 0, sd = 1)
    },
    sd = 1
  ),
  Uniform = list(
    r = function(n) {
      n <- as.integer(n)[1]
      if (is.na(n) || n < 0L) stop("Invalid n in Uniform generator: ", n)
      runif(n, min = 0, max = 1)
    },
    sd = sqrt(1 / 12)
  ),
  Exponential = list(
    r = function(n) {
      n <- as.integer(n)[1]
      if (is.na(n) || n < 0L) stop("Invalid n in Exponential generator: ", n)
      rexp(n, rate = 1)
    },
    sd = 1
  ),
  Gamma = list(
    r = function(n) {
      n <- as.integer(n)[1]
      if (is.na(n) || n < 0L) stop("Invalid n in Gamma generator: ", n)
      rgamma(n, shape = 2, rate = 1)
    },
    sd = sqrt(2)
  ),
  Lognormal = list(
    r = function(n) {
      n <- as.integer(n)[1]
      if (is.na(n) || n < 0L) stop("Invalid n in Lognormal generator: ", n)
      rlnorm(n, meanlog = 0, sdlog = 1)
    },
    sd = sqrt((exp(1) - 1) * exp(1))
  )
)

## -----------------------------------------------------------------------------
## 4. Fixed BR supports
## -----------------------------------------------------------------------------

get_br_support_16_8_bibd <- function() {
  list(
    c(1,6,7,8,9,14,15,16), c(2,3,4,5,10,11,12,13),
    c(1,4,5,8,9,12,13,16), c(2,3,6,7,10,11,14,15),
    c(1,2,3,8,9,10,11,16), c(4,5,6,7,12,13,14,15),
    c(1,3,5,7,9,11,13,15), c(2,4,6,8,10,12,14,16),
    c(1,2,4,7,9,10,12,15), c(3,5,6,8,11,13,14,16),
    c(1,2,5,6,9,10,13,14), c(3,4,7,8,11,12,15,16),
    c(1,3,4,6,9,11,12,14), c(2,5,7,8,10,13,15,16),
    c(1,10,11,12,13,14,15,16), c(2,3,4,5,6,7,8,9),
    c(1,2,3,4,5,14,15,16), c(6,7,8,9,10,11,12,13),
    c(1,2,3,6,7,12,13,16), c(4,5,8,9,10,11,14,15),
    c(1,4,5,6,7,10,11,16), c(2,3,8,9,12,13,14,15),
    c(1,2,4,6,8,11,13,15), c(3,5,7,9,10,12,14,16),
    c(1,3,5,6,8,10,12,15), c(2,4,7,9,11,13,14,16),
    c(1,3,4,7,8,10,13,14), c(2,5,6,9,11,12,15,16),
    c(1,2,5,7,8,11,12,14), c(3,4,6,9,10,13,15,16)
  )
}

get_br_support_26_13_fixed <- function() {
  list(
    c(23,4,1,9,5,3,18,22,6,7,11,25,20),
    c(12,7,22,16,17,18,15,4,8,25,23,6,5),
    c(3,24,2,12,14,11,17,7,13,20,25,4,9),
    c(17,19,13,5,18,12,16,26,2,11,20,8,7),
    c(5,16,25,11,26,14,24,6,17,22,15,19,4),
    c(19,14,21,24,5,15,23,2,7,11,3,18,8),
    c(4,5,25,8,14,1,12,18,21,2,24,9,10),
    c(6,12,26,16,22,21,5,25,8,9,19,20,2),
    c(13,26,11,21,10,14,25,2,17,23,7,16,6),
    c(17,2,19,23,7,24,16,22,5,6,9,10,1),
    c(26,8,24,10,3,17,22,21,20,16,7,25,15),
    c(23,21,8,6,1,4,17,10,20,5,14,11,24),
    c(7,23,6,18,3,5,20,24,14,17,26,12,10),
    c(1,17,15,3,13,9,6,11,21,12,26,5,7),
    c(18,7,20,25,24,8,10,11,22,1,9,26,15),
    c(22,1,4,19,15,24,3,7,14,13,2,6,12),
    c(9,13,12,4,16,3,7,20,10,23,6,24,8),
    c(20,11,8,5,4,26,13,1,12,15,23,10,22),
    c(26,23,4,21,7,6,13,9,15,8,18,11,19),
    c(24,15,11,4,18,23,19,3,2,20,22,17,26),
    c(14,19,6,20,16,22,18,12,1,8,21,7,3),
    c(21,17,16,1,2,11,4,15,24,18,8,12,9),
    c(25,24,23,15,12,10,16,19,18,3,5,13,21),
    c(25,13,23,18,1,7,15,14,4,21,19,20,16),
    c(16,20,2,14,24,9,23,21,12,15,11,26,6),
    c(15,3,10,7,9,5,19,26,24,16,8,14,4),
    c(11,18,10,26,19,2,1,12,16,6,4,3,20),
    c(21,12,7,11,4,20,25,5,19,17,24,16,1),
    c(24,14,19,26,9,18,12,23,22,1,13,5,17),
    c(14,26,9,3,10,7,22,19,11,4,25,21,12),
    c(1,20,12,14,21,26,5,3,15,25,2,23,22),
    c(4,5,24,6,22,2,7,9,25,26,21,13,18),
    c(15,1,14,25,20,19,2,8,6,9,17,7,23),
    c(8,11,15,19,6,10,9,5,25,2,12,3,23),
    c(5,22,16,2,13,8,14,20,3,4,15,9,17),
    c(7,10,5,17,11,21,2,8,3,19,22,1,13),
    c(3,6,18,17,11,16,21,15,9,22,10,1,14),
    c(10,2,17,23,21,12,8,4,7,26,14,22,18),
    c(10,25,18,20,2,6,21,17,13,5,4,15,3),
    c(8,4,3,13,6,22,26,16,23,21,1,2,24),
    c(2,6,13,22,25,19,8,24,10,18,20,14,11),
    c(2,15,7,24,23,25,10,13,11,12,1,22,16),
    c(13,8,20,9,23,1,26,17,19,24,3,21,25),
    c(19,16,3,1,17,4,9,25,26,10,18,23,2),
    c(16,9,22,8,25,13,14,18,23,3,12,17,11),
    c(12,25,14,15,26,17,4,10,1,13,6,8,19),
    c(9,18,26,7,20,15,1,14,5,10,16,2,13),
    c(6,3,5,13,8,25,11,1,26,24,16,18,14),
    c(20,9,17,22,12,13,24,6,18,19,10,15,21),
    c(22,21,9,10,19,20,11,23,16,14,13,4,5)
  )
}

construct_br_support_safe <- function(N, n, m = 300L, seed = NULL) {
  N <- as.integer(N)[1]
  n <- as.integer(n)[1]
  m <- as.integer(m)[1]

  if (!is.null(seed)) set.seed(seed)

  total_alloc <- choose(N, n)

  if (total_alloc <= m) {
    cmb <- combn(N, n)
    return(split(cmb, col(cmb)))
  }

  support_list <- vector("list", m)
  seen <- new.env(hash = TRUE, parent = emptyenv())
  k <- 1L

  while (k <= m) {
    alloc <- sort(sample.int(N, n, replace = FALSE))
    key <- paste(alloc, collapse = ",")

    if (!exists(key, envir = seen, inherits = FALSE)) {
      assign(key, TRUE, envir = seen)
      support_list[[k]] <- alloc
      k <- k + 1L
    }
  }

  support_list
}

get_br_support_fixed <- function(N, n) {
  N <- as.integer(N)[1]
  n <- as.integer(n)[1]

  if (N == 16L && n == 8L) {
    get_br_support_16_8_bibd()
  } else if (N == 26L && n == 13L) {
    get_br_support_26_13_fixed()
  } else if (N == 36L && n == 18L) {
    construct_br_support_safe(N, n, m = 300L, seed = 42L + N)
  } else {
    stop("No fixed BR support defined for N = ", N, ", n = ", n)
  }
}

## -----------------------------------------------------------------------------
## 5. Helper functions
## -----------------------------------------------------------------------------

make_key <- function(N, n) paste0(as.integer(N)[1], "_", as.integer(n)[1])

support_list_to_matrix <- function(support_list) {
  out <- do.call(cbind, lapply(support_list, function(x) as.integer(x)))
  storage.mode(out) <- "integer"
  out
}

sample_alloc_matrix <- function(N, n, S) {
  N <- as.integer(N)[1]
  n <- as.integer(n)[1]
  S <- as.integer(S)[1]
  replicate(S, sort(sample.int(N, n, replace = FALSE)))
}

wilcox_stats_matrix <- function(ranks, alloc_mat) {
  colSums(matrix(ranks[alloc_mat], nrow = nrow(alloc_mat)))
}

perm_pvalue_centered <- function(stats, obs_stat) {
  center <- mean(stats)
  mean(abs(stats - center) >= abs(obs_stat - center) - 1e-15)
}

rmse_fun  <- function(a, b) sqrt(mean((a - b)^2))
mae_fun   <- function(a, b) mean(abs(a - b))
agree_fun <- function(a, b, alpha = 0.05) mean((a <= alpha) == (b <= alpha))

## -----------------------------------------------------------------------------
## 6. Bootstrap p-value engine
## -----------------------------------------------------------------------------

bootstrap_pvalue_wilcox <- function(x1, x2, S_ba, n1, n2) {
  S_ba <- as.integer(S_ba)[1]
  n1 <- as.integer(n1)[1]
  n2 <- as.integer(n2)[1]

  x_all <- c(x1, x2)
  ranks_obs <- rank(x_all, ties.method = "average")
  obs_stat <- sum(ranks_obs[seq_len(n1)])

  W_ba <- numeric(S_ba)

  for (s in seq_len(S_ba)) {
    xb1 <- sample(x1, n1, replace = TRUE)
    xb2 <- sample(x2, n2, replace = TRUE)
    W_ba[s] <- sum(rank(c(xb1, xb2), ties.method = "average")[seq_len(n1)])
  }

  center <- mean(W_ba)
  mean(abs(W_ba - center) >= abs(obs_stat - center) - 1e-15)
}

## -----------------------------------------------------------------------------
## 7. One simulation replicate
## -----------------------------------------------------------------------------

simulate_one_wilcox <- function(N, n1, dist_name, hyp, br_entry,
                                S_mc, S_benchmark, S_ba, max_exact_alloc) {
  N  <- as.integer(N)[1]
  n1 <- as.integer(n1)[1]
  dist_name <- as.character(dist_name)[1]
  hyp <- as.character(hyp)[1]

  if (is.na(N) || is.na(n1) || n1 <= 0L || N <= n1) {
    stop("Invalid N/n1 passed to simulate_one_wilcox(): N=", N, ", n1=", n1)
  }

  if (!dist_name %in% names(distributions)) {
    stop("Unknown distribution: ", dist_name)
  }

  if (!hyp %in% c("H0", "H1")) {
    stop("Unknown hypothesis label: ", hyp)
  }

  n2 <- N - n1
  dist_obj <- distributions[[dist_name]]
  sigma <- as.numeric(dist_obj$sd)[1]

  x1 <- dist_obj$r(n1)
  x2 <- dist_obj$r(n2)

  if (hyp == "H1") {
    x2 <- x2 + 0.5 * sigma
  }

  x <- c(x1, x2)
  ranks <- rank(x, ties.method = "average")
  obs_stat <- sum(ranks[seq_len(n1)])

  p_br <- perm_pvalue_centered(
    wilcox_stats_matrix(ranks, br_entry$br_mat),
    obs_stat
  )

  mc_mat <- sample_alloc_matrix(N, n1, S_mc)
  p_mc <- perm_pvalue_centered(
    wilcox_stats_matrix(ranks, mc_mat),
    obs_stat
  )

  p_ba <- bootstrap_pvalue_wilcox(x1, x2, S_ba, n1, n2)

  total_comb <- choose(N, n1)

  if (total_comb <= max_exact_alloc) {
    bench_mat <- combn(N, n1)
    benchmark_type <- "exact"
  } else {
    bench_mat <- sample_alloc_matrix(N, n1, S_benchmark)
    benchmark_type <- "high_precision_mc"
  }

  p_ex <- perm_pvalue_centered(
    wilcox_stats_matrix(ranks, bench_mat),
    obs_stat
  )

  list(
    p_ex = p_ex,
    p_br = p_br,
    p_mc = p_mc,
    p_ba = p_ba,
    benchmark_type = benchmark_type,
    benchmark_size = ncol(bench_mat),
    br_size = br_entry$br_size
  )
}

## -----------------------------------------------------------------------------
## 8. Build fixed BR cache
## -----------------------------------------------------------------------------

support_cache <- list()
diagnostics_list <- list()

for (reg_name in names(regimes)) {
  N <- as.integer(regimes[[reg_name]]$N)[1]
  n <- as.integer(regimes[[reg_name]]$n)[1]
  key <- make_key(N, n)

  support_list <- get_br_support_fixed(N, n)
  br_mat <- support_list_to_matrix(support_list)

  support_cache[[key]] <- list(
    br_list = support_list,
    br_mat = br_mat,
    br_size = ncol(br_mat),
    regime = reg_name
  )

  diagnostics_list[[key]] <- data.table(
    Regime = reg_name,
    N = N,
    n = n,
    BR_Size = ncol(br_mat)
  )
}

diagnostics_dt <- rbindlist(diagnostics_list)

## -----------------------------------------------------------------------------
## 9. Simulation grid
## -----------------------------------------------------------------------------

grid <- rbindlist(lapply(names(regimes), function(reg_name) {
  N <- as.integer(regimes[[reg_name]]$N)[1]
  n <- as.integer(regimes[[reg_name]]$n)[1]

  as.data.table(expand.grid(
    Distribution = distribution_names,
    Hypothesis = c("H0", "H1"),
    Replicate = seq_len(M_reps),
    stringsAsFactors = FALSE
  ))[, `:=`(
    Regime = reg_name,
    N = N,
    n = n,
    key_id = make_key(N, n)
  )]
}))

setcolorder(grid, c("Distribution", "Regime", "N", "n", "Hypothesis", "Replicate", "key_id"))

## -----------------------------------------------------------------------------
## 10. Run study
## -----------------------------------------------------------------------------

t_start <- Sys.time()

if (use_parallel) {
  progressr::with_progress({
    p <- progressr::progressor(steps = nrow(grid))

    results_list <- future.apply::future_lapply(seq_len(nrow(grid)), function(i) {
      p()

      N_i    <- as.integer(grid[i, "N"])
      n_i    <- as.integer(grid[i, "n"])
      dist_i <- as.character(grid[i, "Distribution"])
      hyp_i  <- as.character(grid[i, "Hypothesis"])
      key_i  <- as.character(grid[i, "key_id"])

      simulate_one_wilcox(
        N = N_i,
        n1 = n_i,
        dist_name = dist_i,
        hyp = hyp_i,
        br_entry = support_cache[[key_i]],
        S_mc = S_mc,
        S_benchmark = S_benchmark,
        S_ba = S_ba,
        max_exact_alloc = max_exact_alloc
      )
    }, future.seed = TRUE)
  })
} else {
  results_list <- lapply(seq_len(nrow(grid)), function(i) {
    N_i    <- as.integer(grid[i, "N"])
    n_i    <- as.integer(grid[i, "n"])
    dist_i <- as.character(grid[i, "Distribution"])
    hyp_i  <- as.character(grid[i, "Hypothesis"])
    key_i  <- as.character(grid[i, "key_id"])

    simulate_one_wilcox(
      N = N_i,
      n1 = n_i,
      dist_name = dist_i,
      hyp = hyp_i,
      br_entry = support_cache[[key_i]],
      S_mc = S_mc,
      S_benchmark = S_benchmark,
      S_ba = S_ba,
      max_exact_alloc = max_exact_alloc
    )
  })
}

t_end <- Sys.time()

## -----------------------------------------------------------------------------
## 11. Final results table
## -----------------------------------------------------------------------------

results_dt <- cbind(grid, rbindlist(results_list))

results_dt[, `:=`(
  rej_br = as.integer(p_br <= alpha),
  rej_mc = as.integer(p_mc <= alpha),
  rej_ba = as.integer(p_ba <= alpha),
  rej_ex = as.integer(p_ex <= alpha)
)]

## -----------------------------------------------------------------------------
## 12. Summary tables
## -----------------------------------------------------------------------------

summary_by_dist <- results_dt[, .(
  BR_RMSE  = rmse_fun(p_br, p_ex),
  MC_RMSE  = rmse_fun(p_mc, p_ex),
  BA_RMSE  = rmse_fun(p_ba, p_ex),
  BR_MAE   = mae_fun(p_br, p_ex),
  MC_MAE   = mae_fun(p_mc, p_ex),
  BA_MAE   = mae_fun(p_ba, p_ex),
  BR_Agree = agree_fun(p_br, p_ex, alpha),
  MC_Agree = agree_fun(p_mc, p_ex, alpha),
  BA_Agree = agree_fun(p_ba, p_ex, alpha),
  Power_BR = mean(rej_br[Hypothesis == "H1"]),
  Power_MC = mean(rej_mc[Hypothesis == "H1"]),
  Power_BA = mean(rej_ba[Hypothesis == "H1"]),
  Power_EX = mean(rej_ex[Hypothesis == "H1"]),
  TypeI_BR = mean(rej_br[Hypothesis == "H0"]),
  TypeI_MC = mean(rej_mc[Hypothesis == "H0"]),
  TypeI_BA = mean(rej_ba[Hypothesis == "H0"]),
  TypeI_EX = mean(rej_ex[Hypothesis == "H0"])
), by = .(Distribution)][order(Distribution)]

summary_by_regime <- results_dt[, .(
  BR_RMSE  = rmse_fun(p_br, p_ex),
  MC_RMSE  = rmse_fun(p_mc, p_ex),
  BA_RMSE  = rmse_fun(p_ba, p_ex),
  BR_MAE   = mae_fun(p_br, p_ex),
  MC_MAE   = mae_fun(p_mc, p_ex),
  BA_MAE   = mae_fun(p_ba, p_ex),
  BR_Agree = agree_fun(p_br, p_ex, alpha),
  MC_Agree = agree_fun(p_mc, p_ex, alpha),
  BA_Agree = agree_fun(p_ba, p_ex, alpha),
  Power_BR = mean(rej_br[Hypothesis == "H1"]),
  Power_MC = mean(rej_mc[Hypothesis == "H1"]),
  Power_BA = mean(rej_ba[Hypothesis == "H1"]),
  Power_EX = mean(rej_ex[Hypothesis == "H1"]),
  TypeI_BR = mean(rej_br[Hypothesis == "H0"]),
  TypeI_MC = mean(rej_mc[Hypothesis == "H0"]),
  TypeI_BA = mean(rej_ba[Hypothesis == "H0"]),
  TypeI_EX = mean(rej_ex[Hypothesis == "H0"])
), by = .(Regime)][order(Regime)]

summary_overall <- results_dt[, .(
  Distribution = "Overall",
  BR_RMSE  = rmse_fun(p_br, p_ex),
  MC_RMSE  = rmse_fun(p_mc, p_ex),
  BA_RMSE  = rmse_fun(p_ba, p_ex),
  BR_MAE   = mae_fun(p_br, p_ex),
  MC_MAE   = mae_fun(p_mc, p_ex),
  BA_MAE   = mae_fun(p_ba, p_ex),
  BR_Agree = agree_fun(p_br, p_ex, alpha),
  MC_Agree = agree_fun(p_mc, p_ex, alpha),
  BA_Agree = agree_fun(p_ba, p_ex, alpha),
  Power_BR = mean(rej_br[Hypothesis == "H1"]),
  Power_MC = mean(rej_mc[Hypothesis == "H1"]),
  Power_BA = mean(rej_ba[Hypothesis == "H1"]),
  Power_EX = mean(rej_ex[Hypothesis == "H1"]),
  TypeI_BR = mean(rej_br[Hypothesis == "H0"]),
  TypeI_MC = mean(rej_mc[Hypothesis == "H0"]),
  TypeI_BA = mean(rej_ba[Hypothesis == "H0"]),
  TypeI_EX = mean(rej_ex[Hypothesis == "H0"])
)]

summary_manuscript <- rbind(summary_by_dist, summary_overall, fill = TRUE)

## -----------------------------------------------------------------------------
## 13. Export
## -----------------------------------------------------------------------------

dir.create("results/wilcoxon", recursive = TRUE, showWarnings = FALSE)

fwrite(results_dt, "results/wilcoxon/raw_data.csv")
fwrite(summary_by_dist, "results/wilcoxon/summary_by_distribution.csv")
fwrite(summary_by_regime, "results/wilcoxon/summary_by_regime.csv")
fwrite(summary_overall, "results/wilcoxon/summary_overall.csv")
fwrite(summary_manuscript, "results/wilcoxon/summary_manuscript.csv")
fwrite(diagnostics_dt, "results/wilcoxon/design_diagnostics.csv")

writeLines(
  c(
    "Wilcoxon BR robustness study runtime log",
    paste("Start:", as.character(t_start)),
    paste("End:", as.character(t_end)),
    paste("Runtime (mins):", round(as.numeric(difftime(t_end, t_start, units = "mins")), 2)),
    paste("Parallel:", use_parallel),
    paste("Workers:", if (use_parallel) n_cores else 1),
    paste("Quick_test:", quick_test),
    paste("M_reps:", M_reps),
    paste("S_mc:", S_mc),
    paste("S_benchmark:", S_benchmark),
    paste("S_ba:", S_ba),
    paste("max_exact_alloc:", max_exact_alloc)
  ),
  con = "results/wilcoxon/runtime.txt"
)

cat("\nSUCCESS: All Wilcoxon files saved to 'results/wilcoxon/'.\n\n")
cat("Summary by distribution:\n")
print(summary_by_dist)

cat("\nSummary overall:\n")
print(summary_overall)