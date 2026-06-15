## scripts/run_simulations.R
## Main two-sample t-statistic simulation study

source("R/br_functions.R")

set.seed(123)

# -----------------------------------------------------------------------------
# Simulation settings
# -----------------------------------------------------------------------------

distributions <- c("normal", "uniform", "exponential", "gamma", "lognormal")
regimes <- list(
  small = 11:20,
  moderate = 21:30,
  large = 31:40
)

quick_test <- FALSE

if (quick_test) {
  M <- 5
  b_prime <- 50
  S_mc <- 1000
  S_ba <- 1000
  S_mc_bench <- 10000
  max_swaps_br <- 500
} else {
  M <- 200
  b_prime <- 300
  S_mc <- 10000
  S_ba <- 10000
  S_mc_bench <- 1000000
  max_swaps_br <- 5000
}

fp_exact_max <- 50000
n_restarts_br <- 10
n_candidates_br <- 20

# -----------------------------------------------------------------------------
# Data-generating models
# -----------------------------------------------------------------------------

rgen <- function(dist, n, shift = 0) {
  switch(
    dist,
    normal = rnorm(n, mean = shift, sd = 1),
    uniform = runif(n, min = 0, max = 1) + shift,
    exponential = rexp(n, rate = 1) + shift,
    gamma = rgamma(n, shape = 2, rate = 1) + shift,
    lognormal = rlnorm(n, meanlog = 0, sdlog = 1) + shift
  )
}

sd_pop <- function(dist) {
  switch(
    dist,
    normal = 1,
    uniform = sqrt(1 / 12),
    exponential = 1,
    gamma = sqrt(2),
    lognormal = sqrt(exp(1) * (exp(1) - 1))
  )
}

# -----------------------------------------------------------------------------
# Cached BR support construction
# -----------------------------------------------------------------------------

br_cache <- new.env(parent = emptyenv())

get_br_design_cached <- function(N, n, b_prime,
                                 max_swaps_br,
                                 n_restarts_br,
                                 n_candidates_br,
                                 seed_base = 1000) {
  key <- paste(
    N, n, b_prime, max_swaps_br, n_restarts_br, n_candidates_br,
    sep = "_"
  )

  if (!exists(key, envir = br_cache, inherits = FALSE)) {
    message(sprintf("Constructing BR support for N=%d, n=%d", N, n))

    des <- generate_near_balanced_design(
      N = N,
      n = n,
      b_prime = b_prime,
      max_swaps = max_swaps_br,
      n_restarts = n_restarts_br,
      n_candidates = n_candidates_br,
      seed = seed_base + N + n
    )

    assign(key, des, envir = br_cache)
  }

  get(key, envir = br_cache, inherits = FALSE)
}

# -----------------------------------------------------------------------------
# Main simulation loop
# -----------------------------------------------------------------------------

results <- vector("list", length = 0L)
row_counter <- 0L

for (dist in distributions) {
  sigma <- sd_pop(dist)

  for (regime_name in names(regimes)) {
    Ns <- regimes[[regime_name]]

    for (hyp in c("null", "alt")) {
      shift <- if (hyp == "alt") 0.5 * sigma else 0
      is_null <- hyp == "null"

      for (m in seq_len(M)) {
        message(sprintf(
          "Running dist=%s | regime=%s | hyp=%s | replicate=%d/%d",
          dist, regime_name, hyp, m, M
        ))

        N <- sample(Ns, 1)
        n <- floor(N / 2)

        x1 <- rgen(dist, n, shift = 0)
        x2 <- rgen(dist, N - n, shift = shift)
        x <- c(x1, x2)

        des <- get_br_design_cached(
          N = N,
          n = n,
          b_prime = b_prime,
          max_swaps_br = max_swaps_br,
          n_restarts_br = n_restarts_br,
          n_candidates_br = n_candidates_br,
          seed_base = 1000
        )
        br_support <- des$support

        total_alloc <- choose_safe(N, n)
        if (is.finite(total_alloc) && total_alloc <= fp_exact_max) {
          ref_method <- "FP_EXACT"
          ref_S <- NULL
        } else {
          ref_method <- "MC_BENCH"
          ref_S <- S_mc_bench
        }

        seed_base <- 1000000 +
          100000 * match(dist, distributions) +
          10000 * match(regime_name, names(regimes)) +
          1000 * match(hyp, c("null", "alt")) +
          m

        seed_mc <- seed_base + 11
        seed_ba <- seed_base + 22
        seed_ref <- seed_base + 33

        t_br <- system.time(
          p_br <- compute_pvalue(
            data = x,
            n = n,
            design = br_support,
            method = "BR",
            statistic = allocation_t_stat
          )
        )[3]

        t_mc <- system.time(
          p_mc <- compute_pvalue(
            data = x,
            n = n,
            method = "MC",
            S = S_mc,
            statistic = allocation_t_stat,
            seed = seed_mc
          )
        )[3]

        t_ba <- system.time(
          p_ba <- compute_pvalue(
            data = x,
            n = n,
            method = "BA",
            S = S_ba,
            statistic = allocation_t_stat,
            seed = seed_ba
          )
        )[3]

        t_ref <- system.time(
          p_ref <- compute_pvalue(
            data = x,
            n = n,
            method = ref_method,
            S = ref_S,
            statistic = allocation_t_stat,
            seed = seed_ref
          )
        )[3]

        row_counter <- row_counter + 1L
        results[[row_counter]] <- data.frame(
          dist = dist,
          regime = regime_name,
          N = N,
          n = n,
          hyp = hyp,
          is_null = is_null,
          truth_shift = shift,
          replicate = m,
          p_br = p_br,
          p_mc = p_mc,
          p_ba = p_ba,
          p_ref = p_ref,
          ref_method = ref_method,
          time_br = unname(t_br),
          time_mc = unname(t_mc),
          time_ba = unname(t_ba),
          time_ref = unname(t_ref),
          discrepancy = des$discrepancy,
          max_first_dev = des$max_first_dev,
          max_second_dev = des$max_second_dev,
          br_support_size = des$support_size,
          total_allocations = total_alloc,
          seed_mc = seed_mc,
          seed_ba = seed_ba,
          seed_ref = seed_ref,
          stringsAsFactors = FALSE
        )
      }
    }
  }
}

# -----------------------------------------------------------------------------
# Save outputs
# -----------------------------------------------------------------------------

results_df <- do.call(rbind, results)

dir.create("results/main", recursive = TRUE, showWarnings = FALSE)

saveRDS(results_df, file = "results/main/simulation_results.rds")
write.csv(results_df, file = "results/main/simulation_results.csv", row.names = FALSE)

cat("Simulation study completed. Results saved in results/main/.\n")
