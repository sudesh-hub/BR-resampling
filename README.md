# BR-Resampling

Code and data for the paper

> **Balanced Randomization for Efficient Approximation of Exact Permutation Tests**  
> Sudesh Srivastav and Oumar Thiero

This repository contains the simulation code, BR design construction scripts, summary tables,
figures, and manuscript source files used to produce the results in the JNS submission.

---

## Repository structure

```text
BR-Resampling/
  code/
    run_simulations.R
    summarize_results.R
    highdim_illustration.R
    Wilcoxon_BR_Study_Long.R
    make_manuscript_tables.R
    ssession_info.R
  data/
    table_allocation_space_summary.csv
    table_agreement.csv
    table_power.csv
    table_design_diagnostics.csv
    table_reference_method_counts.csv
    table_rmse.csv
    table_runtime.csv
  results/
    wilcoxon/
      raw_data.csv
      summary_by_distribution.csv
      summary_by_regime.csv
      summary_manuscript.csv
      summary_overall.csv
      design_diagnostics.csv
      runtime.txt
  figures/
    figure2_br_mc_fp_density.pdf
  manuscript/
    ms_jns.tex
    ms_jns.pdf
    figure2_br_mc_fp_density.tex
  README.md
```

- **`code/`** – R scripts for simulations, BR construction, and table/figure generation:
  - `run_simulations.R`: main driver for the t-statistic simulation study.
  - `summarize_results.R`: aggregates simulation outputs (RMSE, power, agreement, runtimes).
  - `highdim_illustration.R`: high-dimensional BR illustration (large-scale permutation testing).
  - `Wilcoxon_BR_Study_Long.R`: long-form Wilcoxon rank-sum robustness study (BR vs MC vs BA vs benchmark).
  - `make_manuscript_tables.R`: formats CSV outputs into manuscript-ready tables.
  - `ssession_info.R`: records R session information used to generate results.

- **`data/`** – CSV tables used directly in the manuscript:
  - `table_rmse.csv`, `table_power.csv`, `table_agreement.csv`, `table_runtime.csv`:
    summary metrics by distribution and regime.
  - `table_allocation_space_summary.csv`, `table_design_diagnostics.csv`,
    `table_reference_method_counts.csv`: allocation-space and BR support diagnostics.
  - 

- **`results/wilcoxon/`** – replicate-level Wilcoxon outputs and diagnostics:
  - `raw_data.csv`: replicate-level p-values and indicators for BR, MC, BA, and the benchmark.
  - `summary_by_distribution.csv`, `summary_by_regime.csv`,
    `summary_manuscript.csv`, `summary_overall.csv`: Wilcoxon robustness summaries.
  - `design_diagnostics.csv`: BR support sizes and regime information.
  - `runtime.txt`: runtime and configuration log for the Wilcoxon long run.

- **`figures/`** – final figures:
  - `figure2_br_mc_fp_density.pdf`: Figure 2 in the manuscript, showing FP (benchmark) vs
    BR and MC approximation of the permutation distribution for a representative configuration.

- **`manuscript/`** – JNS submission files:
  - `ms_jns.tex`: LaTeX source for the manuscript.
  - `ms_jns.pdf`: compiled version of the manuscript.
  - `figure2_br_mc_fp_density.tex`: LaTeX wrapper for including Figure 2.

---

## Reproducibility instructions

All results in the manuscript can be regenerated from the R scripts in `code/`. The data are simulated; no external datasets are required.

### 1. Setup

1. Install a recent version of R (≥ 4.0.0) and RStudio (optional but recommended).
2. Install required R packages (from an R console):

   ```r
   install.packages(c(
     "data.table",
     "ggplot2",
     "future",
     "future.apply",
     "progressr"
   ))
   ```

3. Set the working directory to the root of this repository, for example:

   ```r
   setwd("path/to/BR-Resampling")
   ```

### 2. Main t-statistic simulation study

1. **Run the simulations**:

   ```r
   source("code/run_simulations.R")
   ```

   This script generates the main t-statistic simulation results and writes
   intermediate outputs and summaries into `data/` (and, if configured, into
   a `results/` subfolder).

2. **Summarize results**:

   ```r
   source("code/summarize_results.R")
   ```

   This script computes RMSE, power, agreement, and runtime summaries used in
   the manuscript tables (e.g., `table_rmse.csv`, `table_power.csv`,
   `table_agreement.csv`, `table_runtime.csv`).

3. **Create manuscript tables**:

   ```r
   source("code/make_manuscript_tables.R")
   ```

   This script formats the CSV summaries into the exact table structures that
   appear in the paper and saves them in `data/` (and/or prints them to the console).

### 3. Wilcoxon rank-sum robustness study

Run the long-form Wilcoxon study:

```r
source("code/Wilcoxon_BR_Study_Long.R")
```

This script:

- Generates replicate-level Wilcoxon p-values for BR, MC, BA, and the benchmark.
- Writes raw results and summaries into `results/wilcoxon/` and `data/`:
  - `results/wilcoxon/raw_data.csv`
  - `results/wilcoxon/summary_*.csv`, `design_diagnostics.csv`, `runtime.txt`
  - `data/summary_by_distribution.csv`, `data/summary_by_regime.csv`,
    `data/summary_manuscript.csv`, `data/summary_overall.csv`

These outputs correspond to the Wilcoxon tables and descriptions in the manuscript.

### 4. High-dimensional illustration

To reproduce the high-dimensional BR illustration:

```r
source("code/highdim_illustration.R")
```

This script runs the large-scale permutation experiments and produces the
high-dimensional tables and figures used in the application section of the paper.

### 5. Figure 2: BR vs MC vs permutation benchmark

Figure 2 is generated by running:

```r
source("code/figure2_br_mc_fp_density.R")  # or equivalent code segment
```

(if you keep the figure code in its own script; otherwise, use the code in
`figure2_br_mc_fp_density.tex`/manuscript Section 4.3 to regenerate it).

The output is:

- `figures/figure2_br_mc_fp_density.pdf`

which is included in the manuscript via `manuscript/figure2_br_mc_fp_density.tex`.

### 6. Session info

For reproducibility, the script:

```r
source("code/ssession_info.R")
```

can be used to write the R version and package versions used to generate the results.

---

## License and acknowledgments

All data in this repository are simulated. The code is provided for research and
reproducibility purposes. Please cite the associated manuscript if you use the
BR methods or scripts in your own work.

The authors thank Tulane University and USTTB for institutional support.
