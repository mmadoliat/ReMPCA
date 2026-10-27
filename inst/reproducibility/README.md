# Reproducibility materials

This folder contains the R scripts that reproduce the results in

> Pourmoshir, M. and Maadooliat, M. *Regularized Principal Component Analysis
> for Multivariate Two-way Hybrid Data.* Journal of Data Science (in revision).

After installing the package, the folder is available at
`system.file("reproducibility", package = "ReMPCA")`.

| Script | Data | Reproduces |
|---|---|---|
| `simulation_study.R` | simulated (seeded) | Section 3: Table 2, Figure 2 (`fig2_reconstruction_error_boxplot.pdf`), Figure 3 (`fig3_representative_run.pdf`); Section 3.1: Table 3 (benchmark against pre-smoothing hybrid PCA, HFV-PCA, PMD and the two-way regularized SVD); Supplementary Section S2: Table S1 and Figure S4 (`figS_convergence_traces.pdf`); Supplementary Section S3: Tables S2 and S3 |
| `ReMPCA_bike_example.R` | `bike_day`, `bike_hour` (UCI Bike Sharing, shipped with the package) | Section 4: Figure 4 (`bike_compare_u_and_v_HPCA_ReHPCA.pdf`), Figure 5 (`heatmap_HPCA_vs_ReHPCA_fitted.pdf`), Figure 6 (`bike_mean_ReHPCA.pdf`) |
| `Appliances_Energy_Example.R` | `energydata` (UCI Appliances Energy Prediction, shipped with the package) | Supplementary Section S1: Figures S1-S3 |
| `figure1/fig1_plot.R` | `figure1/fig1_*_from_pdf.csv.gz` | Figure 1 (chessboard illustration; see below) |

## Running

```r
# install.packages(c("tidyverse", "lubridate", "patchwork", "scales",   # data analyses
#                    "MFPCA", "funData", "PMA", "fda"))                 # simulation benchmark
Sys.setenv(REMPCA_OUTPUT_DIR = "ReMPCA_output")   # where figures and tables are written
rep_dir <- system.file("reproducibility", package = "ReMPCA")
source(file.path(rep_dir, "ReMPCA_bike_example.R"))          # about 15 s
source(file.path(rep_dir, "Appliances_Energy_Example.R"))    # about 30 s
source(file.path(rep_dir, "figure1", "fig1_plot.R"))         # a few seconds

Sys.setenv(REMPCA_CORES = "8")                                # parallel workers (forked; 1 on Windows)
source(file.path(rep_dir, "simulation_study.R"))
```

Each script writes its figures (PDF, at the final print width of 6.3 in), tables
and `sessionInfo()` to `REMPCA_OUTPUT_DIR`.

### Simulation study

`simulation_study.R` fits 900 simulated data sets (9 noise-by-sparsity scenarios
x 100 runs) with six ReMPCA configurations and four competing methods, and then
runs the convergence diagnostics (20 data sets per scenario, 20 random starts
each). Each data set is generated from its own seed, so the results do not
depend on the number of workers. With 64 workers the full study takes about
15 minutes; the run time scales inversely with `REMPCA_CORES`. Further options
(environment variables): `REMPCA_NREP` (runs per scenario, default 100),
`REMPCA_SCENARIOS` (subset of scenarios, e.g. `"5"`), and `REMPCA_STAGE`
(`"all"`, `"fit"`, `"conv"`, `"traces"` or `"figs"`). Outputs include
`simulation_results.csv` (one row per data set and method), `table2_rows.tex`,
`table3_benchmark_rows.tex`, `tableS1_convergence_summary.csv`,
`tableS_benchmark_sinu_rows.tex`, `tableS_benchmark_sinv_rows.tex` and
`tableS_benchmark_time.csv`.

All tuning settings used in the paper (smoothing grids for GCV, sparsity grids
and the minimum-error rule for 5-fold CV, tuning smoothness before sparsity,
soft thresholding) are fixed in the script.

The folder `reference_output/` contains the outputs of the run reported in the
paper (R 4.6.1, ReMPCA 0.1.0, fda 6.3.0, MFPCA 1.3-11, funData 1.3-9, PMA 1.2-4;
see `sessionInfo_simulation.txt`): the per-data-set results
(`simulation_results.csv.gz`), the convergence diagnostics
(`convergence_diagnostics.csv.gz`), the summary table and the LaTeX rows of
Tables 2, 3, S2 and S3. A rerun can be checked against these files, for example

```r
ref <- read.csv(file.path(rep_dir, "reference_output", "simulation_results.csv.gz"))
new <- read.csv(file.path(Sys.getenv("REMPCA_OUTPUT_DIR"), "simulation_results.csv"))
key <- c("scenario", "rep", "method")
m <- merge(ref, new, by = key, suffixes = c(".ref", ".new"))
max(abs(m$rel_err.ref - m$rel_err.new))
```

While tuning, R reports many warnings of the form "Algorithm did not converge
within the maximum number of iterations". They come from the fits on the
cross-validation training folds at candidate tuning values, where thresholded
iterations can oscillate; the last iterate is then used to compute the
cross-validation error. They do not indicate a
failure of the script. Convergence of the final fits at the selected tuning
values is assessed separately in the convergence diagnostics
(`convergence_diagnostics.csv`, Table S1).

The representative run shown in Figure 3 (run 73 of the medium-noise,
medium-sparsity scenario) is chosen by the script as the run with the median
error of the two-way smoothness and sparsity configuration.

### Bike-sharing application

The non-functional block contains the covariates listed in the paper's variable
table (daily temperature, apparent temperature, humidity and wind speed, season
indicators, and weekday/weekend indicators). All tuning parameters are fixed in
the script (one candidate value per parameter), and the cross-validation folds
use a fixed internal seed, so the results are deterministic.

### Figure 1

Figure 1 is a conceptual illustration of two-way regularization on a noisy
chessboard. The folder `figure1/` contains the plotted values (the six
reconstructed matrices and the estimated singular vectors), recovered from the
vector graphics of the figure in the earlier version of the manuscript by
`fig1_extract_from_pdf.py`, and the script `fig1_plot.R` that draws the figure.
