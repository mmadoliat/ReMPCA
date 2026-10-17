# Reproducibility materials

This folder contains the R scripts that reproduce the real-data analyses in

> Pourmoshir, M. and Maadooliat, M. *Regularized Principal Component Analysis
> for Multivariate Two-way Hybrid Data.* Journal of Data Science (in revision).

After installing the package, the folder is available at
`system.file("reproducibility", package = "ReMPCA")`.

| Script | Data (shipped with the package) | Reproduces |
|---|---|---|
| `ReMPCA_bike_example.R` | `bike_day`, `bike_hour` (UCI Bike Sharing) | Bike-sharing application: HPCA vs ReHPCA components (`bike_compare_u_and_v_HPCA_ReHPCA.pdf`), daily-hour reconstructed heatmaps (`heatmap_HPCA_vs_ReHPCA_fitted.pdf`), mean +/- ReHPCA PCs (`bike_mean_ReHPCA.pdf`) |
| `Appliances_Energy_Example.R` | `energydata` (UCI Appliances Energy Prediction) | Supplementary Figures S1-S3 (`appliances_compare_u_and_v_HPCA_ReHPCA_PC_specific.pdf`, `appliances_ReHPCA_mean_plus_minus_PC1_PC2_symbols.pdf`, `appliances_u_scores_HPCA_vs_ReHPCA.pdf`) |

## Running

```r
# install.packages(c("tidyverse", "lubridate", "patchwork", "scales"))
Sys.setenv(REMPCA_OUTPUT_DIR = "ReMPCA_output")   # where figures are written
source(system.file("reproducibility", "ReMPCA_bike_example.R", package = "ReMPCA"))
source(system.file("reproducibility", "Appliances_Energy_Example.R", package = "ReMPCA"))
```

Each script writes the figures (PDF), the fitted `ReMPCA` objects (`*_fits.rds`)
and the `sessionInfo()` to `REMPCA_OUTPUT_DIR`. The bike-sharing script runs in
about 15 seconds and the appliances script in about 30 seconds on a laptop.

All tuning parameters used in the paper are fixed in the scripts (one
candidate value per parameter), and the cross-validation folds use a fixed
internal seed, so the results are deterministic.
