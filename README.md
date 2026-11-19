# R code: VIC first-difference model, GAM peak detection and depth profiles

Two versions of the same analysis:

- `code/` – seven scripts, each runs on its own, in any order.
- `combined/ALL_analysis_combined.R` – the whole analysis in one script.

Figures are shown in the RStudio Plots panel and tables are printed in the console.
Nothing is saved to disk.

## Folder layout

```
data folder/
├── Treatments.csv
├── 2023/TAREC/*.csv
├── 2024/TAREC/*.csv
├── 2025/ESAREC/*.csv
├── 2025/TAREC/*.csv
├── 2026/ESAREC/*.csv
└── Journal_submission/
    ├── code/
    └── combined/
```

The scripts find the data folder by looking for `Treatments.csv` in the current folder,
the script's folder, and one or two levels above. Otherwise set
`main_folder <- "path/to/data folder"` at the top of the script.

## How to run

Open a script in RStudio and click **Source**. Missing packages are installed from CRAN.

## Scripts

| Script | Shows |
|---|---|
| 01_FD_coefficients.R | Table S1 – FD coefficients for every sensor and depth |
| 02_FD_coefficient_boxplots.R | Fig. 1, Fig. 2 – FD coefficients by depth |
| 03_FD_significance_table.R | Table 1 – % of FD models with a significant coefficient |
| 04_GAM_example_sensor.R | Fig. 3, Fig. 4 – GAM of residuals, sensor 01209A5A (TAREC 2024, 10 cm) |
| 05_GAM_peak_table.R | Table 2 – highest GAM peak and Dixon Q test, every sensor and depth |
| 06_Peak_timing_boxplots.R | Fig. 5 – peak timing (days after planting); white, blue, grey |
| 07_Sig_peak_depth_profile.R | Fig. 6 – % of sensors with a significant peak by depth |
| combined/ALL_analysis_combined.R | all of the above (about 5 min) |

## Methods in short

- FD model: ΔVIC = b1·ΔVWC + b2·ΔTemp (no intercept), `glm()`, each sensor and depth.
- Residual = VIC − (b0 + b1·VWC + b2·Temp), b0 = mean(VIC) − b1·mean(VWC) − b2·mean(Temp).
- GAM: `mgcv::gam(resid ~ s(index, k = 20), method = "REML")`.
- Peak: highest turning point of the GAM, first and last points included.
- Dixon Q test (`outliers::dixon.test()`) on the GAM turning points; significant when p < 0.05 and Q ≥ 0.5.
- Planting dates, sensor list and N groups from `Treatments.csv` (55 sensors). ESAREC 2024 not used.

## Software

R 4.3.3; lubridate 1.9.3, ggplot2 3.4.4, mgcv 1.9.1, outliers 0.15, patchwork 1.2.0 (ggplot2 ≥ 3.4.0 needed).
