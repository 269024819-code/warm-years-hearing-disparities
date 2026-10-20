# Unusually warm years and sex disparities in hearing-loss disability

Analysis code for the manuscript **"Unusually warm years are associated with greater sex disparities in hearing-loss disability among older adults"** (160-country panel, 2003–2019).

The code constructs the warm-year exposure, aggregates population, estimates age-specific Hot×Female interactions, calculates burden, performs the counterfactual (Shapley) decomposition, propagates uncertainty by Monte Carlo simulation, and reproduces all main and Extended Data figures.

## Requirements

- **R** (≥ 4.3): `dplyr`, `tidyr`, `fixest`, `ggplot2`, `patchwork`, `MASS`, `sandwich`, `plm`
- **Python** (≥ 3.10): `pandas`, `geopandas`, `matplotlib`

## Data (`data/`)

| File | Content | Source |
|---|---|---|
| `Standard_Gender_Hearing_Panel_2003_2019.csv.gz` | Country–year–age–sex–severity hearing-loss YLD panel | GBD 2023 (IHME) |
| `hearing_temperature_analysis_panel_2003_2019.csv.gz` | Country–year temperature + hearing panel | CRU TS 4.08 / ERA5 + GBD 2023 |
| `country_age_sex_female_population_early_late_2003_2007_vs_2015_2019.csv` | Female population by country, age group and period | WorldPop |
| `country_zh_en_map.csv` | Country-name harmonization table | — |
| `table1_modsevere_160.csv` | Age-specific Hot×Female estimates with 95% CIs (Table 1) | this study |
| `ne_110m_admin_0_countries.*` | Country boundaries for maps | Natural Earth (public domain) |

Decompress the two `.csv.gz` files into `data/` before running (`gunzip *.gz`, or R reads `.gz` directly). GBD-derived estimates are subject to the IHME terms of use; raw data are available from the original providers.

## Running

Run scripts in numerical order from the repository root:

```bash
Rscript code/01_build_panel.R          # reconstruct analysis panels, validate Table 1
Rscript code/02_se_robustness.R        # clustered SEs, region×year FE, leave-one-region-out
Rscript code/03_se_extended.R          # two-way clustering, wild bootstrap, pooled/gradient tests
Rscript code/04_fdr_and_climate_trends.R  # 36-test FDR; country temperature trends
Rscript code/05_mc_rho_sensitivity.R   # Monte Carlo, cross-age correlation rho ladder
Rscript code/06_mc_gbd_uncertainty.R   # Monte Carlo, GBD outcome uncertainty (two-component)
Rscript code/07_prep_fig3_data.R       # country-level values for maps
Rscript code/08_draw_fig1.R            # Figure 1
python3 code/09_draw_fig3.py           # Figure 3 (world maps)
Rscript code/10_draw_fig4.R            # Figure 4 (decomposition)
python3 code/11_draw_edfig5.py         # Extended Data Fig. 5 (warming-trend map)
```

Outputs are written to `out/` (tables) and the repository root (figures).

## License

Code: MIT (see `LICENSE`). Natural Earth boundaries: public domain. GBD-derived estimates: subject to IHME terms of use.
