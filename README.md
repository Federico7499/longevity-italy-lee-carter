# Longevity in Italy: Lombardy, Lazio and Sardinia (1974–2024)

Mortality analysis and Lee-Carter projection for three Italian regions chosen for their different levels of air pollution (PM10, NO₂):

- **Lombardy**: high pollution and population density
- **Lazio**: intermediate case with a large urban area
- **Sardinia**: low pollution benchmark

## What the project does

1. Builds age-specific mortality rates and complete life tables from ISTAT data, for total population, females and males.
2. Computes period life expectancy at birth and at age 65, and measures the COVID-19 shock (2019 → 2020) and the recovery.
3. Fits a Lee-Carter model, checks the residuals (age × year heatmaps) and projects the mortality index *k_t* 30 years ahead with a random walk with drift.

## Repository structure

```
data/raw/                 ISTAT regional life tables, one CSV per year
Progetto_demografia.R     full pipeline: data → models → tables → figures
output/                   result tables (CSV)
output/figures/           figures (PNG)
```

## How to run

Open the project folder in RStudio (or set it as the working directory) and run `Progetto_demografia.R`. Missing packages (`demography`, `tidyverse`, `ggprism`) are installed automatically.

## Data

ISTAT regional life tables, 1974–2024 (`datiregionalicompleti<year>-2.csv`). Ages 100 and over are grouped into an open-ended 100+ class.

## Roadmap

- [ ] Compare stochastic mortality models (Lee-Carter, Renshaw-Haberman, CBD, APC) with out-of-sample backtesting (`StMoMo`), using observed deaths and exposures
- [ ] Include air quality (PM10, PM2.5, NO₂) as a covariate in a Poisson GLM, at the provincial level
- [ ] Price life annuities with projected cohort tables, by region and sex
- [ ] Simulate longevity risk and compare the resulting capital with the Solvency II standard-formula shock

*Author: Federico Cerri, MSc in Statistical, Financial and Actuarial Sciences, University of Bologna*
