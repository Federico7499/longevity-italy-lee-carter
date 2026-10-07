# Longevity in Italy: Lombardy, Lazio and Sardinia (1974–2024)

Mortality analysis and Lee-Carter projection for three Italian regions chosen for their different levels of air pollution (PM10, NO₂):

- **Lombardy**: high pollution and population density
- **Lazio**: intermediate case with a large urban area
- **Sardinia**: low pollution benchmark

![Life expectancy at birth](output/figures/e0.png)

## Key results

| | Lombardy | Lazio | Sardinia |
|---|---|---|---|
| Life expectancy at birth, 1974 | 71.7 | 73.5 | 73.8 |
| Life expectancy at birth, 2024 | 84.1 | 83.4 | 83.0 |
| COVID-19 shock, 2019 → 2020 | **−2.2 years** | −0.6 years | −0.6 years |
| 2024 vs 2019 | +0.5 years | +0.2 years | +0.0 years |
| Variance explained by Lee-Carter | 92.3% | 91.0% | 87.1% |
| Projected life expectancy, 2054 | 89.1 | 87.8 | 87.1 |

*Total population, period life expectancy.*

- Lombardy started with the **lowest** life expectancy of the three regions in 1974 and had the **highest** by 2024. The more polluted region shows the fastest long-term improvement, which suggests that air quality cannot be isolated without controlling for other regional factors (income, health care).
- The 2020 shock was almost **four times larger in Lombardy** than in Lazio and Sardinia (−2.5 years for males). By 2024 all three regions were back above their pre-pandemic level.
- Lee-Carter fits Lombardy and Lazio well (around 91–92% of variance explained). The fit is weaker for Sardinia, especially for females (77%), because a smaller population produces noisier rates.

![k_t forecast](output/figures/lc_kt_forecast.png)

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
