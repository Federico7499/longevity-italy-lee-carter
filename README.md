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

## Phase 2: choosing a stochastic mortality model (Italy, ages 55–100)

Five models fitted with `StMoMo` on Human Mortality Database data for Italy (deaths and exposures, 1974–2023), by sex, and compared out of sample on two backtest windows (fit to 2004 → forecast 2005–2014; fit to 2009 → forecast 2010–2019).

| Model | Average rank in backtest (F / M) | Average error on log q (F / M) | Error on e65 at the end of the window |
|---|---|---|---|
| **Lee-Carter** | **1 / 1** | **0.067 / 0.053** | +0.5 / −0.2 years |
| APC | 2 / 4.5 | 0.108 / 0.088 | +0.7 / +0.6 years |
| CBD | 4 / 4.5 | 0.217 / 0.091 | +0.4 / −0.2 years |
| Renshaw-Haberman | 4 / 6 | 0.206 / 0.587 | +2.4 / from +4.8 to −9.0 years |
| M7 | 6 / 2 | 0.238 / 0.057 | −1.3 / −0.5 years |

- **The best model in sample is not the best forecaster.** Renshaw-Haberman has the lowest BIC for both sexes but the least stable forecasts: for males its error on life expectancy at 65 goes from +4.8 years in one window to −9.0 in the other. Projecting the cohort effect without drift does not fix it.
- **Lee-Carter ranks first in both windows and for both sexes** (forecast error on death probabilities of 4–6%), even though its residuals show an unexplained cohort effect. It is adopted as the reference model, with **APC as a sensitivity**: APC projects higher longevity, which is the prudent direction for annuity pricing.
- **M7 is excluded**: its female projection is implausible (life expectancy at 65 falling to about 6 years by 2050) and its male projection is almost flat, despite good 10-year backtests.
- **COVID years matter**: fitting on 1974–2023 instead of 1974–2019 lowers projected life expectancy at 65 in 2050 by 0.4–0.7 years (Lee-Carter, APC, CBD).

![Life expectancy at 65 by model](output/fase2/figures/e65_proiezioni.png)

## Phase 2b: Lee-Carter and APC by region (ages 55–100, 2002–2024)

The two models selected in phase 2 are fitted separately for each region and sex, using ISTAT data: population on 1 January by single age (reconstructed series 2002–2018, current series 2019–2025) and mortality rates from the regional life tables.

| Life expectancy at 65 (Lee-Carter) | Lombardy | Lazio | Sardinia |
|---|---|---|---|
| Females, 2024 → 2050 | 23.0 → 25.3 | 22.5 → 24.9 | 22.9 → 25.1 |
| Males, 2024 → 2050 | 20.2 → **23.9** (+3.7) | 19.8 → 22.9 (+3.1) | 19.6 → 22.2 (+2.6) |

- **Out-of-sample errors are 4–8%** (fit 2002–2014, forecast 2015–2019), with no clear winner between Lee-Carter and APC at regional level.
- **Lombardy leads and is projected to pull ahead**, especially for males: again, the most polluted region shows the fastest improvement.
- **Caveat: the projected speed depends on the end points.** With a random walk with drift, the slope of *k_t* is set by its first and last values. The sharp fall for Lombardy males in 2023–2024 is likely a post-COVID rebound (frail individuals died earlier), which may inflate the projection. This also explains why excluding the COVID years barely changes Lombardy's male projection (−0.2 years, against about −1 year elsewhere).
- **APC is not always the prudent choice**: for Italy it projects higher longevity than Lee-Carter, but for Sardinia it projects lower. A prudent longevity assumption should take the higher of the two models region by region.
- **Data note:** deaths are reconstructed as ISTAT life-table rates × exposure, because observed deaths by single age and region are not readily available. Life-table rates are slightly smoothed at the oldest ages, so backtest errors may be somewhat understated.

![Life expectancy at 65 by region](output/fase2b/figures/e65_regioni.png)

## Phase 4: pricing a life annuity by region

Single premium for a life annuity of **10,000 euro a year**, paid in advance from age 65, issued in 2025 (cohort born in 1960). Technical rate 2%, no expense loadings. The premium is computed with the **period table 2024** (today's mortality, no future improvement) and with the **projected cohort table** from phase 2b.

| Single premium, Lee-Carter, 2% | Lombardy | Lazio | Sardinia |
|---|---|---|---|
| Females: period → cohort | 185,800 → 195,600 € | 182,300 → 192,500 € | 185,000 → 194,800 € |
| Males: period → cohort | 167,000 → **180,800 €** | 163,800 → 175,100 € | 162,300 → 171,800 € |
| Cost of ignoring improvements (F / M) | +5.2% / **+8.3%** | +5.6% / +6.9% | +5.3% / +5.8% |

- **Ignoring future mortality improvements underprices the annuity by 5–8%**: up to about 13,800 euro per contract for a Lombard man.
- **Region matters for men, much less for women.** A 65-year-old Lombard man costs 5.3% more than a Sardinian one (+9,000 euro); for women the three regions are within about 1.5%.
- **Model risk is as large as regional risk.** For Lombard men, APC gives a premium of 188,300 euro, 4.2% above Lee-Carter. A prudent basis taking the higher of the two models per region selects APC for Lombardy (both sexes) and Lazio men.
- **COVID years lower the premium by 0.7–3.6%** compared with fitting on 2002–2019. Lombard men are again the exception (−0.7% with Lee-Carter), because of the post-COVID rebound discussed in phase 2b.
- **Limitations:** general population mortality, not annuitant mortality (annuitants live longer, so real premiums would be higher); flat technical rate; table closed at age 100; no expense loadings; no uncertainty on the projection (next phase).

![Single premium by region](output/fase4/figures/premi_rendita.png)

## What the project does

1. Builds age-specific mortality rates and complete life tables from ISTAT data, for total population, females and males.
2. Computes period life expectancy at birth and at age 65, and measures the COVID-19 shock (2019 → 2020) and the recovery.
3. Fits a Lee-Carter model, checks the residuals (age × year heatmaps) and projects the mortality index *k_t* 30 years ahead with a random walk with drift.

## Repository structure

```
data/raw/                 ISTAT regional life tables, one CSV per year
data/hmd/                 HMD Italy deaths and exposures (not included, see Data)
Progetto_demografia.R     phase 1: ISTAT life tables → Lee-Carter by region
Fase2_modelli_stocastici.R  phase 2: model comparison and backtesting (StMoMo)
Fase2b_regioni.R          phase 2b: Lee-Carter and APC by region
Fase4_rendite.R           phase 4: life annuity pricing by region
data/istat_pop/           ISTAT population on 1 January by single age, 2002–2025
output/                   phase 1 tables and figures
output/fase2/             phase 2 tables and figures
output/fase2b/            phase 2b tables and figures
output/fase4/             phase 4 tables and figures
```

## How to run

Open the project in RStudio from the `.Rproj` file and run the script of each phase: `Progetto_demografia.R`, `Fase2_modelli_stocastici.R`, `Fase2b_regioni.R`, `Fase4_rendite.R` (run phase 2b first: it saves the fitted models used by phase 4). Missing packages (`demography`, `StMoMo`, `tidyverse`, `ggprism`) are installed automatically.

## Data

- **Phase 1:** ISTAT regional life tables, 1974–2024 (`datiregionalicompleti<year>-2.csv`). Ages 100 and over are grouped into an open-ended 100+ class.
- **Phase 2:** Human Mortality Database, Italy, period `Deaths_1x1.txt` and `Exposures_1x1.txt`. HMD data are not redistributed here: download them from [mortality.org](https://www.mortality.org) (free registration) and save them in `data/hmd/`.

- **Phase 2b:** ISTAT population on 1 January by sex and single age: reconstructed population 2002–2019 (`pop_ricostruita_<region>.csv`) and current series 2019–2025 (`pop_<region>_<year>.csv`), from demo.istat.it. The two series match exactly in 2019.

  *HMD. Human Mortality Database. Max Planck Institute for Demographic Research (Germany), University of California, Berkeley (USA), and French Institute for Demographic Studies (France). Available at www.mortality.org (data downloaded in October 2026).*

## Roadmap

- [x] Compare stochastic mortality models (LC, RH, APC, CBD, M7) with out-of-sample backtesting (`StMoMo`), Italy
- [x] Apply the selected models to the three regions (ISTAT population and life tables, 2002–2024)
- [ ] Include air quality (PM10, PM2.5, NO₂) as a covariate in a Poisson GLM, at the provincial level
- [x] Price life annuities with projected cohort tables, by region and sex
- [ ] Simulate longevity risk and compare the resulting capital with the Solvency II standard-formula shock

*Author: Federico Cerri, MSc in Statistical, Financial and Actuarial Sciences, University of Bologna*
