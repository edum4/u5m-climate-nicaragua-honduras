# Under-5 mortality and climate variability in Nicaragua and Honduras

## Overview

This project examines whether precipitation variability during the first year of life is
associated with under-5 mortality in Nicaragua (1998, 2001) and Honduras
(2005, 2011). Using four waves of Demographic and Health Survey (DHS) data linked
to high-resolution WorldClim precipitation records, I construct an individual-level
climate exposure measure and evaluate its association with mortality through a series
of nested Cox proportional-hazards models.

The results provide no statistically significant evidence that precipitation anomalies
during the first year of life are associated with child mortality in any of the specifications.

This project aims to illustrate the workflow for a survival analysis linking individual-level
data and climate exposure, from constructing the climate exposure model to estimating
and robustness analysis. This exercise does not represent formal research oriented
toward a particular result.

## Data

| Source | Coverage | Unit of observation |
|---|---|---|
| IPUMS-DHS | Nicaragua 1998, 2001; Honduras 2005, 2011 | Child birth record |
| WorldClim CRU TS4.09 (5 arc-min) | 1980–2012 | Monthly raster, ~10 km resolution |
| Administrative boundaries | Honduras (18 dep.), Nicaragua (17 dep.) | GeoJSON polygons |

## Pipeline

```
    Regional monthly precipitation (from WorldClim rasters)
                       ↓
    [Stata] → Rolling 15-year regional baselines
                       ↓
    [Stata] → Child-level 12-month z-score exposure
                       ↓
    [Stata/R] → Cox models + sensitivity analyses
```

1. **Climate baselines (Stata):** Rolling 15-year region × calendar-month means and
   standard deviations, anchored to the year before each birth. Regional monthly
   precipitation was previously extracted from WorldClim rasters using area-weighted
   methods.
2. **Child exposure (Stata):** `z_precip_12m` — mean of 12 monthly z-scores over the
   child's first year of life.
3. **Survival analysis (Stata + R):** Five nested Cox models (M0–M4) plus five
   sensitivity specifications (S1–S5).

## Skills demonstrated

- **Spatial-temporal data linkage:** DHS survey records linked to
  regional climate data by administrative region and birth date.
- **Survival analysis in Stata and R** `stcox` / `coxph`, Efron ties, survey weights,
  time-varying coefficients (`tvc()` / `tt()`), Schoenfeld residual diagnostics.
- **Reproducible pipeline:** modular scripts with a master do-file/R script.
- **DHS microdata:** experience with multi-wave survey data, including data cleaning and construction of demographic indicators.

## Replication

Each script defines its own root path at the top of the file (`local root` in Stata, `root <-` in R). To run the full pipeline or any individual script, update that line to match your local directory. No other changes are required.

## Repository structure

```
├── Scripts/
│   ├── Python/          
│   ├── Stata/           # Data construction, exploration, and Cox models (0_master.do)
│   └── R/               # Adjusted models and sensitivity analyses (0_master.R)
├── Data/
│   ├── Raw climate data/
│   ├── GIS maps/
│   ├── IPUMS-DHS/
│   └── processed/
└── Output/
    ├── figures/{Stata,R}/
    └── tables/{Stata,R}/
```
