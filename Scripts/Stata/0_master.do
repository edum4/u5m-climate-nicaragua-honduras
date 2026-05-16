/*==================================================
project:       Survival Analysis – DHS + WorldClim Climate Data
               Cox PH Model: Under-5 Child Mortality
               Nicaragua (1998, 2001) & Honduras (2005, 2011)
Author:        Eduardo Pacheco
----------------------------------------------------
Creation Date:    16 March 2026
Modification Date: 29 April 2026
Do-file version:    02
==================================================*/

/*==================================================
          Master do-file

Runs the full pipeline in order:
  1a — Rolling precipitation baselines
  1b — Per-child climate exposure
  2a — DHS data inspection
  2b — Merge DHS + climate exposure
  3  — Kaplan-Meier exploration
  4a — Unadjusted Cox model
  4b — Adjusted Cox models
  4c — Sensitivity analyses
==================================================*/

clear all
set more off

/*==================================================
          0: Path configuration
==================================================*/

local root    `""'  // ← INSERT YOUR ROOT PATH HERE (all other paths build from this)
local scripts `"`root'\Child mortality and climate\Scripts\Stata"'

/*==================================================
          1: Run pipeline
==================================================*/

do `"`scripts'\1a_precipitation baselines.do"'
do `"`scripts'\1b_child_climate_exposure.do"'
do `"`scripts'\2a_dhs_inspection.do"'
do `"`scripts'\2b_merge_dhs_climate.do"'
do `"`scripts'\3_exploration_kmc.do"'
do `"`scripts'\4a_coxmodels_unadjusted.do"'
do `"`scripts'\4b_coxmodels_adjusted.do"'
do `"`scripts'\4c_sensitivity_analysis.do"'
