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
					Cox models

Four nested models assess whether the crude association between
z_precip_12m and under-5 mortality persists after covariate adjustment:
  Model 1 —>  Child characteristics:       kidsex
  Model 2 —>  + Maternal characteristics:  age, edyrtotal
  Model 3 —>  + Socioeconomic status:      wealthq, urban
  Model 4 —>  Full (+ context):            country_n, year
==================================================*/

clear all
set more off
set scheme s2color

/*==================================================
          0: Path configuration
==================================================*/

local root    `""'  // ← INSERT YOUR ROOT PATH HERE (all other paths build from this)
local in_csv  `"`root'\Child mortality and climate\Data\processed\dhs_final.csv"'
local tbl_dir `"`root'\Child mortality and climate\Output\tables\Stata"'
local fig_dir `"`root'\Child mortality and climate\Output\figures\Stata"'

/*==================================================
          1: Load and prepare data
==================================================*/

import delimited using `"`in_csv'"', varnames(1) encoding(UTF-8) clear

gen byte country_n = (country_str == "HN")
label define country_lbl 0 "Nicaragua" 1 "Honduras", replace
label values country_n country_lbl

foreach v in kidsex urban wealthq {
	encode `v', gen(`v'_n)
	drop `v'
	rename `v'_n `v'
}

recode wealthq (1=3) (3=1)
label define wealthq 1 "Poorest" 2 "Poorer" 3 "Middle" 4 "Richer" 5 "Richest", replace
label values wealthq wealthq

stset surv_time, failure(event==1) scale(1)

/*==================================================
          2: Cox models

Efron method for tied event times. 
M0 = unadjusted (baseline),
M1–M4 = progressive covariate adjustment.
==================================================*/

local covs0 ""
local covs1 "i.kidsex"
local covs2 "i.kidsex age edyrtotal"
local covs3 "i.kidsex age edyrtotal i.wealthq i.urban"
local covs4 "i.kidsex age edyrtotal i.wealthq i.urban country_n i.year"

forvalues m = 0/4 {
	stcox z_precip_12m `covs`m'', efron
	estimates store M`m'
	estat concordance
}

/*==================================================
          3: Summary table — z_precip_12m HR across models
==================================================*/

estimates table M0 M1 M2 M3 M4, ///
	eform keep(z_precip_12m) ///
	b(%6.4f) se(%6.4f) p(%6.4f) ///
	stats(N N_fail)

/*==================================================
          4: Export results table to Excel
==================================================*/

quietly etable, estimates(M0 M1 M2 M3 M4) ///
	keep(z_precip_12m) ///
	cstat(_r_b,  label("HR") nformat(%6.3f)) ///
	cstat(_r_se, label("SE") nformat(%6.3f)) ///
	mstat(N) mstat(N_fail) ///
	column(estimates) ///
	showstars stars(0.05 "*" 0.01 "**" 0.001 "***") ///
	name(cox_models)
collect preview
collect export `"`tbl_dir'/cox_models.xlsx"', name(cox_models) replace

/*==================================================
          5: PH assumption test — full model (M4)

Schoenfeld residuals test on the fully-adjusted model.
A significant p-value (p<0.05) suggests non-proportional hazards.
==================================================*/

stcox z_precip_12m i.kidsex age edyrtotal i.wealthq i.urban country_n i.year, efron

log using `"`tbl_dir'\cox_phtest_ad.txt"', text replace
estat phtest, detail
log close

/*==================================================
          6: Forest plot — z_precip_12m HR across models
==================================================*/

coefplot M4 M3 M2 M1 M0 ///
	, eform keep(z_precip_12m) ///
	xline(1, lpattern(dash) lcolor(red)) ///
	xlabel(0.9(.1)1.5, format(%4.1f)) ///
	xtitle("Hazard ratio (z_precip_12m)") ///
	title("HR for 12-month precipitation z-score") ///
	subtitle("Sequential covariate adjustment — under-5 mortality") ///
	note("Dashed line = null (HR=1). Error bars = 95% CI.") ///
	mlabel format(%6.4f) mlabposition(12) ///
	xsize(6) ysize(4)
graph export `"`fig_dir'\forestplot.png"', replace width(1200)
