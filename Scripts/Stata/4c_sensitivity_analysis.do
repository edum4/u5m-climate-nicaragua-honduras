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
          Sensitivity analyses

Five alternative specifications test the robustness of the
z_precip_12m–mortality association from Model 4:
  S1 — Country×wave stratified Cox: separate baseline hazard per survey
  S2 — Nicaragua only (1998 & 2001)
  S3 — Honduras only (2005 & 2011)
  S4 — Survey-weighted Cox using DHS person weights
  S5 — Time-varying coefficients for z_precip_12m, edyrtotal, and wealthq
       (all three violated the PH test in M4)
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

* Wave variable for S1 stratification
gen byte wave4 = cond(country_str=="NIC" & year==1998, 1, ///
                 cond(country_str=="NIC" & year==2001, 2, ///
                 cond(country_str=="HN"  & year==2005, 3, ///
                 cond(country_str=="HN"  & year==2011, 4, .))))
label define wave4_lbl 1 "NIC-1998" 2 "NIC-2001" 3 "HN-2005" 4 "HN-2011"
label values wave4 wave4_lbl
tab wave4, miss

stset surv_time, failure(event==1) scale(1)

* Covariate list for M3 specification: child + maternal + SES, no country/wave
local covs_m3 "i.kidsex age edyrtotal i.wealthq i.urban"

/*==================================================
          2: S1–S4 sensitivity models
==================================================*/

local covs_s1 "`covs_m3'"
local covs_s2 "`covs_m3' i.year"
local covs_s3 "`covs_m3' i.year"

local cond_s1 ""
local cond_s2 `"if country_str == "NIC""'
local cond_s3 `"if country_str == "HN""'

local opts_s1 "efron strata(wave4)"
local opts_s2 "efron"
local opts_s3 "efron"

forvalues s = 1/3 {
	stcox z_precip_12m `covs_s`s'' `cond_s`s'', `opts_s`s''
	estimates store S`s'
}

* S4: survey-weighted Cox
stset surv_time [pweight=perweight], failure(event==1) scale(1)
stcox z_precip_12m `covs_m3' country_n i.year, vce(robust)
estimates store S4
stset surv_time, failure(event==1) scale(1)   // restore unweighted stset for S5

log using `"`tbl_dir'\cox_phtest_s.txt"', text replace
estat phtest, d
log close

/*==================================================
          3: S5 — Time-varying coefficients for z_precip_12m,
                  edyrtotal, and wealthq (all violated PH test)
==================================================*/

stcox z_precip_12m i.kidsex age edyrtotal i.wealthq i.urban country_n i.year, ///
	efron tvc(z_precip_12m edyrtotal wealthq) texp(log(_t))
estimates store S5

/*==================================================
          4: Summary table — z_precip_12m HR across sensitivity models
==================================================*/

estimates table S1 S2 S3 S4, ///
	eform keep(z_precip_12m) ///
	b(%6.4f) se(%6.4f) p(%6.4f) ///
	stats(N N_fail) // S5 (TVC model) is excluded from this table

* S5 reported separately
estimates table S5, eform keep(z_precip_12m) b(%6.4f) se(%6.4f) p(%6.4f) stats(N N_fail)

/*==================================================
          5: Export results table to Excel
==================================================*/

quietly etable, estimates(S1 S2 S3 S4) ///
	keep(z_precip_12m) ///
	cstat(_r_b,  label("HR") nformat(%6.3f)) ///
	cstat(_r_se, label("SE") nformat(%6.3f)) ///
	mstat(N) mstat(N_fail) ///
	column(estimates) ///
	showstars stars(0.05 "*" 0.01 "**" 0.001 "***") ///
	name(cox_4c)
collect preview
collect export `"`tbl_dir'/cox_sensitivity.xlsx"', name(cox_4c) replace

/*==================================================
          6: Forest plot — S1–S5
==================================================*/

coefplot S5 S4 S3 S2 S1 ///
	, eform keep(z_precip_12m) ///
	xline(1, lpattern(dash) lcolor(red)) ///
	xlabel(0.7(.2)1.8, format(%4.1f)) ///
	xtitle("Hazard ratio (z_precip_12m)") ///
	title("Sensitivity analyses: z_precip_12m HR") ///
	subtitle("Alternative model specifications") ///
	note("Dashed line = null (HR=1). Error bars = 95% CI.") ///
	mlabel format(%6.4f) mlabposition(12) ///
	xsize(6) ysize(4)
graph export `"`fig_dir'\forest_s.png"', replace width(1200)
