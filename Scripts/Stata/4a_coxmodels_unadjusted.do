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
          Unadjusted Cox model

Fits the crude Cox model with z_precip_12m as predictor.
Reports hazard ratio, 95% CI, and p-value. Tests the
proportional hazards (PH) assumption using Schoenfeld residuals
and saves diagnostics.
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
          2: Unadjusted Cox model

Efron method for tied event times
==================================================*/

stcox z_precip_12m, efron nohr   // log-hazard scale
stcox z_precip_12m, efron        // hazard ratio scale

/*==================================================
          3: Proportional hazards assumption test

estat phtest uses Schoenfeld residuals to test whether the
hazard ratio is constant over time. A significant p-value (p<0.05)
suggests that the effect of z_precip_12m varies with survival time
==================================================*/

* Re-fit storing Schoenfeld residuals for the diagnostic plot
stcox z_precip_12m, efron schoenfeld(sch_z) scaledsch(ssch_z)

log using `"`tbl_dir'\cox_phtest_un.txt"', text replace
estat phtest, detail
log close

/*==================================================
          4: Schoenfeld residuals plot

Scaled Schoenfeld residuals vs ln(failure time) with a LOWESS smoother
Note: A flat trend supports the proportional hazards assumption
==================================================*/

gen double lnt = log(_t) if _d == 1   // ln(time) at event observations only

twoway ///
	(scatter ssch_z lnt, msize(vsmall) mcolor(gs10%40)) ///
	(lowess  ssch_z lnt, lcolor(navy) lwidth(medthick)) ///
	, yline(0, lpattern(dash) lcolor(red)) ///
	xtitle("ln(failure time)") ///
	ytitle("Scaled Schoenfeld residual") ///
	title("Schoenfeld Residuals: z_precip_12m") ///
	legend(order(2 "LOWESS smoother") pos(1) ring(0)) ///
	note("Unadjusted Cox model")
graph export `"`fig_dir'\schoenfeld_un.png"', replace width(1200)

/*==================================================
          5: Baseline cumulative hazard
==================================================*/

stcox z_precip_12m, efron basehc(H0)

stcurve, cumhaz ///
	title("Baseline Cumulative Hazard") ///
	subtitle("Unadjusted Cox - z_precip_12m") ///
	xtitle("Age (months)") ytitle("Cumulative hazard H(t)") ///
	xlabel(0(12)60)
graph export `"`fig_dir'\cox_basehazard.png"', replace width(1200)
