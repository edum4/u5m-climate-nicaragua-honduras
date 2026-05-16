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
          Exploratory analyisis and Kaplan-Meier curves

Descriptive statistics by country and wave, Kaplan-Meier
curves by key subgroups, log-rank tests, and log-log plots to inspect 
the proportional hazards assumption before computing Cox models.
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
          1: Load data and prepare variables
==================================================*/

import delimited using `"`in_csv'"', varnames(1) encoding(UTF-8) clear

gen byte country_n = (country_str == "HN")
label define country_lbl 0 "Nicaragua" 1 "Honduras", replace
label values country_n country_lbl

xtile precip_q4 = z_precip_12m, nquantiles(4)
label define q4lbl 1 "Q1 (driest)" 2 "Q2" 3 "Q3" 4 "Q4 (wettest)", replace
label values precip_q4 q4lbl

foreach v in kidsex urban wealthq {
	encode `v', gen(`v'_n)
	drop `v'
	rename `v'_n `v'
}

recode wealthq (1=3) (3=1)
label define wealthq 1 "Poorest" 2 "Poorer" 3 "Middle" 4 "Richer" 5 "Richest", replace
label values wealthq wealthq

/*==================================================
          2: Declare survival time
==================================================*/

* failure(event==1): died within 60 months; scale(1): time in months
stset surv_time, failure(event==1) scale(1)
local N_st = r(N_sub)

/*==================================================
          3: Descriptive statistics
==================================================*/

tabstat surv_time event z_precip_12m precip_12m age edyrtotal ///
, stat(n mean sd min max) by(country_str) nototal long

foreach v of varlist kidsex wealthq urban {
	tab country_str `v', row nofreq
}

tabstat event, by(year) stat(n sum mean) nototal
tabstat z_precip_12m, by(year) stat(n mean sd p25 p50 p75) nototal

* Table: key variables
quietly table (country_str year), ///
statistic(count event)                          ///
statistic(sum event)                            ///
statistic(mean event age edyrtotal precip_12m) ///
nformat(%9.3f) name(table1)
collect preview
collect export `"`tbl_dir'/table_descriptive_stata.xlsx"', name(table1) replace

/*==================================================
          4: Kaplan-Meier curves
          Cumulative failure (1 - survival)
==================================================*/

* Overall
sts graph, failure ///
title("Under-5 survival, pooled sample") ///
xtitle("Age (months)") ytitle("Cumulative mortality") ///
ylabel(0(.02).10, format(%4.2f)) xlabel(0(12)60) ///
note("DHS Nicaragua (1998,2001) & Honduras (2005,2011)")
graph export `"`fig_dir'\km_overall_stata.png"', replace width(1200)

local title1 "Country"
local title2 "Precipitation quartile"
local title3 "Wealth quintile"
local title4 "Area of residence"
local ymax1  .10
local ymax2  .10
local ymax3  .14
local ymax4  .10
local sub2   `"subtitle("Q1 = driest, Q4 = wettest (z_precip_12m quartiles)")"'
local fname1 country
local fname2 precip_quartile
local fname3 wealth
local fname4 urban
local leg1   `"legend(order(1 "Nicaragua" 2 "Honduras") pos(11) ring(0))"'
local leg2   `"legend(order(1 "Q1 (driest)" 2 "Q2" 3 "Q3" 4 "Q4 (wettest)") pos(11) ring(0) cols(1))"'
local leg3   `"legend(order(1 "Poorest" 2 "Poorer" 3 "Middle" 4 "Richer" 5 "Richest") pos(11) ring(0) cols(1))"'
local leg4   `"legend(order(1 "Rural" 2 "Urban") pos(11) ring(0))"'

local i = 1
foreach v in country_n precip_q4 wealthq urban {
sts graph, by(`v') failure ///
title("Under-5 cumulative mortality by `title`i''") ///
`sub`i'' ///
xtitle("Age (months)") ytitle("Cumulative mortality") ///
ylabel(0(.02)`ymax`i'', format(%4.2f)) xlabel(0(12)60) ///
`leg`i''
graph export `"`fig_dir'\km_by_`fname`i''_stata.png"', replace width(1200)
local i = `i' + 1
}

/*==================================================
          5: Log-rank tests
==================================================*/

log using `"`tbl_dir'\logrank_tests_stata.txt"', text replace
foreach v in country_n precip_q4 wealthq urban kidsex {
	sts test `v'
}
log close

/*==================================================
          6: Log-log plots (PH assumption visual check)
==================================================*/

local ll_title1 "Survival by Country"
local ll_title2 "Survival by Precipitation Quartile"
local ll_leg1   `"legend(order(1 "Nicaragua" 2 "Honduras") pos(1) ring(0))"'
local ll_leg2   `"legend(order(1 "Q1 (driest)" 2 "Q2" 3 "Q3" 4 "Q4 (wettest)") pos(1) ring(0) cols(1))"'
local ll_fname1 country
local ll_fname2 precip_quartile

local i = 1
foreach v in country_n precip_q4 {
	stphplot, by(`v') ///
	title("Log-Log Plot: `ll_title`i''") ///
	subtitle("Parallel lines support PH assumption") ///
	xtitle("ln(time)") ytitle("ln[-ln(S(t))]") ///
	`ll_leg`i''
	graph export `"`fig_dir'\loglog_`ll_fname`i''_stata.png"', replace width(1200)
	local i = `i' + 1
}
