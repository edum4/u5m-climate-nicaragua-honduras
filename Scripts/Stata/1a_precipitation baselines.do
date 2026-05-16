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
          Rolling 15-year precipitation baselines

For each child's birth in year Y, the individual climate baseline is the
per-region × per-calendar-month mean and SD of precipitation computed
over the 15 years prior (Y-15 to Y-1), or all available years back to
1980 when fewer than 15 are found. With WorldClim data starting in 1980,
partial windows affect only the earliest cohorts: birth_year=1993 uses
1980–1992 (13 years) and birth_year=1994 uses 1980–1993 (14 years).
Full 15-year baselines are available from birth_year=1995 onward.
==================================================*/

clear all
set more off
set type double

/*==================================================
          0: Path configuration
==================================================*/

local root `""'  // ← INSERT YOUR ROOT PATH HERE (all other paths build from this)

local in_csv  `"`root'\Child mortality and climate\Data\processed\worldclim_regional_monthly.csv"'
local out_csv `"`root'\Child mortality and climate\Data\processed\worldclim_historical_baselines.csv"'
local out_dta `"`root'\Child mortality and climate\Data\processed\worldclim_historical_baselines.dta"'

* Birth year range (from DHS kidbirthyr) and earliest WorldClim year available in the files
local y_min      = 1993
local y_max      = 2012
local data_start = 1980

/*==================================================
          1: Load WorldClim monthly data
==================================================*/

import delimited using `"`in_csv'"', varnames(1) encoding(UTF-8) clear

tab country
sum precip_mm

/*==================================================
          2: Rolling baseline

For each birth in year Y:
- Compute window boundaries: start = max(data_start, Y-15), end = Y-1
- Preserve the full dataset
- Keep months within the window, collapse data to region × month stats
- Tag with birth-year data and append to a tempfile
- Restore for the next iteration
==================================================*/

tempfile rolling
preserve
	clear
	save `"`rolling'"', emptyok replace
restore

forvalues Y = `y_min'/`y_max' {
	local y_end   = `Y' - 1
	local y_start = max(`data_start', `Y' - 15)
	local n_yrs   = `y_end' - `y_start' + 1
	preserve
		keep if inrange(year, `y_start', `y_end')
		collapse ///
		(mean) baseline_mean = precip_mm ///
		(sd)   baseline_sd   = precip_mm ///
		, by(country dhs_code month)
		gen int birth_year          = `Y'
		gen int n_years_used        = `n_yrs'
		gen int baseline_year_start = `y_start'
		gen int baseline_year_end   = `y_end'
		append using `"`rolling'"'
		save `"`rolling'"', replace
	restore
}

/*==================================================
          3: Check results
==================================================*/

use `"`rolling'"', clear
sort country dhs_code birth_year month

tab birth_year n_years_used, missing // Checking windows, particularly for early cohorts
sum baseline_mean baseline_sd, d

list country dhs_code birth_year month baseline_mean baseline_sd n_years_used baseline_year_start baseline_year_end ///
     if country=="HN" & dhs_code==5 & month==1 & birth_year==2005, noobs abbreviate(25)

/*==================================================
          5: Save new file
==================================================*/

order country dhs_code birth_year month baseline_mean baseline_sd n_years_used baseline_year_start baseline_year_end
save `"`out_dta'"', replace                                          
export delimited using `"`out_csv'"', delimiter(",") replace         