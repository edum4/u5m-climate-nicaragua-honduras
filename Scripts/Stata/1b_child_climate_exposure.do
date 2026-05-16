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
          Climate exposure per-child

For each child in DHS data, monthly precipitation is extracted for the
12 months following birth and standardizes each month against its
calendar-month-specific rolling baseline. The output variable
z_precip_12m is the mean of 12 monthly z-scores: an index of how different
the child's first year was relative to recent regional precipitation.

Steps:
- Prepare DHS data for merge (identifiers + birth date + region)
- Expand x12: one row per child × month of life (offset 0–11)
- Derive calendar (year, month) for each offset
- Merge with worldclim_regional_monthly on (country_str, dhs_code, year, month)
- Merge with worldclim_historical_baselines on (country_str, dhs_code, birth_year, month)
- z_m = (precip_mm - baseline_mean) / baseline_sd
- Collapse: precip_12m = mean(precip_mm), z_precip_12m = mean(z_m)
==================================================*/

clear all
set more off
set type double

/*==================================================
          0: Path configuration
==================================================*/

local root `""'  // ← INSERT YOUR ROOT PATH HERE (all other paths build from this)

local dhs_dta  `"`root'\Child mortality and climate\Data\IPUMS-DHS\dhs_nic_hn_4waves.dta"'
local clim_csv `"`root'\Child mortality and climate\Data\processed\worldclim_regional_monthly.csv"'
local base_dta `"`root'\Child mortality and climate\Data\processed\worldclim_historical_baselines.dta"'
local out_csv  `"`root'\Child mortality and climate\Data\processed\dhs_climate_exposure.csv"'

/*==================================================
          1: Prepare DHS data for merge
==================================================*/

use `"`dhs_dta'"', clear
keep idhspid bidx country kidbirthmo kidbirthyr geo_hn2005_2011 geo_nc1998_2001
gen str3 country_str = cond(country == 340, "HN", "NIC") // Country string
gen int dhs_code = round(cond(country == 340, geo_hn2005_2011, geo_nc1998_2001)) // Unified region code
keep idhspid bidx country_str dhs_code kidbirthmo kidbirthyr
tempfile dhs_m
save `"`dhs_m'"'

/*==================================================
          2: Load WorldClim monthly data
==================================================*/

import delimited using `"`clim_csv'"', varnames(1) encoding(UTF-8) clear
keep country dhs_code year month precip_mm
rename country country_str
destring dhs_code year month precip_mm, replace
sort country_str dhs_code year month
tempfile clim
save `"`clim'"'

/*==================================================
          3: Expand to 12 rows per child and compute window dates
==================================================*/

use `"`dhs_m'"', clear
expand 12

bysort idhspid bidx: gen int offset = _n - 1  // 0 = birth month, 11 = last month

gen int month = mod(kidbirthmo - 1 + offset, 12) + 1 // mod() wraps the month back to 1 after December
gen int year  = kidbirthyr + floor((kidbirthmo - 1 + offset) / 12) // floor() counts full years elapsed

/*==================================================
          4: Merge with WorldClim monthly precipitation
==================================================*/

sort country_str dhs_code year month
merge m:1 country_str dhs_code year month using `"`clim'"', keepusing(precip_mm) keep(match master)
tab _merge
drop _merge

/*==================================================
          5: Merge with rolling baselines
==================================================*/

gen int birth_year = kidbirthyr // birth_year is merge key
rename country_str country // baseline_dta stores the country variable as "country" 
sort country dhs_code birth_year month
merge m:1 country dhs_code birth_year month using `"`base_dta'"', ///
keepusing(baseline_mean baseline_sd n_years_used) keep(match master)
tab _merge
drop _merge birth_year
rename country country_str

/*==================================================
          6: Monthly z-scores
==================================================*/

gen double z_monthly = (precip_mm - baseline_mean) / baseline_sd
replace z_monthly = 0 if baseline_sd == 0 | missing(baseline_sd) // if baseline_sd == 0 (no variance over baseline window), z is set to 0

/*==================================================
          7: Collapse to one row per child
==================================================*/

collapse ///
(mean)  precip_12m       = precip_mm    ///
(mean)  z_precip_12m     = z_monthly    ///
(mean)  n_years_avg      = n_years_used ///
(count) n_months_matched = precip_mm    ///
(first) country_str                     ///
(first) dhs_code                        ///
(first) kidbirthmo                      ///
(first) kidbirthyr                      ///
, by(idhspid bidx)

/*==================================================
          8: Checking results
==================================================*/

sum z_precip_12m precip_12m, d
tabstat z_precip_12m precip_12m, by(country_str) stat(n mean sd min max) nototal

/*==================================================
          9: Save new file
==================================================*/

order idhspid bidx country_str dhs_code kidbirthmo kidbirthyr precip_12m z_precip_12m n_months_matched n_years_avg

export delimited using `"`out_csv'"', delimiter(",") replace
