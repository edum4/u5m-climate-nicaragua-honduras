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
          Merge DHS + climate exposure

Merges the full DHS data with the per-child climate exposure data, 
constructs survival outcome variables, recodes sentinel codes,
adds region names, and saves the dataset for subsequent tasks.
==================================================*/

clear all
set more off
set type double

/*==================================================
          0: Path configuration
==================================================*/

local root    `""'  // ← INSERT YOUR ROOT PATH HERE (all other paths build from this)
local dhs_dta `"`root'\Child mortality and climate\Data\IPUMS-DHS\dhs_nic_hn_4waves.dta"'
local exp_csv `"`root'\Child mortality and climate\Data\processed\dhs_climate_exposure.csv"'
local reg_csv `"`root'\Child mortality and climate\Data\processed\region_mapping.csv"'
local out_csv `"`root'\Child mortality and climate\Data\processed\dhs_final.csv"'

/*==================================================
          1: Load DHS data and keep key variables
==================================================*/

use `"`dhs_dta'"', clear
keep idhspid bidx country year ///
     kidalive kidagediedimp kidagemo ///
     kidsex kidbirthmo kidbirthyr ///
     age edyrtotal wealthq urban ///
     geo_hn2005_2011 geo_nc1998_2001 ///
     perweight

replace idhspid = strtrim(idhspid) // .dta pads string variables to fixed width; trim for consistent merging with CSVs

/*==================================================
          2: Survival outcome variables

surv_time = min(age at death or current age, 60)
event     = died before 60 months
Note: The cond() handles sentinel values: kidagediedimp=999 (NIU for alive children)
==================================================*/

gen double surv_time = min(cond(kidalive==0, kidagediedimp, kidagemo), 60)
gen byte event = (kidalive==0 & kidagediedimp<=60)

* Midpoint convention: DHS records age at death in whole months, so a child
* coded as dying at month 0 lived on average 0.5 months.
replace surv_time = 0.5 if surv_time == 0

label variable surv_time "Survival time in months (0 recoded to 0.5), max 60."
label variable event "Event: 1 = died before 60 months, 0 = censored"

sum surv_time, d
tab event

/*==================================================
          3: Recode sentinel codes in edyrtotal
==================================================*/

tab edyrtotal if inlist(edyrtotal, 96, 97, 98), miss
recode edyrtotal (96=.) (97=.) (98=.)

/*==================================================
          4: New variables
==================================================*/

gen int dhs_code = round(cond(country==340, geo_hn2005_2011, geo_nc1998_2001)) // region code: round() needed because geo vars are stored as float
gen str3 country_str = cond(country==340, "HN", "NIC")
label variable dhs_code "Unified DHS region code"
label variable country_str "Country string (HN/NIC)"

/*==================================================
          5: Add region names from region_mapping.csv
==================================================*/

tempfile dhs_rm
save `"`dhs_rm'"'

import delimited using `"`reg_csv'"', varnames(1) encoding(UTF-8) clear
keep country dhs_code geojson_name
rename country country_str
rename geojson_name region_name
destring dhs_code, replace
duplicates drop country_str dhs_code, force
tempfile region_lkp
save `"`region_lkp'"'

use `"`dhs_rm'"', clear
sort country_str dhs_code
merge m:1 country_str dhs_code using `"`region_lkp'"', keepusing(region_name) keep(match master)
tab _merge
drop _merge
label variable region_name "Administrative region name"

/*==================================================
          6: Merge with climate exposure
==================================================*/

tempfile dhs_master
save `"`dhs_master'"'

import delimited using `"`exp_csv'"', varnames(1) encoding(UTF-8) clear
replace idhspid = strtrim(idhspid)
keep idhspid bidx precip_12m z_precip_12m n_years_avg n_months_matched
destring bidx precip_12m z_precip_12m n_years_avg n_months_matched, replace
tempfile climate_exp
save `"`climate_exp'"'

use `"`dhs_master'"', clear
sort idhspid bidx
merge 1:1 idhspid bidx using `"`climate_exp'"'
tab _merge
drop _merge

/*==================================================
          7: Variable labels and column order
==================================================*/

label variable idhspid          "Unique cross-sample respondent ID"
label variable bidx             "Birth order within respondent"
label variable country          "Country ISO-3166 (340=HN, 558=NIC)"
label variable year             "DHS survey year"
label variable kidalive         "Child alive at interview (0=died, 1=alive)"
label variable kidsex           "Child sex (1=male, 2=female)"
label variable kidbirthmo       "Month of birth (1-12)"
label variable kidbirthyr       "Year of birth"
label variable age              "Maternal age at interview (years)"
label variable edyrtotal        "Maternal years of education"
label variable wealthq          "Wealth quintile (1=poorest, 5=richest)"
label variable urban            "Urban/rural (1=urban, 2=rural)"
label variable perweight        "DHS individual weight"
label variable precip_12m       "Mean monthly precipitation in 12m post-birth (mm)"
label variable z_precip_12m     "Standardized 12m precipitation (mean of monthly z-scores)"
label variable n_years_avg      "Average years in rolling baseline (max 15)"
label variable n_months_matched "Climate months matched"

order idhspid bidx country country_str year dhs_code region_name ///
      kidalive surv_time event ///
      kidsex kidbirthmo kidbirthyr ///
      age edyrtotal wealthq urban perweight ///
      precip_12m z_precip_12m n_years_avg n_months_matched

drop geo_hn2005_2011 geo_nc1998_2001 kidagediedimp kidagemo

/*==================================================
          8: Check results
==================================================*/

sum surv_time z_precip_12m, d
tabstat event, by(country_str) stat(n sum mean) nototal

/*==================================================
          9: Save new file
==================================================*/

export delimited using `"`out_csv'"', delimiter(",") replace
