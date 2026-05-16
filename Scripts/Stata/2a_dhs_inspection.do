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
          DHS data inspection

This script documents the structure, distributions, and missing 
rates of all variables entering the survival model.

==================================================*/

clear all
set more off

/*==================================================
          0: Path configuration
==================================================*/

local root `""'  // ← INSERT YOUR ROOT PATH HERE (all other paths build from this)
local dhs  `"`root'\Child mortality and climate\Data\IPUMS-DHS\dhs_nic_hn_4waves.dta"'

/*==================================================
          1: Overview
==================================================*/

use `"`dhs'"', clear
describe
tab country   // 340 = Honduras, 558 = Nicaragua (ISO-3166)
tab year
tab country year

/*==================================================
          2: Survival outcome variables
==================================================*/

* kidalive: 0 = died, 1 = alive at interview
tab kidalive, miss

* kidagediedimp: age at death in months. 999 = NIU for alive children
sum kidagediedimp if kidalive == 0, d

* kidagemo: current age in months. 99 = NIU for dead children
sum kidagemo if kidalive == 1, d

* Survival time and event indicator (as constructed in the final model)
* surv_time = min(age at death or current age, 60); event = died before 60 months
tempvar surv_time event
gen `surv_time' = min(cond(kidalive == 0, kidagediedimp, kidagemo), 60)
gen `event'     = (kidalive == 0 & kidagediedimp <= 60)

sum `surv_time', d
tab `event'
tabstat `event', by(country) stat(n sum mean) nototal
tabstat `event', by(year)    stat(n sum mean) nototal

/*==================================================
          3: Child and maternal characteristics
==================================================*/

tab kidsex, miss
tab kidbirthmo, miss
tab kidbirthyr, miss

sum age, d // maternal age at interview

* edyrtotal: years of education. Sentinel codes: 96 = inconsistent, 97 = don't know, 98 = missing
sum edyrtotal if edyrtotal < 96, d
tab edyrtotal if edyrtotal >= 95, miss

/*==================================================
          4: Household and geographic variables
==================================================*/

tab wealthq, miss
tab urban,   miss

tab geo_hn2005_2011 if country == 340, miss   // Honduras regions
tab geo_nc1998_2001 if country == 558, miss   // Nicaragua regions

/*==================================================
          5: Missing values
==================================================*/

misstable summarize idhspid bidx country year kidalive kidagediedimp kidagemo ///
          kidsex kidbirthmo kidbirthyr edyrtotal wealthq age urban          ///
          geo_hn2005_2011 geo_nc1998_2001 perweight

/*==================================================
          6: Sentinel values
==================================================*/

tab kidagediedimp if kidagediedimp >= 990, miss // 999: NIU (alive children) 
tab kidagemo      if kidagemo      >= 95,  miss // 99:  NIU (dead children)  
tab edyrtotal     if edyrtotal     >= 95,  miss // 96/97/98: inconsistent / don't know / missing 
duplicates report idhspid bidx   // Confirming the are 0 duplicates
