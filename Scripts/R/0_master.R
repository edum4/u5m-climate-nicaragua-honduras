# ==================================================
# project:       Survival Analysis – DHS + WorldClim Climate Data
#                Cox PH Model: Under-5 Child Mortality
#                Nicaragua (1998, 2001) & Honduras (2005, 2011)
# Author:        Eduardo Pacheco
# --------------------------------------------------
# Creation Date:     16 March 2026
# Modification Date: 29 April 2026
# Script version:    02
# ==================================================

# ==================================================
#         Master script
#
# Runs the modelling pipeline in order:
#   4b — Adjusted Cox models
#   4c — Sensitivity analyses
#
# Note: Data construction (scripts 1a–2b), exploratory
# analysis (script 3) and unadjusted Cox model (script 4a) 
# are conducted in Stata. This master covers only the modelling.
# ==================================================

root    <- ""  # ← INSERT YOUR ROOT PATH HERE (all other paths build from this)
scripts <- file.path(root, "Child mortality and climate", "Scripts", "R")

source(file.path(scripts, "4b_coxmodels_adjusted.R"))
source(file.path(scripts, "4c_sensitivity_analysis.R"))
