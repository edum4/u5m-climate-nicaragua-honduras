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
#         Cox models
#
# Four nested models assess whether the crude association
# between z_precip_12m and under-5 mortality persists after
# covariate adjustment:
#   M0 --> Unadjusted (baseline)
#   M1 --> + Child characteristics:      kidsex
#   M2 --> + Maternal characteristics:   age, edyrtotal
#   M3 --> + Socioeconomic status:       wealthq, urban
#   M4 --> Full (+ context):             country_n, year
# ==================================================

library(survival)
library(broom)
library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)
library(writexl)

# ==================================================
#         0: Path configuration
# ==================================================

root    <- ""  # ← INSERT YOUR ROOT PATH HERE (all other paths build from this)
in_csv  <- file.path(root, "Child mortality and climate", "Data", "processed", "dhs_final.csv")
tbl_dir <- file.path(root, "Child mortality and climate", "Output", "tables", "R")
fig_dir <- file.path(root, "Child mortality and climate", "Output", "figures", "R")

# ==================================================
#         1: Load and prepare data
# ==================================================

df <- read.csv(in_csv, stringsAsFactors = FALSE)

df <- df |>
  mutate(
    country_n  = factor(country_str,
                        levels = c("NIC", "HN"),
                        labels = c("Nicaragua", "Honduras")),
    year       = factor(year),
    kidsex     = factor(tolower(kidsex),
                        levels = c("male", "female"),
                        labels = c("Male", "Female")),
    urban      = factor(tolower(urban),
                        levels = c("urban", "rural"),
                        labels = c("Urban", "Rural")),
    wealthq    = factor(tolower(wealthq),
                        levels = c("poorest", "poorer", "middle", "richer", "richest"),
                        labels = c("Poorest", "Poorer", "Middle", "Richer", "Richest"))
  )

# ==================================================
#         2: Cox models
#
# Efron method for tied event times. 
# M0 = unadjusted (baseline),
# M1–M4 = progressive covariate adjustment.
# ==================================================

rhs <- list(
  M0 = "z_precip_12m",
  M1 = "z_precip_12m + kidsex",
  M2 = "z_precip_12m + kidsex + age + edyrtotal",
  M3 = "z_precip_12m + kidsex + age + edyrtotal + wealthq + urban",
  M4 = "z_precip_12m + kidsex + age + edyrtotal + wealthq + urban + country_n + year"
)

fits <- lapply(rhs, function(r) {
  coxph(as.formula(paste("Surv(surv_time, event == 1) ~", r)),
        data = df, ties = "efron")
})

# Print summary and concordance for each model
for (m in names(fits)) {
  cat(sprintf("\n=== %s ===\n", m))
  print(summary(fits[[m]]))
}

# ==================================================
#         3: Summary table — z_precip_12m HR across models
# ==================================================

cat("\n--- HR for z_precip_12m across models ---\n")
imap(fits, function(fit, name) {
  tidy(fit, exponentiate = TRUE, conf.int = TRUE) |>
    filter(term == "z_precip_12m") |>
    mutate(
      model    = name,
      n        = fit$n,
      n_events = fit$nevent,
      concordance = summary(fit)$concordance["C"]
    ) |>
    select(model, n, n_events, concordance,
           HR = estimate, conf.low, conf.high, p.value)
}) |>
  bind_rows() |>
  mutate(across(where(is.double), \(x) round(x, 4))) |>
  print()

# ==================================================
#         4: Export results table to Excel
# ==================================================

tbl_models <- imap(fits, function(fit, name) {
  tidy(fit, exponentiate = TRUE, conf.int = TRUE) |>
    filter(term == "z_precip_12m") |>
    mutate(
      model       = name,
      n           = fit$n,
      n_events    = fit$nevent,
      concordance = round(summary(fit)$concordance["C"], 3),
      se_log_hr   = std.error,              
      stars       = case_when(
        p.value < 0.001 ~ "***",
        p.value < 0.01  ~ "**",
        p.value < 0.05  ~ "*",
        TRUE            ~ ""
      )
    ) |>
    select(model, n, n_events, concordance,
           HR = estimate, SE = se_log_hr,
           conf.low, conf.high, p.value, stars)
}) |>
  bind_rows() |>
  mutate(across(c(HR, SE, conf.low, conf.high), \(x) round(x, 3)),
         p.value = round(p.value, 4))

write_xlsx(tbl_models, file.path(tbl_dir, "cox_models.xlsx"))

# ==================================================
#         5: PH assumption test — full model (M4)
#
# Schoenfeld residuals test on the fully-adjusted model.
# A significant p-value (p < 0.05) suggests
# non-proportional hazards.
# ==================================================

zph_M4 <- cox.zph(fits[["M4"]], transform = "identity")

sink(file.path(tbl_dir, "cox_phtest_ad.txt"))
cat("=== PH assumption test: full model M4 ===\n\n")
print(zph_M4)
sink()

cat("\n--- PH test: full model M4 ---\n")
print(zph_M4)

# ==================================================
#         6: Forest plot — z_precip_12m HR across models
# ==================================================

forest_df <- imap(fits, function(fit, name) {
  tidy(fit, exponentiate = TRUE, conf.int = TRUE) |>
    filter(term == "z_precip_12m") |>
    mutate(model = name)
}) |>
  bind_rows() |>
  mutate(
    model  = factor(model, levels = rev(names(fits))),  # M4 on top, M0 at bottom
    label  = sprintf("%.4f", estimate)
  )

p_forest <- ggplot(forest_df,
                   aes(x = estimate, y = model,
                       xmin = conf.low, xmax = conf.high)) +
  geom_pointrange(linewidth = 0.7, size = 0.5, colour = "steelblue4") +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "red") +
  geom_text(aes(label = label), vjust = -0.7, size = 3.2, colour = "grey30") +
  scale_x_continuous(breaks = seq(0.9, 1.5, 0.1),
                     labels = scales::label_number(accuracy = 0.1)) +
  labs(
    title    = "HR for 12-month precipitation z-score",
    subtitle = "Sequential covariate adjustment — under-5 mortality",
    x        = "Hazard ratio (z_precip_12m)",
    y        = NULL,
    caption  = "Dashed line = null (HR = 1). Error bars = 95% CI."
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title    = element_text(face = "bold"),
    plot.subtitle = element_text(colour = "grey40", size = 10),
    plot.caption  = element_text(colour = "grey50", size = 9),
    axis.text.y   = element_text(face = "bold")
  )

ggsave(file.path(fig_dir, "forestplot.png"), p_forest,
       width = 8, height = 5, dpi = 150)
