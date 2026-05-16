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
#         Sensitivity analyses
#
# Five alternative specifications test the robustness of the
# z_precip_12m–mortality association from Model 4:
#   S1 — Country×wave stratified Cox: separate baseline hazard per survey
#   S2 — Nicaragua only (1998 & 2001)
#   S3 — Honduras only (2005 & 2011)
#   S4 — Survey-weighted Cox using DHS person weights
#   S5 — Time-varying coefficients for z_precip_12m, edyrtotal, and wealthq
#        (all three violated the PH test in M4)
# ==================================================

library(survival)
library(broom)
library(dplyr)
library(purrr)
library(ggplot2)
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
                        labels = c("Poorest", "Poorer", "Middle", "Richer", "Richest")),
    # Wave variable for S1 stratification
    wave4      = factor(
      case_when(
        country_str == "NIC" & year == "1998" ~ "NIC-1998",
        country_str == "NIC" & year == "2001" ~ "NIC-2001",
        country_str == "HN"  & year == "2005" ~ "HN-2005",
        country_str == "HN"  & year == "2011" ~ "HN-2011"
      ),
      levels = c("NIC-1998", "NIC-2001", "HN-2005", "HN-2011")
    )
  )

cat("\n--- Wave distribution ---\n")
print(table(df$wave4, useNA = "ifany"))

# ==================================================
#         2: S1–S4 sensitivity models
#
# covs_m3: child + maternal + SES covariates, no country/wave
# ==================================================

covs_m3 <- "kidsex + age + edyrtotal + wealthq + urban"

# S1: stratified on wave4 — separate baseline hazard per survey wave
S1 <- coxph(
  as.formula(paste("Surv(surv_time, event == 1) ~ z_precip_12m +",
                   covs_m3, "+ strata(wave4)")),
  data = df, ties = "efron"
)

# S2: Nicaragua only
S2 <- coxph(
  as.formula(paste("Surv(surv_time, event == 1) ~ z_precip_12m +",
                   covs_m3, "+ year")),
  data = filter(df, country_str == "NIC"), ties = "efron"
)

# S3: Honduras only
S3 <- coxph(
  as.formula(paste("Surv(surv_time, event == 1) ~ z_precip_12m +",
                   covs_m3, "+ year")),
  data = filter(df, country_str == "HN"), ties = "efron"
)

# S4: Survey-weighted Cox using DHS person weights.
S4 <- coxph(
  as.formula(paste("Surv(surv_time, event == 1) ~ z_precip_12m +",
                   covs_m3, "+ country_n + year")),
  data    = df,
  ties    = "efron",
  weights = df$perweight,
  robust  = TRUE
)

# Print summaries
for (s in list(S1 = S1, S2 = S2, S3 = S3, S4 = S4)) {
  print(summary(s))
}

zph_S4 <- cox.zph(S4, transform = "identity")

sink(file.path(tbl_dir, "cox_phtest_s.txt"))
cat("=== PH assumption test: sensitivity model S4 (weighted) ===\n\n")
print(zph_S4)
sink()

cat("\n--- PH test: S4 ---\n")
print(zph_S4)

# ==================================================
#         4: S5 — Time-varying coefficients for z_precip_12m,
#                 edyrtotal, and wealthq (all violated PH test)
# ==================================================

df <- df |> mutate(wealthq_num = as.integer(wealthq))

S5 <- coxph(
  as.formula(paste("Surv(surv_time, event == 1) ~ z_precip_12m +",
                   covs_m3, "+ tt(z_precip_12m) + tt(edyrtotal) + tt(wealthq_num) + country_n + year")),
  data = df,
  ties = "efron",
  tt   = function(x, t, ...) x * log(pmax(t, 0.5))  # x * ln(t)
)

cat("\n--- S5: Time-varying coefficients for z_precip_12m, edyrtotal, wealthq ---\n")
print(summary(S5))

# ==================================================
#         5: Summary table — z_precip_12m HR across S1–S4
#            S5 reported separately (TVC model)
# ==================================================

cat("\n--- HR for z_precip_12m: S1–S4 ---\n")
imap(list(S1 = S1, S2 = S2, S3 = S3, S4 = S4),
     function(fit, name) {
       tidy(fit, exponentiate = TRUE, conf.int = TRUE) |>
         filter(term == "z_precip_12m") |>
         mutate(model = name, n = fit$n, n_events = fit$nevent) |>
         select(model, n, n_events, HR = estimate,
                conf.low, conf.high, p.value)
     }) |>
  bind_rows() |>
  mutate(across(where(is.double), \(x) round(x, 4))) |>
  print()

cat("\n--- HR for z_precip_12m: S5 (TVC) ---\n")
tidy(S5, exponentiate = TRUE, conf.int = TRUE) |>
  filter(term == "z_precip_12m") |>
  mutate(across(where(is.double), \(x) round(x, 4))) |>
  print()

# ==================================================
#         6: Export results table to Excel
# ==================================================

tbl_sens <- imap(list(S1 = S1, S2 = S2, S3 = S3, S4 = S4),
                 function(fit, name) {
                   tidy(fit, exponentiate = TRUE, conf.int = TRUE) |>
                     filter(term == "z_precip_12m") |>
                     mutate(
                       model    = name,
                       n        = fit$n,
                       n_events = fit$nevent,
                       stars    = case_when(
                         p.value < 0.001 ~ "***",
                         p.value < 0.01  ~ "**",
                         p.value < 0.05  ~ "*",
                         TRUE            ~ ""
                       )
                     ) |>
                     select(model, n, n_events,
                            HR = estimate, SE = std.error,
                            conf.low, conf.high, p.value, stars)
                 }) |>
  bind_rows() |>
  mutate(across(c(HR, SE, conf.low, conf.high), \(x) round(x, 3)),
         p.value = round(p.value, 4))

# S5 appended as a separate block
tbl_S5 <- tidy(S5, exponentiate = TRUE, conf.int = TRUE) |>
  filter(term == "z_precip_12m") |>
  mutate(
    model    = "S5 (TVC)",
    n        = S5$n,
    n_events = S5$nevent,
    stars    = case_when(
      p.value < 0.001 ~ "***",
      p.value < 0.01  ~ "**",
      p.value < 0.05  ~ "*",
      TRUE            ~ ""
    )
  ) |>
  select(model, n, n_events,
         HR = estimate, SE = std.error,
         conf.low, conf.high, p.value, stars) |>
  mutate(across(c(HR, SE, conf.low, conf.high), \(x) round(x, 3)),
         p.value = round(p.value, 4))

write_xlsx(
  list("S1-S4" = tbl_sens, "S5_TVC" = tbl_S5),
  file.path(tbl_dir, "cox_sensitivity.xlsx")
)

# ==================================================
#         7: Forest plot — S1–S5
# ==================================================

forest_sens <- imap(
  list(S1 = S1, S2 = S2, S3 = S3, S4 = S4, S5 = S5),
  function(fit, name) {
    tidy(fit, exponentiate = TRUE, conf.int = TRUE) |>
      filter(term == "z_precip_12m") |>
      mutate(model = name)
  }) |>
  bind_rows() |>
  mutate(
    model = factor(model, levels = rev(c("S1", "S2", "S3", "S4", "S5"))),
    label = sprintf("%.4f", estimate)
  )

p_forest_s <- ggplot(forest_sens,
                     aes(x = estimate, y = model,
                         xmin = conf.low, xmax = conf.high)) +
  geom_pointrange(linewidth = 0.7, size = 0.5, colour = "steelblue4") +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "red") +
  geom_text(aes(label = label), vjust = -0.7, size = 3.2, colour = "grey30") +
  scale_x_continuous(breaks = seq(0.7, 1.8, 0.2),
                     labels = scales::label_number(accuracy = 0.1)) +
  labs(
    title    = "Sensitivity analyses: z_precip_12m HR",
    subtitle = "Alternative model specifications",
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

ggsave(file.path(fig_dir, "forest_s.png"), p_forest_s,
       width = 8, height = 5, dpi = 150)
