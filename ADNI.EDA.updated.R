# ADNI Exploratory Data Analysis (EDA)
# Group 7 (Brain Bugs)
# Author: Paa Kwesi Danso.
# Date: Feb 2026
#
# Purpose:
#   1) Analysis 1: CSF Aβ42 vs MMSE (regression + diagnostics)
#   2) Analysis 2: APOE ε4 status vs MMSE trajectories (0–60 months)
#      (mean plot + spaghetti + boxplot + t-test + mixed-effects)

library(tidyverse)
library(lme4)
library(lmerTest)

# Set working directory
setwd("C:/Users/quacy/Downloads/ADNI_EDA_data")

# Load Data (All 6 Tables)
adas       <- read_csv("All_Subjects_ADAS_17Feb2026.csv", show_col_types = FALSE)
apoe       <- read_csv("All_Subjects_APOERES_17Feb2026.csv", show_col_types = FALSE)
cdr        <- read_csv("All_Subjects_CDR_17Feb2026.csv", show_col_types = FALSE)
mmse       <- read_csv("All_Subjects_MMSE_17Feb2026.csv", show_col_types = FALSE)
demog      <- read_csv("All_Subjects_PTDEMOG_17Feb2026.csv", show_col_types = FALSE)
biomarkers <- read_csv("All_Subjects_UPENNBIOMK_ROCHE_ELECSYS_17Feb2026.csv", show_col_types = FALSE)


mmse <- mmse %>%
  mutate(MMSCORE = if_else(MMSCORE == -1, NA_real_, MMSCORE))

# Helper: VISCODE2 -> months
# bl -> 0, m06 -> 6, m12 -> 12
viscode_to_months <- function(v) {
  case_when(
    v == "bl" ~ 0,
    str_detect(v, "^m\\d+$") ~ as.numeric(str_remove(v, "m")),
    TRUE ~ NA_real_
  )
}


# Analysis 1: CSF Aβ42 and Tau vs MMSE (matched by RID + VISCODE2)

analysis1 <- mmse %>%
  select(RID, VISCODE2, MMSCORE) %>%
  inner_join(
    biomarkers %>% select(RID, VISCODE2, ABETA42, TAU),
    by = c("RID", "VISCODE2")
  ) %>%
  drop_na(MMSCORE, ABETA42, TAU)

cat("Analysis 1 rows:", nrow(analysis1), "\n")
cat("Analysis 1 unique participants:", n_distinct(analysis1$RID), "\n\n")

# Plot: scatter + regression line (Aβ42 vs MMSE)
plot_analysis1 <- ggplot(analysis1, aes(x = ABETA42, y = MMSCORE)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm", se = TRUE) +
  labs(
    title = "CSF Aβ42 vs MMSE",
    x = "CSF Aβ42",
    y = "MMSE Score"
  ) +
  theme_minimal()

print(plot_analysis1)

# Regression models
# Model 1: Aβ42 only
model1_abeta <- lm(MMSCORE ~ ABETA42, data = analysis1)

# Model 2: Aβ42 + Tau
model2_abeta_tau <- lm(MMSCORE ~ ABETA42 + TAU, data = analysis1)

cat("\n--- Model 1: MMSE ~ ABETA42 ---\n")
print(summary(model1_abeta))

cat("\n--- Model 2: MMSE ~ ABETA42 + TAU ---\n")
print(summary(model2_abeta_tau))

# Compare overall fit (useful for your manuscript tables)
cat("\nAIC comparison:\n")
print(AIC(model1_abeta, model2_abeta_tau))

cat("\nNested model comparison (ANOVA):\n")
print(anova(model1_abeta, model2_abeta_tau))

# Extra checks requested by feedback (simple + useful)
cat("\nCorrelation (ABETA42 vs MMSCORE): ",
    round(cor(analysis1$ABETA42, analysis1$MMSCORE, use = "complete.obs"), 3), "\n")

cat("Proportion of MMSE == 30 (ceiling effect): ",
    round(mean(analysis1$MMSCORE == 30, na.rm = TRUE), 3), "\n\n")

# Regression diagnostics (assumption checks)
par(mfrow = c(2, 2))
plot(model1_abeta)
par(mfrow = c(2, 2))
plot(model2_abeta_tau)
par(mfrow = c(1, 1))


# Analysis 2: APOE ε4 status vs MMSE trajectories (0–60 months)


# APOE genotype does not change across visits, so keep one row per RID
apoe_status <- apoe %>%
  select(RID, GENOTYPE) %>%
  distinct() %>%
  mutate(apoe4_status = if_else(str_detect(GENOTYPE, "4"), "Carrier", "Non-carrier"))

analysis2 <- mmse %>%
  select(RID, VISCODE2, MMSCORE) %>%
  inner_join(apoe_status, by = "RID") %>%
  mutate(
    VISCODE2 = trimws(tolower(VISCODE2)),
    months = suppressWarnings(as.numeric(sub("^m", "", VISCODE2))),
    MMSCORE = ifelse(MMSCORE == -1, NA, MMSCORE)
  ) %>%
  drop_na(MMSCORE, months) %>%
  filter(months <= 60)

# Set "Carrier" as the reference group
analysis2$apoe4_status <- factor(analysis2$apoe4_status, levels = c("Carrier", "Non-carrier"))

cat("Analysis 2 rows (<=60 months):", nrow(analysis2), "\n")
cat("Analysis 2 unique participants:", n_distinct(analysis2$RID), "\n\n")

# Mean trajectory plot
plot_analysis2_mean <- ggplot(analysis2, aes(x = months, y = MMSCORE, color = apoe4_status)) +
  stat_summary(fun = mean, geom = "line", linewidth = 1.2) +
  stat_summary(fun = mean, geom = "point", size = 2) +
  labs(
    title = "MMSE Trajectory by APOE ε4 Status (First 60 Months)",
    x = "Months Since Baseline",
    y = "MMSE Score",
    color = "APOE ε4 Status"
  ) +
  scale_y_continuous(limits = c(0, 30), breaks = seq(0, 30, by = 5)) +
  theme_minimal()

print(plot_analysis2_mean)

# Spaghetti plot (sample participants so the plot is readable)
set.seed(7)
sample_rids <- sample(unique(analysis2$RID), size = min(150, n_distinct(analysis2$RID)))
analysis2_sample <- analysis2 %>% filter(RID %in% sample_rids)

plot_analysis2_spaghetti <- ggplot(analysis2_sample, aes(x = months, y = MMSCORE, group = RID, color = apoe4_status)) +
  geom_line(alpha = 0.35, linewidth = 0.5) +
  labs(
    title = "Individual MMSE Trajectories (Sample) by APOE ε4 Status",
    x = "Months Since Baseline",
    y = "MMSE Score",
    color = "APOE ε4 Status"
  ) +
  theme_minimal()

print(plot_analysis2_spaghetti)

# Boxplot (distribution by group)
plot_analysis2_box <- ggplot(analysis2, aes(x = apoe4_status, y = MMSCORE, fill = apoe4_status)) +
  geom_boxplot(outlier.alpha = 0.6) +
  labs(
    title = "Distribution of MMSE Scores by APOE ε4 Status (0–60 Months)",
    x = "APOE ε4 Status",
    y = "MMSE Score"
  ) +
  theme_minimal() +
  theme(legend.position = "none")

print(plot_analysis2_box)

# Summary + t-test (exploratory group difference)
apoe_summary <- analysis2 %>%
  group_by(apoe4_status) %>%
  summarise(
    mean_mmse = mean(MMSCORE, na.rm = TRUE),
    sd_mmse   = sd(MMSCORE, na.rm = TRUE),
    n         = n(),
    .groups = "drop"
  )
print(apoe_summary)

print(t.test(MMSCORE ~ apoe4_status, data = analysis2))

# Mixed-effects model (handles repeated measures)
mixed_model1 <- lmer(
  MMSCORE ~ months + apoe4_status + months:apoe4_status + (1 | RID),
  data = analysis2,
  REML = TRUE
)
print(summary(mixed_model1))
print(fixef(mixed_model1))