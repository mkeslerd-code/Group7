# =============================================================================
# ADNI Manuscript Analysis Script
# Group 7 (Brain Bugs) | Arizona State University
# Date: March 2026
#
# Structure:
#   SECTION 0: Libraries + Data Loading
#   SECTION 1: Cross-Sectional Analysis (OLS Models 1–4)
#   SECTION 2: Longitudinal EDA (APOE ε4 Trajectories)
#   SECTION 3: Linear Mixed-Effects Models (LMM)
#   SECTION 4: Secondary Outcome — ADAS-Cog Validation
#   SECTION 5: Model Validation — Train/Test Split
#   SECTION 6: Publication Tables (gtsummary)
#
# Manuscript Figures:
#   Fig 1: plot_analysis1          — CSF Aβ42 vs MMSE (basic, unstratified)
#   Fig 2: plot_apoe_scatter       — CSF Aβ42 vs MMSE by APOE ε4 Status
#   Fig 3: plot_analysis2_mean     — Mean MMSE trajectories by APOE ε4
#   Fig 4: plot_analysis2_spaghetti — Individual MMSE trajectories by APOE ε4
#   Fig 5: plot_analysis2_box      — MMSE distribution by APOE ε4 (boxplot)
#   Fig 6: plot_random_slopes      — Predicted trajectories (random slopes model)
#   Fig 7: plot_ceiling            — MMSE vs ADAS-Cog ceiling effect
#
# Manuscript Tables:
#   Table 1: Dataset summary (reported in text)
#   Table 2: tbl_F1 — Participant characteristics by APOE ε4 status
#   Table 3: tbl_F2 — OLS Models 1–4 (cross-sectional)
#   Table 4: tbl_F3 — LMM full model with demographic covariates
# =============================================================================


# =============================================================================
# SECTION 0: Libraries + Data Loading
# =============================================================================


library(tidyverse)
library(lme4)
library(lmerTest)
library(gtsummary)
library(gt)
library(broom.mixed)   # required for tbl_regression() on lmer objects

CommonPTID <- Reduce(intersect,list(All_Subjects_ADAS_17Feb2026$PTID,
                                    All_Subjects_APOERES_17Feb2026$PTID,
                                    All_Subjects_CDR_17Feb2026$PTID,
                                    All_Subjects_MMSE_17Feb2026$PTID,
                                    All_Subjects_PTDEMOG_17Feb2026$PTID,
                                    All_Subjects_UPENNBIOMK_ROCHE_ELECSYS_17Feb2026$PTID
))
length(unique(CommonPTID))

#filtered datasets

adas <- All_Subjects_ADAS_17Feb2026[All_Subjects_ADAS_17Feb2026$PTID 
                                            %in% CommonPTID,]

apoe <- All_Subjects_APOERES_17Feb2026[All_Subjects_APOERES_17Feb2026$PTID 
                                                  %in% CommonPTID,]

cdr  <- All_Subjects_CDR_17Feb2026[All_Subjects_CDR_17Feb2026$PTID 
                                          %in% CommonPTID,]

mmse <- All_Subjects_MMSE_17Feb2026[All_Subjects_MMSE_17Feb2026$PTID 
                                            %in% CommonPTID,]

demog <- All_Subjects_PTDEMOG_17Feb2026[All_Subjects_PTDEMOG_17Feb2026$PTID 
                                                  %in% CommonPTID,]

biomarkers <- All_Subjects_UPENNBIOMK_ROCHE_ELECSYS_17Feb2026[All_Subjects_UPENNBIOMK_ROCHE_ELECSYS_17Feb2026$PTID 
                                                                      %in% CommonPTID,]

# Set working directory — update to your local path
#setwd("./outputs/ADNI_EDA_data")

# Load all 6 ADNI data tables
#adas       <- read_csv("All_Subjects_ADAS_17Feb2026.csv",                      show_col_types = FALSE)
#apoe       <- read_csv("All_Subjects_APOERES_17Feb2026.csv",                   show_col_types = FALSE)
#cdr        <- read_csv("All_Subjects_CDR_17Feb2026.csv",                       show_col_types = FALSE)
#mmse       <- read_csv("All_Subjects_MMSE_17Feb2026.csv",                      show_col_types = FALSE)
#demog      <- read_csv("All_Subjects_PTDEMOG_17Feb2026.csv",                   show_col_types = FALSE)
#biomarkers <- read_csv("All_Subjects_UPENNBIOMK_ROCHE_ELECSYS_17Feb2026.csv",  show_col_types = FALSE)

# Recode MMSE sentinel value (-1 = missing/invalid)
mmse <- mmse %>%
  mutate(MMSCORE = if_else(MMSCORE == -1, NA_real_, MMSCORE))

# --- Shared lookup tables (used across multiple sections) ---

# APOE ε4 carrier status: one row per participant
apoe_status <- apoe %>%
  select(RID, GENOTYPE) %>%
  distinct() %>%
  mutate(apoe4_status = if_else(str_detect(GENOTYPE, "4"), "Carrier", "Non-carrier"))

# Demographic variables: one row per participant
demog_clean <- demog %>%
  select(RID, PTGENDER, PTEDUCAT, PTDOB) %>%
  distinct(RID, .keep_all = TRUE) %>%
  rename(SEX = PTGENDER, EDUCATION = PTEDUCAT, DOB = PTDOB) %>%
  mutate(SEX = factor(SEX, levels = c(1, 2), labels = c("Male", "Female")))


# =============================================================================
# SECTION 1: Cross-Sectional Analysis — OLS Models 1–4
#
#   Model 1: MMSE ~ Aβ42
#   Model 2: MMSE ~ Aβ42 + Tau
#   Model 3: MMSE ~ Aβ42 + Tau + APOE ε4
#   Model 4: MMSE ~ Aβ42 + Tau + APOE ε4 + Age + Sex + Education
#
#   Story: Each model adds a layer — biomarkers → genetic risk → demographics.
#          R² progression shows incremental improvement in explanatory power.
# =============================================================================

cat("=== SECTION 1: Cross-Sectional Analysis ===\n\n")

# --- 1a. Missing data report (pre-filter) ---
mmse_raw_n       <- nrow(mmse %>% select(RID, VISCODE2, MMSCORE))
mmse_missing_n   <- sum(is.na(mmse$MMSCORE))
biomarker_raw_n  <- nrow(biomarkers %>% select(RID, VISCODE2, ABETA42, TAU))
biomarker_ab_na  <- sum(is.na(biomarkers$ABETA42))
biomarker_tau_na <- sum(is.na(biomarkers$TAU))

cat("Pre-filter counts:\n")
cat("  MMSE table: total rows =", mmse_raw_n,
    "| missing MMSCORE =", mmse_missing_n, "\n")
cat("  Biomarker table: total rows =", biomarker_raw_n,
    "| missing ABETA42 =", biomarker_ab_na,
    "| missing TAU =", biomarker_tau_na, "\n\n")

# --- 1b. Build base cross-sectional dataset (complete-case) ---
analysis1 <- mmse %>%
  select(RID, VISCODE2, MMSCORE) %>%
  mutate(VISCODE2 = trimws(tolower(VISCODE2))) %>%
  inner_join(
    biomarkers %>%
      select(RID, VISCODE2, ABETA42, TAU) %>%
      mutate(VISCODE2 = trimws(tolower(VISCODE2))),
    by = c("RID", "VISCODE2")
  ) %>%
  drop_na(MMSCORE, ABETA42, TAU)

cat("Analysis 1 — after complete-case filter:\n")
cat("  Rows (visit-observations):", nrow(analysis1), "\n")
cat("  Unique participants:       ", n_distinct(analysis1$RID), "\n")
cat("  Pearson r (Aβ42 vs MMSE):",
    round(cor(analysis1$ABETA42, analysis1$MMSCORE, use = "complete.obs"), 3), "\n")
cat("  Proportion MMSE == 30 (ceiling):",
    round(mean(analysis1$MMSCORE == 30, na.rm = TRUE), 3), "\n\n")

# --- 1c. Figure 1: CSF Aβ42 vs MMSE (basic, unstratified) ---
plot_analysis1 <- ggplot(analysis1, aes(x = ABETA42, y = MMSCORE)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm", se = TRUE, fill = "steelblue", alpha = 0.2) +
  scale_y_continuous(limits = c(0, 30), breaks = seq(0, 30, 5)) +
  labs(
    title = "CSF Aβ42 vs MMSE",
    x     = "CSF Aβ42 (pg/mL)",
    y     = "MMSE Score"
  ) +
  theme_minimal()

print(plot_analysis1)
ggsave("./outputs/Fig1_ABeta42_vs_MMSE.png",
       plot = plot_analysis1, width = 8, height = 6, dpi = 300)

# --- 1d. OLS Models 1 & 2 ---
model1_abeta     <- lm(MMSCORE ~ ABETA42,       data = analysis1)
model2_abeta_tau <- lm(MMSCORE ~ ABETA42 + TAU, data = analysis1)

cat("--- Model 1: MMSE ~ Aβ42 ---\n")
print(summary(model1_abeta))

cat("\n--- Model 2: MMSE ~ Aβ42 + Tau ---\n")
print(summary(model2_abeta_tau))

cat("\nAIC comparison (Models 1 vs 2):\n")
print(AIC(model1_abeta, model2_abeta_tau))

cat("\nNested F-test (Model 1 vs Model 2):\n")
print(anova(model1_abeta, model2_abeta_tau))

# --- 1e. Add APOE ε4 status to cross-sectional dataset ---
analysis1_apoe <- analysis1 %>%
  inner_join(apoe_status %>% select(RID, apoe4_status), by = "RID") %>%
  mutate(apoe4_status = factor(apoe4_status, levels = c("Carrier", "Non-carrier")))

# --- 1f. OLS Model 3: Biomarkers + APOE ε4 ---
model1_apoe_refit <- lm(MMSCORE ~ ABETA42,                       data = analysis1_apoe)
model2_apoe_refit <- lm(MMSCORE ~ ABETA42 + TAU,                 data = analysis1_apoe)
model3_apoe       <- lm(MMSCORE ~ ABETA42 + TAU + apoe4_status,  data = analysis1_apoe)

cat("\n--- Model 3: MMSE ~ Aβ42 + Tau + APOE ε4 ---\n")
print(summary(model3_apoe))

cat("\nNested F-test (Model 2 vs Model 3):\n")
print(anova(model2_apoe_refit, model3_apoe))

# --- 1g. Formal interaction test: does the Aβ42–MMSE slope differ by APOE status? ---
# This tests whether APOE ε4 moderates the biomarker-cognition relationship
model_interaction <- lm(MMSCORE ~ ABETA42 * apoe4_status, data = analysis1_apoe)

cat("\n--- Interaction Model: MMSE ~ Aβ42 * APOE ε4 ---\n")
print(summary(model_interaction))
cat("\n  Key: p-value for ABETA42:apoe4_statusNon-carrier tests whether the\n")
cat("  Aβ42-MMSE slope is significantly different between carriers and non-carriers.\n\n")

# --- 1h. Figure 2: CSF Aβ42 vs MMSE stratified by APOE ε4 (y-axis capped at 30) ---
plot_apoe_scatter <- ggplot(
  analysis1_apoe,
  aes(x = ABETA42, y = MMSCORE, color = apoe4_status)
) +
  geom_point(alpha = 0.35, size = 0.9) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 1.1) +
  scale_y_continuous(limits = c(0, 30), breaks = seq(0, 30, 5)) +
  labs(
    title = "CSF Aβ42 vs MMSE by APOE ε4 Status",
    x     = "CSF Aβ42 (pg/mL)",
    y     = "MMSE Score",
    color = "APOE ε4 Status"
  ) +
  theme_minimal()

print(plot_apoe_scatter)
ggsave("./outputs/Fig2_ABeta42_vs_MMSE_by_APOE.png",
       plot = plot_apoe_scatter, width = 8, height = 6, dpi = 300)

# --- 1i. Add demographic covariates ---
# Note: ADNI MMSE table stores visit date as USERDATE (not EXAMDATE)
analysis1_demog <- analysis1_apoe %>%
  inner_join(demog_clean, by = "RID") %>%
  inner_join(
    mmse %>%
      select(RID, VISCODE2, EXAMDATE = USERDATE) %>%
      mutate(VISCODE2 = trimws(tolower(VISCODE2))),
    by = c("RID", "VISCODE2")
  ) %>%
  mutate(
    DOB      = as.Date(paste0("01/", DOB), format = "%d/%m/%Y"),
    EXAMDATE = as.Date(EXAMDATE),
    AGE      = as.numeric(difftime(EXAMDATE, DOB, units = "days")) / 365.25
  ) %>%
  drop_na(AGE, SEX, EDUCATION)

cat("Demographic-adjusted dataset: rows =", nrow(analysis1_demog),
    "| unique participants =", n_distinct(analysis1_demog$RID), "\n\n")

# --- 1j. OLS Model 4: Biomarkers + APOE ε4 + Demographics ---
model4_demog <- lm(
  MMSCORE ~ ABETA42 + TAU + apoe4_status + AGE + SEX + EDUCATION,
  data = analysis1_demog
)

cat("--- Model 4: MMSE ~ Aβ42 + Tau + APOE + Age + Sex + Education ---\n")
print(summary(model4_demog))

# Refit Models 1–3 on same demographic sample for valid AIC/R² comparison
m1_demog <- lm(MMSCORE ~ ABETA42,                              data = analysis1_demog)
m2_demog <- lm(MMSCORE ~ ABETA42 + TAU,                        data = analysis1_demog)
m3_demog <- lm(MMSCORE ~ ABETA42 + TAU + apoe4_status,         data = analysis1_demog)

cat("\nNested F-test (Model 3 vs Model 4):\n")
print(anova(m3_demog, model4_demog))

cat("\nAIC comparison — Models 1–4 (same sample):\n")
print(AIC(m1_demog, m2_demog, m3_demog, model4_demog))

cat("\nR-squared progression across Models 1–4:\n")
cat("  Model 1 R²:", round(summary(m1_demog)$r.squared,    3), "\n")
cat("  Model 2 R²:", round(summary(m2_demog)$r.squared,    3), "\n")
cat("  Model 3 R²:", round(summary(m3_demog)$r.squared,    3), "\n")
cat("  Model 4 R²:", round(summary(model4_demog)$r.squared, 3), "\n\n")
cat("  (R² progression quantifies incremental gain from each added predictor)\n")


# =============================================================================
# SECTION 2: Longitudinal EDA — APOE ε4 Cognitive Trajectories
#
#   Visualises MMSE trajectories over the first 60 months by APOE ε4 status.
#   Sets up the rationale for mixed-effects modelling in Section 3.
# =============================================================================

cat("\n=== SECTION 2: Longitudinal EDA ===\n\n")

# Build longitudinal dataset (0–60 months, valid MMSE only)
analysis2 <- mmse %>%
  select(RID, VISCODE2, MMSCORE) %>%
  inner_join(apoe_status, by = "RID") %>%
  mutate(
    VISCODE2     = trimws(tolower(VISCODE2)),
    months       = suppressWarnings(as.numeric(sub("^m", "", VISCODE2))),
    apoe4_status = factor(apoe4_status, levels = c("Carrier", "Non-carrier"))
  ) %>%
  drop_na(MMSCORE, months) %>%
  filter(months <= 60)

cat("Longitudinal dataset (0–60 months):\n")
cat("  Rows (visit-observations):", nrow(analysis2), "\n")
cat("  Unique participants:       ", n_distinct(analysis2$RID), "\n\n")

# Group-level MMSE descriptives
apoe_summary <- analysis2 %>%
  group_by(apoe4_status) %>%
  summarise(
    mean_mmse = round(mean(MMSCORE, na.rm = TRUE), 2),
    sd_mmse   = round(sd(MMSCORE,   na.rm = TRUE), 2),
    n         = n(),
    .groups   = "drop"
  )

cat("Group-level MMSE descriptives by APOE ε4 status:\n")
print(apoe_summary)

cat("\nWelch two-sample t-test (MMSE ~ APOE ε4 status):\n")
print(t.test(MMSCORE ~ apoe4_status, data = analysis2))

# --- Figure 3: Mean MMSE trajectories by APOE ε4 ---
plot_analysis2_mean <- ggplot(
  analysis2,
  aes(x = months, y = MMSCORE, color = apoe4_status)
) +
  stat_summary(fun = mean, geom = "line",  linewidth = 1.2) +
  stat_summary(fun = mean, geom = "point", size = 2) +
  scale_y_continuous(limits = c(0, 30), breaks = seq(0, 30, by = 5)) +
  labs(
    title = "Mean MMSE Trajectory by APOE ε4 Status (0–60 Months)",
    x     = "Months Since Baseline",
    y     = "Mean MMSE Score",
    color = "APOE ε4 Status"
  ) +
  theme_minimal()

print(plot_analysis2_mean)
ggsave("./outputs/Fig3_Mean_MMSE_Trajectory_by_APOE.png",
       plot = plot_analysis2_mean, width = 8, height = 6, dpi = 300)

# --- Figure 4: Spaghetti plot — individual trajectories by APOE ε4 ---
set.seed(7)
sample_rids      <- sample(unique(analysis2$RID), size = min(150, n_distinct(analysis2$RID)))
analysis2_sample <- analysis2 %>% filter(RID %in% sample_rids)

plot_analysis2_spaghetti <- ggplot(
  analysis2_sample,
  aes(x = months, y = MMSCORE, group = RID, color = apoe4_status)
) +
  geom_line(alpha = 0.35, linewidth = 0.5) +
  labs(
    title = "Individual MMSE Trajectories (n = 150 sample) by APOE ε4 Status",
    x     = "Months Since Baseline",
    y     = "MMSE Score",
    color = "APOE ε4 Status"
  ) +
  theme_minimal()

print(plot_analysis2_spaghetti)
ggsave("./outputs/Fig4_Individual_MMSE_Trajectories.png",
       plot = plot_analysis2_spaghetti, width = 9, height = 6, dpi = 300)

# --- Figure 5: MMSE distribution by APOE ε4 (boxplot) ---
plot_analysis2_box <- ggplot(
  analysis2,
  aes(x = apoe4_status, y = MMSCORE, fill = apoe4_status)
) +
  geom_boxplot(outlier.alpha = 0.6) +
  labs(
    title = "Distribution of MMSE Scores by APOE ε4 Status (0–60 Months)",
    x     = "APOE ε4 Status",
    y     = "MMSE Score"
  ) +
  theme_minimal() +
  theme(legend.position = "none")

print(plot_analysis2_box)
ggsave("./outputs/Fig5_MMSE_Boxplot_by_APOE.png",
       plot = plot_analysis2_box, width = 7, height = 6, dpi = 300)


# =============================================================================
# SECTION 3: Linear Mixed-Effects Models
#
#   Model A: Random intercept — baseline model for LRT comparison
#   Model B: Random slopes — each participant's decline rate varies freely
#   LRT:     Confirms random slopes model is significantly better fit
#   Model C: Random slopes + demographic covariates — FINAL model for tbl_F3
#
#   All models estimated with ML for LRT; REML for final parameter estimates.
# =============================================================================

cat("\n=== SECTION 3: Linear Mixed-Effects Models ===\n\n")

# --- Model A: Random intercept (REML for estimates) ---
mixed_model_ri <- lmer(
  MMSCORE ~ months + apoe4_status + months:apoe4_status + (1 | RID),
  data    = analysis2,
  REML    = TRUE,
  control = lmerControl(optimizer = "bobyqa")
)

cat("--- LMM Model A (Random Intercept): MMSE ~ months * APOE + (1|RID) ---\n")
print(summary(mixed_model_ri))

# --- Model B: Random slopes (REML for estimates) ---
mixed_model_rs <- lmer(
  MMSCORE ~ months + apoe4_status + months:apoe4_status + (1 + months | RID),
  data    = analysis2,
  REML    = TRUE,
  control = lmerControl(optimizer = "bobyqa")
)

cat("\n--- LMM Model B (Random Slopes): MMSE ~ months * APOE + (1 + months|RID) ---\n")
print(summary(mixed_model_rs))
cat("\nFixed effects:\n")
print(fixef(mixed_model_rs))

# --- Likelihood Ratio Test: Random Intercept vs Random Slopes ---
# Models must be refit with ML (not REML) for valid LRT comparison
cat("\n--- Likelihood Ratio Test: Model A vs Model B ---\n")
mixed_ri_ml <- lmer(
  MMSCORE ~ months + apoe4_status + months:apoe4_status + (1 | RID),
  data    = analysis2,
  REML    = FALSE,
  control = lmerControl(optimizer = "bobyqa")
)
mixed_rs_ml <- lmer(
  MMSCORE ~ months + apoe4_status + months:apoe4_status + (1 + months | RID),
  data    = analysis2,
  REML    = FALSE,
  control = lmerControl(optimizer = "bobyqa")
)

print(anova(mixed_ri_ml, mixed_rs_ml))
cat("\n  (Significant LRT p-value confirms random slopes model is the better fit)\n")

# --- Figure 6: Predicted MMSE trajectories — random slopes model ---
set.seed(7)
sample_rids2   <- sample(unique(analysis2$RID), size = min(80, n_distinct(analysis2$RID)))
analysis2_sub  <- analysis2 %>% filter(RID %in% sample_rids2)
analysis2_pred <- analysis2_sub %>%
  mutate(predicted = predict(mixed_model_rs, newdata = analysis2_sub))

plot_random_slopes <- ggplot(
  analysis2_pred,
  aes(x = months, y = predicted, group = RID, color = apoe4_status)
) +
  geom_line(alpha = 0.4, linewidth = 0.5) +
  scale_y_continuous(limits = c(0, 30), breaks = seq(0, 30, 5)) +
  labs(
    title = "Predicted MMSE Trajectories — Random Slopes Model (n = 80 sample)",
    x     = "Months Since Baseline",
    y     = "Predicted MMSE Score",
    color = "APOE ε4 Status"
  ) +
  theme_minimal()

print(plot_random_slopes)
ggsave("./outputs/Fig6_Predicted_Trajectories_Random_Slopes.png",
       plot = plot_random_slopes, width = 9, height = 6, dpi = 300)

# --- Model C: Random slopes + demographic covariates (full model) ---
# Build longitudinal dataset with age, sex, education joined in
analysis2_demog <- analysis2 %>%
  inner_join(demog_clean, by = "RID") %>%
  inner_join(
    mmse %>%
      select(RID, VISCODE2, EXAMDATE = USERDATE) %>%   # USERDATE is the correct column
      mutate(VISCODE2 = trimws(tolower(VISCODE2))),
    by = c("RID", "VISCODE2")
  ) %>%
  mutate(
    DOB      = as.Date(paste0("01/", DOB), format = "%d/%m/%Y"),
    EXAMDATE = as.Date(EXAMDATE),
    AGE      = as.numeric(difftime(EXAMDATE, DOB, units = "days")) / 365.25
  ) %>%
  drop_na(AGE, SEX, EDUCATION)

cat("\nDemographic-adjusted longitudinal dataset: rows =", nrow(analysis2_demog),
    "| unique participants =", n_distinct(analysis2_demog$RID), "\n\n")

mixed_model_full <- lmer(
  MMSCORE ~ months + apoe4_status + months:apoe4_status +
    AGE + SEX + EDUCATION + (1 + months | RID),
  data    = analysis2_demog,
  REML    = TRUE,
  control = lmerControl(optimizer = "bobyqa")
)

cat("--- LMM Model C (Full Model): Random Slopes + Age + Sex + Education ---\n")
print(summary(mixed_model_full))
cat("\nFixed effects:\n")
print(fixef(mixed_model_full))


# =============================================================================
# SECTION 4: Secondary Outcome — ADAS-Cog Validation
#
#   ADAS-Cog Total Score ~ Aβ42 + Tau
#   Higher ADAS-Cog = worse cognition (inverse direction vs MMSE)
#   Purpose 1: Validates the biomarker-cognition associations across a second measure
#   Purpose 2: Score distribution comparison illustrates MMSE ceiling effect (Fig 7)
# =============================================================================

cat("\n=== SECTION 4: ADAS-Cog Secondary Outcome ===\n\n")

adas_clean <- adas %>%
  select(RID, VISCODE2, ADAS_TOTAL = TOTSCORE) %>%
  mutate(VISCODE2 = trimws(tolower(VISCODE2))) %>%
  filter(!is.na(ADAS_TOTAL), ADAS_TOTAL >= 0)

analysis_adas <- biomarkers %>%
  select(RID, VISCODE2, ABETA42, TAU) %>%
  mutate(VISCODE2 = trimws(tolower(VISCODE2))) %>%
  inner_join(adas_clean, by = c("RID", "VISCODE2")) %>%
  drop_na(ADAS_TOTAL, ABETA42, TAU)

cat("ADAS-Cog dataset: rows =", nrow(analysis_adas),
    "| unique participants =", n_distinct(analysis_adas$RID), "\n")
cat("\nADAS-Cog score summary (higher score = worse cognition):\n")
print(summary(analysis_adas$ADAS_TOTAL))
cat("\nCeiling check — proportion MMSE == 30 in cross-sectional sample:",
    round(mean(analysis1$MMSCORE == 30, na.rm = TRUE), 3), "\n\n")

adas_model <- lm(ADAS_TOTAL ~ ABETA42 + TAU, data = analysis_adas)
cat("--- ADAS-Cog Model: ADAS ~ Aβ42 + Tau ---\n")
print(summary(adas_model))

# --- Figure 7: MMSE vs ADAS-Cog score distribution (ceiling effect) ---
mmse_dist <- analysis1 %>%
  select(MMSCORE) %>%
  mutate(Scale = "MMSE",     score = MMSCORE / 30)

adas_dist <- analysis_adas %>%
  select(ADAS_TOTAL) %>%
  mutate(Scale = "ADAS-Cog", score = ADAS_TOTAL / 70)

ceiling_compare <- bind_rows(mmse_dist, adas_dist)

plot_ceiling <- ggplot(ceiling_compare, aes(x = score, fill = Scale)) +
  geom_histogram(bins = 30, alpha = 0.65, position = "identity") +
  facet_wrap(~ Scale, scales = "free_x") +
  labs(
    title = "Score Distribution: MMSE vs ADAS-Cog (Normalized 0–1)",
    x     = "Normalized Score",
    y     = "Count"
  ) +
  theme_minimal() +
  theme(legend.position = "none")

print(plot_ceiling)
ggsave("./outputs/Fig7_MMSE_vs_ADAS_Distribution.png",
       plot = plot_ceiling, width = 9, height = 5, dpi = 300)


# =============================================================================
# SECTION 5: Model Validation — Train/Test Split (OLS)
#
#   Addresses reviewer feedback: "use test data to make predictions and assess
#   model accuracy before drawing conclusions."
#
#   Approach: 80/20 random split of cross-sectional dataset.
#     - Model 4 is fit on the 80% training set only
#     - Predictions are made on the held-out 20% test set
#     - RMSE and R² on the test set measure out-of-sample performance
#     - If train R² ≈ test R², the model generalises (not overfitting)
# =============================================================================

cat("\n=== SECTION 5: Train/Test Validation ===\n\n")

set.seed(42)
n_obs      <- nrow(analysis1_demog)
train_idx  <- sample(seq_len(n_obs), size = floor(0.8 * n_obs))
train_data <- analysis1_demog[train_idx, ]
test_data  <- analysis1_demog[-train_idx, ]

cat("Train set: n =", nrow(train_data), "\n")
cat("Test set:  n =", nrow(test_data),  "\n\n")

# Fit Model 4 on training data only
model4_train <- lm(
  MMSCORE ~ ABETA42 + TAU + apoe4_status + AGE + SEX + EDUCATION,
  data = train_data
)

cat("--- Model 4 trained on 80% of data ---\n")
print(summary(model4_train))

# Predict on held-out 20% test set
test_preds <- predict(model4_train, newdata = test_data)

# Compute out-of-sample RMSE and R²
rmse_test <- sqrt(mean((test_data$MMSCORE - test_preds)^2, na.rm = TRUE))
ss_res    <- sum((test_data$MMSCORE - test_preds)^2, na.rm = TRUE)
ss_tot    <- sum((test_data$MMSCORE - mean(test_data$MMSCORE, na.rm = TRUE))^2, na.rm = TRUE)
r2_test   <- 1 - (ss_res / ss_tot)

cat("\n--- Out-of-Sample Performance (20% Test Set) ---\n")
cat("  RMSE (test):              ", round(rmse_test, 3), "\n")
cat("  R²   (test):              ", round(r2_test,   3), "\n")
cat("  R²   (train, comparison): ", round(summary(model4_train)$r.squared, 3), "\n\n")
cat("  Interpretation: Similar train/test R² indicates the model generalises\n")
cat("  to unseen data without overfitting.\n")


# =============================================================================
# SECTION 6: Publication Tables (gtsummary)
#
#   tbl_F1: Participant characteristics by APOE ε4 status
#   tbl_F2: Cross-sectional OLS model comparison (Models 1–4)
#   tbl_F3: Full longitudinal LMM (random slopes + demographic covariates)
# =============================================================================

cat("\n=== SECTION 6: Publication Tables ===\n\n")
if (!requireNamespace("webshot2", quietly = TRUE)) install.packages("webshot2")
library(webshot2)
dir.create("./outputs", showWarnings = FALSE)
# --- tbl_F1: Participant Characteristics by APOE ε4 Status ---
cat("--- Table 1: Participant Characteristics ---\n")

char_table_data <- analysis1_demog %>%
  select(MMSCORE, ABETA42, TAU, AGE, SEX, EDUCATION, apoe4_status) %>%
  rename(
    `MMSE Score`        = MMSCORE,
    `CSF Aβ42`       = ABETA42,
    `CSF Tau (pg/mL)`   = TAU,
    `Age (years)`       = AGE,
    `Sex`               = SEX,
    `Education (years)` = EDUCATION,
    `APOE ε4 Status`    = apoe4_status
  )

tbl_F1 <- tbl_summary(
  char_table_data,
  by        = `APOE ε4 Status`,
  statistic = list(
    all_continuous()  ~ "{mean} ({sd})",
    all_categorical() ~ "{n} ({p}%)"
  ),
  digits = all_continuous() ~ 2
) %>%
  add_p() %>%
  bold_p(t = 0.05) %>%
  bold_labels() %>%
  modify_caption("**Table 1. Participant Characteristics by APOE ε4 Status**")

print(tbl_F1)
tbl_F1 %>%
  as_gt() %>%
  gt::gtsave("./outputs/Table1_characteristics.png", vwidth = 1400, vheight = 900)

cat("--- Table 1 saved as PNG ---\n")

# --- tbl_F2: Cross-Sectional OLS Model Comparison (Models 1–4) ---
cat("\n--- Table 3: Cross-Sectional OLS Models 1–4 ---\n")

# All four models refitted on same demographic-filtered sample for comparability
tbl_m1 <- tbl_regression(m1_demog,
  label        = list(ABETA42 ~ "CSF Aβ42"),
  estimate_fun = ~style_sigfig(., digits = 4)) %>%
  bold_p(t = 0.05)

tbl_m2 <- tbl_regression(m2_demog,
  label        = list(ABETA42 ~ "CSF Aβ42",
                      TAU     ~ "CSF Tau"),
  estimate_fun = ~style_sigfig(., digits = 4)) %>%
  bold_p(t = 0.05)

tbl_m3 <- tbl_regression(m3_demog,
  label        = list(ABETA42      ~ "CSF Aβ42",
                      TAU          ~ "CSF Tau",
                      apoe4_status ~ "APOE ε4 Status"),
  estimate_fun = ~style_sigfig(., digits = 4)) %>%
  bold_p(t = 0.05)

tbl_m4 <- tbl_regression(model4_demog,
  label        = list(ABETA42      ~ "CSF Aβ42",
                      TAU          ~ "CSF Tau",
                      apoe4_status ~ "APOE ε4 Status",
                      AGE          ~ "Age (years)",
                      SEX          ~ "Sex",
                      EDUCATION    ~ "Education (years)"),
  estimate_fun = ~style_sigfig(., digits = 4)) %>%
  bold_p(t = 0.05)

tbl_F2 <- tbl_merge(tbls= list(tbl_m1, tbl_m2, tbl_m3, tbl_m4), 
                    tab_spanner = c("**Model 1**", "**Model 2**", "**Model 3**", "**Model 4**")
                    ) %>% modify_caption("**Table 2. Cross-Sectional OLS Regression Models (Outcome: MMSE Score)**")

print(tbl_F2)
tbl_F2 %>%
  as_gt() %>%
  gt::gtsave("./outputs/Table2_crosssectional_models.png", vwidth = 1600, vheight = 1000)
cat("--- Table 2 saved as PNG ---\n")

# --- tbl_F3: Longitudinal LMM — Full Model (Random Slopes + Demographics) ---
cat("\n--- Table 4: Longitudinal LMM Full Model ---\n")

tbl_F3 <- tbl_regression(
  mixed_model_full,
  label = list(
    months                ~ "Time (months)",
    apoe4_status          ~ "APOE ε4 Status",
    `months:apoe4_status` ~ "Time x APOE ε4 (Non-carrier vs Carrier)",
    AGE                   ~ "Age (years)",
    SEX                   ~ "Sex",
    EDUCATION             ~ "Education (years)"
  )
) %>%
  bold_p(t = 0.05) %>%
  bold_labels() %>%
  modify_caption("**Table 3. Linear Mixed-Effects Model — Random Slopes + Demographics (Outcome: MMSE Score)**")

print(tbl_F3)
tbl_F3 %>%
  as_gt() %>%
  gt::gtsave("./outputs/Table3_longitudinal_LMM.png", vwidth = 1400, vheight = 900)
cat("--- Table 3 saved as PNG ---\n")

cat("\n=== Script complete. All manuscript figures and tables generated. ===\n")
