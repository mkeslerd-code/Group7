# Install packages (only runs if missing)
if (!requireNamespace("tidyverse", quietly = TRUE)) install.packages("tidyverse")
if (!requireNamespace("broom", quietly = TRUE)) install.packages("broom")
if (!requireNamespace("car", quietly = TRUE)) install.packages("car")

# Load library packages
library(tidyverse)
library(broom)
library(car)

# Locate Working Folder
file.choose()
data_folder <- "C:/Users/natal/OneDrive/Documents/ASU Course/ADNI_DatasetWorkingFolder"

# List files in folder + read files including column names
list.files(data_folder)

mmse_df <- readr::read_csv(file.path(data_folder,
                                     "All_Subjects_MMSE_17Feb2026.csv"))

bio_df  <- readr::read_csv(file.path(data_folder,
                                     "All_Subjects_UPENNBIOMK_ROCHE_ELECSYS_17Feb2026.csv"))

names(mmse_df)

names(bio_df)

# Create mini files for easy working
mmse_small <- mmse_df %>%
  select(RID, VISCODE2, MMSCORE)

bio_small <- bio_df %>%
  select(RID, VISCODE2, ABETA42, TAU)

# Merge the minis!
df_merged <- inner_join(mmse_small, bio_small, by = c("RID", "VISCODE2"))
nrow(df_merged)
head(df_merged)

# Clean up the merged data
df_clean <- df_merged %>%
  filter(!is.na(MMSCORE),
         !is.na(ABETA42),
         !is.na(TAU))

nrow(df_clean)

# Amyloid-Only Model
model1 <- lm(MMSCORE ~ ABETA42, data = df_clean)
summary(model1)

# Combined Model
model2 <- lm(MMSCORE ~ ABETA42 + TAU, data = df_clean)
summary(model2)

# Model 1 and Model 2 Comparison
anova(model1, model2)
