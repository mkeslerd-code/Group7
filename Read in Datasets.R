#this code will read in all datasets

library(readr)
All_Subjects_ADAS_17Feb2026 <- read_csv("All_Subjects_ADAS_17Feb2026.csv")
View(All_Subjects_ADAS_17Feb2026)
library(readr)
All_Subjects_APOERES_17Feb2026 <- read_csv("All_Subjects_APOERES_17Feb2026.csv")
View(All_Subjects_APOERES_17Feb2026)
library(readr)
All_Subjects_CDR_17Feb2026 <- read_csv("All_Subjects_CDR_17Feb2026.csv")
View(All_Subjects_CDR_17Feb2026)
library(readr)
All_Subjects_MMSE_17Feb2026 <- read_csv("All_Subjects_MMSE_17Feb2026.csv")
View(All_Subjects_MMSE_17Feb2026)
library(readr)
All_Subjects_PTDEMOG_17Feb2026 <- read_csv("All_Subjects_PTDEMOG_17Feb2026.csv")
View(All_Subjects_PTDEMOG_17Feb2026)
library(readr)
All_Subjects_UPENNBIOMK_ROCHE_ELECSYS_17Feb2026 <- read_csv("All_Subjects_UPENNBIOMK_ROCHE_ELECSYS_17Feb2026.csv")
View(All_Subjects_UPENNBIOMK_ROCHE_ELECSYS_17Feb2026)
library(readr)
DATADIC_19Feb2026 <- read_csv("DATADIC_19Feb2026.csv")
View(DATADIC_19Feb2026)
library(readr)
AgeAndSexDist <- read_csv("PTDEMOG_22Feb2026.csv")
View(AgeAndSexDist)




# Code for Sex distribution 
AgeAndSex_long <- AgeAndSexDist %>% pivot_longer(cols = c(Female, Male), names_to = "Gender", values_to = "Count")

library(dplyr)
library(ggplot2)
library(scales)

gender_summary <- AgeAndSexDist %>%
  mutate(
    Gender = case_when(
      PTGENDER == 1 ~ "Male",
      PTGENDER == 2 ~ "Female",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(Gender)) %>%              
  count(Gender, name = "n") %>%
  mutate(
    pct = n / sum(n),
    pct_label = percent(pct, accuracy = 0.1) 
  )

gender_summary


ggplot(gender_summary, aes(x = Gender, y = n, fill = Gender)) +
  geom_col(width = 0.6) +
  geom_text(aes(label = pct_label),
            vjust = -0.4, size = 5) +
  scale_fill_manual(values = c("Male" = "coral", "Female" = "maroon"))+
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(
    title = "Gender Distribution",
    x = "Gender",
    y = "Count"
  ) +
  theme_minimal(base_size = 14)