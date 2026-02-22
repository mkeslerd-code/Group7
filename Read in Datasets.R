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
AgeAndSexDist <- read_csv("ADNI_Participant_Age_Distribution_02-21-2026.csv")
View(AgeAndSexDist )




# Code for Age and Sex distribution 
AgeAndSex_long <- AgeAndSexDist %>% pivot_longer(cols = c(Female, Male), names_to = "Gender", values_to = "Count")

ggplot(AgeAndSex_long, aes(x= `Age Group`, y= Count, fill = Gender))+
  geom_col(position = "dodge")+ 
  labs(
    title = "Distribution of Age and Gender of ADNI Participants",
    x= "Age Group",
    y= "Number of Participants"
  ) +
  theme_minimal(base_size = 13)+
  theme(plot.title = element_text(hjust = 0.5))

#ggplot(AgeAndSex_long, aes(x= `Age Group`, y= Count, fill = Gender))+
#  geom_col(position = "fill")+ 
#  scale_y_continuous(labels = scales::percent)+
#  labs(
#    title = "Proportional Distribution of Age and Gender of ADNI Participants",
#    x= "Age Group",
#    y= "Percentage"
#  ) +
#  theme_minimal()

Totals <- AgeAndSex_long %>%
  group_by(Gender)%>%
  summarise(Total=sum(Count,na.rm = TRUE))

ggplot(Totals, aes(x=Gender, y= Total, fill=Gender))+
  geom_col()+
  labs(
    title = "Total Participant Counts by Gender",
    x = "",
    y = "Total Number of Participants"
  )+
  theme_minimal(base_size = 13)+
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position = "none")


Totals2 <- Totals %>%
  mutate(Percent = Total / sum(Total)*100)

ggplot(Totals2, aes(x=Gender, y= Total, fill=Gender))+
  geom_col()+
  geom_text(
    aes(label = paste0(round(Percent, 1),"%")),
    vjust = -0.5,
    size = 5
  )+
  scale_fill_manual(values = c(
    "Female" = "coral", 
    "Male" = "maroon"
  ))+
  labs(
    title = "Total AD Participant Counts by Gender",
    x = "",
    y = "Total Number of Participants"
  )+
  theme_minimal(base_size = 13)+
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position = "none")




