library(readr)
All_Subjects_ADAS_17Feb2026 <- read_csv("All_Subjects_ADAS_17Feb2026.csv")
#View(All_Subjects_ADAS_17Feb2026)
All_Subjects_APOERES_17Feb2026 <- read_csv("All_Subjects_APOERES_17Feb2026.csv")
#View(All_Subjects_APOERES_17Feb2026)
All_Subjects_CDR_17Feb2026 <- read_csv("All_Subjects_CDR_17Feb2026.csv")
#View(All_Subjects_CDR_17Feb2026)
All_Subjects_MMSE_17Feb2026 <- read_csv("All_Subjects_MMSE_17Feb2026.csv")
#View(All_Subjects_MMSE_17Feb2026)
All_Subjects_PTDEMOG_17Feb2026 <- read_csv("All_Subjects_PTDEMOG_17Feb2026.csv")
#View(All_Subjects_PTDEMOG_17Feb2026)
All_Subjects_UPENNBIOMK_ROCHE_ELECSYS_17Feb2026 <- read_csv("All_Subjects_UPENNBIOMK_ROCHE_ELECSYS_17Feb2026.csv")
#View(All_Subjects_UPENNBIOMK_ROCHE_ELECSYS_17Feb2026)
DATADIC_19Feb2026 <- read_csv("DATADIC_19Feb2026.csv")
#View(DATADIC_19Feb2026)
#AgeAndSexDist <- read_csv("PTDEMOG_22Feb2026.csv")
#View(AgeAndSexDist)


# count of participants in Demographics Data 
#CountDemographicsParticipants  <- length(unique(All_Subjects_PTDEMOG_17Feb2026$PTID))

#count of participants in ADAS-Cognitive Behavior data 
#CountADASCognativeBehaviorParticipants  <- length(unique(All_Subjects_ADAS_17Feb2026$PTID))

#count of participants in UPENN Biomarker data 
#CountUPENNBiomarkerParticipants <- length(unique(All_Subjects_UPENNBIOMK_ROCHE_ELECSYS_17Feb2026$PTID))

#count of participants in ApoE Genotyping - Results
#CountApoEGenotypingParticipants <- length(unique(All_Subjects_APOERES_17Feb2026$PTID))

#count of participants in Clinical Dementia Rating
#CountClinicalDementialParticipants <- length(unique(All_Subjects_CDR_17Feb2026$PTID))

#count of participants in Mini Mental State Exam
#CountMMSEParticipants <- length(unique(All_Subjects_MMSE_17Feb2026$PTID))


# find common PTID throughout all datasets 
CommonPTID <- Reduce(intersect,list(All_Subjects_ADAS_17Feb2026$PTID,
                                    All_Subjects_APOERES_17Feb2026$PTID,
                                    All_Subjects_CDR_17Feb2026$PTID,
                                    All_Subjects_MMSE_17Feb2026$PTID,
                                    All_Subjects_PTDEMOG_17Feb2026$PTID,
                                    All_Subjects_UPENNBIOMK_ROCHE_ELECSYS_17Feb2026$PTID
                                    ))
length(unique(CommonPTID))

#filtered datasets

FilteredADAS <- All_Subjects_ADAS_17Feb2026[All_Subjects_ADAS_17Feb2026$PTID 
                                            %in% CommonPTID,]

FilteredAPOERES <- All_Subjects_APOERES_17Feb2026[All_Subjects_APOERES_17Feb2026$PTID 
                                            %in% CommonPTID,]

FilteredCDR <- All_Subjects_CDR_17Feb2026[All_Subjects_CDR_17Feb2026$PTID 
                                            %in% CommonPTID,]

FilteredMMSE <- All_Subjects_MMSE_17Feb2026[All_Subjects_MMSE_17Feb2026$PTID 
                                            %in% CommonPTID,]

FilteredPTDEMOG <- All_Subjects_PTDEMOG_17Feb2026[All_Subjects_PTDEMOG_17Feb2026$PTID 
                                            %in% CommonPTID,]

FilteredUPENNBIOMK <- All_Subjects_UPENNBIOMK_ROCHE_ELECSYS_17Feb2026[All_Subjects_UPENNBIOMK_ROCHE_ELECSYS_17Feb2026$PTID 
                                            %in% CommonPTID,]


#To check for correct cleaning
# count of participants in ADAS-Cognitive Behavior data 
#CountFilteredADAS  <- length(unique(FilteredADAS$PTID))

#count of participants in ApoE Genotyping - Results
#CountFilteredAPOERES  <- length(unique(FilteredAPOERES$PTID))

#count of participants in Clinical Dementia Rating
#CountFilteredCDR <- length(unique(FilteredCDR$PTID))

#count of participants in Mini Mental State Exam
#CountFilteredMMSE <- length(unique(FilteredMMSE$PTID))

#count of participants in Demographics Data 
#CountFilteredPTDEMOG <- length(unique(FilteredPTDEMOG$PTID))

#count of participants in UPENN Biomarker data 
#CountFilteredUPENNBIOMK <- length(unique(FilteredUPENNBIOMK$PTID))



## must run ADNI_EDA_Manuscript first 
library(dplyr)
library(ggplot2)
demog$PTRACCAT_label <- dplyr::case_when(
      demog$PTRACCAT == 1 ~ "American Indian or Alaska Native",
      demog$PTRACCAT == 2 ~ "Asian",
      demog$PTRACCAT == 3 ~ "Native Hawaiian or Other Pacific Islander",
      demog$PTRACCAT == 4 ~ "Black or African American",
      demog$PTRACCAT == 5 ~ "White",
      demog$PTRACCAT == 7 ~ "Unknown",
      demog$PTRACCAT == 8 ~ "Native Hawaiian",
      demog$PTRACCAT == 9 ~ "Other Pacific Islander",
      TRUE ~ "Missing"
    )
head(demog$PTRACCAT_label)
table(demog$PTRACCAT_label)

demographics_summary <- demog %>%
  count(PTRACCAT_label)%>%
  mutate(percent = round(n/sum(n)* 100, 2))


ggplot(demographics_summary, aes(x = reorder(PTRACCAT_label, -n), y = n)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(x = "Race/Ethnicity", y = "Count") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


