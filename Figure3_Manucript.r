library(dplyr)
library(ggplot2)

# Randomly select participants
set.seed(1)

sample_ids <- ALL_Subjects_MMSE %>%
  distinct(PTID) %>%
  slice_sample(n = 120)

# Plot trajectories
ALL_Subjects_MMSE %>%
  filter(PTID %in% sample_ids$PTID) %>%
  ggplot(aes(x = months_since_bl,
             y = MMSCORE,
             group = PTID)) +
  
  geom_line(alpha = 0.35, linewidth = 0.4, color = "darkblue") +
  
  labs(
    title = "Individual MMSE Trajectories (Sample of ADNI Participants)",
    x = "Months Since Baseline",
    y = "MMSE Score"
  ) +
  
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )