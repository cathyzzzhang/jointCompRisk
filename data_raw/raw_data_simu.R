#############
# Load the data
set.seed(123)
actt_ch <- read.csv("actt_ch.csv")

# Define the number of unique IDs to sample
sample_size <- 30
sampled_ids <- sample(unique(actt_ch$USUBJID), size = sample_size, replace = FALSE)

#sampled_data <- actt_ch%>%
#  filter(USUBJID %in% sampled_ids)%>%
#select(USUBJID,TTRECOV,TTDEATH,RECCNSR,DTHCNSR,ordscr_bs,trt)

sampled_data <- actt_ch %>%
  filter(USUBJID %in% sampled_ids) %>%
  rename(
    ID = USUBJID,  # Rename unique subject identifier
    TimeToRecovery = TTRECOV,  # Match recovery time column
    TimetoDeath = TTDEATH,  # Match death time column
    RecoveryCensoringIndicator = RECCNSR,  # Match recovery censoring column
    DeathCensoringIndicator = DTHCNSR,  # Match death censoring column
    BaselineScore = ordscr_bs,  # Match baseline score column
    Treatment = trt  # Match treatment column
  )%>%
  select(ID, TimeToRecovery, TimetoDeath, RecoveryCensoringIndicator,
         DeathCensoringIndicator, BaselineScore, Treatment)


print(sampled_data)
write.csv(sampled_data, "data_raw/DUMMY_CH.csv", row.names = FALSE)

#####
# Load actt_long dataset
actt_long <- read.csv("actt_long.csv")

# Extract all rows where USUBJID matches sampled_data$ID
longitudinal <- actt_long %>%
  filter(USUBJID %in% sampled_data$ID)

# Print the extracted data
print(longitudinal)

# Save to CSV if needed
write.csv(longitudinal, "data_raw/DUMMY_LONG.csv", row.names = FALSE)



dummy_long<- read.csv("vignettes/dummy_long.csv")
head(dummy_long)
dummy_long <- dummy_long%>%
  rename(PersonID = USUBJID,
         OrdinalScore = ORDSCOR,
         RelativeDay = ADYC) %>%
  select(PersonID, OrdinalScore, RelativeDay)
write.csv(dummy_long, "data_raw/DUMMY_LONG.csv", row.names = FALSE)
