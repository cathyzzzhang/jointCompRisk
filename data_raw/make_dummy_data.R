# file: data-raw/make_dummy_data.R

# 1) Load usethis for saving data
# install.packages("usethis") # if not already installed
set.seed(123) # For reproducibility

# ------------------------------------------------------------------
# 2) Create a dummy actt.ch dataset
#    (based on the real structure you've shown)
# ------------------------------------------------------------------
N <- 15  # small example number of patients
dummy_ch <- data.frame(
  USUBJID  = paste0("COV.DUMMY", 700 + seq_len(N)),
  TRTP     = sample(c("Remdesivir", "Placebo"), N, replace=TRUE),
  SEX      = sample(c("M", "F"), N, replace=TRUE),
  RACE     = sample(c("WHITE","BLACK OR AFRICAN AMERICAN","ASIAN"), N, replace=TRUE),
  ETHNIC   = sample(c("HISPANIC OR LATINO","NOT HISPANIC OR LATINO"), N, replace=TRUE),
  REGION   = sample(c("North America","Asia","Europe"), N, replace=TRUE),
  STRATUM  = sample(c("Severe Disease","Mild-Moderate Disease"), N, replace=TRUE),
  ORDSCRG  = sample(c("Baseline Clinical Status Score 4",
                      "Baseline Clinical Status Score 5",
                      "Baseline Clinical Status Score 6",
                      "Baseline Clinical Status Score 7"),
                    N, replace=TRUE),
  HYPFL    = sample(c("N","Y"), N, replace=TRUE),
  CADFL    = sample(c("N","Y"), N, replace=TRUE),
  CHFFL    = sample(c("N","Y"), N, replace=TRUE),
  CRDFL    = sample(c("N","Y"), N, replace=TRUE),
  CORFL    = sample(c("N","Y"), N, replace=TRUE),
  CLDFL    = sample(c("N","Y"), N, replace=TRUE),
  CKDFL    = sample(c("N","Y"), N, replace=TRUE),
  DIAB1FL  = sample(c("N","Y"), N, replace=TRUE),
  DIAB2FL  = sample(c("N","Y"), N, replace=TRUE),
  OBESIFL  = sample(c("N","Y"), N, replace=TRUE),
  CANCERFL = sample(c("N","Y"), N, replace=TRUE),
  IMMDFL   = sample(c("N","Y"), N, replace=TRUE),
  ASTHMAFL = sample(c("N","Y"), N, replace=TRUE),
  COMORB1  = sample(c("No Comorbidities","Any Comorbidities"), N, replace=TRUE),
  COMORB2  = sample(c("No Comorbidities","1 Comorbidity","2 or more Comorbidities"), N, replace=TRUE),
  AGE      = sample(30:90, N, replace=TRUE),
  BMI      = round(runif(N, 18, 45), 1),
  BDURSYMP = sample(1:15, N, replace=TRUE),
  ORDSCR15 = sample(4:7, N, replace=TRUE),   # baseline ordinal score
  TTRECOV  = sample(1:28, N, replace=TRUE),
  RECCNSR  = sample(c(0,1), N, replace=TRUE), # 0=event, 1=censored for recovery
  TTDEATH  = sample(1:28, N, replace=TRUE),
  DTHCNSR  = sample(c(0,1), N, replace=TRUE), # 0=event, 1=censored for death
  cn       = seq_len(N),                      # random ID column
  trt      = sample(c(0,1), N, replace=TRUE), # 0=control, 1=treatment
  sexn     = sample(c(0,1), N, replace=TRUE),
  regionn  = sample(1:3, N, replace=TRUE),
  strat    = sample(c(0,1), N, replace=TRUE),
  ordscr_bs= sample(4:7, N, replace=TRUE),
  com_bs1  = sample(c(0,1), N, replace=TRUE),
  com_bs2  = sample(c(1,2,3), N, replace=TRUE)
)

# ------------------------------------------------------------------
# 3) Create a dummy actt.long dataset
# ------------------------------------------------------------------
# We'll pretend each patient has some repeated measurements over time
M <- 50  # total long-format rows
dummy_long <- data.frame(
  USUBJID    = sample(dummy_actt_ch$USUBJID, M, replace=TRUE),
  ARM        = sample(c("Placebo", "Remdesivir"), M, replace=TRUE),
  ACTARM     = sample(c("Placebo", "Remdesivir"), M, replace=TRUE),
  BCSOSN     = sample(4:7, M, replace=TRUE),
  D29DTHE0   = sample(c(NA,5,6,7), M, replace=TRUE),
  D29DTHE1   = sample(c(NA,26,28), M, replace=TRUE),
  TTRECOV0   = sample(1:28, M, replace=TRUE),
  TTRECOV1   = sample(c(1,5,10,20,28), M, replace=TRUE),
  ADYC       = sample(1:28, M, replace=TRUE),
  ORDSCOR    = sample(4:7, M, replace=TRUE),
  OR15SCOR   = sample(1:7, M, replace=TRUE),
  AGEC       = sample(30:90, M, replace=TRUE),
  SEX        = sample(c("M","F"), M, replace=TRUE),
  BDURSYMP   = sample(1:15, M, replace=TRUE),
  COMORB2    = sample(c("No Comorbidities","1 Comorbidity","2 or more Comorbidities"), M, replace=TRUE)
)

# ------------------------------------------------------------------
# 4) Save them in the data/ folder as .rda files
# ------------------------------------------------------------------
usethis::use_data(dummy_ch, dummy_long, overwrite = TRUE)
