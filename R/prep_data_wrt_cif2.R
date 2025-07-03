#' @title Prepare Data for Weighted CIF (Mirroring Original Part II, Dashboard-Style Merge)
#' @description This function replicates the original Part II logic:
#'   1) Exclude patients with zero survival or missing baseline, handle discharge-to-die.
#'   2) Merge main data and long data, create 'data.l'.
#'   3) Compute weighted times for death vs. discharge using a 'count' array and a loop.
#'   4) Return final subsets for Weighted CIF analysis (\code{Treatment.death, Control.death}, etc.).
#'
#' @param data_main A data.frame with ID, TTR, TTD, RECCNSR, DTHCNSR, baseline score, trt, etc.
#' @param data_long A data.frame with repeated clinical scores over time
#'   (e.g. ADYC, ORDSCOR).
#' @param wID_main Name of the patient ID column in the main dataset (default "USUBJID").
#' @param wTimeToRecovery_main Name of the time-to-recovery column (default "TTRECOV").
#' @param wTimeToDeath_main Name of the time-to-death column (default "TTDEATH").
#' @param wRecov_Censoring_main Name of the recovery-censor column (default "RECCNSR").
#' @param wDeath_Censoring_main Name of the death-censor column (default "DTHCNSR").
#' @param wBaselineScore_main Name of the baseline ordinal column (default "ordscr_bs").
#' @param wTreatment_main Name of the treatment indicator column (0=control,1=treatment). Default "trt".
#'
#' @param wID_long Name of the patient ID column in the long dataset (default "USUBJID").
#' @param wADY_long Name of the day-since-treatment column in the long dataset (default "ADYC").
#' @param wScore_long Name of the ordinal score column in the long dataset (default "ORDSCOR").
#'
#' @param wStates_death Vector of ordinal states for death weighting (default c(4,5,6,7)).
#' @param wWeights_death Numeric weights, same length as wStates_death (default c(2,1.5,1,0.5)).
#' @param wStates_discharge Vector of states for discharge weighting (default c(4,5,6,7)).
#' @param wWeights_discharge Numeric weights, same length as wStates_discharge
#'   (default c(0.5,1,1.5,2)).
#' @importFrom survival survfit Surv
#' @importFrom stats var cov pnorm
#' @return A list containing:
#'   \itemize{
#'     \item \code{data.ws.death} and \code{data.ws.discharge}: Full merged datasets
#'           with an added \code{wU} column for (death) or (discharge).
#'     \item \code{Treatment.death} and \code{Control.death}: Subsets for Weighted RMLT (death).
#'     \item \code{Treatment.discharge} and \code{Control.discharge}: Subsets for Weighted RMGT (discharge).
#'   }
#' @export
prep_data_weighted_cif2 <- function(
    data_main,
    data_long,

    wID_main              = "USUBJID",
    wTimeToRecovery_main  = "TTRECOV",
    wTimeToDeath_main     = "TTDEATH",
    wRecov_Censoring_main = "RECCNSR",
    wDeath_Censoring_main = "DTHCNSR",
    wBaselineScore_main   = "ordscr_bs",
    wTreatment_main       = "trt",

    wID_long              = "USUBJID",
    wADY_long             = "ADYC",
    wScore_long           = "ORDSCOR",

    wStates_death         = c(4,5,6,7),
    wWeights_death        = c(2,1.5,1,0.5),
    wStates_discharge     = c(4,5,6,7),
    wWeights_discharge    = c(0.5,1,1.5,2)
){

  # Quick validation
  if (length(wWeights_death) != length(wStates_death)) {
    stop("wWeights_death must be same length as wStates_death")
  }
  if (length(wWeights_discharge) != length(wStates_discharge)) {
    stop("wWeights_discharge must be same length as wStates_discharge")
  }

  # ============================
  # Part A) Exclude zero survival, missing baseline, discharge-to-die
  # ============================

  # 1) Exclude zero survival time
  cn.t0 <- which(
    data_main[[wTimeToRecovery_main]] == 0 |
      data_main[[wTimeToDeath_main]] == 0
  )

  # 2) Patients missing baseline score
  cn.noob <- which(is.na(data_main[[wBaselineScore_main]]))

  # 3) Discharge-to-die => set recov_censor=1, push TTR to max
  cn.dtd <- which(
    data_main[[wRecov_Censoring_main]] == 0 &
      data_main[[wDeath_Censoring_main]] == 0
  )
  if(length(cn.dtd) > 0){
    id.dtd <- data_main[[wID_main]][cn.dtd]
    myMaxTime <- max(data_main[[wTimeToRecovery_main]], na.rm=TRUE)
    data_main[data_main[[wID_main]] %in% id.dtd, wRecov_Censoring_main] <- 1
    data_main[data_main[[wID_main]] %in% id.dtd, wTimeToRecovery_main]  <- myMaxTime
  }

  # 4) Create etime
  data_main$etime <- pmin(
    data_main[[wTimeToRecovery_main]],
    data_main[[wTimeToDeath_main]]
  )

  # Exclude zero survival or missing baseline
  cn.discard <- c(cn.t0, cn.noob)
  if(length(cn.discard) > 0){
    data_main <- data_main[-cn.discard, ]
  }

  # Censoring indicator => 1 if event occurred, else 0
  data_main$estatus <- 1 - (
    (data_main[[wRecov_Censoring_main]] == 1) *
      (data_main[[wDeath_Censoring_main]] == 1)
  )

  # etype2 => 1 = recovery, 2 = death
  data_main$etype2 <-
    1 * (data_main[[wRecov_Censoring_main]]==0 &
           data_main[[wDeath_Censoring_main]]==1) +
    2 * (data_main[[wDeath_Censoring_main]]==0)

  # final data.w => select (ID, etime, estatus, etype2, trt)
  data.w <- data_main %>%
    dplyr::select(
      all_of(wID_main),
      "etime",
      "estatus",
      "etype2",
      all_of(wTreatment_main)
    )
  colnames(data.w) <- c("USUBJID","etime","estatus","etype2","trt")

  # ============================
  # Part B) Merge with long data
  # ============================

  # convert ID col to character
  data_long[[wID_long]] <- as.character(data_long[[wID_long]])

  # keep only matching IDs
  data_long <- data_long[data_long[[wID_long]] %in% data.w$USUBJID, ]

  # rename that column => USUBJID
  data_long <- data_long %>%
    dplyr::rename(USUBJID = all_of(wID_long))

  # left_join => data.long
  data.long <- dplyr::left_join(data_long, data.w, by="USUBJID")

  # replicate original "data.l"
  data.l <- data.long %>%
    dplyr::select("USUBJID","etime","estatus",
                  all_of(wADY_long),
                  all_of(wScore_long))

  # ensure ADYC is numeric, fill NA with 0
  data.l[[wADY_long]] <- suppressWarnings(as.numeric(data.l[[wADY_long]]))
  data.l[[wADY_long]][is.na(data.l[[wADY_long]])] <- 0

  data.l <- data.l %>%
    dplyr::arrange(USUBJID, !!rlang::sym(wADY_long))

  # keep only rows where etime >= ADYC
  data.l <- data.l[data.l$etime >= data.l[[wADY_long]], ]

  # rename => c("USUBJID","D_time","D_status","resp_time","resp")
  colnames(data.l) <- c("USUBJID","D_time","D_status","resp_time","resp")

  # ============================
  # Part C) FLEXIBLE VERSION: Replace hardcoded 4:7 with general lookup
  # ============================

  idList <- unique(data.l$USUBJID)
  wU <- data.frame(
    USUBJID = idList,
    death_w = 0,
    disc_w  = 0,
    stringsAsFactors = FALSE
  )

  for(i in seq_along(idList)){
    data.id <- data.l[data.l$USUBJID == idList[i], ]

    # FLEXIBLE: Initialize count arrays for any state range
    count_death <- rep(0, length(wStates_death))
    count_discharge <- rep(0, length(wStates_discharge))

    if(nrow(data.id) == 1){
      # ORIGINAL LOGIC: single row => add 1 day
      current_score <- data.id$resp[1]

      # FLEXIBLE: Use match() instead of hardcoded arithmetic
      death_idx <- match(current_score, wStates_death)
      if (!is.na(death_idx)) {
        count_death[death_idx] <- 1
      }

      discharge_idx <- match(current_score, wStates_discharge)
      if (!is.na(discharge_idx)) {
        count_discharge[discharge_idx] <- 1
      }

    } else {
      # ORIGINAL LOGIC: multi row => sum intervals (excluding final stretch)
      for(j in seq_len(nrow(data.id) - 1)){
        current_score <- data.id$resp[j]
        interval_length <- data.id$resp_time[j+1] - data.id$resp_time[j]

        # FLEXIBLE: Use match() instead of hardcoded arithmetic
        death_idx <- match(current_score, wStates_death)
        if (!is.na(death_idx)) {
          count_death[death_idx] <- count_death[death_idx] + interval_length
        }

        discharge_idx <- match(current_score, wStates_discharge)
        if (!is.na(discharge_idx)) {
          count_discharge[discharge_idx] <- count_discharge[discharge_idx] + interval_length
        }
      }
    }

    # ORIGINAL LOGIC: Direct weighted sum
    wU$death_w[i] <- sum(wWeights_death * count_death)
    wU$disc_w[i] <- sum(wWeights_discharge * count_discharge)
  }

  # ============================
  # Part D) Merge wU back, create treatment/control subsets
  # ============================

  data.ws.death <- dplyr::left_join(
    data.w, wU[, c("USUBJID","death_w")], by="USUBJID"
  ) %>%
    dplyr::rename(wU = death_w)

  data.ws.discharge <- dplyr::left_join(
    data.w, wU[, c("USUBJID","disc_w")], by="USUBJID"
  ) %>%
    dplyr::rename(wU = disc_w)

  # -- Subset: Treatment.death
  Treatment.death <- data.ws.death[data.ws.death$trt == 1,
                                   c("USUBJID","etime","estatus","etype2","wU")]
  colnames(Treatment.death) <- c("cn","D_time","D_status","etype","wU")

  # -- Subset: Control.death
  Control.death <- data.ws.death[data.ws.death$trt == 0,
                                 c("USUBJID","etime","estatus","etype2","wU")]
  colnames(Control.death) <- c("cn","D_time","D_status","etype","wU")

  # -- Subset: Treatment.discharge
  Treatment.discharge <- data.ws.discharge[data.ws.discharge$trt == 1,
                                           c("USUBJID","etime","estatus","etype2","wU")]
  colnames(Treatment.discharge) <- c("cn","D_time","D_status","etype","wU")

  # -- Subset: Control.discharge
  Control.discharge <- data.ws.discharge[data.ws.discharge$trt == 0,
                                         c("USUBJID","etime","estatus","etype2","wU")]
  colnames(Control.discharge) <- c("cn","D_time","D_status","etype","wU")

  # Final return
  list(
    data.ws.death     = data.ws.death,
    data.ws.discharge = data.ws.discharge,

    Treatment.death    = Treatment.death,
    Control.death      = Control.death,
    Treatment.discharge = Treatment.discharge,
    Control.discharge   = Control.discharge
  )
}

# ============================
# Analysis wrapper function
# ============================

do_weighted_cif_analysis <- function(
    data_main,
    data_long,

    # Data preparation parameters
    wID_main              = "USUBJID",
    wTimeToRecovery_main  = "TTRECOV",
    wTimeToDeath_main     = "TTDEATH",
    wRecov_Censoring_main = "RECCNSR",
    wDeath_Censoring_main = "DTHCNSR",
    wBaselineScore_main   = "ordscr_bs",
    wTreatment_main       = "trt",
    wID_long              = "USUBJID",
    wADY_long             = "ADYC",
    wScore_long           = "ORDSCOR",

    # Flexible state and weight parameters
    wStates_death         = c(4,5,6,7),
    wWeights_death        = c(2,1.5,1,0.5),
    wStates_discharge     = c(4,5,6,7),
    wWeights_discharge    = c(0.5,1,1.5,2),

    # Analysis parameters
    tau                   = c(15, 29)
){

  # Prepare data
  prep_result <- prep_data_weighted_cif2(
    data_main = data_main,
    data_long = data_long,
    wID_main = wID_main,
    wTimeToRecovery_main = wTimeToRecovery_main,
    wTimeToDeath_main = wTimeToDeath_main,
    wRecov_Censoring_main = wRecov_Censoring_main,
    wDeath_Censoring_main = wDeath_Censoring_main,
    wBaselineScore_main = wBaselineScore_main,
    wTreatment_main = wTreatment_main,
    wID_long = wID_long,
    wADY_long = wADY_long,
    wScore_long = wScore_long,
    wStates_death = wStates_death,
    wWeights_death = wWeights_death,
    wStates_discharge = wStates_discharge,
    wWeights_discharge = wWeights_discharge
  )

  # Run analysis for each tau
  results <- lapply(tau, function(tt) {
    list(
      WRMLT = table_weighted(prep_result$Treatment.death,
                             prep_result$Control.death,
                             eta = 2, tau = tt),  # eta=2 for death
      WRMGT = table_weighted(prep_result$Treatment.discharge,
                             prep_result$Control.discharge,
                             eta = 1, tau = tt)   # eta=1 for recovery/discharge
    )
  })

  # Add tau labels for clarity
  names(results) <- paste0("tau_", tau)

  return(results)
}

# ============================
# Helper function to create common weight patterns
# ============================

create_weights <- function(states, pattern = "decreasing") {
  n <- length(states)
  weights <- switch(pattern,
                    "decreasing" = seq(2, 0.5, length.out = n),        # Worse states get higher weight
                    "increasing" = seq(0.5, 2, length.out = n),        # Better states get higher weight
                    "equal" = rep(1, n),                               # All states equal weight
                    "binary_high_low" = if(n == 2) c(2, 0.5) else stop("binary_high_low only for 2 states"),
                    "binary_low_high" = if(n == 2) c(0.5, 2) else stop("binary_low_high only for 2 states"),
                    stop("Unknown pattern. Use: decreasing, increasing, equal, binary_high_low, binary_low_high")
  )

  return(weights)
}

###
# prepped_w <- prep_data_weighted_cif2(
#   data_main = main_df,
#   data_long = long_df,
#
#   wID_main              = "ID",
#   wTimeToRecovery_main  = "TimeToRecovery",
#   wTimeToDeath_main     = "TimeToDeath",
#   wRecov_Censoring_main = "RecoveryCensoringIndicator",
#   wDeath_Censoring_main = "DeathCensoringIndicator",
#   wTreatment_main       = "Treatment",
#   wBaselineScore_main   = "BaselineScore",
#
#   wID_long              = "PersonID",
#   wADY_long             = "RelativeDay",
#   wScore_long           = "OrdinalScore",
#
#   wStates_death         = c(4,5,6,7),
#   wWeights_death        = c(2,1.5,1,0.5),
#   wStates_discharge     = c(4,5,6,7),
#   wWeights_discharge    = c(0.5,1,1.5,2)
# )
