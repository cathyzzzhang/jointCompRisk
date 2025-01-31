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
#'
#' @return A list containing:
#'   \itemize{
#'     \item \code{data.ws.death} and \code{data.ws.discharge}: Full merged datasets
#'           with an added \code{wU} column for (death) or (discharge).
#'     \item \code{Treatment.death} and \code{Control.death}: Subsets for Weighted RMLT (death).
#'     \item \code{Treatment.discharge} and \code{Control.discharge}: Subsets for Weighted RMGT (discharge).
#'   }
#' @export
prep_data_weighted_cif <- function(
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
  # Part C) Loop to compute weighted times
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

    # "count" for states 4,5,6,7
    count <- rep(0, 4)  # index => [1]=state4, [2]=state5, [3]=state6, [4]=state7

    if(nrow(data.id) == 1){
      # single row => add 1 day
      count[data.id$resp - 3] <- 1
    } else {
      # multi row => sum intervals
      for(j in seq_len(nrow(data.id) - 1)){
        count[data.id$resp[j] - 3] <-
          count[data.id$resp[j] - 3] +
          data.id$resp_time[j+1] - data.id$resp_time[j]
      }
    }

    # Weighted time for death
    wU$death_w[i] <-
      wWeights_death[1] * count[1] +
      wWeights_death[2] * count[2] +
      wWeights_death[3] * count[3] +
      wWeights_death[4] * count[4]

    # Weighted time for discharge
    wU$disc_w[i] <-
      wWeights_discharge[1] * count[1] +
      wWeights_discharge[2] * count[2] +
      wWeights_discharge[3] * count[3] +
      wWeights_discharge[4] * count[4]
  }

  # ============================
  # Part D) Merge wU back, now separate (Treatment vs Control) for death and discharge
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

