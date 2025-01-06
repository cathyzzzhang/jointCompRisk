#' Prepare ACTT Trial Data (Using Dummy Data)
#'
#' @description
#' Demonstrates the data cleaning steps for the ACTT trial analysis,
#' both standard CIF (Part I) and weighted CIF (Part II),
#' using the \code{dummy_actt_ch} and \code{dummy_actt_long} datasets.
#'
#' @details
#' Steps included:
#' \enumerate{
#'   \item Combine and exclude invalid patients or missing observations as needed.
#'   \item Create competing-risks variables (\code{etime, estatus, etype2}).
#'   \item Split data by treatment group for standard CIF (Part I).
#'   \item Merge long-format data to compute weighted times for "death" and "discharge" (Part II).
#' }
#'
#' \return A named list with the following elements:
#' \describe{
#'   \item{\code{data.w}}{Cleaned dataset (wide) for CIF analysis (Part I).}
#'   \item{\code{trt1}, \code{trt0}}{Treatment vs Control subsets of \code{data.w}.}
#'   \item{\code{data.ws.death}, \code{data.ws.discharge}}{Wide data frames with weighted times.}
#'   \item{\code{trt1.death}, \code{trt0.death}}{Treatment vs Control subsets for Weighted CIF (death).}
#'   \item{\code{trt1.discharge}, \code{trt0.discharge}}{Treatment vs Control subsets for Weighted CIF (discharge).}
#' }
#'
#' @export
#'
#' @examples
#' data("dummy_actt_ch")
#' data("dummy_actt_long")
#'
#' result <- prepare_actt_data()
#' names(result)
#' head(result$data.w)
prepare_actt_data <- function() {
  #---------------------------------------------------------------------------
  # 1. Load dummy data from package environment
  #---------------------------------------------------------------------------
  data("dummy_actt_ch", envir = environment())
  data("dummy_actt_long", envir = environment())

  actt <- dummy_actt_ch
  data_long <- dummy_actt_long

  #---------------------------------------------------------------------------
  # 2. Potential exclusions or transformations
  #    (In real analysis, you might remove patients with zero survival time,
  #     or missing baseline scores, etc. We'll keep it simpler here.)
  #---------------------------------------------------------------------------

  # Example: define "zero survival time" as TTRECOV == 0 or TTDEATH == 0
  exclude_zero_surv <- which(actt$TTRECOV == 0 | actt$TTDEATH == 0)
  if(length(exclude_zero_surv) > 0) {
    actt <- actt[-exclude_zero_surv, ]
  }

  # Example: define missing baseline ordinal score => exclude
  exclude_missing_score <- which(is.na(actt$ordscr_bs))
  if(length(exclude_missing_score) > 0) {
    actt <- actt[-exclude_missing_score, ]
  }

  # Maybe handle "discharge-to-die" logic: if both RECCNSR and DTHCNSR are 0 => treat them as death
  # This is just an example consistent with your original code
  dtd_idx <- which(actt$RECCNSR == 0 & actt$DTHCNSR == 0)
  actt[dtd_idx, "RECCNSR"]  <- 1
  actt[dtd_idx, "TTRECOV"]  <- 28

  #---------------------------------------------------------------------------
  # 3. Create standard CIF data (Part I)
  #---------------------------------------------------------------------------
  actt$etime <- pmin(actt$TTRECOV, actt$TTDEATH)

  # estatus: 1 if event observed, 0 if censored
  # Example: "event" if at least one is uncensored
  actt$estatus <- 1 - ((actt$RECCNSR == 1) & (actt$DTHCNSR == 1))

  # etype2: 1 = recovery, 2 = death
  actt$etype2 <- 1*(actt$RECCNSR == 0 & actt$DTHCNSR == 1) +
    2*(actt$DTHCNSR == 0)

  # The final wide dataset used for standard CIF
  data.w <- data.frame(
    USUBJID = actt$USUBJID,
    etime   = actt$etime,
    estatus = actt$estatus,
    etype2  = actt$etype2,
    trt     = actt$trt
  )

  # Split by treatment
  trt1 <- subset(data.w, trt == 1)
  trt0 <- subset(data.w, trt == 0)

  #---------------------------------------------------------------------------
  # 4. Weighted CIF data (Part II)
  #    We'll join data_long with the relevant columns from actt
  #---------------------------------------------------------------------------

  # Minimal columns from 'actt' to merge on USUBJID
  data_for_merge <- data.frame(
    USUBJID = actt$USUBJID,
    etime   = actt$etime,
    estatus = actt$estatus,
    etype2  = actt$etype2,
    trt     = actt$trt
  )

  # Keep only matching USUBJID in data_long
  data_long <- subset(data_long, data_long$USUBJID %in% data_for_merge$USUBJID)

  # Join
  library(dplyr)
  data.l <- left_join(data_long, data_for_merge, by = "USUBJID")

  # Possibly rename columns if needed:
  # Suppose ADYC is the "day" or "resp_time," ORDSCOR is the "resp" score
  data.l <- data.l %>%
    rename(D_time_original = etime) %>%  # keep the old name for reference
    mutate(
      # We'll preserve 'etime' from the merged data_for_merge as 'D_time_original'
      # but keep ADYC as 'resp_time', ORDSCOR as 'resp'
      resp_time = ADYC,
      resp      = ORDSCOR
    )

  # In real code, you might do something like ensuring resp_time <= D_time_original, etc.

  #---------------------------------------------------------------------------
  # 4a. Compute Weighted Times
  #    Example weighting:
  #    - time to death:     score 4->2, 5->1.5, 6->1, 7->0.5
  #    - time to discharge: score 4->0.5, 5->1, 6->1.5, 7->2
  #---------------------------------------------------------------------------

  # For demonstration, let's create a function that calculates weighted times for each ID
  compute_weighted_times <- function(subdf) {
    # subdf is all rows for a single USUBJID
    # we’ll accumulate time spent in each ordinal score
    # NOTE: This is a simplified approach
    scores <- subdf$resp
    times  <- subdf$resp_time

    # Sort by resp_time
    ord <- order(times)
    scores <- scores[ord]
    times  <- times[ord]

    # For consecutive intervals
    wtime_death      <- 0
    wtime_discharge  <- 0

    # Weighted values for each score
    w_death_map <- c(`4`=2, `5`=1.5, `6`=1, `7`=0.5)
    w_disc_map  <- c(`4`=0.5, `5`=1, `6`=1.5, `7`=2)

    if(length(times) > 1) {
      for(i in seq_len(length(times) - 1)) {
        score_i <- as.character(scores[i])
        interval <- (times[i+1] - times[i])
        wtime_death     <- wtime_death     + w_death_map[[score_i]] * interval
        wtime_discharge <- wtime_discharge + w_disc_map[[score_i]]  * interval
      }
    }
    # Possibly add last interval up to 'D_time_original' if you do that in your analysis
    # Or stop at the last measurement time

    data.frame(
      wU_death     = wtime_death,
      wU_discharge = wtime_discharge
    )
  }

  # Apply per USUBJID
  library(dplyr)
  wU_matrix <- data.l %>%
    group_by(USUBJID) %>%
    do( compute_weighted_times(.) ) %>%
    ungroup()

  # Merge these weighted times back
  data.ws.death <- left_join(data_for_merge, wU_matrix, by="USUBJID")
  data.ws.discharge <- data.ws.death  # copy structure

  # Rename columns for clarity
  colnames(data.ws.death)[ colnames(data.ws.death) == "wU_death" ] <- "wU"
  colnames(data.ws.discharge)[ colnames(data.ws.discharge) == "wU_discharge" ] <- "wU"

  # Split by treatment
  trt1.death <- subset(data.ws.death, trt == 1,
                       select=c("USUBJID","etime","estatus","etype2","wU"))
  trt0.death <- subset(data.ws.death, trt == 0,
                       select=c("USUBJID","etime","estatus","etype2","wU"))

  trt1.discharge <- subset(data.ws.discharge, trt == 1,
                           select=c("USUBJID","etime","estatus","etype2","wU"))
  trt0.discharge <- subset(data.ws.discharge, trt == 0,
                           select=c("USUBJID","etime","estatus","etype2","wU"))

  # Rename columns for usage in your analysis
  colnames(trt1.death)      <- c("cn","D_time","D_status","etype","wU")
  colnames(trt0.death)      <- c("cn","D_time","D_status","etype","wU")
  colnames(trt1.discharge)  <- c("cn","D_time","D_status","etype","wU")
  colnames(trt0.discharge)  <- c("cn","D_time","D_status","etype","wU")

  #---------------------------------------------------------------------------
  # Return everything
  #---------------------------------------------------------------------------
  list(
    data.w            = data.w,
    trt1              = trt1,
    trt0              = trt0,
    data.ws.death     = data.ws.death,
    data.ws.discharge = data.ws.discharge,
    trt1.death        = trt1.death,
    trt0.death        = trt0.death,
    trt1.discharge    = trt1.discharge,
    trt0.discharge    = trt0.discharge
  )
}
