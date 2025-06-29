#' @title Prepare Data for Standard CIF
#' @description Cleans and prepares a single dataset for standard (competing risks) CIF analysis.
#'
#' @param data A data frame with columns for ID, time to recovery, time to death,
#'   recovery censor, death censor, and treatment indicator.
#' @param ID Name of the patient ID column. Default is "USUBJID".
#' @param TimeToRecovery Name of the time-to-recovery column. Default "TTRECOV".
#' @param TimeToDeath Name of the time-to-death column. Default "TTDEATH".
#' @param Recov_Censoring Name of the recovery-censor column. Default "RECCNSR" (0=event,1=censor).
#' @param Death_Censoring Name of the death-censor column. Default "DTHCNSR" (0=event,1=censor).
#' @param Treatment Name of the treatment indicator column (0=control,1=treatment). Default "trt".
#'
#' @return A list with:
#'   \itemize{
#'     \item \code{data.w}: The processed data frame with columns \code{cn, etime, estatus, etype2, Treatment}.
#'     \item \code{Treatment}: Subset of \code{data.w} where \code{Treatment==1}.
#'     \item \code{Control}: Subset of \code{data.w} where \code{Treatment==0}.
#'   }
#' @import survival
#' @export
prep_data_cif <- function(
    data,
    ID = "USUBJID",
    TimeToRecovery = "TTRECOV",
    TimeToDeath = "TTDEATH",
    Recov_Censoring = "RECCNSR",
    Death_Censoring = "DTHCNSR",
    Treatment = "trt"
){
  # Exclude zero survival times
  cn.t0 <- which(data[[TimeToRecovery]] == 0 | data[[TimeToDeath]] == 0)

  # Handle discharge-to-die cases
  cn.dtd <- which(data[[Recov_Censoring]] == 0 & data[[Death_Censoring]] == 0)
  if (length(cn.dtd) > 0) {
    id.dtd <- data[[ID]][cn.dtd]
    myMaxTime <- max(data[[TimeToRecovery]], na.rm = TRUE)
    data[data[[ID]] %in% id.dtd, Recov_Censoring] <- 1
    data[data[[ID]] %in% id.dtd, TimeToRecovery]  <- myMaxTime
  }

  # Create event time
  data$etime <- pmin(data[[TimeToRecovery]], data[[TimeToDeath]])

  # Only exclude zero survival times
  if (length(cn.t0) > 0) {
    data <- data[-cn.t0, ]
  }

  # Create event status and type
  data$estatus <- 1 - ((data[[Recov_Censoring]] == 1) * (data[[Death_Censoring]] == 1))
  data$etype2 <- 1 * (data[[Recov_Censoring]] == 0 & data[[Death_Censoring]] == 1) +
    2 * (data[[Death_Censoring]] == 0)

  # Select final columns
  data.w <- data %>%
    dplyr::select(
      all_of(ID),
      "etime",
      "estatus",
      "etype2",
      all_of(Treatment)
    )
  colnames(data.w) <- c("cn", "etime", "estatus", "etype2", "Treatment")

  # Split by treatment
  Treatment_sub <- data.w[data.w$Treatment == 1, ]
  Control_sub   <- data.w[data.w$Treatment == 0, ]

  list(
    data.w    = data.w,
    Treatment = Treatment_sub,
    Control   = Control_sub
  )
}
