#' Prepare Data (Part II: Weighted CIF)
#'
#' @description
#' Extends the cleaned wide data (from \code{\link{prep_data_part1}}) by merging
#' a long-format dataset and computing weighted times for death and discharge.
#' This is only relevant if the user has repeated measurements (e.g., ordinal scores).
#'
#' @param part1_output  The list returned by \code{prep_data_part1()}, containing
#'   \code{df_ch} among other elements.
#' @param path_long  (Optional) Path to the long-format CSV. If \code{NULL}, uses
#'   \code{dummy_actt_long} instead (for demonstration).
#' @param col_map_long A named list specifying columns in the long data:
#'   \describe{
#'     \item{\code{adyc}}{Measurement day/time.}
#'     \item{\code{ordscor}}{Ordinal score.}
#'   }
#'
#' @return A list with the Weighted CIF elements:
#' \itemize{
#'   \item \code{data.ws.death}, \code{data.ws.discharge}: full wide+long merges.
#'   \item \code{trt1.death}, \code{trt0.death}: for event=death subsets.
#'   \item \code{trt1.discharge}, \code{trt0.discharge}: for event=discharge subsets.
#' }
#'
#' @export
prep_data_wrt_cif <- function(
    part1_output,
    path_long = NULL,
    col_map_long = list(
      adyc    = "ADYC",
      ordscor = "ORDSCOR"
    )
) {
  # 1) Extract the cleaned wide data from part1
  df_ch <- part1_output$df_ch
  if (is.null(df_ch)) {
    stop("part1_output$df_ch is NULL. Did you provide the correct object from prep_data_part1()?")
  }

  # 2) Load or read the long data
  if (is.null(path_long)) {
    message("No path_long provided, using dummy_actt_long for Part II.")
    data("dummy_long", package = "jointCompRisk", envir = environment())
    df_long <- dummy_long
  } else {
    message("Reading user-supplied long data from: ", path_long)
    df_long <- utils::read.csv(path_long, stringsAsFactors = FALSE)
  }

  # 3) Validate needed columns in df_long
  needed_cols <- c("USUBJID", col_map_long$adyc, col_map_long$ordscor)
  missing_cols <- setdiff(needed_cols, colnames(df_long))
  if (length(missing_cols) > 0) {
    stop("Missing columns in long data: ", paste(missing_cols, collapse=", "))
  }

  # 4) Minimal wide columns for merging
  minimal_merge <- data.frame(
    USUBJID = df_ch$USUBJID,
    etime   = df_ch$etime,
    estatus = df_ch$estatus,
    etype2  = df_ch$etype2,
    trt     = df_ch$trt
  )

  # Keep only those in the long data who appear in wide
  df_long <- subset(df_long, df_long$USUBJID %in% minimal_merge$USUBJID)

  # Merge
  suppressWarnings(requireNamespace("dplyr", quietly = TRUE))
  library(dplyr)
  data.l <- left_join(df_long, minimal_merge, by = "USUBJID")

  # Weighted times function
  compute_weighted_times <- function(subdf) {
    times  <- subdf[[col_map_long$adyc]]
    scores <- subdf[[col_map_long$ordscor]]

    ord_idx <- order(times)
    times   <- times[ord_idx]
    scores  <- scores[ord_idx]

    # example weighting
    w_death_map <- c(`4`=2, `5`=1.5, `6`=1, `7`=0.5)
    w_disc_map  <- c(`4`=0.5, `5`=1, `6`=1.5, `7`=2)

    wtime_death     <- 0
    wtime_discharge <- 0

    if (length(times) > 1) {
      for (i in seq_len(length(times) - 1)) {
        s_char  <- as.character(scores[i])
        interval <- times[i+1] - times[i]
        wtime_death     <- wtime_death + w_death_map[s_char] * interval
        wtime_discharge <- wtime_discharge + w_disc_map[s_char] * interval
      }
    }
    data.frame(wU_death = wtime_death, wU_discharge = wtime_discharge)
  }

  wU_matrix <- data.l %>%
    group_by(USUBJID) %>%
    do( compute_weighted_times(.) ) %>%
    ungroup()

  # Combine with minimal_merge
  data.ws.death     <- left_join(minimal_merge, wU_matrix, by="USUBJID")
  data.ws.discharge <- data.ws.death

  colnames(data.ws.death)[     colnames(data.ws.death)     == "wU_death"     ] <- "wU"
  colnames(data.ws.discharge)[ colnames(data.ws.discharge) == "wU_discharge" ] <- "wU"

  # Split by treatment, rename columns for Weighted CIF
  trt1.death <- subset(data.ws.death, trt == 1,
                       select=c("USUBJID","etime","estatus","etype2","wU"))
  trt0.death <- subset(data.ws.death, trt == 0,
                       select=c("USUBJID","etime","estatus","etype2","wU"))

  trt1.discharge <- subset(data.ws.discharge, trt == 1,
                           select=c("USUBJID","etime","estatus","etype2","wU"))
  trt0.discharge <- subset(data.ws.discharge, trt == 0,
                           select=c("USUBJID","etime","estatus","etype2","wU"))

  colnames(trt1.death)     <- c("cn","D_time","D_status","etype","wU")
  colnames(trt0.death)     <- c("cn","D_time","D_status","etype","wU")
  colnames(trt1.discharge) <- c("cn","D_time","D_status","etype","wU")
  colnames(trt0.discharge) <- c("cn","D_time","D_status","etype","wU")

  # Return Part II results
  list(
    data.ws.death     = data.ws.death,
    data.ws.discharge = data.ws.discharge,
    trt1.death        = trt1.death,
    trt0.death        = trt0.death,
    trt1.discharge    = trt1.discharge,
    trt0.discharge    = trt0.discharge
  )
}
