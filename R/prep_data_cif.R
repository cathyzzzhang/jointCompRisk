#' Prepare Data (Part I: Standard CIF Only)
#'
#' @description
#' This function reads in (or uses built-in) wide-format data and
#' prepares it for standard CIF analysis (Part I).
#'
#' @param path_ch  (Optional) Path to the wide-format CSV. If \code{NULL},
#'   the built-in dummy data (\code{dummy_actt_ch}) is used.
#' @param col_map  Named list mapping columns in the wide dataset to needed variables:
#'   \describe{
#'     \item{\code{ttrecov}}{Time to recovery.}
#'     \item{\code{ttdeath}}{Time to death.}
#'     \item{\code{reccnsr}}{Indicator for censored/not for recovery.}
#'     \item{\code{dthcnsr}}{Indicator for censored/not for death.}
#'     \item{\code{ordscr_bs}}{Baseline ordinal score.}
#'     \item{\code{trt}}{Treatment group. (0=control, 1=treatment)}
#'   }
#'
#' @return A list with:
#' \itemize{
#'   \item \code{df_ch}: the cleaned wide dataframe.
#'   \item \code{data.w}: the final dataset for standard CIF with columns
#'         (\code{USUBJID, etime, estatus, etype2, trt}).
#'   \item \code{trt1}: subset of \code{data.w} where \code{trt==1}.
#'   \item \code{trt0}: subset of \code{data.w} where \code{trt==0}.
#' }
#'
#' @export
prep_data_cif <- function(
    path_ch = NULL,
    col_map = list(
      ttrecov   = "TTRECOV",
      ttdeath   = "TTDEATH",
      reccnsr   = "RECCNSR",
      dthcnsr   = "DTHCNSR",
      ordscr_bs = "ordscr_bs",
      trt       = "trt"
    )
) {
  # 1) Load wide-format data
  if (is.null(path_ch)) {
    message("No path_ch provided, using built-in dummy_actt_ch for Part I.")
    data("dummy_ch", package = "jointCompRisk", envir = environment())
    df_ch <- dummy_ch
  } else {
    message("Reading user-supplied wide data from: ", path_ch)
    df_ch <- utils::read.csv(path_ch, stringsAsFactors = FALSE)
  }

  # 2) Validate columns
  req_cols <- c(
    col_map$ttrecov,
    col_map$ttdeath,
    col_map$reccnsr,
    col_map$dthcnsr,
    col_map$ordscr_bs,
    col_map$trt
  )
  missing_cols <- setdiff(req_cols, colnames(df_ch))
  if (length(missing_cols) > 0) {
    stop("Missing required columns in wide data: ",
         paste(missing_cols, collapse=", "))
  }

  # 3) Basic Exclusions & Adjustments
  idx_zero <- which(df_ch[[col_map$ttrecov]] == 0 | df_ch[[col_map$ttdeath]] == 0)
  if (length(idx_zero) > 0) {
    df_ch <- df_ch[-idx_zero, ]
  }

  idx_miss <- which(is.na(df_ch[[col_map$ordscr_bs]]))
  if (length(idx_miss) > 0) {
    df_ch <- df_ch[-idx_miss, ]
  }

  # discharge-to-die logic
  dtd_idx <- which(df_ch[[col_map$reccnsr]] == 0 & df_ch[[col_map$dthcnsr]] == 0)
  if (length(dtd_idx) > 0) {
    df_ch[dtd_idx, col_map$reccnsr] <- 1
    df_ch[dtd_idx, col_map$ttrecov] <- 28
  }

  # 4) Create standard CIF variables
  df_ch$etime <- pmin(df_ch[[col_map$ttrecov]], df_ch[[col_map$ttdeath]])
  df_ch$estatus <- 1 - ((df_ch[[col_map$reccnsr]] == 1) & (df_ch[[col_map$dthcnsr]] == 1))
  df_ch$etype2  <- 1*(df_ch[[col_map$reccnsr]] == 0 & df_ch[[col_map$dthcnsr]] == 1) +
    2*(df_ch[[col_map$dthcnsr]] == 0)

  # The final dataset for standard CIF
  data.w <- data.frame(
    USUBJID = df_ch$USUBJID,
    etime   = df_ch$etime,
    estatus = df_ch$estatus,
    etype2  = df_ch$etype2,
    trt     = df_ch[[col_map$trt]]
  )

  trt1 <- subset(data.w, trt == 1)
  trt0 <- subset(data.w, trt == 0)

  # Return Part I results
  list(
    df_ch  = df_ch,   # cleaned wide data if user wants to see or pass on
    data.w = data.w,
    trt1   = trt1,
    trt0   = trt0
  )
}

