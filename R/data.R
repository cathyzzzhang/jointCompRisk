#' Dummy ACTT Trial Data (De-Identified)
#'
#' @description
#' These two synthetic datasets mimic the structure of the real ACTT trial data,
#' but use randomly generated values for demonstration purposes.
#'
#' \itemize{
#'   \item \code{dummy_actt_ch}: The "wide" data containing patient-level variables.
#'   \item \code{dummy_actt_long}: The "long" data containing repeated measurements over time.
#' }
#'
#' @format
#' \describe{
#'   \item{\code{dummy_actt_ch}}{A data frame with 15 rows (example) and 40 columns.}
#'   \item{\code{dummy_actt_long}}{A data frame with 50 rows (example) and 15 columns.}
#' }
#'
#' @usage
#' data("dummy_actt_ch")
#' data("dummy_actt_long")
#'
#' @examples
#' data("dummy_actt_ch")
#' data("dummy_actt_long")
#' head(dummy_actt_ch)
#' head(dummy_actt_long)
"dummy_actt_ch"

"dummy_actt_long"
