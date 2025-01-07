#' Dummy Data for CIF Analysis (De-Identified)
#'
#' @description
#' These synthetic datasets mimic the structure required for CIF (Cumulative Incidence Function)
#' and weighted CIF analyses, but use randomly generated values for demonstration purposes.
#'
#' \itemize{
#'   \item \code{dummy_ch}: The "wide" dataset containing patient-level variables.
#'   \item \code{dummy_long}: The "long" dataset containing repeated measurements over time.
#' }
#'
#' @format
#' \describe{
#'   \item{\code{dummy_ch}}{A data frame with 15 rows (example) and 40 columns, used for standard CIF analysis.}
#'   \item{\code{dummy_long}}{A data frame with 50 rows (example) and 15 columns, used for weighted CIF analysis.}
#' }
#'
#' @usage
#' data("dummy_ch")
#' data("dummy_long")
#'
#' @examples
#' # Load the dummy datasets
#' data("dummy_ch")
#' data("dummy_long")
#'
#' # Preview data
#' head(dummy_ch)
#' head(dummy_long)
"dummy_ch"

#' Dummy Data for Weighted CIF Analysis (Long Format)
#'
#' @description
#' A synthetic dataset simulating the long-format structure used for weighted CIF analysis.
#'
#' @format A data frame with 50 rows and 15 columns:
#' \describe{
#'   \item{USUBJID}{Patient identifier.}
#'   \item{ARM}{Treatment arm (e.g., Placebo or Treatment).}
#'   \item{ACTARM}{Actual treatment arm.}
#'   \item{BCSOSN}{Baseline ordinal score (4–7).}
#'   \item{D29DTHE0, D29DTHE1}{Day 29 ordinal score values.}
#'   \item{TTRECOV0, TTRECOV1}{Time to recovery (days).}
#'   \item{ADYC}{Assessment day.}
#'   \item{ORDSCOR}{Ordinal score during follow-up.}
#'   \item{AGEC}{Age of patient.}
#'   \item{SEX}{Gender.}
#'   \item{BDURSYMP}{Baseline symptom duration.}
#'   \item{COMORB2}{Comorbidity categories.}
#' }
#'
#' @examples
#' # Load the long-format dataset
#' data(dummy_long)
#' head(dummy_long)
"dummy_long"

