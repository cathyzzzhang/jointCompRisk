#' Main Competing Risks Dataset
#'
#' Simulated clinical trial data with competing risks survival outcomes.
#' This dataset follows the structure of Adaptive COVID-19 Treatment Trials (ACTT).
#'
#' @format A data frame with 30 rows and 7 variables:
#' \describe{
#'   \item{ID}{Patient identifier (character)}
#'   \item{TimeToRecovery}{Time to recovery event in days (numeric)}
#'   \item{TimeToDeath}{Time to death event in days (numeric)}
#'   \item{RecoveryCensoringIndicator}{Recovery censoring indicator (0=event observed, 1=censored)}
#'   \item{DeathCensoringIndicator}{Death censoring indicator (0=event observed, 1=censored)}
#'   \item{BaselineScore}{Baseline severity score, range 4-8 (numeric)}
#'   \item{Treatment}{Treatment arm indicator (0=control, 1=treatment)}
#' }
#' @details
#' This is a simulated dataset created for demonstration purposes. The data represents
#' a clinical trial with competing risks where patients can either recover or die,
#' with administrative censoring at 30 days.
#' @source Simulated data based on Weibull distributions for event times
#' @examples
#' data(main_df)
#' head(main_df)
#' summary(main_df)
"main_df"

#' Longitudinal Severity Scores Dataset
#'
#' Repeated measurements of ordinal severity scores over time for the same patients
#' in the main_df dataset.
#'
#' @format A data frame with variable rows per patient:
#' \describe{
#'   \item{PersonID}{Patient identifier matching ID in main_df (character)}
#'   \item{OrdinalScore}{Severity score on 1-8 scale (numeric)}
#'   \item{RelativeDay}{Study day (numeric) or "Baseline" (character)}
#' }
#' @details
#' Measurements are taken at scheduled visits: days 1, 2, 3, 4, 5, 7, 10, 14, 21, 28.
#' The trajectory follows a random walk pattern around the baseline score.
#' @source Simulated data using random walk with boundaries
#' @examples
#' data(long_df)
#' head(long_df)
#' # See data for first patient
#' subset(long_df, PersonID == "DUMMY.01")
"long_df"
