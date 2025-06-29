#' Main Competing Risks Dataset
#'
#' Simulated clinical trial data with competing risks survival outcomes.
#' This dataset follows the structure of Adaptive COVID-19 Treatment Trials (ACTT)
#' with built-in treatment effects for demonstration purposes.
#'
#' @format A data frame with 150 rows and 7 variables:
#' \describe{
#'   \item{ID}{Patient identifier (character)}
#'   \item{TimeToRecovery}{Time to recovery event in days (numeric)}
#'   \item{TimeToDeath}{Time to death event in days (numeric)}
#'   \item{RecoveryCensoringIndicator}{Recovery censoring indicator (0=event observed, 1=censored)}
#'   \item{DeathCensoringIndicator}{Death censoring indicator (0=event observed, 1=censored)}
#'   \item{BaselineScore}{Baseline severity score, range 4-7 (numeric)}
#'   \item{Treatment}{Treatment arm indicator (0=control, 1=treatment)}
#' }
#' @details
#' This is a simulated dataset created for demonstration purposes with realistic
#' treatment effects built in: treatment group has 1.5× faster recovery times and
#' 1.8× improved survival compared to control. The data represents a clinical trial
#' with competing risks where patients can either recover or die, with administrative
#' censoring at 30 days.
#' @source Simulated data based on Weibull distributions with treatment-specific parameters
#' @examples
#' data(main_df)
#' head(main_df)
#' summary(main_df)
#' # Compare outcomes by treatment
#' tapply(main_df$TimeToRecovery, main_df$Treatment, summary)
#' tapply(main_df$TimeToDeath, main_df$Treatment, summary)
"main_df"

#' Longitudinal Severity Scores Dataset
#'
#' Repeated measurements of ordinal severity scores over time for the same patients
#' in the main_df dataset, with treatment-specific trajectory patterns.
#'
#' @format A data frame with variable rows per patient:
#' \describe{
#'   \item{PersonID}{Patient identifier matching ID in main_df (character)}
#'   \item{OrdinalScore}{Severity score on 1-8 scale (numeric)}
#'   \item{RelativeDay}{Study day (numeric) starting from day 0 (baseline)}
#' }
#' @details
#' Measurements are taken at scheduled visits: days 0 (baseline), 1, 3, 5, 7, 10, 14, 18, 21, 25, 28.
#' The trajectory follows treatment-specific probabilities: treatment patients have 45% improvement
#' and 15% worsening probability per visit, while control patients have 30% improvement and 25%
#' worsening probability, creating realistic differential clinical progression patterns.
#' @source Simulated data using treatment-specific random walk with boundaries
#' @examples
#' data(long_df)
#' head(long_df)
#' # See data for first patient
#' subset(long_df, PersonID == "Patient_001")
#' # Compare average scores by treatment
#' library(dplyr)
#' long_df %>%
#'   left_join(main_df[,c("ID","Treatment")], by=c("PersonID"="ID")) %>%
#'   group_by(Treatment) %>%
#'   summarise(mean_score = mean(OrdinalScore))
"long_df"
