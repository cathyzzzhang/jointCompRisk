% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{main_df}
\alias{main_df}
\title{Main Competing Risks Dataset
Simulated clinical trial data with competing risks survival outcomes.
This dataset follows the structure of Adaptive COVID-19 Treatment Trials (ACTT)
with built-in treatment effects for demonstration purposes.}
\format{
A data frame with 150 rows and 7 variables:
\describe{
  \item{ID}{Patient identifier (character)}
  \item{TimeToRecovery}{Time to recovery event in days (numeric)}
  \item{TimeToDeath}{Time to death event in days (numeric)}
  \item{RecoveryCensoringIndicator}{Recovery censoring indicator (0=event observed, 1=censored)}
  \item{DeathCensoringIndicator}{Death censoring indicator (0=event observed, 1=censored)}
  \item{BaselineScore}{Baseline severity score, range 4-7 (numeric)}
  \item{Treatment}{Treatment arm indicator (0=control, 1=treatment)}
}
}
\source{
Simulated data based on Weibull distributions with treatment-specific parameters
}
\usage{
main_df
}
\description{
Main Competing Risks Dataset
Simulated clinical trial data with competing risks survival outcomes.
This dataset follows the structure of Adaptive COVID-19 Treatment Trials (ACTT)
with built-in treatment effects for demonstration purposes.
}
\details{
This is a simulated dataset created for demonstration purposes with realistic
treatment effects built in: treatment group has 1.5× faster recovery times and
1.8× improved survival compared to control. The data represents a clinical trial
with competing risks where patients can either recover or die, with administrative
censoring at 30 days.
}
\examples{
data(main_df)
head(main_df)
summary(main_df)
# Compare outcomes by treatment
tapply(main_df$TimeToRecovery, main_df$Treatment, summary)
tapply(main_df$TimeToDeath, main_df$Treatment, summary)
}
\keyword{datasets}
