% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cif_inference.R
\name{auc.func}
\alias{auc.func}
\title{Compute AUC under a cumulative incidence function}
\usage{
auc.func(table, start.time, tau)
}
\arguments{
\item{table}{A two-column matrix/data frame where the first column is time and the second column is the CIF estimate.}

\item{start.time}{The lower time bound.}

\item{tau}{The upper time bound.}
}
\value{
A numeric value of the AUC.
}
\description{
Computes the area under the curve (AUC) from a \code{survfit} output for a specific
cumulative incidence function (CIF) between a given time window [\code{start.time}, \code{tau}].
}
