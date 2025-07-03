#' @title Run Standard CIF Analysis
#' @description Given a prepped data list from \code{\link{prep_data_cif}}, run the standard CIF analysis.
#' @importFrom survival survfit Surv
#' @importFrom stats var cov pnorm
#' @param prepped A list returned by \code{prep_data_cif()}, containing \code{Treatment} and \code{Control}.
#' @param tau Numeric, time horizon (e.g. 15 or 29).
#'
#' @details
#' - RMGT uses parameters \code{(a,b,c) = (0,1,0)}.
#' - RMLT uses \code{(a,b,c) = (0,0,1)}.
#' - Net time uses \code{(a,b,c) = (0,1,-1)}.
#'
#' @return A list with elements \code{RMGT}, \code{RMLT}, and \code{Net}, each a 3x4 matrix.
#' @export
do_cif_analysis <- function(prepped, tau=15){
  Treatment <- prepped$Treatment
  Control <- prepped$Control

  # RMGT
  rmgt_tab <- table1_cif(Treatment, Control, tau, 0, 1, 0)
  # RMLT
  rmlt_tab <- table1_cif(Treatment, Control, tau, 0, 0, 1)
  # Net
  net_tab  <- table1_cif(Treatment, Control, tau, 0, 1, -1)

  list(
    RMGT = rmgt_tab,
    RMLT = rmlt_tab,
    Net  = net_tab
  )
}
