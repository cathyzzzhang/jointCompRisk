#' @title Run Weighted CIF Analysis
#' @description Given the list from \code{\link{prep_data_weighted_cif}}, which contains
#'   \code{Treatment.death, Control.death, Treatment.discharge, Control.discharge}, run Weighted RMLT (death)
#'   and Weighted RMGT (discharge) at a user-specified time horizon \code{tau}.
#'
#' @param prepped A list returned by \code{prep_data_weighted_cif()}.
#' @param tau Numeric time horizon (e.g., 15 or 29).
#' @importFrom survival survfit Surv
#' @importFrom stats var cov pnorm
#' @details
#' - Weighted RMLT uses \code{eta=2}.
#' - Weighted RMGT uses \code{eta=1}.
#'
#' @return A list with \code{$WRMLT} and \code{$WRMGT}, each a 3x4 matrix.
#' @export
do_weighted_cif_analysis <- function(prepped, tau){
  # Pull subsets named (Treatment/Control).(death/discharge)
  Treatment_death    <- prepped$Treatment.death
  Control_death      <- prepped$Control.death
  Treatment_discharge <- prepped$Treatment.discharge
  Control_discharge   <- prepped$Control.discharge

  # Weighted RMLT => eta=2, Weighted RMGT => eta=1
  WRMLT <- table_weighted(Treatment_death, Control_death, eta=2, tau)
  WRMGT <- table_weighted(Treatment_discharge, Control_discharge, eta=1, tau)

  list(
    WRMLT = WRMLT,
    WRMGT = WRMGT
  )
}
