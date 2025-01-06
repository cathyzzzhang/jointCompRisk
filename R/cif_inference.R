#' Compute AUC under a cumulative incidence function
#'
#' @description
#' Computes the area under the curve (AUC) from a \code{survfit} output for a specific
#' cumulative incidence function (CIF) between a given time window [\code{start.time}, \code{tau}].
#'
#' @param table A two-column matrix/data frame where the first column is time and the second column is the CIF estimate.
#' @param start.time The lower time bound.
#' @param tau The upper time bound.
#'
#' @return A numeric value of the AUC.
#'
#' @export
auc.func <- function(table, start.time, tau) {
  loc1 <- which(table[,1] >= start.time)[1]
  loc2 <- max(which(table[,1] <= tau))
  length_interval <- loc2 - loc1 + 1
  auc <- 0
  if (length_interval == 1) {
    auc <- auc + table[loc2,2]*(tau - table[loc2,1])
  } else {
    for(i in seq_len(length_interval - 1)) {
      auc <- auc + table[loc1 + i - 1, 2] * (table[loc1 + i, 1] - table[loc1 + i - 1, 1])
    }
    auc <- auc + table[loc2, 2] * (tau - table[loc2, 1])
  }
  return(auc)
}

#' Compute estimate and variance for joint CIF measure
#'
#' @description
#' Calculates the value of the measure \eqn{\Psi} (defined in Section 4 of the paper) and
#' its standard error for competing risks data.
#'
#' @param data.w A data frame containing time to event (\code{etime}), censoring indicator (\code{estatus}),
#'   event type (\code{etype2}), and treatment group indicator (\code{trt} if relevant).
#' @param tau The time window upper bound.
#' @param a,b,c Numeric coefficients for the measure \eqn{\Psi = a + b \cdot AUC(\text{CIF}_1) + c \cdot AUC(\text{CIF}_2)}.
#'
#' @return A numeric vector of length two:
#'   \itemize{
#'     \item \code{[1]} = estimated value of \eqn{\Psi}
#'     \item \code{[2]} = standard error
#'   }
#'
#' @details
#' Implementation follows equations from Web Appendix D3 in the referenced article.
#'
#' @export
auc.var.joint <- function(data.w, tau, a, b, c) {
  fit <- survival::survfit(survival::Surv(etime, estatus) ~ 1, etype=etype2, data=data.w)
  n.risk <- fit$n.risk[,1]
  n.event <- fit$n.event
  # combine event 1 and event 2 for 1st column
  n.event[,1] <- n.event[,2] + n.event[,3]
  survival.matrix <- fit$pstate
  time <- fit$time
  length_idx <- sum(time <= tau)

  var <- 0
  for(i in seq_len(length_idx)) {
    auc.cif1 <- auc.func(cbind(time, survival.matrix[,2]), time[i], tau)
    cif1     <- survival.matrix[i,2]
    auc.cif2 <- auc.func(cbind(time, survival.matrix[,3]), time[i], tau)
    cif2     <- survival.matrix[i,3]
    surv     <- survival.matrix[i,1]

    var <- var + ( b^2*auc.cif1^2 - 2*b^2*(tau - time[i])*cif1*auc.cif1 + b^2*(tau - time[i])^2*cif1^2 +
                     2*b*c*auc.cif1*auc.cif2 - 2*b*c*(tau - time[i])*(cif2*auc.cif1 + cif1*auc.cif2) +
                     2*b*c*(tau - time[i])^2*cif1*cif2 +
                     c^2*auc.cif2^2 - 2*c^2*(tau - time[i])*cif2*auc.cif2 + c^2*(tau - time[i])^2*cif2^2
    ) * n.event[i,1] / n.risk[i]^2

    var <- var + ( b^2*(tau - time[i])^2*(surv^2 + 2*surv*cif1) - 2*b^2*(tau - time[i])*surv*auc.cif1 -
                     2*b*c*(tau - time[i])*surv*auc.cif2 + 2*b*c*(tau - time[i])^2*surv*cif2
    ) * n.event[i,2] / n.risk[i]^2

    var <- var + ( c^2*(tau - time[i])^2*(surv^2 + 2*surv*cif2) - 2*c^2*(tau - time[i])*surv*auc.cif2 -
                     2*b*c*(tau - time[i])*surv*auc.cif1 + 2*b*c*(tau - time[i])^2*surv*cif1
    ) * n.event[i,3] / n.risk[i]^2
  }

  est <- a + b * auc.func(cbind(time, survival.matrix[,2]), 1, tau) +
    c * auc.func(cbind(time, survival.matrix[,3]), 1, tau)
  se <- sqrt(var)
  return(c(est, se))
}

#' Create summary table for the difference of joint CIF measures
#'
#' @description
#' Creates a 3x4 matrix that contains the estimates and confidence intervals of \eqn{\Psi}
#' for two groups (e.g., treatment vs. control) and their difference.
#'
#' @param data1 A data frame for group 1 (treatment) with columns \code{etime}, \code{estatus}, \code{etype2}.
#' @param data2 A data frame for group 2 (control) with columns \code{etime}, \code{estatus}, \code{etype2}.
#' @param tau The time window upper bound.
#' @param a,b,c Numeric coefficients defining the \eqn{\Psi} measure.
#'
#' @return A 3x4 matrix with rows:
#'   \enumerate{
#'     \item \eqn{\Psi} in group 1
#'     \item \eqn{\Psi} in group 2
#'     \item Difference (\eqn{| \Psi_1 - \Psi_2 |})
#'   }
#'
#' @export
table1 <- function(data1, data2, tau, a, b, c){
  diff1 <- auc.var.joint(data1, tau, a, b, c)
  diff2 <- auc.var.joint(data2, tau, a, b, c)
  table_out <- NULL

  # line1 for \Psi in group 1
  line1 <- c(diff1[1], diff1[1] - 1.96*diff1[2], diff1[1] + 1.96*diff1[2], 0)
  # line2 for \Psi in group 2
  line2 <- c(diff2[1], diff2[1] - 1.96*diff2[2], diff2[1] + 1.96*diff2[2], 0)

  # line3 for difference
  psi_diff <- abs(diff1[1] - diff2[1])
  psi_se   <- sqrt(diff1[2]^2 + diff2[2]^2)
  line3 <- c(psi_diff,
             psi_diff - 1.96 * psi_se,
             psi_diff + 1.96 * psi_se,
             2*(1 - pnorm(psi_diff / psi_se, 0, 1)) )

  table_out <- rbind(line1, line2, line3)
  return(table_out)
}
