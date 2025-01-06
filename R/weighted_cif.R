#' Indicator function
#'
#' A small helper that returns 1 if expression is TRUE, otherwise 0.
#' Used internally for clarity.
#'
#' @param x A logical value.
#' @export
I <- function(x) ifelse(x, 1, 0)

#' Weighted survival function for competing risks
#'
#' @description
#' Estimates \eqn{E[I(U \le t, \eta = j, X \le m)]} given a data set with
#' weighted survival times.
#'
#' @param data A data frame that must include:
#'   \describe{
#'     \item{D_time}{Time to event}
#'     \item{D_status}{Censoring indicator (1 = event, 0 = censored)}
#'     \item{etype}{Event type}
#'     \item{wU}{Weighted survival time}
#'   }
#' @param t A numeric value, the cutoff for \code{wU}.
#' @param eta The event type of interest (1 or 2).
#' @param m The upper bound on \code{D_time}.
#'
#' @return A single numeric value of the weighted survival function at the given cutoff.
#'
#' @export
w_surv <- function(data, t, eta, m){
  n <- nrow(data)
  surv_time_all <- sort(unique(data$D_time), decreasing = FALSE)
  surv_time     <- surv_time_all[surv_time_all <= m]

  # K_t[i] = # of obs with (wU<=t, etype=eta, D_status=1, D_time<=surv_time[i]) / n
  K_t <- numeric(length(surv_time))
  for(i in seq_along(surv_time)){
    K_t[i] <- sum(I(data$wU <= t & data$etype == eta & data$D_status == 1 &
                      data$D_time <= surv_time[i])) / n
  }
  K_t <- c(0, K_t)

  fit <- survival::survfit(survival::Surv(D_time, D_status) ~ 1, data=data)
  fit_time <- fit$time
  fit_surv <- fit$surv
  cif <- 0

  for(i in seq_along(surv_time)){
    # Empirical S_X = proportion of people who have D_time >= surv_time[i]
    Sx <- sum(I(data$D_time >= surv_time[i])) / n
    # Overall survival function at surv_time[i]
    idx <- max(which(fit_time <= surv_time[i]))
    ST <- fit_surv[idx]
    cif <- cif + ST * (K_t[i+1] - K_t[i]) / Sx
  }
  return(cif)
}

#' Influence Function for Weighted Survival
#'
#' @description
#' Calculates the influence function (a_j) for each individual, used in variance estimation.
#'
#' @param data Same structure as \code{\link{w_surv}}.
#' @param t Cutoff for \code{wU}.
#' @param eta Event type of interest.
#'
#' @return A numeric vector of length = number of rows in \code{data}, containing the influence function for each subject.
#'
#' @export
var_est <- function(data, t, eta){
  fit.c <- survival::survfit(survival::Surv(D_time, 1 - D_status) ~ 1, data=data)
  fit.ctime <- fit.c$time
  fit.csurv <- fit.c$surv
  cumhaz <- c(0, fit.c$cumhaz)
  n <- nrow(data)
  surv_time <- fit.ctime

  # CIF at t for X up to the max
  cif <- w_surv(data, t, eta, max(data$D_time) + 1)
  exp_u <- numeric(length(surv_time))
  for(i in seq_along(surv_time)){
    exp_u[i] <- w_surv(data, t, eta, surv_time[i])
  }

  aj <- numeric(n)
  for(i in seq_len(n)){
    int1 <- 0
    int2 <- 0
    for(j in seq_along(surv_time)){
      Sx <- sum(I(data$D_time >= surv_time[j])) / n
      int1 <- int1 + ( sum(I(data$D_time[i] == surv_time[j] & data$D_status[i] == 0)) -
                         sum(I(data$D_time[i] >= surv_time[j])) * (cumhaz[j+1] - cumhaz[j]) ) / Sx
      int2 <- int2 + exp_u[j] *
        ( sum(I(data$D_time[i] == surv_time[j] & data$D_status[i] == 0)) -
            sum(I(data$D_time[i] >= surv_time[j])) * (cumhaz[j+1] - cumhaz[j]) ) / Sx
    }
    if(I(data$wU[i] <= t & data$etype[i] == eta & data$D_status[i] == 1)){
      aj[i] <- 1/fit.csurv[max(which(fit.ctime <= data$D_time[i]))] - cif + cif*int1 - int2
    } else {
      aj[i] <- - cif + cif*int1 - int2
    }
  }
  return(aj)
}

#' AUC for Weighted CIF
#'
#' @description
#' Computes the area under the weighted survival function curve for event type \code{eta}.
#'
#' @param data A data frame with structure described in \code{\link{w_surv}}.
#' @param eta The event type of interest.
#' @param tau The cutoff on the weighted time scale (\code{wU}).
#'
#' @return Numeric value of the weighted AUC.
#'
#' @export
auc.ws <- function(data, eta, tau){
  ws.time <- sort(unique(data$wU), decreasing = FALSE)
  ws.est  <- numeric(length(ws.time))

  for(i in seq_along(ws.time)){
    ws.est[i] <- w_surv(data, ws.time[i], eta, max(data$D_time) + 1)
  }
  table_tmp <- cbind(ws.time, ws.est)
  loc <- suppressWarnings(max(which(ws.time <= tau)))
  if(loc == -Inf) return(0)

  auc_val <- 0
  if(loc == 1){
    auc_val <- auc_val + table_tmp[loc,2]*(tau - table_tmp[loc,1])
  } else {
    for(i in seq_len(loc-1)){
      auc_val <- auc_val + table_tmp[i,2]*(table_tmp[i+1,1] - table_tmp[i,1])
    }
    auc_val <- auc_val + table_tmp[loc,2]*(tau - table_tmp[loc,1])
  }
  return(as.numeric(auc_val))
}

#' Standard Error for \code{auc.ws}
#'
#' @description
#' Computes the standard error of the weighted AUC using influence functions.
#'
#' @param data Data frame with columns \code{D_time}, \code{D_status}, \code{etype}, and \code{wU}.
#' @param eta Event type of interest.
#' @param tau Upper bound on weighted time \code{wU}.
#'
#' @return Standard error (numeric scalar).
#'
#' @export
sd.auc.ws <- function(data, eta, tau){
  ws.time <- sort(unique(data$wU), decreasing = FALSE)
  loc <- suppressWarnings(max(which(ws.time <= tau)))
  aj.table <- matrix(0, nrow=nrow(data), ncol=loc+1)

  # We need the influence function at each ws.time + an extra index
  for(i in seq_len(loc+1)){
    aj.table[,i] <- var_est(data, ws.time[i], eta)
  }

  var_val <- 0
  if(loc == -Inf) {
    var_val <- 0
  } else {
    # Summation region by region (Web Appendix D3 detail)
    if(loc == 1){
      var_val <- 2*(tau - ws.time[1])*ws.time[1]*cov(aj.table[,1], aj.table[,2]) +
        (tau - ws.time[1])^2*var(aj.table[,2])
    } else {
      # region expansions
      for(j in 2:loc){
        var_val <- var_val + 2*ws.time[1]*(ws.time[j] - ws.time[j-1])*cov(aj.table[,1], aj.table[,j])
      }
      var_val <- var_val + 2*ws.time[1]*(tau - ws.time[loc])*cov(aj.table[,1], aj.table[,loc+1])

      for(i in 2:loc){
        for(j in 2:loc){
          var_val <- var_val + (ws.time[i] - ws.time[i-1])*(ws.time[j] - ws.time[j-1])*cov(aj.table[,i], aj.table[,j])
        }
        var_val <- var_val + (ws.time[i] - ws.time[i-1])*(tau - ws.time[loc])*cov(aj.table[,i], aj.table[,loc+1])
      }
      var_val <- var_val + (tau - ws.time[loc])^2*var(aj.table[,loc+1])
    }
  }
  return(sqrt(var_val / nrow(data)))
}

#' Create Table for Weighted CIF AUC
#'
#' @description
#' Summarizes weighted AUC for two groups, including point estimates and 95% confidence intervals,
#' plus the absolute difference between groups and p-value.
#'
#' @param data1 Data frame for group 1.
#' @param data2 Data frame for group 2.
#' @param eta Event type (1 = discharge, 2 = death).
#' @param tau Upper bound on the weighted time scale.
#'
#' @return A 3x4 numeric matrix:
#'   \enumerate{
#'     \item Group 1 estimate, lower CI, upper CI
#'     \item Group 2 estimate, lower CI, upper CI
#'     \item Absolute difference, lower CI, upper CI, p-value
#'   }
#'
#' @export
table_wr <- function(data1, data2, eta, tau){
  out <- matrix(0, 3, 4)
  sd1 <- sd.auc.ws(data1, eta, tau)
  sd2 <- sd.auc.ws(data2, eta, tau)
  sd.com <- sqrt(sd1^2 + sd2^2)

  est1 <- auc.ws(data1, eta, tau)
  est2 <- auc.ws(data2, eta, tau)

  out[1,1] <- est1
  out[1,2] <- est1 - 1.96*sd1
  out[1,3] <- est1 + 1.96*sd1

  out[2,1] <- est2
  out[2,2] <- est2 - 1.96*sd2
  out[2,3] <- est2 + 1.96*sd2

  diff_val <- abs(est1 - est2)
  out[3,1] <- diff_val
  out[3,2] <- diff_val - 1.96*sd.com
  out[3,3] <- diff_val + 1.96*sd.com
  out[3,4] <- 2*(1 - pnorm(diff_val/sd.com, 0, 1))

  return(out)
}

