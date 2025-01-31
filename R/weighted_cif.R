#' @keywords internal
I <- function(x) as.numeric(x)

w_surv <- function(data,t,eta,m){
  n <- nrow(data)
  # unique observed survival time
  surv_time_all <- sort(unique(data$D_time),decreasing = F)
  # observed survival time < m
  surv_time <- surv_time_all[surv_time_all<=m]
  # generate empirical estimate for E[I(U<=t,eta=j,X<=m,\Delta=1)]
  K_t <- rep(0,length(surv_time))
  for(i in 1:length(surv_time)){
    K_t[i] <- sum(I(data$wU<=t & data$etype==eta & data$D_status==1 & data$D_time<=surv_time[i]))/n
  }
  K_t <- c(0,K_t)
  fit <- survfit(Surv(D_time,D_status) ~ 1, data=data)
  fit.time <- fit$time
  fit.surv <- fit$surv
  cif <- 0
  for(i in 1:length(surv_time)){
    # empirical estimate for S_X
    Sx <- sum(I(data$D_time>=surv_time[i]))/n
    # overall survival function estimate
    ST <- fit.surv[max(which(fit.time<=surv_time[i]))]
    # sum up in the integral
    cif <- cif+ST*(K_t[i+1]-K_t[i])/Sx
  }
  return(cif)
}

# function for calculating the aj term in Theorem 2
var_est <- function(data, t, eta) {
  # Generate competing risks output
  fit.c <- survfit(Surv(D_time, 1 - D_status) ~ 1, data = data)
  fit.ctime <- fit.c$time
  fit.csurv <- fit.c$surv
  cumhaz <- c(0, fit.c$cumhaz)

  # Number of observations
  n <- nrow(data)
  surv_time <- fit.c$time

  # Compute CIF
  cif <- w_surv(data, t, eta, max(data$D_time) + 1)
  exp_u <- rep(0, length(surv_time))
  for (k in seq_along(surv_time)) {
    exp_u[k] <- w_surv(data, t, eta, surv_time[k])
  }
  aj <- rep(0, n)

  for (i in 1:n) {
    int1 <- 0
    int2 <- 0

    # Loop over survival times
    for (j in seq_along(surv_time)) {
      Sx <- sum(data$D_time >= surv_time[j]) / n
      Sx <- max(Sx, 1e-6)  # Avoid division by zero

      # Cumulative hazard difference with finite check
      delta_cumhaz <- cumhaz[j + 1] - cumhaz[j]
      delta_cumhaz <- ifelse(is.finite(delta_cumhaz), delta_cumhaz, 0)

      term1_num <- sum(data$D_time[i] == surv_time[j] & data$D_status[i] == 0) -
        sum(data$D_time[i] >= surv_time[j]) * delta_cumhaz
      int1 <- int1 + term1_num / Sx

      term2_num <- exp_u[j] * term1_num
      int2 <- int2 + term2_num / Sx
    }

    # Validity checks for row data
    wU_valid <- !is.na(data$wU[i]) && is.finite(data$wU[i])
    etype_valid <- !is.na(data$etype[i]) && is.finite(data$etype[i])
    D_status_valid <- !is.na(data$D_status[i]) && is.finite(data$D_status[i])
    D_time_valid <- !is.na(data$D_time[i]) && is.finite(data$D_time[i])

    # Logical conditions
    cond1 <- ifelse(!is.na(data$wU[i] <= t), data$wU[i] <= t, FALSE)
    cond2 <- ifelse(!is.na(data$etype[i] == eta), data$etype[i] == eta, FALSE)
    cond3 <- ifelse(!is.na(data$D_status[i] == 1), data$D_status[i] == 1, FALSE)

    if (wU_valid && etype_valid && D_status_valid && D_time_valid && cond1 && cond2 && cond3) {
      # Handle survival probability
      fit_index <- which(fit.ctime <= data$D_time[i])
      if (length(fit_index) > 0) {
        survival_prob <- fit.csurv[max(fit_index)]
        survival_prob <- max(survival_prob, 1e-6)  # Cap small probabilities
      } else {
        survival_prob <- 1  # Default to 1 if no valid index
      }

      # Compute aj value with finite check
      aj_value <- 1 / survival_prob - cif + cif * int1 - int2
      aj[i] <- ifelse(is.finite(aj_value), aj_value, 0)

    } else {
      aj_value <- -cif + cif * int1 - int2
      aj[i] <- ifelse(is.finite(aj_value), aj_value, 0)
    }
  }
  return(aj)
}

# function for calculating the AUC for weighted CIF
auc.ws <- function(data,eta,tau){
  # idea is the same above
  ws.time <- sort(unique(data$wU),decreasing = F)
  ws.est <- rep(0,length(ws.time))
  for(i in 1:length(ws.time)){
    ws.est[i] <- w_surv(data,ws.time[i],eta,max(data$D_time)+1)
  }
  table <- cbind(ws.time,ws.est)
  loc <- max(which(ws.time<=tau))
  auc=0
  if(loc==-Inf){ auc <- 0 }
  else{
    if(loc==1){ auc <- auc+table[loc,2]*(tau-table[loc,1]) }
    else{
      for(i in 1:(loc-1)){
        auc <- auc+table[i,2]*(table[i+1,1]-table[i,1])
      }
      auc <- auc+table[loc,2]*(tau-table[loc,1])
    }
  }
  return(as.numeric(auc))
}

# function for standard error estimates for auc.ws
sd.auc.ws <- function(data,eta,tau){
  ws.time <- sort(unique(data$wU),decreasing = F)
  loc <- max(which(ws.time<=tau))
  aj.table <- matrix(0,nrow=nrow(data),ncol=loc+1)
  for(i in 1:(loc+1)){
    aj.table[,i] <- var_est(data,ws.time[i],eta)
  }
  # calculate the variance region by region
  # details found in Web Appendix D3
  var <- 0
  if(loc==-Inf) { var <- 0 }
  else{
    if(loc==1){
      var <- 2*(tau-ws.time[1])*ws.time[1]*cov(aj.table[,1],aj.table[,2])+
        (tau-ws.time[1])^2*var(aj.table[,2])
    }
    else{
      for(j in 2:loc){
        var <- var+2*ws.time[1]*(ws.time[j]-ws.time[j-1])*cov(aj.table[,1],aj.table[,j])
      }
      var <- var+2*ws.time[1]*(tau-ws.time[loc])*cov(aj.table[,1],aj.table[,loc+1])
      for(i in 2:loc){
        for(j in 2:loc){
          var <- var+(ws.time[i]-ws.time[i-1])*(ws.time[j]-ws.time[j-1])*cov(aj.table[,i],aj.table[,j])
        }
        var <- var+(ws.time[i]-ws.time[i-1])*(tau-ws.time[loc])*cov(aj.table[,i],aj.table[,loc+1])
      }
      var <- var+(tau-ws.time[loc])^2*var(aj.table[,loc+1])
    }
  }
  return(sqrt(var/nrow(data)))
}

# create tables for output
table_weighted <- function(data1,data2,eta,tau){
  table <- matrix(0,3,4)
  sd1 <- sd.auc.ws(data1,eta,tau)
  sd2 <- sd.auc.ws(data2,eta,tau)
  sd.com <- sqrt(sd1^2+sd2^2)
  # WRMLT/WRMGT in data1 and its confidence interval
  table[1,1] <- auc.ws(data1,eta,tau)
  table[1,2] <- table[1,1]-1.96*sd1
  table[1,3] <- table[1,1]+1.96*sd1
  # WRMLT/WRMGT in data2 and its confidence interval
  table[2,1] <- auc.ws(data2,eta,tau)
  table[2,2] <- table[2,1]-1.96*sd2
  table[2,3] <- table[2,1]+1.96*sd2
  # WRMLT/WRMGT difference between two groups and its confidence interval
  table[3,1] <- abs(table[1,1]-table[2,1])
  table[3,2] <- table[3,1]-1.96*sd.com
  table[3,3] <- table[3,1]+1.96*sd.com
  table[3,4] <- 2*(1-pnorm(table[3,1]/sd.com,0,1))
  return(table)
}
