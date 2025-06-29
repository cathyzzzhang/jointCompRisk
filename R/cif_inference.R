#' @importFrom survival survfit Surv
#' @import survival
#' @importFrom dplyr %>% select all_of rename left_join arrange

auc.func <- function(table, start.time, tau){
  loc1 <- which(table[,1] >= start.time)[1]
  loc2 <- max(which(table[,1] <= tau))
  length_int <- loc2 - loc1 + 1
  auc <- 0
  if(length_int == 1){
    auc <- auc + table[loc2,2] * (tau - table[loc2,1])
  } else {
    for(i in 1:(length_int - 1)){
      auc <- auc + table[loc1 + i - 1,2] * (table[loc1 + i, 1] - table[loc1 + i - 1, 1])
    }
    auc <- auc + table[loc2,2] * (tau - table[loc2,1])
  }
  return(auc)
}

auc.var.joint <- function(data.w, tau, a, b, c){
  fit <- survival::survfit(Surv(etime, estatus) ~ 1, etype=etype2, data=data.w)
  n.risk <- fit$n.risk[,1]
  n.event <- fit$n.event
  # Combine events 2 & 3 into 1st column => total events
  n.event[,1] <- n.event[,2] + n.event[,3]
  survival.matrix <- fit$pstate
  time <- fit$time
  length_idx <- sum(time <= tau)
  var <- 0
  for(i in 1:length_idx){
    auc.cif1 <- auc.func(cbind(time, survival.matrix[,2]), time[i], tau)
    cif1     <- survival.matrix[i,2]
    auc.cif2 <- auc.func(cbind(time, survival.matrix[,3]), time[i], tau)
    cif2     <- survival.matrix[i,3]
    surv     <- survival.matrix[i,1]

    var <- var + (
      b^2*auc.cif1^2 - 2*b^2*(tau - time[i])*cif1*auc.cif1 + b^2*(tau - time[i])^2*cif1^2 +
        2*b*c*auc.cif1*auc.cif2 - 2*b*c*(tau - time[i])*(cif2*auc.cif1 + cif1*auc.cif2) +
        2*b*c*(tau - time[i])^2*cif1*cif2 +
        c^2*auc.cif2^2 - 2*c^2*(tau - time[i])*cif2*auc.cif2 + c^2*(tau - time[i])^2*cif2^2
    ) * n.event[i,1]/n.risk[i]^2

    var <- var + (
      b^2*(tau - time[i])^2*(surv^2 + 2*surv*cif1) -
        2*b^2*(tau - time[i])*surv*auc.cif1 -
        2*b*c*(tau - time[i])*surv*auc.cif2 +
        2*b*c*(tau - time[i])^2*surv*cif2
    ) * n.event[i,2]/n.risk[i]^2

    var <- var + (
      c^2*(tau - time[i])^2*(surv^2 + 2*surv*cif2) -
        2*c^2*(tau - time[i])*surv*auc.cif2 -
        2*b*c*(tau - time[i])*surv*auc.cif1 +
        2*b*c*(tau - time[i])^2*surv*cif1
    ) * n.event[i,3]/n.risk[i]^2
  }
  est <- a + b*auc.func(cbind(time, survival.matrix[,2]), 1, tau) +
    c*auc.func(cbind(time, survival.matrix[,3]), 1, tau)
  return(c(est, sqrt(var)))
}

table1_cif <- function(data1, data2, tau, a, b, c){
  diff1 <- auc.var.joint(data1, tau, a, b, c)
  diff2 <- auc.var.joint(data2, tau, a, b, c)

  line1 <- c(diff1[1], diff1[1] - 1.96*diff1[2], diff1[1] + 1.96*diff1[2], 0)
  line2 <- c(diff2[1], diff2[1] - 1.96*diff2[2], diff2[1] + 1.96*diff2[2], 0)
  psi_diff <- abs(diff1[1] - diff2[1])
  psi_se   <- sqrt(diff1[2]^2 + diff2[2]^2)
  line3 <- c(psi_diff,
             psi_diff - 1.96*psi_se,
             psi_diff + 1.96*psi_se,
             2*(1 - pnorm(psi_diff/psi_se, 0, 1)) )
  out <- rbind(line1, line2, line3)
  rownames(out) <- c("Group1 (trt=1)", "Group2 (trt=0)", "Difference")
  colnames(out) <- c("Estimate", "Lower95", "Upper95", "p-value")
  return(out)
}

