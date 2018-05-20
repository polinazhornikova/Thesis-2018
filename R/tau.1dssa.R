# tau
angle.fun <- function(P,Q){
  angle <- function(P1,P2,Q1,Q2){
    acos((P1*P2 + Q1*Q2)/sqrt(P1^2+Q1^2)/sqrt(P2^2+Q2^2))
  }
  vv <- var(angle(P[-length(P)],P[-1],Q[-length(Q)],Q[-1]))
  mm <- mean(angle(P[-length(P)],P[-1],Q[-length(Q)],Q[-1]))
  # print(print(c(vv, mm, vv/min(1,mm^2))))
  vv/min(1,mm^2)
}

draft.grouping.auto.tau.1dssa <- function(x, threshold = 0.01, numcomp=0, groups,...){
  if (missing(groups)) {
    groups <- 1:nu(x)
  }
  
  tau <- numeric(length(groups) - 1)
  tau_prev <- 10000
  for (j in seq_along(groups)[-length(groups)]){
    tau[j] <- angle.fun(x$U[,groups[j]], s$U[,groups[j+1]])
    names(tau)[j] <- paste0(j, ", " ,j+1)
    if (tau[j] > tau_prev) {
      tau[j] <- 10000
    }
    tau_prev <- tau[j]
  }
  
  if (numcomp > 0) {
    tau_sort <- sort(tau)[1:numcomp]
    idx_final <- groups[which(tau %in% tau_sort)]
  }
  else {
    idx_final <- groups[which(tau <= threshold)]
  }
  
  list(idx=unique(sort(c(idx_final, idx_final + 1))), tau = tau)
}