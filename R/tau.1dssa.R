# tau
angle.fun <- function(P,Q){
  angle <- function(P1,P2,Q1,Q2){
    acos((P1*P2 + Q1*Q2)/sqrt(P1^2+Q1^2)/sqrt(P2^2+Q2^2))
  }
  var(angle(P[-length(P)],P[-1],Q[-length(Q)],Q[-1]))
}

draft.grouping.auto.tau.1dssa <- function(x, threshold = 0.01, numcomp=0, groups,...){
  if (missing(groups)) {
    groups <- 1:nu(x)
  }
  
  tau <- numeric(length(groups) - 1)
  tau_prev <- 1000
  for (j in seq_along(groups)[-length(groups)]){
    tau[j] <- angle.fun(x$U[,groups[j]], s$U[,groups[j+1]])
    if (tau[j] > tau_prev) {
      tau[j] <- 1000
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
  
  unique(sort(c(idx_final, idx_final + 1)))
}