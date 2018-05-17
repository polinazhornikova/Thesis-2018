library(Rssa)

# tau
angle.fun <- function(P,Q){
  angle <- function(P1,P2,Q1,Q2){
    acos((P1*P2 + Q1*Q2)/sqrt(P1^2+Q1^2)/sqrt(P2^2+Q2^2))
  }
  var(angle(P[-length(P)],P[-1],Q[-length(Q)],Q[-1]))
}

AI_em_garm_tau_1dssa <- function(s, threshold = 0.01, numcomp=0, index){
  if (missing(index)) {
    index <- 1:s$window
  }

  tau <- numeric(length(index) - 1)
  tau_prev <- 1000
  for (j in seq_along(index)[-length(index)]){
    tau[j] <- angle.fun(s$U[,index[j]], s$U[,index[j+1]])
    if (tau[j] > tau_prev) {
      tau[j] <- 1000
    }
    tau_prev <- tau[j]
  }
  
  if (numcomp > 0) {
    tau_sort <- sort(tau)[1:numcomp]
    idx_itog <- index[which(tau %in% tau_sort)]
  }
  else {
    idx_itog <- index[which(tau <= threshold)]
  }
  
  unique(sort(c(idx_itog, idx_itog + 1)))
}