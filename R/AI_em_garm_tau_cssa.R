library(Rssa)

# tau
angle.fun <- function(P,Q){
  temp <- function(x,y){
    x[which(y!=0)]
  }
  
  angle <- function(P1,P2,Q1,Q2){
    t <- acos((P1*P2 + Q1*Q2)/sqrt(P1^2+Q1^2)/sqrt(P2^2+Q2^2))
    t[!is.nan(t)]
  }
  var(angle(P[-length(P)],P[-1],Q[-length(Q)],Q[-1]))
}

# The calculation of the optimal tau and t for two vectors
optim.turn <- function(P, Q, fun){
  f_optim <- function(t) {
    fun(Re(P), Re(Q*exp(1.0i*2*pi*t)))+fun(Im(P), Im(Q*exp(1.0i*2*pi*t)))
  }
  
  opt <- optimize(f =  f_optim, interval = c(0,0.5))
  list(tau = opt$objective, t_min = opt$minimum)
}

AI_em_garm_tau_cssa <- function(s, threshold = 0.01, numcomp1=0, d1=TRUE, d2=TRUE, 
                                numcomp2=0, index, all.pairs=FALSE, f=angle.fun){
  if (missing(index)) {
    index <- 1:nu(s)
  }
  
  idx_itog_1 <- numeric(0)
  idx_itog_2 <- numeric(0)
  
  # d = 1
  if (d1){
    if (numcomp1 > 0) {  
      tau <- numeric(length(index))
      for (i in seq_along(index)){
        tau[i]  <- f(Re(s$U[,index[i]]),Im(s$U[,index[i]]))
      }
      tau_sort <- sort(tau)[(1:numcomp1)]
      idx_itog_1 <- index[which(tau %in% tau_sort)]
      index <- index[-which(index %in% idx_itog_1)]
    }
    else {
      tau <- numeric(length(index))
      for (i in seq_along(index)){
        tau[i]  <- f(Re(s$U[,index[i]]),Im(s$U[,index[i]]))
      }
      idx_itog_1 <- index[which(tau <= threshold)]
      index <- index[-which(index %in% idx_itog_1)]
    }
  }
  
  # d = 2
  if (d2){
    if (numcomp2 > 0) {  
      # if (!all.pairs){
      tau <- numeric(length(index)-1)
      for (i in seq_along(index)[-length(index)]){
        tau[i]  <- optim.turn(s$U[,index[i]], s$U[,index[i+1]], f)$tau
      }
      tau_sort <- sort(tau)[1:numcomp2]
      idx_itog_2 <- index[which(tau %in% tau_sort)]
      # }
      # 
      # else {
      #   idx_itog_2 <- NULL
      # }
    }
    else {
      # if (!all.pairs){
      tau <- numeric(length(index)-1)
      for (i in seq_along(index)[-length(index)]){
        tau[i]  <- optim.turn(s$U[,index[i]], s$U[,index[i+1]], f)$tau
      }
      idx_itog_2 <- index[which(tau <= threshold)]
      # }
      # 
      # else {
      #   idx_itog_2 <- NULL
      # }
    }
  }
  
  unique(sort(c(c(idx_itog_2, idx_itog_2 + 1), idx_itog_1)))
}