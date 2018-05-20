# tau
angle.fun <- function(P,Q){
  angle <- function(P1,P2,Q1,Q2){
    t <- acos((P1*P2 + Q1*Q2)/sqrt(P1^2+Q1^2)/sqrt(P2^2+Q2^2))
    t[!is.nan(t)]
  }
  vv <- var(angle(P[-length(P)],P[-1],Q[-length(Q)],Q[-1]))
  mm <- mean(angle(P[-length(P)],P[-1],Q[-length(Q)],Q[-1]))
  # print(print(c(vv, mm, vv/min(1,mm^2))))
  vv/min(1,mm^2)
}

# The calculation of the optimal tau and t for two vectors
optim.turn <- function(P, Q, fun){
  f_optim <- function(t) {
    fun(Re(P), Re(Q*exp(1.0i*2*pi*t)))+fun(Im(P), Im(Q*exp(1.0i*2*pi*t)))
  }
  
  opt <- optimize(f =  f_optim, interval = c(0,0.5))
  list(tau = opt$objective, t_min = opt$minimum)
}

draft.grouping.auto.tau.cssa <- function(x, threshold = 0.01, numcomp1=0, numcomp2=0, 
                                         groups, all.pairs=FALSE, f=angle.fun,...){
  if (missing(groups)) {
    groups <- 1:nu(x)
  }
  
  idx_final_1 <- numeric(0)
  idx_final_2 <- numeric(0)
  
  # d = 1
  
  tau_d1 <- numeric(length(groups))
  for (i in seq_along(groups)){
    tau_d1[i]  <- f(Re(x$U[,groups[i]]),Im(x$U[,groups[i]]))
    names(tau_d1)[i] <- paste0(i)
  }
  if (numcomp1 > 0) {
    tau_sort <- sort(tau_d1)[(1:numcomp1)]
    idx_final_1 <- groups[which(tau_d1 %in% tau_sort)]
    # groups <- groups[-which(groups %in% idx_final_1)]
  }
  else {
    idx_final_1 <- groups[which(tau_d1 <= threshold)]
    # groups <- groups[-which(groups %in% idx_final_1)]
  }
  
  
  # d = 2
  
  if (numcomp2 > 0) {  
    if (!all.pairs){
      tau_d2 <- numeric(length(groups)-1)
      tau_prev <- 10000
      for (i in seq_along(groups)[-length(groups)]){
        tau_d2[i]  <- optim.turn(x$U[,groups[i]], x$U[,groups[i+1]], f)$tau
        names(tau_d2)[i] <- paste0(i, ", " ,i+1)
        if (tau_d2[i] > tau_prev) {
          tau_d2[i] <- 10000
        }
        tau_prev <- tau_d2[i]
      }
      tau_sort <- sort(tau_d2)[1:numcomp2]
      idx_final_2_0 <- groups[which(tau_d2 %in% tau_sort)]
      idx_final_2 <- c(idx_final_2_0, idx_final_2_0+1)
    }
    
    else {
      res <- data.frame(i = 0, j = 0, tau = 0)
      for (ii in seq_along(groups))
        for (jj in (ii + seq_len(length(groups) - ii))){
          tau <- optim.turn(x$U[,groups[ii]], x$U[,groups[jj]], f)$tau
          res <- rbind(res, c(ii,jj,tau)) 
        }
      res <- res[-1,]
      tau_d2 <- res$tau
      names(tau_d2) <- paste0(res$i,', ', res$j)
      tau_sort <- sort(res$tau)[1:numcomp2]
      idx_final_2 <- c(res$i[which(res$tau %in% tau_sort)], res$j[which(res$tau %in% tau_sort)])
    }
  }
  else {
    if (!all.pairs){
      tau_d2 <- numeric(length(groups)-1)
      tau_prev <- 1000
      for (i in seq_along(groups)[-length(groups)]){
        tau_d2[i]  <- optim.turn(x$U[,groups[i]], x$U[,groups[i+1]], f)$tau
        if (tau_d2[i] > tau_prev) {
          tau_d2[i] <- 1000
        }
        tau_prev <- tau_d2[i]
      }
      idx_final_2_0 <- groups[which(tau_d2 <= threshold)]
      idx_final_2 <- c(idx_final_2_0, idx_final_2_0 + 1)
    }
    
    else {
      res <- data.frame(i = 0, j = 0, tau = 0)
      for (ii in seq_along(groups))
        for (jj in (ii + seq_len(length(groups) - ii))){
          tau <- optim.turn(x$U[,groups[ii]], x$U[,groups[jj]], f)$tau
          res <- rbind(res, c(ii,jj,tau)) 
        }
      res <- res[-1,]
      tau_d2 <- res$tau
      names(tau_d2) <- paste0(res$i,', ', res$j)
      idx_final_2 <- c(res$i[which(res$tau <= threshold)], res$j[which(res$tau <= threshold)])
    }
  }
  
  idx_final_2 <- idx_final_2[-which(idx_final_2 %in% idx_final_1)]
  
  list(d1_idx = unique(sort(idx_final_1)), d2_idx = unique(sort(idx_final_2)),
       tau_d1 = tau_d1, tau_d2 = tau_d2)
}