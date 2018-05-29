general.grouping.auto.tau.mssa <- function(x, threshold = 0.01, numcomp=0, groups,
                                           base = c("eigen", "factor"),...){
  if (missing(groups)) {
    groups <- 1:nu(x)
  }
  
  base <- match.arg(base)
  
  if (identical(base, "eigen")) {
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
  }
  else if (identical(base, "factor")) {
    Fs <- calc.v(x, groups, ...)
    Fs_list <- list()
    K <- cumsum(c(1,x$length - x$window + 1))
    for (i in 1:(length(K)-1)){
      Fs_list[[i]] <- Fs[K[i]:(K[i+1] -1),]
    }
    K_main <- x$length - x$window + 1
    max_k <- 1:length(K_main)
    for (i in (1:length(max_k))){
      max_k[i] <- length((0:(K_main[i] %/% 2)) / K_main[i])
    }
    
    tau_main <- list()
    for (i in (1:length(Fs_list))){
      V <- Fs_list[[i]]
      tau <- numeric(length(groups) - 1)
      tau_prev <- 10000
      for (j in seq_along(groups)[-length(groups)]){
        tau[j] <- angle.fun(V[,groups[j]], V[,groups[j+1]])
        names(tau)[j] <- paste0(j, ", " ,j+1)
        if (tau[j] > tau_prev) {
          tau[j] <- 10000
        }
        tau_prev <- tau[j]
      }
      tau_main[[i]] <- tau
    }
    
    get.element <- function(x, i){
      unlist(lapply(x, function(y) y[j]))
    }
    
    tau_fin <- numeric(length(tau_main[[1]]))
    for (j in 1:length(tau_fin)){
      tau_fin[j] <- min(get.element(tau_main,j))
    }
    
    if (numcomp > 0) {
      tau_sort <- sort(tau_fin)[1:numcomp]
      idx_final <- groups[which(tau_fin %in% tau_sort)]
    }
    else {
      idx_final <- groups[which(tau_fin <= threshold)]
    }
  }
  
  
  return <- list(idx=unique(sort(c(idx_final, idx_final + 1))), tau = tau)
}