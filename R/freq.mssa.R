general.grouping.auto.freq.mssa <- function(x, groups, s_0 = 1, rho_0 = 0.9,  
                                               base = c("eigen", "factor"),...){
  
  base <- match.arg(base)
  
  L <- x$window
  n <- nu(x)
  max_k <- length((0:(L %/% 2)) / L)
  
  if (missing(groups))
    groups <- as.list(1:min(nsigma(x),n))
  
  groups <- sort(unique(unlist(groups)))
  groups2 <- groups[-length(groups)]
  
  if (identical(base, "eigen")) {
    # Fs <- .U(x)[, groups, drop = FALSE]
    Fs <- s$U[, groups]
    
    # periodogram
    pgs <- pgram_1d(Fs)
    pgs$spec <- pgs$spec / L
    
    ### part one
    max_pgram <- apply(pgs$spec, 2, function(x) pgs$freq[which.max(x)])
    
    I_1 <- groups2[L * abs(max_pgram[1:(n-1)] - max_pgram[2:n]) <= s_0]
    I_2 <- groups[L * abs(max_pgram - 0.5) <= s_0]
    
    
    ### part two
    if (length(I_1) != 0){
      rho_I_1 <- (pgs$spec[,I_1]  + pgs$spec[,(I_1 + 1)])/2
      I_1_prefinal <- I_1[apply(rho_I_1, 2, function(x) max(x[1:(max_k-1)] + x[2:max_k])) >= rho_0 ]
      
      if (length(I_1_prefinal) != 0){
        I_1_final <- sort(c(I_1_prefinal, I_1_prefinal+1))
      }
      else {I_1_final <- numeric(0)}
    }
    else {I_1_final <- numeric(0)}
    
    if (length(I_2) != 0){
      rho_I_2 <- pgs$spec[,I_2] 
      r1 <- rho_I_2[which(pgs$freq == (L %/% 2-1) / L)]
      r2 <- rho_I_2[which(pgs$freq == L %/% 2 / L)]
      I_2_final <- I_2[r1+r2 >= rho_0]
    }
    else {I_2_final <- numeric(0)}
    
    ### final
    return <- list(I_1 = unique(sort(I_1_final)), I_2 = unique(sort(I_2_final)))
    
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
    
    #normalization
    Fs_n <- lapply(Fs_list, function(z) apply(z,2, function(y) y/(sqrt(sum(y^2)))))
    
    # periodogram
    pgs_all <- list()
    for (i in (1:length(Fs_n))){
      pgs_all[[i]] <- pgram_1d(Fs_n[[i]])
      pgs_all[[i]]$spec <- pgs_all[[i]]$spec / K_main[i]
    }
    
    ### part one
    max_pgram <- lapply(pgs_all, function(pgs) {
      apply(pgs$spec, 2, function(y) pgs$freq[which.max(y)])
    })
    
    I_1_0 <- list()
    for (i in (1:length(max_pgram))){
      pgram <- max_pgram[[i]]
      I_1_0[[i]] <- groups2[K_main[i] * abs(pgram[1:(n-1)] - pgram[2:n]) <= s_0]
    }
    
    I_2_0 <- list()
    for (i in (1:length(max_pgram))){
      pgram <- max_pgram[[i]]
      I_2_0[[i]] <- groups[K_main[i] * abs(pgram - 0.5) <= s_0]
    }
    
    I_1 <- sort(unique(unlist(I_1_0)))
    I_2 <- sort(unique(unlist(I_2_0)))
    
    rho_1 <- list()
    for (i in (1:length(pgs_all))){
      pgs <- pgs_all[[i]]
      if (length(I_1) != 0){
        rho_I_1 <- (pgs$spec[,I_1]  + pgs$spec[,(I_1 + 1)])/2
        rho_1[[i]] <- apply(rho_I_1, 2, function(x) max(x[1:(max_k[i]-1)] + x[2:max_k[i]]))
      }
      else {rho_1[[i]] <- numeric(0)}
    }
    
    get.element <- function(x, i){
      unlist(lapply(x, function(y) y[j]))
    }
    
    rho_fin <- numeric(length(rho_1[[1]]))
    for (j in 1:length(rho_fin)){
      rho_fin[j] <- max(get.element(rho_1,j))
    }
    
    I_1_prefinal <- I_1[rho_fin >= rho_0]
    
    if (length(I_1_prefinal) != 0){
      I_1_final <- sort(c(I_1_prefinal, I_1_prefinal+1))
    }
    else {I_1_final <- numeric(0)}
    
    ## I_2
    if (length(I_2) != 0){
      rho_2 <- list()
      for (i in (1:length(pgs_all))){
        pgs <- pgs_all[[i]]
        rho_I_2 <- pgs$spec[,I_2]
        r1 <- rho_I_2[which(round(pgs$freq,6) == round((K_main[i] %/% 2-1) / K_main[i],6))]
        r2 <- rho_I_2[which(round(pgs$freq,6) == round((K_main[i] %/% 2) / K_main[i],6))]
        rho_2[[i]] <- r1+r2
      }
      
      if (length(rho_2) !=0){
        rho_fin_2 <- numeric(length(rho_2[[1]]))
        for (j in 1:length(rho_fin_2)){
          rho_fin_2[j] <- max(get.element(rho_2,j))
        }
        I_2_final <- I_2[rho_fin_2 >= rho_0]
      }
      else {I_2_final <- numeric(0)}
      
    }
    else {I_2_final <- numeric(0)}
    
    ### final
    return <- list(I_1 = unique(sort(I_1_final)), I_2 = unique(sort(I_2_final)))
  }
}


