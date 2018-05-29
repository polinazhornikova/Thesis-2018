general.grouping.auto.freq.cssa <- function(x, groups,  s_0 = 1, rho_0 = 0.9,...){
  L <- x$window
  n <- nu(x)
  max_k <- length((0:(L %/% 2)) / L)
  
  if (missing(groups))
    groups <- as.list(1:min(nsigma(x),n))
  
  groups <- sort(unique(unlist(groups)))
  
  Fs <- x$U[, groups, drop = FALSE]
  
  #normalization
  Fs_re_n <- apply(Re(Fs), 2, function(v) v/(sqrt(sum(v^2))))
  Fs_im_n <- apply(Im(Fs), 2, function(v) v/(sqrt(sum(v^2))))
  
  # periodogram
  pgs_re <- pgram_1d(Fs_re_n)
  pgs_im <- pgram_1d(Fs_im_n)
  pgs_re$spec <- pgs_re$spec / L
  pgs_im$spec <- pgs_im$spec / L
  
  ### part one
  max_pgram_re <- apply(pgs_re$spec, 2, function(x) pgs_re$freq[which.max(x)])
  max_pgram_im <- apply(pgs_im$spec, 2, function(x) pgs_im$freq[which.max(x)])
  
  I_1 <- groups[L * abs(max_pgram_re - max_pgram_im) <= s_0]

  ### part two
  if (length(I_1) != 0){
    rho_I_1 <- (pgs_re$spec[,I_1]  + pgs_im$spec[,I_1])/2
    I_1_final <- I_1[apply(rho_I_1, 2, function(x) max(x[1:(max_k-1)] + x[2:max_k])) >= rho_0 ]
  }
  else {I_1_final <- numeric(0)}
  
  ### final
  unique(sort(I_1_final))
}


