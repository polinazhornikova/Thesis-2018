pgram1 <- function(x) {
  if (!is.matrix(x)) x <- as.matrix(x)
  stopifnot(all(is.finite(x)))
  
  X <- mvfft(x)
  n <- nrow(x)
  
  N <- n %/% 2 + 1
  spec <- abs(X[seq_len(N),, drop = FALSE])^2
  
  if (n %% 2 == 0) {
    if (N > 2) spec[2:(N-1), ] <- 2 * spec[2:(N-1), ]
  } else {
    if (N >= 2) spec[2:N, ] <- 2 * spec[2:N, ]
  }
  
  freq <- seq(0, 1, length.out = n + 1)[seq_len(N)]
  
  cumspecfuns <- lapply(seq_len(ncol(x)),
                        function(j)
                          approxfun(c(0, freq[-N] + 1/(2*n), 0.5),
                                    c(0, cumsum(spec[, j])),
                                    rule = 2))
  
  list(spec = spec, freq = freq, cumspecfuns = cumspecfuns)
}

# returns a list of two vectors I_1 and I_2
# I1 are the indices of the components with omega = 0.5
# I2 are the indices of other components (with omega != 0.5)

draft.grouping.auto.pair.freq.1dssa <- function(x, groups,  s_0 = 1, rho_0 = 0.9,...){
  L <- x$window
  n <- nu(x)
  max_k <- length((0:(L %/% 2)) / L)
  
  if (missing(groups))
    groups <- as.list(1:min(nsigma(x),n))
  
  groups <- sort(unique(unlist(groups)))
  groups2 <- groups[-length(groups)]
  
  Fs <- x$U[, groups, drop = FALSE]
  # periodogram
  pgs <- pgram1(Fs)
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
  list(I_1 = unique(sort(I_1_final)), I_2 = unique(sort(I_2_final)))
}


