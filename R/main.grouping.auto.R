library(Rssa)

# common functions
pgram_1d <- function(x) {
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

pgram_2d <- function(x,n) {
  
  N <- nrow(x[[1]])
  M <- ncol(x[[1]])
  
  X <- list()
  
  for (i in 1:n){
    NN <- dim(x[[i]])
    shift.exp <- exp(2i * pi * floor(NN/2) / NN)
    shift1 <- shift.exp[1]^(0:(NN[1] - 1))
    shift2 <- shift.exp[2]^(0:(NN[2] - 1))
    X[[i]] <- Mod(t(mvfft(t(mvfft(outer(shift1, shift2) * x[[i]])))))
  }
  
  spec <- list()
  
  for (i in 1:n){
    spec[[i]] <- X[[i]]
  }
  
  freq1 <- seq(-0.5, 0.5, length.out = N) 
  freq2 <- seq(-0.5, 0.5, length.out = M)
  
  list(spec = spec, freq1 = freq1, freq2 = freq2)
}

# tau
angle.fun <- function(P,Q){
  angle <- function(P1,P2,Q1,Q2){
    t <- acos((P1*P2 + Q1*Q2)/sqrt(P1^2+Q1^2)/sqrt(P2^2+Q2^2))
    t[!is.nan(t)]
  }
  vv <- var(angle(P[-length(P)],P[-1],Q[-length(Q)],Q[-1]))
  mm <- mean(angle(P[-length(P)],P[-1],Q[-length(Q)],Q[-1]))
  vv/min(1,mm^2)
}

source('pair.freq.1dssa.R')
source('pair.freq.cssa.R')
source('low.freq.2dssa.R')
source('tau.1dssa.R')
source('tau.cssa.R')
source('low.freq.cssa.R')
source('low.freq.mssa.R')


draft.grouping.auto <- function(x, ...,
                          grouping.method = c("pair.freq.1dssa", "tau.1dssa","tau.cssa",
                                              "low.freq.2dssa","low.freq.cssa","low.freq.mssa",
                                              "pair.freq.cssa")) {
  switch(match.arg(grouping.method),
         pair.freq.1dssa = draft.grouping.auto.pair.freq.1dssa(x, ...),
         tau.1dssa  = draft.grouping.auto.tau.1dssa(x, ...),
         tau.cssa  = draft.grouping.auto.tau.cssa(x, ...),
         low.freq.2dssa  = draft.grouping.auto.low.freq.2dssa(x, ...),
         low.freq.cssa  = draft.grouping.auto.low.freq.cssa(x, ...),
         low.freq.mssa  = draft.grouping.auto.low.freq.mssa(x, ...),
         pair.freq.cssa  = draft.grouping.auto.pair.freq.cssa(x, ...))
}
