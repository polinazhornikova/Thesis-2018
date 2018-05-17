library(Rssa)
library(lattice)
library(latticeExtra)
library(ggplot2)
library(xtable)


angle.fun <- function(P,Q){
  angle <- function(P1,P2,Q1,Q2){
    acos((P1*P2 + Q1*Q2)/sqrt(P1^2+Q1^2)/sqrt(P2^2+Q2^2))
  }
  var(angle(P[-length(P)],P[-1],Q[-length(Q)],Q[-1]))
}

pgram_fun <- function(x) {
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

get_rho <- function(x, groups=(1:2),  s_0 = 1){
  L <- x$window
  n <- nu(x)
  max_k <- length((0:(L %/% 2)) / L)
  
  if (missing(groups)) {groups <- as.list(1:min(nsigma(x),n))}
  
  groups <- sort(unique(unlist(groups)))
  groups2 <- groups[-length(groups)]
  
  Fs <- x$U[, groups, drop = FALSE]
  # periodogram
  pgs <- pgram_fun(Fs)
  pgs$spec <- pgs$spec / n
  
  ### part one
  
  max_pgram <- apply(pgs$spec, 2, function(x) pgs$freq[which.max(x)])
  I_1 <- groups2[L * abs(max_pgram[1:(length(groups)-1)] - max_pgram[2:length(groups)]) <= s_0]
  
  ### part two
  
  if (length(I_1) != 0){
    I_1_temp <- which(groups2 == I_1)
    rho_I_1 <- as.matrix((pgs$spec[,I_1_temp]  + pgs$spec[,I_1_temp + 1])/2)
    rho_final <- apply(rho_I_1, 2, function(x) max(x[1:(max_k-1)] + x[2:max_k]))
  }
  else {
    rho_final <- 0
  }
 
  return(rho_final)
}


n <- 1000

N <- 99
L <- 50
alpha <- 0.01
sigma_l <- c(0.2,0.4,0.6,0.8,1)

name_sigma <- numeric(length(sigma_l))
i <- 1
for (sigma in sigma_l){
  name_sigma[i] <- paste0('sigma = ',sigma)
  i <- i + 1
}


tau1_itog_signal <- data.frame(matrix(nrow = n, ncol = length(sigma_l)))
tau1_itog_noise <- data.frame(matrix(nrow = n, ncol = length(sigma_l)))

pgram_itog_signal <- data.frame(matrix(nrow = n, ncol = length(sigma_l)))
pgram_itog_noise <- data.frame(matrix(nrow = n, ncol = length(sigma_l)))


colnames(pgram_itog_signal) <- colnames(pgram_itog_noise) <- 
  colnames(tau1_itog_noise)  <- colnames(tau1_itog_signal) <- name_sigma


### signal + noise, L*omega --- not integer 
omega <- 1/7
for (sigma in sigma_l){
  name <- paste0('sigma = ',sigma)
  pgram_noise <- pgram_signal <- numeric(n)
  tau1_noise <- tau1_signal <- numeric(n)
  for (i in 1:n){
    signal <- exp(alpha * (1:N)) * cos(2 * pi * omega * (1:N))
    noise <- exp(alpha * (1:N)) * sigma * rnorm(N)
    x <- signal + noise 
    s <- ssa(x, L = L)
    pgram_signal[i] <- get_rho(s,groups=1:2)
    pgram_noise[i] <- get_rho(s,groups=3:4)
    tau1_signal[i] <- angle.fun(s$U[,1], s$U[,2])
    tau1_noise[i] <- angle.fun(s$U[,3], s$U[,4])
  }
  pgram_itog_signal[,name] <- pgram_signal
  pgram_itog_noise[,name] <- pgram_noise
  tau1_itog_signal[,name] <- tau1_signal
  tau1_itog_noise[,name] <- tau1_noise
}


res_signal_tau1 <- res_noise_tau1 <- summary(tau1_itog_signal[,1])
res_signal_pgram <- res_noise_pgram <- summary(pgram_itog_signal[,1])
overlap_tau1 <- overlap_pgram <- numeric(ncol(pgram_itog_signal))

for (j in 1:ncol(pgram_itog_signal)){
  # print(paste0('sigma = ',sigma_l[j]))
  # print('signal')
  # print(summary(pgram_itog_signal[,j]))
  
  res_signal_tau1 <- rbind(res_signal_tau1, summary(tau1_itog_signal[,j]))
  res_noise_tau1 <- rbind(res_noise_tau1, summary(tau1_itog_noise[,j]))
  
  res_signal_pgram <- rbind(res_signal_pgram, summary(pgram_itog_signal[,j]))
  res_noise_pgram <- rbind(res_noise_pgram, summary(pgram_itog_noise[,j]))
  
  # print('noise')
  # print(summary(pgram_itog_noise[,j]))
  
  a <- pgram_itog_signal[,j]
  b <- pgram_itog_noise[,j]
  
  lower <- min(c(a, b)) - 1 
  upper <- max(c(a, b)) + 1
  
  # generate kernel densities
  da <- density(a, from=lower, to=upper)
  db <- density(b, from=lower, to=upper)
  d <- data.frame(x=da$x, a=da$y, b=db$y)
  
  # calculate intersection densities
  d$w <- pmin(d$a, d$b)
  
  # integrate areas under curves
  library(sfsmisc)
  total <- integrate.xy(d$x, d$a) + integrate.xy(d$x, d$b)
  intersection <- integrate.xy(d$x, d$w)
  
  # compute overlap coefficient
  overlap_pgram[j] <- 2 * intersection / total
  
  a <- tau1_itog_signal[,j]
  b <- tau1_itog_noise[,j]
  
  lower <- min(c(a, b)) - 1 
  upper <- max(c(a, b)) + 1
  
  # generate kernel densities
  da <- density(a, from=lower, to=upper)
  db <- density(b, from=lower, to=upper)
  d <- data.frame(x=da$x, a=da$y, b=db$y)
  
  # calculate intersection densities
  d$w <- pmin(d$a, d$b)
  
  # integrate areas under curves
  library(sfsmisc)
  total <- integrate.xy(d$x, d$a) + integrate.xy(d$x, d$b)
  intersection <- integrate.xy(d$x, d$w)
  
  # compute overlap coefficient
  overlap_tau1[j] <- 2 * intersection / total
}

res_signal_tau1 <- res_signal_tau1[-1,]
res_noise_tau1 <- res_noise_tau1[-1,]
res_signal_pgram <- res_signal_pgram[-1,]
res_noise_pgram <- res_noise_pgram[-1,]
rownames(res_signal_tau1) <- rownames(res_noise_tau1) <- name_sigma

xtable(res_signal_tau1, digits = 5)
xtable(res_noise_tau1, digits = 4)
xtable(res_signal_pgram, digits = 3)
xtable(res_noise_pgram, digits = 3)

overlap_tau1 <- as.data.frame(overlap_tau1)
rownames(overlap_tau1) <- name_sigma

overlap_pgram <- as.data.frame(overlap_pgram)
rownames(overlap_pgram) <- name_sigma

xtable(overlap_tau1, digits = 4)
xtable(overlap_pgram, digits = 3)

tau_1 <- overlap_tau1
pgram <- overlap_pgram
xtable(cbind(tau_1,pgram), digits = c(2,4,3))



######################################################

### only noise and only signal, L*omega --- not integer 
omega <- 1/7
for (sigma in sigma_l){
  name <- paste0('sigma = ',sigma)
  pgram_noise <- pgram_signal <- numeric(n)
  tau1_noise <- tau1_signal <- numeric(n)
  for (i in 1:n){
    signal <- exp(alpha * (1:N)) * cos(2 * pi * omega * (1:N))
    noise <- exp(alpha * (1:N)) * sigma * rnorm(N)
    x <- signal + noise 
    s <- ssa(signal, L = L)
    s2 <- ssa(noise, L = L)
    # plot(s, type = 'paired')
    pgram_signal[i] <- get_rho(s,groups=(1:2))
    pgram_noise[i] <- get_rho(s2,groups=(1:2))
    tau1_signal[i] <- angle.fun(s$U[,1], s$U[,2])
    tau1_noise[i] <- angle.fun(s2$U[,1], s2$U[,2])
  }
  pgram_itog_signal[,name] <- pgram_signal
  pgram_itog_noise[,name] <- pgram_noise
  tau1_itog_signal[,name] <- tau1_signal
  tau1_itog_noise[,name] <- tau1_noise
}


res_signal_tau1 <- summary(tau1_itog_signal[,1])
res_noise_tau1 <- summary(tau1_itog_noise[,1])
res_signal_pgram <- summary(pgram_itog_signal[,1])
res_noise_pgram <- summary(pgram_itog_noise[,1])

overlap_tau1 <- overlap_pgram <- numeric(ncol(pgram_itog_signal))

for (j in 1:ncol(pgram_itog_signal)){
  # print(paste0('sigma = ',sigma_l[j]))
  # print('signal')
  # print(summary(pgram_itog_signal[,j]))
  
  res_signal_tau1 <- rbind(res_signal_tau1, summary(tau1_itog_signal[,j]))
  res_noise_tau1 <- rbind(res_noise_tau1, summary(tau1_itog_noise[,j]))
  
  res_signal_pgram <- rbind(res_signal_pgram, summary(pgram_itog_signal[,j]))
  res_noise_pgram <- rbind(res_noise_pgram, summary(pgram_itog_noise[,j]))
  
  # print('noise')
  # print(summary(pgram_itog_noise[,j]))
  
  a <- pgram_itog_signal[,j]
  b <- pgram_itog_noise[,j]
  
  lower <- min(c(a, b)) - 1 
  upper <- max(c(a, b)) + 1
  
  # generate kernel densities
  da <- density(a, from=lower, to=upper)
  db <- density(b, from=lower, to=upper)
  d <- data.frame(x=da$x, a=da$y, b=db$y)
  
  # calculate intersection densities
  d$w <- pmin(d$a, d$b)
  
  # integrate areas under curves
  library(sfsmisc)
  total <- integrate.xy(d$x, d$a) + integrate.xy(d$x, d$b)
  intersection <- integrate.xy(d$x, d$w)
  
  # compute overlap coefficient
  # overlap_pgram[j] <- 2 * intersection / total
  overlap_pgram[j] <- intersection
  
  a <- tau1_itog_signal[,j]
  b <- tau1_itog_noise[,j]
  
  lower <- min(c(a, b)) - 1 
  upper <- max(c(a, b)) + 1
  
  # generate kernel densities
  da <- density(a, from=lower, to=upper)
  db <- density(b, from=lower, to=upper)
  d <- data.frame(x=da$x, a=da$y, b=db$y)
  
  # calculate intersection densities
  d$w <- pmin(d$a, d$b)
  
  # integrate areas under curves
  library(sfsmisc)
  total <- integrate.xy(d$x, d$a) + integrate.xy(d$x, d$b)
  intersection <- integrate.xy(d$x, d$w)
  
  # compute overlap coefficient
  # overlap_tau1[j] <- 2 * intersection / total
  overlap_tau1[j] <- intersection 
}

res_signal_tau1 <- res_signal_tau1[-1,]
res_noise_tau1 <- res_noise_tau1[-1,]
res_signal_pgram <- res_signal_pgram[-1,]
res_noise_pgram <- res_noise_pgram[-1,]
rownames(res_signal_tau1) <- rownames(res_noise_tau1) <- name_sigma

xtable(res_signal_tau1, digits = 5)
xtable(res_noise_tau1, digits = 3)
xtable(res_signal_pgram, digits = 3)
xtable(res_noise_pgram, digits = 3)

overlap_tau1 <- as.data.frame(overlap_tau1)
rownames(overlap_tau1) <- name_sigma

overlap_pgram <- as.data.frame(overlap_pgram)
rownames(overlap_pgram) <- name_sigma

xtable(overlap_tau1, digits = 4)
xtable(overlap_pgram, digits = 2)

tau_1 <- overlap_tau1
pgram <- overlap_pgram
xtable(cbind(tau_1,pgram), digits = c(2,5,2))





