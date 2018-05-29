source('tau.1dssa.R')
source('pair.freq.1dssa.R')

library(Rssa)
library(lattice)
library(xtable)


###########################

error <- function(x, y, alpha){
  N <- length(x)
  1/N * sum(exp(-2*alpha * (1:N)) * (x - y)^2)
}

#number of simulations
n_rep <- 200

N <- 99
L <- 50
# alpha <- 0
alpha <- 0.02
sigma_l <- c(0.2,0.4,0.6,0.8,1)
omega <- 1/7

# for threshold 
left <- 0
right <- 1
delta <- 0.01
M <- (right - left) / delta

name_sigma <- numeric(length(sigma_l))
i <- 1
for (sigma in sigma_l){
  name_sigma[i] <- paste0("sigma = ",sigma)
  i <- i + 1
}


### tau1 vs pgram 

for (sigma in sigma_l){
  print(sigma)
  
  err_tau1 <- numeric(n_rep)
  err_pgram <- numeric(n_rep)
  
  for (i in (1:n_rep)){
    print(i)
    signal <- exp(alpha * (1:N)) * cos(2 * pi * omega * (1:N))
    noise <- exp(alpha * (1:N)) * sigma * rnorm(N)
    x <- signal + noise
    s <- ssa(x, L = L)
    # plot(s, type = 'paired')
    rec_V <- reconstruct(s, groups = list(1:2))
    x_recon_V <- rec_V$F1
    
    err_tau1_k <- numeric(n_rep)
    err_pgram_k <- numeric(n_rep)
    
    for (k in (1:M)){
      idx_tau1 <- general.grouping.auto.tau.1dssa(s, threshold = delta*k)$idx
      rec_tau1 <- reconstruct(s, groups = list(idx_tau1))
      x_recon_tau1 <- rec_tau1$F1
      
      idx_pgram <- general.grouping.auto.pair.freq.1dssa(s, rho_0 = delta*k)$I_1
      rec_pgram <- reconstruct(s, groups = list(idx_pgram))
      x_recon_pgram <- rec_pgram$F1
      
      err_tau1_k[k] <- error(x_recon_tau1, x_recon_V, alpha = alpha)
      err_pgram_k[k] <- error(x_recon_pgram, x_recon_V, alpha = alpha)
    }
    k_pgram <- max(which.min(err_pgram_k))
    k_tau1 <- min(which.min(err_tau1_k))
    
    rho_0 <- delta * k_pgram
    threshold <- delta * k_tau1
    
    signal <- exp(alpha * (1:N)) * cos(2 * pi * omega * (1:N))
    noise <- exp(alpha * (1:N)) * sigma * rnorm(N)
    x <- signal + noise
    s <- ssa(x, L = L)
    # plot(s, type = 'paired')
    rec_V <- reconstruct(s, groups = list(1:2))
    x_recon_V <- rec_V$F1
    
    idx_tau1 <- general.grouping.auto.tau.1dssa(s, threshold = threshold)$idx
    rec_tau1 <- reconstruct(s, groups = list(idx_tau1))
    x_recon_tau1 <- rec_tau1$F1
    
    idx_pgram <- general.grouping.auto.pair.freq.1dssa(s, rho_0 = rho_0)$I_1
    rec_pgram <- reconstruct(s, groups = list(idx_pgram))
    x_recon_pgram <- rec_pgram$F1
    
    err_tau1[i] <- error(x_recon_V, x_recon_tau1, alpha = alpha)
    err_pgram[i] <- error(x_recon_V, x_recon_pgram, alpha = alpha)
  }
  
  mean(err_tau1)
  mean(err_pgram)
  
  median(err_tau1)
  median(err_pgram)
  
  d <-  data.frame(err_tau1, err_pgram)
  
  write.table(d,paste0('error_tau_pgram\\err_alpha',alpha,'_omega',omega,'_sigma',sigma,'.txt'))
}


results <- data.frame(mean_tau_1 = numeric(length(sigma_l)),
                      mean_pgram = numeric(length(sigma_l)),
                      median_tau_1 =  numeric(length(sigma_l)),
                      median_pgram =  numeric(length(sigma_l)))
rownames(results) <- name_sigma


i <- 1
for (sigma in sigma_l){
  d <- read.table(paste0('error_tau_pgram\\err_alpha',alpha,'_omega',omega,'_sigma',sigma,'.txt'))
  results$mean_tau_1[i] <- mean(d[,1])
  results$median_tau_1[i] <- median(d[,1])
  results$mean_pgram[i] <- mean(d[,2])
  results$median_pgram[i] <- median(d[,2])
  i <- i + 1
}

results

xtable(results, digits = c(0,4,3,0,4))
