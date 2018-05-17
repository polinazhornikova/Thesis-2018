SSA <- function(x, L = 100){
  ssa(x, kind = "1d-ssa", svd.method = "svd", L = L)
}

recon <- function(s, dat, trend = 1){
  r <- reconstruct(s, groups = list(trend = trend))
  dat$T <- r$trend
  dat$R <- dat$S - dat$T
  dat
}

recon.groups <- function(s, dat, g){
  r <- reconstruct(s, groups = g)
  dat$T <- r$F1
  dat$R <- dat$S - dat$T
  dat
}

deriv1 <- function(f,x){
  n <- length(f)
  df <- numeric(n)
  for (i in (2:(n-1))){
    df[i] <- (f[i + 1] - f[i-1])/(x[i+1]-x[i-1]) 
  }
  df
}


# IRLS
it.proc.mult <- function(dat, init = c(runif(1), runif(1)), eps = 10^(-8),
                         max_count = 200, regular = TRUE){
  res <- init
  prev_res <- 1:2
  variance <- 2*(res[2]^2) * (1 + 2*res[1])^2 * dat$D2^2 + 
    2*(res[1]^2) * dat$T2^2 + 
    4*res[2] * res[1] * (1 + 3 * res[1]) * dat$D2 * dat$T2
  if (regular){
    dat$variance <- variance + mean(variance)
  }
  else{
    dat$variance <- variance 
  }
  
  dat <- dat[which(dat$variance > 0),]
  count <- 0
  
  while (((abs(res[1] - prev_res[1]) > eps) | (abs(res[2] - prev_res[2]) > eps)) & count < max_count){
    count <- count + 1
    model <- summary(lm(R2 ~ 0 +T2 + D2, weights = 1 / variance, data = dat))
    result <- model$coefficients[,1]
    prev_res <- res
    res[1] <- result[1]
    res[2] <- result[2] / (1 + result[1])
    variance <- 2*(res[2]^2) * (1 + 2*res[1])^2 * (dat$D2)^2 + 
      2*(res[1]^2) * dat$T2^2 +
      4*res[2] * res[1] * (1 + 3 * res[1]) * dat$D2 * dat$T2
    if (regular){
      dat$variance <- variance + mean(variance)
    }
    else{
      dat$variance <- variance 
    }
    dat <- dat[which(dat$variance > 0),]
  }
  
  list(model, dat)
}

it.proc.add <- function(dat, init = c(runif(1), runif(1)), eps = 10^(-8),
                        max_count = 200, regular = FALSE){
  res <- init
  prev_res <- 1:2
  variance <- 2 * (res[2] * dat$D2 + res[1])^2
  if (regular){
    dat$variance <- variance + mean(variance)
  }
  else{
    dat$variance <- variance 
  }
  count <- 0
  
  while (((abs(res[1] - prev_res[1]) > eps) | (abs(res[2] - prev_res[2]) > eps)) & count < max_count){
    count <- count + 1
    model <- summary(lm(R2 ~ D2, weights = 1 / variance, data = dat))
    result <- model$coefficients[,1]
    prev_res <- res
    res <- result
    variance <- 2 * (res[2] * dat$D2 + res[1])^2
    if (regular){
      dat$variance <- variance + mean(variance)
    }
    else{
      dat$variance <- variance 
    }
  }
  
  list(model,dat)
}

lm_mult <- function(dat){
  summary(lm(R2 ~ 0 + T2 + D2, data = dat))
}

lm_add <- function(dat){
  summary(lm(R2 ~ D2, data = dat))
}

bias <- function(true_param, res){
  sqrt((true_param - res)^2)
} 

get_conf_bound_med <- function(x){
  n <- length(x)
  p_est <- approxfun(density(x))
  error <- qnorm(0.975) * 1/(2 * sqrt(n) * p_est(median(x)))
  left <- median(x) - error
  right <- median(x) + error
  if (0 < left){
    return <- left
  }
  else if (0 > right){
    return <- right
  }
  else {
    return <- 0
  }
}

get_conf_bound_mean <- function(x){
  n <- length(x)
  error <- qnorm(0.975) * sd(x) / sqrt(n)
  left <- mean(x) - error
  right <- mean(x) + error
  if (0 < left){
    return <- abs(left)
  }
  else if (0 > right){
    return <- abs(right)
  }
  else {
    return <- 0
  }
}

conf_interval_mean <- function(x){
  error <- qnorm(0.975)*sd(x)/sqrt(length(x))
  left <- mean(x) - error
  right <- mean(x) + error
  c(left, right) 
}

conf_interval_med <- function(x){
  n <- length(x)
  p_est <- approxfun(density(x))
  error <- qnorm(0.975) * 1/(2 * sqrt(n) * p_est(median(x)))
  left <- median(x) - error
  right <- median(x) + error
  c(left, right) 
}

cosxy <- function(x,y){
  x.new <- x[!(is.na(x) | is.na(y))]
  y.new <- y[!(is.na(x) | is.na(y))]
  sum(x.new * y.new) / (sqrt(sum(x.new^2)) * sqrt(sum(y.new^2)))
}
