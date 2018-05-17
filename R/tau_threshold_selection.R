source('AI_em_garm_tau_1d.R')


conf_interval_mean95 <- function(x){
  error <- qnorm(0.975)*sd(x)/sqrt(length(x))
  left <- mean(x) - error
  right <- mean(x) + error
  c(left, right) 
}

get_conf_bound_mean <- function(x){
  n <- length(x)
  error <- qnorm(0.975) * sd(x) / sqrt(n)
  left <- mean(x) - error
  right <- mean(x) + error
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

n <- 1000

N <- 99
L <- 50
alpha <- 0.01
sigma_list <- c(0, 0.2,0.4,0.6,0.8,1,1.2,1.4)
omega <- 1/7
A <- 1


## 1 гармоника

tau1_q <- numeric(length(sigma_list))
tau1_q2 <- numeric(length(sigma_list))

for (sigma in sigma_list){
  signal <- A*exp(alpha * (1:N)) * cos(2 * pi * omega * (1:N))
  tau1 <- numeric(n)
  tau1_2 <- numeric(n)
  for (k in (1:n)){
    noise <- A*exp(alpha * (1:N)) * sigma * rnorm(N)
    x <- signal + noise 
    # print(plot(x, type = 'l'))
    s <- ssa(x, L = L)
    # print(plot(s, type = 'paired'))
    # norm <- numeric(ncol(s$U) - 1)
   
    # idx_final <- numeric(0)
    # idx_final2 <- numeric(0)
    
    tau1[k] <- angle.fun(s$U[,1], s$U[,2])
  }
  tau1_q[which(sigma_list == sigma)] <- conf_interval_mean95(tau1)[2]
  tau1_q2[which(sigma_list == sigma)] <- conf_interval_mean95(tau1_2)[2]
}

res <- data.frame(quantile95_tau = tau1_q, quantile_tau_sd = tau1_q2, sigma = sigma_list)

xyplot(quantile95_tau + quantile_tau_sd ~ sigma, data=res, type = "l",auto.key = list(points = FALSE, lines = TRUE))
xyplot(quantile95_tau ~ sigma, data=res, type = "l",auto.key = list(points = FALSE, lines = TRUE))



    # r <- reconstruct(s, groups = list(t = (1:2)))
    # norm[1] <- sum((r$t)^2)
    # order_tau1 <- numeric(ncol(s$U) - 1)
    # for (j in (2:(ncol(s$U) - 1))){
    #   
    #   tau1[j] <- angle.fun(fs$U[,j], fs$U[,j+1])
    #   tau1_2[j] <- angle.fun(fs$U[,j], fs$U[,j+1])
    #   if (j %% 2 == 0){
    #     if (tau1[j] < tau1[j-1]){
    #       idx_final <- c(idx_final, j)
    #     }
    #     else {
    #       idx_final <- c(idx_final, j-1)
    #     }
    #     if (tau1_2[j] < tau1_2[j-1]){
    #       idx_final2 <- c(idx_final2, j)
    #     }
    #     else {
    #       idx_final2 <- c(idx_final2, j-1)
    #     }
    #   }
    # }
    # for (i in idx_final[-1]){
    #   if (i == i+1){
    #     if (tau1[i] > tau1[i+1]){
    #       idx_final <- idx_final[-which(idx_final == i)]
    #     }
    #     else {
    #       idx_final <- idx_final[-which(idx_final == i + 1)]
    #     }
    #   }
    # }
    # 
    # for (i in idx_final2[-1]){
    #   if (i == i+1){
    #     if (tau1_2[i] > tau1_2[i+1]){
    #       idx_final2 <- idx_final2[-which(idx_final2 == i)]
    #     }
    #     else {
    #       idx_final2 <- idx_final2[-which(idx_final2 == i + 1)]
    #     }
    #   }

  # tau1_2 <- tau1[idx]
  # tau1_2_sort <- sort(tau1_2)
  # 
  # tau1_sort <- sort(tau1)
  # norm_sort <- sort(norm, decreasing = TRUE)
  # norm2 <- norm[idx_final]
  # tau1_2 <- tau1[idx_final]
  # idx_final_sort <- idx_final[order(tau1_2)]
  # # tau2_sort <- sort(tau2)
  # dat <- data.frame(tau1 = sort(tau1_2), norm = norm2[order(tau1_2)]/max(norm2),
  #                   norm_1 = norm2[order(tau1_2)]/max(norm2),
  #                   norm_2 = (norm2[order(tau1_2)]/max(norm2))^2)
  # x <- seq(0,1,by=0.001)
  # dat2 <- data.frame(x = x, x2 = exp(x^(for_exp))-1)
  # print(ggplot()+
  #         geom_point(data = dat, aes(norm,tau1))+
  #         geom_line(data = dat2, aes(x,x2), col='red')+
  #         labs(title = paste0('sigma = ', sigma)))
  # print(paste0('sigma = ', sigma))
  # ii <- which(dat$tau1 < exp((dat$norm)^(for_exp))-1)
  # print('How much to identify')
  # print(sum(dat$tau1 < exp((dat$norm)^(for_exp))-1))
  # print('Which components')
  # print(idx_final_sort[ii])
}



## 2 гармоники

omega2 <- 1/5

for (sigma in sigma_list){
  signal <- A*exp(alpha * (1:N)) * cos(2 * pi * omega * (1:N)) +
    A*2*exp(alpha * (1:N)) * cos(2 * pi * omega2 * (1:N))
  noise <- A*exp(alpha * (1:N)) * sigma * rnorm(N)
  x <- signal + noise 
  print(plot(x, type = 'l'))
  s <- ssa(x, L = L)
  print(plot(s, type = 'paired'))
  
  norm <- numeric(ncol(s$U) - 1)
  tau1 <- numeric(ncol(s$U) - 1)
  idx_final <- numeric(0)
  tau1[1] <- angle.fun(s$U[,1], s$U[,2])
  r <- reconstruct(s, groups = list(t = (1:2)))
  norm[1] <- sum((r$t)^2)
  # order_tau1 <- numeric(ncol(s$U) - 1)
  for (j in (2:(ncol(s$U) - 1))){
    tau1[j] <- angle.fun(s$U[,j], s$U[,j+1])
    if (j %% 2 == 0){
      if (tau1[j] < tau1[j-1]){
        idx_final <- c(idx_final, j)
      }
      else {
        idx_final <- c(idx_final, j-1)
      }
    }
    for (i in idx_final[-1]){
      if (idx_final[i] == idx_final[i+1]){
        if (tau1[i] > tau1[i+1]){
          idx_final <- idx_final[-which(idx_final == i)]
        }
        else {
          idx_final <- idx_final[-which(idx_final == i + 1)]
        }
      }
    }
    r <- reconstruct(s, groups = list(t = (j:(j+1))))
    norm[j] <- sum((r$t)^2)
    # tau2[j] <- sqrt(sum(s$U[,j]^2 + s$U[,j+1]^2))
    # temp <-  tau1[j]
    # i_order <- 0
    # while (temp < 1){
    #   temp <- temp * 10
    #   i_order <- i_order + 1
    # }
    # order_tau1[j] <- i_order
  }
  
  tau1_sort <- sort(tau1)
  norm_sort <- sort(norm, decreasing = TRUE)
  norm2 <- norm[idx_final]
  tau1_2 <- tau1[idx_final]
  idx_final_sort <- idx_final[order(tau1_2)]
  # tau2_sort <- sort(tau2)
  dat <- data.frame(tau1 = sort(tau1_2), norm = norm2[order(tau1_2)]/max(norm2),
                    norm_1 = norm2[order(tau1_2)]/max(norm2),
                    norm_2 = (norm2[order(tau1_2)]/max(norm2))^2)
  x <- seq(0,1,by=0.001)
  dat2 <- data.frame(x = x, x2 = exp(x^(for_exp))-1)
  print(ggplot()+
          geom_point(data = dat, aes(norm,tau1))+
          geom_line(data = dat2, aes(x,x2), col='red')+
          labs(title = paste0('sigma = ', sigma)))
  print(paste0('sigma = ', sigma))
  ii <- which(dat$tau1 < exp((dat$norm)^(for_exp))-1)
  print('How much to identify')
  print(sum(dat$tau1 < exp((dat$norm)^(for_exp))-1))
  print('Which components')
  print(idx_final_sort[ii])
}







## real data
for_exp <- 1.9

for (k in (c(1:4,6:7))){
  x <- na.omit(AustralianWine[,k])
  plot(x, type = 'l', main = colnames(AustralianWine)[k])
  s <- ssa(x, L = 84)
  r <- reconstruct(s, groups = list(trend = 1))
  x <- x - r$trend
  s <- ssa(x, L = 84)
  print(plot(s, type = 'paired',idx=1:20))
  
  norm <- numeric(ncol(s$U) - 1)
  tau1 <- numeric(ncol(s$U) - 1)
  idx_final <- numeric(0)
  tau1[1] <- angle.fun(s$U[,1], s$U[,2])
  r <- reconstruct(s, groups = list(t = (1:2)))
  norm[1] <- sum((r$t)^2)
  # order_tau1 <- numeric(ncol(s$U) - 1)
  for (j in (2:(ncol(s$U) - 1))){
    tau1[j] <- angle.fun(s$U[,j], s$U[,j+1])
    if (j %% 2 == 0){
      if (tau1[j] < tau1[j-1]){
        idx_final <- c(idx_final, j)
      }
      else {
        idx_final <- c(idx_final, j-1)
      }
    }
    for (i in idx_final[-1]){
      if (idx_final[i] == idx_final[i+1]){
        if (tau1[i] > tau1[i+1]){
          idx_final <- idx_final[-which(idx_final == i)]
        }
        else {
          idx_final <- idx_final[-which(idx_final == i + 1)]
        }
      }
    }
    r <- reconstruct(s, groups = list(t = (j:(j+1))))
    norm[j] <- sum((r$t)^2)
    # tau2[j] <- sqrt(sum(s$U[,j]^2 + s$U[,j+1]^2))
    # temp <-  tau1[j]
    # i_order <- 0
    # while (temp < 1){
    #   temp <- temp * 10
    #   i_order <- i_order + 1
    # }
    # order_tau1[j] <- i_order
  }
  
  tau1_sort <- sort(tau1)
  plot(tau1_sort)
  norm_sort <- sort(norm, decreasing = TRUE)
  norm2 <- norm[idx_final]
  tau1_2 <- tau1[idx_final]
  idx_final_sort <- idx_final[order(tau1_2)]
  # tau2_sort <- sort(tau2)
  dat <- data.frame(tau1 = sort(tau1_2), norm = norm2[order(tau1_2)]/max(norm2),
                    norm_1 = norm2[order(tau1_2)]/max(norm2),
                    norm_2 = (norm2[order(tau1_2)]/max(norm2))^2)
  x <- seq(0,1,by=0.001)
  dat2 <- data.frame(x = x, x2 = exp(x^(for_exp))-1)
  print(ggplot()+
          geom_point(data = dat, aes(norm,tau1))+
          geom_line(data = dat2, aes(x,x2), col='red')+
          labs(title = colnames(AustralianWine)[k]))
  print(colnames(AustralianWine)[k])
  ii <- which(dat$tau1 < exp((dat$norm)^(for_exp))-1)
  print('How much to identify')
  print(sum(dat$tau1 < exp((dat$norm)^(for_exp))-1))
  print('Which components')
  print(idx_final_sort[ii])
}




#### noise

N <- 199
L <- 100

n <- 500
rr <- 24
all_tau1 <- numeric(n*rr)
all_norm <- numeric(n*rr)
alpha <- 0.01
sigma_list <- c(0.2, 0.4, 0.6, 0.8, 1)
# sigma_list <- c(0.2)

omega <- 1/10
A <- 1
signal <- A*exp(alpha * (1:N)) * cos(2 * pi * omega * (1:N))
for_exp <- 5

for (sigma in sigma_list){
  for (i in (1:n)){
    noise <- A*exp(alpha * (1:N)) * sigma * rnorm(N)
    x <- noise + signal
    # print(plot(x, type = 'l'))
    s <- ssa(x, L = L, svd.method = 'svd')
    # print(plot(s, type = 'paired'))
    norm <- numeric(ncol(s$U) - 1)
    tau1 <- numeric(ncol(s$U) - 1)
    idx_final <- numeric(0)
    tau1[1] <- angle.fun(s$U[,1], s$U[,2])
    r <- reconstruct(s, groups = list(t = (1:2)))
    norm[1] <- sum((r$t)^2)
    # order_tau1 <- numeric(ncol(s$U) - 1)
    for (j in (2:(ncol(s$U) - 1))){
      tau1[j] <- angle.fun(s$U[,j], s$U[,j+1])
      if (j %% 2 == 0){
        if (tau1[j] < tau1[j-1]){
          idx_final <- c(idx_final, j)
        }
        else {
          idx_final <- c(idx_final, j-1)
        }
      }
      for (i in idx_final[-1]){
        if (idx_final[i] == idx_final[i+1]){
          if (tau1[i] > tau1[i+1]){
            idx_final <- idx_final[-which(idx_final == i)]
          }
          else {
            idx_final <- idx_final[-which(idx_final == i + 1)]
          }
        }
      }
      r <- reconstruct(s, groups = list(t = (j:(j+1))))
      norm[j] <- sum((r$t)^2)
      # tau2[j] <- sqrt(sum(s$U[,j]^2 + s$U[,j+1]^2))
      # temp <-  tau1[j]
      # i_order <- 0
      # while (temp < 1){
      #   temp <- temp * 10
      #   i_order <- i_order + 1
      # }
      # order_tau1[j] <- i_order
    }
    
    tau1_sort <- sort(tau1)
    norm_sort <- sort(norm, decreasing = TRUE)
    norm2 <- norm[idx_final]
    tau1_2 <- tau1[idx_final]
    idx_final_sort <- idx_final[order(tau1_2)]
    # tau2_sort <- sort(tau2)
    all_tau1[((i-1)*rr+1):(i*rr)] <- sort(tau1_2)
    all_norm[((i-1)*rr+1):(i*rr)] <- norm2[order(tau1_2)]/max(norm2)
  }
  dat <- data.frame(tau1 = all_tau1, norm = all_norm)
  dat <- dat[-which(dat$norm == 1),]
  h <- 10
  step <- max(dat$norm)/h
  cat <- (1:h)*step
  dat$group  <- numeric(nrow(dat))
  # dat$left <- dat$right <- numeric(nrow(dat))
  # dat$left <- 
  # dat$right <- conf_interval_mean(dat$tau1)[2]
  
  for (i in (1:nrow(dat))){
    dat$group[i] <- dat$norm[i] %/% step + 1
    # dat$left[i] <- conf_interval_mean(dat[which(dat$group == dat$group[i]),'tau1'])[1] 
  }
  
  dat$group <-factor(dat$group, levels = c(1:h))
  
  # library(gridExtra)
  
  for_exp <- 3
  max_x <- min(1, log(2)^(1/for_exp))
  x <- seq(0,max_x,by=max_x/1000)
  dat2 <- data.frame(x = x, x2 = exp(x^(for_exp))-1)
  
  dat <- na.omit(dat)
  
  
  print(ggplot(dat,aes(norm,tau1))+
          geom_point(aes(col = group))+
          labs(title = paste0('sigma = ', sigma)))
  # geom_ribbon(aes(ymin=left,ymax=right),alpha=0.3))
  print(ggplot(dat, aes(group, tau1,col = group))+
          # geom_line(data = dat2, aes(x,x2), col='red')+
          labs(title = paste0('sigma = ', sigma))+
          geom_boxplot())
  
  print(paste0('sigma = ', sigma))
  print(table(dat$group))
  # grid.arrange(plot1, plot2, ncol=2)
  
  # print(paste0('sigma = ', sigma))
  # ii <- which(dat$tau1 < exp((dat$norm)^(for_exp))-1)
  # print('How much to identify')
  # print(sum(dat$tau1 < exp((dat$norm)^(for_exp))-1))
  # print('Which components')
  # print(idx_final_sort[ii])
}




#### noise without sin 

N <- 199
L <- 100

n <- 500
rr <- 24
all_tau1 <- numeric(n*rr)
all_norm <- numeric(n*rr)
alpha <- 0.01
sigma_list <- c(0.2, 0.4, 0.6, 0.8, 1)
# sigma_list <- c(0.2)

omega <- 1/10
A <- 1
# signal <- A*exp(alpha * (1:N)) * cos(2 * pi * omega * (1:N))
for_exp <- 5

for (sigma in sigma_list){
  for (i in (1:n)){
    noise <- A*exp(alpha * (1:N)) * sigma * rnorm(N)
    x <- noise 
    # print(plot(x, type = 'l'))
    s <- ssa(x, L = L, svd.method = 'svd')
    # print(plot(s, type = 'paired'))
    norm <- numeric(ncol(s$U) - 1)
    tau1 <- numeric(ncol(s$U) - 1)
    idx_final <- numeric(0)
    tau1[1] <- angle.fun(s$U[,1], s$U[,2])
    r <- reconstruct(s, groups = list(t = (1:2)))
    norm[1] <- sum((r$t)^2)
    # order_tau1 <- numeric(ncol(s$U) - 1)
    for (j in (2:(ncol(s$U) - 1))){
      tau1[j] <- angle.fun(s$U[,j], s$U[,j+1])
      if (j %% 2 == 0){
        if (tau1[j] < tau1[j-1]){
          idx_final <- c(idx_final, j)
        }
        else {
          idx_final <- c(idx_final, j-1)
        }
      }
      for (i in idx_final[-1]){
        if (idx_final[i] == idx_final[i+1]){
          if (tau1[i] > tau1[i+1]){
            idx_final <- idx_final[-which(idx_final == i)]
          }
          else {
            idx_final <- idx_final[-which(idx_final == i + 1)]
          }
        }
      }
      r <- reconstruct(s, groups = list(t = (j:(j+1))))
      norm[j] <- sum((r$t)^2)
      # tau2[j] <- sqrt(sum(s$U[,j]^2 + s$U[,j+1]^2))
      # temp <-  tau1[j]
      # i_order <- 0
      # while (temp < 1){
      #   temp <- temp * 10
      #   i_order <- i_order + 1
      # }
      # order_tau1[j] <- i_order
    }
    
    tau1_sort <- sort(tau1)
    norm_sort <- sort(norm, decreasing = TRUE)
    norm2 <- norm[idx_final]
    tau1_2 <- tau1[idx_final]
    idx_final_sort <- idx_final[order(tau1_2)]
    # tau2_sort <- sort(tau2)
    all_tau1[((i-1)*rr+1):(i*rr)] <- sort(tau1_2)
    all_norm[((i-1)*rr+1):(i*rr)] <- norm2[order(tau1_2)]/max(norm2)
  }
  dat <- data.frame(tau1 = all_tau1, norm = all_norm)
  dat <- dat[-which(dat$norm == 1),]
  h <- 10
  step <- max(dat$norm)/h
  cat <- (1:h)*step
  dat$group  <- numeric(nrow(dat))
  # dat$left <- dat$right <- numeric(nrow(dat))
  # dat$left <- 
  # dat$right <- conf_interval_mean(dat$tau1)[2]
  
  for (i in (1:nrow(dat))){
    dat$group[i] <- dat$norm[i] %/% step + 1
    # dat$left[i] <- conf_interval_mean(dat[which(dat$group == dat$group[i]),'tau1'])[1] 
  }
  
  dat$group <-factor(dat$group, levels = c(1:h))
  
  # library(gridExtra)
  
  for_exp <- 3
  max_x <- min(1, log(2)^(1/for_exp))
  x <- seq(0,max_x,by=max_x/1000)
  dat2 <- data.frame(x = x, x2 = exp(x^(for_exp))-1)
  
  dat <- na.omit(dat)
  
  
  print(ggplot(dat,aes(norm,tau1))+
          geom_point(aes(col = group))+
          labs(title = paste0('sigma = ', sigma)))
  # geom_ribbon(aes(ymin=left,ymax=right),alpha=0.3))
  print(ggplot(dat, aes(group, tau1,col = group))+
          # geom_line(data = dat2, aes(x,x2), col='red')+
          labs(title = paste0('sigma = ', sigma))+
          geom_boxplot())
  
  print(paste0('sigma = ', sigma))
  print(table(dat$group))
  # grid.arrange(plot1, plot2, ncol=2)
  
  # print(paste0('sigma = ', sigma))
  # ii <- which(dat$tau1 < exp((dat$norm)^(for_exp))-1)
  # print('How much to identify')
  # print(sum(dat$tau1 < exp((dat$norm)^(for_exp))-1))
  # print('Which components')
  # print(idx_final_sort[ii])
}



###example 
N <- 1000

omega1 <- 0.12
omega2 <- 0.121
omega3 <- 0.122

omega1 <- 0.1
omega2 <- 1/7
omega3 <- 1/13

omega1 <- 0.12
omega2 <- 0.13
omega3 <- 0.14

s1 <- cos(2 * pi * omega1 * (1:N))
s2 <- cos(2 * pi * omega2 * (1:N))
s3 <- cos(2 * pi * omega3 * (1:N))

signal <- s1+s2+s3
s <- ssa(signal)
plot(s, type = 'paired')

tau1 <- numeric(ncol(s$U) - 1)

for (j in (2:(ncol(s$U) - 1))){
  tau1[j] <- angle.fun(s$U[,j], s$U[,j+1])
}
tau1[1:6]

###

s1 <- exp(alpha * (1:N)) * cos(2 * pi * omega1 * (1:N))
s2 <- exp(alpha * (1:N)) * cos(2 * pi * omega2 * (1:N))
s3 <- exp(alpha * (1:N)) * cos(2 * pi * omega3 * (1:N))

signal <- s1+s2+s3
s <- ssa(signal)
plot(s, type = 'paired')

tau1 <- numeric(ncol(s$U) - 1)

for (j in (2:(ncol(s$U) - 1))){
  tau1[j] <- angle.fun(s$U[,j], s$U[,j+1])
}
tau1[1:6]





######### other method 

## model data
## 2 гармоники

eps <- 0.05

omega2 <- 1/5

for (sigma in sigma_list){
  signal <- A*exp(alpha * (1:N)) * cos(2 * pi * omega * (1:N)) +
    A*2*exp(alpha * (1:N)) * cos(2 * pi * omega2 * (1:N))
  noise <- A*exp(alpha * (1:N)) * sigma * rnorm(N)
  x <- signal + noise 
  print(plot(x, type = 'l'))
  s <- ssa(x, L = L)
  print(plot(s, type = 'paired'))
  
  norm <- numeric(ncol(s$U) - 1)
  tau1 <- numeric(ncol(s$U) - 1)
  idx_final <- numeric(0)
  tau1[1] <- angle.fun(s$U[,1], s$U[,2])
  r <- reconstruct(s, groups = list(t = (1:2)))
  norm[1] <- sum((r$t)^2)
  
  for (j in (2:(ncol(s$U) - 1))){
    r <- reconstruct(s, groups = list(t = (j:(j+1))))
    norm[j] <- sum((r$t)^2)
    
    tau1[j] <- angle.fun(s$U[,j], s$U[,j+1])
    if (j %% 2 == 0){
      if (tau1[j] < tau1[j-1]){
        idx_final <- c(idx_final, j)
      }
      else {
        idx_final <- c(idx_final, j-1)
      }
    }
    # for (i in idx_final[-1]){
    #   if (idx_final[i] == idx_final[i+1]){
    #     if (tau1[i] > tau1[i+1]){
    #       idx_final <- idx_final[-which(idx_final == i)]
    #     }
    #     else {
    #       idx_final <- idx_final[-which(idx_final == i + 1)]
    #     }
    #   }
    # }
    
  }
  
  norm2_2 <- norm/max(norm)
  
  for ( i in idx_final){
    if (norm2_2[i] < eps){
      idx_final <- idx_final[-which(idx_final == i)]
    }
  }
  
  norm2 <- norm[idx_final]
  tau1_2 <- tau1[idx_final]
  
  idx_final_sort <- idx_final[order(tau1_2)]
  
  
  # tau2_sort <- sort(tau2)
  dat <- data.frame(tau1 = sort(tau1_2), norm = norm2[order(tau1_2)]/max(norm2),
                    norm_1 = norm2[order(tau1_2)]/max(norm2),
                    norm_2 = (norm2[order(tau1_2)]/max(norm2))^2)
  
  
  tau1_diff <- numeric(length(dat$tau1)-1)
  for (i in (1:length(tau1_diff))){
    tau1_diff[i] <- dat$tau1[i+1] - dat$tau1[i]
  }
  print(tau1_diff)
  
  x <- seq(0,1,by=0.001)
  dat2 <- data.frame(x = x, x2 = exp(x^(for_exp))-1)
  print(ggplot()+
          geom_point(data = dat, aes(norm,tau1))+
          # geom_line(data = dat2, aes(x,x2), col='red')+
          labs(title = paste0('sigma = ', sigma)))
  print(paste0('sigma = ', sigma))
  plot(dat$tau1)
  print(idx_final_sort)
  plot(dat$tau1)
  plot(dat$norm)
  # ii <- which(dat$tau1 < exp((dat$norm)^(for_exp))-1)
  # print('How much to identify')
  # print(sum(dat$tau1 < exp((dat$norm)^(for_exp))-1))
  # print('Which components')
  # print(idx_final_sort[ii])
}



### real data
eps <- 0.05

idx <- c(1:4,6:7)

for (k in (idx)){
  x <- na.omit(AustralianWine[,k])
  plot(x, type = 'l', main = colnames(AustralianWine)[k])
  s <- ssa(x, L = 84)
  r <- reconstruct(s, groups = list(trend = 1))
  x <- x - r$trend
  s <- ssa(x, L = 84)
  print(plot(s, type = 'paired',idx=1:20))
  
  norm <- numeric(ncol(s$U) - 1)
  tau1 <- numeric(ncol(s$U) - 1)
  idx_final <- numeric(0)
  tau1[1] <- angle.fun(s$U[,1], s$U[,2])
  r <- reconstruct(s, groups = list(t = (1:2)))
  norm[1] <- sum((r$t)^2)
  
  for (j in (2:(ncol(s$U) - 1))){
    r <- reconstruct(s, groups = list(t = (j:(j+1))))
    norm[j] <- sum((r$t)^2)
    
    tau1[j] <- angle.fun(s$U[,j], s$U[,j+1])
    if (j %% 2 == 0){
      if (tau1[j] < tau1[j-1]){
        idx_final <- c(idx_final, j)
      }
      else {
        idx_final <- c(idx_final, j-1)
      }
      for (i in idx_final[-1]){
        if (idx_final[i] == idx_final[i+1]){
          if (tau1[i] > tau1[i+1]){
            idx_final <- idx_final[-which(idx_final == i)]
          }
          else {
            idx_final <- idx_final[-which(idx_final == i + 1)]
          }
        }
      }
    }
  }
  
  norm2_2 <- norm/max(norm)
  
  # print(norm2_2[1:10])
  
  for ( i in idx_final){
    if (norm2_2[i] < eps){
      idx_final <- idx_final[-which(idx_final == i)]
    }
  }
  
  norm2 <- norm[idx_final]
  tau1_2 <- tau1[idx_final]
  
  idx_final_sort <- idx_final[order(tau1_2)]
  
  
  # tau2_sort <- sort(tau2)
  dat <- data.frame(tau1 = sort(tau1_2), norm = norm2[order(tau1_2)]/max(norm2),
                    norm_1 = norm2[order(tau1_2)]/max(norm2),
                    norm_2 = (norm2[order(tau1_2)]/max(norm2))^2)
  
  plot(dat$tau1)
  tau1_diff <- numeric(length(dat$tau1)-1)
  for (i in (1:length(tau1_diff))){
    tau1_diff[i] <- dat$tau1[i+1] - dat$tau1[i]
  }
  # print(tau1_diff)
  
  x <- seq(0,1,by=0.001)
  dat2 <- data.frame(x = x, x2 = exp(x^(for_exp))-1)
  print(ggplot()+
          geom_point(data = dat, aes(norm,tau1))+
          # geom_line(data = dat2, aes(x,x2), col='red')+
          labs(title = colnames(AustralianWine)[k]))
  print(colnames(AustralianWine)[k])
  print(idx_final_sort)
  # print(norm2_2[idx_final_sort])
  # print(tau1_2)
  # ii <- which(dat$tau1 < exp((dat$norm)^(for_exp))-1)
  # print('How much to identify')
  # print(sum(dat$tau1 < exp((dat$norm)^(for_exp))-1))
  # print('Which components')
  # print(idx_final_sort[ii])
}


