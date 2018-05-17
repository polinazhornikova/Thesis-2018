### We use the pattern, which is extracted from the embryo bk1 (age 8) and simulate
### 200 1D profiles with random noise fitted to either additive or multiplicative models

#########
#cur_dir <- getwd() #run one time
#########


## !!!! Model data is created in 2noise1D_modeldata_generate.R

##########################
# setting the parameters #
##########################

setwd(cur_dir)

### from the file
source('params.R')

### or set the values yourself
params_1D <- list(
  # what percentage is cut in the middle
  cut_per = 0.2,
  # for SSA
  L = 30,
  # for automatic identification
  threshold = 0.4,
  freq = 0.04,
  # for LOESS
  degree = 2,
  alpha = 0.1,
  # for moving median
  k_med=51
)

setwd(cur_dir)

# setting the values for the parameters that we will evaluate
sigma_y_mult_sq <- 0.0009
sigma_y_add_sq <- 1
sigma_x_sq <- 0.1
####################################################################
# reading data and get estimates of the parameters alpha and sigma #
####################################################################

# vectors for saving results
coef_mult_x <- numeric(M)
coef_mult_y <- numeric(M)
coef_add_x <- numeric(M)
coef_add_y <- numeric(M)

results_confbound_med_mult <- numeric(M)
results_confbound_med_add <- numeric(M)

true_more_than_false <- numeric(M)

dir <- "simul_1D_mult"
setwd(dir)

for (i in 1:M) {
  
  # read data
  file <- paste0(i,".txt")
  dat <- read.table(file = file, header = TRUE)  
  
  #######################################################
  # obtaining estimates of the trend and the derivative #
  #######################################################
  
  # applying SSA
  s <- ssa(dat$S, L = params_1D$L)
  
  print(plot(s, type='vectors'))
  dat <- recon(s, dat, 1:3)
  
  # applying LOESS for trend
  temp <-
    loess(T ~ X,
          data = dat,
          degree = params_1D$degree,
          span = params_1D$alpha)
  dat$T <- predict(temp, dat$X)
  dat$R <- dat$S - dat$T
  
  # derivative
  df <- deriv1(dat$T, dat$X)
  # applying LOESS for derivative
  model <-
    loess(df ~ dat$X, degree = params_1D$degree, span = params_1D$alpha)
  dat$D <- predict(model, dat$X)
  
  dat$T2 <- dat$T ^ 2
  dat$D2 <- dat$D ^ 2
  dat$R2 <- dat$R ^ 2
  
  #################################################################
  # obtaining estimates of the parameters sigma_x^2 and sigma_y^2 #
  #################################################################
  
  # if the multiplicative model is true
  
  # applying IRLS
  res <- it.proc.mult(dat, regular = TRUE)
  dat_main <- res[[2]]
  res <- res[[1]]
  coef.mult <- res$coefficients[, 1]
  # noise variance estimate
  dat$var_mult <-
    dat$T2 * coef.mult[1]  + dat$D2 * coef.mult[2] * (1 + coef.mult[1])
  # residuals after subtracting the variance from the noise
  dat$R_off_var_mult <- dat$R2 - dat$var_mult
  # this is the estimation of the parameters sigma_x^2 and sigma_y^2
  coef_mult_x[i] <- coef.mult[2]
  coef_mult_y[i] <- coef.mult[1]
  
  
  # if the additive model is true
  
  # applying IRLS
  res <- it.proc.add(dat)
  dat_main <- res[[2]]
  res <- res[[1]]
  coef.add <- res$coefficients[, 1]
  # noise variance estimate
  dat$var_add <- coef.add[1]  + dat$D2 * coef.add[2]
  # residuals after subtracting the variance from the noise
  dat$R_off_var_add <- dat$R2 - dat$var_add
  coef_add_x[i] <- coef.add[2]
  coef_add_y[i] <- coef.add[1]
  
  
  ############################################################
  # checking which model is true, additive or multiplicative #
  ############################################################
  
  # applying the moving median with window of length 3, to remove outliers
  R_off_var_mult_med3 <- rollmedian(dat$R_off_var_mult, 3)
  R_off_var_add_med3 <- rollmedian(dat$R_off_var_add, 3)
  
  nn3 <- length(R_off_var_mult_med3) - params_1D$k_med + 1
  
  # obtaining the limits of confidence intervals
  
  conf_bound_mult_med3 <- numeric(nn3)
  conf_bound_add_med3 <- numeric(nn3)
  
  for (j in 1:nn3) {
    conf_bound_mult_med3[j] <-
      get_conf_bound_mean(R_off_var_mult_med3[j:(j + params_1D$k_med - 1)])
    conf_bound_add_med3[j] <-
      get_conf_bound_mean(R_off_var_add_med3[j:(j + params_1D$k_med - 1)])
  }
  
  results_confbound_med_add[i] <- median(conf_bound_add_med3)
  results_confbound_med_mult[i] <- median(conf_bound_mult_med3)
  true_more_than_false[i] <- median(conf_bound_mult_med3) > median(conf_bound_add_med3)
}


########################
# printing the results #
########################

print(paste0('measure_add = ' , median(results_confbound_med_add)))
print(paste0('measure_mult = ' , median(results_confbound_med_mult)))
print(paste0('measure_add = ' , mean(results_confbound_med_add)))
print(paste0('measure_mult = ' , mean(results_confbound_med_mult)))
print(paste0('num of genes = ' , length(true_more_than_false)))
print(paste0('num mult more than add = ' , sum(true_more_than_false)))
print("*****")
print(sigma_x_sq)
print(median(coef_mult_x))
print(mean(coef_mult_x))
bias(sigma_x_sq, mean(coef_mult_x))
print(sigma_y_mult_sq)
print(median(coef_mult_y))
print(mean(coef_mult_y))
bias(sigma_y_mult_sq, mean(coef_mult_y))

print(c(t.test(coef_mult_y, mu = sigma_y_mult_sq)$p.value, t.test(coef_mult_x, mu = sigma_x_sq)$p.value))



print(c(t.test(coef_add_y, mu = sigma_y_mult_sq)$p.value, t.test(coef_add_x, mu = sigma_x_sq)$p.value))


###############################
# writing the results to file #
###############################

results1D <- data.frame(true_mult = c(median(results_confbound_med_add),
                                      mean(results_confbound_med_add),
                                      median(results_confbound_med_mult),
                                      mean(results_confbound_med_mult),
                                      length(true_more_than_false),
                                      sum(true_more_than_false),
                                      sigma_x_sq,
                                      median(coef_mult_x),
                                      mean(coef_mult_x),
                                      sigma_y_mult_sq,
                                      median(coef_mult_y),
                                      mean(coef_mult_y),
                                      mean(sqrt(coef_mult_y)),
                                      sd(sqrt(coef_mult_y))))
rownames(results1D) <- c('measure_add_med', 'measure_add_mean',
                         'measure_mult_med', 'measure_mult_mean',
                         "number_of_genes",
                         'num_true_more_than_false',
                         'true_x',
                         'coef_x_med', 'coef_x_mean',
                         'true_y',
                         'coef_y_med', 'coef_y_mean',
                         'coef_y_sqrt_mean', 'coef_y_sqrt_sd')


setwd(cur_dir)
write.table(results1D, paste0("modeldata_results/1D_",M,"repeats_mult_sig_eq_", sigma_y_mult_sq, ".txt"))

results.data <- data.frame(coef_y = sqrt(coef_mult_y), coef_x=sqrt(coef_mult_x) )
write.table(results.data, paste0("modeldata_results/data/1D_",M,"repeats_mult_sig_eq_", sigma_y_mult_sq, ".txt"))

