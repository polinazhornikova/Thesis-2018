#############################
# 1D estimate for real data #
#############################

#########
#cur_dir <- getwd() #run one time
#########


##########################
# setting the parameters #
##########################

### from the file
source('params.R')

### or set the values yourself
params_1D <- list(
  # what percentage is cut in the middle
  cut_per = 0.1,
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


#############
# main part #
#############

for (dir in c('age_4', 'age_8')){
  for (gene in c('Kr', 'gt')){
    
    # vectors for saving results
    results_confbound_med_mult <- numeric(0)
    results_confbound_med_add <- numeric(0)
    coef_mult_x <- numeric(0)
    coef_mult_y <- numeric(0)
    mult_more_than_add <- numeric(0)
    sigma.mult1 <- numeric(0)
    sigma.mult2 <- numeric(0)
    sigma.mult3 <- numeric(0)
    sigma.mult4 <- numeric(0)
    
    # preparing the directory 
    setwd(cur_dir)
    files <- list.files(dir, pattern = "*.txt")
    setwd(dir)
    
    # just a loop counter
    i <- 1
    
    all_embyo <- character(0)
    
    for (file in files) {
      
      #########################
      # preparing the 1D data #
      #########################
      
      # reading 2D data
      data_main <- read.emb.data(file)
      gr <- grep(gene, names(data_main))
      if (isempty(gr)) {print("oops"); next}
      dat <- cbind(data_main[, c(2, 3)], data_main[, gene])
      colnames(dat) <- c('X', 'Y', 'S')
      
      # making one-dimensional (1D) data from two-dimensional (2D) data
      # cut out a strip along the AP (X) axis and then omit the DV (Y) coordinate
      temp <- max(dat$Y) - min(dat$Y)
      cut <- round(temp * ((1 - params_1D$cut_per) / 2))
      dat <- dat[which(dat$Y > cut & dat$Y < round(temp - cut)),]
      dat <- dat[order(dat$X),]
      # making the mean zero
      dat$S <- dat$S - min(dat$S)
      # cut the ends
      dat <- dat[which(dat$X > 20 & dat$X < 80),]
      
    
      #######################################################
      # obtaining estimates of the trend and the derivative #
      #######################################################
      
      # applying SSA
      s <- ssa(dat$S, L = params_1D$L)
      # using automatic grouping
      g <- grouping.auto(
        s,
        grouping.method = 'pgram',
        freq.bins = list(params_1D$freq)
        ,
        threshold = params_1D$threshold
      )
      dat <- recon.groups(s, dat, g)
      
      # applying LOESS for trend
      temp <-
        loess(T ~ X,
              data = dat,
              degree = params_1D$degree,
              span = params_1D$alpha)
      
      # trend and residual estimates
      dat$T <- predict(temp, dat$X)
      dat$R <- dat$S - dat$T
      print(xyplot(T + R + S ~ X, data = dat, type = "l", main = file))
      
      # derivative
      df <- deriv1(dat$T, dat$X)
      # applying LOESS for derivative
      model <-
        loess(df ~ dat$X, degree = params_1D$degree, span = params_1D$alpha)
      temp <- predict(model, dat$X)
      # estimation of derivative
      dat$D <- temp
      dat <- dat[2:(nrow(dat) - 1), ]
      
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
      coef_mult_x <- c(coef_mult_x, coef.mult[2])
      coef_mult_y <- c(coef_mult_y, coef.mult[1])
      
      
      # if the additive model is true
      
      # applying IRLS
      res <- it.proc.add(dat)
      dat_main <- res[[2]]
      res <- res[[1]]
      coef <- res$coefficients[, 1]
      # noise variance estimate
      dat$var_add <- coef[1]  + dat$D2 * coef[2]
      # residuals after subtracting the variance from the noise
      dat$R_off_var_add <- dat$R2 - dat$var_add
      
      
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
      
      results_confbound_med_add <-c(results_confbound_med_add,median(conf_bound_add_med3))
      results_confbound_med_mult <- c(results_confbound_med_mult,median(conf_bound_mult_med3))
      mult_more_than_add <- c(mult_more_than_add,
                              median(conf_bound_mult_med3) > median(conf_bound_add_med3))
      
      all_embyo <- c(all_embyo, file)
      
      
    }
    
    
    ########################
    # printing the results #
    ########################
    
    print(paste0('measure_add = ' , median(results_confbound_med_add)))
    print(paste0('measure_mult = ' , median(results_confbound_med_mult)))
    print(paste0('measure_add = ' , mean(results_confbound_med_add)))
    print(paste0('measure_mult = ' , mean(results_confbound_med_mult)))
    print(paste0('num of genes = ' , length(mult_more_than_add)))
    print(paste0('num mult more than add = ' , sum(mult_more_than_add)))
    print(sqrt(median(coef_mult_x)))
    print(sqrt(mean(coef_mult_x)))
    print(sqrt(median(coef_mult_y)))
    print(sqrt(mean(coef_mult_y)))
    
    print(mean(sqrt(coef_mult_x)))
    print(mean(sqrt(coef_mult_x)))
    print(median(sqrt(coef_mult_y)))
    
    print("*****")
    print(sqrt(mean(coef_mult_y)))
    print(mean(sqrt(coef_mult_y)))
    
    
    ###############################
    # writing the results to file #
    ###############################
    
    results1D_real <- data.frame(mult = c(median(results_confbound_med_add),
                                          mean(results_confbound_med_add),
                                          median(results_confbound_med_mult),
                                          mean(results_confbound_med_mult),
                                          length(mult_more_than_add),
                                          sum(mult_more_than_add),
                                          median(sqrt(coef_mult_x)),
                                          mean(sqrt(coef_mult_x)),
                                          median(sqrt(coef_mult_y)),
                                          mean(sqrt(coef_mult_y)),
                                          mean(coef_mult_x),
                                          mean(coef_mult_y),
                                          sd(sqrt(coef_mult_x)),
                                          sd(sqrt(coef_mult_y))))
    rownames(results1D_real) <- c('measure_add_med', 'measure_add_mean',
                                  'measure_mult_med', 'measure_mult_mean',
                                  'number_of_genes',
                                  'num_mult_more_than_add',
                                  'coef_x_med', 'coef_x_mean',
                                  'coef_y_med', 'coef_y_mean',
                                  'coef_x_sq_mean', 'coef_y_sq_mean',
                                  'coef_x_sd', 'coef_y_sd')
    
    setwd(cur_dir)
    write.table(results1D_real, paste0("realdata_results_auto/1D_", dir, "_",gene, ".txt"))
    
    results.data <- data.frame(coef_y = sqrt(coef_mult_y), coef_x=sqrt(coef_mult_x) )
    row.names(results.data) <- all_embyo
    write.table(results.data, paste0("realdata_results_auto/data/1D_", dir, "_",gene, ".txt"))
  }
}