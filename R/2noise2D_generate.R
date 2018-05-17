######################################################################
# check the quality of sigma estimate and generate model data for 2D #
######################################################################

#########
#cur_dir <- getwd() #run one time
#########


######################
# setting parameters #
######################

### from the file
setwd(cur_dir)
source('params.R')

### or set the values yourself
xlim <- c(20, 80)
ylim <- c(20, 80)
# for SSA
step <- 0.5
size = 5
L <- c(round(size/step), round(size/step))
# number of components
r_4 <- 6  # for age 4
r_8 <- 6  # for age 8

# number of repeats
M <- 200

setwd(cur_dir)


###########################################################################
# obtaining the real trend of the embryo bk1 of the Kruppel gene of age 8 #
###########################################################################

dir <- "age_8"
gene <- "Kr"
file <- "bk1.txt"

# reading data
setwd(dir)
df_bk1 <- read.emb.data(file)

# applying SSA
bss_bk1 <- BioSSA(get(gene) ~ AP + DV,  data = df_bk1, 
                  L = L,
                  step = step,
                  xlim = xlim, ylim = ylim)
rec_bk1 <- reconstruct(bss_bk1, groups = list(1:r_8))

# checking of decomposition quality
plot(residuals(bss_bk1, 1:r_8), type = "nuclei-section",
     coord = "x", at = 50, tol = 5,
     ref = TRUE, col = "blue")
plot(residuals(bss_bk1, 1:r_8), type = "nuclei-section",
     coord = "y", at = 50, tol = 5,
     ref = TRUE, col = "blue")
plot(noise.model(bss_bk1, groups = 1:r_8, model = "mult"))
### it's ok

r_bk1 <- residuals(bss_bk1, 1:r_8)
#plot(rec)
#print(plot(rec$F1, main = file))
print(plot(residuals(rec_bk1), main = file))

#residual <- abs(residuals(rec_bk1)$field$z)
#trend <- abs(rec_bk1$F1$field$z)   
#sigma3 <- mean((residual/trend)^2, na.rm = TRUE)
#standard <- sqrt(sigma3)
#standard

residual <- abs(residuals(rec_bk1)$values)
trend <- abs(rec_bk1$F1$values)   
sigma3 <- mean((residual/trend)^2, na.rm = TRUE)
# standard <- sqrt(sigma3)
standard <- sqrt(sigma3)
standard #correct version

############################################
# making model data and writing it to file #
############################################

# setting the values of sigma and alpha for model data (multiplicative model)
alpha_true <- 1
# sigma_true <- standard # sigma ^ 2 = 0.0001 (similarly as in 1D)
sigma_true <- 0.03 # sigma ^ 2 = 0.0009 (similarly as in 1D)


# trend of the embryo bk1
values <- rec_bk1$F1$values

length(values)

# making model data and writing it to file
# M = sample size
set.seed(23)
setwd(cur_dir)
for (i in 1:M){
  newsignal <- values +  values ^ alpha_true * sigma_true * rnorm(length(values))
  newdata <- data.frame(N=(0:(length(values)-1)), AP=rec_bk1$F1$x2d, DV=rec_bk1$F1$y2d, Kr=newsignal)
  #newdata <- na.omit(newdata)
  setwd(cur_dir)
  write.table(newdata,  paste0("simul_2D_mult/",i,".txt"))
}



##########################################################################
# reading data and obtaining estimates of the parameters alpha and sigma #
##########################################################################

# vectors for saving results
alpha1 <- numeric(M)
alpha2 <- numeric(M)
sigma.mult1 <- numeric(M)
sigma.mult2 <- numeric(M)
sigma.mult3 <- numeric(M)
sigma.mult4 <- numeric(M)
sig <- numeric(M)
prop.corr <- numeric(M)
prop.mean <- numeric(M)
prop.sd <- numeric(M)
prop.lm.coef <- numeric(M)
prop.lm.pvalue <- numeric(M)
prrp.cos <- numeric(M)

setwd(cur_dir)
dir <- "simul_2D_mult"
setwd(dir)

for (i in 1:M) {
  
  # reading data
  file <- paste0(i,".txt")
  newdata <- read.table(file = file, header = TRUE) 
  
  # applying SSA
  bss <- BioSSA(get(gene) ~ AP + DV,  data = newdata, 
                L = L,
                step = step,
                xlim = c(0,100), ylim = c(0,100))
  rec <- reconstruct(bss, groups = list(1:r_8))
  print(plot(residuals(rec), main = file))
  
  print(plot(rec))
  #print(plot(rec_bk1))
  
  # checking of decomposition quality
  print(plot(residuals(bss, 1:r_8), type = "nuclei-section",
             coord = "x", at = 50, tol = 5,
             ref = TRUE, col = "blue"))
  print(plot(residuals(bss, 1:r_8), type = "nuclei-section",
             coord = "y", at = 50, tol = 5,
             ref = TRUE, col = "blue"))
  print(plot(noise.model(bss, groups = 1:r_8, model = "mult")))
  #print(plot(noise.model(bss, groups = 1:r_8)))
  ### this is not good :(
  
  
  # obtaining different estimates of sigma and alpha
  
  #residual <- abs(residuals(rec)$field$z)
  #trend <- abs(rec$F1$field$z)   
  residual <- abs(residuals(rec)$values)
  trend <- abs(rec$F1$values)   
  sigma3 <- mean((residual/trend)^2, na.rm = TRUE)
  sigma.mult1[i] <- exp(mean(log(residual)-log(trend), na.rm = TRUE))
  sigma.mult2[i] <- mean(exp(log(residual)-log(trend)), na.rm = TRUE)
  sigma.mult3[i] <- mean((residual/trend)^2, na.rm = TRUE)
  sigma.mult4[i] <- mean((residual/trend), na.rm = TRUE)
  
  R <- lm(as.vector(log(residual)) ~ as.vector(log(trend)), na.action=na.omit)
  #plot(as.vector(residuals(rec)$values) ~ as.vector(rec$F1$values))
  #plot(as.vector(residuals(rec_bk1)$values) ~ as.vector(rec_bk1$F1$values))
  alpha2[i] <- as.numeric(coef(R)[2])
  
  noise <- noise.model(bss, 1:r_8, averaging.type = "none", offset = 0)
  alpha1[i] <- noise$alpha
  sig[i] <- noise$sigma
  
  # correlation
  x <- rec$F1$field$xperc
  y <- rec$F1$field$yperc
  z <- rec$F1$field$z
  middle <- 50
  lower <- middle+5
  upper <- middle-5
  pmiddle.ssa <- z[round(middle/step),]
  pupper.ssa <- z[round(upper/step),]
  plower.ssa <- z[round(lower/step),]
  pdiff.y.ssa <- ((pupper.ssa - plower.ssa)/(upper-middle))[2:(100/step-1)]
  pdiff.x.ssa <- diff(pmiddle.ssa, lag = 2)/2
  #I am not sure that correlation is good. We should check proportionality
  prop.corr[i] <- cor(abs(pdiff.y.ssa),abs(pdiff.x.ssa), use =  "complete.obs")
  
  prop.mean[i] <- mean(abs(pdiff.y.ssa)/abs(pdiff.x.ssa), na.rm = TRUE)
  prop.sd[i] <- sd(abs(pdiff.y.ssa)/abs(pdiff.x.ssa), na.rm = TRUE)
  
  prop.lm.coef[i] <- summary(lm(abs(pdiff.y.ssa) ~ 0+ abs(pdiff.x.ssa)))$coefficients[,1]
  prop.lm.pvalue[i] <- summary(lm(abs(pdiff.y.ssa) ~ 0+ abs(pdiff.x.ssa)))$coefficients[,4]
  
  prop.cos[i] <- cosxy(abs(pdiff.y.ssa), abs(pdiff.x.ssa))
  #print(corr[i])
}




########################
# printing the results #
########################

print(mean(alpha1))
print(mean(alpha2))
print(mean(corr))
print("******")
print(mean((sigma.mult1)))
print(mean((sigma.mult2)))
print("******")
print(mean(sqrt(sigma.mult3)))
print(sqrt(mean(sigma.mult3)))

print(standart)

print(mean((sigma.mult4)))
print("******")
print(mean(sig))


###############################
# writing the results to file #
###############################

results2D <- data.frame(true_mult = c(mean(alpha1),
                                      sigma_true,
                                      sqrt(mean(sigma.mult3)),
                                      mean(sqrt(sigma.mult3)),
                                      sd(sqrt(sigma.mult3))))
rownames(results2D) <- c('alpha_mean',
                         'true_sigma',
                         'sqrt_mean_sigma3',
                         'mean_sqrt_sigma3',
                         'sd_sqrt_sigma3')

setwd(cur_dir)
write.table(results2D, paste0("modeldata_results/2D_",M,"repeats_sigma_eq_",sigma_true,".txt"))

results.data <- data.frame(coef_y_2D = sqrt(sigma.mult3))
write.table(results.data, paste0("modeldata_results/data/2D_",M,"repeats_sigma_eq_",sigma_true,".txt"))

results.data.prop <- data.frame(prop.corr,prop.mean, prop.sd, prop.lm.coef, prop.lm.pvalue, prop.cos)
write.table(results.data.prop, paste0("modeldata_results/data/prop_2D_",M,"repeats_sigma_eq_",sigma_true,".txt"))
