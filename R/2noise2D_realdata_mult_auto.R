#############################
# 2D estimate for real data #
#############################


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
# params for auto SSA
freq1 <- 0.08 * step
freq2 <- 0.1 * step
th <- 0.2

setwd(cur_dir)
i <- 0


for (dir in c('age_4', 'age_8')){
  for (gene in c('Kr','gt')){
    i <- i + 1
    
    setwd(cur_dir)
    
    files <- list.files(dir, pattern = "*.txt")
    setwd(dir)
    
    alpha1 <- numeric(0)
    alpha2 <- numeric(0)
    sig <- numeric(0)
    sigma.mult1 <- numeric(0)
    sigma.mult2 <- numeric(0)
    sigma.mult3 <- numeric(0)
    sigma.mult4 <- numeric(0)
    prop.corr <- numeric(0)
    prop.mean <- numeric(0)
    prop.sd <- numeric(0)
    prop.lm.coef <- numeric(0)
    prop.lm.pvalue <- numeric(0)
    prop.cos <- numeric(0)
    
    all_embyo <- character(0)
    
    for (file in files){
      
      df <- read.emb.data(file)
      #df[gene] <- df[gene] - min(df[gene]) 
      gr <- grep(gene, names(df))
      if (isempty(gr)) next;
      
      bss <- BioSSA(get(gene) ~ AP + DV,  data = df, 
                    L = L,
                    step = step,
                    xlim = xlim, ylim = ylim)
      
      g <- grouping.auto.pgram.2d.ssa_my_bio(bss$ssa,freq.bins1 = freq1,freq.bins2 =freq2,
                                             threshold = th)$g
      
      rec <- reconstruct(bss, groups = g)
      
      print(plot(residuals(bss, g), type = "nuclei-section",
           coord = "x", at = 50, tol = 5,
           ref = TRUE, col = "blue", main=file))
      print(plot(residuals(bss, g), type = "nuclei-section",
           coord = "y", at = 50, tol = 5,
           ref = TRUE, col = "blue", main=file))
      print(plot(noise.model(bss, groups = g, model = "mult"), main=file))
      
      residual <- abs(residuals(rec)$field$z)
      trend <- abs(rec$F1$field$z)   
      sigma1 <- exp(mean(log(residual)-log(trend), na.rm = TRUE))
      sigma2 <- mean(exp(log(residual)-log(trend)), na.rm = TRUE)
      
      sigma3 <- mean((residual/trend)^2, na.rm = TRUE)
      sigma4 <- mean((residual/trend), na.rm = TRUE)
      
      sigma.mult1 <- c(sigma.mult1, sigma1)
      sigma.mult2 <- c(sigma.mult2, sigma2)
      sigma.mult3 <- c(sigma.mult3, sigma3)
      sigma.mult4 <- c(sigma.mult4, sigma4)
      
      R <- lm(as.vector(log(residual)) ~ as.vector(log(trend)), na.action=na.omit)
      alpha <- as.numeric(coef(R)[2])
      
      noise <- noise.model(bss, g, averaging.type = "none", offset = 0)
      alpha1 <- c(alpha1, noise$alpha)
      alpha2 <- c(alpha2, alpha)
      sig <- c(sig, noise$sigma)
      
      x <- rec$F1$field$xperc
      y <- rec$F1$field$yperc
      z <- rec$F1$field$z
      
      #dy <- rec$F1$field$yperc[2]-rec$F1$field$yperc[1]
      #dx <- rec$F1$field$xperc[2]-rec$F1$field$xperc[1]
      #rec$F1$field$xperc[round(50/dx)]
      #rec$F1$field$yperc[round(50/dy)]
      middle <- 50
      lower <- middle+5
      upper <- middle-5
      pmiddle.ssa <- z[round(middle/step),]
      pupper.ssa <- z[round(upper/step),]
      plower.ssa <- z[round(lower/step),]
      
      pdiff.y.ssa <- ((pupper.ssa - plower.ssa)/(upper-middle))[2:(100/step-1)]
      
      pdiff.x.ssa <- diff(pmiddle.ssa, lag = 2)/2
      
      #plot(pdiff.y.ssa ~ pdiff.x.ssa)
      #plot(abs(pdiff.y.ssa)~abs(pdiff.x.ssa))
      
      #I am not sure that correlation is good. We should check proportionality
      c <- cor(abs(pdiff.y.ssa),abs(pdiff.x.ssa), use =  "complete.obs")
      prop.corr <- c(prop.corr, c)
      prop.mean <- c(prop.mean, mean(abs(pdiff.y.ssa)/abs(pdiff.x.ssa), na.rm = TRUE))
      prop.sd <- c(prop.sd ,sd(abs(pdiff.y.ssa)/abs(pdiff.x.ssa), na.rm = TRUE))
      
      prop.lm.coef <- c(prop.lm.coef, summary(lm(abs(pdiff.y.ssa) ~ 0+ abs(pdiff.x.ssa)))$coefficients[,1])
      prop.lm.pvalue <- c(prop.lm.pvalue, summary(lm(abs(pdiff.y.ssa) ~ 0+ abs(pdiff.x.ssa)))$coefficients[,4])
      
      prop.cos <- c(prop.cos, cosxy(abs(pdiff.y.ssa), abs(pdiff.x.ssa)))
        
      all_embyo <- c(all_embyo, file)
      #print(c)
    }
    
    print(mean(alpha1))
    print(mean(alpha2))
    print(mean(prop.corr))
    
    
    print(mean((sigma.mult1)))
    print(mean((sigma.mult2)))
    print("******")
    print(mean(sqrt(sigma.mult3)))
    print(sqrt(mean(sigma.mult3)))
    print(mean((sigma.mult4)))
    
    
    results2D_real <- data.frame(true_mult = c(mean(alpha1),
                                               sqrt(mean(sigma.mult3)),
                                               mean(sqrt(sigma.mult3)),
                                               sd(sqrt(sigma.mult3))))
    rownames(results2D_real) <- c('alpha_mean',
                                  'sqrt_mean_sigma3',
                                  'mean_sqrt_sigma3',
                                  'sd_sqrt_sigma3')

    setwd(cur_dir)
    write.table(results2D_real, paste0("realdata_results_auto/2D_", dir, "_",gene, ".txt"))
    results.data <- data.frame(coef_y = sqrt(sigma.mult3))
    row.names(results.data) <- all_embyo
    write.table(results.data, paste0("realdata_results_auto/data/2D_", dir, "_",gene, ".txt"))
    
    results.data.prop <- data.frame(prop.corr,prop.mean, prop.sd, prop.lm.coef, prop.lm.pvalue, prop.cos)
    row.names(results.data.prop) <- all_embyo
    write.table(results.data.prop, paste0("realdata_results_auto/data/prop_2D_", dir, "_",gene, ".txt"))
    
  }
}


