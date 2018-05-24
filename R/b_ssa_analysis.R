###Fragment initial
fragmentStart("fragments_ch1/initial.tex")
source('main.grouping.auto.R')
N <- 199
omega1 <- 0.1
omega2 <- 0.25
fragmentStop()
###end 

###Fragment 1dssa_series
fragmentStart("fragments_ch1/1dssa_series.tex")
d <- data.frame(x = exp(0.01 * (1:N)) + 2*cos(2 * pi * omega1 * (1:N)) +  
  exp(0.009 * (1:N)) * cos(2 * pi * omega2 * (1:N)) + 
  2*cos(2 * pi * 0.5 * (1:N))+ rnorm(N,1.2), N = 1:N)
fragmentSkip(pdf("img_ch1/1dssa_series.pdf", paper = "special", width = 4, height = 3))
xyplot(x ~ N, data=d, type="l")
fragmentSkip(dev.off())
s <- ssa(x)
fragmentStop()
###end 
