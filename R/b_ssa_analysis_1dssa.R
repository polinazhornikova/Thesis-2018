###Fragment initial
fragmentStart("fragments_ch1/1dssa_init.tex")
source('main.grouping.auto.R')
library(lattice)
library(Rssa)
N <- 199
omega1 <- 0.1
omega2 <- 0.3
fragmentStop()
###end 

###Fragment 1dssa_series
fragmentStart("fragments_ch1/1dssa_series.tex")
x <- exp(0.01 * (1:N)) + 2*cos(2 * pi * omega1 * (1:N)) +  
  exp(0.009 * (1:N)) * cos(2 * pi * omega2 * (1:N)) + 
  2*cos(2 * pi * 0.5 * (1:N))+ rnorm(N,1.5)
s <- ssa(x, L = 100)
fragmentSkip(pdf("img_ch1/1dssa_vectors.pdf", paper = "special", width = 5, height = 2))
plot(s, type="vectors", layout = c(5, 2))
fragmentSkip(dev.off())
fragmentSkip(pdf("img_ch1/1dssa_vectors_pair.pdf", paper = "special", width = 5, height = 2))
plot(s, type="paired", layout = c(5, 2))
fragmentSkip(dev.off())
fragmentStop()
###end 

###Fragment 1dssa_trend
fragmentStart("fragments_ch1/1dssa_trend.tex")
g_trend <- grouping.auto(s, grouping.method = 'pgram',
                         freq.bins = list(0.01), threshold = 0.9)
print(g_trend$F1)
fragmentStop()
###end 

###Fragment 1dssa_em_freq
fragmentStart("fragments_ch1/1dssa_em_freq.tex")
g_em_freq <- general.grouping.auto(s, grouping.method = "freq.1dssa", s_0 = 1, 
                                   rho_0 = 0.9)
print(g_em_freq)
fragmentStop()
###end 

###Fragment 1dssa_em_tau
fragmentStart("fragments_ch1/1dssa_em_tau.tex")
g_em_tau <- general.grouping.auto(s, grouping.method = "tau.1dssa", treshold = 0.01)
print(g_em_tau$idx)
fragmentStop()
###end 

###Fragment 1dssa_recon
fragmentStart("fragments_ch1/1dssa_rec.tex")
r <- reconstruct(s, groups = list(T = g_trend, P = c(g_em_freq$I_1, 
                                                     g_em_freq$I_2)))
fragmentSkip(pdf("img_ch1/1dssa_rec.pdf", paper = "special", width = 4, height = 3))
d <- data.frame(X=x, N=1:N, T=r$T, P=r$P)
xyplot(T + P +X ~ N, data = d, type ='l', ylab = '',
       auto.key = list(points = FALSE, lines = TRUE))
fragmentSkip(dev.off())
fragmentStop()
###end 
