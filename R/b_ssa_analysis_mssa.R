###Fragment initial
fragmentStart("fragments_ch1/mssa_init.tex")
source('main.grouping.auto.R')
library(lattice)
library(Rssa)
N.A <- 150
N.B <- 120
omega1 <- 0.1
fragmentStop()
###end 

###Fragment mssa_series
fragmentStart("fragments_ch1/mssa_series.tex")
tt.A <- 1:N.A
tt.B <- 1:N.B
X1 <- list(A = 2 * sin(2*pi * omega1 * tt.A), B = cos(2*pi * omega1 * tt.B))
X2 <- list(A = sin(2*pi * 0.5 * tt.A), B = cos(2*pi * 0.5 * tt.B))
X3 <- list(A = exp(0.01 * tt.A))
X4 <- list(A=rnorm(tt.A, 1), B=rnorm(tt.B, 1))
X <- list(A = X1$A + X2$A + X3$A+ X4$A, B = X1$B + X2$B + X4$B)
s <- ssa(X, kind = "mssa")
fragmentSkip(pdf("img_ch1/mssa_vectors.pdf", paper = "special", width = 5, height = 2))
plot(s,type='vectors')
fragmentSkip(dev.off())
fragmentStop()
###end 


###Fragment mssa_trend
fragmentStart("fragments_ch1/mssa_trend.tex")
g_trend <- general.grouping.auto(s, grouping.method = 'low.freq.mssa', 
                                 base='factor',
                                 freq.bins = list(0.01), threshold = 0.9)
print(g_trend$F1)
fragmentStop()
###end 


###Fragment mssa_em_freq
fragmentStart("fragments_ch1/mssa_em_freq.tex")
g_em_freq <- general.grouping.auto(s, grouping.method = "freq.mssa", base='factor',
                                   s_0 = 1, rho_0 = 0.9)
print(g_em_freq)
fragmentStop()
###end 

###Fragment mssa_em_tau
fragmentStart("fragments_ch1/mssa_em_tau.tex")
g_em_tau <- general.grouping.auto(s, grouping.method = "tau.mssa", treshold = 0.02)
print(g_em_tau$idx)
fragmentStop()
###end 

###Fragment mssa_recon
fragmentStart("fragments_ch1/mssa_rec.tex")
r <- reconstruct(s, groups = list(T = g_trend$F1, P = g_em_tau$idx))
fragmentSkip(pdf("img_ch1/mssa_rec.pdf", paper = "special", width = 4, height = 3))
d_A <- data.frame(T_A = r$T$A, P_A = r$P$A,  X_A = X$A, N = 1:N.A)
xyplot(T_A + P_A + X_A  ~ N, data = d_A, type ='l', ylab = '',
       auto.key = list(points = FALSE, lines = TRUE))
fragmentSkip(dev.off())
fragmentStop()
###end 
