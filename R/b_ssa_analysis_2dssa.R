###Fragment initial
fragmentStart("fragments_ch1/2dssa_init.tex")
source('main.grouping.auto.R')
library(lattice)
library(Rssa)
N_x <- 99
N_y <- 149
omega1 <- 0.1
fragmentStop()
###end 

###Fragment 2dssa_series
fragmentStart("fragments_ch1/2dssa_series.tex")
matr <- matrix(1,N_x,N_y)
for (i in (1:N_x)){
  for (j in (1:N_y)){
    matr[i,j] <- exp(0.01 * i) + cos(2 * pi * omega1 * (i + j))  + 
      rnorm(1,sd=0.5)
  }
}
fragmentSkip(pdf("img_ch1/2dssa_series.pdf", paper = "special", width = 9, height = 6))
plot2d(matr)
fragmentSkip(dev.off())
s <- ssa(matr, kind = "2d-ssa", L = c(50, 50))
fragmentSkip(pdf("img_ch1/2dssa_vectors.pdf", paper = "special", width = 4, height = 2))
plot(s, type = "vectors", cuts = 255, layout = c(4, 2))
fragmentSkip(dev.off())
fragmentStop()
###end 


###Fragment 2dssa_trend
fragmentStart("fragments_ch1/2dssa_trend.tex")
g_trend <- draft.grouping.auto(s, grouping.method="low.freq.2dssa", 
                               freq.bins1 = 0.1, freq.bins2 = 0.1,
                               threshold = 0.7, base='series')
print(g_trend$g)
fragmentStop()
###end 


###Fragment 2dssa_recon
fragmentStart("fragments_ch1/2dssa_rec.tex")
r <- reconstruct(s, groups = list(T = g_trend$g))
fragmentSkip(pdf("img_ch1/2dssa_rec.pdf", paper = "special", width = 9, height = 6))
plot2d(r$T)
fragmentSkip(dev.off())
fragmentStop()
###end 
