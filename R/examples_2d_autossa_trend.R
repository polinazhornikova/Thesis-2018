library("lattice")
library('Rssa')

plot2d <- function(x) {
  regions <- list(col = colorRampPalette(grey(c(0, 1))));
  # m <- t(x[seq(nrow(x), 1, -1), ])
  m <- t(x)
  d = data.frame(x=rep(seq(0, nrow(m), length=nrow(m)), ncol(m)), 
                 y=rep(seq(0, ncol(m), length=ncol(m)), each=nrow(m)), 
                 z=c(m))
  levelplot(z~x*y, data=d, aspect = "iso",
            par.settings = list(regions = regions), colorkey = TRUE,
            xlab = "", ylab = "", scales=list(x=list(at=c(0, 50, 100, 150) 
                                                     # ,labels=c('0','1/3','2/3','1')
                                                     ),
                                              y=list(at=c(0, 50, 100)
                                                     # ,labels=c('0','1/2','1')
                                                     )))
}

centered.mod.fft.2d <- function(x) {
  N <- dim(x)
  shift.exp <- exp(2i * pi * floor(N/2) / N)
  shift1 <- shift.exp[1]^(0:(N[1] - 1))
  shift2 <- shift.exp[2]^(0:(N[2] - 1))
  Mod(t(mvfft(t(mvfft(outer(shift1, shift2) * x)))))
}

mod.fft.2d <- function(x) {
  Mod(t(mvfft(t(mvfft(x)))))
}

N <- 100
M <- 150
alpha <- 0.01
matr <- matrix(1,N,M)

omega1 <- 0.1
omega2 <- 0.25


for (i in (1:N)){
  for (j in (1:M)){
    matr[i,j] <- cos(2 * pi * omega1 * i) * cos(2 * pi * omega2 *j) + sin(2 * pi * omega1 * i) * sin(2 * pi * omega2 *j) + rnorm(1,sd=0.01) 
  }
}
matr_cos_rank2 <- matr

for (i in (1:N)){
  for (j in (1:M)){
    matr[i,j] <- cos(2 * pi * omega1 * i) * cos(2 * pi * omega2 *j) + rnorm(1,sd=0.05)
  }
}

matr_cos_rank4 <- matr

for (i in (1:N)){
  for (j in (1:M)){
    matr[i,j] <- exp(0.01 * i) + cos(2 * pi * omega3 * i) * cos(2 * pi * omega4 *j)  + rnorm(1,sd=0.5)
      # cos(2 * pi * omega1 * i) + 
      
  }
}

matr_cos_rank8 <- matr

plot2d(matr_cos_rank8)
plot2d(matr_cos_rank4)
plot2d(matr_cos_rank2)


s_rank2 <- ssa(matr_cos_rank2, kind = "2d-ssa", L = c(8, 8))
plot(s_rank2, type = "vectors", idx = 1:50, cuts = 255)
plot(wcor(s_rank2, groups= 1:50), scales = list(at = c(10, 20, 30)))


s_rank4 <- ssa(matr_cos_rank4, kind = "2d-ssa", L = c(8, 8))
plot(s_rank4, type = "vectors", idx = 1:32, cuts = 255, layout = c(8, 4))
plot(wcor(s_rank4, groups= 1:8), scales = list(at = c(10, 20, 30)))


s_rank8 <- ssa(matr_cos_rank8, kind = "2d-ssa", L = c(8, 8))
plot(s_rank8, type = "vectors", idx = 1:32, cuts = 255, layout = c(8, 4))
plot(wcor(s_rank8, groups= 1:15), scales = list(at = c(10, 20, 30)))

plot2d(centered.mod.fft.2d(matr_cos_rank8))
plot2d(centered.mod.fft.2d(matr_cos_rank4))
plot2d(centered.mod.fft.2d(matr_cos_rank2))
# plot2d(mod.fft.2d(matr_cos_rank2))

g <- grouping.auto.pgram.2d.ssa_my(s_rank8,freq.bins1 = 0.05,freq.bins2 = 0.05,
                              threshold = 0.6)


d2 <- data.frame(T=sort(g$contr, decreasing = TRUE), N = 1:length(g$contr))

xyplot(T ~ N, data = d2, xlab = 'Component')

# grouping.auto.pgram.2d.ssa_my(s_rank4,freq.bins1 = 0.1,freq.bins2 = 0.1,
#                               threshold = 0.4)
# 
# grouping.auto.pgram.2d.ssa_my(s_rank2,freq.bins1 = 0.1,freq.bins2 = 0.1,
#                               threshold = 0.4)



N <- 100
M <- 150
alpha <- 0.01
matr <- matrix(1,N,M)

omega1 <- 0.1
omega2 <- 0.25


for (i in (1:N)){
  for (j in (1:M)){
    matr[i,j] <- exp(0.01 * i)+ cos(2 * pi * omega1 * i) * cos(2 * pi * omega2 *j) + rnorm(1,sd=0.5)
  }
}

matr_cos_rank4 <- matr

plot2d(matr_cos_rank4)

s_rank4 <- ssa(matr_cos_rank4, kind = "2d-ssa", L = c(50, 50))
plot(s_rank4, type = "vectors", idx = 1:32, cuts = 255, layout = c(8, 4))
plot(wcor(s_rank4, groups= 1:8), scales = list(at = c(10, 20, 30)))

g <- AI_trend_freq_2d(s_rank4,freq.bins1 = 0.1,freq.bins2 = 0.1,
                                   threshold = 0.2, groups=(1:7))
g$g

r <- reconstruct(s_rank4, groups=list(g$g))

plot2d(r$F1)
