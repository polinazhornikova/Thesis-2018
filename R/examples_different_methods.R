library(lattice)
source('main.grouping.auto.R')

N <- 199
omega1 <- 0.1
omega2 <- 0.25

### 1D-SSA
x <- exp(0.01 * (1:N)) + 2*cos(2 * pi * omega1 * (1:N)) +  exp(0.009 * (1:N)) * cos(2 * pi * omega2 * (1:N)) + 2*cos(2 * pi * 0.5 * (1:N))+ rnorm(N,1.2)
plot(x, type="l")
s <- ssa(x)
plot(s, type="vectors")
plot(s, type="paired")

# trend
g <- grouping.auto(
  s,
  grouping.method = 'pgram',
  freq.bins = list(0.01)
  ,
  threshold = 0.9
)
r <- reconstruct(s, groups = list(T = g))
d <- data.frame(T = r$T, X = x, N = 1:N)
xyplot(X + T ~ N, data = d, type ='l')

# e-m garm
# tau 
g <- draft.grouping.auto(s, grouping.method = "tau.1dssa")
print(g$idx)
# print(g$tau)
r <- reconstruct(s, groups = list(S = g$idx))
d <- data.frame(S = r$S, X = x, N = 1:N)
xyplot(X + S ~ N, data = d, type ='l')

# e-m garm
# paired frequency
g <- draft.grouping.auto(s, grouping.method = "pair.freq.1dssa")
print(g)
r <- reconstruct(s, groups = list(S=c(g$I_1, g$I_2)))
d <- data.frame(S = r$S, X = x, N = 1:N)
xyplot(X + S ~ N, data = d, type ='l')


### CSSA
x <- exp(0.01 * (1:N)) + cos(2 * pi * omega1 * (1:N)) +1.0i*sin(2 * pi * omega1 * (1:N)) + 1.0i*cos(2 * pi * omega2 * (1:N)) +rnorm(N, 0.2) 
s <- ssa(x, kind = 'cssa')
plot(s, type="vectors")

# trend
g <- draft.grouping.auto(s, grouping.method = 'low.freq.cssa', freq.bins = list(0.01), threshold = 0.9)
print(g$F1)
r <- reconstruct(s, groups = list(T=g))
d_re <- data.frame(T_re = Re(r$T),  X_re = Re(x),N = 1:N)
xyplot(X_re  + T_re  ~ N, data = d_re, type ='l')
# d_im <- data.frame(T_im = Im(r$F1), X_im = Im(x), N = 1:N)
# xyplot(X_im + T_im ~ N, data = d_im, type ='l')


# e-m garm
# tau 

# Drawing scatterplot for d = 1
plot.d1 <- function(s,index=1:8){
  index <- unique(index)
  u <- s$U[,index]
  data <- data.frame(z = as.vector(u), 
                     g = factor(rep(seq_len(ncol(u)),                     
                                    each = nrow(u)),
                                labels = paste("U", seq_len(ncol(u)))))
  xyplot(Im(z) ~ Re(z) | g, data = data, type = "l", 
         as.table = TRUE,  
         scales = list(relation = "free", draw = FALSE))
}


plot.d1(s)
plot(s, type="paired")
g <- draft.grouping.auto(s, grouping.method = "tau.cssa",  treshold=0.05)
print(g$d1_idx)
# print(g$tau_d1)
print(g$d2_idx)
# print(g$tau_d2)
r <- reconstruct(s, groups = list(S=c(g$d1_idx, g$d2_idx)))
d_re <- data.frame(S_re = Re(r$S),  X_re = Re(x),N = 1:N)
xyplot(X_re  + S_re  ~ N, data = d_re, type ='l')
d_im <- data.frame(S_im = Im(r$S), X_im = Im(x), N = 1:N)
xyplot(X_im + S_im ~ N, data = d_im, type ='l')


### MSSA
N.A <- 150
N.B <- 120
omega1 <- 0.1
omega2 <- 0.05
tt.A <- 1:N.A
tt.B <- 1:N.B
F1 <- list(A = 2 * sin(2*pi * omega1 * tt.A), B = cos(2*pi * omega1 * tt.B))
F2 <- list(A = 1 * sin(2*pi * omega2 * tt.A), B = cos(2*pi * omega2 * tt.B))
F3 <- list(A = exp(0.01 * tt.A))
F4 <- list(A=rnorm(tt.A), B=rnorm(tt.B))
F <- list(A = F1$A + F2$A + F3$A+ F4$A, B = F1$B + F2$B + F4$B)

s <- ssa(F, kind = "mssa")
plot(s,type='vectors')

# trend
# series
g <- draft.grouping.auto(s, grouping.method = 'low.freq.mssa', base='series', freq.bins = list(0.01), threshold = 0.9)
print(g$F1)
r <- reconstruct(s, groups = list(T = g$F1))

d_A <- data.frame(T_A = r$T$A,  X_A = F$A, N = 1:N.A)
xyplot(X_A  + T_A  ~ N, data = d_A, type ='l')
d_B <- data.frame(T_B = r$T$B,  X_B = F$B, N = 1:N.B)
xyplot(X_B  + T_B  ~ N, data = d_B, type ='l')

# trend
# factor
g <- draft.grouping.auto(s, grouping.method = 'low.freq.mssa', base='factor', freq.bins = list(0.01), threshold = 0.9)
print(g$F1)
r <- reconstruct(s, groups = list(T = g$F1))

d_A <- data.frame(T_A = r$T$A,  X_A = F$A, N = 1:N.A)
xyplot(X_A  + T_A  ~ N, data = d_A, type ='l')
d_B <- data.frame(T_B = r$T$B,  X_B = F$B, N = 1:N.B)
xyplot(X_B  + T_B  ~ N, data = d_B, type ='l')

# trend
# eigen
g <- draft.grouping.auto(s, grouping.method = 'low.freq.mssa', base='eigen', freq.bins = list(0.01), threshold = 0.9)
print(g$F1)
r <- reconstruct(s, groups = list(T = g$F1))

d_A <- data.frame(T_A = r$T$A,  X_A = F$A, N = 1:N.A)
xyplot(X_A  + T_A  ~ N, data = d_A, type ='l')
d_B <- data.frame(T_B = r$T$B,  X_B = F$B, N = 1:N.B)
xyplot(X_B  + T_B  ~ N, data = d_B, type ='l')


### 2D-SSA
plot2d <- function(x) {
  regions <- list(col = colorRampPalette(grey(c(0, 1))));
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

N <- 99
M <- 149
alpha <- 0.01
omega1 <- 0.001
omega2 <- 0.0025
matr <- matrix(1,N,M)


for (i in (1:N)){
  for (j in (1:M)){
    matr[i,j] <- cos(2 * pi * omega1 * i) * cos(2 * pi * omega2 *j) + sin(2 * pi * omega1 * i) * sin(2 * pi * omega2 *j) + rnorm(1,sd=0.01) 
  }
}
matr_cos_rank2 <- matr

for (i in (1:N)){
  for (j in (1:M)){
    matr[i,j] <- cos(2 * pi * omega1 * i) * cos(2 * pi * omega2 *j) + rnorm(1,sd=0.01)
  }
}

matr_cos_rank4 <- matr

plot2d(matr_cos_rank2)
plot2d(matr_cos_rank4)

s_rank4 <- ssa(matr_cos_rank4, kind = "2d-ssa", L = c(50, 50))
plot(s_rank4, type = "vectors", cuts = 255, layout = c(8, 4))

g <- draft.grouping.auto(s_rank4, grouping.method="low.freq.2dssa", freq.bins1 = 0.15,freq.bins2 = 0.15,
                         threshold = 0.7)
print(g$g)
r <- reconstruct(s_rank4, groups=list(g$g))
plot2d(r$F1)


s_rank2 <- ssa(matr_cos_rank2, kind = "2d-ssa", L = c(50, 50))
plot(s_rank2, type = "vectors", cuts = 255)

g <- draft.grouping.auto(s_rank2, grouping.method="low.freq.2dssa", freq.bins1 = 0.1,freq.bins2 = 0.1,
                         threshold = 0.4)
print(g$g)
r <- reconstruct(s_rank4, groups=list(g$g))
plot2d(r$F1)

