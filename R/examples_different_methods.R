library(lattice)
library(Rssa)

N <- 199
omega1 <- 0.1
omega2 <- 0.25

set.seed(11)
X <- exp(0.01 * (1:N)) + 2*cos(2 * pi * omega1 * (1:N)) +  cos(2 * pi * omega2 * (1:N)) + rnorm(N,1)

plot(X, type="l")

d <- data.frame(X = X, N = 1:N)
xyplot(X ~ N, data = d, type ='l')

s <- ssa(X)
plot(s, type="vectors")

g <- grouping.auto(
  s,
  grouping.method = 'pgram',
  freq.bins = list(0.01)
  ,
  threshold = 0.9
)

plot(g)

r <- reconstruct(s, groups = g)
d2 <- data.frame(T = r$F1, X = X, N = 1:N)
xyplot(X + T ~ N, data = d2, type ='l')



lst <- grouping.auto(s, grouping.method = "wcor", nclust=2)

print(lst)
plot(lst)
# Check separability
w <- wcor(s)
plot(w)



aa <- auto.idenf.e_m.pgram(s, rho_0 = 0.82, s_0=0.99)

d2 <- data.frame(rho=sort(aa$rho, decreasing = TRUE), N = 1:length(aa$rho))

xyplot(rho ~ N, data = d2, xlab = 'Component')

r <- reconstruct(s, groups = list(2:5))
d3 <- data.frame(P = r$F1, N = 1:N)
xyplot(P ~ N, data = d3, type ='l')

angle.fun <- function(P,Q){
  angle <- function(P1,P2,Q1,Q2){
    acos((P1*P2 + Q1*Q2)/sqrt(P1^2+Q1^2)/sqrt(P2^2+Q2^2))
  }
  var(angle(P[-length(P)],P[-1],Q[-length(Q)],Q[-1]))
}

tau1 <- numeric(ncol(s$U) - 1)
for (j in (1:(ncol(s$U) - 1))){
  tau1[j] <- angle.fun(s$U[,j], s$U[,j+1])
}

d3 <- data.frame(tau=sort(tau1), N = 1:length(tau1))

xyplot(tau ~ N, data = d3, xlab = 'Component')


plot(wcor(s, groups = 1:30), scales = list(at = c(10, 20, 30)))




###### CSSA
X <-  2*cos(2 * pi * omega1 * (1:N)) + 1.0i * cos(2 * pi * omega2 * (1:N)) + rnorm(N,1)
s <- ssa(X, kind='cssa')
plot(s, type='paired')
AI_em_garm_tau_cssa(s, numcomp2 = 3, d1=FALSE)






