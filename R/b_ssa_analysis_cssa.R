###Fragment cssa_series
fragmentStart("fragments_ch1/cssa_series.tex")
x <- exp(0.01 * (1:N)) + cos(2 * pi * omega1 * (1:N)) + 
  1.0i*sin(2 * pi * omega1 * (1:N)) + 
  1.0i*exp(0.007 * (1:N))*cos(2 * pi * omega2 * (1:N)) + rnorm(N, 1) 
s <- ssa(x, kind = 'cssa', L = 100)
fragmentSkip(pdf("img_ch1/cssa_vectors.pdf", paper = "special", width = 5, height = 2))
plot(s, type="vectors", layout = c(5, 2))
fragmentSkip(dev.off())
fragmentStop()
###end 


###Fragment cssa_trend
fragmentStart("fragments_ch1/cssa_trend.tex")
g_trend <- general.grouping.auto(s, grouping.method = 'low.freq.cssa', 
                                 freq.bins = list(0.01), threshold = 0.9)
print(g_trend$F1)
fragmentStop()
###end 


###Fragment cssa_em_freq
fragmentStart("fragments_ch1/cssa_em_freq.tex")
g_em_freq <- general.grouping.auto(s, grouping.method = "freq.cssa", s_0 = 1, 
                                   rho_0 = 0.95)
print(g_em_freq)
fragmentStop()
###end 

###Fragment cssa_em_tau
fragmentStart("fragments_ch1/cssa_em_tau.tex")
g_em_tau <- general.grouping.auto(s, grouping.method = "tau.cssa", treshold = 0.01)
print(g_em_tau$idx)
fragmentStop()
###end 

###Fragment cssa_recon
fragmentStart("fragments_ch1/cssa_rec.tex")
r <- reconstruct(s, groups = list(T = g_trend, P = g_em_tau$idx))
fragmentSkip(pdf("img_ch1/cssa_rec.pdf", paper = "special", width = 4, height = 3))
d_re <- data.frame(T_re = Re(r$T), P_re = Re(r$P),  X_re = Re(x), N = 1:N)
xyplot(T_re + P_re + X_re  ~ N, data = d_re, type ='l', ylab = '',
       auto.key = list(points = FALSE, lines = TRUE))
fragmentSkip(dev.off())
fragmentStop()
###end 
