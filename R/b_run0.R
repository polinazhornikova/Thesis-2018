source("ecruos.R")

library(lattice)
# Setup various margins
lattice.options(default.theme =
                  modifyList(standard.theme("pdf", color = TRUE),
                             list(
                               fontsize = list(text = 7, points = 5),
                               layout.heights =
                                 list(top.padding = 0.1,
                                      key.axis.padding = 0.25,
                                      axis.xlab.padding = 0.1,
                                      xlab.key.padding = 0.1,
                                      key.sub.padding = 0.25,
                                      bottom.padding = 0.1),
                               layout.widths =
                                 list(left.padding = 0.1,
                                      key.ylab.padding = 0.1,
                                      ylab.axis.padding = 0.1,
                                      axis.key.padding = 0.1,
                                      right.padding = 1))))
library(ssabook)
library(Rssa)
plot.d1 <- function(s,index=1:10){
  index <- unique(index)
  u <- s$U[,index]
  data <- data.frame(z = as.vector(u), 
                     g = factor(rep(seq_len(ncol(u)),                     
                                    each = nrow(u)),
                                labels = paste("U", seq_len(ncol(u)))))
  xyplot(Im(z) ~ Re(z) | g, data = data, type = "l", 
         as.table = TRUE,  
         scales = list(relation = "free", draw = FALSE), layout = c(5, 2))
}

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

trellis.par.get("fontsize")
ecruos("b_ssa_analysis_1dssa.R")
