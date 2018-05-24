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
trellis.par.get("fontsize")
ecruos("b_ssa_analysis.R")
