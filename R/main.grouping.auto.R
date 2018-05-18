library(Rssa)

source('pair.freq.1dssa.R')
source('low.freq.2dssa.R')
source('tau.1dssa.R')
source('tau.cssa.R')
source('low.freq.cssa.R')


draft.grouping.auto <- function(x, ...,
                          grouping.method = c("pair.freq.1dssa", "tau.1dssa","tau.cssa",
                                              "low.freq.2dssa","low.freq.cssa")) {
  switch(match.arg(grouping.method),
         pair.freq.1dssa = draft.grouping.auto.pair.freq.1dssa(x, ...),
         tau.1dssa  = draft.grouping.auto.tau.1dssa(x, ...),
         tau.cssa  = draft.grouping.auto.tau.cssa(x, ...),
         low.freq.2dssa  = draft.grouping.auto.low.freq.2dssa(x, ...),
         low.freq.cssa  = draft.grouping.auto.low.freq.cssa(x, ...))
}
