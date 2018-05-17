#####################################################################
# load libraries and set the values of the parameters that are used #
#####################################################################

source('2noise_all_function_1D.R')

library(Rssa)
library(lattice)
library(zoo)
library(dplyr)
library(ggplot2)
library(BioSSA)
library(pracma)
library(fields)
library(pracma)


# number of simulation repeats
M <- 200

# params for 1D
params_1D <- list(
  # what percentage is cut in the middle
  cut_per = 0.1,
  # for SSA
  L = 30,
  # for automatic identification
  threshold = 0.4,
  freq = 0.04,
  # for LOESS
  degree = 2,
  alpha = 0.1,
  # for moving median
  k_med=51
)

# params for 2D
xlim <- c(20, 80)
ylim <- c(20, 80)
# for SSA
step <- 0.5
size = 5
L <- c(round(size/step), round(size/step))
# number of components
r_4 <- 6  # for age 4
r_8 <- 6  # for age 8


