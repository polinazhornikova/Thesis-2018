library("Rssa")
library("lattice")
source('main.grouping.auto.R')

# Create a zero matrix
new.matrix <- function(N_t = 100, N_c = 100){
  matrix(0, nrow = N_c, ncol = N_t)
}

# Adding lines in matrix
add.line <- function(m, a = 1 , b = 0){
  N_c <- nrow(m)
  N_t <- ncol(m)
  k <- 1
  for (k in 1:N_t){
    if (a*k + b >= 1 && a*k + b <= N_c){
      m[a*k + b, k] <- 1
      k <- k+1
    }
  }
  m
}

# Adding noise in matrix
add.noise <- function(m, sigma = 0.2){
  N_c <- nrow(m)
  N_t <- ncol(m)
  for (i in (1:N_t)){
    m[,i] <- m[,i] + rnorm(N_c, sd = sigma)
  }
  m
}

# Step DFT
dft <- function(m){
  N_c <- nrow(m)
  N_t <- ncol(m)
  for (i in (1:N_c)){
    m[i,] <- fft(m[i,])
  }
  m
}

# Step Complex SSA by rows
# You can use the Threshold: step = TRUE
# eps -- value threshold
# functional -- used functional
cssa.row <- function(m, num.line = 1, step = FALSE, eps = 0.001, 
                    functional = angle.fun){
  N_c <- nrow(m)
  N_t <- ncol(m)
  for (i in (1:N_c)){
    s <- ssa(m[i,], kind = "cssa", svd.method = "svd")
    if (!step){
      r <- reconstruct(s, groups = list(Seasonality = 1:num.line))
      m[i,] <- r$Seasonality
    } 
    else{
      v <- numeric(N_t)
      for (j in 1:num.line){
        if (functional(s$U[,j]) < eps){  
        r <- reconstruct(s, groups = list(Seasonality = j))
        v <- v + r$Seasonality
        }
        m[i,] <- v
      }
    }
  }
  m
}


# Step IDFT by rows
idft.row <- function(m){
  N_c <- nrow(m)
  N_t <- ncol(m)
  for (i in (1:N_c)){
    m[i,] <- fft(m[i,], inverse=TRUE)/length(m[i,])
  }
  m
}

# Step IDFT by columns
idft.col <- function(m){
  N_c <- nrow(m)
  N_t <- ncol(m)
  for (i in (1:N_t)){
    m[,i] <- fft(m[,i], inverse=TRUE)/length(m[,i])
  }
  m
}


# Drawing matrix
# If from.0.to.1 = TRUE, all values less than 0, replaced by 0,
# all values greater than 1, are replaced by 1

plot.matrix <- function(m, from.0.to.1 = FALSE){
  rgb.palette <- colorRampPalette(c("white", "black"), space = "rgb")
  if (!from.0.to.1){ 
    levelplot(m, xlab="", ylab="", col.regions=rgb.palette,  
              scales = list(draw=FALSE), auto.key= FALSE, 
              colorkey = FALSE)
  }
  else {
    for (i in 1:nrow(m)){
      for (j in 1:ncol(m)){
        if (m[i,j] > 1) {
          m[i,j] <- 1
        }
      }
    }
    levelplot(m, xlab="", ylab="", col.regions=rgb.palette,  
              scales = list(draw=FALSE), auto.key= FALSE, 
              at=seq(0,1,0.01), colorkey =FALSE)
  }
}

# Step Complex SSA by columns
cssa.col <- function(m, num.line = 1, auto=1){
  N_c <- nrow(m)
  N_t <- ncol(m)
  for (i in (1:N_t)){
    s <- ssa(m[,i], kind = "cssa", svd.method = "svd")
    if (auto){
      idx <- draft.grouping.auto(x=s, grouping.method = "tau.cssa",numcomp1 = num.line)$d1
      r <- reconstruct(s, groups = list(Seasonality = (idx)))
    }
    else{
      r <- reconstruct(s, groups = list(Seasonality = (1:num.line)))
    }
    m[,i] <- r$Seasonality
  }
  m
}

# Rotation matrix - for the correct display of matrices
rotate <- function(x) t(apply(x, 2, rev))

# Example
N_t <- 120  # length of series
N_c <- 100  # number of series

m.in <- new.matrix(N_t, N_c)
# m.in <- add.line(m.in, b = 10)
m.in <- add.line(m.in, b=-10, a = 1)
m.in <- add.line(m.in, b=110, a = -1)
# m.in <- add.line(m.in, b= 70, a = 0)
m.noise <- add.noise(m.in, 0.2)

plot.matrix(rotate(Re(m.noise)))

m.dft <- dft(m.noise)
m.cssa.col <- cssa.col(m.dft, num.line = 2)
m.col.row <- idft.row(m.cssa.col)
m.cssa.col2 <- cssa.col(m.dft, num.line = 2, auto=0)
m.col.row2 <- idft.row(m.cssa.col2)

plot.matrix(rotate(Re(m.noise)))
plot.matrix(rotate(Re(m.col.row)), from.0.to.1 = TRUE)
plot.matrix(rotate(Re(m.col.row2)), from.0.to.1 = TRUE)
