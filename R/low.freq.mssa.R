draft.grouping.auto.low.freq.mssa <- function(x, groups,
                                              base = c("series", "eigen", "factor"),
                                              freq.bins = list(0.01),
                                              threshold = 0.9,
                                              method = c("constant", "linear"),
                                              ...,
                                              drop = TRUE) {
  
  base <- match.arg(base)
  method <- match.arg(method)
  
  if (missing(groups))
    groups <- as.list(1:min(nsigma(x), nu(x)))
  
  groups <- sort(unique(unlist(groups)))
  
  if (!is.list(freq.bins)) {
    if (length(freq.bins) == 1 && freq.bins >= 2) {
      freq.bins <- seq(0, 0.5, length.out = freq.bins + 1)[-1]
    }
    
    freq.lower.bound <- c(-Inf, head(freq.bins, -1))
    freq.upper.bound <- freq.bins
  } else {
    freq.bins <- lapply(freq.bins, function(x) if (length(x) == 1) c(-Inf, x) else x)
    freq.lower.bound <- lapply(freq.bins, function(x) x[1])
    freq.upper.bound <- lapply(freq.bins, function(x) x[2])
  }
  
  nresult <- max(length(threshold),
                 length(freq.bins))
  
  if (length(freq.lower.bound) < nresult) {
    freq.lower.bound <- rep_len(freq.lower.bound, nresult)
  }
  
  if (length(freq.upper.bound) < nresult) {
    freq.upper.bound <- rep_len(freq.upper.bound, nresult)
  }
  
  if (length(threshold) < nresult) {
    threshold <- rep_len(threshold, nresult)
  }
  
  
  if (identical(base, "eigen")) {
    # Fs <- .U(x)[, groups, drop = FALSE]
    Fs <- s$U[, groups]
    pgs <- pgram_1d(Fs)
    
    norms <- colSums(pgs$spec)
    contributions <- matrix(NA, length(groups), nresult)
    for (i in seq_len(nresult)) {
      contributions[, i] <- if (identical(method, "constant"))
        colSums(pgs$spec[pgs$freq < freq.upper.bound[i] & pgs$freq >= freq.lower.bound[i],, drop = FALSE]) / norms
      else if (identical(method, "linear"))
        sapply(pgs$cumspecfuns, function(f) diff(f(c(freq.lower.bound[i], freq.upper.bound[i])))) / norms
    }
    
  } else if (identical(base, "factor")) {
    Fs <- calc.v(x, groups, ...)
    Fs_list <- list()
    K <- cumsum(c(1,x$length - x$window + 1))
    for (i in 1:(length(K)-1)){
      Fs_list[[i]] <- Fs[K[i]:(K[i+1] -1),]
    }
    
    pgram.all <- lapply(Fs_list, pgram_1d)
    
    contrib <- function(pgs){
      norms <- colSums(pgs$spec)
      contributions <- matrix(NA, length(groups), nresult)
      for (i in seq_len(nresult)) {
        contributions[, i] <- if (identical(method, "constant"))
          colSums(pgs$spec[pgs$freq < freq.upper.bound[i] & pgs$freq >= freq.lower.bound[i],, drop = FALSE]) / norms
        else if (identical(method, "linear"))
          sapply(pgs$cumspecfuns, function(f) diff(f(c(freq.lower.bound[i], freq.upper.bound[i])))) / norms
      }
      contributions
    }
    
    contr.all <- lapply(pgram.all, contrib)
    
    get.element <- function(x, i, j){
      unlist(lapply(x, function(y) y[i,j]))
    }
    
    n <- nrow(contr.all[[1]])
    m <- ncol(contr.all[[1]])
    contributions <- matrix(0, nrow = n, ncol = m)
    for (i in (1:n)){
      for (j in (1:m)){
        contributions[i,j] <- max(get.element(contr.all, i, j))
      }
    }
    
  }
  else if (identical(base, "series")) {
    Fs <- reconstruct(x, groups = as.list(groups), ...)
    
    Fs_list <- list()
    for (j in (1:length(Fs[[1]]))){
      Fs_list[[j]] <- (Fs[[1]][[j]])
    }
    
    for (k in (2:length(Fs))){
      for (j in (1:length(Fs[[1]]))){
        Fs_list[[j]] <- cbind(Fs_list[[j]], Fs[[k]][[j]])
      }
    }
    
    pgram.all <- lapply(Fs_list, pgram_1d)
    
    contrib <- function(pgs){
      norms <- colSums(pgs$spec)
      contributions <- matrix(NA, length(groups), nresult)
      for (i in seq_len(nresult)) {
        contributions[, i] <- if (identical(method, "constant"))
          colSums(pgs$spec[pgs$freq < freq.upper.bound[i] & pgs$freq >= freq.lower.bound[i],, drop = FALSE]) / norms
        else if (identical(method, "linear"))
          sapply(pgs$cumspecfuns, function(f) diff(f(c(freq.lower.bound[i], freq.upper.bound[i])))) / norms
      }
      contributions
    }
    
    contr.all <- lapply(pgram.all, contrib)
    
    get.element <- function(x, i, j){
      unlist(lapply(x, function(y) y[i,j]))
    }
    
    n <- nrow(contr.all[[1]])
    m <- ncol(contr.all[[1]])
    contributions <- matrix(0, nrow = n, ncol = m)
    for (i in (1:n)){
      for (j in (1:m)){
        contributions[i,j] <- max(get.element(contr.all, i, j))
      }
    }
  }
  
  type <- if (all(threshold <= 0)) "splitting" else "independent"
  
  if (identical(type, "splitting")) {
    gi <- max.col(contributions, ties.method = "first")
    result <- lapply(seq_len(nresult), function(i) groups[gi == i])
  } else if (identical(type, "independent")) {
    result <- lapply(seq_len(nresult), function(i) groups[contributions[, i] >= threshold[i]])
  } else {
    stop("Unknown type for pgrouping.auto.pgram")
  }
  
  # names(result) <- if (!is.null(names(freq.bins))) .group.names(freq.bins) else .group.names(threshold)
  # colnames(contributions) <- names(result)
  # rownames(contributions) <- as.character(groups)
  
  names(result) <- "F1"
  colnames(contributions) <- names(result)
  rownames(contributions) <- as.character(groups)
  
  if (drop) {
    result <- result[sapply(result, length) > 0]
  }
  
  attr(result, "contributions") <- contributions
  attr(result, "type") <- type
  attr(result, "threshold") <- threshold
  
  class(result) <- "grouping.auto.pgram"
  
  return <- result
}
