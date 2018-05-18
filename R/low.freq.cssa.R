pgram_c <- function(x) {
  if (!is.matrix(x)) x <- as.matrix(x)
  stopifnot(all(is.finite(x)))
  
  X <- mvfft(x)
  n <- nrow(x)
  
  N <- n %/% 2 + 1
  spec <- abs(X[seq_len(N),, drop = FALSE])^2
  
  if (n %% 2 == 0) {
    if (N > 2) spec[2:(N-1), ] <- 2 * spec[2:(N-1), ]
  } else {
    if (N >= 2) spec[2:N, ] <- 2 * spec[2:N, ]
  }
  
  freq <- seq(0, 1, length.out = n + 1)[seq_len(N)]
  
  cumspecfuns <- lapply(seq_len(ncol(x)),
                        function(j)
                          approxfun(c(0, freq[-N] + 1/(2*n), 0.5),
                                    c(0, cumsum(spec[, j])),
                                    rule = 2))
  
  list(spec = spec, freq = freq, cumspecfuns = cumspecfuns)
}

draft.grouping.auto.low.freq.cssa <- function(x, groups,
                          base = c("series", "eigen", "factor"),
                          freq.bins = 2,
                          threshold = 0,
                          method = c("constant", "linear"),
                          ...,
                          drop = TRUE) {
  
  base <- match.arg(base)
  method <- match.arg(method)
  
  if (missing(groups))
    groups <- as.list(1:min(nsigma(x), nu(x)))
  
  groups <- sort(unique(unlist(groups)))
  
  if (identical(base, "eigen")) {
    Fs <- .U(x)[, groups, drop = FALSE]
  } else if (identical(base, "factor")) {
    Fs <- calc.v(x, groups, ...)
  } else if (identical(base, "series")) {
    N <- x$length
    Fs <- matrix(unlist(reconstruct(x, groups = as.list(groups), ...)), nrow = N)
  }
  
  pgs_re <- pgram_c(Re(Fs))
  pgs_im <- pgram_c(Im(Fs))
  
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
  
  norms_re <- colSums(pgs_re$spec)
  contributions_re <- matrix(NA, length(groups), nresult)
  for (i in seq_len(nresult)) {
    contributions_re[, i] <- if (identical(method, "constant"))
      colSums(pgs_re$spec[pgs_re$freq < freq.upper.bound[i] & pgs_re$freq >= freq.lower.bound[i],, drop = FALSE]) / norms_re
    else if (identical(method, "linear"))
      sapply(pgs_re$cumspecfuns, function(f) diff(f(c(freq.lower.bound[i], freq.upper.bound[i])))) / norms_re
  }
  
  norms_im <- colSums(pgs_im$spec)
  contributions_im <- matrix(NA, length(groups), nresult)
  for (i in seq_len(nresult)) {
    contributions_im[, i] <- if (identical(method, "constant"))
      colSums(pgs_im$spec[pgs_im$freq < freq.upper.bound[i] & pgs_im$freq >= freq.lower.bound[i],, drop = FALSE]) / norms_im
    else if (identical(method, "linear"))
      sapply(pgs_im$cumspecfuns, function(f) diff(f(c(freq.lower.bound[i], freq.upper.bound[i])))) / norms_im
  }
  
  contributions <- matrix(0, nrow = nrow(contributions_re), ncol = ncol(contributions_re))
  
  for (i in (1:nrow(contributions_re))){
    for (j in (1:ncol(contributions_re))){
      contributions[i,j] <- max(contributions_re[i,j], contributions_im[i,j])
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
  
  result
}
