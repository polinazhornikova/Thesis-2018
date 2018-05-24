library(methods)

fragmentStart <- function(file) {
  cat(sprintf("@@@fragmentStart\n%s\n", file))
}

fragmentStop <- function() {
  cat(sprintf("@@@fragmentStop\n"))
}

fragmentSkip <- function(exp) {
  capture.output(eval(substitute(exp), envir = parent.frame()))
  cat(sprintf("@@@fragmentSkip\n"))
}

ecruos <- function(x, ..., remove.comment.lines = FALSE) {
  lines <- capture.output(source(x, ..., echo = TRUE, max.deparse.length = Inf, keep.source = TRUE))
  lines <- lines[lines != ""]

  skip.start <- grep("fragmentSkip\\s*\\(", lines)
  skip.stop <- grep("^@@@fragmentSkip$", lines)

  stopifnot(length(skip.start) == length(skip.stop))
  if (length(skip.start) > 0) {
    skips <- lapply(seq_along(skip.start), function(i) seq(skip.start[i], skip.stop[i]))
    skips <- do.call("c", skips)
    lines <- lines[-skips]
  }

  fragmentStart.calls <- grep("^@@@fragmentStart$", lines) - 1
  fragmentStop.calls <- grep("^@@@fragmentStop$", lines) - 1
  stopifnot(length(fragmentStart.calls) == length(fragmentStop.calls))
  if (length(fragmentStart.calls) > 0) {
    lines <- lines[-c(fragmentStart.calls, fragmentStop.calls)]

    fragment.start <- grep("^@@@fragmentStart$", lines)
    fragment.stop <- grep("^@@@fragmentStop$", lines)

    fragment.names <- lines[fragment.start + 1]
    for (i in seq_along(fragment.names)) {
      cur.lines <- lines[seq(fragment.start[i] + 2, fragment.stop[i] - 1)]
      if (remove.comment.lines) {
        comment.lines <- grepl("^\\s*#.*$", cur.lines)
        cur.lines <- cur.lines[!comment.lines]
      }

      cur.lines <- c("\\begin{CodeChunk}",
                     "\\begin{CodeInput}",
                     "\n",
                     cur.lines,
                     "\\end{CodeInput}\n",
                     "\\end{CodeChunk}")
      cat(file = fragment.names[i], paste(cur.lines, collapse = "\n"))
    }
  }
}

ecruos.all <- function(..., ignore.list = c("ecruos.r", "r2txt_frag.r", "run.r")) {
  R.names <- dir(pattern = "^.*\\.[rR]$")
  R.names <- R.names[!R.names %in% ignore.list]
  for (name in R.names) {
    cat("Processing file: ", name, "\n")
    ecruos(name)
  }
}

regex.extract.group <- function(pattern, x, group = 1) {
  match <- gregexpr(pattern, x, perl = TRUE)
  match <- lapply(match, function(x) {
      x[] <-attr(x, "capture.start")[group, ]
      attr(x, "match.length") <- attr(x, "capture.length")[group, ]
      attr(x, "capture.start") <- attr(x, "capture.names") <- attr(x, "capture.names") <- NULL
      x
    })
  regmatches(x, match)
}

extract.fragments <- function(infile, outfile = "",
                              dropEmptyLines = TRUE,
                              dropComments = FALSE,
                              drop3Comments = TRUE,
                              dropOuter = FALSE,
                              dropSkips = TRUE) {
  lines <- readLines(infile)

  if (dropEmptyLines)
    lines <- lines[lines != ""]

  if (dropComments)
    lines <- lines[!grepl("^\\s*#", lines)]

  if (drop3Comments)
    lines <- lines[!grepl("^\\s*###", lines)]

  skips <- grepl("\\s*fragmentSkip\\s*\\(.*)", lines)
  if (dropSkips) {
    lines <- lines[!skips]
  } else {
    lines[skips] <- unlist(regex.extract.group("\\s*fragmentSkip\\s\\(\\s*(.*)\\s*\\)", lines))
  }

  fragment.start.regex <- "\\s*fragmentStart\\s*\\(\\s*\"(.*)\"\\s*\\)"
  fragment.start <- grep(fragment.start.regex, lines)
  fragment.stop <- grep("\\s*fragmentStop\\s*\\(*\\)", lines)

  fragment.names <- unlist(regex.extract.group(fragment.start.regex, lines))
  fragment.names <- unlist(regex.extract.group("\\.\\./fragments/(.*)\\.tex", lines))

  stopifnot(length(fragment.start) == length(fragment.stop))
  stopifnot(length(fragment.start) == length(fragment.names))

  for (i in seq_along(fragment.names)) {
    lines[fragment.start[i]] <- sprintf("### Begin fragment: %s", fragment.names[i])
    lines[fragment.stop[i]] <- sprintf("### End fragment: %s\n", fragment.names[i])
  }

  fragment.idx <- unlist(lapply(seq_along(fragment.names), function(i) {
      seq(fragment.start[i], fragment.stop[i])
    }))

  if (dropOuter)
    lines <- lines[fragment.idx]

  res <- paste0(lines, collapse = "\n")

  cat(res, file = outfile)
}
