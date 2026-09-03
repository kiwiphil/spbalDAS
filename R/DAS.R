
#' @name DAS
#'
#' @title Dynamic Assignment Sampling (DAS).
#'
#' @description DAS draws well-spread master samples over an auxiliary space
#' (Robertson, Price and Reale 2024). Each row of the returned matrix is a
#' candidate sample of size \code{n}; the first \code{k} columns of any row
#' are a size-\code{k} DAS sample.
#'
#' @author This function was first written by Blair Robertson in Matlab and later re-written
#' in R/C++ by Phil Davies.
#'
#' @param pop Numeric matrix of population coordinates, one unit per row.
#' @param n The number of sample points to draw from the population.
#' @param J1 Maximum number of candidate samples, \code{2 <= J1 <= floor(N/2)}.
#'   Default is \code{min(200, floor(N/2))}. Values above 500 are allowed but
#'   trigger a warning (large assignment problems dominate run time).
#' @param n_threads Threads used when filling the cost matrix. \code{0} uses
#'   every logical CPU. Parallel fill runs when \eqn{J^2 i} is at least 2048
#'   (so \code{J1 = 250} does use extra cores). Default is 1.
#' @param verbose Print step timings. Default is FALSE.
#' @param cache_W If \code{TRUE}, build one \eqn{N \times N} inverse-distance
#'   matrix in C++. If \code{FALSE}, compute inverse distances on the fly.
#'   Default is \code{TRUE} when \code{N <= 2500}, otherwise \code{FALSE}.
#'   At large \code{N} the cached matrix is RAM-bound and usually slower.
#'
#' @return An integer matrix with one candidate sample per row (1-based
#'   population indices).
#'
#' @examples
#' set.seed(511)
#' pop <- matrix(runif(200), ncol = 2)
#' sampMx <- DAS(pop = pop, n = 15)
#' sampMx[1:3, ]
#'
#' @export
DAS <- function(pop,
                n,
                J1 = NULL,
                n_threads = 1,
                verbose = FALSE,
                cache_W = NULL) {

  if (is.null(pop) || !is.matrix(pop) || !is.numeric(pop)) {
    stop("spbalDAS(DAS) pop must be a numeric matrix with one unit per row.")
  }
  if (ncol(pop) < 1) {
    stop("spbalDAS(DAS) pop must have at least one coordinate column.")
  }

  N <- nrow(pop)
  if (n < 1 || n >= N) {
    stop("spbalDAS(DAS) n must satisfy 1 <= n < nrow(pop).")
  }

  j_max <- max(2L, as.integer(N %/% 2L))
  if (is.null(J1)) {
    J1 <- min(200L, j_max)
  } else {
    J1 <- as.integer(J1)
  }
  if (J1 < 2L) {
    stop("spbalDAS(DAS) J1 must be >= 2.")
  }
  if (J1 > j_max) {
    #if (verbose) {
      message(sprintf("spbalDAS(DAS) J1=%d exceeds floor(N/2)=%d; using %d.",
                      J1, j_max, j_max))
    #}
    J1 <- j_max
  }
  if (J1 > 500L) {
    warning("spbalDAS(DAS) J1 > 500; large assignment problems dominate run time. ",
            "See Robertson, Price and Reale (2025, Environmetrics).")
  }

  if (is.null(cache_W)) {
    cache_W <- (N <= 2500L)
  }
  cache_W <- isTRUE(cache_W)
  if (cache_W && N > 4000L) {
    warning(sprintf(
      "spbalDAS(DAS) cache_W=TRUE at N=%d stores an N x N matrix (~%.0f MiB) and is usually slower than cache_W=FALSE.",
      N, 8 * N * N / (1024 * 1024)
    ), call. = FALSE)
  }

  #if (verbose) {
  message(sprintf("[%s] Starting DAS N=%d n=%d J1=%d cache_W=%s n_threads=%s",
                    format(Sys.time(), "%Y-%m-%d %H:%M:%OS3"),
                    N, n, J1, cache_W, n_threads))
  #}

  SampleMatrix <- LinearAssignmentProblem_cpp(
    pop = pop,
    target_n = as.integer(n),
    initial_J = as.integer(J1),
    n_threads = as.integer(n_threads),
    verbose = verbose,
    cache_W = cache_W
  )

  #if (verbose) {
  message(sprintf("[%s] Finished DAS.", format(Sys.time(), "%Y-%m-%d %H:%M:%OS3")))
  #}

  attr(SampleMatrix, "spbal") <- "DAS"
  SampleMatrix
}
