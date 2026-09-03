#' @name doubleDAS
#'
#' @title Doubly balanced Dynamic Assignment Sampling (doubleDAS).
#'
#' @description Port of the MATLAB \code{doubleDAS} sampler. Candidate samples
#' grow by successive linear assignment with cost
#' \deqn{C = \alpha\,SB + (1-\alpha)\,B}{C = alpha * SB + (1 - alpha) * B}
#' where \eqn{SB} is the inverse-distance spread used by \code{\link{DAS}}
#' (MATLAB \code{SBMoranIvec}) and \eqn{B} is the squared sum of auxiliary
#' variables after appending a candidate unit (MATLAB \code{LinearAssignment}).
#'
#' \code{alpha = 1} is spatially balanced (same objective as \code{DAS}).
#' \code{alpha = 0} is approximately balanced on \code{aux}. Values in
#' between are doubly balanced.
#'
#' Reuses \code{arma_dist_al} weights, the DAS candidate-growth schedule, and
#' \code{native_lapjv} (Jonker-Volgenant) in place of MATLAB \code{matchpairs}.
#'
#' @author Blair Robertson (Matlab); R/C++ port by Phil Davies.
#'
#' @param pop Numeric matrix of population coordinates, one unit per row.
#' @param aux Numeric matrix of auxiliary variables, one row per unit.
#' @param n Sample size.
#' @param alpha Mix in \eqn{[0, 1]}. Default 0.5.
#' @param J1 Maximum number of candidate samples,
#'   \code{2 <= J1 <= floor(N/2)}. Default \code{min(500, floor(N/2))},
#'   matching the MATLAB.
#' @param n_threads Threads for the cost-matrix fill. See \code{\link{DAS}}.
#' @param verbose Print step timings. Default FALSE.
#' @param cache_W Cache the \eqn{N \times N} inverse-distance matrix when
#'   \code{alpha > 0}. Default \code{N <= 2500}. At large \code{N} this is
#'   usually slower than on-the-fly distances.
#'
#' @return Integer matrix of candidate samples (1-based population indices).
#'   Each row is a sample of size \code{n}; the first \code{k} columns of any
#'   row are a size-\code{k} doubleDAS sample.
#'
#' @examples
#' set.seed(511)
#' pop <- matrix(runif(200), ncol = 2)
#' aux <- matrix(runif(300), ncol = 3)
#' samp <- doubleDAS(pop, aux, n = 12, alpha = 0.5)
#' samp[1, ]
#'
#' @export
doubleDAS <- function(pop,
                      aux,
                      n,
                      alpha = 0.5,
                      J1 = NULL,
                      n_threads = 1,
                      verbose = FALSE,
                      cache_W = NULL) {

  if (is.null(pop) || !is.matrix(pop) || !is.numeric(pop)) {
    stop("spbalDAS(doubleDAS) pop must be a numeric matrix with one unit per row.")
  }
  if (ncol(pop) < 1L) {
    stop("spbalDAS(doubleDAS) pop must have at least one coordinate column.")
  }
  if (is.null(aux) || !is.matrix(aux) || !is.numeric(aux)) {
    stop("spbalDAS(doubleDAS) aux must be a numeric matrix with one row per unit.")
  }
  if (nrow(aux) != nrow(pop)) {
    stop("spbalDAS(doubleDAS) aux must have the same number of rows as pop.")
  }
  if (ncol(aux) < 1L) {
    stop("spbalDAS(doubleDAS) aux must have at least one column.")
  }

  N <- nrow(pop)
  if (n < 1 || n >= N) {
    stop("spbalDAS(doubleDAS) n must satisfy 1 <= n < nrow(pop).")
  }
  if (!is.numeric(alpha) || length(alpha) != 1L || alpha < 0 || alpha > 1) {
    stop("spbalDAS(doubleDAS) alpha must be a single number in [0, 1].")
  }

  # MATLAB: J1 = min(floor(N/2), 500)
  j_max <- max(2L, as.integer(N %/% 2L))
  if (is.null(J1)) {
    J1 <- min(500L, j_max)
  } else {
    J1 <- as.integer(J1)
  }
  if (J1 < 2L) {
    stop("spbalDAS(doubleDAS) J1 must be >= 2.")
  }
  if (J1 > j_max) {
    #if (verbose) {
      message(sprintf("spbalDAS(doubleDAS) J1=%d exceeds floor(N/2)=%d; using %d.",
                      J1, j_max, j_max))
    #}
    J1 <- j_max
  }
  if (J1 > 500L) {
    warning("spbalDAS(doubleDAS) J1 > 500; large assignment problems dominate run time. ",
            "MATLAB caps J1 at 500. See Robertson, Price and Reale (2025, Environmetrics).")
  }

  if (is.null(cache_W)) {
    cache_W <- (N <= 2500L)
  }
  cache_W <- isTRUE(cache_W)
  if (cache_W && N > 4000L) {
    warning(sprintf(
      "spbalDAS(doubleDAS) cache_W=TRUE at N=%d stores an N x N matrix (~%.0f MiB) and is usually slower than cache_W=FALSE.",
      N, 8 * N * N / (1024 * 1024)
    ), call. = FALSE)
  }

  #if (verbose) {
  message(sprintf(
      "[%s] Starting doubleDAS N=%d n=%d J1=%d alpha=%.3f cache_W=%s",
      format(Sys.time(), "%Y-%m-%d %H:%M:%OS3"),
      N, n, J1, alpha, cache_W
    ))
  #}

  SampleMatrix <- DoubleAssignmentProblem_cpp(
    pop = pop,
    aux = aux,
    alpha = as.numeric(alpha),
    target_n = as.integer(n),
    initial_J = as.integer(J1),
    n_threads = as.integer(n_threads),
    verbose = verbose,
    cache_W = cache_W
  )

  #if (verbose) {
  message(sprintf("[%s] Finished doubleDAS.",
                    format(Sys.time(), "%Y-%m-%d %H:%M:%OS3")))
  #}

  attr(SampleMatrix, "spbal") <- "doubleDAS"
  SampleMatrix
}

#' Inverse-distance spread used by MATLAB SBMoranIvec.
#'
#' \code{2 * sum(W[A, C])} for each candidate in \code{A}. Kept for tests and
#' for comparing against the C++ cost fill.
#'
#' @keywords internal
SBMoranIvec <- function(C, A, W) {
  as.numeric(2 * rowSums(W[A, C, drop = FALSE]))
}

#' MATLAB LinearAssignment balance term for one candidate sample.
#'
#' @keywords internal
doubleDAS_balance_row <- function(Ci, Ai, X) {
  b <- colSums(X[Ci, , drop = FALSE])
  XA <- X[Ai, , drop = FALSE]
  as.numeric(rowSums(sweep(XA, 2L, b, "+")^2))
}
