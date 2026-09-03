# MATLAB-faithful doubleDAS tests.
# test-spbal-DAS-5.R

testthat::test_that("MATLAB SBMoranIvec is 2 * rowSums(W[A, C])", {
  W <- matrix(c(
    0, 1, 2, 3,
    1, 0, 4, 5,
    2, 4, 0, 6,
    3, 5, 6, 0
  ), 4, 4, byrow = TRUE)
  C <- c(1L, 3L)
  A <- c(2L, 4L)
  # MATLAB: 2 * sum(W(A,C), 2)'
  expected <- 2 * c(sum(W[2, C]), sum(W[4, C]))
  testthat::expect_equal(SBMoranIvec(C, A, W), expected)
})

testthat::test_that("MATLAB balance row matches (b + X[Ai,])^2 summed over columns", {
  X <- matrix(c(
    1, 0,
    2, 1,
    0, 3,
    4, 1
  ), 4, 2, byrow = TRUE)
  Ci <- c(1L, 3L)
  Ai <- c(2L, 4L)
  b <- colSums(X[Ci, , drop = FALSE])
  expected <- c(
    sum((b + X[2, ])^2),
    sum((b + X[4, ])^2)
  )
  testthat::expect_equal(doubleDAS_balance_row(Ci, Ai, X), expected)
})

testthat::test_that("doubleDAS returns unique 1-based indices per row", {
  set.seed(511)
  pop <- matrix(runif(80), ncol = 2)
  aux <- matrix(runif(120), ncol = 3)
  samp <- doubleDAS(pop, aux, n = 8L, alpha = 0.5, J1 = 10L, verbose = FALSE)
  testthat::expect_true(is.matrix(samp))
  testthat::expect_equal(ncol(samp), 8L)
  testthat::expect_true(all(samp >= 1L & samp <= nrow(pop)))
  for (r in seq_len(nrow(samp))) {
    testthat::expect_equal(length(unique(samp[r, ])), 8L)
  }
  testthat::expect_equal(attr(samp, "spbal"), "doubleDAS")
})

testthat::test_that("doubleDAS alpha 0, 0.5 and 1 all run", {
  set.seed(2)
  pop <- matrix(runif(60), ncol = 2)
  aux <- cbind(pop[, 1], runif(nrow(pop)))
  for (a in c(0, 0.5, 1)) {
    samp <- doubleDAS(pop, aux, n = 6L, alpha = a, J1 = 8L, cache_W = TRUE)
    testthat::expect_equal(ncol(samp), 6L)
  }
  onthefly <- doubleDAS(pop, aux, n = 6L, alpha = 0.5, J1 = 8L, cache_W = FALSE)
  testthat::expect_equal(ncol(onthefly), 6L)
})

testthat::test_that("doubleDAS rejects mismatched aux and alpha outside [0,1]", {
  pop <- matrix(runif(20), ncol = 2)
  testthat::expect_error(doubleDAS(pop, aux = matrix(1, 3, 1), n = 3, alpha = 0.5),
                         "same number of rows")
  testthat::expect_error(doubleDAS(pop, aux = matrix(runif(20), ncol = 2), n = 3, alpha = 1.5),
                         "alpha")
})

# MATLAB LinearAssignment.m: C = alpha*SB + (1-alpha)*B, then matchpairs.
testthat::test_that("one MATLAB LinearAssignment step is a min-cost bijection", {
  pop <- matrix(c(
    0, 0,
    1, 0,
    0, 1,
    1, 1,
    0.5, 0.5,
    2, 2,
    3, 0,
    0, 3
  ), ncol = 2, byrow = TRUE)
  X <- cbind(pop, seq_len(nrow(pop)))
  D <- as.matrix(dist(pop))
  W <- 1 / D
  diag(W) <- 0
  W[!is.finite(W)] <- 0

  # Complement of union(Ci) is Ai, as in MATLAB setdiff(1:N, Ci(:)).
  Ci <- rbind(c(1L, 2L), c(3L, 4L), c(1L, 3L))
  Ai <- c(5L, 6L, 7L)
  alpha <- 0.4
  n_samples <- nrow(Ci)

  SB <- matrix(0, n_samples, n_samples)
  B <- matrix(0, n_samples, n_samples)
  for (i in seq_len(n_samples)) {
    SB[i, ] <- SBMoranIvec(Ci[i, ], Ai, W)
    B[i, ] <- doubleDAS_balance_row(Ci[i, ], Ai, X)
  }
  C <- alpha * SB + (1 - alpha) * B
  res <- native_lapjv(C)
  asg <- res$assignment
  assigned <- Ai[asg + 1L]
  testthat::expect_equal(sort(assigned), sort(Ai))

  # Brute-force minimum over the 3! permutations of columns.
  perms <- matrix(c(
    0, 1, 2,
    0, 2, 1,
    1, 0, 2,
    1, 2, 0,
    2, 0, 1,
    2, 1, 0
  ), ncol = 3, byrow = TRUE)
  brute <- apply(perms, 1L, function(p) {
    C[1, p[1] + 1] + C[2, p[2] + 1] + C[3, p[3] + 1]
  })
  testthat::expect_equal(res$cost, min(brute))
})
