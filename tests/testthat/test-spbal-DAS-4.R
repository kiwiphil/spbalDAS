# DAS driver and cost-path smoke tests.
# test-spbal-DAS-4.R

testthat::test_that("DAS runs without an explicit J1 and returns unique units per row", {
  set.seed(511)
  pop <- matrix(runif(80), ncol = 2)
  samp <- DAS(pop = pop, n = 8L, n_threads = 1L, verbose = FALSE)
  testthat::expect_true(is.matrix(samp))
  testthat::expect_equal(ncol(samp), 8L)
  testthat::expect_gt(nrow(samp), 0L)
  testthat::expect_true(all(samp >= 1L & samp <= nrow(pop)))
  for (r in seq_len(nrow(samp))) {
    testthat::expect_equal(length(unique(samp[r, ])), 8L)
  }
})

testthat::test_that("cache_W TRUE and FALSE both produce valid candidate matrices", {
  set.seed(17)
  pop <- matrix(runif(60), ncol = 2)
  cached <- DAS(pop = pop, n = 6L, J1 = 8L, cache_W = TRUE, verbose = FALSE)
  onthefly <- DAS(pop = pop, n = 6L, J1 = 8L, cache_W = FALSE, verbose = FALSE)
  testthat::expect_equal(ncol(cached), 6L)
  testthat::expect_equal(ncol(onthefly), 6L)
  testthat::expect_true(all(cached >= 1L & cached <= nrow(pop)))
  testthat::expect_true(all(onthefly >= 1L & onthefly <= nrow(pop)))
})

testthat::test_that("DAS rejects n >= N", {
  pop <- matrix(runif(20), ncol = 2)
  testthat::expect_error(DAS(pop = pop, n = 10L), "1 <= n < nrow")
})
