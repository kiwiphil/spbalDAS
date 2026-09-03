# Validate DAS functions, features and parameter validation.
# test-spbal-DAS-1.R

testthat::test_that("1. Verify native_lapjv #1.", {
  # humans are rows (3), dogs are columns (3)
  cost <- rbind(c(1, 2, 0),
                c(2, 0, 1),
                c(1, 4, 19))
  res <- native_lapjv(cost)
  # Throw an error if solver computes incorrect cost
  testthat::expect_equal(res$cost, 1)
  # Throw an error if solver computes incorrect pairs
  testthat::expect_equal(res$assignment, c(2, 1, 0))
})


testthat::test_that("1. Verify native_lapjv #2.", {
  # Agents are rows (3), Tasks are columns (3)
  cost <- rbind(c(4, 2, 8),
                c(2, 3, 7),
                c(3, 1, 6))
  res <- native_lapjv(cost)
  # Throw an error if solver computes incorrect cost
  testthat::expect_equal(res$cost, 10)
  # Throw an error if solver computes incorrect pairs
  testthat::expect_equal(res$assignment, c(1, 0, 2))
})

