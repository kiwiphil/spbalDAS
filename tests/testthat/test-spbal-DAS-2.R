# Validate DAS functions, features and parameter validation.
# test-spbal-DAS-2.R

testthat::test_that("2. Verify native_lapjv #2.", {
  # humans are rows (3), dogs are columns (3)
  cost <- rbind(c(1, 2, 3),
                c(2, 3, 4),
                c(0, 1, 2))
  res <- native_lapjv(cost)
  # Throw an error if solver computes incorrect cost
  testthat::expect_equal(res$cost, 6)
  # Throw an error if solver computes incorrect pairs
  testthat::expect_equal(res$assignment, c(1, 2, 0))
})
