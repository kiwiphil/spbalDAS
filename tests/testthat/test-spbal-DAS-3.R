# Validate DAS functions, features and parameter validation.
# test-spbal-DAS-3.R

testthat::test_that("3. Verify arma_dist_al #1.", {
  # calculate euclidean distance between points.
  A <- rbind(c(1, 2, 3),
             c(2, 3, 4),
             c(0, 1, 2))
  A_dist <- rbind(c(0.00000000, 1.73205081, 1.73205081),
                  c(1.73205081, 0.00000000, 3.46410162),
                  c(1.73205081, 3.46410162, 0.00000000))
  dist <- as.matrix(arma_dist_al(A))
  # Throw an error if solver computes incorrect distance matrix.
  testthat::expect_equal(dist, A_dist)
})

testthat::test_that("3. Verify arma_dist_al #2.", {
  # calculate euclidean distance between points.
  B <- rbind(c(1, 2, 3),
             c(4, 3, 2))
  B_dist <- rbind(c(0.00000000, 3.31662479),
                  c(3.31662479, 0.00000000))
  dist <- as.matrix(arma_dist_al(B))
  # Throw an error if solver computes incorrect distance matrix.
  testthat::expect_equal(dist, B_dist)
})

testthat::test_that("3. Verify arma_dist_al #3.", {
  # calculate euclidean distance between points.
  C <- matrix(c(2, 3, 0, 9, 4, 5), nrow = 3, byrow = TRUE)
  C_dist <- rbind(c(0.00000000, 6.32455532, 2.82842712),
                  c(6.32455532, 0.00000000, 5.65685425),
                  c(2.82842712, 5.65685425, 0.00000000))
  dist <- as.matrix(arma_dist_al(C))
  # Throw an error if solver computes incorrect distance matrix.
  testthat::expect_equal(dist, C_dist)
})



