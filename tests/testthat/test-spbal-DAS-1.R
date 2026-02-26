# Validate DAS functions, features and parameter validation.
# test-spbal-DAS-1.R

testthat::test_that("1. Verify HungarianSolver.", {
  # humans are rows (5), dogs are columns (4)
  cost <- rbind(c(1, 5, 2, 19),
                c(4, 0, 3, 10),
                c(6, 9, 6, 1),
                c(9, 8, 1, 3),
                c(6, 1, 1, 1))
  res <- HungarianSolver(cost)
  # Throw an error if solver computes incorrect cost
  expect_equal(res$cost, 3)
  # Throw an error if solver computes incorrect pairs
  expect_equal(res$pairs, matrix(c(1, 2, 3, 4, 5, 1, 2, 4, 3, 0), nrow = 5, ncol = 2))
})

testthat::test_that("2. Verify rSBMoranIvec.", {
  # Create a 6x6 weighted matrix (floats, symmetric)
  W <- matrix(0, nrow = 6, ncol = 6)

  # Add weighted edges (node pairs with weights)
  edges <- rbind(
    c(1,2, 2.5),   # edge 1-2 weight 2.5
    c(1,3, 1.0),
    c(2,3, 3.0),
    c(3,4, 4.5),
    c(4,5, 10.0),  # strong connection to outside
    c(4,6, 1.5),
    c(5,6, 0.5)    # weak between outside nodes
  )

  for (i in 1:nrow(edges)) {
    a <- edges[i,1]
    b <- edges[i,2]
    w <- edges[i,3]
    W[a,b] <- w
    W[b,a] <- w  # symmetric (undirected graph)
  }

  # Cluster: nodes 1,2,3,4
  C <- c(1, 2, 3, 4)

  # Outside: nodes 5 and 6
  A <- c(5, 6)

  result <- spbalDAS::rSBMoranIvec(C, A, W)
  expect_equal(result, c(42.0, 25.0))
  # Expected: [1] 33 15
})
