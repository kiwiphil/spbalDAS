# doubleDAS.R


#' @name doubleDAS
#'
#' @title Double Dynamic Assignment Sampling (doubleDAS).
#'
#' @description DAS
#'
#' @author This function was first written by Blair Robertson in Matlab and later re-written in R by Phil Davies.
#'
#' @param pop Population points.
#' @param aux Auxilliary variables for each population unit.
#' @param n Sample size.
#' @param alpha The objective function parameter which allows the user to control properties of a DAS sample:
#'  alpha = 1 to achieve spatially balanced samples. alpha = 0 < alpha < 1 adds double balance to the sample, setting alpha to 0 forces
#'  the sample to be approximately balanced.
#' @param verbose Boolean if you want to see any output printed to screen. Helpful if taking a
#' long time. Default is FALSE i.e. no informational messages are displayed.
#'
#' @return A list containing one variable, \code{$SampleMatrix} Jn by n matrix containing
#'  population indices, where each row is a candidate sample.
#'

#' @examples
#' # doubleDAS sample ----------------------------------------------
#'
#' my_matrix <- matrix(sample(0:500, 10000, replace = TRUE), nrow = 100, ncol = 100)
#' newpop <- matrix(my_matrix, ncol = 2, byrow = TRUE)
#' n_samples <- 50

#' sampMx <- spbalDAS::doubleDAS(pop = newpop, n = n_samples, verbose = FALSE)
#' # display first twenty sample points.
#' sampMx[1:20]
#'

#' @export
doubleDAS <- function(pop, aux, n, alpha, verbose = FALSE) {
  # parm checks:
  # spbal::validate_parameters("pop", pop) # check - if pop all numeric.
  #       - pop 2 or more columns
  # spbal::validate_parameters("n", n)     #       - n > 0
  # n < nrow(pop)
  if(n >= nrow(pop)){
    base::stop(base::c("spbal(doubleDAS) Number of samples must not exceed number of observations."))
  }

  # alpha needs to be between 0 and 1.
  # Initialisation
  N <- nrow(pop)
  J1 <- floor(N / 2)
  J <- c(J1, J1)
  I <- sample(N)
  SampleMatrix <- matrix(I[1:J1], nrow = J1, ncol = 1)
  print(SampleMatrix[1:20])

  X = aux

  # Compute weight matrix - this takes a fair chunk of time.
  #W <- 1. / as.matrix(arma_dist_al(pop))
  W <- 1. / as.matrix(dist_al(pop))
  W[is.infinite(W)] <- 0
  print(W)

  message("before 2:n loop")

  for (i in 2:n) {
    # Define Ci and Ai
    r <- sample(J[i-1])
    print("r:")
    print(r)
    print("SampleMatrix:")
    print(SampleMatrix)
    print("r[1:J[i]]:")
    print(r[1:J[i]])
    Ci <- SampleMatrix[r[1:J[i]], , drop = FALSE]
    print("Ci:")
    print(Ci)
    OmegaNotCi <- setdiff(1:N, as.vector(Ci))
    print("OmegaNotCi:")
    print(OmegaNotCi)
    r <- sample(length(OmegaNotCi))
    Ai <- OmegaNotCi[r[1:J[i]]]
    print("Ai:")
    print(Ai)

    #message("Calling cppLinearAssignment.")

    CAi <- cbind(Ci, Ai)
    print("dim CAi")
    print(dim(CAi))

    if(!is.numeric(Ci)){
      message(Ci)
      stop("Ci not numeric.")
    }
    if(!is.numeric(Ai)){
      message(Ai)
      stop("Ai not numeric.")
    }
    # Linear Assignment (assign Ai to Ci)
    SampleMatrix <- pfdLinearAssignment(CAi, W, X, alpha, verbose)

    # Update J, the number of points in C_i+1 and A_i+1
    JJ <- max(1, min(J[i], floor(N / (i + 1))))
    J <- c(J, JJ)
  }
  return(SampleMatrix)
}

# Linear Assignment function
pfdLinearAssignment <- function(SampleMatrix, W, X, alpha, verbose = FALSE) {
  nSamples <- nrow(SampleMatrix)
  n <- ncol(SampleMatrix)
  Ci <- SampleMatrix[, 1:(n-1), drop = FALSE]
  Ai <- SampleMatrix[, n]
  # Cost matrix for assignment problem
  SB <- matrix(0, nSamples, nSamples)
  #C <- matrix(0, nSamples, nSamples)
  for (i in 1:nSamples) {
    #message("b4 ci:")
    #message(C[i, ])
    #message("b4 ai:")
    #message(Ai)
    #fred <- Ci[i, ] # save
    SB[i, ] <- SBMoranIvec(as.vector(Ci[i, ]), as.vector(Ai), W)
    #message("cppSBMoranIvec:")
    #message(C[i, ])
    #xxxx <- SBMoranIvec(Ci[i, ], Ai, W)
    #message("SBMoranIvec:")
    #message(xxxx)
  }
  B <- matrix(0, nSamples, nSamples)
  for (i in 1:nSamples) {
    #for (ii in 1:length(X, 2)
    #      b = sum(X(Ci[i]), ii)
    for (ii in 1:ncol(X)) {
      b <- sum(X[Ci[i, ], ii])
      B[i, ] <- B[i, ] + (b + X[Ai, ii])^2
    }
  }
  C = alpha*SB + (1-alpha)*B

  # Solve assignment problem ===============================================
  result = HungarianSolver(C)
  #browser()
  M = result$pairs
  for (i in 1:nSamples)
    #SampleMatrix(i,n) = Ai(M(M(:,1) == i,2));
    SampleMatrix[i, n] <- Ai[ M[M[,1] == i, 2] ]


  # Solve assignment problem
  #M <- solve_LSAP(C, maximum = FALSE)
  #M2 <- spbal:::LAPJV(C)
  #M2 <- XXLAPJV(C)
  #M <- M2$matching
  #browser()

  #for (i in 1:nSamples) {
  #  SampleMatrix[i, n] <- Ai[M[i]]
  #print("M[i]: ")
  #print(M[i])
  #}
  return(SampleMatrix)
}

#SBMoranIvec function
SBMoranIvec <- function(C, A, W) {
  n <- length(C) + 1
  # Compute measure
  #message("W[C, C]:")
  #message(W[C, C])
  WCij <- sum(W[C, C])
  #message("WCij:")
  #message(WCij)

  #print("W[A, C]:")
  #print(W[A, C])

  #print("W[C, A]:")
  #print(W[C, A])

  if (n == 2) {
    MoranVec <- WCij + W[A, C] + W[C, A]
  } else {
    MoranVec <- WCij + rowSums(W[A, C]) + colSums(W[C, A])
  }
  #message("MoranVec:")
  #message(MoranVec)
  return(MoranVec)
}

#Load necessary libraries

#library(pracma) # For pdist2 function
library(clue) # For solving linear assignment problem
# Prepare for plotting
library(ggplot2)
library(reshape2)

#Main function

dist_al <- function(X) {
  n <- nrow(X)
  # Compute squared norms of rows: ||x_i||^2
  G <- rowSums(X^2)  # Sum of squares along rows
  # Compute D^2 = G + G^T - 2 * X * X^T
  D <- -2 * (X %*% t(X))  # -2 * dot products
  D <- sweep(D, 2, G, "+")  # Add G to each column
  D <- sweep(D, 1, G, "+")  # Add G to each row
  # Take square root and ensure non-negative values
  D <- sqrt(abs(D))
  return(D)
}

