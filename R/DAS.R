
#' @name DAS
#'
#' @title Dynamic Assignment Sampling (DAS).
#'
#' @description DAS draws spatially balanced samples from areal resources. Needs updating!!!
#'
#' @author This function was first written by Blair Robertson in Matlab and later re-written
#' in R/C++ by Phil Davies.
#'
#' @param pop Population point pairs. What data types will be supported?
#' @param n The number of sample points to draw from the population.
#' @param J1 maximum number of candidate samples 2<= J1 <= [N/2].
#' @param n_threads The number of threads to use, when calculation the cost matrix, a value of 0
#' will result in a thread being dispatched on each available logical CPU. The default is n_threads=1.
#' This parameter takes affect only when the current J value is greater than or equal to 256 and
#' n_threads is greater than 1.
#' @param verbose Boolean if you want to see any output printed to screen. Helpful if taking a
#' long time. Default is FALSE i.e. no informational messages are displayed.
#'
#' @return A list containing one variable, \code{$SampleMatrix} Jn by n matrix containing
#' population indices, where each row is a candidate sample. I could also return the cost
#' and pairs.
#'

#' @examples
#' # DAS sample ----------------------------------------------
#'
#' my_matrix <- matrix(sample(0:500, 10000, replace = TRUE), nrow = 100, ncol = 100)
#' newpop <- matrix(my_matrix, ncol = 2, byrow = TRUE)
#' n_samples <- 50

#' sampMx <- spbalDAS::DAS(pop = newpop, n = n_samples, verbose = FALSE)
#' # display first twenty sample points.
#' sampMx[1:20]
#'

#' @export
DAS <- function(pop,
                n,
                J1,
                n_threads = 1,
                verbose = FALSE){

  options(digits.secs = 3)

  if (verbose){
    current_time <- Sys.time()
    timestamp_message <- paste(format(current_time, "[%Y-%m-%d %H:%M:%OS3]"), " Starting...")
    message(timestamp_message)
  }

  # Validate parameters.
  #xvalidate_parameters("pop", pop) # ensure pop values are all numeric.
  #       - pop 2 or more columns
  #xvalidate_parameters("n", n)     # ensure n larger than 0.

  if(n >= nrow(pop)){
    base::stop(base::c("spbal(DAS) Number of samples to take must not exceed number of observations."))
  }

  # Initialisation
  N <- nrow(pop)
  # maximum number of candidate samples 2<= J1 <= [N/2]
  #J1 <- floor(N / 2)
  J <- c(J1, J1)
  I <- sample(N)
  SampleMatrix <- matrix(I[1:J1], nrow = J1, ncol = 1)
  if(verbose){message(SampleMatrix[1:20])}


  # Compute weight matrix - this can take a fair chunk of time.
  W <- 1. / as.matrix(arma_dist_al(pop))

  if (verbose){
    current_time <- Sys.time()
    timestamp_message <- paste(format(current_time, "[%Y-%m-%d %H:%M:%OS3]"), " After arma_dist_al")
    message(timestamp_message)
  }

  W[is.infinite(W)] <- 0
  #if (verbose){message(W)}

  # this:
  SampleMatrix <- LinearAssignmentProblem_cpp(W, n, N, J[1], n_threads = 6, verbose = FALSE)
  # replaces this:
  #for (i in 2:n) {
  #  message('2:n')
    # Define Ci and Ai
  #  r <- sample(J[i-1])
  #  Ci <- SampleMatrix[r[1:J[i]], , drop = FALSE]
  #  OmegaNotCi <- setdiff(1:N, as.vector(Ci))
    # Draw an SRS of Ji units from U, denoted Ai={a1,...,aji}
  #  r <- sample(length(OmegaNotCi))
  #  Ai <- OmegaNotCi[r[1:J[i]]]

  #  CAi <- cbind(Ci, Ai)

   # if(!is.numeric(Ci)){
    #  message(Ci)
    #  base::stop(base::c("spbal(DAS) Ci not numeric."))
    #}
    #if(!is.numeric(Ai)){
    #  message(Ai)
    #  base::stop(base::c("spbal(DAS) Ai not numeric."))
    #}
    # Linear Assignment (assign Ai to Ci)
    # Solve assignment problem to find the optimal permutation sigma of {1,...,Ji}
    #result <- SolveLinearAssignment(CAi, W, verbose)
    #SampleMatrix <- result #$SampleMatrix

    # Update J, the number of points in C_i+1 and A_i+1
    #JJ <- max(1, min(J[i], floor(N / (i + 1))))
    #J <- c(J, JJ)
  #}

  if (verbose){
    current_time <- Sys.time()
    timestamp_message <- paste(format(current_time, "[%Y-%m-%d %H:%M:%OS3]"), " Exiting DAS()...")
    message(timestamp_message)
  }

  # assign the spbal attribute to the sample being returned, i.e. the function that created it.
  base::attr(sample, "spbal") <- "DAS"
  result <- base::list(SampleMatrix = SampleMatrix)
  return(SampleMatrix)
}

#' @name SolveLinearAssignment
#'
#' @title SolveLinearAssignment.
#'
#' @description SolveLinearAssignment
#'
#' @author This function was first written by Blair Robertson in Matlab and later re-written in R by Phil Davies.
#'
#' @param SampleMatrix Population points.
#' @param W Weight matrix.
#' @param verbose Boolean if you want to see any output printed to screen. Helpful if taking a
#' long time. Default is FALSE i.e. no informational messages are displayed.
#'
#' @return A list containing one variable, \code{$SampleMatrix} Jn by n matrix containing
#'  population indices, where each row is a candidate sample.

SolveLinearAssignment <- function(SampleMatrix, W, verbose = FALSE) {

  current_time <- Sys.time()
  timestamp_message <- paste(format(current_time, "[%Y-%m-%d %H:%M:%OS3]"), " Enter SolveLinearAssignment")
  print(timestamp_message)

  nSamples <- nrow(SampleMatrix)
  n <- ncol(SampleMatrix)
  Ci <- SampleMatrix[, 1:(n-1), drop = FALSE]
  Ai <- SampleMatrix[, n]
  # Cost matrix for assignment problem
  C <- matrix(0, nSamples, nSamples)
  for (i in 1:nSamples) {
    C_i <- as.vector(SampleMatrix[i, 1:(n-1), drop = FALSE])     # cluster for sample i
    A_all <- as.vector(SampleMatrix[, n])                 # all outside nodes (column 10 across all samples)

    #C[i, ] <- rSBMoranIvecBatch(C_i, A_all, W)
    C[i, ] <- SBMoranIvec(C_i, A_all, W)
  }

  # Solve assignment problem ===============================================
  # The RcppHungarian package solves one problem well. Namely, it solves the minimum cost
  # bipartite matching problem.
  current_time <- Sys.time()
  timestamp_message <- paste(format(current_time, "[%Y-%m-%d %H:%M:%OS3]"), " calling h_arma")
  print(timestamp_message)

  #result = hungarian_arma(C)
  result = my_hungarian(C)

  current_time <- Sys.time()
  timestamp_message <- paste(format(current_time, "[%Y-%m-%d %H:%M:%OS3]"), " back from h_arma")
  print(timestamp_message)
  # returns $cost and $pairs
  M = result$pairs

  for (i in 1:nSamples) {
    #SampleMatrix[i, n] <- Ai[M[i]]
    SampleMatrix[i, n] <- Ai[ M[M[,1] == i, 2] ]
  }
  #result <- base::list(SampleMatrix = SampleMatrix,
  #                     cost = cost,
  #                     pairs = pairs)
  return(SampleMatrix)
}


#' Sample spread using the equal probability modified Moran's I measure
#'
#' @param C Current sample indices (integer vector)
#' @param A Candidate indices to evaluate (integer vector).
#'          The function returns one value per element of `A`.
#' @param W Spatial weight matrix (square numeric matrix)
#'
#' @return Numeric vector of length `length(A)` containing the modified
#'         Moran's I spread measure for each candidate in `A`.
SBMoranIvec <- function(C, A, W) {
  # Initialize
  n <- length(C) + 1L

  # Sum of all pairwise weights inside the current sample
  WCij <- sum(W[C, C])

  # Contribution from each candidate in A to the current sample C,
  # plus contribution from C back to each candidate in A.
  # (These are identical when W is symmetric, but we keep both to
  #  match the original MATLAB logic exactly.)
  contrib_AC <- rowSums(W[A, C, drop = FALSE])
  contrib_CA <- colSums(W[C, A, drop = FALSE])

  # Combine
  MoranVec <- WCij + contrib_AC + contrib_CA
  return(MoranVec)
}




#' @name rSBMoranIvec
#'
#' @title rSBMoranIvec.
#'
#' @description Moran’s I is a measure of spatial autocorrelation–how related the values of a variable are based
#' on the locations where they were measured.
#' To calculate Moran’s I, we will need to generate a matrix of inverse distance weights.
#' In the matrix, entries for pairs of points that are close together are higher than for pairs of points that are far apart.
#' For simplicity, we will treat the latitude and longitude as values on a plane rather than on a sphere–our locations are
#' close together and far from the poles. When using latitude and longitude coordinates from more distant locations,
#' it’s wise to calculate distances based on spherical coordinates (the geosphere package can be used).
#'
#' @author This function was first written by Blair Robertson in Matlab and later re-written in R by Phil Davies.
#'
#' @param C bla.
#' @param A bla.
#' @param W bla.
#'
#' @return To be filled in.

rSBMoranIvec <- function(C, A, W) {
  browser()
  n <- length(C) + 1
  WCij <- sum(W[C, C])

  if (n == 2) {
    MoranVec <- WCij + t(W[A, C]) + W[C, A]
  } else {
    MoranVec <- WCij + rowSums(W[A, C]) + colSums(W[C, A])
  }
  return(MoranVec)
}

rSBMoranIvecNew <- function(C, A, W) {
  internal <- sum(W[C, C])
  cut_out  <- rowSums(as.matrix(W[A, C]))   # A → C
  cut_in   <- colSums(as.matrix(W[C, A]))   # C → A
  return(internal + cut_out + cut_in)
}

rSBMoranIvecBatch <- function(C, A_vec, W) {
  # C: vector of nodes in the cluster (fixed)
  # A_vec: vector of many outside nodes (one per sample, length nSamples)
  # W: weight matrix

  #current_time <- Sys.time()
  #timestamp_message <- paste(format(current_time, "[%Y-%m-%d %H:%M:%OS3]"), " Enter rSBMoranIvecBatch")
  #print(timestamp_message)

  internal <- sum(W[C, C])  # scalar, same for all

  # W[A_vec, C]: rows = each a in A_vec, columns = nodes in C
  out_sum <- rowSums(as.matrix(W[A_vec, C]))   # A_vec → C

  # W[C, A_vec]: rows = C, columns = each a
  in_sum <- colSums(as.matrix(W[C, A_vec]))    # C → A_vec

  result <- internal + out_sum + in_sum

  #current_time <- Sys.time()
  #timestamp_message <- paste(format(current_time, "[%Y-%m-%d %H:%M:%OS3]"), " Exit rSBMoranIvecBatch")
  #print(timestamp_message)
  return(result)
}
