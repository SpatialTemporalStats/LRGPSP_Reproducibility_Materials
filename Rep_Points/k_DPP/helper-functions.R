# K-Determinantal Point Processes (k-DPPs)
## Helper functions

#' @title Spectral decomposition of a matrix
#' @description This function is a simple wrapper of base::eigen. It decomposes the input matrix and
#'   outputs a list with the real eigenvalues and eigenvectors.
#' @param L (matrix) Square matrix about which spectral decomposition will be obtained.
#' @inheritDotParams base::eigen
#' @return A list containing the real eigenvalues denoted as `lambda` as well as the eigenvectors
#'   denoted as `LV`.
decompose.matrix <- function(L, ...) {
  decomposed.matrix <- eigen(L, ...)
  real.eigenvalues <- Re(decomposed.matrix[["values"]])
  return(list(lambda = real.eigenvalues,
              LV = decomposed.matrix[["vectors"]]))
}

#' @title Gaussian kernel
#' @description This function obtains a N by N matrix where each element was obtained
#'   based on a Gaussian kernel. It is parametrized in such way that \code{\tau^2} represents the
#'   bandwidth of the kernel.
#' @param X (matrix) N by D numeric matrix, where N is the total number of location coordinates.
#'   Each row must represent the coordinates of a given location.
#' @param tau.sq (numeric) Bandwidth of the Gaussian kernel.
#' @return N by N matrix where each element was obtained based on a Gaussian kernel.
gaussian.kernel <- function(X, tau.sq) {
  assertthat::assert_that(tau.sq > 0)
  dist.X <- fields::rdist(X)
  gauss.kernel <- exp(- as.matrix(dist.X ^ 2) / tau.sq)
  return(gauss.kernel)
}

#' @title Normalization constant for k-DPPs
#' @description This function compute k-th elementary symmetric polynomials
#' @param lambda (numeric) Eigenvalues of the L-ensemble matrix.
#' @param k (integer) Fixed cardinality of the DPP.
#' @return Matrix of elementary symmetric polynomials
compute.elementary.symmetric.polynomial <- function(lambda, k) {
  if ((k %% 1 != 0) || (k <= 0)) {
    stop(simpleError(paste("`k` must be a positive integer.")))
  }
  n.eigen <- length(lambda)
  assertthat::assert_that(n.eigen > 0)
  E <- matrix(0, nrow = k + 1, ncol = n.eigen + 1)
  E[1, ] <- 1
  for (l in seq(2, k + 1)) {
    for (n in seq(2, n.eigen + 1)) {
      E[l, n] <- E[l, n - 1] + lambda[n - 1] * E[l - 1, n - 1]
    }
  }
  return(E)
}

#' @title Normalization constant for k-DPPs
#' @description This function compute the logrithm for the k-th elementary symmetric polynomials
#' @param lambda (numeric) Eigenvalues of the L-ensemble matrix.
#' @param k (integer) Fixed cardinality of the DPP.
#' @return Matrix of log(elementary symmetric polynomials)
log.elementary.symmetric.polynomial <- function(lambda, k) {
  if ((k %% 1 != 0) || (k <= 0)) {
    stop(simpleError(paste("`k` must be a positive integer.")))
  }
  n.eigen <- length(lambda)
  assertthat::assert_that(n.eigen > 0)
  E <- matrix(log(0), nrow = k + 1, ncol = n.eigen + 1)
  E[1, ] <- log(1)
  for (l in seq(2, k + 1)) {
    for (n in seq(2, n.eigen + 1)) {
      if(E[l-1,n-1]==-Inf){
        E[l,n]=E[l,n-1]
      }
      else{
        E[l, n] <- E[l-1,n-1]+log(lambda[n-1]+exp(E[l,n-1]-E[l-1,n-1]))
      }
    }
  }
  return(E)
}


#' @title Sample k eigenvectors
#' @description This function is responsible for sampling eigenvectors, where the probability depends
#'   on its associated eigenvalue.
#' @inheritParams compute.elementary.symmetric.polynomial
#' @return k-dimensional array indicating the indices of the sampled eigenvector
sample.eigenvector <- function(lambda, k) {
  # Initialize iterator as the total number of eigenvalues
  i <- length(lambda)
  if (k > i) {
    stop(simpleError(paste("Input 'k' must be less than or equal to", i, "i.e. the total number of",
                           "eigenvalues as indicated in the 'lambda' input.")))
  }
  # Compute elementary symmetric polynomial
  E <- compute.elementary.symmetric.polynomial(lambda, k)
  # Initialize object to store the k selected eigenvalues
  selected.eigenvalues <- array(0, dim = k)
  # Iterator indicating how many eigenvectors are yet to be sampled
  remaining <- k
  # This loop will stop once `k` eigenvectors have been selected
  while (remaining > 0) {
    if (i == remaining) {
      marginal.p <- 1
    } else {
      marginal.p <- lambda[i] * E[remaining, i] / E[remaining + 1, i + 1]
    }
    if (is.na(marginal.p)) {
       stop(simpleError("NaN detected. Marginal is undefined."))
    }
    if (runif(1) < marginal.p) {
      selected.eigenvalues[remaining] <- i
      remaining <- remaining - 1
    }
    i <- i - 1
  }
  return(selected.eigenvalues)
}


#' @title Sample k eigenvectors
#' @description This function is responsible for sampling eigenvectors, where the probability depends
#'   on its associated eigenvalue.
#' @inheritParams log.elementary.symmetric.polynomial
#' @return k-dimensional array indicating the indices of the sampled eigenvector
log.sample.eigenvector <- function(lambda, k) {
  # Initialize iterator as the total number of eigenvalues
  i <- length(lambda)
  if (k > i) {
    stop(simpleError(paste("Input 'k' must be less than or equal to", i, "i.e. the total number of",
                           "eigenvalues as indicated in the 'lambda' input.")))
  }
  # Compute log of elementary symmetric polynomial
  E <- log.elementary.symmetric.polynomial(lambda, k)
  # Initialize object to store the k selected eigenvalues
  selected.eigenvalues <- array(0, dim = k)
  # Iterator indicating how many eigenvectors are yet to be sampled
  remaining <- k
  # This loop will stop once `k` eigenvectors have been selected
  while (remaining > 0) {
    if (i == remaining) {
      marginal.p <- 1
    } else {
      marginal.p <- lambda[i]*exp(E[remaining, i]-E[remaining + 1, i + 1])
    }
    if (is.na(marginal.p)) {
      stop(simpleError("NaN detected. Marginal is undefined."))
    }
    if (runif(1) < marginal.p) {
      selected.eigenvalues[remaining] <- i
      remaining <- remaining - 1
    }
    i <- i - 1
  }
  return(selected.eigenvalues)
}
  
