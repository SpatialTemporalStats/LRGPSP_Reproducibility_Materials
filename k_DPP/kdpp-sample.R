# K-Determinantal Point Processes (k-DPPs)
## Sampling

library(magrittr)
library(foreach)

#' @title k-DPP sample
#' @description Obtains a single realization of a k-DPP.
#' @param L (list) Named list. The names of the corresponding elements are `eigenval.name` and
#'   `eigenvec.name`, indicating the eigenvalues and eigenvectors respectively.
#' @param k (numeric) Fixed cardinality of the DPP.
#' @param eigenval.name (character) Name of element of L where eigenvalues are stored.
#'   Defaults to "lambda".
#' @param eigenvec.name (character) Name of element of L where eigenvectors are stored.
#'   Defaults to "LV".
#' @return Sorted vector of indices resulting from one realization of a k-DPP.
k.dpp.sample <- function(L, k, eigenval.name = "lambda", eigenvec.name = "LV") {
  if ((k %% 1 != 0) || (k <= 0)) {
    stop(simpleError(paste("`k` must be a positive integer.")))
  }
  if (!all(c(eigenval.name, eigenvec.name) %in% names(L))) {
    stop(simpleError(paste("Input `L` must represent decomposed matrix by a named list with",
                           eigenval.name, "indicating the eigenvalues and", eigenvec.name,
                           "indicating the eigenvectors.")))
  }

  lambda <- L[[eigenval.name]]
  LV <- L[[eigenvec.name]]

  if (k > length(lambda)) {
    stop(simpleError(paste("Input 'k' must be less than or equal to", length(lambda), "i.e. the",
                           "total number of eigenvalues of 'L'.")))
  }

  # Initializing vector to save k-DPP sample
  out <- array(0, dim = k)

  # Selecting eigenvalues for k-DPP
  # v <-sample.eigenvector(lambda, k)
  v <- log.sample.eigenvector(lambda, k)

  # Corresponding columns in LV matrix
  V <- as.matrix(LV[, v])

  for (i in seq(k, 1)) {
    # Compute probabilities for each item
    sq.V <- rowSums(V ^ 2)
    P <- sq.V / sum(sq.V)

    # Choose a new item to include
    out[i] <- which.max(runif(1) <= cumsum(P))

    # If the sample should be of size 1, then stop loop
    if (k == 1) break

    # Choose a vector to eliminate
    j <- which.max(V[out[i], ] != 0) # first nonzero element

    # Update V
    Vj <- as.matrix(V[, j])
    V <- as.matrix(V[, -j])
    V <- V - Vj %*% (V[out[i], ] / Vj[out[i]])

    # Orthonormalize
    if (i > 1) {
      for (a in seq(1, i - 1)) {
        if (a > 1) {
          for (b in seq(1, a - 1)) {
            V[, a] <- V[, a] - as.numeric(t(V[, a]) %*% V[, b]) * V[, b]
          }
        }
        V[, a] <- V[, a] / max(svd(V[, a])[["d"]])
      }
    }
  }
  return(sort(out))
}

#' @title Empirically maximize k-DPP
#' @description This function is a wrapper of `k.dpp.sample`. Essentially, `n.sim` realizations of
#'   a k-DPP is obtained and the sample yielding the maximum sub-determinant is outputted.
#' @inheritParams k.dpp.sample
#' @param n.sim (integer) Number of realizations of the k-DPP.
#' @param parallelize (logical) Indicates whether the k-DPP realizations should be run in parallel.
#'   Defaults to FALSE (i.e. run sequentially). When it is set to TRUE, it assumes that the parallel
#'   backend was previously set up (e.g. `doMC::registerDoMC()`), otherwise the function will be run
#'   sequentially.
#' @param return.samples (logical) Indicates whether the individual samples should be outputted.
#'   Defaults to FALSE.
#' @return k-DPP realization with the maximum sub-determinant.
empirically.maximize.k.dpp <- function(L, k, n.sim,
                                       parallelize = FALSE,
                                       return.samples = FALSE,
                                       eigenval.name = "lambda",
                                       eigenvec.name = "LV") {
  if (!all(c(eigenval.name, eigenvec.name) %in% names(L))) {
    stop(simpleError(paste("Input `L` must represent decomposed matrix by a named list with",
                           eigenval.name, "indicating the eigenvalues and", eigenvec.name,
                           "indicating the eigenvectors.")))
  }
  if ((k %% 1 != 0) || (k <= 0)) {
    stop(simpleError(paste("`k` must be a positive integer.")))
  }
  `%DO%` <- if (parallelize) `%dopar%` else `%do%`
  LV <- L[[eigenvec.name]]
  dpp.sample.output.list <- foreach::foreach(sim = seq(n.sim)) %DO% {
    determinant <- NULL
    one.sample <- tryCatch({
      k.dpp.sample(L, k, eigenval.name, eigenvec.name)
    }, error = function(e) NULL)
    determinant <- if (!is.null(one.sample)) det(LV[one.sample, one.sample])
    return(list(sample = one.sample, determinant = determinant))
  }
  determinants <- dpp.sample.output.list %>%
    purrr::map_dbl("determinant")
  if (length(determinants) != n.sim) {
    warning("One or more samples were not produced due to numerical instability.")
  }
  max.det.idx <- which.max(determinants)
  if (length(max.det.idx) == 0) {
    stop(simpleError("No sample produced due to numerical instability."))
  }
  # Return the indices of with the maximum sub-determinant
  output.list <- list(max.dpp.sample = dpp.sample.output.list[[max.det.idx]][["sample"]])
  if (return.samples) {
    # If requested, append info about the samples
    sample.matrix <- dpp.sample.output.list %>%
      purrr::map("sample") %>%
      rlist::list.rbind()
    output.list <- c(output.list, list(samples = sample.matrix, determinants = determinants))
  }
  return(output.list)
}