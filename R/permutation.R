#' Generate null distribution via sign-flipping
#'
#' Generates null scores for regions by randomly sign-flipping the
#' Z-score map. This is valid when the original data are symmetric
#' around zero under the null (e.g., group difference maps).
#'
#' @param z_vec Numeric vector of Z-scores (mask-space)
#' @param pi_vec Numeric vector of prior weights (mask-space)
#' @param regions List of index vectors defining regions
#' @param n_perm Number of permutations
#' @param kappa_grid Numeric vector of positive kappa values
#' @param seed Random seed for reproducibility (optional)
#'
#' @return Matrix of dimension (n_perm x n_regions) with null scores
#'
#' @details
#' For each permutation, a random sign (+1 or -1) is drawn independently
#' for each voxel and multiplied with the Z-scores. The omnibus score
#' (max of U_0 and best S_kappa) is computed for each region.
#'
#' @examples
#' \dontrun{
#' z_vec <- rnorm(1000)
#' pi_vec <- rep(1, 1000)
#' regions <- list(1:100, 101:200, 201:300)
#' null_matrix <- generate_null_scores(z_vec, pi_vec, regions,
#'                                     n_perm = 1000)
#' }
#'
#' @export
generate_null_scores <- function(z_vec, pi_vec, regions, n_perm,
                                 kappa_grid = c(0.5, 1.0, 2.0),
                                 seed = NULL) {
  if (!is.null(seed)) {
    set.seed(seed)
  }

  n_regions <- length(regions)
  n_voxels <- length(z_vec)

  null_matrix <- matrix(NA_real_, nrow = n_perm, ncol = n_regions)

  for (b in seq_len(n_perm)) {
    # Random sign-flip
    signs <- sample(c(-1, 1), n_voxels, replace = TRUE)
    z_perm <- signs * z_vec

    # Score each region
    for (r in seq_along(regions)) {
      idx <- regions[[r]]
      result <- score_set_omnibus(idx, z_perm, pi_vec, kappa_grid)
      null_matrix[b, r] <- result$T_omnibus
    }
  }

  null_matrix
}

#' Generate null scores for parcel-level analysis
#'
#' Vectorized null score generation for parcels using efficient
#' group-wise operations.
#'
#' @param z_vec Numeric vector of Z-scores (mask-space)
#' @param pi_vec Numeric vector of prior weights (mask-space)
#' @param lab_vec Integer vector of parcel labels (mask-space)
#' @param n_perm Number of permutations
#' @param kappa_grid Numeric vector of positive kappa values
#' @param seed Random seed for reproducibility (optional)
#'
#' @return Matrix of dimension (n_perm x n_parcels) with null scores
#'
#' @details
#' This is more efficient than \code{\link{generate_null_scores}} for
#' parcel-level analysis because it uses vectorized group operations
#' (rowsum) instead of looping over regions.
#'
#' @export
generate_null_scores_parcels <- function(z_vec, pi_vec, lab_vec, n_perm,
                                         kappa_grid = c(0.5, 1.0, 2.0),
                                         seed = NULL) {
  if (!is.null(seed)) {
    set.seed(seed)
  }

  n_voxels <- length(z_vec)

  # Precompute parcel constants
  pi2_vec <- pi_vec * pi_vec
  den1 <- as.numeric(rowsum(pi_vec, group = lab_vec, reorder = FALSE))
  den2 <- as.numeric(rowsum(pi2_vec, group = lab_vec, reorder = FALSE))
  log_den1 <- log(den1)
  inv_sqrt_den2 <- 1 / sqrt(den2)

  parcels <- sort(unique(lab_vec))
  n_parcels <- length(parcels)

  null_matrix <- matrix(NA_real_, nrow = n_perm, ncol = n_parcels)

  for (b in seq_len(n_perm)) {
    # Random sign-flip
    signs <- sample(c(-1, 1), n_voxels, replace = TRUE)
    z_perm <- signs * z_vec

    # Vectorized parcel scoring
    scores <- .score_parcels_vec(z_perm, pi_vec, lab_vec,
                                 log_den1, inv_sqrt_den2, kappa_grid)
    null_matrix[b, ] <- scores
  }

  null_matrix
}

#' Vectorized parcel scoring
#'
#' @keywords internal
.score_parcels_vec <- function(z_vec, pi_vec, lab_vec,
                               log_den1, inv_sqrt_den2, kappa_grid) {
  # U_0: variance-stabilized diffuse score
  num0 <- as.numeric(rowsum(pi_vec * z_vec, group = lab_vec, reorder = FALSE))
  U0 <- num0 * inv_sqrt_den2

  # S_kappa: soft-max scores
  best_soft <- rep(-Inf, length(U0))

  for (kappa in kappa_grid) {
    if (kappa <= 0) next

    a <- kappa * z_vec
    a_max <- max(a)
    numk <- as.numeric(rowsum(pi_vec * exp(a - a_max),
                              group = lab_vec, reorder = FALSE))
    Sk <- (a_max + log(numk) - log_den1) / kappa
    best_soft <- pmax(best_soft, Sk)
  }

  pmax(U0, best_soft)
}

#' Create permutation function for hierarchical analysis
#'
#' Creates a function that returns sign-flipped Z vectors for use
#' in hierarchical descent.
#'
#' @param z_vec Numeric vector of Z-scores (mask-space)
#' @param n_perm Number of permutations
#' @param seed Random seed for reproducibility (optional)
#'
#' @return A function that takes permutation index b (1 to n_perm)
#'   and returns a sign-flipped Z vector
#'
#' @details
#' The sign matrix is pre-generated and stored, so calling the

#' returned function multiple times with the same index gives
#' the same result.
#'
#' @export
make_perm_fun <- function(z_vec, n_perm, seed = NULL) {
  if (!is.null(seed)) {
    set.seed(seed)
  }

  n <- length(z_vec)

  # Pre-generate all sign matrices
  sign_matrix <- matrix(sample(c(-1, 1), n * n_perm, replace = TRUE),
                        nrow = n_perm, ncol = n)

  function(b) {
    if (b < 1 || b > n_perm) {
      stop("Permutation index out of range")
    }
    sign_matrix[b, ] * z_vec
  }
}
