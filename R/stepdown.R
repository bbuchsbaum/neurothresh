#' Westfall-Young Step-Down Procedure
#'
#' Implements the Westfall-Young step-down maxT procedure for
#' strong control of the family-wise error rate (FWER).
#'
#' @param observed Numeric vector of observed test statistics
#' @param null_matrix Matrix of null statistics (n_perm x n_tests)
#' @param alpha Significance level (default 0.05)
#'
#' @return A data.frame with columns:
#'   \describe{
#'     \item{region}{Index of the region (1 to n_tests)}
#'     \item{score}{Observed test statistic}
#'     \item{p_adj}{Adjusted p-value (monotonically non-decreasing)}
#'     \item{rejected}{Logical, whether the null is rejected}
#'   }
#'
#' @details
#' The Westfall-Young step-down procedure provides strong FWER control
#' with higher power than single-step methods by:
#' 1. Sorting regions by observed score (descending)
#' 2. Computing successive maxima from the null distribution
#' 3. Ensuring adjusted p-values are monotonically non-decreasing
#'
#' This is the recommended method for hierarchical LR-MFT as it
#' preserves the logical ordering of rejections.
#'
#' @references
#' Westfall, P. H., & Young, S. S. (1993). Resampling-based multiple
#' testing: Examples and methods for p-value adjustment. Wiley.
#'
#' @examples
#' \dontrun{
#' # Observed scores for 5 regions
#' observed <- c(3.5, 2.1, 4.2, 1.8, 2.9)
#'
#' # Null distribution from 1000 permutations
#' null_matrix <- matrix(rnorm(5000), nrow = 1000, ncol = 5)
#'
#' # Run step-down procedure
#' result <- wy_stepdown(observed, null_matrix, alpha = 0.05)
#' }
#'
#' @export
wy_stepdown <- function(observed, null_matrix, alpha = 0.05) {
  m <- length(observed)
  B <- nrow(null_matrix)

  if (ncol(null_matrix) != m) {
    stop("Number of columns in null_matrix must match length of observed")
  }

  if (m == 0) {
    return(data.frame(
      region = integer(0),
      score = numeric(0),
      p_adj = numeric(0),
      rejected = logical(0)
    ))
  }

  # Sort regions by observed score (descending)
  ord <- order(observed, decreasing = TRUE)
  obs_sorted <- observed[ord]

  # Reorder null matrix columns to match
  null_sorted <- null_matrix[, ord, drop = FALSE]

  # Compute successive maxima from the right
  # max_null[b, j] = max(null_sorted[b, j:m])
  successive_max <- matrix(NA_real_, nrow = B, ncol = m)
  successive_max[, m] <- null_sorted[, m]

  for (j in (m - 1):1) {
    successive_max[, j] <- pmax(null_sorted[, j], successive_max[, j + 1])
  }

  # Compute raw p-values
  p_raw <- numeric(m)
  for (j in seq_len(m)) {
    p_raw[j] <- (1 + sum(successive_max[, j] >= obs_sorted[j])) / (B + 1)
  }

  # Enforce monotonicity (p_adj[j] >= p_adj[j-1])
  p_adj <- p_raw
  for (j in 2:m) {
    p_adj[j] <- max(p_adj[j], p_adj[j - 1])
  }

  # Create output in original order
  result <- data.frame(
    region = seq_len(m),
    score = observed,
    p_adj = numeric(m),
    rejected = logical(m)
  )

  # Map back to original order
  result$p_adj[ord] <- p_adj
  result$rejected <- result$p_adj <= alpha

  result
}

#' Single-step MaxT Procedure
#'
#' Implements the single-step maxT procedure for FWER control.
#' Less powerful than step-down but simpler.
#'
#' @param observed Numeric vector of observed test statistics
#' @param null_matrix Matrix of null statistics (n_perm x n_tests)
#' @param alpha Significance level (default 0.05)
#'
#' @return A data.frame with columns:
#'   \describe{
#'     \item{region}{Index of the region}
#'     \item{score}{Observed test statistic}
#'     \item{p_adj}{Adjusted p-value}
#'     \item{rejected}{Logical, whether the null is rejected}
#'   }
#'
#' @export
maxT_singlestep <- function(observed, null_matrix, alpha = 0.05) {
  m <- length(observed)
  B <- nrow(null_matrix)

  if (m == 0) {
    return(data.frame(
      region = integer(0),
      score = numeric(0),
      p_adj = numeric(0),
      rejected = logical(0)
    ))
  }

  # Compute max statistic across all tests for each permutation
  max_null <- apply(null_matrix, 1, max)

  # P-value for each test
  p_adj <- vapply(observed, function(t) {
    (1 + sum(max_null >= t)) / (B + 1)
  }, numeric(1))

  data.frame(
    region = seq_len(m),
    score = observed,
    p_adj = p_adj,
    rejected = p_adj <= alpha
  )
}

#' Compute adjusted p-values from max null distribution
#'
#' Helper function to compute FWER-corrected p-values given
#' a max null distribution.
#'
#' @param observed Numeric vector of observed values
#' @param max_null Numeric vector of max statistics under null
#'
#' @return Numeric vector of adjusted p-values
#'
#' @keywords internal
fwer_p_from_maxnull <- function(observed, max_null) {
  B <- length(max_null)
  vapply(observed, function(t) {
    (1 + sum(max_null >= t)) / (B + 1)
  }, numeric(1))
}
