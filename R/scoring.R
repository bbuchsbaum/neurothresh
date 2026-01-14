#' Compute prior-weighted soft evidence statistic T_kappa(R)
#'
#' Computes the prior-weighted log-sum-exp statistic for a region R:
#' \deqn{T_\kappa(R) = \log \sum_{v \in R} \pi(v) \cdot \exp(\kappa \cdot Z(v))}
#'
#' @param indices Integer vector of 1-based mask-space indices defining the region
#' @param z_vec Numeric vector of Z-scores in mask-space
#' @param pi_vec Numeric vector of prior weights in mask-space
#' @param kappa Numeric scalar, temperature parameter (default 1.0).
#'   Higher kappa emphasizes peaks; kappa -> 0 gives diffuse detection.
#'
#' @return Numeric scalar, the T_kappa score for the region
#'
#' @details
#' This is the core statistic for LR-MFT. For kappa = 0, use
#' \code{\link{score_set_stabilized}} instead for variance-stabilized
#' diffuse detection.
#'
#' @seealso \code{\link{score_set_stabilized}} for kappa = 0 case
#'
#' @examples
#' \dontrun{
#' # Score a region with 100 voxels
#' z_vec <- rnorm(1000)
#' pi_vec <- rep(1, 1000)
#' indices <- 1:100
#' score <- score_set(indices, z_vec, pi_vec, kappa = 1.0)
#' }
#'
#' @export
score_set <- function(indices, z_vec, pi_vec, kappa = 1.0) {
  if (length(indices) == 0) {
    return(-Inf)
  }

  z <- z_vec[indices]
  w <- pi_vec[indices]

  # Handle NA values
  valid <- !is.na(z) & !is.na(w)
  if (!any(valid)) {
    return(-Inf)
  }

  z <- z[valid]
  w <- w[valid]

  if (kappa <= 0) {
    # For kappa = 0, return the weighted sum (not normalized)
    return(sum(w * z))
  }

  # Log-sum-exp with numerical stability
  a <- kappa * z
  a_max <- max(a)
  log_sum <- a_max + log(sum(w * exp(a - a_max)))

  return(log_sum)
}

#' Compute variance-stabilized score U_0(R)
#'
#' Computes the variance-stabilized diffuse detection score:
#' \deqn{U_0(R) = \frac{\sum_{v \in R} \pi(v) \cdot Z(v)}{\sqrt{\sum_{v \in R} \pi(v)^2}}}
#'
#' This provides fair comparison across regions of different sizes by
#' normalizing by the effective sample size.
#'
#' @param indices Integer vector of 1-based mask-space indices defining the region
#' @param z_vec Numeric vector of Z-scores in mask-space
#' @param pi_vec Numeric vector of prior weights in mask-space
#'
#' @return A list with components:
#'   \describe{
#'     \item{U0}{The variance-stabilized score}
#'     \item{n_eff}{The effective sample size: (sum(pi))^2 / sum(pi^2)}
#'   }
#'
#' @details
#' Under the null hypothesis, U_0(R) has approximately unit variance
#' regardless of region size, enabling fair comparison across scales.
#'
#' The effective sample size \code{n_eff} represents the equivalent
#' number of independent observations.
#'
#' @seealso \code{\link{score_set}} for kappa > 0 focal detection
#'
#' @examples
#' \dontrun{
#' # Compare scores across regions of different sizes
#' z_vec <- rnorm(1000)
#' pi_vec <- rep(1, 1000)
#'
#' # Small region
#' small <- score_set_stabilized(1:10, z_vec, pi_vec)
#'
#' # Large region
#' large <- score_set_stabilized(1:100, z_vec, pi_vec)
#'
#' # Both have comparable variance under null
#' }
#'
#' @export
score_set_stabilized <- function(indices, z_vec, pi_vec) {
  if (length(indices) == 0) {
    return(list(U0 = -Inf, n_eff = 0))
  }

  z <- z_vec[indices]
  w <- pi_vec[indices]

  # Handle NA values
  valid <- !is.na(z) & !is.na(w)
  if (!any(valid)) {
    return(list(U0 = -Inf, n_eff = 0))
  }

  z <- z[valid]
  w <- w[valid]

  sum_wz <- sum(w * z)
  sum_w <- sum(w)
  sum_w2 <- sum(w * w)

  if (sum_w2 <= 0) {
    return(list(U0 = 0, n_eff = 0))
  }

  U0 <- sum_wz / sqrt(sum_w2)
  n_eff <- (sum_w * sum_w) / sum_w2

  list(U0 = U0, n_eff = n_eff)
}

#' Compute omnibus score combining diffuse and focal detection
#'
#' Computes the maximum of the variance-stabilized diffuse score U_0
#' and the best soft-max score S_kappa across a grid of kappa values.
#'
#' @param indices Integer vector of 1-based mask-space indices
#' @param z_vec Numeric vector of Z-scores in mask-space
#' @param pi_vec Numeric vector of prior weights in mask-space
#' @param kappa_grid Numeric vector of positive kappa values to try
#'
#' @return A list with components:
#'   \describe{
#'     \item{T_omnibus}{The maximum score across all statistics}
#'     \item{U0}{The variance-stabilized diffuse score}
#'     \item{S_best}{The best soft-max score}
#'     \item{kappa_best}{The kappa value that gave S_best}
#'   }
#'
#' @details
#' The omnibus statistic is:
#' \deqn{T(R) = \max(U_0(R), \max_{\kappa} S_\kappa(R))}
#'
#' where S_kappa is the normalized soft-max:
#' \deqn{S_\kappa(R) = \frac{1}{\kappa} \left[ \log \sum \pi \exp(\kappa Z) - \log \sum \pi \right]}
#'
#' @export
score_set_omnibus <- function(indices, z_vec, pi_vec,
                              kappa_grid = c(0.5, 1.0, 2.0)) {
  if (length(indices) == 0) {
    return(list(T_omnibus = -Inf, U0 = -Inf, S_best = -Inf, kappa_best = NA))
  }

  z <- z_vec[indices]
  w <- pi_vec[indices]

  # Handle NA values
  valid <- !is.na(z) & !is.na(w)
  if (!any(valid)) {
    return(list(T_omnibus = -Inf, U0 = -Inf, S_best = -Inf, kappa_best = NA))
  }

  z <- z[valid]
  w <- w[valid]

  # U_0: variance-stabilized diffuse score
  sum_wz <- sum(w * z)
  sum_w <- sum(w)
  sum_w2 <- sum(w * w)

  if (sum_w2 <= 0 || sum_w <= 0) {
    return(list(T_omnibus = -Inf, U0 = -Inf, S_best = -Inf, kappa_best = NA))
  }

  U0 <- sum_wz / sqrt(sum_w2)
  log_den1 <- log(sum_w)

  # S_kappa: normalized soft-max scores
  S_best <- -Inf
  kappa_best <- NA

  for (kappa in kappa_grid) {
    if (kappa <= 0) next

    a <- kappa * z
    a_max <- max(a)
    log_sum <- a_max + log(sum(w * exp(a - a_max)))
    S_k <- (log_sum - log_den1) / kappa

    if (S_k > S_best) {
      S_best <- S_k
      kappa_best <- kappa
    }
  }

  T_omnibus <- max(U0, S_best)

  list(T_omnibus = T_omnibus, U0 = U0, S_best = S_best, kappa_best = kappa_best)
}
