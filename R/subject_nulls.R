#' Subject-level sign-flip null generator (one-sample)
#'
#' Constructs an observed Z-equivalent group map and a \code{null_fun(b)} that
#' generates sign-flip randomizations by flipping each subject map by a single
#' \eqn{\pm 1} (constant across voxels). This preserves the spatial covariance
#' structure within each subject map and is the recommended way to drive
#' resampling-based inference for regional statistics.
#'
#' @param subject_vols List of subject-level 3D volumes (arrays or \code{NeuroVol}s)
#'   with identical dimensions.
#' @param mask Optional logical mask aligned to the volumes. If NULL, uses all
#'   finite voxels common to all subjects.
#' @param seed Optional RNG seed.
#' @param z_equiv If TRUE (default), returns Z-equivalent vectors using a normal
#'   approximation to the one-sample t statistic with df = n_subjects - 1.
#'
#' @return A list with components:
#'   \describe{
#'     \item{z_vol}{Observed group Z-equivalent volume (array-like, same dims).}
#'     \item{mask}{Logical mask used.}
#'     \item{null_fun}{Function \code{function(b)} returning a length-\code{n_mask}
#'       Z-equivalent vector for permutation index \code{b}.}
#'     \item{df}{Degrees of freedom used for the t-to-Z conversion.}
#'   }
#'
#' @export
#' @examples
#' subj <- replicate(10, array(rnorm(20 * 20 * 20), dim = c(20, 20, 20)), simplify = FALSE)
#' sf <- make_null_fun_subject_signflip(subj, seed = 1)
#' res <- octree_scan_fwer(sf$z_vol, mask = sf$mask, n_perm = 200, null_fun = sf$null_fun)
make_null_fun_subject_signflip <- function(subject_vols,
                                           mask = NULL,
                                           seed = NULL,
                                           z_equiv = TRUE) {
  if (!is.null(seed)) set.seed(seed)
  if (!is.list(subject_vols) || length(subject_vols) < 2) {
    stop("subject_vols must be a list of length >= 2")
  }

  dims <- dim(subject_vols[[1]])
  if (length(dims) != 3) stop("subject_vols[[1]] must be 3D")

  n_subj <- length(subject_vols)
  subj_vecs <- vector("list", n_subj)

  for (i in seq_len(n_subj)) {
    vi <- as.numeric(subject_vols[[i]])
    if (!identical(dim(subject_vols[[i]]), dims)) stop("All subject_vols must have identical dimensions")
    subj_vecs[[i]] <- vi
  }

  if (is.null(mask)) {
    ok <- rep(TRUE, prod(dims))
    for (i in seq_len(n_subj)) {
      ok <- ok & is.finite(subj_vecs[[i]])
    }
    mask <- array(ok, dim = dims)
  } else {
    mask <- as.logical(mask)
  }

  mask_idx <- which(mask)
  n_mask <- length(mask_idx)
  if (n_mask == 0) stop("Mask is empty")

  X <- matrix(NA_real_, nrow = n_subj, ncol = n_mask)
  for (i in seq_len(n_subj)) {
    X[i, ] <- subj_vecs[[i]][mask_idx]
  }

  sumsq <- colSums(X * X)
  df <- n_subj - 1

  t_from_signs <- function(signs) {
    s <- as.numeric(crossprod(signs, X)) # length n_mask
    mean <- s / n_subj
    var <- (sumsq - n_subj * mean * mean) / df
    var[var < 0] <- 0
    se <- sqrt(var / n_subj)
    t <- mean / se
    t[!is.finite(t)] <- 0
    t
  }

  z_from_t <- function(t) {
    if (!z_equiv) return(t)
    p <- stats::pt(t, df = df, lower.tail = FALSE)
    p <- pmin(pmax(p, .Machine$double.xmin), 1 - .Machine$double.eps)
    stats::qnorm(p, lower.tail = FALSE)
  }

  # Observed: all +1
  t_obs <- t_from_signs(rep(1, n_subj))
  z_obs <- z_from_t(t_obs)

  z_vol <- subject_vols[[1]]
  z_vol[] <- NA_real_
  z_vol[mask] <- z_obs

  null_fun <- function(b) {
    signs <- sample(c(-1, 1), n_subj, replace = TRUE)
    z_from_t(t_from_signs(signs))
  }

  list(
    z_vol = z_vol,
    mask = mask,
    null_fun = null_fun,
    df = df
  )
}
