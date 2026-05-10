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
#' @param n_perm Number of pre-computed sign-flip draws (default 5000); the
#'   resulting `null_fun(b)` is deterministic in `b`.
#'
#' @export
#' @examples
#' subj <- replicate(5, array(rnorm(4 * 4 * 4), dim = c(4, 4, 4)), simplify = FALSE)
#' sf <- make_null_fun_subject_signflip(subj, seed = 1)
#' res <- octree_scan_fwer(sf$z_vol, mask = sf$mask, n_perm = 5, null_fun = sf$null_fun)
make_null_fun_subject_signflip <- function(subject_vols,
                                           mask = NULL,
                                           n_perm = 5000,
                                           seed = NULL,
                                           z_equiv = TRUE) {
  ns <- null_one_sample(subject_vols, mask = mask, n_perm = n_perm,
                        seed = seed, z_equiv = z_equiv)
  list(
    z_vol = ns$stat_vol,
    mask = ns$mask,
    null_fun = ns$draw,
    df = ns$df
  )
}
