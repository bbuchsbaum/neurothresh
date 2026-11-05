#' Empirical FDR for LISA filtered maps
#'
#' Estimates local FDR/q-values in the filtered-statistic space using the
#' empirical filtered null distribution and the empirical observed-map mixture
#' distribution.
#'
#' @param filtered_obs Numeric observed filtered statistic vector or volume.
#' @param filtered_nulls Null filtered statistics. Accepts a matrix with one row
#'   per permutation and one column per voxel, or a list of numeric vectors.
#' @param mask Optional logical mask when `filtered_obs` is a volume.
#' @param tail One of `"pos"`, `"neg"`, or `"two"`.
#' @param q Target FDR level used to compute `threshold` and `reject`.
#' @param p0 Conservative null proportion. Defaults to 1.
#' @param monotone Whether to enforce monotone q-values over evidence strength.
#'
#' @return List with vector `fdr`, scalar `threshold`, logical `reject`, and an
#'   empirical `curve` data frame.
#' @export
lisa_fdr_map <- function(filtered_obs,
                         filtered_nulls,
                         mask = NULL,
                         tail = c("pos", "neg", "two"),
                         q = 0.05,
                         p0 = 1,
                         monotone = TRUE) {
  tail <- match.arg(tail)
  if (!is.numeric(q) || length(q) != 1L || q <= 0 || q >= 1) stop("q must be in (0, 1)")
  if (!is.numeric(p0) || length(p0) != 1L || p0 <= 0 || p0 > 1) stop("p0 must be in (0, 1]")

  obs <- as.numeric(filtered_obs)
  if (!is.null(mask) && length(mask) == length(obs)) {
    obs <- obs[as.logical(mask)]
  }
  if (any(!is.finite(obs))) stop("filtered_obs contains non-finite values")

  null_mat <- .lisa_nulls_to_matrix(filtered_nulls, n_vox = length(obs))
  if (any(!is.finite(null_mat))) stop("filtered_nulls contain non-finite values")

  obs_score <- .lisa_tail_score(obs, tail)
  null_score <- .lisa_tail_score(as.numeric(null_mat), tail)

  ord <- order(obs_score, decreasing = TRUE)
  score_ord <- obs_score[ord]
  n0 <- length(null_score)
  nmix <- length(obs_score)
  p_null <- .lisa_count_ge(null_score, score_ord) / n0
  p_mix <- .lisa_count_ge(obs_score, score_ord) / nmix
  raw_ord <- pmin(1, p0 * p_null / pmax(p_mix, .Machine$double.eps))
  q_ord <- if (isTRUE(monotone)) rev(cummin(rev(raw_ord))) else raw_ord
  fdr <- numeric(length(obs_score))
  fdr[ord] <- pmin(1, pmax(0, q_ord))

  reject <- fdr <= q
  threshold <- if (any(reject)) min(obs_score[reject]) else Inf
  curve <- data.frame(
    score = score_ord,
    fdr = q_ord,
    raw_fdr = raw_ord
  )
  list(fdr = fdr, threshold = threshold, reject = reject, curve = curve)
}

#' @keywords internal
.lisa_nulls_to_matrix <- function(filtered_nulls, n_vox) {
  if (is.matrix(filtered_nulls)) {
    if (ncol(filtered_nulls) != n_vox) stop("filtered_nulls matrix must have one column per observed voxel")
    return(filtered_nulls)
  }
  if (is.list(filtered_nulls)) {
    if (!length(filtered_nulls)) stop("filtered_nulls list must be non-empty")
    mat <- matrix(NA_real_, nrow = length(filtered_nulls), ncol = n_vox)
    for (i in seq_along(filtered_nulls)) {
      v <- as.numeric(filtered_nulls[[i]])
      if (length(v) != n_vox) stop("each filtered_nulls list element must have length n_vox")
      mat[i, ] <- v
    }
    return(mat)
  }
  stop("filtered_nulls must be a matrix or list")
}

#' @keywords internal
.lisa_tail_score <- function(x, tail) {
  switch(tail,
         pos = x,
         neg = -x,
         two = abs(x))
}

#' @noRd
.lisa_count_ge <- function(x, thresholds) {
  # left.open counts values strictly below each threshold, retaining every tie
  # in the >= tail. Sort once rather than scanning all voxels per threshold.
  length(x) - findInterval(thresholds, sort(x), left.open = TRUE)
}

#' @keywords internal
.lisa_fdr_counter_init <- function(filtered_obs, tail) {
  obs <- as.numeric(filtered_obs)
  if (any(!is.finite(obs))) stop("filtered_obs contains non-finite values")
  obs_score <- .lisa_tail_score(obs, tail)
  ord <- order(obs_score, decreasing = TRUE)
  score_ord <- obs_score[ord]
  mix_counts <- .lisa_count_ge(obs_score, score_ord)
  list(
    obs_score = obs_score,
    ord = ord,
    score_ord = score_ord,
    mix_counts = mix_counts,
    null_counts = numeric(length(score_ord)),
    null_total = 0
  )
}

#' @keywords internal
.lisa_fdr_counter_update <- function(counter, filtered_null, tail) {
  null_score <- .lisa_tail_score(as.numeric(filtered_null), tail)
  if (any(!is.finite(null_score))) stop("filtered null map contains non-finite values")
  counter$null_counts <- counter$null_counts +
    .lisa_count_ge(null_score, counter$score_ord)
  counter$null_total <- counter$null_total + length(null_score)
  counter
}

#' @keywords internal
.lisa_fdr_counter_finalize <- function(counter, q, p0, monotone = TRUE) {
  if (!is.numeric(q) || length(q) != 1L || q <= 0 || q >= 1) stop("q must be in (0, 1)")
  if (!is.numeric(p0) || length(p0) != 1L || p0 <= 0 || p0 > 1) stop("p0 must be in (0, 1]")
  if (counter$null_total < 1L) stop("no filtered null values accumulated")

  p_null <- counter$null_counts / counter$null_total
  p_mix <- counter$mix_counts / length(counter$obs_score)
  raw_ord <- pmin(1, p0 * p_null / pmax(p_mix, .Machine$double.eps))
  q_ord <- if (isTRUE(monotone)) rev(cummin(rev(raw_ord))) else raw_ord

  fdr <- numeric(length(counter$obs_score))
  fdr[counter$ord] <- pmin(1, pmax(0, q_ord))
  reject <- fdr <= q
  threshold <- if (any(reject)) min(counter$obs_score[reject]) else Inf
  curve <- data.frame(
    score = counter$score_ord,
    fdr = q_ord,
    raw_fdr = raw_ord
  )
  list(fdr = fdr, threshold = threshold, reject = reject, curve = curve)
}
