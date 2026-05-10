#' Threshold-Free Cluster Enhancement (TFCE)
#'
#' Implementation of FSL-style TFCE with permutation-based FWER control.
#'
#' @name tfce-methods
NULL

#' TFCE Transform
#'
#' Computes the Threshold-Free Cluster Enhancement transform:
#' \deqn{TFCE(v) = \int_0^{h_v} e(h)^E \cdot h^H \, dh}
#'
#' where \eqn{e(h)} is the extent of the cluster containing voxel \eqn{v}
#' at threshold \eqn{h}.
#'
#' @param stat_vol 3D Z-equivalent statistic map (array or `NeuroVol`).
#' @param mask Optional logical 3D mask. Defaults to finite voxels in
#'   `stat_vol`.
#' @param H Height exponent (default 2.0).
#' @param E Extent exponent (default 0.5).
#' @param dh Threshold increment (default 0.1).
#' @param tail One of `"pos"`, `"neg"`, `"two"`. For `"two"` the result is
#'   the per-voxel maximum of the positive- and negative-tail TFCE maps.
#'
#' @return Volume (same type as `stat_vol`) with TFCE-transformed values.
#'
#' @references Smith & Nichols (2009), NeuroImage 44(1), 83-98.
#' @seealso [tfce_fwer()].
#' @export
tfce_transform <- function(stat_vol, mask = NULL,
                           H = 2.0, E = 0.5, dh = 0.1,
                           tail = c("pos", "neg", "two")) {
  tail <- match.arg(tail)

  if (is.null(mask)) {
    mask <- is.finite(as.numeric(stat_vol))
  }

  dims <- dim(stat_vol)

  if (tail == "two") {
    pos_vec <- .tfce_from_vec(as.numeric(stat_vol), dims, mask, H, E, dh, "pos")
    neg_vec <- .tfce_from_vec(as.numeric(stat_vol), dims, mask, H, E, dh, "neg")
    return(.array_to_vol(pmax(pos_vec, neg_vec), dims, stat_vol))
  }

  vals <- .tfce_from_vec(as.numeric(stat_vol), dims, mask, H, E, dh, tail)
  .array_to_vol(vals, dims, stat_vol)
}

#' TFCE with Permutation FWER
#'
#' Computes TFCE with permutation-based FWER control using the maxT
#' procedure. The null distribution is supplied via a [`null_spec`]; for
#' subject-level inference build one with [null_one_sample()] /
#' [null_two_sample()] / [null_glm()]. For a fast voxel-level fallback use
#' [null_voxel_signflip()] (note: voxel sign-flips are calibrated only when
#' the field is uncorrelated under the null).
#'
#' @param stat_vol 3D Z-equivalent statistic map.
#' @param mask Optional logical 3D mask.
#' @param alpha Significance level (default 0.05).
#' @param tail One of `"pos"`, `"neg"`, `"two"`.
#' @param H Height exponent (default 2.0).
#' @param E Extent exponent (default 0.5).
#' @param dh Threshold increment (default 0.1).
#' @param null_spec A [`null_spec`] supplying `n_perm` deterministic null
#'   draws.
#' @param null_fun Optional function \code{function(b)} returning a null draw
#'   for permutation \code{b}. Kept for compatibility with the package's older
#'   permutation API; prefer \code{null_spec} for new code.
#' @param n_perm Number of permutations to consume from `null_spec`. Defaults
#'   to `null_spec$n_perm`.
#' @param allow_voxel_signflip Logical; if \code{TRUE}, allows a fallback to
#'   independent voxel sign-flips when neither \code{null_spec} nor
#'   \code{null_fun} is supplied. This fallback is fast but can be
#'   anti-conservative for spatially correlated fields.
#' @param seed Optional random seed; affects only RNG state outside `null_spec`
#'   (the spec itself is deterministic by index).
#'
#' @return A [`neurothresh_result`] with `method = "tfce_fwer"`. `pvals` is the
#'   voxelwise FWER-corrected p-value volume; `extras$tfce` holds the observed
#'   TFCE map; `null_max` is the vector of per-permutation max-TFCE values.
#'
#' @seealso [tfce_transform()], [null_voxel_signflip()], [null_one_sample()].
#' @export
tfce_fwer <- function(stat_vol, mask = NULL,
                      alpha = 0.05, tail = c("pos", "neg", "two"),
                      H = 2.0, E = 0.5, dh = 0.1,
                      null_spec = NULL, null_fun = NULL, n_perm = NULL,
                      allow_voxel_signflip = FALSE,
                      seed = NULL) {
  tail <- match.arg(tail)
  call_capture <- match.call()

  if (!is.null(seed)) set.seed(seed)

  if (is.null(mask)) {
    mask <- is.finite(as.numeric(stat_vol))
  }

  dims <- dim(stat_vol)
  n <- prod(dims)
  mask_idx <- which(mask)

  if (!is.null(null_spec) && !is.null(null_fun)) {
    stop("provide only one of null_spec or null_fun")
  }
  if (!is.null(null_spec)) {
    if (!is_null_spec(null_spec)) {
      stop("null_spec must be a neurothresh_null_spec (see ?null_spec)")
    }
    if (is.null(n_perm)) n_perm <- null_spec$n_perm
    if (n_perm > null_spec$n_perm) {
      stop(sprintf("n_perm = %d exceeds null_spec$n_perm = %d",
                   n_perm, null_spec$n_perm))
    }
    draw_fun <- null_spec$draw
    null_kind <- null_spec$kind
  } else {
    if (is.null(n_perm)) n_perm <- 5000L
    if (is.null(null_fun)) {
      if (!isTRUE(allow_voxel_signflip)) {
        stop(
          "tfce_fwer() requires null_spec for calibrated inference; provide null_fun or null_spec. ",
          "Use make_null_fun_subject_signflip(... )$null_fun or null_one_sample(subjects) ",
          "(recommended), or set allow_voxel_signflip=TRUE for a fast but potentially ",
          "anti-conservative fallback."
        )
      }
      null_fun <- null_voxel_signflip(stat_vol, mask = mask, n_perm = n_perm, seed = seed)$draw
      null_kind <- "voxel_signflip"
    } else {
      null_kind <- "custom_null_fun"
    }
    draw_fun <- null_fun
  }

  if (tail == "two") {
    return(.tfce_fwer_twosided(stat_vol, mask, n_perm, H, E, dh, alpha,
                                draw_fun, null_kind, call_capture))
  }

  tfce_obs <- tfce_transform(stat_vol, mask, H, E, dh, tail = tail)
  tfce_obs_vec <- as.numeric(tfce_obs)

  max_null <- numeric(n_perm)
  for (b in seq_len(n_perm)) {
    z_perm <- .coerce_null_draw(draw_fun(b), n_full = n, mask_idx = mask_idx)
    tfce_perm <- .tfce_from_vec(z_perm, dims, mask, H, E, dh, tail)
    max_null[b] <- max(tfce_perm[mask_idx], na.rm = TRUE)
  }

  threshold <- fwer_threshold_from_maxnull(max_null, alpha = alpha)
  reject_vec <- tfce_obs_vec > threshold & mask
  reject <- array(reject_vec, dim = dims)

  p_corr <- rep(1, n)
  for (v in mask_idx) {
    p_corr[v] <- (1 + sum(max_null >= tfce_obs_vec[v])) / (n_perm + 1)
  }
  pvals <- .array_to_vol(p_corr, dims, stat_vol)

  out <- new_neurothresh_result(
    method = "tfce_fwer",
    stat_vol = stat_vol,
    mask = mask,
    threshold = threshold,
    reject = reject,
    pvals = pvals,
    null_max = max_null,
    params = list(alpha = alpha, tail = tail, H = H, E = E, dh = dh,
                  n_perm = n_perm, null_kind = null_kind,
                  allow_voxel_signflip = allow_voxel_signflip),
    call = call_capture,
    extras = list(tfce = tfce_obs)
  )
  .tfce_add_legacy_aliases(out, max_null, reject)
}

#' Two-sided TFCE FWER
#' @keywords internal
.tfce_fwer_twosided <- function(stat_vol, mask, n_perm, H, E, dh, alpha,
                                draw_fun, null_kind, call_capture) {
  dims <- dim(stat_vol)
  n <- prod(dims)
  mask_idx <- which(mask)

  tfce_pos <- as.numeric(tfce_transform(stat_vol, mask, H, E, dh, tail = "pos"))
  tfce_neg <- as.numeric(tfce_transform(stat_vol, mask, H, E, dh, tail = "neg"))
  tfce_2   <- pmax(tfce_pos, tfce_neg)

  max_null <- numeric(n_perm)
  for (b in seq_len(n_perm)) {
    z_perm <- .coerce_null_draw(draw_fun(b), n_full = n, mask_idx = mask_idx)
    tfce_p <- .tfce_from_vec(z_perm, dims, mask, H, E, dh, "pos")
    tfce_n <- .tfce_from_vec(-z_perm, dims, mask, H, E, dh, "pos")
    max_null[b] <- max(pmax(tfce_p, tfce_n)[mask_idx], na.rm = TRUE)
  }

  threshold <- fwer_threshold_from_maxnull(max_null, alpha = alpha)
  reject_vec <- tfce_2 > threshold & mask
  reject <- array(reject_vec, dim = dims)

  p_corr <- rep(1, n)
  for (v in mask_idx) {
    p_corr[v] <- (1 + sum(max_null >= tfce_2[v])) / (n_perm + 1)
  }
  pvals <- .array_to_vol(p_corr, dims, stat_vol)

  reject_pos <- .array_to_vol(tfce_pos > threshold & tfce_pos >= tfce_neg & mask, dims, stat_vol)
  reject_neg <- .array_to_vol(tfce_neg > threshold & tfce_neg > tfce_pos & mask, dims, stat_vol)

  out <- new_neurothresh_result(
    method = "tfce_fwer",
    stat_vol = stat_vol,
    mask = mask,
    threshold = threshold,
    reject = reject,
    pvals = pvals,
    null_max = max_null,
    params = list(alpha = alpha, tail = "two", H = H, E = E, dh = dh,
                  n_perm = n_perm, null_kind = null_kind),
    call = call_capture,
    extras = list(
      tfce_pos = .array_to_vol(tfce_pos, dims, stat_vol),
      tfce_neg = .array_to_vol(tfce_neg, dims, stat_vol),
      tfce_two = .array_to_vol(tfce_2, dims, stat_vol),
      reject_pos = reject_pos,
      reject_neg = reject_neg
    )
  )
  .tfce_add_legacy_aliases(out, max_null, reject,
                           sig_pos = reject_pos, sig_neg = reject_neg)
}

#' Add compatibility fields used by earlier list-style TFCE callers.
#' @keywords internal
#' @noRd
.tfce_add_legacy_aliases <- function(out, max_null, sig_mask,
                                     sig_pos = NULL, sig_neg = NULL) {
  out$max_null <- max_null
  out$sig_mask <- sig_mask
  if (!is.null(sig_pos)) out$sig_pos <- sig_pos
  if (!is.null(sig_neg)) out$sig_neg <- sig_neg
  out
}

#' TFCE from vector (internal helper)
#' @keywords internal
.tfce_from_vec <- function(z_vec, dims, mask, H, E, dh, tail) {
  z <- z_vec
  if (tail == "neg") z <- -z

  n <- prod(dims)
  z[!mask] <- 0
  z[z < 0] <- 0

  h_max <- max(z, na.rm = TRUE)
  if (h_max <= 0) return(rep(0, n))

  hs <- seq(0, h_max, by = dh)
  tfce_vals <- rep(0, n)

  for (h in hs) {
    if (h == 0) next
    supra <- z >= h
    clusters <- .find_clusters(supra, dims)
    if (length(clusters$sizes) == 0) next

    for (cid in seq_along(clusters$sizes)) {
      extent <- clusters$sizes[cid]
      voxels <- which(clusters$labels == cid)
      tfce_vals[voxels] <- tfce_vals[voxels] + (extent^E) * (h^H) * dh
    }
  }

  tfce_vals
}

#' Convert array to vol-like object (placeholder)
#' @keywords internal
.array_to_vol <- function(vals, dims, template) {
  if (length(vals) != prod(dims)) {
    stop("vals length does not match dims")
  }

  out <- template
  ok <- TRUE
  tryCatch({
    out[] <- vals
  }, error = function(e) {
    ok <<- FALSE
  })

  if (ok) {
    return(out)
  }

  array(vals, dim = dims)
}
