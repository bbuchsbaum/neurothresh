#' Cluster-Level FDR Control
#'
#' Implements topological FDR control at the cluster level using
#' Benjamini-Hochberg procedure applied to cluster p-values.
#'
#' @name cluster-fdr
NULL

#' Cluster-FDR with RFT or Permutation P-values
#'
#' Controls the False Discovery Rate at the cluster level by applying
#' the Benjamini-Hochberg procedure to cluster-extent p-values.
#'
#' @param stat_vol 3D Z-equivalent statistic map (array or `NeuroVol`).
#' @param mask Optional logical 3D mask. Defaults to finite voxels in
#'   `stat_vol`.
#' @param q FDR level (default 0.05).
#' @param tail One of `"pos"`, `"neg"`, `"two"`. For `"two"`, both tails are
#'   pooled into a single BH ordering.
#' @param cluster_thresh Cluster-forming threshold (Z-score).
#' @param p_method Method for computing cluster p-values: `"rft"` (Random
#'   Field Theory) or `"perm"` (permutation).
#' @param smoothness A [`smoothness_fwhm()`] (or numeric mm) — required when
#'   `p_method = "rft"`. Ignored otherwise.
#' @param null_spec A [`null_spec`] supplying null draws — required when
#'   `p_method = "perm"`.
#' @param n_perm Number of permutations to consume from `null_spec`. Defaults
#'   to `null_spec$n_perm`.
#' @param seed Optional random seed.
#' @param fwhm_mm Deprecated alias for `smoothness`.
#' @param null_fun Deprecated bare-function null source; prefer `null_spec`.
#' @param allow_voxel_signflip Deprecated; if `TRUE` and `null_spec`/`null_fun`
#'   are missing, falls back to independent voxel sign-flips.
#' @param two_sided_policy Deprecated. The `"two"` tail now uses BH on the
#'   pooled set of clusters (the previous `"BH_all"` behaviour). Pass
#'   `"split_q"` only if you need the legacy q/2 split per tail.
#'
#' @return A [`neurothresh_result`] with `method = "cluster_fdr"` and the
#'   cluster table in `$clusters`. `$threshold` is `cluster_thresh`.
#'
#' @references Chumbley & Friston (2009), NeuroImage 44(1), 62-70.
#' @seealso [rft_cluster_fwer()], [tfce_fwer()], [hier_scan()].
#' @export
cluster_fdr <- function(stat_vol, mask = NULL,
                        q = 0.05, tail = c("pos", "neg", "two"),
                        cluster_thresh = 3.0,
                        p_method = c("rft", "perm"),
                        smoothness = NULL,
                        null_spec = NULL,
                        n_perm = NULL,
                        seed = NULL,
                        fwhm_mm = NULL,
                        null_fun = NULL,
                        allow_voxel_signflip = FALSE,
                        two_sided_policy = c("BH_all", "split_q")) {

  p_method <- match.arg(p_method)
  tail <- match.arg(tail)
  two_sided_policy <- match.arg(two_sided_policy)
  call_capture <- match.call()

  if (is.null(smoothness) && !is.null(fwhm_mm)) smoothness <- fwhm_mm

  if (p_method == "rft") {
    if (is.null(smoothness)) {
      stop("p_method='rft' requires smoothness (a smoothness_fwhm() or numeric mm)")
    }
    fwhm_mm_used <- .resolve_smoothness_fwhm_mm(smoothness)
  } else {
    fwhm_mm_used <- NULL
  }

  if (is.null(mask)) {
    mask <- is.finite(as.numeric(stat_vol))
  }
  dims <- dim(stat_vol)
  if (is.null(dims) || length(dims) != 3L) {
    stop("stat_vol must be a 3D volume-like object")
  }
  n <- prod(dims)
  mask_idx <- which(mask)

  draw_fun <- NULL
  if (p_method == "perm") {
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
    } else {
      if (is.null(n_perm)) n_perm <- 5000L
      if (is.null(null_fun)) {
        if (!isTRUE(allow_voxel_signflip)) {
          stop(
            "cluster_fdr(p_method='perm') requires null_spec (or null_fun) for ",
            "calibrated inference. Use null_one_sample(subjects) (recommended), ",
            "or set allow_voxel_signflip=TRUE for a fast but potentially ",
            "anti-conservative voxel-sign-flip fallback."
          )
        }
        z_full <- as.numeric(stat_vol)
        draw_fun <- null_voxel_signflip(stat_vol, mask = mask, n_perm = n_perm,
                                        seed = seed)$draw
      } else {
        draw_fun <- null_fun
      }
    }
  }

  if (tail == "two") {
    return(.cluster_fdr_twosided(stat_vol, mask, cluster_thresh, q,
                                 p_method, fwhm_mm_used, n_perm,
                                 two_sided_policy, seed, draw_fun, call_capture))
  }

  z <- as.numeric(stat_vol)
  if (tail == "neg") z <- -z

  supra <- (z >= cluster_thresh) & mask
  clusters <- .find_clusters(supra, dims)

  empty <- length(clusters$sizes) == 0L

  if (empty) {
    cluster_table <- data.frame(
      cluster_id = integer(0), size_voxels = integer(0),
      p_unc = numeric(0), significant = logical(0), tail = character(0)
    )
    out <- new_neurothresh_result(
      method = "cluster_fdr",
      stat_vol = stat_vol, mask = mask,
      threshold = cluster_thresh,
      reject = array(FALSE, dim = dims),
      pvals = array(NA_real_, dim = dims),
      clusters = cluster_table,
      params = list(q = q, tail = tail, p_method = p_method,
                    cluster_thresh = cluster_thresh,
                    smoothness = smoothness, fwhm_mm = fwhm_mm_used,
                    n_perm = n_perm),
      call = call_capture
    )
    return(.cfdr_add_legacy_aliases(out, clusters, cluster_table))
  }

  if (p_method == "rft") {
    vox_size <- .vox_size_mm(stat_vol)
    p_unc <- .cluster_pvals_rft(clusters$sizes, cluster_thresh,
                                fwhm_mm_used, vox_size)
  } else {
    p_unc <- .cluster_pvals_perm(
      dims, mask, cluster_thresh, clusters, n_perm, seed,
      draw_fun, n, mask_idx, tail
    )
  }

  rejected <- .bh_reject(p_unc, q)

  cluster_table <- data.frame(
    cluster_id = seq_along(clusters$sizes),
    size_voxels = clusters$sizes,
    p_unc = p_unc,
    significant = rejected,
    tail = tail
  )

  reject_arr <- mask & FALSE
  for (id in which(rejected)) {
    reject_arr[clusters$labels == id] <- TRUE
  }
  reject <- array(reject_arr, dim = dims)

  out <- new_neurothresh_result(
    method = "cluster_fdr",
    stat_vol = stat_vol, mask = mask,
    threshold = cluster_thresh,
    reject = reject,
    pvals = array(NA_real_, dim = dims),
    clusters = cluster_table,
    params = list(q = q, tail = tail, p_method = p_method,
                  cluster_thresh = cluster_thresh,
                  smoothness = smoothness, fwhm_mm = fwhm_mm_used,
                  n_perm = n_perm),
    call = call_capture
  )
  .cfdr_add_legacy_aliases(out, clusters, cluster_table)
}

#' Two-sided cluster FDR
#' @keywords internal
.cluster_fdr_twosided <- function(stat_vol, mask, cluster_thresh, q,
                                  p_method, fwhm_mm, n_perm,
                                  two_sided_policy, seed, draw_fun,
                                  call_capture) {
  q_tail <- if (two_sided_policy == "split_q") q / 2 else q
  dims <- dim(stat_vol)

  pos <- cluster_fdr(stat_vol, mask, q = q_tail, tail = "pos",
                     cluster_thresh = cluster_thresh,
                     p_method = p_method,
                     smoothness = if (p_method == "rft") fwhm_mm else NULL,
                     null_spec = NULL,
                     null_fun = draw_fun,
                     n_perm = n_perm, seed = seed)
  neg <- cluster_fdr(stat_vol, mask, q = q_tail, tail = "neg",
                     cluster_thresh = cluster_thresh,
                     p_method = p_method,
                     smoothness = if (p_method == "rft") fwhm_mm else NULL,
                     null_spec = NULL,
                     null_fun = draw_fun,
                     n_perm = n_perm, seed = seed)

  pos_table <- pos$clusters
  neg_table <- neg$clusters

  if (two_sided_policy == "split_q") {
    sig_pos <- as.logical(pos$reject)
    sig_neg <- as.logical(neg$reject)
    sig_arr <- sig_pos | sig_neg
  } else {
    p_all <- c(pos_table$p_unc, neg_table$p_unc)
    if (length(p_all) > 0L) {
      rej_all <- .bh_reject(p_all, q)
      n_pos <- nrow(pos_table)
      n_all <- length(rej_all)

      rej_pos <- if (n_pos > 0L) rej_all[seq_len(n_pos)] else logical(0)
      rej_neg <- if (n_pos < n_all) rej_all[(n_pos + 1L):n_all] else logical(0)

      pos_table$significant <- rej_pos
      neg_table$significant <- rej_neg

      pos_labels <- pos$cluster_labels
      neg_labels <- neg$cluster_labels

      sig_pos <- mask & FALSE
      if (!is.null(pos_labels)) {
        for (id in which(rej_pos)) sig_pos[pos_labels == id] <- TRUE
      }
      sig_neg <- mask & FALSE
      if (!is.null(neg_labels)) {
        for (id in which(rej_neg)) sig_neg[neg_labels == id] <- TRUE
      }
      sig_arr <- sig_pos | sig_neg
    } else {
      sig_pos <- mask & FALSE
      sig_neg <- mask & FALSE
      sig_arr <- mask & FALSE
    }
  }

  table <- rbind(pos_table, neg_table)
  reject <- array(sig_arr, dim = dims)

  out <- new_neurothresh_result(
    method = "cluster_fdr",
    stat_vol = stat_vol, mask = mask,
    threshold = cluster_thresh,
    reject = reject,
    pvals = array(NA_real_, dim = dims),
    clusters = table,
    params = list(q = q, tail = "two", p_method = p_method,
                  cluster_thresh = cluster_thresh,
                  fwhm_mm = fwhm_mm, n_perm = n_perm,
                  two_sided_policy = two_sided_policy),
    call = call_capture,
    extras = list(
      reject_pos = array(sig_pos, dim = dims),
      reject_neg = array(sig_neg, dim = dims)
    )
  )
  out$sig_mask <- reject
  out$sig_pos <- array(sig_pos, dim = dims)
  out$sig_neg <- array(sig_neg, dim = dims)
  out$table <- table
  out$n_clusters <- sum(table$significant, na.rm = TRUE)
  out
}

#' RFT cluster p-values
#' @keywords internal
.cluster_pvals_rft <- function(sizes, u0, fwhm_mm, vox_size) {
  if (length(fwhm_mm) == 1) fwhm_mm <- rep(fwhm_mm, 3)

  s_resel <- sizes * prod(vox_size) / prod(fwhm_mm)
  .rft_cluster_extent_tail(s_resel, u0)
}

#' Permutation cluster p-values
#' @keywords internal
.cluster_pvals_perm <- function(dims, mask, u0, obs_clusters, n_perm, seed,
                                null_fun, n_full, mask_idx, tail = "pos") {
  if (!is.null(seed)) set.seed(seed)
  obs_sizes <- obs_clusters$sizes

  null_sizes <- numeric(0)
  for (b in seq_len(n_perm)) {
    z_perm <- .coerce_null_draw(null_fun(b), n_full = n_full, mask_idx = mask_idx)
    if (tail == "neg") z_perm <- -z_perm
    supra <- (z_perm >= u0) & mask
    perm_clusters <- .find_clusters(supra, dims)
    if (length(perm_clusters$sizes) > 0) {
      null_sizes <- c(null_sizes, perm_clusters$sizes)
    }
  }

  if (length(null_sizes) == 0) {
    p_unc <- rep(1 / (n_perm + 1), length(obs_sizes))
  } else {
    p_unc <- vapply(obs_sizes, function(s) {
      (1 + sum(null_sizes >= s)) / (1 + length(null_sizes))
    }, numeric(1))
  }
  p_unc
}

#' Benjamini-Hochberg rejection
#' @keywords internal
.bh_reject <- function(p, q) {
  m <- length(p)
  if (m == 0) return(logical(0))

  o <- order(p)
  ps <- p[o]
  thr <- (seq_len(m) / m) * q
  k <- max(which(ps <= thr), 0L)

  rej <- rep(FALSE, m)
  if (k > 0) rej[o[seq_len(k)]] <- TRUE
  rej
}

#' Attach legacy aliases to a cluster_fdr result
#'
#' Allows older callers using `$sig_mask`, `$table`, `$n_clusters`, and
#' `$cluster_labels` to keep working.
#' @keywords internal
.cfdr_add_legacy_aliases <- function(out, clusters, table) {
  out$sig_mask <- out$reject
  out$table <- table
  out$n_clusters <- sum(table$significant, na.rm = TRUE)
  out$cluster_labels <- clusters$labels
  out
}
