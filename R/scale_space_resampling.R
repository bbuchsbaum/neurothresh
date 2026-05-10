#' Subject-Level One-Sample Sign-Flip Null Generator
#'
#' Builds an observed one-sample statistic map and deterministic subject-level
#' sign-flip null draws for scale-space and other resampling methods.
#'
#' @param subject_images List of subject-level 3D images, or a subject-by-voxel
#'   matrix with \code{template}/\code{mask}.
#' @param mask Optional logical mask.
#' @param n_perm Number of sign-flip draws.
#' @param seed Optional random seed.
#' @param z_equiv If TRUE, return Z-equivalent maps; otherwise return t maps.
#' @param template Optional template for matrix inputs.
#'
#' @return List with \code{stat_vol}, \code{mask}, deterministic
#'   \code{null_fun}, and metadata.
#'
#' @export
make_null_fun_onesample <- function(subject_images,
                                    mask = NULL,
                                    n_perm = 1000,
                                    seed = NULL,
                                    z_equiv = TRUE,
                                    template = NULL) {
  prep <- as_subject_matrix(subject_images, mask = mask, template = template)
  X <- prep$X
  n <- nrow(X)
  if (n_perm < 1L) stop("n_perm must be >= 1")

  if (!is.null(seed)) set.seed(seed)
  signs <- matrix(sample(c(-1, 1), n * n_perm, replace = TRUE),
                  nrow = n_perm, ncol = n)

  stat_from_signs <- function(s) {
    t <- .ss_one_sample_t(X * as.numeric(s))
    if (isTRUE(z_equiv)) .ss_t_to_z(t, df = n - 1) else t
  }

  z_obs <- stat_from_signs(rep(1, n))
  stat_vol <- .ss_stat_template(prep, z_obs)

  null_fun <- function(b) {
    if (b < 1L || b > n_perm) stop("Permutation index out of range")
    stat_from_signs(signs[b, ])
  }

  list(
    stat_vol = stat_vol,
    z_vol = stat_vol,
    mask = prep$mask,
    null_fun = null_fun,
    prep = prep,
    df = n - 1,
    signs = signs,
    statistic = if (isTRUE(z_equiv)) "z" else "t",
    design = "one_sample"
  )
}

#' Two-Group Label-Permutation Null Generator
#'
#' Builds an observed two-group statistic map and deterministic fixed-size
#' label-permutation null draws.
#'
#' @param subject_images List of subject-level 3D images, or a subject-by-voxel
#'   matrix with \code{template}/\code{mask}.
#' @param group Two-level group label.
#' @param mask Optional logical mask.
#' @param n_perm Number of label permutations.
#' @param seed Optional random seed.
#' @param variance "welch" (default) or "pooled".
#' @param z_equiv If TRUE, return Z-equivalent maps; otherwise return t maps.
#' @param template Optional template for matrix inputs.
#'
#' @return List with \code{stat_vol}, \code{mask}, deterministic
#'   \code{null_fun}, and metadata.
#'
#' @export
make_null_fun_twosample <- function(subject_images,
                                    group,
                                    mask = NULL,
                                    n_perm = 1000,
                                    seed = NULL,
                                    variance = c("welch", "pooled"),
                                    z_equiv = TRUE,
                                    template = NULL) {
  variance <- match.arg(variance)
  prep <- as_subject_matrix(subject_images, mask = mask, template = template)
  X <- prep$X
  group <- as.factor(group)
  if (length(group) != nrow(X)) stop("length(group) must equal number of images")
  if (anyNA(group)) stop("group contains missing labels")
  if (nlevels(group) != 2L) stop("group must have exactly two levels")
  if (n_perm < 1L) stop("n_perm must be >= 1")

  if (!is.null(seed)) set.seed(seed)
  perm_index <- replicate(n_perm, sample.int(nrow(X)), simplify = FALSE)

  stat_from_group <- function(g) {
    t <- .ss_two_sample_t(X, g, variance = variance)
    if (!isTRUE(z_equiv)) return(t)
    df <- .ss_two_sample_df(X, g, variance = variance)
    .ss_t_to_z(t, df = df)
  }

  z_obs <- stat_from_group(group)
  stat_vol <- .ss_stat_template(prep, z_obs)

  null_fun <- function(b) {
    if (b < 1L || b > n_perm) stop("Permutation index out of range")
    stat_from_group(group[perm_index[[b]]])
  }

  list(
    stat_vol = stat_vol,
    z_vol = stat_vol,
    mask = prep$mask,
    null_fun = null_fun,
    prep = prep,
    group = group,
    perm_index = perm_index,
    variance = variance,
    statistic = if (isTRUE(z_equiv)) "z" else "t",
    design = "two_sample"
  )
}

#' Scale-Space One-Sample Wrapper
#'
#' Thin raw-image wrapper around \code{\link{scale_space_energy}} using
#' subject-level sign-flip resampling.
#'
#' @inheritParams make_null_fun_onesample
#' @param scales Scale-space specification.
#' @param statistics List of scale-space channel objects.
#' @param alpha FWER level.
#' @param studentize Whether to studentize scale-space fields.
#' @param ... Additional arguments passed to \code{\link{scale_space_energy}}.
#'
#' @export
scale_space_onesample <- function(subject_images,
                                  mask = NULL,
                                  scales = scale_space(c(6, 10, 14, 18, 24, 32, 48)),
                                  statistics = list(local_energy(), signed_coherence(tail = "two")),
                                  n_perm = 1000,
                                  seed = NULL,
                                  z_equiv = TRUE,
                                  template = NULL,
                                  alpha = 0.05,
                                  studentize = TRUE,
                                  ...) {
  nf <- make_null_fun_onesample(subject_images, mask = mask, n_perm = n_perm,
                                seed = seed, z_equiv = z_equiv, template = template)
  fit <- scale_space_energy(nf$stat_vol, mask = nf$mask, scales = scales,
                            statistics = statistics, null_fun = nf$null_fun,
                            n_perm = n_perm, studentize = studentize,
                            alpha = alpha, ...)
  fit$design <- nf
  fit
}

#' Scale-Space Two-Sample Wrapper
#'
#' Thin raw-image wrapper around \code{\link{scale_space_energy}} using
#' fixed-size two-group label permutations.
#'
#' @inheritParams make_null_fun_twosample
#' @param scales Scale-space specification.
#' @param statistics List of scale-space channel objects.
#' @param alpha FWER level.
#' @param studentize Whether to studentize scale-space fields.
#' @param ... Additional arguments passed to \code{\link{scale_space_energy}}.
#'
#' @export
scale_space_twosample <- function(subject_images,
                                  group,
                                  mask = NULL,
                                  scales = scale_space(c(6, 10, 14, 18, 24, 32, 48)),
                                  statistics = list(local_energy(), signed_coherence(tail = "two")),
                                  n_perm = 1000,
                                  seed = NULL,
                                  variance = c("welch", "pooled"),
                                  z_equiv = TRUE,
                                  template = NULL,
                                  alpha = 0.05,
                                  studentize = TRUE,
                                  ...) {
  variance <- match.arg(variance)
  nf <- make_null_fun_twosample(subject_images, group = group, mask = mask,
                                n_perm = n_perm, seed = seed,
                                variance = variance, z_equiv = z_equiv,
                                template = template)
  fit <- scale_space_energy(nf$stat_vol, mask = nf$mask, scales = scales,
                            statistics = statistics, null_fun = nf$null_fun,
                            n_perm = n_perm, studentize = studentize,
                            alpha = alpha, ...)
  fit$design <- nf
  fit
}

#' Construct Within-Subject Contrast Images
#'
#' Builds one planned repeated-measures contrast image per subject. The first
#' implementation supports a complete 2 x 2 interaction contrast.
#'
#' @param images List of 3D images, one per row of the design vectors.
#' @param subject Subject identifier for each image.
#' @param within Data frame or list containing within-subject factor columns.
#' @param contrast Contrast string. Currently supports "factor1:factor2" for
#'   two 2-level within-subject factors.
#'
#' @return List of contrast images, one per subject.
#'
#' @export
make_within_subject_contrast_images <- function(images,
                                                subject,
                                                within,
                                                contrast) {
  if (!is.list(images) || length(images) < 1L) stop("images must be a non-empty list")
  if (length(subject) != length(images)) stop("subject length must match images")
  within <- as.data.frame(within)
  if (nrow(within) != length(images)) stop("within rows must match images")
  parts <- strsplit(contrast, ":", fixed = TRUE)[[1L]]
  if (length(parts) != 2L) stop("only two-factor interaction contrasts are currently supported")
  if (!all(parts %in% names(within))) stop("contrast factors must be columns in within")

  f1 <- factor(within[[parts[1L]]], levels = unique(within[[parts[1L]]]))
  f2 <- factor(within[[parts[2L]]], levels = unique(within[[parts[2L]]]))
  if (nlevels(f1) != 2L || nlevels(f2) != 2L) {
    stop("interaction contrast currently requires two 2-level factors")
  }

  ids <- unique(subject)
  out <- vector("list", length(ids))
  names(out) <- as.character(ids)
  dims <- dim(images[[1L]])

  for (s in seq_along(ids)) {
    idx_s <- which(subject == ids[s])
    if (length(idx_s) != 4L) stop("each subject must have exactly four cells")
    cell <- list()
    for (a in levels(f1)) {
      for (b in levels(f2)) {
        idx <- idx_s[f1[idx_s] == a & f2[idx_s] == b]
        if (length(idx) != 1L) stop("each subject must have one image per cell")
        if (!identical(dim(images[[idx]]), dims)) stop("all images must share dimensions")
        cell[[paste(a, b, sep = "___")]] <- as.numeric(images[[idx]])
      }
    }
    l1 <- levels(f1)
    l2 <- levels(f2)
    v <- (cell[[paste(l1[2L], l2[2L], sep = "___")]] -
            cell[[paste(l1[2L], l2[1L], sep = "___")]]) -
      (cell[[paste(l1[1L], l2[2L], sep = "___")]] -
         cell[[paste(l1[1L], l2[1L], sep = "___")]])
    img <- images[[1L]]
    img[] <- v
    out[[s]] <- img
  }
  out
}
