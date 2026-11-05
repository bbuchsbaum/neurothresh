#' Coerce subject-level images to a canonical (subject by voxel) matrix
#'
#' Converts heterogeneous subject-level image inputs into a shared canonical
#' representation used by every subject-level null factory. Accepts a list of
#' 3D arrays or `NeuroVol`s, a 4D array, a `NeuroVec`, a `SparseNeuroVec`, or
#' a numeric `matrix` already in subject-by-voxel form.
#'
#' @param x Input images. One of: list of 3D arrays/`NeuroVol`s, 4D numeric
#'   array, `NeuroVec`, `SparseNeuroVec`, or a numeric `matrix` with one row
#'   per subject and one column per mask voxel.
#' @param mask Optional logical 3D mask. When `NULL`, the mask is inferred
#'   from the finite intersection across subjects (or, for `SparseNeuroVec`,
#'   from the intrinsic mask). When supplied, the mask is used exactly: any
#'   voxel inside it that is non-finite in the input is an error.
#' @param template Optional 3D template carrying spatial metadata. Required
#'   only for the `matrix` branch when `mask` does not already supply
#'   dimensions.
#' @param ... Additional arguments (currently unused; reserved for future
#'   per-branch options).
#'
#' @return An object of S3 class `"neurothresh_image_matrix"` with components:
#'   \describe{
#'     \item{X}{`n_subjects` x `n_mask` numeric matrix.}
#'     \item{mask}{3D logical array, or `LogicalNeuroVol` when input carried
#'       `NeuroSpace` metadata.}
#'     \item{mask_idx}{Integer vector of mask voxel positions in column-major
#'       order.}
#'     \item{template}{3D template volume (array or `NeuroVol`).}
#'     \item{dims}{Integer length-3 dimensions.}
#'     \item{coords}{`n_mask` x 3 numeric matrix of array indices.}
#'     \item{space}{`NeuroSpace` if available, else `NULL`.}
#'     \item{source}{Character: `"list"`, `"neurovec"`, `"sparseneurovec"`,
#'       `"array4d"`, or `"matrix"`.}
#'   }
#'
#' @export
as_subject_matrix <- function(x, mask = NULL, template = NULL, ...) {
  if (is.matrix(x)) {
    return(.asm_matrix(x, mask = mask, template = template))
  }
  if (inherits(x, "SparseNeuroVec")) {
    return(.asm_neurovec(x, mask = mask, sparse = TRUE))
  }
  if (inherits(x, "NeuroVec")) {
    return(.asm_neurovec(x, mask = mask, sparse = FALSE))
  }
  if (is.array(x) && length(dim(x)) == 4L) {
    return(.asm_array4d(x, mask = mask))
  }
  if (is.list(x)) {
    return(.asm_list(x, mask = mask))
  }
  stop(sprintf(
    "as_subject_matrix(): unsupported input class '%s'",
    paste(class(x), collapse = "/")
  ))
}

#' @keywords internal
.asm_finalize <- function(X, final_mask, dims, template, space, source) {
  mask_idx <- which(final_mask)
  if (length(mask_idx) == 0L) stop("Mask is empty")

  if (any(!is.finite(X))) {
    stop("Image data contain non-finite values inside the requested mask")
  }

  mask_arr <- array(as.logical(final_mask), dim = dims)
  mask_obj <- if (!is.null(space)) {
    tryCatch(
      neuroim2::LogicalNeuroVol(mask_arr, space),
      error = function(e) mask_arr
    )
  } else {
    mask_arr
  }

  grid <- which(mask_arr, arr.ind = TRUE)
  storage.mode(grid) <- "double"

  structure(
    list(
      X = X,
      mask = mask_obj,
      mask_idx = mask_idx,
      template = template,
      dims = as.integer(dims),
      coords = grid,
      space = space,
      source = source
    ),
    class = "neurothresh_image_matrix"
  )
}

#' @keywords internal
.asm_resolve_mask <- function(supplied, finite_intersection, dims, intrinsic = NULL) {
  if (is.null(supplied)) {
    if (!is.null(intrinsic)) {
      return(as.logical(intrinsic))
    }
    return(as.logical(finite_intersection))
  }
  user <- as.logical(supplied)
  if (length(user) != prod(dims)) {
    stop("mask length must match volume dimensions")
  }
  if (!is.null(intrinsic) && any(user & !as.logical(intrinsic))) {
    stop("Supplied mask includes voxels outside the SparseNeuroVec intrinsic mask")
  }
  if (any(user & !as.logical(finite_intersection))) {
    stop("Supplied mask includes non-finite voxels in input data")
  }
  user
}

#' @keywords internal
.asm_space <- function(x) {
  tryCatch(neuroim2::space(x), error = function(e) NULL)
}

#' @keywords internal
.asm_list <- function(images, mask = NULL) {
  if (length(images) < 2L) stop("images must be a list of at least two volumes")

  template <- images[[1L]]
  dims <- dim(template)
  if (is.null(dims) || length(dims) != 3L) {
    stop("images[[1]] must be 3D (array or NeuroVol)")
  }

  space <- .asm_space(template)

  vecs <- vector("list", length(images))
  ok <- rep(TRUE, prod(dims))
  for (i in seq_along(images)) {
    if (!identical(dim(images[[i]]), dims)) {
      stop("all images must have identical dimensions")
    }
    vecs[[i]] <- as.numeric(images[[i]])
    ok <- ok & is.finite(vecs[[i]])
  }

  final_mask <- .asm_resolve_mask(mask, ok, dims, intrinsic = NULL)

  mask_idx <- which(final_mask)
  if (length(mask_idx) == 0L) stop("Mask is empty")

  X <- matrix(NA_real_, nrow = length(images), ncol = length(mask_idx))
  for (i in seq_along(images)) {
    X[i, ] <- vecs[[i]][mask_idx]
  }

  .asm_finalize(X, final_mask, dims, template, space, source = "list")
}

#' @keywords internal
.asm_array4d <- function(arr, mask = NULL) {
  dims4 <- dim(arr)
  if (length(dims4) != 4L) stop("4D array expected")
  dims3 <- dims4[1:3]
  n_subj <- dims4[4]
  if (n_subj < 2L) stop("4D array must have at least two subjects on the 4th dim")

  X_full <- t(matrix(as.numeric(arr), nrow = prod(dims3), ncol = n_subj))
  finite_intersection <- apply(is.finite(X_full), 2, all)

  final_mask <- .asm_resolve_mask(mask, finite_intersection, dims3, intrinsic = NULL)
  mask_idx <- which(final_mask)
  if (length(mask_idx) == 0L) stop("Mask is empty")

  X <- X_full[, mask_idx, drop = FALSE]
  template <- array(NA_real_, dim = dims3)

  .asm_finalize(X, final_mask, dims3, template, space = NULL, source = "array4d")
}

#' @keywords internal
.asm_neurovec <- function(x, mask = NULL, sparse = FALSE) {
  dims4 <- dim(x)
  if (is.null(dims4) || length(dims4) != 4L) stop("NeuroVec must be 4D")
  dims3 <- as.integer(dims4[1:3])
  n_subj <- as.integer(dims4[4])
  if (n_subj < 2L) stop("NeuroVec must have at least two subjects on the 4th dim")

  X_full <- matrix(NA_real_, nrow = n_subj, ncol = prod(dims3))
  for (s in seq_len(n_subj)) {
    vol <- as.numeric(x[[s]])
    if (length(vol) != prod(dims3)) {
      stop(sprintf(
        "Subject %d extracted to length %d (expected %d)",
        s, length(vol), prod(dims3)
      ))
    }
    X_full[s, ] <- vol
  }
  finite_intersection <- apply(is.finite(X_full), 2, all)

  intrinsic <- NULL
  if (sparse) {
    intrinsic <- tryCatch(
      as.logical(as.array(neuroim2::mask(x))),
      error = function(e) NULL
    )
    if (is.null(intrinsic)) {
      intrinsic <- finite_intersection
    }
  }

  final_mask <- .asm_resolve_mask(mask, finite_intersection, dims3, intrinsic = intrinsic)
  mask_idx <- which(final_mask)
  if (length(mask_idx) == 0L) stop("Mask is empty")

  X <- X_full[, mask_idx, drop = FALSE]

  space <- .asm_space(x)
  template <- if (!is.null(space)) {
    tryCatch(
      neuroim2::DenseNeuroVol(array(NA_real_, dim = dims3), space),
      error = function(e) array(NA_real_, dim = dims3)
    )
  } else {
    array(NA_real_, dim = dims3)
  }

  .asm_finalize(X, final_mask, dims3, template, space,
                source = if (sparse) "sparseneurovec" else "neurovec")
}

#' @keywords internal
.asm_matrix <- function(images, mask = NULL, template = NULL) {
  if (is.null(template) && is.null(mask)) {
    stop("matrix inputs require template or mask")
  }

  if (is.null(template)) {
    dims <- dim(mask)
    if (is.null(dims) || length(dims) != 3L) {
      stop("mask must be a 3D logical array for matrix inputs")
    }
    template <- array(NA_real_, dim = dims)
    # No NeuroSpace metadata is available, so downstream voxel-spacing lookups
    # fall back to 1 mm isotropic. Warn the caller: any FWHM-in-mm setting
    # then behaves as FWHM-in-voxels.
    warning(
      "as_subject_matrix(): matrix input supplied without a template; ",
      "voxel spacing is unknown and defaults to 1 mm isotropic. Pass a ",
      "NeuroVol/NeuroSpace template to preserve mm-based bandwidths."
    )
  } else {
    dims <- dim(template)
    if (is.null(dims) || length(dims) != 3L) stop("template must be 3D")
  }

  if (is.null(mask)) {
    mask_arr <- array(is.finite(as.numeric(template)), dim = dims)
  } else {
    mask_arr <- array(as.logical(mask), dim = dims)
  }

  mask_idx <- which(mask_arr)
  if (length(mask_idx) == 0L) stop("Mask is empty")
  if (ncol(images) != length(mask_idx)) {
    stop("ncol(images) must equal the number of mask voxels")
  }

  X <- matrix(as.numeric(images), nrow = nrow(images), ncol = ncol(images))
  if (any(!is.finite(X))) {
    stop("Image data contain non-finite values inside the requested mask")
  }

  space <- .asm_space(template)

  .asm_finalize(X, mask_arr, dims, template, space, source = "matrix")
}

#' @keywords internal
.ss_vec_to_vol <- function(z_vec, prep) {
  full <- rep(NA_real_, prod(prep$dims))
  full[prep$mask_idx] <- z_vec
  .array_to_vol(full, prep$dims, prep$template)
}

#' @keywords internal
.ss_stat_template <- function(prep, z_vec) {
  .ss_vec_to_vol(z_vec, prep)
}

#' @keywords internal
.ss_t_to_z <- function(t, df, two_sided = FALSE) {
  # A signed Z-equivalent is the same for either tail convention. Work with
  # the absolute upper log-tail to avoid cancellation at large negative t
  # and underflow at large positive t; t = 0 maps exactly to Z = 0.
  log_p <- stats::pt(abs(t), df = df, lower.tail = FALSE, log.p = TRUE)
  sign(t) * stats::qnorm(log_p, lower.tail = FALSE, log.p = TRUE)
}

#' @keywords internal
.ss_one_sample_t <- function(X) {
  n <- nrow(X)
  if (n < 2L) stop("one-sample statistic requires at least two subjects")
  mu <- colMeans(X)
  s <- apply(X, 2L, stats::sd)
  se <- s / sqrt(n)
  t <- mu / se
  t[!is.finite(t)] <- 0
  t
}

#' @keywords internal
.ss_two_sample_t <- function(X, group, variance = c("welch", "pooled")) {
  variance <- match.arg(variance)
  group <- as.factor(group)
  if (nlevels(group) != 2L) stop("group must have exactly two levels")
  if (length(group) != nrow(X)) stop("length(group) must equal nrow(X)")

  lev <- levels(group)
  X1 <- X[group == lev[1L], , drop = FALSE]
  X2 <- X[group == lev[2L], , drop = FALSE]
  n1 <- nrow(X1)
  n2 <- nrow(X2)
  if (n1 < 2L || n2 < 2L) stop("each group must contain at least two subjects")

  m1 <- colMeans(X1)
  m2 <- colMeans(X2)
  v1 <- apply(X1, 2L, stats::var)
  v2 <- apply(X2, 2L, stats::var)

  if (variance == "pooled") {
    sp2 <- ((n1 - 1) * v1 + (n2 - 1) * v2) / (n1 + n2 - 2)
    se <- sqrt(sp2 * (1 / n1 + 1 / n2))
  } else {
    se <- sqrt(v1 / n1 + v2 / n2)
  }

  t <- (m2 - m1) / se
  t[!is.finite(t)] <- 0
  t
}

#' @keywords internal
.ss_two_sample_df <- function(X, group, variance = c("welch", "pooled")) {
  variance <- match.arg(variance)
  group <- as.factor(group)
  lev <- levels(group)
  X1 <- X[group == lev[1L], , drop = FALSE]
  X2 <- X[group == lev[2L], , drop = FALSE]
  n1 <- nrow(X1)
  n2 <- nrow(X2)
  if (variance == "pooled") return(n1 + n2 - 2)
  v1 <- apply(X1, 2L, stats::var)
  v2 <- apply(X2, 2L, stats::var)
  a <- v1 / n1
  b <- v2 / n2
  df <- (a + b)^2 / (a^2 / (n1 - 1) + b^2 / (n2 - 1))
  df[!is.finite(df)] <- n1 + n2 - 2
  df
}
