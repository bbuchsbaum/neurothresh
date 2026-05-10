#' Null specification S3 class
#'
#' A `neurothresh_null_spec` is the uniform handle every permutation-based
#' inference method consumes. It bundles the observed (canonicalized)
#' statistic map, the mask used to build the null, a deterministic-by-index
#' draw function, and metadata describing the resampling design.
#'
#' Factories that return `neurothresh_null_spec` objects:
#' [null_voxel_signflip()], [null_subject_signflip()], [null_one_sample()],
#' [null_two_sample()], [null_glm()], [null_mc_fwhm()], [null_mc_acf()],
#' [null_parcel_perm()].
#'
#' @section Slots:
#' \describe{
#'   \item{stat_vol}{Observed Z-equivalent statistic map (3D array or
#'     `NeuroVol`).}
#'   \item{mask}{Logical mask (3D array or `LogicalNeuroVol`).}
#'   \item{draw}{`function(b)` returning a length-`sum(mask)` numeric vector
#'     for permutation index `1 <= b <= n_perm`. Deterministic by index.}
#'   \item{n_perm}{Maximum supported draws (integer).}
#'   \item{df}{Degrees of freedom (numeric, may be `NA`).}
#'   \item{kind}{One of `"voxel_signflip"`, `"subject_signflip"`,
#'     `"one_sample"`, `"two_sample"`, `"glm"`, `"mc_fwhm"`, `"mc_acf"`,
#'     `"parcel_perm"`.}
#'   \item{meta}{Free-form named list of method-specific metadata
#'     (smoothness, design matrix, perm indices, etc.).}
#' }
#'
#' @name null_spec
NULL

#' Construct a neurothresh null specification
#'
#' Internal constructor. Factories should call this rather than building
#' the structure by hand so validation runs uniformly.
#'
#' @param stat_vol Observed (canonicalized) statistic map.
#' @param mask Logical mask.
#' @param draw `function(b)` returning length-`sum(mask)` numeric.
#' @param n_perm Integer maximum number of draws.
#' @param df Degrees of freedom (numeric scalar, may be `NA`).
#' @param kind Character tag identifying the resampling design.
#' @param meta Optional named list of extras.
#'
#' @return Object of class `"neurothresh_null_spec"`.
#' @keywords internal
new_null_spec <- function(stat_vol, mask, draw, n_perm,
                          df = NA_real_, kind = "unknown", meta = list()) {
  if (!is.function(draw)) stop("draw must be a function")
  if (!is.numeric(n_perm) || length(n_perm) != 1L || n_perm < 1L) {
    stop("n_perm must be a positive integer scalar")
  }
  if (!is.character(kind) || length(kind) != 1L) {
    stop("kind must be a length-1 character")
  }
  if (!is.list(meta)) stop("meta must be a list")

  structure(
    list(
      stat_vol = stat_vol,
      mask = mask,
      draw = draw,
      n_perm = as.integer(n_perm),
      df = as.numeric(df),
      kind = kind,
      meta = meta
    ),
    class = "neurothresh_null_spec"
  )
}

#' Test whether an object is a neurothresh_null_spec
#'
#' @param x Any R object.
#' @return Logical scalar.
#' @export
is_null_spec <- function(x) inherits(x, "neurothresh_null_spec")

#' Print a neurothresh null specification
#'
#' @param x A `neurothresh_null_spec`.
#' @param ... Unused.
#' @export
print.neurothresh_null_spec <- function(x, ...) {
  cat("<neurothresh_null_spec>\n")
  cat("  kind   :", x$kind, "\n")
  cat("  n_perm :", x$n_perm, "\n")
  if (!is.na(x$df)) cat("  df     :", format(x$df, digits = 4), "\n")
  n_mask <- if (is.null(x$mask)) NA_integer_ else sum(as.logical(x$mask), na.rm = TRUE)
  if (!is.na(n_mask)) cat("  n_mask :", n_mask, "\n")
  if (length(x$meta) > 0L) {
    cat("  meta   :", paste(names(x$meta), collapse = ", "), "\n")
  }
  invisible(x)
}

#' Wrap a bare draw function in a null_spec
#'
#' Convenience for ad-hoc resampling models that already produce a
#' deterministic-by-index draw function. Most users should reach for the
#' high-level factories (e.g. [null_one_sample()]) instead.
#'
#' @param draw `function(b)` returning a length-`sum(mask)` numeric.
#' @param stat_vol Observed statistic map.
#' @param n_perm Number of permutations the draw function supports.
#' @param mask Optional logical 3D mask. Defaults to finite voxels in
#'   `stat_vol`.
#' @param kind Optional kind tag (default `"custom"`).
#' @return A [`null_spec`].
#' @export
null_from_function <- function(draw, stat_vol, n_perm,
                               mask = NULL, kind = "custom") {
  if (is.null(mask)) mask <- is.finite(as.numeric(stat_vol))
  new_null_spec(stat_vol = stat_vol, mask = mask, draw = draw,
                n_perm = n_perm, kind = kind, meta = list())
}

#' Validate a neurothresh_null_spec
#'
#' Checks that the object's slots are internally consistent: `draw(1)`
#' returns a finite numeric vector of length `sum(mask)`. Useful for
#' factories or tests; production code should not need to call this.
#'
#' @param ns A `neurothresh_null_spec`.
#' @return The object invisibly. Throws on failure.
#' @export
validate_null_spec <- function(ns) {
  if (!is_null_spec(ns)) stop("not a neurothresh_null_spec")
  if (!is.function(ns$draw)) stop("draw is not a function")
  n_mask <- sum(as.logical(ns$mask), na.rm = TRUE)
  if (n_mask == 0L) stop("mask is empty")
  test <- tryCatch(ns$draw(1L), error = function(e) e)
  if (inherits(test, "error")) {
    stop("draw(1) errored: ", conditionMessage(test))
  }
  if (!is.numeric(test) || length(test) != n_mask) {
    stop(sprintf(
      "draw(1) returned length %d (expected %d, the mask size)",
      length(test), n_mask
    ))
  }
  invisible(ns)
}
