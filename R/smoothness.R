#' FWHM-based smoothness specification
#'
#' Tagged container for full-width-at-half-maximum (FWHM) smoothness in
#' millimetres. RFT methods accept only this form; Monte-Carlo null factories
#' may also accept it (alongside [smoothness_acf()]).
#'
#' @param fwhm_mm Numeric; either a length-1 isotropic FWHM or a length-3
#'   per-axis FWHM, in millimetres. Must be strictly positive.
#'
#' @return Object of class `c("neurothresh_smoothness_fwhm",
#'   "neurothresh_smoothness")`.
#' @export
#' @examples
#' s1 <- smoothness_fwhm(6)
#' s3 <- smoothness_fwhm(c(6, 6, 8))
smoothness_fwhm <- function(fwhm_mm) {
  if (!is.numeric(fwhm_mm)) stop("fwhm_mm must be numeric")
  if (length(fwhm_mm) == 1L) fwhm_mm <- rep(fwhm_mm, 3L)
  if (length(fwhm_mm) != 3L) {
    stop("fwhm_mm must have length 1 or 3")
  }
  if (any(!is.finite(fwhm_mm)) || any(fwhm_mm <= 0)) {
    stop("fwhm_mm must be strictly positive and finite")
  }
  structure(
    list(fwhm_mm = as.numeric(fwhm_mm)),
    class = c("neurothresh_smoothness_fwhm", "neurothresh_smoothness")
  )
}

#' AFNI-style ACF smoothness specification
#'
#' Tagged container for the three-parameter mixed Gaussian + monoexponential
#' autocorrelation model used by AFNI's `3dFWHMx -acf` and downstream
#' Monte-Carlo null generators.
#'
#' @param a Numeric scalar in (0, 1); mixing weight of the Gaussian component.
#' @param b Numeric scalar > 0; Gaussian width parameter (mm).
#' @param c Numeric scalar > 0; exponential length scale (mm).
#'
#' @return Object of class `c("neurothresh_smoothness_acf",
#'   "neurothresh_smoothness")`.
#' @export
#' @examples
#' s <- smoothness_acf(a = 0.7, b = 3.0, c = 8.0)
smoothness_acf <- function(a, b, c) {
  if (!is.numeric(a) || length(a) != 1L) stop("a must be a length-1 numeric")
  if (!is.numeric(b) || length(b) != 1L) stop("b must be a length-1 numeric")
  if (!is.numeric(c) || length(c) != 1L) stop("c must be a length-1 numeric")
  if (!is.finite(a) || a <= 0 || a >= 1) {
    stop("a must lie strictly between 0 and 1")
  }
  if (!is.finite(b) || b <= 0) stop("b must be strictly positive and finite")
  if (!is.finite(c) || c <= 0) stop("c must be strictly positive and finite")
  structure(
    list(a = as.numeric(a), b = as.numeric(b), c = as.numeric(c)),
    class = c("neurothresh_smoothness_acf", "neurothresh_smoothness")
  )
}

#' Print a smoothness specification
#'
#' @param x A `neurothresh_smoothness`.
#' @param ... Unused.
#' @export
print.neurothresh_smoothness <- function(x, ...) {
  if (inherits(x, "neurothresh_smoothness_fwhm")) {
    cat("<smoothness_fwhm> mm:",
        paste(format(x$fwhm_mm, digits = 4), collapse = " x "), "\n")
  } else if (inherits(x, "neurothresh_smoothness_acf")) {
    cat(sprintf("<smoothness_acf> a = %.4f, b = %.4f, c = %.4f\n",
                x$a, x$b, x$c))
  } else {
    cat("<neurothresh_smoothness>\n")
    utils::str(unclass(x))
  }
  invisible(x)
}

#' Test whether an object is a smoothness specification
#'
#' @param x Any R object.
#' @return Logical scalar.
#' @export
is_smoothness <- function(x) inherits(x, "neurothresh_smoothness")
