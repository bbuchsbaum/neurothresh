#' GLM Permutation Specification
#'
#' @param method Currently only "freedman_lane".
#' @param n_perm Number of permutations.
#' @param seed Optional random seed.
#'
#' @export
permutation_glm <- function(method = c("freedman_lane"), n_perm = 1000, seed = NULL) {
  method <- match.arg(method)
  if (!is.numeric(n_perm) || length(n_perm) != 1L || n_perm < 1L) stop("n_perm must be >= 1")
  structure(list(method = method, n_perm = as.integer(n_perm), seed = seed),
            class = "neurothresh_permutation_glm")
}

#' GLM Freedman-Lane Null Generator
#'
#' Builds an observed GLM contrast statistic map and deterministic
#' Freedman-Lane residual-permutation null draws.
#'
#' @param images List of subject-level contrast images or matrix.
#' @param design Formula, model matrix, or data frame. If a formula is supplied,
#'   \code{data} must contain the variables.
#' @param contrast Named coefficient or numeric contrast vector.
#' @param mask Optional logical mask.
#' @param data Optional data for formula designs.
#' @param permutation A \code{\link{permutation_glm}} object.
#' @param statistic "t" or "z".
#' @param template Optional template for matrix inputs.
#'
#' @export
make_null_fun_glm <- function(images,
                              design,
                              contrast,
                              mask = NULL,
                              data = NULL,
                              permutation = permutation_glm(),
                              statistic = c("t", "z"),
                              template = NULL) {
  statistic <- match.arg(statistic)
  if (!inherits(permutation, "neurothresh_permutation_glm")) {
    stop("permutation must be created by permutation_glm()")
  }
  prep <- as_subject_matrix(images, mask = mask, template = template)
  Y <- prep$X
  X <- .ss_design_matrix(design, data = data)
  if (nrow(X) != nrow(Y)) stop("design rows must match number of images")
  cvec <- .ss_contrast_vector(contrast, X)

  # Restrict the reduced model to the single-coefficient case. A multi-column
  # contrast (e.g. c(0, 1, -1, 0)) requires partialling the tested subspace
  # out of X, not just dropping zero-weight columns. Until that QR-based
  # construction lands, refuse vector contrasts so users do not silently get
  # a reduced model that still spans the contrast direction.
  nz <- which(abs(cvec) > .Machine$double.eps)
  if (length(nz) != 1L) {
    stop("make_null_fun_glm() currently supports single-coefficient ",
         "contrasts only; supply a coefficient name or a contrast vector ",
         "with exactly one non-zero entry. Multi-column contrasts will be ",
         "added once a QR-based partial-out is implemented.")
  }
  keep <- seq_len(ncol(X)) != nz
  X0 <- X[, keep, drop = FALSE]
  if (ncol(X0) == 0L) X0 <- matrix(1, nrow = nrow(X), ncol = 1L)

  if (!is.null(permutation$seed)) set.seed(permutation$seed)
  perm_index <- replicate(permutation$n_perm, sample.int(nrow(Y)), simplify = FALSE)

  fit0 <- .ss_lm_fit_matrix(X0, Y)
  fitted0 <- X0 %*% fit0$beta
  resid0 <- Y - fitted0

  stat_from_Y <- function(Yi) {
    st <- .ss_glm_contrast_stat(X, Yi, cvec)
    if (statistic == "z") .ss_t_to_z(st$t, df = st$df) else st$t
  }

  z_obs <- stat_from_Y(Y)
  stat_vol <- .ss_stat_template(prep, z_obs)

  null_fun <- function(b) {
    if (b < 1L || b > permutation$n_perm) stop("Permutation index out of range")
    Yb <- fitted0 + resid0[perm_index[[b]], , drop = FALSE]
    stat_from_Y(Yb)
  }

  list(
    stat_vol = stat_vol,
    z_vol = stat_vol,
    mask = prep$mask,
    null_fun = null_fun,
    prep = prep,
    design = X,
    reduced_design = X0,
    contrast = cvec,
    perm_index = perm_index,
    statistic = statistic,
    method = "freedman_lane"
  )
}

#' @keywords internal
.ss_design_matrix <- function(design, data = NULL) {
  if (inherits(design, "formula")) {
    if (is.null(data)) stop("data is required when design is a formula")
    return(stats::model.matrix(design, data = data))
  }
  if (is.data.frame(design)) {
    return(stats::model.matrix(~ ., data = design))
  }
  X <- as.matrix(design)
  if (!is.numeric(X)) stop("design matrix must be numeric")
  X
}

#' @keywords internal
.ss_contrast_vector <- function(contrast, X) {
  p <- ncol(X)
  if (is.character(contrast)) {
    if (length(contrast) != 1L) stop("character contrast must name one coefficient")
    nm <- colnames(X)
    if (is.null(nm) || !(contrast %in% nm)) stop("contrast coefficient not found in design")
    cvec <- numeric(p)
    cvec[match(contrast, nm)] <- 1
    return(cvec)
  }
  cvec <- as.numeric(contrast)
  if (length(cvec) != p) stop("numeric contrast length must equal ncol(design)")
  cvec
}

#' @keywords internal
.ss_lm_fit_matrix <- function(X, Y) {
  qrX <- qr(X)
  beta <- qr.coef(qrX, Y)
  beta[!is.finite(beta)] <- 0
  list(beta = beta, qr = qrX, df = nrow(X) - qrX$rank)
}

#' @keywords internal
.ss_glm_contrast_stat <- function(X, Y, cvec) {
  fit <- .ss_lm_fit_matrix(X, Y)
  beta <- fit$beta
  fitted <- X %*% beta
  resid <- Y - fitted
  df <- fit$df
  rss <- colSums(resid * resid)
  sigma2 <- rss / max(df, 1L)
  R <- qr.R(fit$qr)
  XtX_inv <- chol2inv(R)
  cvar <- as.numeric(t(cvec) %*% XtX_inv %*% cvec)
  est <- as.numeric(t(cvec) %*% beta)
  se <- sqrt(pmax(sigma2 * cvar, 0))
  t <- est / se
  t[!is.finite(t)] <- 0
  list(t = t, df = df)
}
