#' Canonicalize statistical maps to Z-scores
#'
#' Converts various statistical map types (Z-scores, t-statistics, or
#' -log10(p) values) to a common Z-score representation for unified
#' thresholding procedures.
#'
#' @param stat_vol A NeuroVol object containing the statistical map
#' @param stat_type Character string specifying the input type:
#'   \describe{
#'     \item{"Z"}{Z-scores (identity transform)}
#'     \item{"t"}{t-statistics (requires \code{df})}
#'     \item{"neglog10p"}{Negative log10 p-values}
#'   }
#' @param df Degrees of freedom, required when \code{stat_type = "t"}
#' @param tail Character string specifying the tail:
#'   \describe{
#'     \item{"pos"}{Positive tail (activations)}
#'     \item{"neg"}{Negative tail (deactivations)}
#'     \item{"two"}{Two-sided test}
#'   }
#' @param p_side Character string for -log10(p) inputs specifying whether
#'   the p-values are one-sided or two-sided: "one" or "two"
#' @param sign_vol Optional NeuroVol providing sign information for -log10(p)
#'   inputs to recover directionality
#'
#' @return A list with components:
#'   \describe{
#'     \item{Zeq}{NeuroVol of equivalent Z-scores}
#'     \item{p_one}{NeuroVol of one-sided p-values}
#'     \item{p_two}{NeuroVol of two-sided p-values}
#'     \item{mask}{Logical mask of valid voxels}
#'   }
#'
#' @details
#' The conversions are:
#' \itemize{
#'   \item Z-scores: Used directly
#'   \item t-statistics: Converted via \code{qnorm(pt(t, df))}
#'   \item -log10(p): Converted via \code{qnorm(1 - 10^(-x))} with optional
#'     sign recovery from \code{sign_vol}
#' }
#'
#' @examples
#' \dontrun{
#' # From t-statistics with 50 degrees of freedom
#' result <- canonicalize_stat(t_map, stat_type = "t", df = 50)
#'
#' # From -log10(p) with sign information
#' result <- canonicalize_stat(logp_map, stat_type = "neglog10p",
#'                             sign_vol = t_map)
#' }
#'
#' @export
canonicalize_stat <- function(stat_vol,
                              stat_type = c("Z", "t", "neglog10p"),
                              df = NULL,
                              tail = c("pos", "neg", "two"),
                              p_side = c("one", "two"),
                              sign_vol = NULL) {

  stat_type <- match.arg(stat_type)
  tail <- match.arg(tail)
  p_side <- match.arg(p_side)

  # Get mask of valid (finite) values
  mask <- is.finite(as.numeric(stat_vol))
  s <- as.numeric(stat_vol)[mask]

  n <- length(s)
  Zeq <- rep(NA_real_, n)
  p_one <- rep(NA_real_, n)
  p_two <- rep(NA_real_, n)

  if (stat_type == "Z") {
    z <- s
    ppos <- 1 - pnorm(z)
    pneg <- pnorm(z)

    if (tail == "pos") {
      p_one <- .clamp01(ppos)
      Zeq <- z
    } else if (tail == "neg") {
      p_one <- .clamp01(pneg)
      Zeq <- z
    } else {
      p_two <- .clamp01(2 * (1 - pnorm(abs(z))))
      Zeq <- z
    }

    if (tail != "two") {
      p_two <- .clamp01(2 * pmin(ppos, pneg))
    }

  } else if (stat_type == "t") {
    if (is.null(df)) stop("df is required for stat_type='t'")

    tval <- s
    p_two <- .clamp01(2 * (1 - pt(abs(tval), df = df)))

    ppos <- 1 - pt(tval, df = df)
    pneg <- pt(tval, df = df)
    p_one <- if (tail == "neg") .clamp01(pneg) else .clamp01(ppos)

    # Convert via probability integral transform
    u <- .clamp01(pt(tval, df = df))
    Zeq <- qnorm(u)

  } else if (stat_type == "neglog10p") {
    lp <- s
    p_raw <- .log10p_to_p(lp)

    if (p_side == "one") {
      p_one <- p_raw
      p_two <- .clamp01(2 * p_raw)
      zabs <- qnorm(1 - p_raw)
    } else {
      p_two <- p_raw
      p_one <- .clamp01(p_raw / 2)
      zabs <- qnorm(1 - p_raw / 2)
    }

    # Recover sign if available
    if (!is.null(sign_vol)) {
      sg <- sign(as.numeric(sign_vol)[mask])
      sg[sg == 0] <- 1
      Zeq <- sg * zabs
    } else {
      if (tail == "pos") {
        Zeq <- zabs
      } else if (tail == "neg") {
        Zeq <- -zabs
      } else {
        Zeq <- zabs
      }
    }
  }

  # Create output volumes with same space as input
  Zeq_vol <- stat_vol
  Zeq_vol[mask] <- Zeq

  pone_vol <- stat_vol
  pone_vol[mask] <- p_one

  ptwo_vol <- stat_vol
  ptwo_vol[mask] <- p_two

  list(Zeq = Zeq_vol, p_one = pone_vol, p_two = ptwo_vol, mask = mask)
}

#' @keywords internal
.clamp01 <- function(u) {
  pmin(pmax(u, .Machine$double.xmin), 1 - .Machine$double.eps)
}

#' @keywords internal
.log10p_to_p <- function(log10p) {
  # p = 10^(-log10p) = exp(-log(10) * log10p)
  p <- exp(-log(10) * log10p)
  .clamp01(p)
}
