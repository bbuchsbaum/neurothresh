#' @keywords internal
"_PACKAGE"

#' neurothresh: Neuroimaging Statistical Thresholding with LR-MFT
#'
#' Statistical thresholding for neuroimaging data using Likelihood-Ratio
#' Matched-Filter Thresholding (LR-MFT) with hierarchical stepdown inference.
#'
#' @section Main Functions:
#' \describe{
#'   \item{\code{\link{hier_scan}}}{Main entry point for hierarchical LR-MFT analysis}
#'   \item{\code{\link{score_set}}}{Compute prior-weighted statistic T_kappa(R)}
#'   \item{\code{\link{canonicalize_stat}}}{Convert Z/t/neglog10p to Z-scores}
#' }
#'
#' @section Baseline Methods:
#' \describe{
#'   \item{\code{\link{rft_peak_fwer}}}{RFT peak-level FWER}
#'   \item{\code{\link{rft_cluster_fwer}}}{RFT cluster-level FWER}
#'   \item{\code{\link{tfce_transform}}}{TFCE transformation}
#'   \item{\code{\link{tfce_fwer}}}{TFCE with permutation FWER}
#'   \item{\code{\link{cluster_fdr}}}{Topological cluster-FDR}
#' }
#'
#' @section Key Concepts:
#' LR-MFT is provably more powerful than traditional GRF/RFT methods for
#' detecting spatially extended activation blobs. It uses a prior-weighted
#' soft evidence statistic:
#'
#' \deqn{T_\kappa(R) = \log \sum_{v \in R} \pi(v) \cdot \exp(\kappa \cdot Z(v))}
#'
#' where \eqn{\pi(v)} is the prior weight and \eqn{\kappa} is the temperature.
#'
#' @useDynLib neurothresh, .registration=TRUE
#' @importFrom Rcpp evalCpp
#' @importFrom methods new
#' @importFrom stats qnorm pt pnorm qt quantile uniroot
#' @name neurothresh-package
NULL
