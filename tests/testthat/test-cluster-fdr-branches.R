# Branch/semantic checks for cluster_fdr

make_two_blob_map <- function(dims = c(18, 18, 18), amp = 7) {
  z <- array(rnorm(prod(dims), sd = 0.5), dim = dims)
  z[4:7, 4:7, 4:7] <- z[4:7, 4:7, 4:7] + amp
  z[12:15, 12:15, 12:15] <- z[12:15, 12:15, 12:15] - amp
  z
}

test_that("cluster_fdr requires fwhm_mm for RFT p-values", {
  z <- array(rnorm(10 * 10 * 10), dim = c(10, 10, 10))
  expect_error(
    cluster_fdr(z, p_method = "rft", fwhm_mm = NULL),
    "requires smoothness"
  )
})

test_that("cluster_fdr permutation branch runs and returns valid p-values", {
  set.seed(3)
  z <- make_two_blob_map()
  res <- cluster_fdr(
    z,
    cluster_thresh = 2.8,
    q = 0.2,
    p_method = "perm",
    n_perm = 40,
    allow_voxel_signflip = TRUE,
    tail = "pos",
    seed = 42
  )

  expect_true(is.data.frame(res$table))
  expect_equal(res$n_clusters, sum(res$table$significant))
  expect_true(all(res$table$p_unc >= 0 & res$table$p_unc <= 1))
})

test_that("split_q two-sided policy matches one-sided q/2 when only one tail has clusters", {
  dims <- c(16, 16, 16)
  z <- array(rnorm(prod(dims), sd = 0.2), dim = dims)
  z[5:8, 5:8, 5:8] <- z[5:8, 5:8, 5:8] + 7

  q <- 0.2
  two_split <- cluster_fdr(
    z,
    cluster_thresh = 2.5,
    q = q,
    p_method = "rft",
    fwhm_mm = 4,
    tail = "two",
    two_sided_policy = "split_q"
  )
  one_pos <- cluster_fdr(
    z,
    cluster_thresh = 2.5,
    q = q / 2,
    p_method = "rft",
    fwhm_mm = 4,
    tail = "pos"
  )

  expect_identical(two_split$sig_pos, one_pos$sig_mask)
})

test_that("BH_all two-sided mask is union of per-tail BH_all decisions", {
  set.seed(9)
  z <- make_two_blob_map(amp = 8)
  q <- 0.15

  two <- cluster_fdr(
    z,
    cluster_thresh = 2.8,
    q = q,
    p_method = "perm",
    n_perm = 30,
    allow_voxel_signflip = TRUE,
    tail = "two",
    two_sided_policy = "BH_all",
    seed = 11
  )

  expect_identical(two$sig_mask, two$sig_pos | two$sig_neg)
})
