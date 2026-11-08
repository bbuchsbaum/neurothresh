ns <- asNamespace("neurothresh")

smooth_field <- function(dims, sd = 1, scale = 3) {
  z <- rnorm(prod(dims))
  z <- ns$gaussian_smooth3d_cpp(z, as.integer(dims), rep(sd, 3), 3L)
  array(z / stats::sd(z) * scale / 3, dims)
}

test_that("compiled TFCE equals the R reference on random masked volumes", {
  set.seed(101)
  for (rep in 1:25) {
    dims <- sample(4:9, 3, TRUE)
    z <- smooth_field(dims)
    mask <- array(runif(prod(dims)) > 0.25, dims)
    conn <- c(6L, 18L, 26L)[(rep %% 3) + 1L]
    H <- c(2, 1.5)[(rep %% 2) + 1L]
    E <- c(0.5, 1)[(rep %% 2) + 1L]
    dh <- c(0.1, 0.05, 0.37)[(rep %% 3) + 1L]
    for (tail in c("pos", "neg")) {
      a <- ns$.tfce_from_vec(as.numeric(z), dims, mask, H, E, dh, tail, conn)
      b <- ns$.tfce_from_vec_r(as.numeric(z), dims, mask, H, E, dh, tail, conn)
      expect_equal(a, b, tolerance = 1e-12)
    }
  }
})

test_that("tfce_transform two-sided is the voxelwise max of the tails", {
  set.seed(102)
  dims <- c(7, 6, 5)
  z <- smooth_field(dims)
  mask <- array(TRUE, dims)
  pos <- ns$.tfce_from_vec_r(as.numeric(z), dims, mask, 2, 0.5, 0.1, "pos")
  neg <- ns$.tfce_from_vec_r(as.numeric(z), dims, mask, 2, 0.5, 0.1, "neg")
  two <- tfce_transform(z, mask, tail = "two")
  expect_equal(as.numeric(two), pmax(pos, neg), tolerance = 1e-12)
})

test_that("TFCE edge cases", {
  dims <- c(5, 5, 5)
  mask <- array(TRUE, dims)
  # all non-positive
  z <- array(-abs(rnorm(125)), dims)
  expect_identical(ns$.tfce_from_vec(as.numeric(z), dims, mask, 2, 0.5, 0.1, "pos"),
                   rep(0, 125))
  # maximum below dh: only the h = 0 step, so zero everywhere
  z <- array(0.05, dims)
  expect_identical(ns$.tfce_from_vec(as.numeric(z), dims, mask, 2, 0.5, 0.1, "pos"),
                   rep(0, 125))
  # single supra-threshold voxel: sum_{i=1}^{n} 1 * (i dh)^2 dh
  z <- array(0, dims)
  z[3, 3, 3] <- 1
  out <- ns$.tfce_from_vec(as.numeric(z), dims, mask, 2, 0.5, 0.1, "pos")
  hs <- seq(0, 1, by = 0.1)[-1]
  expect_equal(out[which(z > 0)], sum(hs^2 * 0.1), tolerance = 1e-14)
  expect_true(all(out[z == 0] == 0))
  # masked-out voxels are zero and do not bridge clusters
  z <- array(0, dims)
  z[1:5, 3, 3] <- 2
  m2 <- mask
  m2[3, 3, 3] <- FALSE
  a <- ns$.tfce_from_vec(as.numeric(z), dims, m2, 2, 0.5, 0.1, "pos")
  b <- ns$.tfce_from_vec_r(as.numeric(z), dims, m2, 2, 0.5, 0.1, "pos")
  expect_equal(a, b, tolerance = 1e-14)
  expect_equal(a[which(!m2)], 0)
  # non-finite values inside the mask never exceed a height
  z <- smooth_field(dims)
  z[2, 2, 2] <- NA
  z[4, 4, 4] <- Inf
  zc <- z
  zc[!is.finite(zc)] <- 0
  expect_equal(ns$.tfce_from_vec(as.numeric(z), dims, mask, 2, 0.5, 0.1, "pos"),
               ns$.tfce_from_vec_r(as.numeric(zc), dims, mask, 2, 0.5, 0.1, "pos"),
               tolerance = 1e-12)
  # empty mask
  expect_identical(ns$.tfce_from_vec(as.numeric(z), dims, array(FALSE, dims),
                                     2, 0.5, 0.1, "pos"), rep(0, 125))
})

test_that("connectivity changes TFCE as expected and is validated", {
  dims <- c(3, 3, 3)
  z <- array(0, dims)
  z[1, 1, 1] <- 2
  z[2, 2, 2] <- 2  # corner neighbour only
  mask <- array(TRUE, dims)
  t6 <- tfce_transform(z, mask, connectivity = 6)
  t26 <- tfce_transform(z, mask, connectivity = 26)
  expect_true(all(t26[z > 0] > t6[z > 0]))
  expect_error(tfce_transform(z, mask, connectivity = 8), "6, 18 or 26")
})

test_that("compiled cluster labelling equals the R reference", {
  set.seed(103)
  for (rep in 1:20) {
    dims <- sample(4:9, 3, TRUE)
    z <- smooth_field(dims)
    supra <- as.numeric(z) > 0.3 & runif(prod(dims)) > 0.1
    conn <- c(6L, 18L, 26L)[(rep %% 3) + 1L]
    a <- ns$.find_clusters(supra, dims, conn)
    b <- ns$.find_clusters_r(supra, dims, conn)
    expect_identical(a$labels, as.integer(b$labels))
    expect_identical(a$sizes, b$sizes)
  }
  empty <- ns$.find_clusters(rep(FALSE, 27), c(3, 3, 3))
  expect_identical(empty$sizes, integer(0))
})

test_that("label_clusters reports extent, mass and peak", {
  z <- array(0, c(6, 6, 6))
  z[2:3, 2:3, 2:3] <- 4
  z[3, 3, 3] <- 6
  z[5, 5, 5] <- -5
  res <- label_clusters(z, threshold = 3)
  expect_equal(nrow(res$table), 1L)
  expect_equal(res$table$size, 8L)
  expect_equal(res$table$mass, 7 * 1 + 3)
  expect_equal(res$table$peak_value, 6)
  expect_equal(unname(unlist(res$table[1, c("peak_i", "peak_j", "peak_k")])), c(3, 3, 3))
  neg <- label_clusters(z, threshold = 3, tail = "neg")
  expect_equal(neg$table$size, 1L)
  expect_equal(neg$table$mass, 2)
  expect_equal(neg$table$peak_value, -5)
  expect_equal(sum(neg$labels > 0), 1L)
  none <- label_clusters(z, threshold = 10)
  expect_equal(nrow(none$table), 0L)
})

test_that("tfce_fwer null maxima match a per-draw R reference", {
  set.seed(104)
  dims <- c(6, 6, 5)
  subj <- replicate(8, smooth_field(dims) + 0.3, simplify = FALSE)
  spec <- null_one_sample(subj, n_perm = 15, seed = 3)
  mask <- array(as.logical(spec$mask), dims)
  mask_idx <- which(mask)
  for (tail in c("pos", "neg", "two")) {
    res <- tfce_fwer(spec$stat_vol, mask = mask, null_spec = spec, tail = tail)
    ref <- vapply(seq_len(15), function(b) {
      d <- numeric(prod(dims))
      d[mask_idx] <- spec$draw(b)
      if (tail == "two") {
        max(pmax(ns$.tfce_from_vec_r(d, dims, mask, 2, 0.5, 0.1, "pos"),
                 ns$.tfce_from_vec_r(d, dims, mask, 2, 0.5, 0.1, "neg")))
      } else {
        max(ns$.tfce_from_vec_r(d, dims, mask, 2, 0.5, 0.1, tail))
      }
    }, numeric(1))
    expect_equal(res$null_max, ref, tolerance = 1e-12)
    obs <- as.numeric(if (tail == "two") res$extras$tfce_two else res$extras$tfce)
    p_ref <- vapply(obs[mask_idx], function(o) (1 + sum(ref >= o)) / 16, numeric(1))
    expect_equal(as.numeric(res$pvals)[mask_idx], p_ref)
  }
})
