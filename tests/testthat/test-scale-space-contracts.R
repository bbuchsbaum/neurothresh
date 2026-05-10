test_that("as_subject_matrix preserves mask and subject dimensions", {
  imgs <- replicate(3, array(rnorm(27), c(3, 3, 3)), simplify = FALSE)
  mask <- array(TRUE, c(3, 3, 3))
  mask[1, 1, 1] <- FALSE

  prep <- as_subject_matrix(imgs, mask = mask)

  expect_s3_class(prep, "neurothresh_image_matrix")
  expect_equal(dim(prep$X), c(3, 26))
  expect_equal(length(prep$mask_idx), 26)
  expect_equal(prep$dims, c(3L, 3L, 3L))
  expect_equal(dim(prep$coords), c(26, 3))
  expect_identical(prep$source, "list")
})

test_that("as_subject_matrix yields identical X across input shapes", {
  set.seed(7)
  arrs <- lapply(seq_len(4), function(i) array(rnorm(2 * 3 * 4), c(2, 3, 4)))
  mask <- array(TRUE, c(2, 3, 4))
  mask[1, 1, 1] <- FALSE

  arr4d <- array(NA_real_, c(2, 3, 4, 4))
  for (i in seq_along(arrs)) arr4d[, , , i] <- arrs[[i]]

  prep_list <- as_subject_matrix(arrs, mask = mask)
  prep_4d   <- as_subject_matrix(arr4d, mask = mask)

  expect_identical(prep_list$X, prep_4d$X)
  expect_identical(prep_list$mask_idx, prep_4d$mask_idx)
  expect_identical(prep_4d$source, "array4d")

  flat_X <- matrix(NA_real_, nrow = length(arrs), ncol = sum(mask))
  for (i in seq_along(arrs)) flat_X[i, ] <- as.numeric(arrs[[i]])[which(mask)]
  template <- array(NA_real_, dim = c(2, 3, 4))
  prep_mat <- as_subject_matrix(flat_X, mask = mask, template = template)
  expect_identical(prep_list$X, prep_mat$X)
  expect_identical(prep_mat$source, "matrix")
})

test_that("as_subject_matrix accepts NeuroVec and SparseNeuroVec", {
  skip_if_not_installed("neuroim2")
  set.seed(11)
  dims3 <- c(3L, 3L, 3L)
  n_subj <- 3L
  mat <- array(rnorm(prod(dims3) * n_subj), c(dims3, n_subj))
  sp <- neuroim2::NeuroSpace(c(dims3, n_subj))

  nv <- neuroim2::NeuroVec(mat, sp)
  prep_nv <- as_subject_matrix(nv)
  expect_equal(dim(prep_nv$X), c(n_subj, prod(dims3)))
  expect_identical(prep_nv$source, "neurovec")
  expect_false(is.null(prep_nv$space))

  mask3 <- array(TRUE, dims3)
  mask3[1, 1, 1] <- FALSE
  msp <- neuroim2::NeuroSpace(dims3)
  mvol <- neuroim2::LogicalNeuroVol(mask3, msp)

  sv <- neuroim2::SparseNeuroVec(mat, neuroim2::NeuroSpace(c(dims3, n_subj)), mask = mvol)
  prep_sv <- as_subject_matrix(sv)
  expect_equal(ncol(prep_sv$X), sum(mask3))
  expect_identical(prep_sv$source, "sparseneurovec")

  expect_error(
    as_subject_matrix(sv, mask = array(TRUE, dims3)),
    "outside the SparseNeuroVec intrinsic mask"
  )
})

test_that("as_subject_matrix enforces strict mask precedence on non-finite values", {
  arrs <- list(
    array(c(1, 2, NA, 4, 5, 6, 7, 8), c(2, 2, 2)),
    array(rnorm(8), c(2, 2, 2))
  )
  bad_mask <- array(TRUE, c(2, 2, 2))
  expect_error(
    as_subject_matrix(arrs, mask = bad_mask),
    "non-finite"
  )

  good_mask <- array(TRUE, c(2, 2, 2))
  good_mask[1, 2, 1] <- FALSE
  prep <- as_subject_matrix(arrs, mask = good_mask)
  expect_equal(ncol(prep$X), 7L)
})

test_that("one-sample statistic kernel matches base t.test", {
  X <- matrix(c(
    1, 2, 3,
    2, 3, 5,
    3, 4, 7,
    4, 5, 9
  ), nrow = 4, byrow = TRUE)

  got <- neurothresh:::.ss_one_sample_t(X)
  ref <- vapply(seq_len(ncol(X)), function(j) {
    unname(stats::t.test(X[, j])$statistic)
  }, numeric(1))

  expect_equal(got, ref, tolerance = 1e-12)
})

test_that("two-sample statistic kernel matches base Welch t.test", {
  X <- matrix(c(
    1, 2,
    2, 4,
    3, 6,
    5, 8,
    7, 9,
    9, 12
  ), nrow = 6, byrow = TRUE)
  group <- factor(c("a", "a", "a", "b", "b", "b"))

  got <- neurothresh:::.ss_two_sample_t(X, group, variance = "welch")
  ref <- vapply(seq_len(ncol(X)), function(j) {
    # Package convention is level 2 minus level 1.
    -unname(stats::t.test(X[, j] ~ group)$statistic)
  }, numeric(1))

  expect_equal(got, ref, tolerance = 1e-12)
})

test_that("deterministic one-sample null generator repeats by index", {
  imgs <- replicate(5, array(rnorm(8), c(2, 2, 2)), simplify = FALSE)
  nf <- make_null_fun_onesample(imgs, n_perm = 10, seed = 1)

  expect_equal(nf$null_fun(3), nf$null_fun(3), tolerance = 0)
  expect_equal(length(nf$null_fun(3)), sum(nf$mask))
  expect_error(nf$null_fun(11), "Permutation index out of range")
})

test_that("within-subject 2x2 contrast helper matches hand contrast", {
  subjects <- rep(c("s1", "s2"), each = 4)
  task <- rep(rep(c("E", "R"), each = 2), 2)
  model <- rep(c("L", "H"), 4)
  vals <- c(1, 4, 2, 8, 10, 16, 12, 20)
  imgs <- lapply(vals, function(v) array(v, c(2, 2, 2)))

  out <- make_within_subject_contrast_images(
    imgs,
    subject = subjects,
    within = data.frame(task = task, model = model),
    contrast = "task:model"
  )

  expect_equal(length(out), 2)
  expect_equal(as.numeric(out[[1]])[1], (8 - 2) - (4 - 1))
  expect_equal(as.numeric(out[[2]])[1], (20 - 12) - (16 - 10))
})
