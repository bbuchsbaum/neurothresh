# Simulation-based sanity checks for FWER/power behavior

make_kernel <- function(radius = 2, sigma = 1) {
  ax <- seq(-radius, radius)
  g <- outer(ax, ax, function(x, y) exp(-(x^2 + y^2) / (2 * sigma^2)))
  g / sum(g)
}

smooth_field <- function(z, kernel) {
  nr <- nrow(z); nc <- ncol(z)
  kr <- nrow(kernel); kc <- ncol(kernel)
  pad_r <- floor(kr / 2)
  pad_c <- floor(kc / 2)

  zpad <- matrix(0, nrow = nr + 2 * pad_r, ncol = nc + 2 * pad_c)
  zpad[(pad_r + 1):(pad_r + nr), (pad_c + 1):(pad_c + nc)] <- z

  out <- matrix(0, nrow = nr, ncol = nc)
  for (i in seq_len(nr)) {
    for (j in seq_len(nc)) {
      block <- zpad[i:(i + kr - 1), j:(j + kc - 1)]
      out[i, j] <- sum(block * kernel)
    }
  }
  out
}

insert_blob <- function(z, center, radius = 2, amp = 1.5) {
  nr <- nrow(z); nc <- ncol(z)
  rr <- seq_len(nr); cc <- seq_len(nc)
  mask <- outer(rr, cc, function(i, j) (i - center[1])^2 + (j - center[2])^2 <= radius^2)
  z + amp * mask
}

softmax_score <- function(z, kappa = 1) {
  a <- kappa * as.vector(z)
  a_max <- max(a)
  a_max + log(mean(exp(a - a_max)))
}

test_that("softmax score increases with a localized blob", {
  set.seed(42)

  n <- 25
  kernel <- make_kernel(radius = 2, sigma = 1)

  z_null <- smooth_field(matrix(rnorm(n * n), nrow = n), kernel)
  z_blob <- insert_blob(z_null, c(13, 13), radius = 2, amp = 1.5)

  score_null <- softmax_score(z_null, kappa = 1)
  score_blob <- softmax_score(z_blob, kappa = 1)

  expect_gt(score_blob, score_null)
})
