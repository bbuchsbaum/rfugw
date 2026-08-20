make_gw_pair <- function(n, seed, d = 2L) {
  set.seed(seed)
  X1 <- matrix(rnorm(n * d), n, d)
  X2 <- matrix(rnorm(n * d), n, d)
  C1 <- as.matrix(dist(X1))
  C2 <- as.matrix(dist(X2))
  C1 <- C1 / max(C1)
  C2 <- C2 / max(C2)
  list(
    X1 = X1,
    X2 = X2,
    C1 = C1,
    C2 = C2,
    p = rep(1 / n, n),
    q = rep(1 / n, n)
  )
}

lowrank_rel_error <- function(T, rank) {
  r <- as.integer(rank)
  sv <- svd(T, nu = r, nv = r)
  s <- pmax(sv$d[seq_len(r)], 0)
  S <- diag(sqrt(s), nrow = r, ncol = r)
  Q <- sv$u[, seq_len(r), drop = FALSE] %*% S
  R <- sv$v[, seq_len(r), drop = FALSE] %*% S
  Tr <- Q %*% t(R)
  normT <- sqrt(sum(T^2))
  sqrt(sum((T - Tr)^2)) / normT
}

test_that("unusable sampled budgets error or warn explicitly", {
  expect_error(rfugw:::.parse_sampled_budget(0L, 8L, 8L), ">= 1")
  expect_error(rfugw:::.parse_sampled_budget(-2L, 8L, 8L), ">= 1")
  expect_error(rfugw:::.parse_sampled_budget(c(2L, 0L), 8L, 8L), "length-2")
  expect_error(rfugw:::.parse_sampled_budget(c(2L, 3L, 4L), 8L, 8L), "length-2")
  expect_error(rfugw:::.parse_sampled_budget(NA_integer_, 8L, 8L), ">= 1")

  expect_warning(
    over <- rfugw:::.parse_sampled_budget(c(20L, 15L), 8L, 10L),
    "clamping"
  )
  expect_equal(over$nb_p, 8L)
  expect_equal(over$nb_q, 10L)
  expect_true(over$clamped)

  scalar <- rfugw:::.parse_sampled_budget(5L, 8L, 8L)
  expect_equal(scalar$nb_p, 5L)
  expect_equal(scalar$nb_q, 1L)
  expect_false(scalar$clamped)

  expect_warning(
    remapped <- rfugw:::.parse_sampled_budget(80L, 8L, 6L),
    "clamping"
  )
  expect_equal(remapped$nb_p, 8L)
  expect_equal(remapped$nb_q, 6L)
  expect_true(remapped$clamped)

  d <- make_gw_pair(6L, 20260816L)
  expect_warning(
    sampled_gromov_wasserstein(
      d$C1, d$C2, d$p, d$q,
      nb_samples_grad = c(20L, 20L),
      epsilon = 0.2,
      max_iter = 5L,
      random_state = 1L
    ),
    "clamping"
  )
  expect_warning(
    sampled_gromov_wasserstein_coords(
      d$X1, d$X2, d$p, d$q,
      nb_samples_grad = c(20L, 20L),
      epsilon = 0.2,
      max_iter = 5L,
      random_state = 1L
    ),
    "clamping"
  )
})

test_that("unusable low-rank ranks error or warn explicitly", {
  expect_error(rfugw:::.parse_lowrank_rank(0L, 8L, 6L), ">= 1")
  expect_equal(rfugw:::.parse_lowrank_rank(NULL, 8L, 6L), 6L)
  expect_equal(rfugw:::.parse_lowrank_rank(4L, 8L, 6L), 4L)
  expect_warning(
    over <- rfugw:::.parse_lowrank_rank(20L, 8L, 6L),
    "clamping"
  )
  expect_equal(over, 6L)

  set.seed(20260816L)
  Xs <- matrix(rnorm(24L), 6L, 4L)
  Xt <- matrix(rnorm(32L), 8L, 4L)
  expect_warning(
    lowrank_gromov_wasserstein_samples(
      Xs, Xt, reg = 0.1, rank = 20L, numItermax = 20L
    ),
    "clamping"
  )
})

test_that("sampled GW quality generally improves with budget", {
  n <- 12L
  d <- make_gw_pair(n, 20260816L)
  dense <- entropic_gromov_wasserstein(
    d$C1, d$C2, d$p, d$q,
    epsilon = 0.1,
    max_iter = 80L,
    tol = 1e-7
  )
  gw_ref <- ot_gw_square(d$C1, d$C2, dense$plan)

  seeds <- 20260816L + seq_len(6L)
  gap_low <- numeric(length(seeds))
  frob_low <- numeric(length(seeds))
  for (i in seq_along(seeds)) {
    low <- sampled_gromov_wasserstein(
      d$C1, d$C2, d$p, d$q,
      nb_samples_grad = c(2L, 1L),
      epsilon = 0.1,
      max_iter = 40L,
      random_state = seeds[[i]],
      log = TRUE
    )
    gap_low[[i]] <- abs(ot_gw_square(d$C1, d$C2, low$plan) - gw_ref)
    frob_low[[i]] <- sqrt(sum((low$plan - dense$plan)^2))
  }
  high <- sampled_gromov_wasserstein(
    d$C1, d$C2, d$p, d$q,
    nb_samples_grad = c(n, n),
    epsilon = 0.1,
    max_iter = 40L,
    random_state = seeds[[1]],
    log = TRUE
  )
  gap_high <- abs(ot_gw_square(d$C1, d$C2, high$plan) - gw_ref)
  frob_high <- sqrt(sum((high$plan - dense$plan)^2))

  expect_lt(gap_high, stats::median(gap_low))
  expect_lt(frob_high, stats::median(frob_low))
})

test_that("low-rank reconstruction error decreases with rank", {
  d <- make_gw_pair(10L, 20260816L)
  dense <- entropic_gromov_wasserstein(
    d$C1, d$C2, d$p, d$q,
    epsilon = 0.08,
    max_iter = 80L,
    tol = 1e-7
  )
  err2 <- lowrank_rel_error(dense$plan, 2L)
  err4 <- lowrank_rel_error(dense$plan, 4L)
  err8 <- lowrank_rel_error(dense$plan, 8L)
  expect_lt(err4, err2)
  expect_lt(err8, err4)
  expect_lt(err8, 0.15)

  out <- lowrank_gromov_wasserstein_samples(
    d$X1, d$X2,
    a = d$p,
    b = d$q,
    reg = 0.08,
    rank = 4L,
    numItermax = 80L,
    log = TRUE
  )
  recon <- out$Q %*% t(out$R)
  rel <- sqrt(sum((out$plan - recon)^2)) / sqrt(sum(out$plan^2))
  expect_lt(rel, err2)
})

test_that("coordinate and graph inputs store less than dense structure costs", {
  skip_if_not_installed("Matrix")
  n <- 40L
  d <- 3L
  k <- 5L
  set.seed(20260816L)
  X <- matrix(rnorm(n * d), n, d)
  C <- as.matrix(dist(X))
  W <- matrix(0, n, n)
  D <- C
  for (i in seq_len(n)) {
    keep <- order(D[i, ])[seq_len(k + 1L)]
    keep <- keep[keep != i]
    W[i, keep] <- 1
  }
  W <- pmax(W, t(W))
  Ws <- Matrix::Matrix(W, sparse = TRUE)
  coords <- matrix(rnorm(n * k), n, k)

  expect_lt(as.numeric(object.size(X)), as.numeric(object.size(C)))
  expect_lt(as.numeric(object.size(coords)), as.numeric(object.size(C)))
  expect_lt(as.numeric(object.size(Ws)), as.numeric(object.size(C)))
  expect_lt(
    as.numeric(object.size(Ws)) + as.numeric(object.size(coords)),
    2 * as.numeric(object.size(C))
  )
})
