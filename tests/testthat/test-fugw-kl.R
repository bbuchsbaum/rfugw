test_that("fugw_kl matches POT fixture (KL + Sinkhorn)", {
  fx <- read_fixture("fugw_kl_sinkhorn_fixture.json")
  pi_samp_ref <- as.matrix(fx$outputs$pi_samp)
  pi_feat_ref <- as.matrix(fx$outputs$pi_feat)

  out <- rfugw::fugw_kl(
    Cx = fx$inputs$Cx,
    Cy = fx$inputs$Cy,
    wx = fx$inputs$wx,
    wy = fx$inputs$wy,
    reg_marginals = unlist(fx$params$reg_marginals),
    epsilon = fx$params$epsilon,
    alpha = fx$params$alpha,
    M = fx$inputs$M,
    max_iter = fx$params$max_iter,
    tol = fx$params$tol,
    max_iter_ot = fx$params$max_iter_ot,
    tol_ot = fx$params$tol_ot
  )

  expect_equal(dim(out$pi_samp), dim(pi_samp_ref))
  expect_equal(dim(out$pi_feat), dim(pi_feat_ref))
  expect_equal(out$pi_samp, pi_samp_ref, tolerance = 3e-4)
  expect_equal(out$pi_feat, pi_feat_ref, tolerance = 3e-4)
  expect_equal(out$fugw_cost, fx$outputs$fugw_cost, tolerance = 5e-4)
})

test_that("fugw_kl mixed precision stays close to double", {
  set.seed(221)
  nx <- 26L
  ny <- 22L
  Xx <- matrix(rnorm(nx * 3L), nx, 3L)
  Xy <- matrix(rnorm(ny * 3L), ny, 3L)
  Fx <- matrix(rnorm(nx * 4L), nx, 4L)
  Fy <- matrix(rnorm(ny * 4L), ny, 4L)

  Cx <- as.matrix(dist(Xx)); Cx <- Cx / max(Cx)
  Cy <- as.matrix(dist(Xy)); Cy <- Cy / max(Cy)
  M <- as.matrix(dist(rbind(Fx, Fy)))[seq_len(nx), nx + seq_len(ny)]
  M <- M / max(M)
  wx <- rep(1 / nx, nx)
  wy <- rep(1 / ny, ny)

  out_d <- rfugw::fugw_kl(
    Cx = Cx,
    Cy = Cy,
    wx = wx,
    wy = wy,
    reg_marginals = c(30, 12),
    epsilon = 1e-2,
    alpha = 0.6,
    M = M,
    max_iter = 60L,
    tol = 1e-8,
    max_iter_ot = 220L,
    tol_ot = 1e-8,
    check_every = 1L,
    precision = "double"
  )

  out_m <- rfugw::fugw_kl(
    Cx = Cx,
    Cy = Cy,
    wx = wx,
    wy = wy,
    reg_marginals = c(30, 12),
    epsilon = 1e-2,
    alpha = 0.6,
    M = M,
    max_iter = 60L,
    tol = 1e-8,
    max_iter_ot = 220L,
    tol_ot = 1e-8,
    check_every = 1L,
    precision = "mixed"
  )

  expect_equal(out_m$fugw_cost, out_d$fugw_cost, tolerance = 8e-4)
  expect_equal(out_m$linear_cost, out_d$linear_cost, tolerance = 8e-4)
  expect_equal(out_m$pi_samp, out_d$pi_samp, tolerance = 4e-4)
  expect_equal(out_m$pi_feat, out_d$pi_feat, tolerance = 4e-4)
})

test_that("fugw_kl exposes inner Sinkhorn telemetry", {
  set.seed(901)
  n <- 18L
  X1 <- matrix(rnorm(n * 3L), n, 3L)
  X2 <- matrix(rnorm(n * 3L), n, 3L)
  Cx <- as.matrix(dist(X1)); Cx <- Cx / max(Cx)
  Cy <- as.matrix(dist(X2)); Cy <- Cy / max(Cy)
  M <- matrix(1, n, n)
  idx <- seq_len(n)
  M[cbind(idx, rev(idx))] <- 0

  out <- rfugw::fugw_kl(
    Cx = Cx,
    Cy = Cy,
    reg_marginals = c(20, 20),
    epsilon = 1e-2,
    alpha = 0.5,
    M = M,
    max_iter = 30L,
    tol = 1e-7,
    max_iter_ot = 120L,
    tol_ot = 1e-7
  )

  expect_equal(length(out$inner_iters_feat), out$iterations)
  expect_equal(length(out$inner_iters_samp), out$iterations)
  expect_true(is.logical(out$inner_warm_feat))
  expect_true(is.logical(out$inner_warm_samp))
  expect_true(is.logical(out$inner_warm_fallback_feat))
  expect_true(is.logical(out$inner_warm_fallback_samp))
  expect_equal(
    out$inner_iters_total,
    sum(out$inner_iters_feat) + sum(out$inner_iters_samp)
  )
})
