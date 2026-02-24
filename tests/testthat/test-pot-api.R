test_that("entropic_gromov_wasserstein matches FGW zero-feature formulation", {
  fx <- read_fixture("fgw_entropic_square_fixture.json")
  M0 <- matrix(0, nrow = nrow(fx$inputs$C1), ncol = nrow(fx$inputs$C2))

  out_ref <- rfugw::fgw_entropic(
    M = M0,
    C1 = fx$inputs$C1,
    C2 = fx$inputs$C2,
    p = fx$inputs$p,
    q = fx$inputs$q,
    alpha = 1,
    epsilon = fx$params$epsilon,
    max_iter = fx$params$max_iter,
    tol = fx$params$tol,
    sinkhorn_max_iter = fx$params$sinkhorn_numItermax,
    sinkhorn_tol = fx$params$sinkhorn_stopThr,
    sinkhorn_method = "scaling",
    precision = "double",
    solver = "PGD"
  )
  out <- rfugw::entropic_gromov_wasserstein(
    C1 = fx$inputs$C1,
    C2 = fx$inputs$C2,
    p = fx$inputs$p,
    q = fx$inputs$q,
    epsilon = fx$params$epsilon,
    max_iter = fx$params$max_iter,
    tol = fx$params$tol,
    sinkhorn_max_iter = fx$params$sinkhorn_numItermax,
    sinkhorn_tol = fx$params$sinkhorn_stopThr,
    sinkhorn_method = "scaling",
    precision = "double",
    solver = "PGD"
  )

  expect_equal(out$plan, out_ref$plan, tolerance = 1e-10)
  expect_equal(out$gw_dist, out_ref$fgw_dist, tolerance = 1e-12)
  expect_equal(
    rfugw::entropic_gromov_wasserstein2(
      C1 = fx$inputs$C1,
      C2 = fx$inputs$C2,
      p = fx$inputs$p,
      q = fx$inputs$q,
      epsilon = fx$params$epsilon,
      max_iter = fx$params$max_iter,
      tol = fx$params$tol,
      sinkhorn_max_iter = fx$params$sinkhorn_numItermax,
      sinkhorn_tol = fx$params$sinkhorn_stopThr,
      sinkhorn_method = "scaling",
      precision = "double",
      solver = "PGD"
    ),
    out$gw_dist,
    tolerance = 1e-12
  )
})

test_that("gromov_wasserstein matches exact FGW zero-feature formulation", {
  set.seed(42)
  n <- 8L
  X <- matrix(rnorm(n * 3L), n, 3L)
  C1 <- as.matrix(dist(X)); C1 <- C1 / max(C1)
  C2 <- C1[rev(seq_len(n)), rev(seq_len(n))]
  p <- rep(1 / n, n)
  q <- rep(1 / n, n)
  M0 <- matrix(0, n, n)

  out_ref <- rfugw::fgw_exact_cg(
    M = M0,
    C1 = C1,
    C2 = C2,
    p = p,
    q = q,
    alpha = 1,
    max_iter = 60L,
    tol_rel = 1e-10,
    tol_abs = 1e-10,
    lp_solver = "cpp_transport"
  )
  out <- rfugw::gromov_wasserstein(
    C1 = C1,
    C2 = C2,
    p = p,
    q = q,
    max_iter = 60L,
    tol_rel = 1e-10,
    tol_abs = 1e-10,
    lp_solver = "cpp_transport"
  )

  expect_equal(out$plan, out_ref$plan, tolerance = 1e-10)
  expect_equal(out$gw_dist, out_ref$fgw_dist, tolerance = 1e-12)
  expect_equal(rfugw::gromov_wasserstein2(
    C1 = C1,
    C2 = C2,
    p = p,
    q = q,
    max_iter = 60L,
    tol_rel = 1e-10,
    tol_abs = 1e-10,
    lp_solver = "cpp_transport"
  ), out$gw_dist, tolerance = 1e-12)
})

test_that("POT-compatible alias names route to existing fast solvers", {
  fx <- read_fixture("fgw_entropic_square_fixture.json")
  out_alias <- rfugw::entropic_fused_gromov_wasserstein(
    M = fx$inputs$M,
    C1 = fx$inputs$C1,
    C2 = fx$inputs$C2,
    p = fx$inputs$p,
    q = fx$inputs$q,
    alpha = fx$params$alpha,
    epsilon = fx$params$epsilon,
    max_iter = fx$params$max_iter,
    tol = fx$params$tol,
    sinkhorn_max_iter = fx$params$sinkhorn_numItermax,
    sinkhorn_tol = fx$params$sinkhorn_stopThr,
    precision = "double"
  )
  out_base <- rfugw::fgw_entropic(
    M = fx$inputs$M,
    C1 = fx$inputs$C1,
    C2 = fx$inputs$C2,
    p = fx$inputs$p,
    q = fx$inputs$q,
    alpha = fx$params$alpha,
    epsilon = fx$params$epsilon,
    max_iter = fx$params$max_iter,
    tol = fx$params$tol,
    sinkhorn_max_iter = fx$params$sinkhorn_numItermax,
    sinkhorn_tol = fx$params$sinkhorn_stopThr,
    precision = "double"
  )
  expect_equal(out_alias$plan, out_base$plan, tolerance = 1e-12)
  expect_equal(out_alias$fgw_dist, out_base$fgw_dist, tolerance = 1e-12)

  fxu <- read_fixture("fugw_kl_sinkhorn_fixture.json")
  out_u_alias <- rfugw::fused_unbalanced_gromov_wasserstein(
    Cx = fxu$inputs$Cx,
    Cy = fxu$inputs$Cy,
    wx = fxu$inputs$wx,
    wy = fxu$inputs$wy,
    reg_marginals = unlist(fxu$params$reg_marginals),
    epsilon = fxu$params$epsilon,
    alpha = fxu$params$alpha,
    M = fxu$inputs$M,
    max_iter = fxu$params$max_iter,
    tol = fxu$params$tol,
    max_iter_ot = fxu$params$max_iter_ot,
    tol_ot = fxu$params$tol_ot
  )
  out_u_base <- rfugw::fugw_kl(
    Cx = fxu$inputs$Cx,
    Cy = fxu$inputs$Cy,
    wx = fxu$inputs$wx,
    wy = fxu$inputs$wy,
    reg_marginals = unlist(fxu$params$reg_marginals),
    epsilon = fxu$params$epsilon,
    alpha = fxu$params$alpha,
    M = fxu$inputs$M,
    max_iter = fxu$params$max_iter,
    tol = fxu$params$tol,
    max_iter_ot = fxu$params$max_iter_ot,
    tol_ot = fxu$params$tol_ot
  )
  expect_equal(out_u_alias$pi_samp, out_u_base$pi_samp, tolerance = 1e-12)
  expect_equal(out_u_alias$fugw_cost, out_u_base$fugw_cost, tolerance = 1e-12)
})

test_that("entropic semirelaxed GW keeps row marginals and finite objective", {
  set.seed(7)
  ns <- 9L
  nt <- 7L
  Xs <- matrix(rnorm(ns * 3L), ns, 3L)
  Xt <- matrix(rnorm(nt * 3L), nt, 3L)
  C1 <- as.matrix(dist(Xs)); C1 <- C1 / max(C1)
  C2 <- as.matrix(dist(Xt)); C2 <- C2 / max(C2)
  p <- rep(1 / ns, ns)

  out <- rfugw::entropic_semirelaxed_gromov_wasserstein(
    C1 = C1,
    C2 = C2,
    p = p,
    epsilon = 0.05,
    max_iter = 300L,
    tol = 1e-10,
    check_every = 5L
  )

  expect_equal(unname(rowSums(out$plan)), p, tolerance = 1e-10)
  expect_true(all(is.finite(out$plan)))
  expect_true(all(out$plan >= 0))
  expect_true(is.finite(out$srgw_dist))
  expect_equal(
    rfugw::entropic_semirelaxed_gromov_wasserstein2(
      C1 = C1,
      C2 = C2,
      p = p,
      epsilon = 0.05,
      max_iter = 300L,
      tol = 1e-10,
      check_every = 5L
    ),
    out$srgw_dist,
    tolerance = 1e-12
  )
})

test_that("entropic semirelaxed FGW is invariant to row-wise M offsets", {
  set.seed(8)
  ns <- 10L
  nt <- 8L
  Xs <- matrix(rnorm(ns * 3L), ns, 3L)
  Xt <- matrix(rnorm(nt * 3L), nt, 3L)
  Fs <- matrix(rnorm(ns * 2L), ns, 2L)
  Ft <- matrix(rnorm(nt * 2L), nt, 2L)
  C1 <- as.matrix(dist(Xs)); C1 <- C1 / max(C1)
  C2 <- as.matrix(dist(Xt)); C2 <- C2 / max(C2)
  M <- as.matrix(dist(rbind(Fs, Ft)))[seq_len(ns), ns + seq_len(nt)]
  M <- M / max(M)
  p <- rep(1 / ns, ns)
  row_offset <- rnorm(ns)
  M_shift <- M + tcrossprod(row_offset, rep(1, nt))

  out0 <- rfugw::entropic_semirelaxed_fused_gromov_wasserstein(
    M = M,
    C1 = C1,
    C2 = C2,
    p = p,
    alpha = 0.6,
    epsilon = 0.05,
    max_iter = 300L,
    tol = 1e-10,
    check_every = 5L
  )
  out1 <- rfugw::entropic_semirelaxed_fused_gromov_wasserstein(
    M = M_shift,
    C1 = C1,
    C2 = C2,
    p = p,
    alpha = 0.6,
    epsilon = 0.05,
    max_iter = 300L,
    tol = 1e-10,
    check_every = 5L
  )

  expect_equal(out0$plan, out1$plan, tolerance = 1e-7)
  expect_equal(unname(rowSums(out0$plan)), p, tolerance = 1e-10)
  expect_equal(
    rfugw::entropic_semirelaxed_fused_gromov_wasserstein2(
      M = M,
      C1 = C1,
      C2 = C2,
      p = p,
      alpha = 0.6,
      epsilon = 0.05,
      max_iter = 300L,
      tol = 1e-10,
      check_every = 5L
    ),
    out0$srfgw_dist,
    tolerance = 1e-12
  )
})

test_that("entropic semirelaxed C++ backend matches R reference backend", {
  set.seed(18)
  ns <- 11L
  nt <- 9L
  Xs <- matrix(rnorm(ns * 3L), ns, 3L)
  Xt <- matrix(rnorm(nt * 3L), nt, 3L)
  Fs <- matrix(rnorm(ns * 2L), ns, 2L)
  Ft <- matrix(rnorm(nt * 2L), nt, 2L)
  C1 <- as.matrix(dist(Xs)); C1 <- C1 / max(C1)
  C2 <- as.matrix(dist(Xt)); C2 <- C2 / max(C2)
  M <- as.matrix(dist(rbind(Fs, Ft)))[seq_len(ns), ns + seq_len(nt)]
  M <- M / max(M)
  p <- rep(1 / ns, ns)

  out_cpp <- rfugw::entropic_semirelaxed_fused_gromov_wasserstein(
    M = M,
    C1 = C1,
    C2 = C2,
    p = p,
    alpha = 0.55,
    epsilon = 0.05,
    precision = "double",
    max_iter = 240L,
    tol = 1e-10,
    check_every = 5L,
    backend = "cpp"
  )
  out_r <- rfugw::entropic_semirelaxed_fused_gromov_wasserstein(
    M = M,
    C1 = C1,
    C2 = C2,
    p = p,
    alpha = 0.55,
    epsilon = 0.05,
    precision = "double",
    max_iter = 240L,
    tol = 1e-10,
    check_every = 5L,
    backend = "r"
  )

  expect_equal(out_cpp$srfgw_dist, out_r$srfgw_dist, tolerance = 5e-8)
  expect_equal(out_cpp$lin_loss, out_r$lin_loss, tolerance = 5e-8)
  expect_equal(out_cpp$quad_loss, out_r$quad_loss, tolerance = 5e-8)
  expect_equal(unname(out_cpp$plan), unname(out_r$plan), tolerance = 5e-7)
  expect_equal(unname(rowSums(out_cpp$plan)), p, tolerance = 1e-10)
})

test_that("fgw_barycenters returns valid fixed-support outputs", {
  set.seed(9)
  make_graph <- function(n) {
    X <- matrix(rnorm(n * 3L), n, 3L)
    F <- matrix(rnorm(n * 2L), n, 2L)
    C <- as.matrix(dist(X))
    C <- C / max(C)
    list(C = C, F = F, p = rep(1 / n, n))
  }
  g1 <- make_graph(12L)
  g2 <- make_graph(10L)
  g3 <- make_graph(11L)

  out <- rfugw::fgw_barycenters(
    N = 7L,
    Ys = list(g1$F, g2$F, g3$F),
    Cs = list(g1$C, g2$C, g3$C),
    ps = list(g1$p, g2$p, g3$p),
    alpha = 0.5,
    epsilon = 0.05,
    max_iter = 5L,
    tol = 1e-7,
    precision = "mixed",
    sinkhorn_max_iter = 200L,
    check_every = 5L
  )

  expect_equal(dim(out$X), c(7L, 2L))
  expect_equal(dim(out$C), c(7L, 7L))
  expect_equal(length(out$couplings), 3L)
  expect_equal(sum(out$p), 1, tolerance = 1e-12)
  expect_true(all(vapply(out$couplings, function(T) all(is.finite(T)), logical(1))))
  expect_equal(vapply(out$couplings, nrow, integer(1)), rep(7L, 3L))
  expect_equal(out$history$iter[[nrow(out$history)]], out$iterations)
  expect_true(is.finite(out$objective))
})
