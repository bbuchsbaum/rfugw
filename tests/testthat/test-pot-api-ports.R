test_that("partial GW satisfies partial-mass constraints", {
  skip_if_not_installed("lpSolve")

  set.seed(101)
  ns <- 8L
  nt <- 7L
  Xs <- matrix(rnorm(ns * 3L), ns, 3L)
  Xt <- matrix(rnorm(nt * 3L), nt, 3L)
  C1 <- as.matrix(dist(Xs)); C1 <- C1 / max(C1)
  C2 <- as.matrix(dist(Xt)); C2 <- C2 / max(C2)
  p <- rep(1 / ns, ns)
  q <- rep(1 / nt, nt)
  m <- 0.65

  out <- rfugw::partial_gromov_wasserstein(
    C1 = C1,
    C2 = C2,
    p = p,
    q = q,
    m = m,
    numItermax = 120L,
    tol = 1e-8,
    log = TRUE,
    lp_solver = "lp_matrix"
  )

  expect_equal(sum(out$plan), m, tolerance = 1e-5)
  expect_true(all(rowSums(out$plan) <= p + 5e-4))
  expect_true(all(colSums(out$plan) <= q + 5e-4))
  expect_true(all(out$plan >= -1e-12))
  expect_true(is.finite(out$partial_gw_dist))
  expect_equal(rfugw::partial_gromov_wasserstein2(
    C1 = C1,
    C2 = C2,
    p = p,
    q = q,
    m = m,
    numItermax = 80L,
    tol = 1e-8,
    lp_solver = "lp_matrix"
  ) > 0, TRUE)
})

test_that("partial FGW objective wrapper agrees with solver output", {
  skip_if_not_installed("lpSolve")

  set.seed(102)
  ns <- 9L
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
  q <- rep(1 / nt, nt)

  out <- rfugw::partial_fused_gromov_wasserstein(
    M = M,
    C1 = C1,
    C2 = C2,
    p = p,
    q = q,
    m = 0.6,
    alpha = 0.5,
    numItermax = 120L,
    tol = 1e-8,
    lp_solver = "lp_matrix",
    log = TRUE
  )
  val <- rfugw::partial_fused_gromov_wasserstein2(
    M = M,
    C1 = C1,
    C2 = C2,
    p = p,
    q = q,
    m = 0.6,
    alpha = 0.5,
    numItermax = 120L,
    tol = 1e-8,
    lp_solver = "lp_matrix"
  )

  expect_equal(val, out$partial_fgw_dist, tolerance = 1e-8)
  expect_equal(sum(out$plan), 0.6, tolerance = 1e-5)
  expect_true(is.finite(out$lin_loss))
  expect_true(is.finite(out$quad_loss))
})

test_that("entropic partial GW/FGW produce valid couplings", {
  set.seed(103)
  ns <- 8L
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
  q <- rep(1 / nt, nt)
  m <- 0.7

  out_gw <- rfugw::entropic_partial_gromov_wasserstein(
    C1 = C1,
    C2 = C2,
    p = p,
    q = q,
    m = m,
    reg = 0.2,
    numItermax = 80L,
    tol = 1e-7,
    log = TRUE
  )

  out_fgw <- rfugw::entropic_partial_fused_gromov_wasserstein(
    M = M,
    C1 = C1,
    C2 = C2,
    p = p,
    q = q,
    m = m,
    reg = 0.2,
    alpha = 0.6,
    numItermax = 80L,
    tol = 1e-7,
    log = TRUE
  )

  expect_equal(sum(out_gw$plan), m, tolerance = 1e-4)
  expect_equal(sum(out_fgw$plan), m, tolerance = 1e-4)
  expect_true(all(rowSums(out_gw$plan) <= p + 1e-6))
  expect_true(all(colSums(out_gw$plan) <= q + 1e-6))
  expect_true(all(is.finite(out_fgw$plan)))
  expect_equal(
    rfugw::entropic_partial_fused_gromov_wasserstein2(
      M = M,
      C1 = C1,
      C2 = C2,
      p = p,
      q = q,
      m = m,
      reg = 0.2,
      alpha = 0.6,
      numItermax = 60L,
      tol = 1e-7
    ) > 0,
    TRUE
  )
})

test_that("non-entropic semirelaxed GW/FGW keep row marginals", {
  set.seed(104)
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

  out_gw <- rfugw::semirelaxed_gromov_wasserstein(
    C1 = C1,
    C2 = C2,
    p = p,
    max_iter = 200L,
    tol_rel = 1e-9,
    tol_abs = 1e-9
  )
  out_fgw <- rfugw::semirelaxed_fused_gromov_wasserstein(
    M = M,
    C1 = C1,
    C2 = C2,
    p = p,
    alpha = 0.55,
    max_iter = 200L,
    tol_rel = 1e-9,
    tol_abs = 1e-9
  )

  expect_equal(unname(rowSums(out_gw$plan)), p, tolerance = 1e-10)
  expect_equal(unname(rowSums(out_fgw$plan)), p, tolerance = 1e-10)
  expect_true(all(out_gw$plan >= -1e-12))
  expect_true(all(out_fgw$plan >= -1e-12))
  expect_true(is.finite(out_gw$srgw_dist))
  expect_true(is.finite(out_fgw$srfgw_dist))
  expect_equal(
    rfugw::semirelaxed_fused_gromov_wasserstein2(
      M = M,
      C1 = C1,
      C2 = C2,
      p = p,
      alpha = 0.55,
      max_iter = 200L,
      tol_rel = 1e-9,
      tol_abs = 1e-9
    ),
    out_fgw$srfgw_dist,
    tolerance = 1e-12
  )
})

test_that("entropic barycenter wrappers return coherent objects", {
  set.seed(105)
  make_graph <- function(n) {
    X <- matrix(rnorm(n * 3L), n, 3L)
    F <- matrix(rnorm(n * 2L), n, 2L)
    C <- as.matrix(dist(X)); C <- C / max(C)
    list(C = C, F = F, p = rep(1 / n, n))
  }
  g1 <- make_graph(10L)
  g2 <- make_graph(11L)
  g3 <- make_graph(12L)

  out_gw <- rfugw::entropic_gromov_barycenters(
    N = 7L,
    Cs = list(g1$C, g2$C, g3$C),
    ps = list(g1$p, g2$p, g3$p),
    epsilon = 0.05,
    max_iter = 5L,
    tol = 1e-8,
    log = TRUE
  )
  out_fgw <- rfugw::entropic_fused_gromov_barycenters(
    N = 7L,
    Ys = list(g1$F, g2$F, g3$F),
    Cs = list(g1$C, g2$C, g3$C),
    ps = list(g1$p, g2$p, g3$p),
    alpha = 0.5,
    epsilon = 0.05,
    max_iter = 5L,
    tol = 1e-8,
    log = TRUE
  )

  expect_equal(dim(out_gw$C), c(7L, 7L))
  expect_equal(dim(out_fgw$C), c(7L, 7L))
  expect_equal(dim(out_fgw$Y), c(7L, 2L))
  expect_equal(length(out_gw$couplings), 3L)
  expect_equal(length(out_fgw$couplings), 3L)
})

test_that("sampled GW returns valid coupling", {
  set.seed(106)
  n <- 18L
  X1 <- matrix(rnorm(n * 3L), n, 3L)
  X2 <- matrix(rnorm(n * 3L), n, 3L)
  C1 <- as.matrix(dist(X1)); C1 <- C1 / max(C1)
  C2 <- as.matrix(dist(X2)); C2 <- C2 / max(C2)
  p <- rep(1 / n, n)
  q <- rep(1 / n, n)

  out <- rfugw::sampled_gromov_wasserstein(
    C1 = C1,
    C2 = C2,
    p = p,
    q = q,
    nb_samples_grad = c(8L, 2L),
    epsilon = 0.1,
    max_iter = 80L,
    random_state = 106,
    log = TRUE
  )

  expect_equal(rowSums(out$plan), p, tolerance = 2e-4)
  expect_equal(colSums(out$plan), q, tolerance = 2e-4)
  expect_true(is.finite(out$gw_dist_estimated))
  expect_gte(out$iterations, 1L)
})

test_that("lowrank GW samples wrapper returns factors and diagnostics", {
  set.seed(107)
  Xs <- matrix(rnorm(40L), 10L, 4L)
  Xt <- matrix(rnorm(48L), 12L, 4L)

  out <- rfugw::lowrank_gromov_wasserstein_samples(
    X_s = Xs,
    X_t = Xt,
    reg = 0.05,
    rank = 5L,
    numItermax = 100L,
    log = TRUE
  )

  expect_equal(dim(out$Q), c(10L, 5L))
  expect_equal(dim(out$R), c(12L, 5L))
  expect_equal(length(out$g), 5L)
  expect_true(is.finite(out$value_quad))
  expect_true(is.finite(out$value))
})

test_that("across-space and UCOOT wrappers produce finite objectives", {
  set.seed(108)
  X <- matrix(rnorm(8L * 5L), 8L, 5L)
  Y <- matrix(rnorm(7L * 6L), 7L, 6L)

  out_joint <- rfugw::fused_unbalanced_across_spaces_divergence(
    X = X,
    Y = Y,
    reg_marginals = c(10, 8),
    epsilon = c(0.05, 0.03),
    reg_type = "joint",
    divergence = "kl",
    unbalanced_solver = "sinkhorn",
    max_iter = 40L,
    tol = 1e-7,
    max_iter_ot = 200L,
    tol_ot = 1e-7,
    log = TRUE
  )

  out_ucoot <- rfugw::unbalanced_co_optimal_transport(
    X = X,
    Y = Y,
    reg_marginals = c(10, 8),
    epsilon = c(0.05, 0.03),
    divergence = "kl",
    unbalanced_solver = "sinkhorn",
    max_iter = 40L,
    tol = 1e-7,
    max_iter_ot = 200L,
    tol_ot = 1e-7,
    log = TRUE
  )

  val_ucoot <- rfugw::unbalanced_co_optimal_transport2(
    X = X,
    Y = Y,
    reg_marginals = c(10, 8),
    epsilon = c(0.05, 0.03),
    divergence = "kl",
    unbalanced_solver = "sinkhorn",
    max_iter = 40L,
    tol = 1e-7,
    max_iter_ot = 200L,
    tol_ot = 1e-7
  )

  expect_equal(dim(out_joint$pi_samp), c(8L, 7L))
  expect_equal(dim(out_joint$pi_feat), c(5L, 6L))
  expect_equal(dim(out_ucoot$pi_samp), c(8L, 7L))
  expect_equal(dim(out_ucoot$pi_feat), c(5L, 6L))
  expect_true(is.finite(out_joint$ucoot_cost))
  expect_true(is.finite(out_ucoot$ucoot_cost))
  expect_true(is.finite(val_ucoot))
})
