test_that("fgw_entropic log Sinkhorn matches scaling Sinkhorn", {
  fx <- read_fixture("fgw_entropic_square_fixture.json")

  out_scaling <- rfugw::fgw_entropic(
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
    sinkhorn_method = "scaling",
    symmetric = TRUE,
    solver = "PGD"
  )

  out_log <- rfugw::fgw_entropic(
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
    sinkhorn_method = "log",
    symmetric = TRUE,
    solver = "PGD"
  )

  expect_equal(out_log$fgw_dist, out_scaling$fgw_dist, tolerance = 1e-6)
  expect_equal(rowSums(out_log$plan), fx$inputs$p, tolerance = 1e-7)
  expect_equal(colSums(out_log$plan), fx$inputs$q, tolerance = 1e-7)
})

test_that("fgw_exact_cg decreases objective and satisfies marginals", {
  set.seed(123)

  n <- 10
  X <- matrix(rnorm(n * 3), n, 3)
  C1 <- as.matrix(dist(X))
  C1 <- C1 / max(C1)
  perm <- rev(seq_len(n))
  C2 <- C1[perm, perm]

  Y <- matrix(rnorm(n * 2), n, 2)
  M <- as.matrix(dist(rbind(Y, Y[perm, ])))[1:n, (n + 1):(2 * n)]
  M <- M / max(M)

  init <- rfugw::fgw_exact_cg(M, C1, C2, alpha = 0.5, max_iter = 0)
  out <- rfugw::fgw_exact_cg(M, C1, C2, alpha = 0.5, max_iter = 80)

  expect_true(is.finite(out$fgw_dist))
  expect_lte(out$fgw_dist, init$fgw_dist + 1e-10)
  expect_equal(rowSums(out$plan), rep(1 / n, n), tolerance = 2e-6)
  expect_equal(colSums(out$plan), rep(1 / n, n), tolerance = 2e-6)
})

test_that("fgw_entropic mixed precision stays close to double precision", {
  fx <- read_fixture("fgw_entropic_square_fixture.json")

  out_double <- rfugw::fgw_entropic(
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
    sinkhorn_method = "scaling",
    precision = "double",
    symmetric = TRUE,
    solver = "PGD"
  )
  out_mixed <- rfugw::fgw_entropic(
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
    sinkhorn_method = "scaling",
    precision = "mixed",
    symmetric = TRUE,
    solver = "PGD"
  )

  expect_equal(out_mixed$fgw_dist, out_double$fgw_dist, tolerance = 5e-5)
  expect_equal(rowSums(out_mixed$plan), fx$inputs$p, tolerance = 2e-6)
  expect_equal(colSums(out_mixed$plan), fx$inputs$q, tolerance = 2e-6)
})

test_that("fgw_entropic low-rank full-rank fallback matches dense", {
  fx <- read_fixture("fgw_entropic_square_fixture.json")
  ns <- nrow(fx$inputs$C1)
  nt <- nrow(fx$inputs$C2)

  out_dense <- rfugw::fgw_entropic(
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
    sinkhorn_method = "scaling",
    precision = "double",
    symmetric = TRUE,
    solver = "PGD",
    structure_rank = 0L
  )
  out_full_rank <- rfugw::fgw_entropic(
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
    sinkhorn_method = "scaling",
    precision = "double",
    symmetric = TRUE,
    solver = "PGD",
    structure_rank = max(ns, nt)
  )

  expect_equal(out_full_rank$fgw_dist, out_dense$fgw_dist, tolerance = 1e-6)
  expect_equal(out_full_rank$plan, out_dense$plan, tolerance = 1e-5)
})

test_that("fgw_entropic supports warm-start init_plan", {
  fx <- read_fixture("fgw_entropic_square_fixture.json")

  out <- rfugw::fgw_entropic(
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
    sinkhorn_method = "scaling",
    precision = "double",
    solver = "PGD",
    check_every = 1L
  )
  out_ws <- rfugw::fgw_entropic(
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
    init_plan = out$plan,
    sinkhorn_method = "scaling",
    precision = "double",
    solver = "PGD",
    check_every = 1L
  )

  expect_equal(out_ws$fgw_dist, out$fgw_dist, tolerance = 1e-10)
  expect_equal(out_ws$plan, out$plan, tolerance = 1e-8)
  expect_lte(out_ws$iterations, out$iterations)
})

test_that("fgw_exact_cg lp_matrix matches lp_transport", {
  skip_if_not_installed("lpSolve")
  set.seed(9)

  n <- 8
  X <- matrix(rnorm(n * 3), n, 3)
  C1 <- as.matrix(dist(X)); C1 <- C1 / max(C1)
  C2 <- C1[rev(seq_len(n)), rev(seq_len(n))]
  Y <- matrix(rnorm(n * 2), n, 2)
  M <- as.matrix(dist(rbind(Y, Y[rev(seq_len(n)), ])))[1:n, (n + 1):(2 * n)]
  M <- M / max(M)

  out_matrix <- rfugw::fgw_exact_cg(
    M, C1, C2,
    alpha = 0.5,
    max_iter = 60,
    lp_solver = "lp_matrix"
  )
  out_transport <- rfugw::fgw_exact_cg(
    M, C1, C2,
    alpha = 0.5,
    max_iter = 60,
    lp_solver = "lp_transport"
  )

  expect_equal(out_matrix$fgw_dist, out_transport$fgw_dist, tolerance = 1e-6)
  expect_equal(out_matrix$plan, out_transport$plan, tolerance = 1e-5)
})

test_that("fgw_exact_cg cpp_transport matches lp_transport", {
  skip_if_not_installed("lpSolve")
  set.seed(21)

  n <- 8
  X <- matrix(rnorm(n * 3), n, 3)
  C1 <- as.matrix(dist(X)); C1 <- C1 / max(C1)
  C2 <- C1[rev(seq_len(n)), rev(seq_len(n))]
  Y <- matrix(rnorm(n * 2), n, 2)
  M <- as.matrix(dist(rbind(Y, Y[rev(seq_len(n)), ])))[1:n, (n + 1):(2 * n)]
  M <- M / max(M)

  out_cpp <- rfugw::fgw_exact_cg(
    M, C1, C2,
    alpha = 0.5,
    max_iter = 60,
    lp_solver = "cpp_transport"
  )
  out_transport <- rfugw::fgw_exact_cg(
    M, C1, C2,
    alpha = 0.5,
    max_iter = 60,
    lp_solver = "lp_transport"
  )

  expect_equal(out_cpp$fgw_dist, out_transport$fgw_dist, tolerance = 1e-6)
  expect_equal(out_cpp$plan, out_transport$plan, tolerance = 1e-5)
})
