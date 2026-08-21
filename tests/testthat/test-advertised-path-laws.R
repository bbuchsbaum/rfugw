.path_problem <- function(seed = 73L, n = 4L) {
  set.seed(seed)
  X1 <- matrix(rnorm(n * 2L), n, 2L)
  X2 <- matrix(rnorm(n * 2L), n, 2L)
  C1 <- as.matrix(dist(X1)); C1 <- C1 / max(C1)
  C2 <- as.matrix(dist(X2)); C2 <- C2 / max(C2)
  M <- matrix(runif(n * n), n, n)
  list(C1 = C1, C2 = C2, M = M, p = rep(1 / n, n), q = rep(1 / n, n))
}

.expect_equivalent_diagnostics <- function(a, b, atol, label) {
  expect_equal(rfugw_plan(a), rfugw_plan(b), tolerance = atol, info = label)
  expect_equal(rfugw_value(a), rfugw_value(b), tolerance = atol, info = label)
  expect_identical(a$formulation, b$formulation, info = label)
  expect_identical(a$converged, b$converged, info = label)
  expect_identical(a$feasible, b$feasible, info = label)
  expect_identical(a$objective_consistent, b$objective_consistent, info = label)
}

test_that("the advertised path matrix covers every required dimension", {
  matrix_path <- trust_test_resource("numerical-path-matrix.csv")
  paths <- utils::read.csv(matrix_path, stringsAsFactors = FALSE)
  expect_identical(anyDuplicated(paths[c("family", "dimension", "path")]), 0L)
  required <- c(
    "symmetry", "regularization", "sinkhorn", "precision", "start",
    "kernel", "threading", "adapter", "backend", "approximation",
    "formulation", "input"
  )
  expect_setequal(unique(paths$dimension), required)
  expect_true(all(paths$scope %in% c("pr", "nightly", "release")))
  expect_true(all(paths$maturity %in% c("supported", "experimental")))
  expect_true(all(nzchar(paths$comparison)))
})

test_that("native and POT balanced FGW adapters are semantically identical", {
  d <- .path_problem()
  native <- fgw_entropic(
    d$M, d$C1, d$C2, d$p, d$q,
    epsilon = 0.1, max_iter = 40L, tol = 1e-7,
    sinkhorn_tol = 1e-8, precision = "strict_double"
  )
  pot <- entropic_fused_gromov_wasserstein(
    d$M, d$C1, d$C2, d$p, d$q,
    epsilon = 0.1, max_iter = 40L, tol = 1e-7,
    sinkhorn_tol = 1e-8, precision = "strict_double"
  )
  .expect_equivalent_diagnostics(native, pot, 1e-12, "native/POT")
})

test_that("warm starts preserve semantics and cannot bypass validation", {
  d <- .path_problem(seed = 74L)
  cold <- fgw_entropic(
    d$M, d$C1, d$C2, d$p, d$q,
    epsilon = 0.1, max_iter = 80L, tol = 1e-7,
    sinkhorn_tol = 1e-8, precision = "strict_double"
  )
  warm <- fgw_entropic(
    d$M, d$C1, d$C2, d$p, d$q,
    epsilon = 0.1, max_iter = 80L, tol = 1e-7,
    sinkhorn_tol = 1e-8, precision = "strict_double",
    init_plan = cold$plan
  )
  expect_equal(warm$fgw_dist, cold$fgw_dist, tolerance = 2e-8)
  expect_true(warm$objective_consistent)
  expect_true(warm$feasible)
  bad <- cold$plan
  bad[1, 1] <- -1
  expect_error(
    fgw_entropic(d$M, d$C1, d$C2, d$p, d$q, init_plan = bad),
    "finite and nonnegative"
  )
  expect_error(
    fgw_exact_cg(d$M, d$C1, d$C2, d$p, d$q, G0 = bad),
    "finite and nonnegative"
  )
})

test_that("semirelaxed C++ and R backends return equivalent diagnostics", {
  d <- .path_problem(seed = 75L)
  args <- list(
    M = d$M, C1 = d$C1, C2 = d$C2, p = d$p,
    epsilon = 0.1, alpha = 0.5, max_iter = 50L, tol = 1e-7,
    precision = "double"
  )
  cpp <- do.call(entropic_semirelaxed_fused_gromov_wasserstein,
                 c(args, list(backend = "cpp")))
  r <- do.call(entropic_semirelaxed_fused_gromov_wasserstein,
               c(args, list(backend = "r")))
  expect_equal(cpp$srfgw_dist, r$srfgw_dist, tolerance = 2e-6)
  expect_equal(unname(cpp$plan), unname(r$plan), tolerance = 2e-5)
  expect_identical(cpp$feasible, r$feasible)
  expect_true(cpp$objective_consistent)
  expect_true(r$objective_consistent)
})

test_that("optional partial LP paths agree with native certified directions", {
  skip_if_not_installed("lpSolve")
  d <- .path_problem(seed = 76L)
  args <- list(
    M = d$M, C1 = d$C1, C2 = d$C2, p = d$p, q = d$q,
    m = 0.7, alpha = 0.5, numItermax = 8L, tol = 1e-7, log = TRUE
  )
  native <- do.call(partial_fused_gromov_wasserstein,
                    c(args, list(lp_solver = "cpp_transport")))
  oracle <- do.call(partial_fused_gromov_wasserstein,
                    c(args, list(lp_solver = "lp_transport", lp_scale = 1e6)))
  expect_equal(native$partial_fgw_dist, oracle$partial_fgw_dist, tolerance = 2e-5)
  expect_equal(native$plan, oracle$plan, tolerance = 2e-5)
  expect_true(native$objective_consistent)
  expect_true(oracle$objective_consistent)
})

test_that("unsupported fallback labels fail before computation", {
  X <- matrix(c(0, 1, 1, 0), 2, 2)
  expect_error(
    unbalanced_co_optimal_transport(X, X, unbalanced_solver = "mm"),
    "Unsupported `unbalanced_solver`"
  )
  expect_error(
    unbalanced_co_optimal_transport(X, X, divergence = "l2"),
    "divergence.*unsupported"
  )
})
