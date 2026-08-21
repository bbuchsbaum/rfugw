.expect_nested_diagnostics <- function(out) {
  fields <- c(
    "inner_residual", "max_inner_residual", "inner_iterations",
    "inner_converged", "inner_status"
  )
  expect_true(all(fields %in% names(out)))
  expect_length(out$inner_converged, 1L)
  expect_type(out$inner_converged, "logical")
  expect_true(length(out$inner_iterations) == 1L)
  if (is.finite(out$inner_residual) && is.finite(out$max_inner_residual)) {
    expect_gte(out$max_inner_residual, out$inner_residual)
  }
  if (isTRUE(out$converged)) {
    expect_true(out$inner_converged)
  }
}

.nested_partial_fixture <- function() {
  C1 <- matrix(c(0, 1, 2, 1, 0, 1, 2, 1, 0), 3, 3, byrow = TRUE)
  C2 <- matrix(c(0, 2, 1, 2, 0, 1, 1, 1, 0), 3, 3, byrow = TRUE)
  list(
    M = matrix(c(0, 1, 2, 1, 0, 1, 2, 1, 0), 3, 3, byrow = TRUE),
    C1 = C1,
    C2 = C2,
    p = c(0.2, 0.3, 0.5),
    q = c(0.4, 0.1, 0.5)
  )
}

test_that("nested exact-partial failures are propagated before an outer update", {
  d <- .nested_partial_fixture()
  failures <- c("max_iter", "numerical_failure", "infeasible", "nonfinite")

  for (failure in failures) {
    out <- cpp_partial_fgw_exact_square(
      M = d$M,
      C1 = d$C1,
      C2 = d$C2,
      p = d$p,
      q = d$q,
      m = 0.7,
      alpha = 0.5,
      symmetric = TRUE,
      init_plan = matrix(numeric(0), 0, 0),
      max_iter = 5L,
      tol = 1e-8,
      nb_dummies = 1L,
      lp_max_iter = 20000L,
      lp_tol = 1e-12,
      test_inner_failure = failure
    )
    expect_false(out$lp_ok, info = failure)
    expect_false(out$inner_converged, info = failure)
    expect_identical(out$iterations, 0L, info = failure)
    expect_true(all(is.finite(out$plan)), info = failure)
    expect_identical(
      out$inner_status,
      if (identical(failure, "max_iter")) "max_iter" else "numerical_failure",
      info = failure
    )
  }
})

test_that("starved balanced inner solves cannot imply outer convergence", {
  d <- .nested_partial_fixture()
  full <- fgw_entropic(
    d$M, d$C1, d$C2, d$p, d$q,
    epsilon = 0.1, max_iter = 8L, tol = 1e-8,
    sinkhorn_max_iter = 500L, sinkhorn_tol = 1e-7,
    precision = "double"
  )
  starved <- fgw_entropic(
    d$M, d$C1, d$C2, d$p, d$q,
    epsilon = 0.1, max_iter = 8L, tol = 1e-8,
    sinkhorn_max_iter = 1L, sinkhorn_tol = 1e-7,
    precision = "double"
  )

  .expect_nested_diagnostics(full)
  .expect_nested_diagnostics(starved)
  expect_false(starved$inner_converged)
  expect_false(starved$converged)
  expect_gt(norm(full$plan - starved$plan, "F"), 1e-6)
})

test_that("partial POT adapters preserve exact and entropic inner status", {
  d <- .nested_partial_fixture()
  exact <- partial_fused_gromov_wasserstein(
    d$M, d$C1, d$C2, d$p, d$q,
    m = 0.7, numItermax = 5L, log = TRUE
  )
  entropic <- entropic_partial_fused_gromov_wasserstein(
    d$M, d$C1, d$C2, d$p, d$q,
    m = 0.7, reg = 0.1, numItermax = 3L,
    inner_max_iter = 1L, inner_tol = 1e-12, log = TRUE
  )

  .expect_nested_diagnostics(exact)
  .expect_nested_diagnostics(entropic)
  expect_true(exact$inner_converged)
  expect_false(entropic$inner_converged)
  expect_false(entropic$converged)
  expect_identical(entropic$status, "inner_failure")
})

test_that("unbalanced families expose final and maximum inner residuals", {
  d <- .nested_partial_fixture()
  uot <- ot_sinkhorn_unbalanced(
    d$M, d$p, d$q, epsilon = 0.1, max_iter = 1L, tol = 1e-12
  )
  fugw <- fugw_kl(
    d$C1, d$C2, d$p, d$q,
    epsilon = 0.05, max_iter = 2L, tol = 1e-12,
    max_iter_ot = 1L, tol_ot = 1e-12
  )
  set.seed(42)
  X <- matrix(0.1 * rnorm(12), 4, 3)
  Y <- matrix(0.1 * rnorm(15), 5, 3)
  ucoot <- unbalanced_co_optimal_transport(
    X, Y, max_iter = 2L, tol = 1e-12,
    max_iter_ot = 1L, tol_ot = 1e-12, log = TRUE
  )

  for (out in list(uot, fugw, ucoot)) {
    .expect_nested_diagnostics(out)
    expect_false(out$inner_converged)
    expect_false(out$converged)
  }
})
