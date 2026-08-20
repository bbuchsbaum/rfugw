audit_gw_value <- function(C1, C2, G) {
  value <- 0
  for (i in seq_len(nrow(C1))) {
    for (j in seq_len(nrow(C1))) {
      for (k in seq_len(nrow(C2))) {
        for (l in seq_len(nrow(C2))) {
          value <- value + (C1[i, j] - C2[k, l])^2 * G[i, k] * G[j, l]
        }
      }
    }
  }
  value
}

audit_gw_gradient <- function(C1, C2, G) {
  gradient <- matrix(0, nrow(G), ncol(G))
  for (a in seq_len(nrow(G))) {
    for (b in seq_len(ncol(G))) {
      for (j in seq_len(nrow(G))) {
        for (l in seq_len(ncol(G))) {
          gradient[a, b] <- gradient[a, b] +
            (C1[a, j] - C2[b, l])^2 * G[j, l]
        }
      }
      for (i in seq_len(nrow(G))) {
        for (k in seq_len(ncol(G))) {
          gradient[a, b] <- gradient[a, b] +
            (C1[i, a] - C2[k, b])^2 * G[i, k]
        }
      }
    }
  }
  gradient
}

test_that("general square-loss gradient matches an O(n^4) oracle", {
  set.seed(20260820)
  C1 <- matrix(runif(16, 0.05, 1.4), 4, 4)
  C2 <- matrix(runif(9, 0.1, 1.7), 3, 3)
  diag(C1) <- 0
  diag(C2) <- 0
  G <- matrix(runif(12, 0.05, 0.25), 4, 3)
  direction <- matrix(rnorm(12), 4, 3)
  direction <- direction / sqrt(sum(direction^2))

  native <- rfugw:::cpp_gw_square_terms_square(C1, C2, G, symmetric = FALSE)
  oracle_gradient <- audit_gw_gradient(C1, C2, G)
  h <- 1e-6
  finite_difference <- (
    audit_gw_value(C1, C2, G + h * direction) -
      audit_gw_value(C1, C2, G - h * direction)
  ) / (2 * h)

  expect_equal(native$loss, audit_gw_value(C1, C2, G), tolerance = 1e-12)
  expect_equal(native$grad, oracle_gradient, tolerance = 1e-10)
  expect_equal(sum(native$grad * direction), finite_difference, tolerance = 1e-7)
})

test_that("mixed precision makes requested and effective tolerances observable", {
  set.seed(1)
  n <- 8L
  C1 <- as.matrix(dist(matrix(rnorm(n * 3), n, 3))); C1 <- C1 / max(C1)
  C2 <- as.matrix(dist(matrix(rnorm(n * 3), n, 3))); C2 <- C2 / max(C2)
  M <- matrix(runif(n * n), n, n)
  out <- fgw_entropic(
    M, C1, C2, epsilon = 0.1, max_iter = 500L,
    tol = 1e-9, sinkhorn_tol = 1e-9,
    precision = "mixed", symmetric = TRUE
  )

  expect_true(
    all(c(
      "requested_tol", "effective_tol", "requested_inner_tol",
      "effective_inner_tol", "requested_precision", "effective_precision",
      "termination_reason"
    ) %in% names(out))
  )
  expect_identical(out$requested_tol, 1e-9)
  expect_identical(out$requested_inner_tol, 1e-9)
  expect_true(out$effective_tol >= out$requested_tol)
  expect_true(out$effective_inner_tol >= out$requested_inner_tol)
  expect_identical(out$converged, out$residual <= out$effective_tol)
  expect_false(out$termination_reason == "max_iter" && out$iterations < out$max_iter)
})

test_that("starved exact partial directions propagate inner failure", {
  set.seed(174)
  ns <- 4L
  nt <- 5L
  M <- matrix(runif(ns * nt), ns, nt)
  C1 <- matrix(0, ns, ns)
  C2 <- matrix(0, nt, nt)
  p <- c(0.1, 0.2, 0.3, 0.4)
  q <- c(0.08, 0.12, 0.2, 0.25, 0.35)
  m <- 0.7
  G0 <- (p %o% q) * m

  limited <- rfugw:::cpp_partial_fgw_exact_square(
    M, C1, C2, p, q, m, 0, TRUE, G0,
    1L, 0, 1L, 1L, 1e-12
  )
  full <- rfugw:::cpp_partial_fgw_exact_square(
    M, C1, C2, p, q, m, 0, TRUE, G0,
    1L, 0, 1L, 20000L, 1e-12
  )

  expect_gt(abs(limited$objective - full$objective), 1e-8)
  expect_true(
    all(c("inner_status", "inner_converged", "inner_iterations", "inner_residual") %in%
      names(limited))
  )
  expect_false(limited$inner_converged)
})

test_that("diagnostic helper cannot certify invalid output", {
  out <- rfugw:::.attach_solver_diagnostics(
    out = list(
      plan = matrix(c(-1, 2, 2, -1), 2, 2),
      ot_dist = Inf,
      objective = Inf
    ),
    residual = 0,
    converged = TRUE,
    iterations = 1L,
    max_iter = 10L,
    p = c(0.5, 0.5),
    q = c(0.5, 0.5),
    inner_residual = Inf
  )

  expect_false(out$converged)
  expect_true(out$status %in% c("numerical_failure", "infeasible", "inner_failure"))
})

test_that("ot_kl implements generalized KL including mass and support terms", {
  independent_kl <- function(plan, p, q) {
    reference <- tcrossprod(p, q)
    if (any(plan > 0 & reference == 0)) return(Inf)
    positive <- plan > 0
    sum(plan[positive] * log(plan[positive] / reference[positive])) -
      sum(plan) + sum(reference)
  }

  plan <- matrix(c(0.30, 0.10, 0.00, 0.20), 2, 2)
  p <- c(0.4, 0.6)
  q <- c(0.7, 0.3)
  expect_equal(ot_kl(plan, p, q), independent_kl(plan, p, q), tolerance = 1e-12)
  expect_equal(ot_kl(matrix(0, 2, 2), p, q), 1, tolerance = 1e-12)

  p_zero <- c(1, 0)
  outside <- matrix(c(0.4, 0.1, 0.4, 0.1), 2, 2)
  expect_identical(ot_kl(outside, p_zero, c(0.5, 0.5)), Inf)
  expect_error(ot_kl(plan, c(1.2, -0.2), q), "nonnegative")
  expect_error(ot_kl(plan, c(NA_real_, 1), q), "finite")
})
