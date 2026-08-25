independent_penalized_direction_lp <- function(M, p, q, discard_penalty) {
  ns <- nrow(M)
  nt <- ncol(M)
  nvar <- ns * nt
  constraint <- matrix(0, ns + nt, nvar)
  index <- matrix(seq_len(nvar), ns, nt, byrow = TRUE)
  for (i in seq_len(ns)) constraint[i, index[i, ]] <- 1
  for (j in seq_len(nt)) constraint[ns + j, index[, j]] <- 1
  linear <- lpSolve::lp(
    direction = "min",
    objective.in = as.vector(t(M - 2 * discard_penalty)),
    const.mat = constraint,
    const.dir = rep("<=", ns + nt),
    const.rhs = c(p, q)
  )
  stopifnot(linear$status == 0L)
  plan <- matrix(linear$solution, ns, nt, byrow = TRUE)
  list(plan = plan, mass = sum(plan))
}

test_that("penalized partial FGW gradient has every advertised term", {
  M <- matrix(c(0.2, 0.7, 0.4, 0.1, 0.8, 0.3), 2L, 3L)
  C1 <- matrix(c(0, 1, 1, 0), 2L)
  C2_symmetric <- matrix(c(0, 2, 1, 2, 0, 3, 1, 3, 0), 3L, 3L, byrow = TRUE)
  C2_asymmetric <- matrix(c(0, 2, 1, 1, 0, 3, 4, 2, 0), 3L, 3L, byrow = TRUE)
  p <- c(0.8, 0.6)
  q <- c(0.3, 0.4, 0.5)
  G <- matrix(c(0.12, 0.08, 0.05, 0.05, 0.13, 0.09), 2L, 3L)
  direction <- matrix(
    c(0.01, -0.02, 0.015, -0.005, 0.007, -0.006), 2L, 3L
  )
  alpha <- 0.4
  penalty <- 0.7
  h <- 1e-6

  for (case in list(
    list(C2 = C2_symmetric, symmetric = TRUE),
    list(C2 = C2_asymmetric, symmetric = FALSE)
  )) {
    terms <- .penalized_partial_fgw_terms(
      M, C1, case$C2, G, p, q, alpha, penalty, case$symmetric
    )
    at <- function(step) .penalized_partial_fgw_terms(
      M, C1, case$C2, G + step * direction, p, q,
      alpha, penalty, case$symmetric
    )
    plus <- at(h)
    minus <- at(-h)
    finite_difference <- function(field) {
      (plus[[field]] - minus[[field]]) / (2 * h)
    }

    expect_equal(
      finite_difference("feature_term"),
      (1 - alpha) * sum(M * direction), tolerance = 1e-6
    )
    expect_equal(
      finite_difference("structure_term"),
      alpha * sum(terms$gw_gradient * direction), tolerance = 1e-6
    )
    expect_equal(
      finite_difference("penalty_term"),
      -2 * penalty * sum(direction), tolerance = 1e-6
    )
    expect_equal(
      finite_difference("objective"),
      sum(terms$gradient * direction), tolerance = 1e-6
    )
  }
})

test_that("line-search polynomial and selected step match direct evaluation", {
  M <- matrix(c(0.1, 1.2, 0.7, 0.2), 2L)
  C1 <- matrix(c(0, 1, 1, 0), 2L)
  C2 <- matrix(c(0, 1.5, 1.5, 0), 2L)
  p <- c(0.6, 0.4)
  q <- c(0.45, 0.55)
  G <- matrix(c(0.08, 0.07, 0.04, 0.06), 2L)
  alpha <- 0.55
  penalty <- 0.35
  terms <- .penalized_partial_fgw_terms(
    M, C1, C2, G, p, q, alpha, penalty, TRUE
  )
  linear_cost <- (1 - alpha) * M + alpha * terms$gw_gradient
  direction <- .solve_penalized_partial_fgw_direction(
    linear_cost, p, q, penalty, 20000L, 1e-12
  )
  expect_true(direction$converged)
  line <- .penalized_partial_fgw_line_search(
    M, C1, C2, G, direction$plan, alpha, penalty, TRUE,
    current_terms = terms
  )
  direct <- function(step) .penalized_partial_fgw_terms(
    M, C1, C2, G + step * (direction$plan - G), p, q,
    alpha, penalty, TRUE
  )$objective

  for (step in c(0, 0.17, 0.5, line$step, 1)) {
    expect_equal(
      direct(step),
      terms$objective + line$quadratic * step^2 + line$linear * step,
      tolerance = 1e-10
    )
  }
  grid <- seq(0, 1, length.out = 10001L)
  expect_lte(direct(line$step), min(vapply(grid, direct, numeric(1))) + 1e-9)
})

test_that("one-cell result matches an independent continuous optimizer", {
  M <- matrix(0.2, 1L, 1L)
  C1 <- matrix(0, 1L, 1L)
  C2 <- matrix(1, 1L, 1L)
  p <- 1
  q <- 0.8
  alpha <- 0.6
  penalty <- 0.3
  direct <- function(mass) {
    (1 - alpha) * 0.2 * mass + alpha * mass^2 +
      penalty * (sum(p) + sum(q) - 2 * mass)
  }
  oracle <- stats::optimize(direct, c(0, min(p, q)), tol = 1e-13)
  out <- penalized_partial_fused_gromov_wasserstein(
    M, C1, C2, penalty, p, q, alpha = alpha, tol = 1e-10
  )

  expect_true(out$converged)
  expect_equal(out$transported_mass, oracle$minimum, tolerance = 1e-7)
  expect_equal(rfugw_value(out), oracle$objective, tolerance = 1e-9)
  expect_lte(out$frank_wolfe_gap, out$frank_wolfe_gap_tolerance)
  expect_true(out$stationarity_consistent)
})

test_that("penalty direction is derived rather than reversed or fixed", {
  M <- matrix(0.2, 1L, 1L)
  C1 <- matrix(0, 1L, 1L)
  C2 <- matrix(1, 1L, 1L)
  penalties <- c(0, 0.1, 0.2, 0.4, 1)
  masses <- vapply(penalties, function(penalty) {
    penalized_partial_fused_gromov_wasserstein(
      M, C1, C2, penalty, p = 1, q = 0.8,
      alpha = 0.6, tol = 1e-10
    )$transported_mass
  }, numeric(1))
  expect_true(all(diff(masses) >= -1e-9))
  expect_gt(max(masses) - min(masses), 0.5)
})

test_that("debug trace retains only feasible variable-mass iterates", {
  M <- matrix(c(0.1, 2, 0.8, 0.4, 1.5, 0.2), 2L, 3L)
  C1 <- matrix(c(0, 1, 1, 0), 2L)
  C2 <- matrix(c(0, 2, 1, 2, 0, 1, 1, 1, 0), 3L, 3L, byrow = TRUE)
  p <- c(0.4, 0.9)
  q <- c(0.2, 0.5, 0.4)
  G0 <- matrix(c(0.05, 0, 0.02, 0, 0.04, 0.03), 2L, 3L)
  out <- penalized_partial_fused_gromov_wasserstein(
    M, C1, C2, 0.5, p, q, G0 = G0, trace = TRUE,
    numItermax = 200L, tol = 1e-9
  )

  expect_true(out$converged)
  expect_true(out$trace_feasible)
  expect_equal(out$mass_trace[[1]], sum(G0))
  expect_gt(diff(range(out$mass_trace)), 0.1)
  for (entry in out$debug_trace) {
    expect_true(entry$feasibility$feasible)
    expect_lte(max(rowSums(entry$plan) - p), entry$feasibility$tolerance)
    expect_lte(max(colSums(entry$plan) - q), entry$feasibility$tolerance)
  }
  expect_true(out$inner_converged)
  expect_identical(out$inner_status, "converged")
  expect_true(all(vapply(
    out$direction_trace, function(x) isTRUE(x$converged), logical(1)
  )))
})

test_that("zero terms, balanced limit, and fixed-mass comparison agree", {
  p <- c(0.5, 0.5)
  q <- c(0.5, 0.5)
  M <- matrix(c(0.1, 2, 2, 0.2), 2L)
  zero_C <- matrix(0, 2L, 2L)
  high <- penalized_partial_fused_gromov_wasserstein(
    M, zero_C, zero_C, 10, p, q, alpha = 0.4, tol = 1e-10
  )
  fixed <- partial_fused_gromov_wasserstein(
    M, zero_C, zero_C, p, q, m = 1, alpha = 0.4,
    numItermax = 100L, tol = 1e-10, log = TRUE
  )
  zero_feature <- penalized_partial_fused_gromov_wasserstein(
    matrix(0, 2L, 2L),
    matrix(c(0, 1, 1, 0), 2L),
    matrix(c(0, 2, 2, 0), 2L),
    1, p, q, alpha = 1, tol = 1e-10
  )

  expect_true(high$converged)
  expect_equal(high$transported_mass, 1, tolerance = 1e-10)
  expect_equal(high$plan, fixed$plan, tolerance = 1e-10)
  expect_equal(
    high$feature_term + high$structure_term,
    fixed$partial_fgw_dist, tolerance = 1e-10
  )
  expect_true(zero_feature$converged)
  expect_equal(zero_feature$feature_term, 0)
  expect_true(is.finite(zero_feature$structure_term))
})

test_that("permutation and compatible cost scaling laws hold", {
  M <- matrix(c(0.2, 1.4, 0.7, 1.1, 0.1, 0.9), 2L, 3L)
  C1 <- matrix(c(0, 1.2, 1.2, 0), 2L)
  C2 <- matrix(c(0, 1, 2, 1, 0, 1.5, 2, 1.5, 0), 3L, 3L, byrow = TRUE)
  p <- c(0.6, 0.7)
  q <- c(0.2, 0.5, 0.4)
  penalty <- 0.8
  base <- penalized_partial_fused_gromov_wasserstein(
    M, C1, C2, penalty, p, q, alpha = 0.45, tol = 1e-9
  )
  ip <- c(2, 1)
  jp <- c(3, 1, 2)
  permuted <- penalized_partial_fused_gromov_wasserstein(
    M[ip, jp], C1[ip, ip], C2[jp, jp], penalty,
    p[ip], q[jp], alpha = 0.45, tol = 1e-9
  )
  scale <- 5
  scaled <- penalized_partial_fused_gromov_wasserstein(
    scale * M, sqrt(scale) * C1, sqrt(scale) * C2,
    scale * penalty, p, q, alpha = 0.45, tol = 5e-9
  )

  expect_true(base$converged)
  expect_true(permuted$converged)
  expect_true(scaled$converged)
  expect_equal(permuted$plan, base$plan[ip, jp], tolerance = 1e-8)
  expect_equal(rfugw_value(permuted), rfugw_value(base), tolerance = 1e-8)
  expect_equal(scaled$plan, base$plan, tolerance = 1e-8)
  expect_equal(rfugw_value(scaled), scale * rfugw_value(base), tolerance = 1e-7)
})

test_that("signed linearized costs retain the exact variable-mass oracle", {
  skip_if_not_installed("lpSolve")
  cost <- matrix(c(-1.2, 0.4, 0.7, -0.3, 1.1, 0.2), 2L, 3L)
  p <- c(0.4, 0.8)
  q <- c(0.3, 0.2, 0.6)
  penalty <- 0.25
  shifted <- .solve_penalized_partial_fgw_direction(
    cost, p, q, penalty, 20000L, 1e-12
  )
  oracle <- independent_penalized_direction_lp(cost, p, q, penalty)

  expect_true(shifted$converged)
  expect_gt(shifted$linear_cost_shift, 0)
  expect_equal(shifted$plan, oracle$plan, tolerance = 1e-9)
  expect_equal(shifted$transported_mass, oracle$mass, tolerance = 1e-9)
})

test_that("objective decomposition, status, and validation fail closed", {
  C <- matrix(c(0, 1, 1, 0), 2L)
  M <- matrix(c(0.1, 2, 2, 0.2), 2L)
  out <- penalized_partial_fused_gromov_wasserstein(
    M, C, C, 0.5, p = c(0.4, 0.8), q = c(0.7, 0.2),
    trace = TRUE, tol = 1e-9
  )

  expect_true(out$converged)
  expect_equal(
    rfugw_value(out),
    out$feature_term + out$structure_term + out$discard_penalty_term,
    tolerance = 1e-10
  )
  expect_true(out$objective_consistent)
  expect_true(out$objective_components_consistent)
  expect_lte(out$line_search_residual, out$line_search_tolerance)
  expect_equal(out$discarded_source_mass, 1.2 - out$transported_mass)
  expect_equal(out$discarded_target_mass, 0.9 - out$transported_mass)
  expect_identical(out$measure_normalization, "finite_measure_preserved")
  expect_identical(out$formulation, "penalized_partial_fgw_square")

  expect_error(
    penalized_partial_fused_gromov_wasserstein(M - 3, C, C, 1),
    "M.*nonnegative"
  )
  expect_error(
    penalized_partial_fused_gromov_wasserstein(M, C, C, -1),
    "discard_penalty.*nonnegative"
  )
  expect_error(
    penalized_partial_fused_gromov_wasserstein(M, C, C, 1, alpha = 2),
    "alpha"
  )
  expect_error(
    penalized_partial_fused_gromov_wasserstein(
      M, C, C, 1, p = c(0.5, 0.5), q = c(0.5, 0.5),
      G0 = matrix(1, 2L, 2L)
    ),
    "rowSums"
  )
})

test_that("manifoldalign fixture accepts only a certified result", {
  fixture <- new.env(parent = globalenv())
  sys.source(
    trust_test_resource(
      "extdata", "clients", "manifoldalign_penalized_fgw_fixture.R"
    ),
    envir = fixture
  )
  C <- matrix(c(0, 1, 1, 0), 2L)
  M <- matrix(c(0.1, 2, 2, 0.2), 2L)
  aligned <- fixture$manifoldalign_penalized_fgw_step(M, C, C, 0.5)
  expect_s3_class(aligned$result, "rfugw_result")
  expect_true(aligned$result$converged)
  expect_equal(aligned$coupling, rfugw_plan(aligned$result))

  rejected <- aligned$result
  rejected$converged <- FALSE
  rejected$status <- "max_iter"
  expect_error(
    fixture$manifoldalign_accept_penalized_fgw_result(rejected),
    "did not certify"
  )
})
