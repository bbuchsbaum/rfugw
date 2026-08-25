independent_penalized_partial_lp <- function(M, p, q, discard_penalty) {
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
  list(
    plan = plan,
    mass = sum(plan),
    objective = sum(M * plan) + discard_penalty *
      (sum(p) + sum(q) - 2 * sum(plan))
  )
}

test_that("one-cell oracle fixes discard-penalty direction and units", {
  low <- ot_partial_penalized(
    matrix(1, 1, 1), p = 2, q = 3, discard_penalty = 0.2
  )
  high <- ot_partial_penalized(
    matrix(1, 1, 1), p = 2, q = 3, discard_penalty = 1
  )
  expect_equal(low$transported_mass, 0, tolerance = 1e-12)
  expect_equal(rfugw_value(low), 0.2 * 5, tolerance = 1e-12)
  expect_equal(high$transported_mass, 2, tolerance = 1e-12)
  expect_equal(high$transport_term, 2, tolerance = 1e-12)
  expect_equal(high$discard_penalty_term, 1, tolerance = 1e-12)
  expect_equal(rfugw_value(high), 3, tolerance = 1e-12)
  expect_identical(
    high$penalty_direction,
    "larger_discard_penalty_weakly_favors_more_transport"
  )
})

test_that("penalized partial OT matches an independent original-variable LP", {
  skip_if_not_installed("lpSolve")
  M <- matrix(c(0.2, 2, 3, 1.5, 0.1, 2.5), 2L, 3L, byrow = TRUE)
  p <- c(0.3, 0.9)
  q <- c(0.2, 0.4, 0.7)
  for (penalty in c(0, 0.4, 1.2, 10)) {
    expected <- independent_penalized_partial_lp(M, p, q, penalty)
    actual <- ot_partial_penalized(M, p, q, discard_penalty = penalty)
    expect_true(actual$converged)
    expect_equal(rfugw_value(actual), expected$objective, tolerance = 1e-9)
    expect_equal(actual$transported_mass, expected$mass, tolerance = 1e-9)
    expect_lte(max(rowSums(actual$plan) - p), actual$feasibility_tolerance)
    expect_lte(max(colSums(actual$plan) - q), actual$feasibility_tolerance)
    expect_lte(abs(actual$duality_gap), actual$duality_gap_tolerance)
  }
})

test_that("zero, large, and all-zero costs have derived mass behavior", {
  positive <- matrix(c(1, 3, 2, 4), 2L)
  zero_penalty <- ot_partial_penalized(positive, discard_penalty = 0)
  large_penalty <- ot_partial_penalized(positive, discard_penalty = 10)
  zero_cost <- ot_partial_penalized(
    matrix(0, 2L, 3L), p = c(0.25, 0.75), q = c(0.2, 0.3, 0.8),
    discard_penalty = 0.1
  )
  expect_equal(zero_penalty$transported_mass, 0, tolerance = 1e-12)
  expect_equal(large_penalty$transported_mass, 1, tolerance = 1e-12)
  expect_equal(zero_cost$transported_mass, 1, tolerance = 1e-12)
  expect_equal(zero_cost$discarded_source_mass, 0, tolerance = 1e-12)
  expect_equal(zero_cost$discarded_target_mass, 0.3, tolerance = 1e-12)
})

test_that("finite unequal masses, zero entries, and rectangles retain provenance", {
  M <- matrix(c(0.1, 5, 2, 4, 0.2, 3), 3L, 2L, byrow = TRUE)
  p <- c(0, 0.4, 1.1)
  q <- c(0.7, 0.2)
  out <- ot_partial_penalized(M, p, q, discard_penalty = 2)
  expect_true(out$converged)
  expect_equal(dim(out$plan), c(3L, 2L))
  expect_equal(out$original_source_mass, 1.5)
  expect_equal(out$original_target_mass, 0.9)
  expect_equal(out$discarded_source_mass, 1.5 - out$transported_mass)
  expect_equal(out$discarded_target_mass, 0.9 - out$transported_mass)
  expect_equal(out$plan[1, ], c(0, 0), tolerance = 1e-12)
  expect_identical(out$measure_normalization, "finite_measure_preserved")
  expect_equal(
    rfugw_value(out),
    out$transport_term + out$discard_penalty_term,
    tolerance = 1e-12
  )
  expect_equal(out$augmented_cost[seq_len(3L), 3L], rep(2, 3L))
  expect_equal(out$augmented_cost[4L, seq_len(2L)], rep(2, 2L))
  expect_equal(out$augmented_cost[4L, 3L], 0)
})

test_that("joint cost and discard-penalty scaling preserves the optimizer", {
  M <- matrix(c(0.2, 2, 3, 0.1), 2L, byrow = TRUE)
  p <- c(0.4, 0.8)
  q <- c(0.5, 0.6)
  base <- ot_partial_penalized(M, p, q, discard_penalty = 0.7)
  scaled <- ot_partial_penalized(5 * M, p, q, discard_penalty = 3.5)
  expect_equal(scaled$plan, base$plan, tolerance = 1e-12)
  expect_equal(scaled$transported_mass, base$transported_mass, tolerance = 1e-12)
  expect_equal(rfugw_value(scaled), 5 * rfugw_value(base), tolerance = 1e-10)
})

test_that("zero measures use the analytic certificate", {
  out <- ot_partial_penalized(
    matrix(c(0, 1, 2, 3), 2L), p = c(0, 0), q = c(0, 0),
    discard_penalty = 4
  )
  expect_true(out$converged)
  expect_identical(out$backend, "analytic_zero_measures")
  expect_equal(out$plan, matrix(0, 2L, 2L))
  expect_equal(rfugw_value(out), 0)
  expect_equal(out$duality_gap, 0)
})

test_that("mutation sentinels reject reversed penalty direction and fixed mass", {
  M <- matrix(c(0.3, 5, 5, 0.4), 2L, byrow = TRUE)
  penalties <- c(0.05, 0.3, 1, 10)
  masses <- vapply(penalties, function(penalty) {
    ot_partial_penalized(M, discard_penalty = penalty)$transported_mass
  }, numeric(1))
  expect_true(all(diff(masses) >= -1e-12))
  expect_gt(max(masses) - min(masses), 0.5)

  reversed_penalty_scores <- vapply(c(0, 1), function(mass) {
    0.3 * mass + 10 * (2 * mass)
  }, numeric(1))
  expect_gt(reversed_penalty_scores[[2]], reversed_penalty_scores[[1]])
  expect_equal(tail(masses, 1), 1, tolerance = 1e-12)
})

test_that("penalized partial validation fails before native computation", {
  M <- matrix(c(0, 1, 1, 0), 2L)
  expect_error(ot_partial_penalized(M, discard_penalty = -1), "nonnegative")
  expect_error(ot_partial_penalized(M, discard_penalty = Inf), "finite")
  expect_error(
    ot_partial_penalized(M - 2, discard_penalty = 1),
    "M.*nonnegative"
  )
  expect_error(
    ot_partial_penalized(M, p = c(1, NA), discard_penalty = 1),
    "p.*finite"
  )
  expect_error(
    ot_partial_penalized(M, q = c(-1, 2), discard_penalty = 1),
    "q.*nonnegative"
  )
})
