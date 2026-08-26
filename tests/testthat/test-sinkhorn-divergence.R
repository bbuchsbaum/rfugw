balanced_2x2_regularized_oracle <- function(M, p, q, epsilon) {
  lower <- max(0, p[[1]] - q[[2]])
  upper <- min(p[[1]], q[[1]])
  objective <- function(t11) {
    plan <- matrix(c(
      t11, p[[1]] - t11,
      q[[1]] - t11, p[[2]] - q[[1]] + t11
    ), 2L, byrow = TRUE)
    sum(M * plan) + epsilon * ot_kl(plan, p, q)
  }
  stats::optimize(objective, c(lower, upper), tol = 1e-14)$objective
}

test_that("balanced Sinkhorn exposes one explicit regularized primal and dual", {
  M <- matrix(c(0, 1.2, 0.7, 0.1), 2L, byrow = TRUE)
  p <- c(0.4, 0.6)
  q <- c(0.5, 0.5)
  epsilon <- 0.3
  scaling <- ot_sinkhorn(
    M, p, q, epsilon = epsilon, method = "scaling",
    max_iter = 5000L, tol = 1e-10
  )
  logarithmic <- ot_sinkhorn(
    M, p, q, epsilon = epsilon, method = "log",
    max_iter = 5000L, tol = 1e-10
  )
  oracle <- balanced_2x2_regularized_oracle(M, p, q, epsilon)

  for (out in list(scaling, logarithmic)) {
    expect_identical(
      out$entropy_convention,
      "epsilon_kl_plan_to_product_measure"
    )
    expect_identical(
      out$entropy_reference,
      "source_weights_tensor_target_weights"
    )
    expect_equal(
      out$regularized_objective,
      out$ot_dist + epsilon * ot_kl(out$plan, p, q),
      tolerance = 1e-12
    )
    expect_equal(out$regularized_objective, oracle, tolerance = 1e-8)
    expect_equal(
      out$regularized_objective,
      out$regularized_dual_objective,
      tolerance = out$regularized_duality_gap_tolerance
    )
    expect_true(out$regularized_dual_consistent)
    expected_offset <- epsilon * (
      1 - sum(p[p > 0] * log(p[p > 0])) -
        sum(q[q > 0] * log(q[q > 0]))
    )
    expect_equal(out$product_reference_offset, expected_offset, tolerance = 1e-9)
  }
  expect_equal(
    scaling$regularized_objective,
    logarithmic$regularized_objective,
    tolerance = 1e-9
  )
  expect_equal(scaling$ot_dist, logarithmic$ot_dist, tolerance = 1e-9)
})

test_that("regularized contract handles zero support and constant offsets", {
  M <- matrix(c(0, 2, 1, 3, 2, 0), 3L, 2L, byrow = TRUE)
  p <- c(0, 0.4, 0.6)
  q <- c(0.75, 0.25)
  out <- ot_sinkhorn(
    M, p, q, epsilon = 0.5, method = "log",
    max_iter = 5000L, tol = 1e-10
  )
  expect_true(out$converged)
  expect_equal(rowSums(out$plan), p, tolerance = 1e-9)
  expect_equal(out$regularized_source_potential[[1]], 0)
  expect_true(is.finite(out$regularized_objective))
  expect_true(out$regularized_dual_consistent)
})

test_that("Sinkhorn divergence is symmetric, nonnegative, and debiased", {
  x <- c(0, 1, 3)
  y <- c(0.5, 2, 4)
  a <- c(0.2, 0.3, 0.5)
  b <- c(0.4, 0.35, 0.25)
  xy <- ot_sinkhorn_divergence(
    x, y, p = 2, source_weights = a, target_weights = b,
    epsilon = 2, method = "log", max_iter = 5000L, tol = 1e-9
  )
  yx <- ot_sinkhorn_divergence(
    y, x, p = 2, source_weights = b, target_weights = a,
    epsilon = 2, method = "log", max_iter = 5000L, tol = 1e-9
  )
  self <- ot_sinkhorn_divergence(
    x, x, p = 2, source_weights = a, target_weights = a,
    epsilon = 2, method = "log", max_iter = 5000L, tol = 1e-9
  )

  expect_true(xy$converged)
  expect_true(xy$nonnegative_within_tolerance)
  expect_gte(rfugw_value(xy), -xy$nonnegative_tolerance)
  expect_equal(rfugw_value(xy), rfugw_value(yx), tolerance = 1e-8)
  expect_equal(rfugw_value(self), 0, tolerance = 1e-10)
  expect_true(all(xy$component_converged))
  expect_setequal(
    names(xy$component_solves),
    c("cross", "source_self", "target_self")
  )
  expect_true(all(xy$component_runtime_seconds >= 0))
  expect_equal(
    rfugw_value(xy),
    xy$component_values[["cross"]] -
      0.5 * xy$component_values[["source_self"]] -
      0.5 * xy$component_values[["target_self"]]
  )
})

test_that("Sinkhorn divergence obeys permutation and joint cost-epsilon scaling", {
  x <- matrix(c(0, 0, 1, 2, 3, 1), 3L, byrow = TRUE)
  y <- matrix(c(1, 0, 2, 2), 2L, byrow = TRUE)
  a <- c(0.2, 0.3, 0.5)
  b <- c(0.6, 0.4)
  base <- ot_sinkhorn_divergence(
    x, y, p = 2, source_weights = a, target_weights = b,
    epsilon = 2, method = "log", max_iter = 5000L, tol = 1e-9
  )
  permuted <- ot_sinkhorn_divergence(
    x[c(3, 1, 2), ], y[c(2, 1), ], p = 2,
    source_weights = a[c(3, 1, 2)], target_weights = b[c(2, 1)],
    epsilon = 2, method = "log", max_iter = 5000L, tol = 1e-9
  )
  scaled <- ot_sinkhorn_divergence(
    sqrt(3) * x, sqrt(3) * y, p = 2,
    source_weights = a, target_weights = b,
    epsilon = 6, method = "log", max_iter = 5000L, tol = 1e-9
  )
  expect_equal(rfugw_value(permuted), rfugw_value(base), tolerance = 1e-8)
  expect_equal(rfugw_value(scaled), 3 * rfugw_value(base), tolerance = 2e-8)
  expect_equal(
    scaled$component_solves$cross$plan,
    base$component_solves$cross$plan,
    tolerance = 1e-8
  )
})

test_that("duplicates, zero weights, tiny epsilon, and constant costs are explicit", {
  duplicate <- ot_sinkhorn_divergence(
    c(0, 0, 1), c(0, 2), p = 1,
    source_weights = c(0, 0.5, 0.5), target_weights = c(0.4, 0.6),
    epsilon = 0.5, method = "log", max_iter = 10000L, tol = 1e-8
  )
  expect_true(duplicate$converged)
  expect_true(is.finite(rfugw_value(duplicate)))

  tiny <- ot_sinkhorn_divergence(
    c(0, 0.1), c(0.02, 0.12), p = 2,
    epsilon = 2e-3, method = "log", max_iter = 20000L, tol = 1e-8
  )
  expect_true(tiny$converged)
  expect_gte(rfugw_value(tiny), -tiny$nonnegative_tolerance)

  constant <- ot_sinkhorn_divergence(
    cost = matrix(1, 2L, 3L),
    source_cost = matrix(0, 2L, 2L),
    target_cost = matrix(0, 3L, 3L),
    cost_power = 1,
    p = 1,
    epsilon = 0.5,
    method = "log",
    max_iter = 5000L,
    tol = 1e-9
  )
  expect_equal(rfugw_value(constant), 1, tolerance = 1e-9)
  expect_false(constant$metric_certified)
})

test_that("component state is reusable without inherited certification", {
  x <- c(0, 1, 3)
  y <- c(0.5, 2, 4)
  cold <- ot_sinkhorn_divergence(
    x, y, p = 2, epsilon = 2, method = "log",
    max_iter = 5000L, tol = 1e-9
  )
  state <- list(
    cross = cold$component_solves$cross,
    source_self = cold$component_solves$source_self,
    target_self = cold$component_solves$target_self
  )
  warm <- ot_sinkhorn_divergence(
    x, y, p = 2, epsilon = 2, method = "scaling",
    max_iter = 5000L, tol = 1e-9, init_state = state
  )
  expect_true(warm$converged)
  expect_equal(rfugw_value(warm), rfugw_value(cold), tolerance = 1e-8)
  expect_true(all(vapply(
    warm$component_solves, function(component) component$warm_started, logical(1)
  )))
  expect_lte(sum(warm$component_iterations), sum(cold$component_iterations))
})

test_that("component failure and invalid self-cost contracts fail closed", {
  expect_error(
    ot_sinkhorn_divergence(
      c(0, 1, 4), c(0.5, 2, 5), p = 2,
      epsilon = 0.01, method = "log", max_iter = 1L, tol = 1e-15
    ),
    "component solve.*failed"
  )
  inspected <- ot_sinkhorn_divergence(
    c(0, 1, 4), c(0.5, 2, 5), p = 2,
    epsilon = 0.01, method = "log", max_iter = 1L, tol = 1e-15,
    allow_uncertified = TRUE
  )
  expect_false(inspected$converged)
  expect_identical(inspected$status, "component_failure")
  expect_false(inspected$value_certified)

  expect_error(
    ot_sinkhorn_divergence(
      cost = matrix(1, 2, 3), source_cost = matrix(0, 2, 2),
      target_cost = NULL, cost_power = 1
    ),
    "requires both"
  )
  expect_error(
    ot_sinkhorn_divergence(
      cost = matrix(1, 2, 2),
      source_cost = matrix(c(0, 1, 2, 0), 2),
      target_cost = matrix(c(0, 1, 1, 0), 2),
      cost_power = 1
    ),
    "symmetric-cost"
  )
  expect_error(
    ot_sinkhorn_divergence(
      cost = matrix(1, 2, 2),
      source_cost = matrix(1, 2, 2),
      target_cost = matrix(1, 2, 2),
      cost_power = 1
    ),
    "zero diagonal"
  )
})
