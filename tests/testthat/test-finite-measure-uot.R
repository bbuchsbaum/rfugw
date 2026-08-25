gkl_vector_reference <- function(x, y) {
  if (any(x > 0 & y == 0)) return(Inf)
  keep <- x > 0
  sum(x[keep] * log(x[keep] / y[keep])) - sum(x) + sum(y)
}

uot_objective_reference <- function(plan, M, p, q, rho, epsilon) {
  transport <- sum(plan * M)
  source_kl <- gkl_vector_reference(rowSums(plan), p)
  target_kl <- gkl_vector_reference(colSums(plan), q)
  product <- tcrossprod(p, q)
  keep <- plan > 0
  plan_kl <- if (any(plan > 0 & product == 0)) {
    Inf
  } else {
    sum(plan[keep] * log(plan[keep] / product[keep])) -
      sum(plan) + sum(product)
  }
  transport + rho[[1L]] * source_kl + rho[[2L]] * target_kl +
    epsilon * plan_kl
}

test_that("finite-measure UOT matches the analytic one-cell oracle", {
  M <- matrix(0.7, 1, 1)
  p <- 2
  q <- 5
  rho <- c(1.3, 2.1)
  epsilon <- 0.4
  expected_plan <- exp(
    ((rho[[1L]] + epsilon) * log(p) +
      (rho[[2L]] + epsilon) * log(q) - M[[1L]]) /
      (sum(rho) + epsilon)
  )

  out <- ot_sinkhorn_unbalanced(
    M, p, q, epsilon = epsilon, rho = rho,
    normalization = "none", method = "log", max_iter = 10000L, tol = 1e-11
  )

  expect_true(out$converged)
  expect_equal(out$plan[[1L]], expected_plan, tolerance = 1e-9)
  expected_objective <- uot_objective_reference(
    out$plan, M, p, q, rho, epsilon
  )
  expect_equal(out$regularized_objective, expected_objective, tolerance = 1e-11)
  expect_equal(out$uot_objective_components$total, expected_objective, tolerance = 1e-11)
  expect_equal(out$original_source_mass, 2)
  expect_equal(out$original_target_mass, 5)
  expect_equal(out$effective_source_mass, 2)
  expect_equal(out$effective_target_mass, 5)
  expect_identical(out$normalization, "none")
  expect_true(out$uot_certificate$converged)
})

test_that("finite-measure UOT obeys its derived common-mass scaling law", {
  M <- matrix(c(0.2, 0.7, 1.1, 0.1, 0.5, 0.9), 2, 3)
  p <- c(0.7, 1.3)
  q <- c(0.4, 0.8, 1.8)
  rho <- c(1.7, 2.4)
  epsilon <- 0.6
  scale <- 7
  exponent <- (sum(rho) + 2 * epsilon) / (sum(rho) + epsilon)

  base <- ot_sinkhorn_unbalanced(
    M, p, q, epsilon = epsilon, rho = rho,
    normalization = "none", method = "log", max_iter = 10000L, tol = 1e-10
  )
  scaled <- ot_sinkhorn_unbalanced(
    M, scale * p, scale * q, epsilon = epsilon, rho = rho,
    normalization = "none", method = "log", max_iter = 10000L, tol = 1e-10
  )

  expect_true(base$converged)
  expect_true(scaled$converged)
  expect_equal(scaled$plan, scale^exponent * base$plan, tolerance = 2e-8)
  expected_objective <- scale^exponent * base$regularized_objective +
    rho[[1L]] * (scale - scale^exponent) * sum(p) +
    rho[[2L]] * (scale - scale^exponent) * sum(q) +
    epsilon * (scale^2 - scale^exponent) * sum(p) * sum(q)
  expect_equal(scaled$regularized_objective, expected_objective, tolerance = 2e-8)
})

test_that("normalization policies are explicit, recorded, and mutation-sensitive", {
  M <- matrix(c(0, 0.4, 0.8, 0.2), 2)
  p <- c(1, 3)
  q <- c(2, 3)
  args <- list(M = M, p = p, q = q, epsilon = 0.3, rho = c(2, 3),
               method = "log", max_iter = 5000L, tol = 1e-9)
  legacy <- do.call(ot_sinkhorn_unbalanced, args)
  separate <- do.call(ot_sinkhorn_unbalanced, c(args, list(normalization = "separate")))
  none <- do.call(ot_sinkhorn_unbalanced, c(args, list(normalization = "none")))
  joint <- do.call(ot_sinkhorn_unbalanced, c(args, list(normalization = "joint")))

  expect_true(legacy$normalization_defaulted)
  expect_identical(legacy$normalization_source, "backward_compatible_0.1_default")
  expect_equal(legacy$plan, separate$plan, tolerance = 1e-12)
  expect_equal(c(separate$effective_source_mass, separate$effective_target_mass), c(1, 1))
  expect_equal(c(none$effective_source_mass, none$effective_target_mass), c(4, 5))
  expect_equal(sum(c(joint$effective_source_mass, joint$effective_target_mass)), 2)
  expect_equal(joint$effective_source_mass / joint$effective_target_mass, 4 / 5)
  expect_false(isTRUE(all.equal(none$plan, separate$plan, tolerance = 1e-6)))
  expect_equal(none$regularized_objective,
               uot_objective_reference(none$plan, M, p, q, c(2, 3), 0.3),
               tolerance = 1e-10)
})

test_that("zero support and zero measures have explicit finite behavior", {
  M <- matrix(c(0, 1, 2, 0, 3, 1), 2, 3)
  supported <- ot_sinkhorn_unbalanced(
    M, c(2, 0), c(0, 3, 1), epsilon = 0.5, rho = c(2, 4),
    normalization = "none", method = "log", max_iter = 5000L, tol = 1e-9
  )
  expect_true(supported$converged)
  expect_equal(supported$plan[2, ], rep(0, 3), tolerance = 0)
  expect_equal(supported$plan[, 1], rep(0, 2), tolerance = 0)
  expect_true(all(is.finite(supported$plan)))

  one_sided <- ot_sinkhorn_unbalanced(
    M, c(0, 0), c(0, 3, 1), epsilon = 0.5, rho = c(2, 4),
    normalization = "none", method = "auto"
  )
  expect_true(one_sided$converged)
  expect_identical(one_sided$termination_reason, "zero_measure_closed_form")
  expect_equal(one_sided$plan, matrix(0, 2, 3), tolerance = 0)
  expect_equal(one_sided$regularized_objective, 4 * 4, tolerance = 0)

  both_zero <- ot_sinkhorn_unbalanced(
    M, c(0, 0), c(0, 0, 0), epsilon = 0.5, rho = c(2, 4),
    normalization = "none"
  )
  expect_equal(both_zero$regularized_objective, 0, tolerance = 0)
  expect_error(
    ot_sinkhorn_unbalanced(M, c(0, 0), c(0, 3, 1), normalization = "separate"),
    "requires positive source and target mass"
  )
})

test_that("scaling and log backends agree for moderate unequal finite measures", {
  M <- matrix(c(0, 0.4, 0.9, 0.1, 0.6, 0.2), 3, 2)
  p <- c(0.2, 0, 1.8)
  q <- c(0.7, 2.3)
  args <- list(M = M, p = p, q = q, epsilon = 0.4, rho = c(2, 5),
               normalization = "none", max_iter = 5000L, tol = 1e-8)
  scaling <- do.call(ot_sinkhorn_unbalanced, c(args, list(method = "scaling")))
  log_domain <- do.call(ot_sinkhorn_unbalanced, c(args, list(method = "log")))

  expect_true(scaling$converged)
  expect_true(log_domain$converged)
  expect_equal(scaling$plan, log_domain$plan, tolerance = 2e-7)
  expect_equal(scaling$regularized_objective,
               log_domain$regularized_objective, tolerance = 2e-7)
})

test_that("log finite-measure UOT remains analytic at an extreme mass ratio", {
  M <- matrix(0.25, 1, 1)
  p <- 1e-12
  q <- 1e12
  rho <- c(2, 5)
  epsilon <- 0.5
  expected <- exp(
    ((rho[[1L]] + epsilon) * log(p) +
      (rho[[2L]] + epsilon) * log(q) - M[[1L]]) /
      (sum(rho) + epsilon)
  )
  out <- ot_sinkhorn_unbalanced(
    M, p, q, epsilon = epsilon, rho = rho,
    normalization = "none", method = "log", max_iter = 10000L, tol = 1e-11
  )
  expect_true(out$converged)
  expect_true(all(is.finite(out$plan)))
  expect_equal(out$plan[[1L]], expected, tolerance = 2e-8)
  expect_equal(c(out$effective_source_mass, out$effective_target_mass), c(p, q))
})

test_that("finite-measure validation fails before solving", {
  M <- matrix(c(0, 1, 1, 0), 2)
  expect_error(ot_sinkhorn_unbalanced(M, c(1, -1), normalization = "none"), "nonnegative")
  expect_error(ot_sinkhorn_unbalanced(M, c(1, Inf), normalization = "none"), "finite")
  expect_error(ot_sinkhorn_unbalanced(M, c(0, 0), c(0, 0), normalization = "joint"),
               "positive combined mass")
  expect_error(
    ot_sinkhorn_unbalanced(M, c(1e308, 1e308), c(1, 1), normalization = "none"),
    "finite total mass"
  )
})
