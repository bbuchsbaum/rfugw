barycenter_fixture <- function() {
  support <- c(-1, 0.2, 1.5)
  source1 <- c(-1.2, -0.1, 1.1)
  source2 <- c(-0.8, 0.5, 1.8, 2.2)
  list(
    support = support,
    costs = list(
      outer(source1, support, function(x, y) (x - y)^2),
      outer(source2, support, function(x, y) (x - y)^2)
    ),
    support_cost = outer(support, support, function(x, y) (x - y)^2),
    measures = list(c(0.2, 0.5, 0.3), c(0.1, 0.25, 0.45, 0.2)),
    coefficients = c(0.35, 0.65)
  )
}

independent_joint_barycenter_lp <- function(costs, measures, coefficients) {
  stopifnot(requireNamespace("lpSolve", quietly = TRUE))
  S <- length(costs)
  m <- ncol(costs[[1L]])
  plan_sizes <- vapply(costs, length, integer(1))
  offsets <- c(0L, cumsum(plan_sizes))
  q_index <- sum(plan_sizes) + seq_len(m)
  nvar <- sum(plan_sizes) + m
  rows <- list()
  rhs <- numeric()
  for (s in seq_len(S)) {
    n <- nrow(costs[[s]])
    for (i in seq_len(n)) {
      row <- numeric(nvar)
      row[offsets[[s]] + i + (seq_len(m) - 1L) * n] <- 1
      rows[[length(rows) + 1L]] <- row
      rhs <- c(rhs, measures[[s]][[i]])
    }
    for (j in seq_len(m)) {
      row <- numeric(nvar)
      row[offsets[[s]] + (j - 1L) * n + seq_len(n)] <- 1
      row[q_index[[j]]] <- -1
      rows[[length(rows) + 1L]] <- row
      rhs <- c(rhs, 0)
    }
  }
  objective <- c(
    unlist(Map(
      function(cost, coefficient) coefficient * as.vector(cost),
      costs, coefficients
    ), use.names = FALSE),
    numeric(m)
  )
  fit <- lpSolve::lp(
    "min", objective, do.call(rbind, rows), rep("=", length(rows)), rhs
  )
  list(
    status = fit$status,
    objective = fit$objval,
    weights = fit$solution[q_index]
  )
}

test_that("exact fixed-support barycenter matches an independent joint LP", {
  skip_if_not_installed("lpSolve")
  d <- barycenter_fixture()
  fit <- ot_barycenter_weights(
    d$costs, d$measures, d$coefficients,
    mode = "exact", tol = 1e-11
  )
  oracle <- independent_joint_barycenter_lp(
    d$costs, d$measures, d$coefficients
  )

  expect_s3_class(fit, "rfugw_barycenter_result")
  expect_identical(oracle$status, 0L)
  expect_true(fit$converged)
  expect_identical(rfugw_status(fit), "converged")
  expect_equal(rfugw_value(fit), oracle$objective, tolerance = 1e-10)
  expect_equal(fit$weights, oracle$weights, tolerance = 1e-10)
  expect_equal(sum(fit$weights), 1, tolerance = 1e-12)
  expect_true(fit$feasible)
  expect_true(fit$objective_consistent)
  expect_true(fit$kkt_consistent)
  expect_lte(fit$kkt_residual, fit$kkt_tolerance)
  expect_lte(abs(fit$duality_gap), fit$objective_tolerance)
  expect_true(all(fit$component_status == "converged"))
  expect_equal(
    fit$barycenter_objective,
    sum(fit$coefficients * vapply(
      fit$component_results, rfugw_value, numeric(1)
    )),
    tolerance = 1e-12
  )
})

test_that("identical measures are fixed points and support permutations commute", {
  skip_if_not_installed("lpSolve")
  support <- c(-1, 0, 2, 3)
  cost <- outer(support, support, function(x, y) (x - y)^2)
  p <- c(0.1, 0.25, 0.4, 0.25)
  exact <- ot_barycenter_weights(
    list(cost, cost, cost), list(p, p, p),
    coefficients = c(0.01, 0.09, 0.9), mode = "exact", tol = 1e-11
  )
  expect_true(exact$converged)
  expect_equal(exact$weights, p, tolerance = 1e-11)
  expect_equal(exact$barycenter_objective, 0, tolerance = 1e-12)

  d <- barycenter_fixture()
  base <- ot_barycenter_weights(
    d$costs, d$measures, d$coefficients, mode = "exact", tol = 1e-11
  )
  support_permutation <- c(3L, 1L, 2L)
  permuted <- ot_barycenter_weights(
    lapply(d$costs, function(x) x[, support_permutation, drop = FALSE]),
    d$measures,
    d$coefficients,
    mode = "exact",
    tol = 1e-11
  )
  restored <- numeric(3L)
  restored[support_permutation] <- permuted$weights
  expect_equal(restored, base$weights, tolerance = 1e-10)
  expect_equal(permuted$barycenter_objective, base$barycenter_objective,
               tolerance = 1e-12)

  reordered <- ot_barycenter_weights(
    rev(d$costs), rev(d$measures), rev(d$coefficients),
    mode = "exact", tol = 1e-11
  )
  expect_equal(reordered$weights, base$weights, tolerance = 1e-10)
  expect_equal(reordered$barycenter_objective, base$barycenter_objective,
               tolerance = 1e-12)
})

test_that("exact mode handles zero weights, duplicate support, and coefficients", {
  skip_if_not_installed("lpSolve")
  support <- c(0, 1, 1, 3)
  source1 <- c(0, 1, 3)
  source2 <- c(0, 2, 3)
  costs <- list(
    outer(source1, support, function(x, y) abs(x - y)),
    outer(source2, support, function(x, y) abs(x - y))
  )
  measures <- list(c(0.4, 0, 0.6), c(0, 0.2, 0.8))
  fit <- ot_barycenter_weights(
    costs, measures, coefficients = c(1 - 1e-10, 1e-10),
    mode = "exact", tol = 1e-10
  )
  expect_true(fit$converged)
  expect_true(all(is.finite(fit$weights)))
  expect_true(all(fit$weights >= 0))
  expect_equal(sum(fit$weights), 1, tolerance = 1e-12)
  expect_true(any(fit$weights == 0))
  expect_true(all(vapply(fit$component_results, function(x) x$converged,
                         logical(1))))

  zero_coefficient <- ot_barycenter_weights(
    costs, measures, coefficients = c(1, 0), mode = "exact", tol = 1e-10
  )
  expect_true(zero_coefficient$converged)
  expect_equal(
    zero_coefficient$barycenter_objective,
    rfugw_value(zero_coefficient$component_results[[1L]]),
    tolerance = 1e-12
  )
})

test_that("regularized dual gradient agrees with finite differences", {
  d <- barycenter_fixture()
  scale <- max(c(unlist(d$costs), d$support_cost))
  costs <- lapply(d$costs, `/`, scale)
  support_cost <- d$support_cost / scale
  problem <- rfugw:::.prepare_barycenter_problem(
    costs, d$measures, d$coefficients
  )
  q <- c(0.22, 0.51, 0.27)
  epsilon <- 0.15
  evaluated <- rfugw:::.evaluate_regularized_barycenter(
    problem, support_cost, q, epsilon, "log", 10000L, 1e-10,
    list(cross = vector("list", 2L), self = NULL)
  )
  expect_true(evaluated$certified)
  objective <- function(weights) {
    rfugw:::.evaluate_regularized_barycenter(
      problem, support_cost, weights, epsilon, "log", 10000L, 1e-10,
      list(cross = vector("list", 2L), self = NULL)
    )$objective
  }
  for (direction in list(c(1, -1, 0), c(0, 1, -1))) {
    h <- 1e-6
    numerical <- (objective(q + h * direction) -
      objective(q - h * direction)) / (2 * h)
    analytic <- sum(evaluated$gradient * direction)
    expect_equal(analytic, numerical, tolerance = 3e-6)
  }
})

test_that("regularized mode has safeguarded descent and certified KKT", {
  support <- c(-1, 0, 1)
  cost <- outer(support, support, function(x, y) (x - y)^2)
  cost <- cost / max(cost)
  p <- c(0.2, 0.5, 0.3)
  fit <- ot_barycenter_weights(
    list(cost, cost), list(p, p),
    mode = "regularized", support_cost = cost,
    epsilon = 0.1, init_weights = c(0, 1, 0), min_weight = 1e-10,
    max_iter = 200L, tol = 3e-6, sinkhorn_method = "auto",
    sinkhorn_max_iter = 10000L, sinkhorn_tol = 1e-9
  )

  expect_true(fit$converged)
  expect_identical(fit$status, "converged")
  expect_equal(fit$weights, p, tolerance = 1e-5)
  expect_equal(sum(fit$weights), 1, tolerance = 1e-12)
  expect_true(all(fit$weights >= fit$min_weight))
  expect_true(fit$objective_monotone)
  expect_true(all(diff(fit$objective_trace$objective) <= 1e-10))
  expect_lte(fit$kkt_residual, fit$kkt_tolerance)
  expect_lte(fit$simplex_residual, fit$feasibility_tolerance)
  expect_true(all(fit$component_status == "converged"))
  expect_true(fit$support_self_result$converged)
  expect_gt(fit$component_solves, length(fit$component_results))
  expect_gt(fit$warm_start_reuse_count, 0)
  expect_match(fit$objective_convention, "source_self_constants_omitted")
})

test_that("regularized solution agrees with an independent convex optimizer", {
  d <- barycenter_fixture()
  scale <- max(c(unlist(d$costs), d$support_cost))
  costs <- lapply(d$costs, `/`, scale)
  support_cost <- d$support_cost / scale
  epsilon <- 0.2
  fit <- ot_barycenter_weights(
    costs, d$measures, d$coefficients,
    mode = "regularized", support_cost = support_cost,
    epsilon = epsilon, max_iter = 300L, tol = 2e-6,
    sinkhorn_method = "log", sinkhorn_max_iter = 10000L,
    sinkhorn_tol = 1e-9
  )

  softmax <- function(z) {
    value <- exp(c(z, 0) - max(c(z, 0)))
    value / sum(value)
  }
  oracle_objective <- function(z) {
    q <- softmax(z)
    cross <- Map(function(cost, measure) {
      ot_sinkhorn(
        cost, measure, q, epsilon = epsilon, method = "log",
        max_iter = 10000L, tol = 1e-9
      )$regularized_objective
    }, costs, d$measures)
    self <- ot_sinkhorn(
      support_cost, q, q, epsilon = epsilon, method = "log",
      max_iter = 10000L, tol = 1e-9
    )$regularized_objective
    sum(d$coefficients * unlist(cross)) - 0.5 * self
  }
  oracle <- stats::optim(
    c(0, 0), oracle_objective, method = "Nelder-Mead",
    control = list(reltol = 1e-10, maxit = 1000L)
  )
  oracle_weights <- softmax(oracle$par)

  expect_true(fit$converged)
  expect_identical(oracle$convergence, 0L)
  expect_equal(fit$weights, oracle_weights, tolerance = 3e-4)
  expect_equal(fit$barycenter_objective, oracle$value, tolerance = 2e-7)
})

test_that("regularized warm state is deterministic and problem-bound", {
  support <- c(0, 0.4, 1)
  cost <- outer(support, support, function(x, y) (x - y)^2)
  p1 <- c(0.2, 0.3, 0.5)
  p2 <- c(0.4, 0.4, 0.2)
  args <- list(
    costs = list(cost, cost), measures = list(p1, p2),
    mode = "regularized", support_cost = cost,
    epsilon = 0.2, max_iter = 200L, tol = 2e-6,
    sinkhorn_method = "auto", sinkhorn_max_iter = 10000L,
    sinkhorn_tol = 1e-9
  )
  first <- do.call(ot_barycenter_weights, args)
  warm <- do.call(ot_barycenter_weights, c(args, list(init_state = first)))
  expect_true(first$converged)
  expect_true(warm$converged)
  expect_true(warm$warm_started)
  expect_true(warm$warm_start_accepted)
  expect_identical(warm$iterations, 0L)
  expect_equal(warm$weights, first$weights, tolerance = 1e-12)
  expect_lte(warm$total_component_iterations, first$total_component_iterations)

  bad <- first$warm_state
  bad$epsilon <- 2 * bad$epsilon
  expect_error(
    do.call(ot_barycenter_weights, c(args, list(init_state = bad))),
    "does not match.*epsilon"
  )
  expect_error(
    ot_barycenter_weights(
      list(cost), list(p1), mode = "exact", init_state = first
    ),
    "does not accept iterative"
  )
})

test_that("regularized support permutations commute", {
  d <- barycenter_fixture()
  scale <- max(c(unlist(d$costs), d$support_cost))
  costs <- lapply(d$costs, `/`, scale)
  support_cost <- d$support_cost / scale
  args <- list(
    measures = d$measures, coefficients = d$coefficients,
    mode = "regularized", epsilon = 0.2, max_iter = 300L, tol = 3e-6,
    sinkhorn_method = "log", sinkhorn_max_iter = 10000L,
    sinkhorn_tol = 1e-9
  )
  base <- do.call(ot_barycenter_weights, c(list(
    costs = costs, support_cost = support_cost
  ), args))
  permutation <- c(3L, 1L, 2L)
  permuted <- do.call(ot_barycenter_weights, c(list(
    costs = lapply(costs, function(x) x[, permutation, drop = FALSE]),
    support_cost = support_cost[permutation, permutation]
  ), args))
  restored <- numeric(3L)
  restored[permutation] <- permuted$weights
  expect_true(base$converged)
  expect_true(permuted$converged)
  expect_equal(restored, base$weights, tolerance = 2e-5)
  expect_equal(permuted$barycenter_objective, base$barycenter_objective,
               tolerance = 2e-7)
})

test_that("modes reject ambiguous controls and tiny entropy fails honestly", {
  d <- barycenter_fixture()
  expect_error(ot_barycenter_weights(matrix(1, 2, 2), d$measures),
               "nonempty list")
  expect_error(ot_barycenter_weights(d$costs, d$measures[1]),
               "one entry")
  expect_error(ot_barycenter_weights(
    d$costs, d$measures, coefficients = c(0, 0)
  ), "positive total")
  expect_error(ot_barycenter_weights(
    d$costs, d$measures, mode = "exact", support_cost = d$support_cost
  ), "only in")
  expect_error(ot_barycenter_weights(
    d$costs, d$measures, mode = "exact", init_weights = rep(1 / 3, 3)
  ), "does not use")
  expect_error(ot_barycenter_weights(
    d$costs, d$measures, mode = "regularized"
  ), "support_cost.*required")
  expect_error(ot_barycenter_weights(
    d$costs, d$measures, mode = "regularized",
    support_cost = d$support_cost[1:2, 1:2]
  ), "common target-support")
  asymmetric <- d$support_cost
  asymmetric[[2L]] <- asymmetric[[2L]] + 0.1
  expect_error(ot_barycenter_weights(
    d$costs, d$measures, mode = "regularized", support_cost = asymmetric
  ), "symmetric")

  scale <- max(c(unlist(d$costs), d$support_cost))
  tiny <- ot_barycenter_weights(
    lapply(d$costs, `/`, scale), d$measures, d$coefficients,
    mode = "regularized", support_cost = d$support_cost / scale,
    epsilon = 1e-6, max_iter = 2L, tol = 1e-12,
    sinkhorn_method = "auto", sinkhorn_max_iter = 2L,
    sinkhorn_tol = 1e-12
  )
  expect_false(tiny$converged)
  expect_true(tiny$status %in% c(
    "inner_failure", "max_iter", "line_search_failure"
  ))
  expect_false(is.null(tiny$warning_payload))
})

test_that("exact and regularized modes report distinct estimands", {
  skip_if_not_installed("lpSolve")
  d <- barycenter_fixture()
  scale <- max(c(unlist(d$costs), d$support_cost))
  costs <- lapply(d$costs, `/`, scale)
  support_cost <- d$support_cost / scale
  exact <- ot_barycenter_weights(
    costs, d$measures, d$coefficients, mode = "exact"
  )
  regularized <- ot_barycenter_weights(
    costs, d$measures, d$coefficients,
    mode = "regularized", support_cost = support_cost,
    epsilon = 0.2, max_iter = 300L, tol = 3e-6,
    sinkhorn_max_iter = 10000L, sinkhorn_tol = 1e-9
  )
  expect_true(exact$converged)
  expect_true(regularized$converged)
  expect_false(identical(exact$formulation, regularized$formulation))
  expect_false(identical(exact$objective_convention,
                         regularized$objective_convention))
  expect_identical(exact$regularization, 0)
  expect_identical(regularized$regularization, 0.2)
  expect_false(isTRUE(all.equal(
    exact$barycenter_objective, regularized$barycenter_objective,
    tolerance = 1e-8
  )))
})

test_that("barycenter benchmark retains quality and continuation evidence", {
  baseline <- utils::read.csv(
    bench_test_resource("wasserstein-barycenter-baseline.csv"),
    stringsAsFactors = FALSE
  )
  required <- c(
    "mode", "status", "converged", "objective", "simplex_residual",
    "kkt_residual", "kkt_tolerance", "outer_iterations",
    "component_solves", "component_iterations", "warm_reuse_count",
    "elapsed_seconds", "allocated_bytes"
  )
  expect_true(all(required %in% names(baseline)))
  expect_true(all(baseline$converged))
  expect_true(all(baseline$status == "converged"))
  expect_true(all(baseline$kkt_residual <= baseline$kkt_tolerance))
  expect_lte(
    baseline$component_iterations[baseline$mode == "regularized_warm"],
    baseline$component_iterations[baseline$mode == "regularized_cold"]
  )
  expect_equal(
    baseline$objective[baseline$mode == "regularized_warm"],
    baseline$objective[baseline$mode == "regularized_cold"],
    tolerance = 1e-12
  )
})
