.attach_certified_example <- function(out, plan, ..., p = NULL, q = NULL) {
  rfugw:::.attach_solver_diagnostics(
    out = out,
    residual = 0,
    converged = TRUE,
    iterations = 2L,
    max_iter = 2L,
    p = p,
    q = q,
    plan = plan,
    ...
  )
}

test_that("balanced convergence requires a finite nonnegative feasible plan", {
  p <- q <- c(0.5, 0.5)
  valid_plan <- diag(p)
  valid <- .attach_certified_example(
    list(plan = valid_plan, ot_dist = 0),
    valid_plan,
    p = p,
    q = q,
    feasibility = "balanced",
    feasibility_tol = 1e-10,
    objective_recomputed = 0
  )
  expect_true(valid$converged)
  expect_true(valid$feasible)
  expect_true(valid$objective_consistent)

  negative_plan <- matrix(c(-0.1, 0.6, 0.6, -0.1), 2, 2)
  negative <- .attach_certified_example(
    list(plan = negative_plan, ot_dist = 0),
    negative_plan,
    p = p,
    q = q,
    feasibility = "balanced",
    feasibility_tol = 1e-10,
    objective_recomputed = 0
  )
  expect_false(negative$converged)
  expect_identical(negative$status, "numerical_failure")

  infeasible_plan <- matrix(0.25, 2, 2)
  infeasible <- .attach_certified_example(
    list(plan = infeasible_plan, ot_dist = 0),
    infeasible_plan,
    p = c(0.8, 0.2),
    q = q,
    feasibility = "balanced",
    feasibility_tol = 1e-10,
    objective_recomputed = 0
  )
  expect_false(infeasible$converged)
  expect_identical(infeasible$status, "infeasible")
  expect_identical(
    rfugw:::.termination_reason_from_result(infeasible, 2L),
    "infeasible"
  )
})

test_that("nonfinite and inconsistent objectives cannot certify convergence", {
  plan <- diag(c(0.5, 0.5))
  nonfinite <- .attach_certified_example(
    list(plan = plan, ot_dist = Inf),
    plan,
    feasibility = "none",
    objective_recomputed = 0
  )
  mismatch <- .attach_certified_example(
    list(plan = plan, ot_dist = 1),
    plan,
    feasibility = "none",
    objective_recomputed = 2,
    objective_tolerance = 1e-12
  )
  component_mismatch <- .attach_certified_example(
    list(plan = plan, fugw_cost = 1, linear_cost = 0.4),
    plan,
    feasibility = "none",
    objective_recomputed = 1,
    objective_components = list(linear_cost = 0.5),
    objective_tolerance = 1e-12
  )

  expect_false(nonfinite$converged)
  expect_identical(nonfinite$status, "numerical_failure")
  expect_false(mismatch$converged)
  expect_identical(mismatch$status, "objective_mismatch")
  expect_false(component_mismatch$converged)
  expect_identical(component_mismatch$status, "objective_mismatch")
  expect_false(component_mismatch$linear_cost_consistent)
})

test_that("nonfinite inner certificates cannot imply unbalanced convergence", {
  plan <- diag(c(0.5, 0.5))
  out <- .attach_certified_example(
    list(plan = plan, ot_dist = 0),
    plan,
    inner_residual = Inf,
    max_inner_residual = Inf,
    inner_iterations = 1L,
    inner_converged = TRUE,
    feasibility = "unbalanced",
    feasibility_tol = 1e-8,
    objective_recomputed = 0
  )
  expect_false(out$converged)
  expect_false(out$feasible)
  expect_identical(out$status, "inner_failure")
})

test_that("partial convergence enforces inequalities and transported mass", {
  p <- c(0.6, 0.4)
  q <- c(0.5, 0.5)
  plan <- matrix(c(0.3, 0.2, 0.1, 0.1), 2, 2, byrow = TRUE)
  valid <- .attach_certified_example(
    list(plan = plan, partial_gw_dist = 0.2),
    plan,
    p = p,
    q = q,
    inner_residual = 0,
    max_inner_residual = 0,
    inner_converged = TRUE,
    feasibility = "partial",
    feasibility_tol = 1e-10,
    mass_target = 0.7,
    objective_recomputed = 0.2
  )
  wrong_mass <- .attach_certified_example(
    list(plan = plan, partial_gw_dist = 0.2),
    plan,
    p = p,
    q = q,
    inner_residual = 0,
    max_inner_residual = 0,
    inner_converged = TRUE,
    feasibility = "partial",
    feasibility_tol = 1e-10,
    mass_target = 0.8,
    objective_recomputed = 0.2
  )

  expect_true(valid$converged)
  expect_equal(valid$mass_residual, 0, tolerance = 1e-15)
  expect_false(wrong_mass$converged)
  expect_identical(wrong_mass$status, "infeasible")
})

test_that("public converged results imply feasibility and objective checks", {
  M <- matrix(c(0, 1, 1, 0), 2, 2)
  emd <- ot_emd(M)
  sink <- ot_sinkhorn(M, epsilon = 0.1, max_iter = 200L, tol = 1e-8)
  C <- M
  partial <- partial_gromov_wasserstein(
    C, C, m = 0.7, numItermax = 10L, log = TRUE
  )

  for (out in list(emd, sink, partial)) {
    expect_true(out$objective_consistent)
    expect_true(is.finite(out$objective_recomputed))
    if (isTRUE(out$converged)) {
      expect_true(out$feasible)
      if (!is.na(out$inner_converged)) {
        expect_true(out$inner_converged)
      }
    }
  }
})

.expect_supported_convergence_law <- function(out, label) {
  expect_true(
    all(c(
      "status", "converged", "termination_reason", "feasible",
      "objective_consistent", "objective_components_consistent"
    ) %in% names(out)),
    info = label
  )
  expect_identical(out$converged, identical(out$status, "converged"), info = label)
  if (isTRUE(out$converged)) {
    plan <- rfugw_plan(out)
    expect_true(all(is.finite(plan)), info = label)
    expect_true(all(plan >= 0), info = label)
    expect_true(out$feasible, info = label)
    expect_true(out$objective_consistent, info = label)
    expect_true(out$objective_components_consistent, info = label)
    if (!is.na(out$inner_converged)) {
      expect_true(out$inner_converged, info = label)
      expect_true(is.finite(out$inner_residual), info = label)
      expect_true(is.finite(out$max_inner_residual), info = label)
    }
    expect_null(out$warning_payload, info = label)
  } else {
    expect_true(is.list(out$warning_payload), info = label)
    expect_identical(out$warning_payload$code, out$status, info = label)
    expect_false(identical(out$termination_reason, "tolerance"), info = label)
  }
}

test_that("the shared convergence law covers supported formulation adapters", {
  M <- matrix(c(0, 1, 1, 0), 2, 2)
  C <- M
  p <- q <- c(0.5, 0.5)
  X <- matrix(c(0, 1, 1, 0, 0.5, 0.5), 3, 2, byrow = TRUE)
  Y <- matrix(c(0, 1, 0.8, 0.2, 0.4, 0.6), 3, 2, byrow = TRUE)
  results <- list(
    balanced_exact = ot_emd(M, p, q),
    balanced_entropic = ot_sinkhorn(
      M, p, q, epsilon = 0.1, max_iter = 200L, tol = 1e-8
    ),
    exact_fgw = fgw_exact_cg(
      M, C, C, p, q, max_iter = 8L, tol_rel = 1e-7, tol_abs = 1e-7
    ),
    entropic_fgw = fgw_entropic(
      M, C, C, p, q, epsilon = 0.1, max_iter = 40L,
      tol = 1e-7, sinkhorn_tol = 1e-8, precision = "strict_double"
    ),
    partial = partial_gromov_wasserstein(
      C, C, p, q, m = 0.7, numItermax = 8L, log = TRUE
    ),
    semirelaxed = semirelaxed_gromov_wasserstein(
      C, C, p, max_iter = 8L, tol_rel = 1e-7, tol_abs = 1e-7
    ),
    semirelaxed_entropic = entropic_semirelaxed_gromov_wasserstein(
      C, C, p, epsilon = 0.1, max_iter = 40L, tol = 1e-7,
      precision = "double"
    ),
    unbalanced_ot = ot_sinkhorn_unbalanced(
      M, p, q, epsilon = 0.1, max_iter = 200L, tol = 1e-7
    ),
    fugw = fugw_kl(
      C, C, p, q, M = M, epsilon = 0.05,
      max_iter = 20L, tol = 1e-6, max_iter_ot = 100L, tol_ot = 1e-6
    ),
    ucoot = unbalanced_co_optimal_transport(
      X, Y, epsilon = 0.05, max_iter = 20L, tol = 1e-6,
      max_iter_ot = 100L, tol_ot = 1e-6, log = TRUE
    )
  )
  for (label in names(results)) {
    .expect_supported_convergence_law(results[[label]], label)
  }
})

test_that("experimental adapters make no convergence claim", {
  set.seed(45)
  X1 <- matrix(rnorm(18), 6, 3)
  X2 <- matrix(rnorm(18), 6, 3)
  C1 <- as.matrix(dist(X1))
  C2 <- as.matrix(dist(X2))
  sampled <- sampled_gromov_wasserstein(
    C1, C2, nb_samples_grad = c(3L, 2L), epsilon = 0.1,
    max_iter = 4L, random_state = 45L, log = TRUE
  )
  approximate <- lowrank_gromov_wasserstein_samples(
    X1, X2, reg = 0.1, rank = 2L, numItermax = 10L, log = TRUE
  )
  for (out in list(sampled = sampled, approximate = approximate)) {
    expect_identical(out$status, "experimental")
    expect_false(out$converged)
    expect_identical(out$certification, "experimental_no_convergence_claim")
    expect_identical(
      out$termination_reason,
      "experimental_no_convergence_certificate"
    )
  }
})

test_that("every injected failure maps to a distinct status and warning payload", {
  plan <- diag(c(0.5, 0.5))
  base <- list(plan = plan, ot_dist = 0)
  max_iter <- rfugw:::.attach_solver_diagnostics(
    base, residual = 1, converged = FALSE, iterations = 2L, max_iter = 2L,
    feasibility = "none", objective_recomputed = 0
  )
  lp_failure <- rfugw:::.attach_solver_diagnostics(
    base, residual = 0, converged = TRUE, iterations = 1L, max_iter = 2L,
    feasibility = "none", objective_recomputed = 0, lp_ok = FALSE
  )
  expect_identical(max_iter$status, "max_iter")
  expect_identical(
    rfugw:::.termination_reason_from_result(max_iter, 2L), "max_iter"
  )
  expect_identical(lp_failure$status, "lp_failure")
  expect_identical(
    rfugw:::.termination_reason_from_result(lp_failure, 2L), "lp_failure"
  )
  expect_identical(max_iter$warning_payload$code, "max_iter")
  expect_identical(lp_failure$warning_payload$code, "lp_failure")

  stagnation <- max_iter
  stagnation$iterations <- 1L
  expect_identical(
    rfugw:::.termination_reason_from_result(stagnation, 2L), "stagnation"
  )
  expect_false(max_iter$converged)
  expect_false(lp_failure$converged)
})
