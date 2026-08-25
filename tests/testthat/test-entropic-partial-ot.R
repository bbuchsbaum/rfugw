independent_partial_entropy_dual <- function(M, p, q, mass, epsilon) {
  ns <- nrow(M)
  nt <- ncol(M)
  evaluate <- function(parameters) {
    source <- parameters[seq_len(ns)]
    target <- parameters[ns + seq_len(nt)]
    log_weight <- -(M + outer(source, target, "+")) / epsilon
    plan <- exp(
      log_weight + log(mass) - .log_sum_exp(as.numeric(log_weight))
    )
    list(
      value = epsilon * mass + sum(source * p) + sum(target * q) +
        epsilon * mass *
          (.log_sum_exp(as.numeric(log_weight)) - log(mass)),
      gradient = c(p - rowSums(plan), q - colSums(plan)),
      plan = plan
    )
  }
  fit <- stats::optim(
    rep(0, ns + nt),
    fn = function(x) evaluate(x)$value,
    gr = function(x) evaluate(x)$gradient,
    method = "L-BFGS-B",
    lower = 0,
    control = list(factr = 1, pgtol = 1e-14, maxit = 10000L)
  )
  stopifnot(fit$convergence == 0L)
  list(plan = evaluate(fit$par)$plan, fit = fit)
}

test_that("scaling and genuine log Dykstra agree in the safe regime", {
  M <- matrix(c(0, 0.4, 1.1, 0.2, 0.8, 0.3), 3L, 2L)
  p <- c(0.2, 0.3, 0.5)
  q <- c(0.65, 0.35)
  scaling <- ot_partial_sinkhorn(
    M, p, q, mass = 0.7, epsilon = 0.2,
    method = "scaling", max_iter = 10000L, tol = 1e-9,
    check_every = 1L
  )
  log_domain <- ot_partial_sinkhorn(
    M, p, q, mass = 0.7, epsilon = 0.2,
    method = "log", max_iter = 10000L, tol = 1e-9,
    check_every = 1L
  )

  expect_true(scaling$converged)
  expect_true(log_domain$converged)
  expect_equal(log_domain$plan, scaling$plan, tolerance = 2e-9)
  expect_equal(rfugw_value(log_domain), rfugw_value(scaling), tolerance = 2e-9)
  expect_identical(scaling$effective_sinkhorn_method, "scaling")
  expect_identical(log_domain$effective_sinkhorn_method, "log")
  expect_lte(abs(scaling$duality_gap), scaling$duality_gap_tolerance)
  expect_lte(abs(log_domain$duality_gap), log_domain$duality_gap_tolerance)
})

test_that("log solution matches an independent constrained dual optimizer", {
  M <- matrix(c(0.1, 0.9, 0.4, 0.2, 1.1, 0.3), 2L, 3L)
  p <- c(0.45, 0.55)
  q <- c(0.2, 0.3, 0.5)
  mass <- 0.72
  epsilon <- 0.17
  oracle <- independent_partial_entropy_dual(M, p, q, mass, epsilon)
  out <- ot_partial_sinkhorn(
    M, p, q, mass, epsilon, method = "log",
    max_iter = 20000L, tol = 1e-11, check_every = 1L
  )

  expect_true(out$converged)
  expect_equal(out$plan, oracle$plan, tolerance = 2e-9)
  expect_lte(max(rowSums(out$plan) - p), out$feasibility_tolerance)
  expect_lte(max(colSums(out$plan) - q), out$feasibility_tolerance)
  expect_equal(sum(out$plan), mass, tolerance = 1e-11)
})

test_that("auto uses log for adversarial range and scaling fails closed", {
  M <- matrix(c(0, 1200, 900, 0), 2L)
  out <- ot_partial_sinkhorn(
    M, mass = 0.7, epsilon = 0.1, method = "auto",
    max_iter = 10000L, tol = 1e-9, check_every = 1L
  )
  expect_true(out$converged)
  expect_identical(out$effective_sinkhorn_method, "log")
  expect_identical(out$sinkhorn_backend_transition, "auto_to_log")
  expect_true(all(is.finite(out$plan)))
  expect_equal(sum(out$plan), 0.7, tolerance = 1e-10)
  expect_error(
    ot_partial_sinkhorn(
      M, mass = 0.7, epsilon = 0.1, method = "scaling"
    ),
    "outside its certified regime"
  )
})

test_that("zero and full mass have the declared boundary behavior", {
  M <- matrix(c(0, 0.7, 1.1, 0.2), 2L)
  zero <- ot_partial_sinkhorn(M, mass = 0, epsilon = 0.2)
  full <- ot_partial_sinkhorn(
    M, mass = 1, epsilon = 0.2, method = "log",
    max_iter = 20000L, tol = 1e-10, check_every = 1L
  )
  balanced <- ot_sinkhorn(
    M, epsilon = 0.2, method = "log", max_iter = 20000L, tol = 1e-10
  )

  expect_true(zero$converged)
  expect_equal(zero$plan, matrix(0, 2L, 2L))
  expect_equal(rfugw_value(zero), 0)
  expect_true(full$converged)
  expect_equal(full$plan, balanced$plan, tolerance = 2e-9)
  expect_equal(rowSums(full$plan), c(0.5, 0.5), tolerance = 2e-9)
  expect_equal(colSums(full$plan), c(0.5, 0.5), tolerance = 2e-9)
})

test_that("zero weights, rectangles, duplicates, constants, and tiny entropy certify", {
  zero_weights <- ot_partial_sinkhorn(
    matrix(c(0, 2, 1, 3, 0.2, 1.5), 2L, 3L),
    p = c(1, 0), q = c(0.2, 0.3, 0.5),
    mass = 0.6, epsilon = 0.08, method = "log",
    max_iter = 20000L, tol = 1e-9, check_every = 1L
  )
  constant <- ot_partial_sinkhorn(
    matrix(1, 2L, 3L), p = c(0.2, 0.8), q = c(0.2, 0.3, 0.5),
    mass = 0.75, epsilon = 0.3, method = "scaling",
    max_iter = 20000L, tol = 1e-9, check_every = 1L
  )
  duplicate <- ot_partial_sinkhorn(
    matrix(c(0, 0, 1, 1, 0, 0), 2L, 3L),
    mass = 0.8, epsilon = 1e-3, method = "log",
    max_iter = 10000L, tol = 1e-8, check_every = 20L
  )

  for (out in list(zero_weights, constant, duplicate)) {
    expect_true(out$converged)
    expect_true(out$feasible)
    expect_true(out$objective_consistent)
    expect_true(out$objective_components_consistent)
    expect_lte(abs(out$duality_gap), out$duality_gap_tolerance)
  }
  expect_equal(zero_weights$plan[2, ], c(0, 0, 0))
  expect_equal(dim(constant$plan), c(2L, 3L))
  expect_true(all(is.finite(duplicate$plan)))
})

test_that("regularized objective and every component are recomputed", {
  M <- matrix(c(0.2, 0.7, 1.2, 0.1, 0.8, 0.4), 2L, 3L)
  epsilon <- 0.25
  out <- ot_partial_sinkhorn(
    M, p = c(0.4, 0.6), q = c(0.2, 0.3, 0.5),
    mass = 0.65, epsilon = epsilon, method = "log",
    max_iter = 20000L, tol = 1e-10, check_every = 1L
  )
  entropy <- sum(ifelse(out$plan > 0, out$plan * log(out$plan), 0))
  expected <- sum(M * out$plan) + epsilon * (entropy - sum(out$plan))

  expect_true(out$converged)
  expect_identical(out$entropy_convention, "counting_measure_entropy_minus_one")
  expect_equal(out$entropy, entropy, tolerance = 1e-12)
  expect_equal(out$entropy_minus_one, entropy - 0.65, tolerance = 1e-12)
  expect_equal(rfugw_value(out), expected, tolerance = 1e-12)
  expect_equal(out$objective_recomputed, expected, tolerance = 1e-12)
  expect_lte(out$stationarity_residual, out$kkt_tolerance)
  expect_lte(out$complementarity_residual, out$kkt_tolerance)
  expect_lte(out$dual_feasibility_residual, out$kkt_tolerance)
})

test_that("certified warm states preserve quality and reject mutations", {
  M <- matrix(c(0, 0.4, 1.1, 0.2, 0.8, 0.3), 3L, 2L)
  p <- c(0.2, 0.3, 0.5)
  q <- c(0.65, 0.35)
  cold <- ot_partial_sinkhorn(
    M, p, q, 0.7, 0.2, method = "log",
    max_iter = 20000L, tol = 1e-9, check_every = 1L
  )
  warm <- ot_partial_sinkhorn(
    M, p, q, 0.7, 0.2, method = "log",
    max_iter = 20000L, tol = 1e-9, check_every = 1L,
    init_state = cold
  )

  expect_true(cold$converged)
  expect_true(warm$converged)
  expect_true(warm$warm_started)
  expect_lte(warm$iterations, cold$iterations)
  expect_equal(warm$plan, cold$plan, tolerance = 2e-8)
  expect_equal(rfugw_value(warm), rfugw_value(cold), tolerance = 2e-8)

  wrong_cost <- cold$warm_state
  wrong_cost$cost[1, 1] <- wrong_cost$cost[1, 1] + 1
  expect_error(
    ot_partial_sinkhorn(M, p, q, 0.7, 0.2, method = "log", init_state = wrong_cost),
    "does not match"
  )
  broken_invariant <- cold$warm_state
  broken_invariant$corrections[[1]][1, 1] <-
    broken_invariant$corrections[[1]][1, 1] + 0.5
  expect_error(
    ot_partial_sinkhorn(
      M, p, q, 0.7, 0.2, method = "log", init_state = broken_invariant
    ),
    "invariant"
  )
  uncertified <- cold
  uncertified$converged <- FALSE
  expect_error(
    ot_partial_sinkhorn(M, p, q, 0.7, 0.2, method = "log", init_state = uncertified),
    "certified and converged"
  )
})

test_that("checked iterations record exact mass and final capacities", {
  out <- ot_partial_sinkhorn(
    matrix(c(0.1, 0.7, 1.3, 0.2, 0.8, 0.4), 2L, 3L),
    p = c(0.4, 0.6), q = c(0.2, 0.3, 0.5),
    mass = 0.68, epsilon = 0.15, method = "log",
    max_iter = 20000L, tol = 1e-10, check_every = 1L
  )
  expect_true(out$converged)
  expect_true(nrow(out$convergence_trace) > 1L)
  expect_true(all(out$convergence_trace$mass_residual <= 1e-12))
  expect_lte(tail(out$convergence_trace$row_violation, 1), out$feasibility_tolerance)
  expect_lte(tail(out$convergence_trace$col_violation, 1), out$feasibility_tolerance)
})

test_that("validation rejects infeasible requests before iteration", {
  M <- matrix(c(0, 1, 1, 0), 2L)
  expect_error(ot_partial_sinkhorn(M, mass = -0.1), "[[]0, 1[]]")
  expect_error(ot_partial_sinkhorn(M, mass = 1.1), "[[]0, 1[]]")
  expect_error(ot_partial_sinkhorn(M, mass = 0.5, epsilon = 0), "epsilon")
  expect_error(ot_partial_sinkhorn(M, mass = 0.5, method = "unknown"), "arg")
  expect_error(ot_partial_sinkhorn(M, mass = 0.5, check_every = 0), "check_every")
})

test_that("public primitive has parity with the former partial-GW projection", {
  set.seed(37)
  ns <- 3L
  nt <- 4L
  C1 <- matrix(runif(ns * ns, 0.05, 1.1), ns, ns); diag(C1) <- 0
  C2 <- matrix(runif(nt * nt, 0.05, 1.2), nt, nt); diag(C2) <- 0
  M <- matrix(runif(ns * nt), ns, nt)
  p <- c(0.2, 0.3, 0.5)
  q <- c(0.1, 0.2, 0.3, 0.4)
  mass <- 0.65
  alpha <- 0.55
  reg <- 0.2
  G0 <- (p %o% q) * mass
  gradient <- rfugw:::.gw_square_terms(C1, C2, G0, symmetric = FALSE)$grad
  linearized_cost <- (1 - alpha) * M + alpha * gradient
  former <- rfugw:::cpp_entropic_partial_wasserstein(
    p, q, linearized_cost, reg, mass, 5000L, 1e-12, FALSE, FALSE
  )
  public <- ot_partial_sinkhorn(
    linearized_cost, p, q, mass, reg, method = "scaling",
    max_iter = 5000L, tol = 1e-12, check_every = 1L
  )
  outer <- entropic_partial_fused_gromov_wasserstein(
    M, C1, C2, p, q, reg = reg, m = mass, alpha = alpha,
    G0 = G0, numItermax = 1L, tol = 1,
    inner_max_iter = 5000L, inner_tol = 1e-12,
    method = "scaling", log = TRUE, symmetric = FALSE
  )

  expect_true(public$converged)
  expect_equal(public$plan, former, tolerance = 5e-10)
  expect_equal(outer$plan, public$plan, tolerance = 5e-10)
  expect_true(outer$inner_converged)
  expect_identical(outer$inner_status, "converged")
  expect_identical(outer$requested_sinkhorn_method, "scaling")
})

test_that("partial Sinkhorn benchmark records dispatch, work, and allocation evidence", {
  baseline <- utils::read.csv(
    bench_test_resource("partial-sinkhorn-baseline.csv"),
    stringsAsFactors = FALSE
  )
  required <- c(
    "case", "requested_method", "effective_method", "warm_started",
    "iterations", "residual", "feasibility_residual",
    "objective_residual", "duality_gap", "dynamic_range",
    "scaling_threshold", "elapsed_ms", "allocation_count",
    "allocation_bytes", "status", "certified", "seed", "r_version",
    "sysname", "machine", "commit"
  )
  expect_true(all(required %in% names(baseline)))
  expect_true(all(baseline$certified))
  expect_true(all(baseline$status == "converged"))
  expect_true(all(is.finite(baseline$iterations)))
  expect_true(all(is.finite(baseline$elapsed_ms)))
  expect_true(all(is.finite(baseline$allocation_bytes)))
  expect_identical(
    baseline$effective_method[baseline$case == "adversarial_auto"],
    "log"
  )
  expect_lte(
    baseline$iterations[baseline$case == "moderate_log_warm"],
    baseline$iterations[baseline$case == "moderate_log"]
  )
})
