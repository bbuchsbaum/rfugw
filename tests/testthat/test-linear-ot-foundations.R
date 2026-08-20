test_that("exact partial OT satisfies analytic mass boundaries", {
  M <- matrix(c(0, 4, 3, 1), 2, 2)
  zero <- ot_partial_emd(M, mass = 0)
  half <- ot_partial_emd(M, mass = 0.5)
  full <- ot_partial_emd(M, mass = 1)
  balanced <- ot_emd(M)

  for (out in list(zero, half, full)) {
    expect_identical(out$status, "converged")
    expect_true(out$feasible)
    expect_true(out$objective_consistent)
    expect_true(out$objective_components_consistent)
    expect_true(out$mass_certified)
    expect_equal(out$duality_gap, 0, tolerance = 1e-12)
  }
  expect_equal(sum(zero$plan), 0, tolerance = 1e-12)
  expect_equal(sum(half$plan), 0.5, tolerance = 1e-12)
  expect_equal(full$plan, balanced$plan, tolerance = 1e-12)
  expect_equal(full$ot_dist, balanced$ot_dist, tolerance = 1e-12)
})

test_that("exact partial OT matches an independent LP oracle", {
  skip_if_not_installed("lpSolve")
  set.seed(20260820L)
  ns <- 3L
  nt <- 4L
  M <- matrix(runif(ns * nt), ns, nt)
  p <- c(0.2, 0.3, 0.5)
  q <- c(0.1, 0.2, 0.3, 0.4)
  mass <- 0.63

  row_constraints <- matrix(0, ns, ns * nt)
  for (i in seq_len(ns)) row_constraints[i, i + ns * (seq_len(nt) - 1L)] <- 1
  col_constraints <- matrix(0, nt, ns * nt)
  for (j in seq_len(nt)) col_constraints[j, seq_len(ns) + ns * (j - 1L)] <- 1
  oracle <- lpSolve::lp(
    direction = "min",
    objective.in = as.vector(M),
    const.mat = rbind(row_constraints, col_constraints, rep(1, ns * nt)),
    const.dir = c(rep("<=", ns + nt), "="),
    const.rhs = c(p, q, mass),
    all.int = FALSE
  )
  expect_equal(oracle$status, 0L)

  out <- ot_partial_emd(M, p, q, mass = mass)
  expect_identical(out$status, "converged")
  expect_equal(out$ot_dist, oracle$objval, tolerance = 1e-10)
  expect_no_error(ot_validate_plan(out, p, q, mass, marginals = "partial", tol = 1e-10))
})

test_that("exact partial OT is permutation equivariant and cost homogeneous", {
  M <- matrix(c(0.1, 1.7, 2.4, 0.4, 1.2, 3.1), 3, 2)
  p <- c(0.2, 0.3, 0.5)
  q <- c(0.6, 0.4)
  mass <- 0.7
  base <- ot_partial_emd(M, p, q, mass)
  ip <- c(3L, 1L, 2L)
  jp <- c(2L, 1L)
  perm <- ot_partial_emd(M[ip, jp], p[ip], q[jp], mass)
  restored <- matrix(0, 3, 2)
  restored[ip, jp] <- perm$plan

  expect_equal(restored, base$plan, tolerance = 1e-12)
  expect_equal(perm$ot_dist, base$ot_dist, tolerance = 1e-12)
  scaled <- ot_partial_emd(7 * M, p, q, mass)
  expect_equal(scaled$plan, base$plan, tolerance = 1e-12)
  expect_equal(scaled$ot_dist, 7 * base$ot_dist, tolerance = 1e-11)
})

test_that("exact partial OT rejects unsupported or infeasible inputs", {
  M <- matrix(c(0, 1, 1, 0), 2)
  expect_error(ot_partial_emd(M, mass = -0.1), "[[]0, 1[]]")
  expect_error(ot_partial_emd(M, mass = 1.1), "[[]0, 1[]]")
  expect_error(ot_partial_emd(M, mass = c(0.2, 0.3)), "one finite")
  expect_error(ot_partial_emd(M - 2, mass = 0.5), "nonnegative")
  expect_error(ot_partial_emd(M, p = c(1, 0, 0)), "length")
})

test_that("genuine log unbalanced Sinkhorn has oracle and permutation laws", {
  M <- matrix(c(0, 0.4, 1.1, 0.2, 0.8, 0.3), 3, 2)
  p <- c(0.2, 0.3, 0.5)
  q <- c(0.65, 0.35)
  scaling <- ot_sinkhorn_unbalanced(
    M, p, q, epsilon = 0.2, rho = c(2, 3), method = "scaling",
    max_iter = 2000L, tol = 1e-7
  )
  log_domain <- ot_sinkhorn_unbalanced(
    M, p, q, epsilon = 0.2, rho = c(2, 3), method = "log",
    max_iter = 2000L, tol = 1e-7
  )
  expect_equal(log_domain$plan, scaling$plan, tolerance = 2e-6)
  expect_equal(log_domain$ot_dist, ot_linear_cost(M, log_domain), tolerance = 1e-12)

  ip <- c(3L, 1L, 2L)
  jp <- c(2L, 1L)
  perm <- ot_sinkhorn_unbalanced(
    M[ip, jp], p[ip], q[jp], epsilon = 0.2, rho = c(2, 3),
    method = "log", max_iter = 2000L, tol = 1e-7
  )
  restored <- matrix(0, 3, 2)
  restored[ip, jp] <- perm$plan
  expect_equal(restored, log_domain$plan, tolerance = 2e-6)
})

test_that("adversarial unbalanced dynamic range routes to genuine log", {
  M <- matrix(c(0, 1200, 900, 0), 2)
  out <- ot_sinkhorn_unbalanced(
    M, epsilon = 0.1, rho = 2, method = "auto",
    max_iter = 3000L, tol = 1e-7
  )
  expect_identical(out$effective_sinkhorn_method, "log")
  expect_true(all(is.finite(out$plan)))
  expect_true(out$objective_consistent)
  expect_error(
    ot_sinkhorn_unbalanced(M, epsilon = 0.1, rho = 2, method = "scaling"),
    "outside its certified regime"
  )
})

test_that("deferred linear foundations are not silently exported", {
  expect_false(exists("ot_partial_sinkhorn", asNamespace("rfugw"), inherits = FALSE))
  expect_false(exists("ot_sinkhorn_divergence", asNamespace("rfugw"), inherits = FALSE))
  expect_false(exists("ot_fixed_support_barycenter", asNamespace("rfugw"), inherits = FALSE))
  decision <- readLines(
    trust_test_resource("linear-ot-foundations.md"),
    warn = FALSE
  )
  expect_true(any(grepl("Deferred", decision, fixed = TRUE)))
  expect_true(any(grepl("log-domain Dykstra", decision, fixed = TRUE)))
  expect_true(any(grepl("regularized primal", decision, fixed = TRUE)))
})

test_that("accepted foundations retain a certificate-first representative baseline", {
  baseline <- utils::read.csv(
    bench_test_resource("linear-foundations-baseline.csv"),
    stringsAsFactors = FALSE
  )
  required <- c(
    "method", "evidence_channel", "certified", "status", "objective",
    "residual", "feasibility_residual", "objective_residual", "mass_residual",
    "setup_ms", "solve_ms", "e2e_ms", "mem_solve_bytes", "mem_e2e_bytes",
    "requested_precision", "effective_precision", "requested_threads",
    "used_threads", "commit", "r_version", "sysname", "machine", "seed"
  )
  expect_true(all(required %in% names(baseline)))
  expect_true(all(baseline$certified))
  expect_true(all(baseline$status == "converged"))
  expect_setequal(
    unique(baseline$method),
    c("ot_partial_emd", "ot_sinkhorn_unbalanced_log")
  )
  expect_true(all(is.finite(baseline$solve_ms)))
  expect_true(all(baseline$evidence_channel == "local"))
})
