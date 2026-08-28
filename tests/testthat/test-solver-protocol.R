protocol_fixture <- function(seed = 104L, ns = 13L, nt = 9L) {
  set.seed(seed)
  source_points <- matrix(rnorm(ns * 3L), ns, 3L)
  target_points <- matrix(rnorm(nt * 3L), nt, 3L)
  cost <- outer(rowSums(source_points^2), rowSums(target_points^2), "+") -
    2 * tcrossprod(source_points, target_points)
  cost[cost < 0 & cost > -1e-12] <- 0
  list(
    cost = cost,
    source = seq_len(ns) / sum(seq_len(ns)),
    target = rev(seq_len(nt)) / sum(seq_len(nt))
  )
}

test_that("protocol problems make estimand, orientation, mass, and units explicit", {
  d <- protocol_fixture()
  expect_error(
    transport_problem_sinkhorn(d$cost, d$source, d$target),
    "mass_policy.*explicit"
  )
  expect_error(
    transport_problem_sinkhorn(
      d$cost, 2 * d$source, d$target, mass_policy = "probability"
    ),
    "requires each supplied measure to sum to one"
  )

  problem <- transport_problem_sinkhorn(
    d$cost, 2 * d$source, 3 * d$target,
    epsilon = 0.6, method = "auto", mass_policy = "normalize"
  )
  expect_s3_class(problem, "rfugw_transport_problem")
  expect_identical(problem$protocol_version, "1.0")
  expect_identical(problem$estimand, "balanced_entropic_linear_ot")
  expect_identical(problem$orientation, "rows_source_columns_target")
  expect_identical(problem$objective$units, "cost_times_probability_mass")
  expect_equal(problem$original_mass, c(2, 3))
  expect_equal(problem$effective_mass, c(1, 1))
  expect_equal(problem$source_measure, d$source)
  expect_equal(problem$target_measure, d$target)

  result <- transport_solve(problem)
  expect_true(result$converged)
  expect_identical(rfugw_problem(result), problem)
  expect_identical(rfugw_status(result), "converged")
  expect_false(result$warm_started)
  expect_false(result$warm_start_accepted)
  expect_identical(result$warm_start_rejection$code, "not_requested")
  expect_equal(rfugw_masses(result)$transported_mass, 1, tolerance = 1e-9)
  expect_identical(
    rfugw_provenance(result)$plan_representation$representation,
    "dense_materialized"
  )
})

test_that("balanced cold and warm protocol solves earn the same certificate", {
  d <- protocol_fixture()
  problem <- transport_problem_sinkhorn(
    d$cost, d$source, d$target, epsilon = 0.7, method = "log",
    max_iter = 5000L, tol = 1e-10, mass_policy = "probability"
  )
  cold <- transport_solve(problem)
  state <- rfugw_state(cold)
  warm <- transport_solve(problem, state)

  expect_true(cold$converged)
  expect_true(warm$converged)
  expect_s3_class(state, "rfugw_solver_state")
  expect_true(warm$warm_started)
  expect_true(warm$warm_start_accepted)
  expect_null(warm$warm_start_rejection)
  expect_identical(warm$status, cold$status)
  expect_equal(rfugw_plan(warm), rfugw_plan(cold), tolerance = 2e-10)
  expect_equal(rfugw_value(warm), rfugw_value(cold), tolerance = 2e-10)
  expect_equal(warm$regularized_objective, cold$regularized_objective,
               tolerance = 2e-10)
  expect_lte(warm$residual, problem$controls$tol)
  expect_lte(warm$iterations, cold$iterations)
})

test_that("tighter tolerance and epsilon continuation reuse portable dual state", {
  d <- protocol_fixture(seed = 9L, ns = 20L, nt = 15L)
  loose_problem <- transport_problem_sinkhorn(
    d$cost, d$source, d$target, epsilon = 0.9, method = "scaling",
    max_iter = 5000L, tol = 1e-6, mass_policy = "probability"
  )
  loose <- transport_solve(loose_problem)
  tight_problem <- transport_problem_sinkhorn(
    d$cost, d$source, d$target, epsilon = 0.55, method = "log",
    max_iter = 5000L, tol = 1e-10, mass_policy = "probability"
  )
  warm <- transport_solve(tight_problem, rfugw_state(loose))
  cold <- transport_solve(tight_problem)

  expect_true(loose$converged)
  expect_true(warm$converged)
  expect_true(cold$converged)
  expect_true(warm$warm_start_accepted)
  expect_identical(warm$effective_sinkhorn_method, "log")
  expect_equal(warm$plan, cold$plan, tolerance = 2e-9)
  expect_equal(warm$regularized_objective, cold$regularized_objective,
               tolerance = 2e-9)
  expect_lte(warm$iterations, cold$iterations)
})

test_that("state validation rejects malformed or incompatible inputs before solve", {
  d <- protocol_fixture()
  problem <- transport_problem_sinkhorn(
    d$cost, d$source, d$target, epsilon = 0.6, method = "log",
    mass_policy = "probability"
  )
  state <- rfugw_state(transport_solve(problem))

  bad_dim <- state
  bad_dim$dimensions[[1L]] <- bad_dim$dimensions[[1L]] + 1L
  expect_error(
    transport_solve(problem, bad_dim),
    "dimension_mismatch"
  )
  bad_finite <- state
  bad_finite$payload$source[[1L]] <- NA_real_
  expect_error(transport_solve(problem, bad_finite), "nonfinite_payload")
  bad_support <- state
  bad_support$source_support[[1L]] <- FALSE
  expect_error(transport_solve(problem, bad_support), "support_mismatch")
  bad_version <- state
  bad_version$protocol_version <- "999"
  expect_error(
    transport_solve(problem, bad_version),
    "unsupported_state_version"
  )

  exact <- transport_problem_emd(
    d$cost, d$source, d$target, mass_policy = "probability"
  )
  expect_error(transport_solve(exact, state), "state_unsupported")

  overflow_problem <- transport_problem_sinkhorn(
    d$cost, d$source, d$target, epsilon = 0.6, method = "scaling",
    mass_policy = "probability"
  )
  overflow <- state
  overflow$payload$source[[1L]] <- 1000
  expect_error(transport_solve(overflow_problem, overflow), "backend_overflow")

  fallback <- transport_solve(problem, bad_dim, state_policy = "cold")
  expect_true(fallback$converged)
  expect_false(fallback$warm_started)
  expect_false(fallback$warm_start_accepted)
  expect_identical(fallback$warm_start_rejection$code, "dimension_mismatch")

  mutated <- problem
  mutated$cost[[1L]] <- NA_real_
  expect_error(transport_solve(mutated), "finite")
  mutated <- problem
  mutated$controls$method <- "not-a-method"
  expect_error(transport_solve(mutated), "should be one of")
})

test_that("protocol state survives serialization without field inspection", {
  d <- protocol_fixture()
  problem <- transport_problem_sinkhorn(
    d$cost, d$source, d$target, epsilon = 0.45, method = "auto",
    max_iter = 5000L, tol = 1e-10, mass_policy = "probability"
  )
  first <- transport_solve(problem)
  path <- tempfile(fileext = ".rds")
  saveRDS(list(problem = problem, state = rfugw_state(first)), path)
  restored <- readRDS(path)
  second <- transport_solve(restored$problem, restored$state)
  expect_true(second$converged)
  expect_true(second$warm_start_accepted)
  expect_equal(second$plan, first$plan, tolerance = 2e-10)
})

test_that("exact and partial estimands use the same boundary without conflation", {
  d <- protocol_fixture(ns = 5L, nt = 4L)
  exact_problem <- transport_problem_emd(
    d$cost, d$source, d$target, mass_policy = "probability"
  )
  exact <- transport_solve(exact_problem)
  direct_exact <- ot_emd(d$cost, d$source, d$target)
  expect_true(exact$converged)
  expect_equal(exact$plan, direct_exact$plan, tolerance = 1e-12)
  expect_equal(rfugw_value(exact), ot_linear_cost(d$cost, exact$plan),
               tolerance = 1e-12)
  expect_error(rfugw_state(exact), "does not expose reusable")

  partial_problem <- transport_problem_partial_emd(
    d$cost, d$source, d$target, mass = 0.65,
    mass_policy = "probability"
  )
  partial <- transport_solve(partial_problem)
  expect_true(partial$converged)
  expect_equal(sum(partial$plan), 0.65, tolerance = 1e-10)
  expect_identical(
    rfugw_provenance(partial)$estimand,
    "exact_fixed_mass_partial_linear_ot"
  )
  expect_false(identical(rfugw_value(partial), rfugw_value(exact)))
})

test_that("partial Dykstra state stays bound to problem and effective backend", {
  d <- protocol_fixture(ns = 4L, nt = 5L)
  scaling_problem <- transport_problem_partial_sinkhorn(
    d$cost, d$source, d$target, mass = 0.7, epsilon = 0.8,
    method = "scaling", max_iter = 10000L, tol = 1e-8,
    check_every = 5L, mass_policy = "probability"
  )
  first <- transport_solve(scaling_problem)
  state <- rfugw_state(first)
  second <- transport_solve(scaling_problem, state)
  expect_true(first$converged)
  expect_true(second$converged)
  expect_true(second$warm_start_accepted)
  expect_equal(second$plan, first$plan, tolerance = 2e-8)

  log_problem <- transport_problem_partial_sinkhorn(
    d$cost, d$source, d$target, mass = 0.7, epsilon = 0.8,
    method = "log", max_iter = 10000L, tol = 1e-8,
    check_every = 5L, mass_policy = "probability"
  )
  expect_error(transport_solve(log_problem, state), "backend_mismatch")
})

test_that("KL-unbalanced protocol preserves explicit finite-measure policy", {
  d <- protocol_fixture(ns = 4L, nt = 3L)
  source <- 2.5 * d$source
  target <- 1.7 * d$target
  finite_problem <- transport_problem_unbalanced(
    d$cost, source, target, epsilon = 0.7, rho = c(1.2, 1.8),
    method = "log", max_iter = 5000L, tol = 1e-8,
    mass_policy = "finite_measure"
  )
  protocol_fit <- transport_solve(finite_problem)
  direct_fit <- ot_sinkhorn_unbalanced(
    d$cost, source, target, epsilon = 0.7, rho = c(1.2, 1.8),
    method = "log", max_iter = 5000L, tol = 1e-8,
    normalization = "none"
  )
  expect_true(protocol_fit$converged)
  expect_equal(protocol_fit$plan, direct_fit$plan, tolerance = 2e-10)
  expect_equal(rfugw_value(protocol_fit), rfugw_value(direct_fit),
               tolerance = 2e-10)
  expect_equal(rfugw_masses(protocol_fit)$source$input_mass, 2.5)
  expect_equal(rfugw_masses(protocol_fit)$target$input_mass, 1.7)

  normalized_a <- transport_solve(transport_problem_unbalanced(
    d$cost, source, target, epsilon = 0.7, rho = c(1.2, 1.8),
    method = "log", max_iter = 5000L, tol = 1e-8,
    mass_policy = "separate_probability"
  ))
  normalized_b <- transport_solve(transport_problem_unbalanced(
    d$cost, 11 * source, 0.2 * target, epsilon = 0.7,
    rho = c(1.2, 1.8), method = "log", max_iter = 5000L, tol = 1e-8,
    mass_policy = "separate_probability"
  ))
  expect_equal(normalized_a$plan, normalized_b$plan, tolerance = 2e-10)
  expect_equal(rfugw_value(normalized_a), rfugw_value(normalized_b),
               tolerance = 2e-10)
  expect_false(isTRUE(all.equal(protocol_fit$plan, normalized_a$plan,
                                tolerance = 1e-7)))
})

test_that("installed generic client uses only the public protocol", {
  d <- protocol_fixture()
  client <- new.env(parent = baseenv())
  sys.source(
    trust_test_resource(
      "extdata", "clients", "generic_solver_protocol_fixture.R"
    ),
    envir = client
  )
  first <- client$generic_transport_client(
    d$cost, d$source, d$target, epsilon = 0.7, tol = 1e-8
  )
  second <- client$generic_transport_client(
    d$cost, d$source, d$target, epsilon = 0.55, tol = 1e-9,
    state = first$state
  )
  exact <- client$generic_exact_transport_client(
    d$cost, d$source, d$target
  )
  expect_identical(first$status, "converged")
  expect_identical(second$status, "converged")
  expect_true(second$result$warm_start_accepted)
  expect_identical(exact$provenance$estimand, "exact_balanced_linear_ot")
  source_text <- readLines(
    trust_test_resource(
      "extdata", "clients", "generic_solver_protocol_fixture.R"
    ),
    warn = FALSE
  )
  expect_false(any(grepl("rfugw:::", source_text, fixed = TRUE)))
})

test_that("DKGE certified log-Sinkhorn is a downstream differential client", {
  skip_if_not_installed("dkge")
  d <- protocol_fixture(seed = 88L, ns = 8L, nt = 6L)
  client <- new.env(parent = baseenv())
  sys.source(
    trust_test_resource("extdata", "clients", "dkge_sinkhorn_fixture.R"),
    envir = client
  )
  dkge::dkge_clear_sinkhorn_cache()
  dkge_solver <- getFromNamespace(".dkge_sinkhorn_plan", "dkge")
  baseline <- dkge_solver(
    d$cost, d$source, d$target, epsilon = 0.8,
    max_iter = 10000L, tol = 1e-9, warm_start = FALSE,
    return_diagnostics = TRUE
  )
  fit <- client$dkge_rfugw_sinkhorn(
    d$cost, d$source, d$target, epsilon = 0.8, tol = 1e-9
  )
  expect_true(baseline$diagnostics$converged)
  expect_true(fit$diagnostics$converged)
  expect_equal(fit$plan, baseline$plan, tolerance = 3e-8)
  expect_equal(
    ot_linear_cost(d$cost, fit$plan),
    ot_linear_cost(d$cost, baseline$plan),
    tolerance = 2e-9
  )
  expect_lte(fit$diagnostics$marginal_error, 1e-9)

  dkge::dkge_clear_sinkhorn_cache()
  cached_first <- dkge_solver(
    d$cost, d$source, d$target, epsilon = 0.8,
    max_iter = 10000L, tol = 1e-9, warm_start = TRUE,
    return_diagnostics = TRUE
  )
  cached_second <- dkge_solver(
    d$cost, d$source, d$target, epsilon = 0.8,
    max_iter = 10000L, tol = 1e-9, warm_start = TRUE,
    return_diagnostics = TRUE
  )
  protocol_reuse <- client$dkge_rfugw_sinkhorn(
    d$cost, d$source, d$target, epsilon = 0.8, tol = 1e-9,
    state = fit$state
  )
  expect_true(cached_second$diagnostics$cache_hit)
  expect_identical(cached_second$diagnostics$iterations, 0L)
  expect_identical(protocol_reuse$diagnostics$iterations, 0L)
  expect_equal(cached_first$plan, cached_second$plan, tolerance = 0)
  expect_equal(protocol_reuse$plan, fit$plan, tolerance = 2e-10)

  expect_error(
    dkge_solver(
      d$cost, d$source, 0.9 * d$target, epsilon = 0.8,
      max_iter = 100L, tol = 1e-8
    ),
    "sum to the same total mass"
  )
  expect_error(
    client$dkge_rfugw_sinkhorn(
      d$cost, d$source, 0.9 * d$target, epsilon = 0.8, tol = 1e-8
    ),
    "requires each supplied measure to sum to one"
  )

  continued <- client$dkge_rfugw_sinkhorn(
    d$cost, d$source, d$target, epsilon = 0.8, tol = 1e-11,
    state = fit$state
  )
  cold <- client$dkge_rfugw_sinkhorn(
    d$cost, d$source, d$target, epsilon = 0.8, tol = 1e-11
  )
  expect_true(continued$diagnostics$converged)
  expect_true(continued$diagnostics$warm_start_accepted)
  expect_lte(continued$diagnostics$iterations, cold$diagnostics$iterations)
  expect_equal(continued$plan, cold$plan, tolerance = 2e-10)
})

test_that("capability matrix keeps scientific formulations distinct", {
  capabilities <- transport_capabilities()
  expect_true(all(c(
    "estimand", "public_solver", "protocol_constructor", "mass_policy",
    "warm_state", "plan_representation", "exact_certificate",
    "relational_costs", "maturity", "coverage_family"
  ) %in% names(capabilities)))
  expect_identical(anyDuplicated(capabilities$estimand), 0L)
  expect_true(all(c(
    "balanced_entropic_linear_ot", "exact_balanced_linear_ot",
    "exact_fixed_mass_partial_linear_ot",
    "entropic_fixed_mass_partial_linear_ot",
    "kl_unbalanced_entropic_linear_ot", "sinkhorn_divergence",
    "fixed_support_wasserstein_barycenter_weights",
    "gromov_wasserstein", "fused_gromov_wasserstein",
    "fused_unbalanced_gromov_wasserstein",
    "matrix_free_fused_unbalanced_gromov_wasserstein",
    "sampled_gromov_wasserstein"
  ) %in% capabilities$estimand))
  expect_true(capabilities$warm_state[
    capabilities$estimand == "balanced_entropic_linear_ot"
  ])
  expect_false(capabilities$protocol_constructor[
    capabilities$estimand == "gromov_wasserstein"
  ])
  expect_identical(
    capabilities$maturity[
      capabilities$estimand == "balanced_entropic_linear_ot"
    ],
    "flagship"
  )
  expect_identical(
    capabilities$maturity[
      capabilities$estimand == "sampled_gromov_wasserstein"
    ],
    "experimental"
  )
  expect_identical(
    capabilities$maturity[
      capabilities$estimand ==
        "matrix_free_fused_unbalanced_gromov_wasserstein"
    ],
    "experimental"
  )
})

test_that("capability evidence fails closed when maturity or ownership drifts", {
  paths <- utils::read.csv(
    trust_test_resource("numerical-path-matrix.csv"),
    stringsAsFactors = FALSE
  )
  capabilities <- transport_capabilities()
  expect_true(rfugw:::.validate_capability_path_matrix(capabilities, paths))

  missing <- paths[paths$family != "sinkhorn_divergence", , drop = FALSE]
  expect_error(
    rfugw:::.validate_capability_path_matrix(capabilities, missing),
    "lack numerical-path evidence.*sinkhorn_divergence"
  )

  drifted <- paths
  drifted$maturity[drifted$family == "balanced_linear_ot"] <- "supported"
  expect_error(
    rfugw:::.validate_capability_path_matrix(capabilities, drifted),
    "Maturity mismatch.*balanced_linear_ot"
  )

  unowned <- rbind(
    paths,
    data.frame(
      family = "invented", dimension = "backend", path = "a",
      comparison = "b", scope = "pr", maturity = "supported"
    )
  )
  expect_error(
    rfugw:::.validate_capability_path_matrix(capabilities, unowned),
    "lack capability ownership.*invented"
  )
})

test_that("release support prose is derived from accepted capability maturity", {
  paths <- utils::read.csv(
    trust_test_resource("numerical-path-matrix.csv"),
    stringsAsFactors = FALSE
  )
  evidence <- rfugw:::.release_support_evidence(
    transport_capabilities(), paths
  )
  verified <- paste(evidence$verified_support, collapse = "\n")
  boundaries <- paste(evidence$experimental_boundaries, collapse = "\n")

  expect_match(verified, "sinkhorn divergence")
  expect_match(verified, "entropic fixed mass partial linear ot")
  expect_match(verified, "fixed support wasserstein barycenter weights")
  expect_false(grepl("sinkhorn divergence", boundaries, fixed = TRUE))
  expect_false(grepl(
    "fixed support wasserstein barycenter weights", boundaries, fixed = TRUE
  ))
  expect_match(boundaries, "sampled gromov wasserstein")
  expect_match(boundaries, "matrix free fused unbalanced gromov wasserstein")
  expect_match(boundaries, "Directed-KL structural GW/FGW loss is deferred")
  expect_match(boundaries, "No end-to-end scalable relational-OT path")
})
