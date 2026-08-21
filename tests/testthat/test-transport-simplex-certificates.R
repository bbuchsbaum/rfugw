.expect_transport_certificate <- function(out, p, q) {
  certificate_fields <- c(
    "termination_reason", "source_potential", "target_potential",
    "primal_objective", "dual_objective", "duality_gap",
    "row_residual", "col_residual", "min_reduced_cost",
    "feasibility_tolerance", "reduced_cost_tolerance",
    "duality_gap_tolerance", "primal_feasible", "dual_feasible",
    "gap_certified"
  )
  expect_true(all(certificate_fields %in% names(out)))
  expect_identical(out$termination_reason, "optimal")
  expect_true(out$lp_ok)
  expect_true(out$primal_feasible)
  expect_true(out$dual_feasible)
  expect_true(out$gap_certified)
  expect_lte(out$row_residual, out$feasibility_tolerance)
  expect_lte(out$col_residual, out$feasibility_tolerance)
  expect_gte(out$min_reduced_cost, -out$reduced_cost_tolerance)
  expect_lte(abs(out$duality_gap), out$duality_gap_tolerance)
  expect_equal(rowSums(out$plan), p, tolerance = out$feasibility_tolerance)
  expect_equal(colSums(out$plan), q, tolerance = out$feasibility_tolerance)
  expect_equal(
    out$dual_objective,
    sum(p * out$source_potential) + sum(q * out$target_potential),
    tolerance = out$duality_gap_tolerance
  )
  expect_equal(out$primal_objective, out$ot_dist, tolerance = 0)
}

test_that("exact transport exposes a checkable optimality certificate", {
  M <- matrix(c(
    3, 0,
    0, 3
  ), 2, 2, byrow = TRUE)
  p <- c(0.4, 0.6)
  q <- c(0.5, 0.5)
  out <- ot_emd(M, p, q)

  .expect_transport_certificate(out, p, q)
  reduced_cost <- M - outer(
    as.numeric(out$source_potential),
    as.numeric(out$target_potential),
    "+"
  )
  expect_gte(min(reduced_cost), -out$reduced_cost_tolerance)
  expect_identical(out$status, "converged")
  expect_true(out$converged)
})

test_that("the uniform assignment fast path satisfies the same certificate", {
  M <- matrix(c(
    0, 4, 2,
    4, 0, 3,
    2, 3, 0
  ), 3, 3, byrow = TRUE)
  p <- q <- rep(1 / 3, 3)
  out <- ot_emd(M, p, q)

  .expect_transport_certificate(out, p, q)
  expect_equal(out$plan, diag(p), tolerance = 1e-14)
})

test_that("every transport-simplex termination branch is explicit", {
  M <- matrix(c(
    3, 0,
    0, 3
  ), 2, 2, byrow = TRUE)
  p <- c(0.4, 0.6)
  q <- c(0.5, 0.5)

  optimal <- cpp_transport_simplex_test(M, p, q, 100L, 1e-12, "none")
  expect_identical(optimal$termination_reason, "optimal")
  expect_true(optimal$lp_ok)

  limited <- cpp_transport_simplex_test(M, p, q, 0L, 1e-12, "none")
  expect_identical(limited$termination_reason, "max_iter")
  expect_false(limited$lp_ok)

  injected <- c(
    "disconnected_basis", "invalid_cycle", "invalid_step",
    "no_leaving_variable", "numerical_failure"
  )
  for (reason in injected) {
    out <- cpp_transport_simplex_test(M, p, q, 100L, 1e-12, reason)
    expect_identical(out$termination_reason, reason, info = reason)
    expect_false(out$lp_ok, info = reason)
  }
  expect_error(
    cpp_transport_simplex_test(M, p, q, 100L, 1e-12, "not_a_reason"),
    "unknown transport-simplex test fault"
  )
})

test_that("rectangular, tied, zero-weight, and degenerate cases agree with an LP oracle", {
  skip_if_not_installed("lpSolve")
  cases <- list(
    rectangular = list(
      M = matrix(c(0, 2, 1, 3, 0, 4), 2, 3, byrow = TRUE),
      p = c(0.45, 0.55),
      q = c(0.2, 0.3, 0.5)
    ),
    ties = list(
      M = matrix(c(0, 0, 2, 0, 0, 2, 2, 2, 0), 3, 3, byrow = TRUE),
      p = c(0.2, 0.3, 0.5),
      q = c(0.4, 0.1, 0.5)
    ),
    zero_weights = list(
      M = matrix(c(2, 0, 1, 1, 3, 0, 0, 2, 2), 3, 3, byrow = TRUE),
      p = c(0, 0.4, 0.6),
      q = c(0.5, 0, 0.5)
    ),
    degenerate = list(
      M = matrix(1, 2, 3),
      p = c(0.5, 0.5),
      q = c(0.2, 0.3, 0.5)
    )
  )

  for (label in names(cases)) {
    case <- cases[[label]]
    out <- ot_emd(case$M, case$p, case$q)
    oracle_scale <- 100L
    oracle <- lpSolve::lp.transport(
      case$M,
      direction = "min",
      row.signs = rep("=", length(case$p)),
      row.rhs = as.integer(round(oracle_scale * case$p)),
      col.signs = rep("=", length(case$q)),
      col.rhs = as.integer(round(oracle_scale * case$q))
    )
    expect_identical(oracle$status, 0L, info = label)
    expect_equal(
      out$ot_dist,
      sum(case$M * oracle$solution) / oracle_scale,
      tolerance = 1e-10,
      info = label
    )
    .expect_transport_certificate(out, case$p, case$q)
  }
})

.make_transport_campaign_case <- function(case_id) {
  seed <- 410000L + case_id
  RNGkind("Mersenne-Twister", "Inversion", "Rejection")
  set.seed(seed)
  mode <- c(
    "rectangular", "ties", "zero_weights", "sparse_support",
    "one_point_source", "one_point_target", "degenerate"
  )[[1L + ((case_id - 1L) %% 7L)]]
  ns <- 1L + (case_id %% 4L)
  nt <- 1L + ((case_id * 3L) %% 4L)
  if (identical(mode, "one_point_source")) ns <- 1L
  if (identical(mode, "one_point_target")) nt <- 1L
  total <- 24L
  source_counts <- tabulate(sample.int(ns, total, replace = TRUE), nbins = ns)
  target_counts <- tabulate(sample.int(nt, total, replace = TRUE), nbins = nt)
  if (identical(mode, "zero_weights") || identical(mode, "sparse_support")) {
    if (ns > 1L) {
      source_counts[[ns]] <- 0L
      source_counts[[1L]] <- source_counts[[1L]] + total - sum(source_counts)
    }
    if (nt > 1L) {
      target_counts[[1L]] <- 0L
      target_counts[[nt]] <- target_counts[[nt]] + total - sum(target_counts)
    }
  }
  M <- switch(
    mode,
    ties = matrix(sample(0:2, ns * nt, replace = TRUE), ns, nt),
    degenerate = matrix(1, ns, nt),
    sparse_support = matrix(sample(c(0, 0, 1, 25), ns * nt, replace = TRUE), ns, nt),
    matrix(round(runif(ns * nt, -3, 5), 3), ns, nt)
  )
  list(
    family = "transport_simplex",
    id = case_id,
    seed = seed,
    mode = mode,
    M = M,
    p = source_counts / total,
    q = target_counts / total,
    source_counts = source_counts,
    target_counts = target_counts,
    total = total
  )
}

.transport_campaign_checks <- function(out, case, oracle_value) {
  reduced <- case$M - outer(
    as.numeric(out$source_potential),
    as.numeric(out$target_potential),
    "+"
  )
  primal <- sum(case$M * out$plan)
  dual <- sum(case$p * out$source_potential) +
    sum(case$q * out$target_potential)
  c(
    optimal = identical(out$termination_reason, "optimal") && isTRUE(out$lp_ok),
    converged = identical(out$status, "converged") && isTRUE(out$converged),
    finite_nonnegative = all(is.finite(out$plan)) && all(out$plan >= -1e-12),
    rows = max(abs(rowSums(out$plan) - case$p)) <= out$feasibility_tolerance,
    cols = max(abs(colSums(out$plan) - case$q)) <= out$feasibility_tolerance,
    reduced_cost = min(reduced) >= -out$reduced_cost_tolerance,
    primal = abs(primal - out$primal_objective) <= out$duality_gap_tolerance,
    dual = abs(dual - out$dual_objective) <= out$duality_gap_tolerance,
    gap = abs(primal - dual) <= out$duality_gap_tolerance,
    oracle = abs(primal - oracle_value) <= 1e-9
  )
}

.transport_failure_fixture <- function(case, checks) {
  root <- Sys.getenv("RFUGW_COUNTEREXAMPLE_DIR", unset = tempdir())
  dir.create(root, recursive = TRUE, showWarnings = FALSE)
  path <- file.path(root, sprintf("transport-seed-%d.rds", case$seed))
  saveRDS(c(case, list(failed_checks = names(checks)[!checks])), path)
  normalizePath(path, winslash = "/", mustWork = FALSE)
}

test_that("hundreds of seeded exact-transport cases have independent certificates", {
  skip_if_not_installed("lpSolve")
  n_cases <- switch(
    Sys.getenv("RFUGW_TRUST_SCOPE", unset = "pr"),
    nightly = 1000L,
    release = 2000L,
    240L
  )
  replay_seed <- suppressWarnings(as.integer(Sys.getenv("RFUGW_REPLAY_SEED", "")))
  case_ids <- seq_len(n_cases)
  if (length(replay_seed) == 1L && is.finite(replay_seed)) {
    case_ids <- replay_seed - 410000L
  }
  for (case_id in case_ids) {
    case <- .make_transport_campaign_case(case_id)
    oracle <- lpSolve::lp.transport(
      case$M,
      direction = "min",
      row.signs = rep("=", length(case$p)),
      row.rhs = case$source_counts,
      col.signs = rep("=", length(case$q)),
      col.rhs = case$target_counts,
      integers = integer(0)
    )
    expect_identical(oracle$status, 0L, info = sprintf("seed=%d", case$seed))
    oracle_value <- sum(case$M * oracle$solution) / case$total
    out <- ot_emd(case$M, case$p, case$q, max_iter = 20000L, tol = 1e-11)
    checks <- .transport_campaign_checks(out, case, oracle_value)
    if (!all(checks)) {
      fixture <- .transport_failure_fixture(case, checks)
      expect_true(
        all(checks),
        info = sprintf(
          paste0(
            "mode=%s seed=%d dims=%dx%d failed=%s fixture=%s ",
            "replay=Rscript tools/numerical-trust/replay-case.R %s"
          ),
          case$mode, case$seed, nrow(case$M), ncol(case$M),
          paste(names(checks)[!checks], collapse = ","), fixture, fixture
        )
      )
    }
  }
  succeed()
})
