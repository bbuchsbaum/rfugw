.fgw_precision_case <- function(
    n,
    precision,
    tol = 1e-6,
    inner_tol = 1e-6,
    sinkhorn_method = "scaling",
    solver = "PGD",
    max_iter = 1L,
    sinkhorn_max_iter = 20L) {
  x <- seq_len(n) / n
  C <- abs(outer(x, x, "-"))
  M <- outer(x, rev(x), function(a, b) (a - b)^2)
  fgw_entropic(
    M, C, C,
    epsilon = 0.1,
    max_iter = max_iter,
    tol = tol,
    sinkhorn_max_iter = sinkhorn_max_iter,
    sinkhorn_tol = inner_tol,
    precision = precision,
    sinkhorn_method = sinkhorn_method,
    solver = solver,
    symmetric = TRUE
  )
}

test_that("mixed tolerance boundaries promote instead of silently flooring", {
  floor <- rfugw:::.precision_tolerance_floor
  below <- floor * (1 - 1e-6)
  above <- floor * (1 + 1e-6)
  cases <- list(
    outer_below = c(tol = below, inner = above),
    outer_at = c(tol = floor, inner = above),
    outer_above = c(tol = above, inner = above),
    inner_below = c(tol = above, inner = below),
    inner_at = c(tol = above, inner = floor),
    inner_above = c(tol = above, inner = above)
  )

  for (label in names(cases)) {
    value <- cases[[label]]
    out <- .fgw_precision_case(
      3L, "mixed", tol = value[["tol"]], inner_tol = value[["inner"]]
    )
    promoted <- grepl("below", label, fixed = TRUE)
    expect_identical(out$requested_tol, unname(value[["tol"]]), info = label)
    expect_identical(out$effective_tol, unname(value[["tol"]]), info = label)
    expect_identical(out$requested_inner_tol, unname(value[["inner"]]), info = label)
    expect_identical(out$effective_inner_tol, unname(value[["inner"]]), info = label)
    expect_identical(
      out$effective_precision,
      if (promoted) "strict_double" else "mixed",
      info = label
    )
    expect_identical(out$used_float_inner, !promoted, info = label)
    expect_identical(out$automatic_backend_transition, promoted, info = label)
  }
})

test_that("strict double blocks the hidden large-problem float dispatch", {
  strict <- .fgw_precision_case(32L, "strict_double")
  expect_identical(strict$requested_precision, "strict_double")
  expect_identical(strict$effective_precision, "strict_double")
  expect_identical(strict$compute_precision, "double")
  expect_false(strict$used_float_inner)
  expect_identical(strict$backend, "cpp_strict_double")
  expect_identical(strict$backend_transition, "none")
})

test_that("double acceleration threshold and transition are explicit", {
  below <- .fgw_precision_case(31L, "double")
  at <- .fgw_precision_case(32L, "double")
  tight <- .fgw_precision_case(32L, "double", tol = 1e-9, inner_tol = 1e-9)
  log_path <- .fgw_precision_case(32L, "double", sinkhorn_method = "log")
  ppa_path <- .fgw_precision_case(32L, "double", solver = "PPA")

  expect_identical(below$effective_precision, "double")
  expect_false(below$used_float_inner)
  expect_identical(at$effective_precision, "mixed_accelerated")
  expect_true(at$used_float_inner)
  expect_identical(at$backend_transition, "double_to_mixed_accelerated")
  expect_identical(tight$effective_precision, "double")
  expect_false(tight$used_float_inner)
  expect_identical(log_path$effective_precision, "double")
  expect_false(log_path$used_float_inner)
  expect_identical(ppa_path$effective_precision, "double")
  expect_false(ppa_path$used_float_inner)
})

test_that("FUGW tight mixed requests are promoted and reported", {
  C1 <- matrix(c(0, 1, 2, 1, 0, 1, 2, 1, 0), 3, 3, byrow = TRUE)
  C2 <- matrix(c(0, 2, 1, 2, 0, 1, 1, 1, 0), 3, 3, byrow = TRUE)
  p <- c(0.2, 0.3, 0.5)
  q <- c(0.4, 0.1, 0.5)
  out <- fugw_kl(
    C1, C2, p, q,
    max_iter = 1L,
    tol = 1e-9,
    max_iter_ot = 10L,
    tol_ot = 1e-9,
    precision = "mixed"
  )
  expect_identical(out$requested_precision, "mixed")
  expect_identical(out$effective_precision, "strict_double")
  expect_identical(out$compute_precision, "double")
  expect_false(out$used_float_inner)
  expect_identical(out$requested_tol, 1e-9)
  expect_identical(out$effective_tol, 1e-9)
  expect_identical(out$backend_transition, "mixed_to_strict_double_for_tight_tolerance")
})

test_that("termination reasons distinguish tolerance, budget, stagnation, and failure", {
  tolerance <- .fgw_precision_case(
    3L, "strict_double", tol = 1e-6, inner_tol = 1e-6,
    max_iter = 10L, sinkhorn_max_iter = 50L
  )
  budget <- .fgw_precision_case(
    3L, "strict_double", tol = 1e-15, inner_tol = 1e-6,
    max_iter = 1L, sinkhorn_max_iter = 50L
  )
  failure <- .fgw_precision_case(
    3L, "strict_double", tol = 1, inner_tol = 1e-12,
    max_iter = 10L, sinkhorn_max_iter = 1L
  )
  stagnation <- rfugw:::.termination_reason_from_result(
    list(converged = FALSE, status = "max_iter", iterations = 2L, residual = 1e-4),
    max_iter = 10L
  )

  expect_identical(tolerance$termination_reason, "tolerance")
  expect_identical(budget$termination_reason, "max_iter")
  expect_identical(stagnation, "stagnation")
  expect_identical(failure$termination_reason, "inner_failure")
  expect_false(failure$converged)
})

test_that("auto Sinkhorn uses the documented dynamic-range dispatch", {
  auto <- .fgw_precision_case(
    5L, "strict_double", sinkhorn_method = "auto",
    tol = 1e-7, inner_tol = 1e-8, max_iter = 20L, sinkhorn_max_iter = 100L
  )
  scaling_path <- .fgw_precision_case(
    5L, "strict_double", sinkhorn_method = "scaling",
    tol = 1e-7, inner_tol = 1e-8, max_iter = 20L, sinkhorn_max_iter = 100L
  )
  expect_identical(auto$requested_sinkhorn_method, "auto")
  expect_identical(auto$effective_sinkhorn_method, "scaling")
  expect_identical(auto$sinkhorn_backend_transition, "auto_to_scaling")
  expect_identical(
    auto$sinkhorn_dispatch_reason,
    "dynamic_range_within_scaling_threshold"
  )
  expect_lte(auto$sinkhorn_dynamic_range, auto$sinkhorn_scaling_threshold)
  expect_identical(scaling_path$sinkhorn_backend_transition, "none")
  expect_equal(auto$plan, scaling_path$plan, tolerance = 1e-12)
  expect_equal(auto$fgw_dist, scaling_path$fgw_dist, tolerance = 1e-12)
})

test_that("log-domain FGW is invariant to balanced feature-cost shifts", {
  set.seed(601)
  n <- 5L
  x <- matrix(rnorm(n * 2L), n, 2L)
  C <- as.matrix(dist(x))
  M <- matrix(runif(n * n, -3, 2), n, n)
  alpha <- 0.4
  shift <- 1e5
  args <- list(
    C1 = C, C2 = C, alpha = alpha, epsilon = 0.05,
    max_iter = 40L, tol = 1e-8, sinkhorn_max_iter = 300L,
    sinkhorn_tol = 1e-9, sinkhorn_method = "log",
    precision = "strict_double", symmetric = TRUE
  )
  base <- do.call(fgw_entropic, c(list(M = M), args))
  shifted <- do.call(fgw_entropic, c(list(M = M + shift), args))
  expect_equal(shifted$plan, base$plan, tolerance = 2e-8)
  expect_equal(
    shifted$fgw_dist - base$fgw_dist,
    (1 - alpha) * shift,
    tolerance = 2e-6
  )
})

test_that("dynamic-range degradation is finite or truthfully non-converged", {
  set.seed(602)
  n <- 4L
  C1 <- matrix(runif(n * n, -2e3, 3e3), n, n)
  C2 <- matrix(runif(n * n, -1e3, 4e3), n, n)
  M <- matrix(runif(n * n, -1e6, 1e6), n, n)
  run <- function(method) fgw_entropic(
    M, C1, C2, epsilon = 1e-5, max_iter = 3L, tol = 1e-9,
    sinkhorn_max_iter = 20L, sinkhorn_tol = 1e-10,
    sinkhorn_method = method, precision = "strict_double", symmetric = FALSE
  )
  auto <- run("auto")
  for (out in list(auto = auto)) {
    expect_true(all(is.finite(out$plan)))
    expect_true(all(out$plan >= 0))
    if (isTRUE(out$converged)) {
      expect_true(out$feasible)
      expect_true(out$objective_consistent)
    } else {
      expect_true(out$status %in% c("inner_failure", "max_iter", "infeasible"))
      expect_true(is.list(out$warning_payload))
    }
  }
  expect_identical(auto$effective_sinkhorn_method, "log")
  expect_identical(
    auto$sinkhorn_dispatch_reason,
    "dynamic_range_exceeds_scaling_threshold"
  )
  expect_error(run("scaling"), "outside its certified regime")
})

test_that("POT aliases preserve auto and strict-double provenance", {
  C <- matrix(c(0, 1, 1, 0), 2, 2)
  alias <- entropic_gromov_wasserstein(
    C, C, epsilon = 0.1, max_iter = 10L,
    sinkhorn_method = "auto", precision = "strict_double"
  )
  expect_identical(alias$requested_sinkhorn_method, "auto")
  expect_identical(alias$effective_sinkhorn_method, "scaling")
  expect_identical(alias$requested_precision, "strict_double")
  expect_identical(alias$compute_precision, "double")
})
