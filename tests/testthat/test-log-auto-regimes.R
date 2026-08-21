test_that("genuine log UOT agrees with scaling in its certified regime", {
  M <- matrix(c(0, 0.4, 1.2, 0.3, 0, 0.8, 1.1, 0.5, 0), 3, 3)
  p <- c(0.2, 0.5, 0.3)
  q <- c(0.4, 0.1, 0.5)
  args <- list(
    M = M, p = p, q = q, epsilon = 0.2, rho = c(2, 3),
    max_iter = 1000L, tol = 1e-9
  )
  scaling <- do.call(ot_sinkhorn_unbalanced, c(args, list(method = "scaling")))
  log_path <- do.call(ot_sinkhorn_unbalanced, c(args, list(method = "log")))

  expect_equal(log_path$plan, scaling$plan, tolerance = 2e-8)
  expect_equal(log_path$mass, scaling$mass, tolerance = 2e-8)
  expect_equal(log_path$ot_dist, scaling$ot_dist, tolerance = 2e-8)
  expect_identical(log_path$backend, "r_log")
  expect_identical(log_path$effective_sinkhorn_method, "log")
})

test_that("UOT auto selects genuine log for tiny epsilon and bad scaling", {
  M <- matrix(c(0, 1e3, 2e3, 1e3, 0, 1e3, 2e3, 1e3, 0), 3, 3)
  p <- c(0, 0.3, 0.7)
  q <- c(0.6, 0, 0.4)
  out <- ot_sinkhorn_unbalanced(
    M, p, q, epsilon = 1e-3, rho = 2, method = "auto",
    max_iter = 2000L, tol = 1e-8
  )
  expect_identical(out$effective_sinkhorn_method, "log")
  expect_identical(out$sinkhorn_backend_transition, "auto_to_log")
  expect_identical(
    out$sinkhorn_dispatch_reason,
    "dynamic_range_exceeds_scaling_threshold"
  )
  expect_true(all(is.finite(out$plan)))
  expect_true(all(out$plan >= 0))
  expect_equal(out$plan[1, ], rep(0, 3), tolerance = 0)
  expect_equal(out$plan[, 2], rep(0, 3), tolerance = 0)
  expect_error(
    ot_sinkhorn_unbalanced(M, p, q, epsilon = 1e-3, method = "scaling"),
    "outside its certified regime"
  )
})

test_that("log UOT obeys joint cost-regularization scaling and handles shifts", {
  M <- matrix(c(0, 0.5, 1.5, 0), 2, 2)
  base <- ot_sinkhorn_unbalanced(
    M, epsilon = 0.05, rho = c(1, 2), method = "log",
    max_iter = 1500L, tol = 1e-9
  )
  scaled <- ot_sinkhorn_unbalanced(
    7 * M, epsilon = 0.35, rho = c(7, 14), method = "log",
    max_iter = 1500L, tol = 1e-9
  )
  shifted <- ot_sinkhorn_unbalanced(
    M + 100, epsilon = 0.05, rho = c(1, 2), method = "log",
    max_iter = 1500L, tol = 1e-9
  )
  expect_equal(scaled$plan, base$plan, tolerance = 1e-8)
  expect_equal(scaled$ot_dist, 7 * base$ot_dist, tolerance = 1e-8)
  expect_true(all(is.finite(shifted$plan)))
  expect_lt(shifted$mass, base$mass)
})

test_that("balanced auto dispatch is shift-aware and never returns clipping", {
  M <- matrix(c(0, 1, 2, 0), 2, 2)
  benign <- ot_sinkhorn(M, epsilon = 0.1, method = "auto", max_iter = 500L)
  shifted <- ot_sinkhorn(
    M + 1e5, epsilon = 0.1, method = "auto", max_iter = 500L
  )
  expect_identical(benign$effective_sinkhorn_method, "scaling")
  expect_identical(shifted$effective_sinkhorn_method, "log")
  expect_equal(shifted$plan, benign$plan, tolerance = 1e-8)
  expect_error(
    ot_sinkhorn(M + 1e5, epsilon = 0.1, method = "scaling"),
    "outside its certified regime"
  )
})

test_that("partial auto is explicitly bounded and unsupported log is truthful", {
  C <- matrix(c(0, 1, 1, 0), 2, 2)
  safe <- entropic_partial_gromov_wasserstein(
    C, C, reg = 1, numItermax = 3L, inner_max_iter = 50L,
    method = "auto", log = TRUE
  )
  expect_identical(safe$effective_sinkhorn_method, "scaling")
  expect_identical(
    safe$sinkhorn_dispatch_reason,
    "dynamic_range_within_scaling_threshold"
  )
  expect_error(
    entropic_partial_gromov_wasserstein(
      1e4 * C, 1e4 * C, reg = 1e-3, method = "auto"
    ),
    "genuine log-domain Dykstra"
  )
  expect_error(
    entropic_partial_gromov_wasserstein(C, C, method = "log"),
    "genuine log-domain Dykstra"
  )
})

test_that("UCOOT no longer aliases sinkhorn_log to scaling", {
  X <- matrix(c(0, 1, 2, 3), 2, 2)
  expect_error(
    unbalanced_co_optimal_transport(
      X, X, unbalanced_solver = "sinkhorn_log"
    ),
    "deprecated and unsupported"
  )
})
