sinkhorn_warm_fixture <- function() {
  set.seed(42)
  X <- matrix(rnorm(60), 20, 3)
  Y <- matrix(rnorm(45), 15, 3)
  Z <- rbind(X, Y)
  M <- as.matrix(stats::dist(Z))[seq_len(nrow(X)), nrow(X) + seq_len(nrow(Y))]^2
  list(M = M, p = rep(1 / nrow(X), nrow(X)), q = rep(1 / nrow(Y), nrow(Y)))
}

test_that("balanced Sinkhorn returns canonical reusable dual potentials", {
  d <- sinkhorn_warm_fixture()
  scaling <- ot_sinkhorn(
    d$M, d$p, d$q, epsilon = 0.4, method = "scaling",
    max_iter = 5000L, tol = 1e-9
  )
  log_domain <- ot_sinkhorn(
    d$M, d$p, d$q, epsilon = 0.4, method = "log",
    max_iter = 5000L, tol = 1e-9
  )
  automatic <- ot_sinkhorn(
    d$M, d$p, d$q, epsilon = 0.4, method = "auto",
    max_iter = 5000L, tol = 1e-9
  )

  for (out in list(scaling, log_domain)) {
    expect_true(out$converged)
    expect_length(out$source_potential, length(d$p))
    expect_length(out$target_potential, length(d$q))
    expect_true(all(is.finite(out$source_potential)))
    expect_true(all(is.finite(out$target_potential)))
    expect_equal(sum(d$p * out$source_potential), 0, tolerance = 1e-12)
    reconstructed <- exp(
      (outer(out$source_potential, out$target_potential, "+") - d$M) / 0.4
    )
    expect_equal(unname(reconstructed), out$plan, tolerance = 2e-10)
    expect_identical(out$potential_gauge, "weighted_source_mean_zero")
  }
  expect_equal(scaling$plan, log_domain$plan, tolerance = 2e-8)
  expect_equal(scaling$source_potential, log_domain$source_potential, tolerance = 2e-8)
  expect_equal(scaling$target_potential, log_domain$target_potential, tolerance = 2e-8)
  expect_equal(automatic$plan, scaling$plan, tolerance = 2e-10)
  expect_equal(automatic$source_potential, scaling$source_potential, tolerance = 2e-10)
  expect_equal(automatic$target_potential, scaling$target_potential, tolerance = 2e-10)
})

test_that("dual and plan warm starts preserve the certified answer", {
  d <- sinkhorn_warm_fixture()
  cold <- ot_sinkhorn(
    d$M, d$p, d$q, epsilon = 0.4, method = "scaling",
    max_iter = 5000L, tol = 1e-9
  )
  dual_warm <- ot_sinkhorn(
    d$M, d$p, d$q, epsilon = 0.4, method = "log",
    max_iter = 5000L, tol = 1e-9, init_duals = cold
  )
  plan_warm <- ot_sinkhorn(
    d$M, d$p, d$q, epsilon = 0.4, method = "scaling",
    max_iter = 5000L, tol = 1e-9, init_plan = cold$plan
  )

  expect_true(dual_warm$converged)
  expect_true(plan_warm$converged)
  expect_equal(dual_warm$plan, cold$plan, tolerance = 2e-9)
  expect_equal(plan_warm$plan, cold$plan, tolerance = 2e-9)
  expect_lte(dual_warm$iterations, cold$iterations)
  expect_lte(plan_warm$iterations, cold$iterations)
  expect_identical(dual_warm$initialization, "duals")
  expect_identical(plan_warm$initialization, "plan")

  precedence <- ot_sinkhorn(
    d$M, d$p, d$q, epsilon = 0.4, method = "log",
    max_iter = 5000L, tol = 1e-9,
    init_plan = cold$plan^0.9, init_duals = cold$dual_state
  )
  expect_identical(precedence$initialization, "duals_over_plan")
  expect_identical(precedence$initialization_precedence, "init_duals_over_init_plan")
  expect_equal(precedence$plan, cold$plan, tolerance = 2e-9)
})

test_that("epsilon continuation uses no more total iterations than cold restarts", {
  d <- sinkhorn_warm_fixture()
  epsilon <- c(1, 0.6, 0.35, 0.2)
  baseline <- utils::read.csv(
    trust_test_resource("bench", "sinkhorn-warm-start-baseline.csv"),
    stringsAsFactors = FALSE
  )
  state <- NULL
  warm_iterations <- cold_iterations <- integer(length(epsilon))
  for (i in seq_along(epsilon)) {
    warm <- ot_sinkhorn(
      d$M, d$p, d$q, epsilon = epsilon[[i]], method = "scaling",
      max_iter = 5000L, tol = 1e-9, init_duals = state
    )
    cold <- ot_sinkhorn(
      d$M, d$p, d$q, epsilon = epsilon[[i]], method = "scaling",
      max_iter = 5000L, tol = 1e-9
    )
    expect_true(warm$converged)
    expect_true(cold$converged)
    expect_lte(max(abs(warm$plan - cold$plan)), 2e-8)
    expect_lte(abs(warm$ot_dist - cold$ot_dist), 2e-8)
    warm_iterations[[i]] <- warm$iterations
    cold_iterations[[i]] <- cold$iterations
    state <- warm$dual_state
  }
  expect_equal(epsilon, baseline$epsilon)
  expect_true(all(warm_iterations <= baseline$warm_iteration_ceiling))
  expect_true(all(cold_iterations <= baseline$cold_iteration_ceiling))
  expect_lte(sum(warm_iterations), sum(cold_iterations))
  expect_lt(sum(warm_iterations), sum(cold_iterations))
})

test_that("invalid balanced warm states fail before computation", {
  d <- sinkhorn_warm_fixture()
  cold <- ot_sinkhorn(d$M, d$p, d$q, epsilon = 0.4, method = "log")
  expect_error(
    ot_sinkhorn(d$M, d$p, d$q, init_plan = matrix(1, 2, 2)),
    "shape 20 x 15"
  )
  bad_plan <- cold$plan
  bad_plan[[1L]] <- -1
  expect_error(ot_sinkhorn(d$M, d$p, d$q, init_plan = bad_plan), "nonnegative")
  bad_plan[[1L]] <- NA_real_
  expect_error(ot_sinkhorn(d$M, d$p, d$q, init_plan = bad_plan), "finite")
  expect_error(
    ot_sinkhorn(d$M, d$p, d$q,
                init_duals = list(source = numeric(2), target = numeric(15))),
    "source.*length 20"
  )
  expect_error(
    ot_sinkhorn(d$M, d$p, d$q,
                init_duals = list(source = rep(0, 20), target = c(NA, rep(0, 14)))),
    "potentials must be finite"
  )
  incompatible_state <- cold$dual_state
  incompatible_state$source_support[[1L]] <- FALSE
  expect_error(
    ot_sinkhorn(d$M, d$p, d$q, init_duals = incompatible_state),
    "source support is incompatible"
  )
  expect_error(
    ot_sinkhorn(
      d$M, d$p, d$q, epsilon = 0.4, method = "scaling",
      init_duals = list(
        source = seq(-1000, 1000, length.out = 20),
        target = seq(1000, -1000, length.out = 15)
      )
    ),
    "overflow the scaling backend"
  )

  p_zero <- d$p
  p_zero[[1L]] <- 0
  p_zero <- p_zero / sum(p_zero)
  outside_support <- cold$plan
  expect_error(
    ot_sinkhorn(d$M, p_zero, d$q, init_plan = outside_support),
    "zero outside"
  )
  missing_active <- cold$plan
  missing_active[[2L]] <- 0
  expect_error(
    ot_sinkhorn(d$M, d$p, d$q, init_plan = missing_active),
    "strictly positive"
  )
})

test_that("a bad warm state cannot inherit converged status", {
  d <- sinkhorn_warm_fixture()
  bad <- list(
    source = seq(-100, 100, length.out = length(d$p)),
    target = seq(80, -80, length.out = length(d$q))
  )
  out <- ot_sinkhorn(
    d$M, d$p, d$q, epsilon = 0.4, method = "log",
    max_iter = 1L, tol = 1e-14, init_duals = bad
  )
  expect_false(out$converged)
  expect_false(out$feasible)
  expect_gt(out$residual, out$feasibility_tolerance)
})

test_that("dual state round-trips and powers a public downstream client", {
  d <- sinkhorn_warm_fixture()
  first <- ot_sinkhorn(
    d$M, d$p, d$q, epsilon = 0.5, method = "auto",
    max_iter = 5000L, tol = 1e-9
  )
  state_file <- tempfile(fileext = ".rds")
  saveRDS(first$dual_state, state_file)
  state <- readRDS(state_file)
  second <- ot_sinkhorn(
    d$M, d$p, d$q, epsilon = 0.4, method = "auto",
    max_iter = 5000L, tol = 1e-9, init_duals = state
  )
  expect_true(second$converged)

  client <- new.env(parent = baseenv())
  sys.source(
    trust_test_resource(
      "extdata", "clients", "manifoldalign_ot_procrustes_fixture.R"
    ),
    envir = client
  )
  set.seed(9)
  X <- matrix(rnorm(36), 12, 3)
  rotation <- qr.Q(qr(matrix(rnorm(9), 3, 3)))
  Y <- X %*% rotation + matrix(rnorm(36, sd = 0.02), 12, 3)
  step1 <- client$manifoldalign_ot_procrustes_step(X, Y, epsilon = 0.8)
  step2 <- client$manifoldalign_ot_procrustes_step(
    X, Y, epsilon = 0.65, state = step1$state
  )
  expect_true(step1$result$converged)
  expect_true(step2$result$converged)
  expect_equal(crossprod(step2$rotation), diag(3), tolerance = 1e-10)
  expect_identical(step2$result$initialization, "duals")
})
