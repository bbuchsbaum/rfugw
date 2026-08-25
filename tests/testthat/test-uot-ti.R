ti_gkl_reference <- function(x, reference) {
  if (any(x > 0 & reference == 0)) return(Inf)
  positive <- x > 0
  sum(x[positive] * log(x[positive] / reference[positive])) -
    sum(x) + sum(reference)
}

ti_direct_reference <- function(cost, p, q, epsilon, rho) {
  n_source <- nrow(cost)
  n_target <- ncol(cost)
  support <- which(tcrossprod(p, q) > 0, arr.ind = TRUE)
  product <- p[support[, 1L]] * q[support[, 2L]]
  objective <- function(log_weight) {
    weight <- exp(pmin(log_weight, 700))
    plan <- matrix(0, n_source, n_target)
    plan[support] <- weight
    transport <- sum(plan * cost)
    plan_kl <- sum(weight * (log_weight - log(product))) - sum(weight) +
      sum(p) * sum(q)
    transport + rho[[1L]] * ti_gkl_reference(rowSums(plan), p) +
      rho[[2L]] * ti_gkl_reference(colSums(plan), q) + epsilon * plan_kl
  }
  gradient <- function(log_weight) {
    weight <- exp(pmin(log_weight, 700))
    plan <- matrix(0, n_source, n_target)
    plan[support] <- weight
    source_mass <- rowSums(plan)
    target_mass <- colSums(plan)
    derivative <- cost[support] +
      rho[[1L]] * log(source_mass[support[, 1L]] /
        p[support[, 1L]]) +
      rho[[2L]] * log(target_mass[support[, 2L]] /
        q[support[, 2L]]) +
      epsilon * (log_weight - log(product))
    weight * derivative
  }
  start <- log(product / max(1, sum(product)))
  fit <- nlminb(
    start, objective, gradient,
    lower = rep(-745, length(start)),
    upper = rep(700, length(start)),
    control = list(iter.max = 50000L, eval.max = 100000L,
                   rel.tol = 1e-13, abs.tol = 1e-13)
  )
  plan <- matrix(0, n_source, n_target)
  plan[support] <- exp(fit$par)
  list(
    plan = plan,
    objective = objective(fit$par),
    convergence = fit$convergence,
    gradient = max(abs(gradient(fit$par)))
  )
}

test_that("TI-UOT matches the analytic one-cell finite-measure oracle", {
  cost <- matrix(0.7, 1, 1)
  p <- 2
  q <- 5
  rho <- c(1.3, 2.1)
  epsilon <- 0.4
  expected <- exp(
    ((rho[[1L]] + epsilon) * log(p) +
       (rho[[2L]] + epsilon) * log(q) - cost[[1L]]) /
      (sum(rho) + epsilon)
  )
  fit <- ot_sinkhorn_unbalanced_ti(
    cost, p, q, epsilon, rho, max_iter = 10000L, tol = 1e-12,
    plan = "dense"
  )
  expect_true(fit$converged)
  expect_equal(rfugw_plan(fit)[[1L]], expected, tolerance = 1e-9)
  expect_lt(abs(fit$primal_dual_gap), 1e-9)
  expect_lt(fit$kkt_residual, fit$kkt_tolerance)
  expect_lt(fit$fixed_point_residual, fit$fixed_point_tolerance)
})

test_that("direct dense optimization covers ordinary and adversarial regimes", {
  cases <- list(
    ordinary_rectangular = list(
      cost = matrix(c(0, .5, 1, .2, .7, .1), 2, 3),
      p = c(.4, 1.1), q = c(.2, .7, 1.4), epsilon = .4,
      rho = c(1.5, 2.2)
    ),
    zero_mass_entries = list(
      cost = matrix(c(0, 4, 1, 2, .3, 3), 3, 2),
      p = c(.2, 0, 1.8), q = c(0, 2.7), epsilon = .3,
      rho = c(3, .8)
    ),
    duplicate_costs = list(
      cost = matrix(c(0, 0, 2, 2, 0, 0), 2, 3),
      p = c(1, 2), q = c(.5, 1, 2), epsilon = .15,
      rho = c(.7, 4)
    ),
    stiff_penalties = list(
      cost = matrix(c(0, .2, 1.5, .1), 2),
      p = c(.1, .9), q = c(1.7, .2), epsilon = .025,
      rho = c(.12, 8)
    ),
    extreme_dynamic_range = list(
      cost = matrix(c(-20, 4, 25, 0), 2),
      p = c(1e-4, 3), q = c(2, 1e-3), epsilon = .08,
      rho = c(2, 5)
    )
  )
  for (case in cases) {
    fit <- do.call(
      ot_sinkhorn_unbalanced_ti,
      c(case, list(max_iter = 20000L, tol = 1e-11, plan = "dense"))
    )
    reference <- ti_direct_reference(
      case$cost, case$p, case$q, case$epsilon, case$rho
    )
    expect_true(reference$convergence %in% c(0L, 1L))
    expect_lt(reference$gradient, 2e-5)
    expect_true(fit$converged)
    expect_equal(
      fit$regularized_objective, reference$objective,
      tolerance = 2e-6, info = deparse(case$cost)
    )
    expect_equal(
      rfugw_plan(fit), reference$plan,
      tolerance = 2e-6, info = deparse(case$cost)
    )
  }
})

test_that("dense, complete sparse, Matrix, and manifoldalign CSR paths agree", {
  set.seed(4201)
  cost <- matrix(runif(30), 6, 5)
  p <- runif(6)
  q <- runif(5)
  args <- list(p = p, q = q, epsilon = .4, rho = c(1.7, .8),
               max_iter = 5000L, tol = 1e-10, plan = "dense")
  dense <- do.call(ot_sinkhorn_unbalanced_ti, c(list(cost = cost), args))
  edges <- data.frame(
    source = rep(seq_len(6), each = 5),
    target = rep(seq_len(5), times = 6),
    cost = as.numeric(t(cost))
  )
  sparse <- do.call(
    ot_sinkhorn_unbalanced_ti,
    c(list(cost = edges, n_source = 6, n_target = 5), args)
  )
  csr <- list(
    row_ptr = as.integer(seq(0, 30, by = 5)),
    col_idx = as.integer(rep(seq_len(5), times = 6)),
    cost = as.numeric(t(cost)),
    n_rows = 6L,
    n_cols = 5L
  )
  csr_fit <- do.call(ot_sinkhorn_unbalanced_ti, c(list(cost = csr), args))
  expect_equal(rfugw_plan(sparse), rfugw_plan(dense), tolerance = 1e-12)
  expect_equal(rfugw_plan(csr_fit), rfugw_plan(dense), tolerance = 1e-12)
  expect_equal(sparse$source_potential, dense$source_potential, tolerance = 1e-12)
  expect_equal(sparse$target_potential, dense$target_potential, tolerance = 1e-12)

  skip_if_not_installed("Matrix")
  positive_cost <- cost + .01
  matrix_fit <- do.call(
    ot_sinkhorn_unbalanced_ti,
    c(list(cost = Matrix::Matrix(positive_cost, sparse = TRUE)), args)
  )
  positive_dense <- do.call(
    ot_sinkhorn_unbalanced_ti,
    c(list(cost = positive_cost), args)
  )
  expect_equal(rfugw_plan(matrix_fit), rfugw_plan(positive_dense), tolerance = 1e-12)
})

test_that("potential gauge shifts leave the plan and objective invariant", {
  cost <- matrix(c(0, .4, .8, .2, 1, .1), 2, 3)
  p <- c(.7, 1.3)
  q <- c(.4, .8, 1.8)
  epsilon <- .25
  fit <- ot_sinkhorn_unbalanced_ti(
    cost, p, q, epsilon, rho = c(2, 3), tol = 1e-11, plan = "dense"
  )
  shift <- 17.3
  shifted <- tcrossprod(p, q) * exp(
    (outer(fit$fbar + shift, rep(1, 3)) +
       outer(rep(1, 2), fit$gbar - shift) - cost) / epsilon
  )
  expect_equal(shifted, rfugw_plan(fit), tolerance = 1e-12)
  expect_true(fit$gauge_invariant)
  expect_lt(fit$gauge_residual, 1e-12)
})

test_that("global max flow diagnoses Hall failures beyond local degree", {
  support <- data.frame(
    source = c(1, 1, 2, 3),
    target = c(1, 2, 3, 3),
    cost = c(0, .1, .2, .3)
  )
  p <- c(.6, .2, .2)
  q <- c(.2, .2, .6)
  fit <- ot_sinkhorn_unbalanced_ti(
    support, p, q, epsilon = .2, rho = 2,
    n_source = 3, n_target = 3, tol = 1e-10, plan = "sparse"
  )
  expect_true(fit$converged)
  expect_true(fit$support_certificate$finite_potential_support)
  expect_length(fit$support_certificate$uncovered_source, 0L)
  expect_length(fit$support_certificate$uncovered_target, 0L)
  expect_false(fit$support_certificate$balanced_marginal_feasible)
  expect_equal(fit$support_certificate$balanced_max_flow, .8, tolerance = 1e-10)
  expect_equal(fit$support_certificate$balanced_flow_deficit, .2, tolerance = 1e-10)

  bad <- subset(support, source != 3)
  expect_error(
    ot_sinkhorn_unbalanced_ti(
      bad, p, q, epsilon = .2, rho = 2, n_source = 3, n_target = 3
    ),
    "finite TI potentials"
  )
})

test_that("operator application, adjoint, sparse extraction, and persistence agree", {
  support <- data.frame(
    source = c(1, 1, 2, 2, 3, 3),
    target = c(1, 2, 2, 3, 1, 3),
    cost = c(0, .7, .1, .5, .8, 0)
  )
  args <- list(
    cost = support, p = c(.2, .5, .7), q = c(.9, .4, .8),
    epsilon = .2, rho = c(1.5, 2.3), n_source = 3, n_target = 3,
    tol = 1e-11, max_iter = 10000L
  )
  operator_fit <- do.call(ot_sinkhorn_unbalanced_ti, c(args, list(plan = "operator")))
  sparse_fit <- do.call(ot_sinkhorn_unbalanced_ti, c(args, list(plan = "sparse")))
  dense <- as.matrix(rfugw_plan(sparse_fit))
  x_target <- matrix(c(1, 2, 3, 4, 5, 6), 3, 2)
  x_source <- matrix(c(6, 5, 4, 3, 2, 1), 3, 2)
  expect_equal(
    transport_plan_apply(rfugw_plan(operator_fit), x_target),
    dense %*% x_target,
    tolerance = 1e-12
  )
  expect_equal(
    transport_plan_adjoint(rfugw_plan(operator_fit), x_source),
    t(dense) %*% x_source,
    tolerance = 1e-12
  )
  expect_equal(
    transport_plan_barycentric(
      rfugw_plan(operator_fit), x_source, "target_to_source"
    ),
    sweep(t(dense) %*% x_source, 1, colSums(dense), "/"),
    tolerance = 1e-12
  )
  expect_identical(
    transport_plan_representation(rfugw_plan(operator_fit))$representation,
    "implicit_operator"
  )
  expect_identical(
    transport_plan_representation(rfugw_plan(sparse_fit))$representation,
    "edge_list"
  )
  restored <- unserialize(serialize(operator_fit, NULL))
  expect_equal(
    transport_plan_apply(rfugw_plan(restored), x_target),
    dense %*% x_target,
    tolerance = 1e-12
  )
})

test_that("uncertified iterations and invalid inputs fail closed", {
  cost <- matrix(c(0, 1, 1, 0), 2)
  short <- ot_sinkhorn_unbalanced_ti(
    cost, epsilon = 1e-3, rho = 10, max_iter = 1L, tol = 1e-15
  )
  expect_false(short$converged)
  expect_false(identical(rfugw_status(short), "converged"))
  expect_error(
    ot_sinkhorn_unbalanced_ti(matrix(c(0, Inf, 1, 0), 2)),
    "finite"
  )
  expect_error(
    ot_sinkhorn_unbalanced_ti(cost, p = c(0, 0)),
    "positive total mass"
  )
  expect_error(
    ot_sinkhorn_unbalanced_ti(
      data.frame(source = c(1, 1), target = c(1, 1), cost = c(0, 1)),
      n_source = 1, n_target = 1
    ),
    "duplicate edges"
  )
})

test_that("large sparse support does not allocate its dense boundary", {
  n <- 10000L
  support <- data.frame(
    source = rep(seq_len(n), each = 2L),
    target = c(seq_len(n), c(seq_len(n - 1L) + 1L, 1L)),
    cost = rep(c(0, .2), n)
  )
  fit <- ot_sinkhorn_unbalanced_ti(
    support, rep(1 / n, n), rep(1 / n, n),
    epsilon = .1, rho = 1, n_source = n, n_target = n,
    max_iter = 1L, tol = 1e-15, plan = "operator"
  )
  expect_equal(fit$support_size, 2L * n)
  expect_lt(as.numeric(object.size(fit)), 20 * 1024^2)
  expect_gt(8 * n * n, 700 * 1024^2)
  expect_identical(
    transport_plan_representation(rfugw_plan(fit))$representation,
    "implicit_operator"
  )
})

test_that("temporary manifoldalign differential fixture remains exact", {
  skip_if_not_installed("manifoldalign")
  set.seed(44)
  cost <- matrix(runif(30), 6, 5)
  p <- runif(6)
  q <- runif(5)
  old <- manifoldalign::uot_ti_sinkhorn_kl(
    cost, p, q, .4, 1.7, .8, max_iter = 5000L, tol = 1e-10
  )
  new <- ot_sinkhorn_unbalanced_ti(
    cost, p, q, .4, c(1.7, .8), max_iter = 5000L, tol = 1e-10,
    plan = "dense"
  )
  old_plan <- tcrossprod(p, q) * exp(
    (outer(as.numeric(old$fbar), rep(1, 5)) +
       outer(rep(1, 6), as.numeric(old$gbar)) - cost) / .4
  )
  expect_equal(rfugw_plan(new), old_plan, tolerance = 1e-12)
  expect_equal(new$f, as.numeric(old$f), tolerance = 1e-12)
  expect_equal(new$g, as.numeric(old$g), tolerance = 1e-12)
})
