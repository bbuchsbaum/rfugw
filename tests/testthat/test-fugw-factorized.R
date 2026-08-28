factorized_fugw_fixture <- function(seed = 1204L, nx = 5L, ny = 6L) {
  set.seed(seed)
  x <- matrix(rnorm(nx * 2L), nx, 2L)
  y <- matrix(rnorm(ny * 2L), ny, 2L)
  source_features <- matrix(rnorm(nx * 3L), nx, 3L)
  target_features <- matrix(rnorm(ny * 3L), ny, 3L)
  wx <- runif(nx)
  wy <- runif(ny)
  wx <- wx / sum(wx)
  wy <- wy / sum(wy)
  list(
    x = x,
    y = y,
    source_features = source_features,
    target_features = target_features,
    Cx = sqeuclidean_cost(x),
    Cy = sqeuclidean_cost(y),
    M = sqeuclidean_cost(source_features, target_features),
    wx = wx,
    wy = wy
  )
}

test_that("factorized TI-UOT matches the dense certified oracle", {
  set.seed(9102)
  source <- matrix(rnorm(24), 6, 4)
  target <- matrix(rnorm(20), 5, 4)
  cost <- sqeuclidean_cost(source, target)
  dense_cost <- cost_block(cost)
  p <- runif(6)
  q <- runif(5)
  p[c(2, 6)] <- 0
  q[4] <- 0
  args <- list(
    p = p, q = q, epsilon = 0.35, rho = c(1.7, 4.2),
    max_iter = 10000L, tol = 1e-11
  )
  implicit <- do.call(
    ot_sinkhorn_unbalanced_ti,
    c(list(cost = cost, plan = "operator", block_size = 3L), args)
  )
  dense <- do.call(
    ot_sinkhorn_unbalanced_ti,
    c(list(cost = dense_cost, plan = "dense"), args)
  )
  implicit_dense <- transport_plan_materialize(rfugw_plan(implicit))

  expect_true(implicit$converged)
  expect_true(implicit$inner_converged)
  expect_equal(implicit_dense, rfugw_plan(dense), tolerance = 2e-11)
  expect_equal(
    implicit$regularized_objective,
    dense$regularized_objective,
    tolerance = 2e-11
  )
  expect_equal(implicit$primal_dual_gap, dense$primal_dual_gap,
               tolerance = 2e-11)
  expect_equal(implicit$source_marginal, rowSums(implicit_dense),
               tolerance = 1e-13)
  expect_equal(implicit$target_marginal, colSums(implicit_dense),
               tolerance = 1e-13)

  source_moment <- matrix(rnorm(12), 6, 2)
  target_moment <- matrix(rnorm(15), 5, 3)
  state <- rfugw:::.factorized_plan_state(implicit$implicit_plan)
  combined <- rfugw:::.factorized_plan_stats_moments(
    state, source_moment, target_moment
  )
  separate_stats <- rfugw:::.factorized_plan_stats(state)
  separate_moments <- rfugw:::.factorized_structure_moments(
    implicit$implicit_plan, source_moment, target_moment
  )
  expect_equal(combined$stats, separate_stats, tolerance = 2e-13)
  expect_equal(combined$moments, separate_moments, tolerance = 2e-13)
  expect_equal(
    combined$workspace_elements,
    3 * 3 + 3 * 3 + 2 * 3 + 2 * (6 + 5)
  )
  expect_true(implicit$source_marginal[[2L]] == 0)
  expect_true(implicit$target_marginal[[4L]] == 0)

  target_values <- matrix(rnorm(15), 5, 3)
  source_values <- matrix(rnorm(12), 6, 2)
  expect_equal(
    transport_plan_apply(rfugw_plan(implicit), target_values),
    implicit_dense %*% target_values,
    tolerance = 2e-13
  )
  expect_equal(
    transport_plan_adjoint(rfugw_plan(implicit), source_values),
    t(implicit_dense) %*% source_values,
    tolerance = 2e-13
  )
  restored <- unserialize(serialize(rfugw_plan(implicit), NULL))
  expect_equal(
    transport_plan_apply(restored, target_values),
    implicit_dense %*% target_values,
    tolerance = 2e-13
  )
})

test_that("factorized FUGW has strict dense-oracle parity", {
  fixture <- factorized_fugw_fixture()
  controls <- list(
    wx = fixture$wx,
    wy = fixture$wy,
    reg_marginals = c(4, 7),
    epsilon = 0.25,
    alpha = 0.55,
    max_iter = 60L,
    tol = 1e-9,
    max_iter_ot = 10000L,
    tol_ot = 1e-11,
    rescale_plan = TRUE
  )
  dense <- do.call(
    fugw_kl,
    c(
      list(
        Cx = cost_block(fixture$Cx),
        Cy = cost_block(fixture$Cy),
        M = cost_block(fixture$M),
        precision = "strict_double"
      ),
      controls
    )
  )
  implicit <- do.call(
    fugw_factorized,
    c(
      list(Cx = fixture$Cx, Cy = fixture$Cy, M = fixture$M,
           block_size = 3L),
      controls
    )
  )
  sample <- transport_plan_materialize(implicit$plans$sample)
  feature <- transport_plan_materialize(implicit$plans$feature)

  expect_true(implicit$converged)
  expect_identical(implicit$status, "converged_stationary")
  expect_true(implicit$inner_uot_certified)
  expect_true(implicit$certificate$outer_stationarity$certified)
  expect_true(implicit$certificate$geometry$exact)
  expect_identical(implicit$certificate$global_optimality, "not_claimed")
  expect_equal(implicit$support_tail_bound, 0)
  expect_true(implicit$outer_plan_residual_exact)
  expect_true(tail(implicit$plan_update_audited_trace, 1L))
  expect_true(
    implicit$runtime_provenance$memory_contract$combined_stats_moments
  )
  expect_true(
    implicit$runtime_provenance$memory_contract$final_kkt_full_support
  )
  expect_identical(
    implicit$runtime_provenance$memory_contract$factor_tile_copy_policy,
    "strided_blas_no_factor_tile_copy"
  )
  expect_true(
    implicit$runtime_provenance$screening_contract$
      stationarity_requires_exact_plan_update
  )
  expect_equal(implicit$fugw_cost, dense$fugw_cost, tolerance = 2e-10)
  expect_equal(
    implicit$objective_decomposition$structure_unweighted,
    dense$objective_decomposition$structure_unweighted,
    tolerance = 2e-9
  )
  expect_equal(
    implicit$objective_decomposition$feature_unweighted,
    dense$objective_decomposition$feature_unweighted,
    tolerance = 2e-9
  )
  expect_equal(
    implicit$objective_decomposition$regularization,
    dense$objective_decomposition$regularization,
    tolerance = 2e-9
  )
  expect_equal(sample, dense$pi_samp, tolerance = 2e-9)
  expect_equal(feature, dense$pi_feat, tolerance = 2e-9)
  expect_lt(implicit$equal_mass_residual, 1e-9)

  source_pair <- rfugw:::.cost_factor_pair(fixture$Cx)
  target_pair <- rfugw:::.cost_factor_pair(fixture$Cy)
  dense_sample_h <- crossprod(
    source_pair$right, dense$pi_samp %*% target_pair$right
  )
  dense_feature_h <- crossprod(
    source_pair$right, dense$pi_feat %*% target_pair$right
  )
  expect_equal(implicit$moments$sample$H, dense_sample_h, tolerance = 3e-9)
  expect_equal(implicit$moments$feature$H, dense_feature_h, tolerance = 3e-9)

  heldout_target <- matrix(rnorm(ncol(feature) * 4L), ncol(feature), 4L)
  heldout_source <- matrix(rnorm(nrow(sample) * 3L), nrow(sample), 3L)
  expect_equal(
    transport_plan_apply(implicit$plans$sample, heldout_target),
    dense$pi_samp %*% heldout_target,
    tolerance = 3e-9
  )
  expect_equal(
    transport_plan_adjoint(implicit$plans$sample, heldout_source),
    t(dense$pi_samp) %*% heldout_source,
    tolerance = 3e-9
  )
})

test_that("factorized FUGW is equivariant to independent node permutations", {
  fixture <- factorized_fugw_fixture(seed = 735L, nx = 4L, ny = 5L)
  source_permutation <- c(3, 1, 4, 2)
  target_permutation <- c(5, 2, 1, 4, 3)
  controls <- list(
    reg_marginals = c(3, 6), epsilon = 0.3, alpha = 0.4,
    max_iter = 50L, tol = 2e-9, max_iter_ot = 8000L,
    tol_ot = 1e-11, block_size = 2L
  )
  original <- do.call(
    fugw_factorized,
    c(
      list(
        Cx = fixture$Cx, Cy = fixture$Cy, M = fixture$M,
        wx = fixture$wx, wy = fixture$wy
      ),
      controls
    )
  )
  permuted <- do.call(
    fugw_factorized,
    c(
      list(
        Cx = sqeuclidean_cost(fixture$x[source_permutation, , drop = FALSE]),
        Cy = sqeuclidean_cost(fixture$y[target_permutation, , drop = FALSE]),
        M = sqeuclidean_cost(
          fixture$source_features[source_permutation, , drop = FALSE],
          fixture$target_features[target_permutation, , drop = FALSE]
        ),
        wx = fixture$wx[source_permutation],
        wy = fixture$wy[target_permutation]
      ),
      controls
    )
  )
  original_plan <- transport_plan_materialize(original$plans$sample)
  permuted_plan <- transport_plan_materialize(permuted$plans$sample)
  restored_plan <- permuted_plan[
    order(source_permutation), order(target_permutation), drop = FALSE
  ]

  expect_true(original$converged)
  expect_true(permuted$converged)
  expect_equal(permuted$fugw_cost, original$fugw_cost, tolerance = 2e-9)
  expect_equal(restored_plan, original_plan, tolerance = 3e-8)
})

test_that("asymmetric factors and zero-weight nodes retain dense parity", {
  set.seed(736)
  nx <- 4L
  ny <- 5L
  source_left <- matrix(runif(nx * 3L), nx, 3L)
  source_right <- matrix(runif(nx * 3L), nx, 3L)
  target_left <- matrix(runif(ny * 2L), ny, 2L)
  target_right <- matrix(runif(ny * 2L), ny, 2L)
  Cx <- factorized_cost(source_left, source_right)
  Cy <- factorized_cost(target_left, target_right)
  M <- factorized_cost(
    matrix(runif(nx * 2L), nx, 2L),
    matrix(runif(ny * 2L), ny, 2L),
    row = runif(nx), column = runif(ny)
  )
  wx <- c(0.2, 0, 0.5, 0.3)
  wy <- c(0, 0.1, 0.2, 0.3, 0.4)
  controls <- list(
    wx = wx, wy = wy, reg_marginals = c(2, 8), epsilon = 0.4,
    alpha = 0.3, max_iter = 50L, tol = 2e-9,
    max_iter_ot = 10000L, tol_ot = 1e-11, rescale_plan = TRUE
  )
  dense <- do.call(
    fugw_kl,
    c(list(Cx = cost_block(Cx), Cy = cost_block(Cy), M = cost_block(M),
           precision = "strict_double"), controls)
  )
  implicit <- do.call(
    fugw_factorized,
    c(list(Cx = Cx, Cy = Cy, M = M, block_size = 2L), controls)
  )
  plan <- transport_plan_materialize(implicit$plans$sample)

  expect_true(implicit$converged)
  expect_equal(plan, dense$pi_samp, tolerance = 3e-8)
  expect_equal(implicit$fugw_cost, dense$fugw_cost, tolerance = 3e-9)
  expect_equal(plan[2, ], rep(0, ny), tolerance = 0)
  expect_equal(plan[, 1], rep(0, nx), tolerance = 0)
})

test_that("matrix-free admission is explicit and callback-only costs fail closed", {
  n <- 1200L
  coordinate <- cbind(seq_len(n) / n, (seq_len(n) %% 13L) / 13)
  cost <- sqeuclidean_cost(coordinate)
  cost$block <- function(...) stop("dense/block evaluation was invoked")
  fit <- ot_sinkhorn_unbalanced_ti(
    cost,
    epsilon = 0.2,
    rho = 2,
    max_iter = 1L,
    tol = 1e-15,
    plan = "operator",
    block_size = 48L
  )
  dense_bytes <- as.double(n) * as.double(n) * 8
  memory <- fit$runtime_provenance$memory_contract

  expect_false(fit$converged)
  expect_false(memory$full_matrix_allocated)
  expect_equal(memory$max_tile_elements, 48^2)
  expect_equal(memory$full_matrix_elements, as.double(n)^2)
  expect_lt(memory$max_tile_elements, memory$full_matrix_elements / 100)
  expect_lt(as.numeric(object.size(fit)), dense_bytes / 2)
  expect_identical(
    transport_plan_representation(rfugw_plan(fit))$representation,
    "implicit_operator"
  )

  callback <- block_cost(3, 3, function(rows, columns) {
    outer(rows, columns, function(i, j) abs(i - j))
  })
  expect_error(
    fugw_factorized(callback, callback, max_iter = 1),
    "no affine-bilinear factors"
  )
  expect_error(
    ot_sinkhorn_unbalanced_ti(cost, plan = "sparse"),
    "cannot return a sparse plan"
  )
  expect_error(
    ot_sinkhorn_unbalanced_ti(cost, n_source = n + 0.5, max_iter = 1),
    "positive integer"
  )
  expect_error(
    ot_sinkhorn_unbalanced_ti(cost, n_target = n + 1L, max_iter = 1),
    "conflicts"
  )
})

test_that("unreferenced approximate geometry stays separate and uncertified", {
  x <- matrix(c(0, 0, 1, 0, 0, 1), 3, 2, byrow = TRUE)
  y <- matrix(c(0, 0, 1, 0, 1, 1, 0, 1), 4, 2, byrow = TRUE)
  source <- sqeuclidean_cost(x)
  target <- sqeuclidean_cost(y)
  source$provenance$exact <- FALSE
  source$provenance$relative_error <- 0.018
  fit <- fugw_factorized(
    source,
    target,
    epsilon = 0.2,
    max_iter = 30L,
    tol = 1e-8,
    max_iter_ot = 5000L,
    tol_ot = 1e-10,
    block_size = 2L
  )

  expect_true(fit$converged)
  expect_identical(fit$status, "converged_uncertified_geometry")
  expect_false(fit$geometry_exact)
  expect_false(fit$geometry_certified)
  expect_identical(
    fit$certificate$geometry$source$status,
    "approximation_unverified_missing_reference"
  )
  expect_equal(fit$geometry_relative_error, 0.018)
  expect_false(fit$certificate$geometry$source$relative_error_verified)
  expect_true(fit$certificate$outer_stationarity$certified)
  expect_true(fit$inner_uot_certified)
  expect_equal(fit$support_tail_bound, 0)
  expect_identical(fit$global_optimality, "not_claimed")
})
