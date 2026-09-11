moment_metamorphic_fixture <- function(seed = 20260827L, nx = 4L, ny = 5L) {
  set.seed(seed)
  x <- matrix(rnorm(nx * 2L), nx, 2L)
  y <- matrix(rnorm(ny * 2L), ny, 2L)
  source_features <- matrix(rnorm(nx * 3L), nx, 3L)
  target_features <- matrix(rnorm(ny * 3L), ny, 3L)
  wx <- runif(nx)
  wy <- runif(ny)
  list(
    x = x,
    y = y,
    source_features = source_features,
    target_features = target_features,
    Cx = sqeuclidean_cost(x),
    Cy = sqeuclidean_cost(y),
    M = sqeuclidean_cost(source_features, target_features),
    wx = wx / sum(wx),
    wy = wy / sum(wy)
  )
}

moment_metamorphic_controls <- function(block_size = 2L) {
  list(
    reg_marginals = c(3, 6),
    epsilon = 0.3,
    feature_weight = 0.4,
    structure_weight = 1,
    max_iter = 80L,
    tol = 1e-8,
    max_iter_ot = 10000L,
    tol_ot = 1e-10,
    rescale_plan = TRUE,
    block_size = block_size
  )
}

moment_metamorphic_fit <- function(
    fixture,
    Cx = fixture$Cx,
    Cy = fixture$Cy,
    M = fixture$M,
    block_size = 2L,
    controls = list()) {
  do.call(
    fugw_factorized,
    c(
      list(
        Cx = Cx,
        Cy = Cy,
        M = M,
        wx = fixture$wx,
        wy = fixture$wy
      ),
      utils::modifyList(moment_metamorphic_controls(block_size), controls)
    )
  )
}

moment_rebuild_cost <- function(cost, left, right) {
  provenance <- cost$provenance
  provenance[c(
    "representation", "rank", "exact", "constant_folded_into_row"
  )] <- NULL
  factorized_cost(
    left = left,
    right = right,
    row = cost$factors$row,
    column = cost$factors$column,
    exact = isTRUE(cost$provenance$exact),
    provenance = provenance,
    geometry_reference = cost$geometry_reference,
    geometry_audit = cost$geometry_audit_spec
  )
}

moment_gauge_cost <- function(cost) {
  rank <- ncol(cost$factors$left)
  gauge <- diag(seq(0.75, 1.75, length.out = rank), rank)
  if (rank > 1L) gauge[1L, rank] <- 0.2
  moment_rebuild_cost(
    cost,
    cost$factors$left %*% gauge,
    cost$factors$right %*% solve(t(gauge))
  )
}

moment_pad_cost <- function(cost, padding = 2L) {
  moment_rebuild_cost(
    cost,
    cbind(cost$factors$left, matrix(0, nrow(cost$factors$left), padding)),
    cbind(cost$factors$right, matrix(0, nrow(cost$factors$right), padding))
  )
}

expect_moment_fits_equivalent <- function(left, right, tolerance = 5e-8) {
  expect_true(left$converged)
  expect_true(right$converged)
  expect_true(left$certificate$outer_stationarity$certified)
  expect_true(right$certificate$outer_stationarity$certified)
  expect_equal(left$fugw_cost, right$fugw_cost, tolerance = tolerance)
  expect_equal(
    left$objective_decomposition$structure_unweighted,
    right$objective_decomposition$structure_unweighted,
    tolerance = tolerance
  )
  expect_equal(
    left$objective_decomposition$feature_unweighted,
    right$objective_decomposition$feature_unweighted,
    tolerance = tolerance
  )
  expect_equal(
    left$objective_decomposition$regularization,
    right$objective_decomposition$regularization,
    tolerance = tolerance
  )
  for (role in c("sample", "feature")) {
    expect_equal(
      transport_plan_materialize(left$plans[[role]]),
      transport_plan_materialize(right$plans[[role]]),
      tolerance = tolerance
    )
  }
  invisible(TRUE)
}

test_that("independent node permutations preserve both Moment-FUGW couplings", {
  fixture <- moment_metamorphic_fixture(seed = 1400L)
  source_permutation <- c(3L, 1L, 4L, 2L)
  target_permutation <- c(5L, 2L, 1L, 4L, 3L)
  permuted <- list(
    x = fixture$x[source_permutation, , drop = FALSE],
    y = fixture$y[target_permutation, , drop = FALSE],
    source_features = fixture$source_features[
      source_permutation, , drop = FALSE
    ],
    target_features = fixture$target_features[
      target_permutation, , drop = FALSE
    ],
    wx = fixture$wx[source_permutation],
    wy = fixture$wy[target_permutation]
  )
  permuted$Cx <- sqeuclidean_cost(permuted$x)
  permuted$Cy <- sqeuclidean_cost(permuted$y)
  permuted$M <- sqeuclidean_cost(
    permuted$source_features,
    permuted$target_features
  )

  original <- moment_metamorphic_fit(fixture)
  reordered <- moment_metamorphic_fit(permuted, block_size = 3L)
  expect_true(original$converged)
  expect_true(reordered$converged)
  expect_equal(reordered$fugw_cost, original$fugw_cost, tolerance = 5e-8)
  for (role in c("sample", "feature")) {
    restored <- transport_plan_materialize(reordered$plans[[role]])[
      order(source_permutation), order(target_permutation), drop = FALSE
    ]
    expect_equal(
      restored,
      transport_plan_materialize(original$plans[[role]]),
      tolerance = 5e-8
    )
  }
})

test_that("coordinate isometries preserve both Moment-FUGW couplings", {
  fixture <- moment_metamorphic_fixture(seed = 1401L)
  rotation <- matrix(c(0, -1, 1, 0), 2L, 2L, byrow = TRUE)
  source_geometry <- sweep(fixture$x %*% rotation, 2L, c(12, -7), "+")
  target_geometry <- sweep(fixture$y %*% rotation, 2L, c(-3, 9), "+")
  feature_rotation <- diag(c(-1, 1, -1))[, c(3, 1, 2)]
  source_features <- sweep(
    fixture$source_features %*% feature_rotation,
    2L,
    c(4, -2, 8),
    "+"
  )
  target_features <- sweep(
    fixture$target_features %*% feature_rotation,
    2L,
    c(4, -2, 8),
    "+"
  )
  transformed <- list(
    Cx = sqeuclidean_cost(source_geometry),
    Cy = sqeuclidean_cost(target_geometry),
    M = sqeuclidean_cost(source_features, target_features)
  )

  expect_equal(cost_block(transformed$Cx), cost_block(fixture$Cx),
               tolerance = 2e-13)
  expect_equal(cost_block(transformed$Cy), cost_block(fixture$Cy),
               tolerance = 2e-13)
  expect_equal(cost_block(transformed$M), cost_block(fixture$M),
               tolerance = 2e-13)

  original <- moment_metamorphic_fit(fixture)
  isometric <- moment_metamorphic_fit(
    fixture,
    transformed$Cx,
    transformed$Cy,
    transformed$M,
    block_size = 3L
  )
  expect_moment_fits_equivalent(original, isometric)
})

test_that("factor gauges, zero padding, composition, and tiling are invariant", {
  fixture <- moment_metamorphic_fixture(seed = 1402L)
  gauged <- lapply(fixture[c("Cx", "Cy", "M")], moment_gauge_cost)
  padded <- lapply(fixture[c("Cx", "Cy", "M")], moment_pad_cost)
  composed <- lapply(fixture[c("Cx", "Cy", "M")], function(cost) {
    sum_costs(scaled_cost(cost, 0.25), scaled_cost(cost, 0.75))
  })

  for (role in c("Cx", "Cy", "M")) {
    reference <- cost_block(fixture[[role]])
    expect_equal(cost_block(gauged[[role]]), reference, tolerance = 2e-13)
    expect_equal(cost_block(padded[[role]]), reference, tolerance = 2e-13)
    expect_equal(cost_block(composed[[role]]), reference, tolerance = 2e-13)
    expect_equal(
      cost_provenance(padded[[role]])$rank,
      cost_provenance(fixture[[role]])$rank + 2L
    )
  }

  reference <- moment_metamorphic_fit(fixture, block_size = 1L)
  gauge_fit <- moment_metamorphic_fit(
    fixture, gauged$Cx, gauged$Cy, gauged$M, block_size = 2L
  )
  padded_fit <- moment_metamorphic_fit(
    fixture, padded$Cx, padded$Cy, padded$M, block_size = 3L
  )
  composed_fit <- moment_metamorphic_fit(
    fixture, composed$Cx, composed$Cy, composed$M, block_size = 64L
  )
  expect_moment_fits_equivalent(reference, gauge_fit)
  expect_moment_fits_equivalent(reference, padded_fit)
  expect_moment_fits_equivalent(reference, composed_fit)
})

test_that("split-weight duplicate nodes aggregate to the original solution", {
  fixture <- moment_metamorphic_fixture(seed = 1403L, nx = 3L, ny = 4L)
  source_index <- c(1L, 2L, 2L, 3L)
  target_index <- c(1L, 2L, 3L, 3L, 4L)
  source_fraction <- c(1, 0.35, 0.65, 1)
  target_fraction <- c(1, 1, 0.4, 0.6, 1)
  duplicated <- list(
    x = fixture$x[source_index, , drop = FALSE],
    y = fixture$y[target_index, , drop = FALSE],
    source_features = fixture$source_features[source_index, , drop = FALSE],
    target_features = fixture$target_features[target_index, , drop = FALSE],
    wx = fixture$wx[source_index] * source_fraction,
    wy = fixture$wy[target_index] * target_fraction
  )
  duplicated$Cx <- sqeuclidean_cost(duplicated$x)
  duplicated$Cy <- sqeuclidean_cost(duplicated$y)
  duplicated$M <- sqeuclidean_cost(
    duplicated$source_features,
    duplicated$target_features
  )

  original <- moment_metamorphic_fit(fixture)
  split <- moment_metamorphic_fit(duplicated, block_size = 3L)
  aggregate_plan <- function(plan) {
    by_source <- rowsum(plan, source_index, reorder = FALSE)
    unname(t(rowsum(t(by_source), target_index, reorder = FALSE)))
  }

  expect_true(original$converged)
  expect_true(split$converged)
  expect_equal(sum(duplicated$wx), sum(fixture$wx), tolerance = 1e-15)
  expect_equal(sum(duplicated$wy), sum(fixture$wy), tolerance = 1e-15)
  expect_equal(split$fugw_cost, original$fugw_cost, tolerance = 5e-8)
  for (role in c("sample", "feature")) {
    expect_equal(
      aggregate_plan(transport_plan_materialize(split$plans[[role]])),
      transport_plan_materialize(original$plans[[role]]),
      tolerance = 8e-8
    )
  }
})

test_that("serialized fits preserve apply and adjoint for both couplings", {
  fixture <- moment_metamorphic_fixture(seed = 1404L)
  fit <- moment_metamorphic_fit(fixture)
  restored <- unserialize(serialize(fit, NULL))
  source_values <- matrix(rnorm(nrow(fixture$x) * 3L), nrow(fixture$x), 3L)
  target_values <- matrix(rnorm(nrow(fixture$y) * 2L), nrow(fixture$y), 2L)

  expect_identical(restored$status, fit$status)
  expect_equal(restored$fugw_cost, fit$fugw_cost, tolerance = 0)
  for (role in c("sample", "feature")) {
    expect_equal(
      transport_plan_apply(restored$plans[[role]], target_values),
      transport_plan_apply(fit$plans[[role]], target_values),
      tolerance = 2e-13
    )
    expect_equal(
      transport_plan_adjoint(restored$plans[[role]], source_values),
      transport_plan_adjoint(fit$plans[[role]], source_values),
      tolerance = 2e-13
    )
  }
})

test_that("thread settings are deterministic and numerical extremes fail closed", {
  fixture <- moment_metamorphic_fixture(seed = 1405L)
  thread_variables <- c(
    "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS", "BLIS_NUM_THREADS"
  )
  previous <- Sys.getenv(thread_variables, unset = NA_character_, names = TRUE)
  on.exit({
    for (name in thread_variables) {
      if (is.na(previous[[name]])) {
        Sys.unsetenv(name)
      } else {
        do.call(Sys.setenv, stats::setNames(list(previous[[name]]), name))
      }
    }
  }, add = TRUE)

  Sys.setenv(
    OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
    MKL_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1",
    BLIS_NUM_THREADS = "1"
  )
  one_thread <- moment_metamorphic_fit(fixture)
  Sys.setenv(OMP_NUM_THREADS = "4")
  four_requested <- moment_metamorphic_fit(fixture)
  expect_moment_fits_equivalent(one_thread, four_requested, tolerance = 2e-12)

  expect_error(
    correlation_cost(matrix(1, 4L, 3L)),
    "zero-norm feature profile"
  )
  expect_error(
    sqeuclidean_cost(matrix(c(1e200, 0, -1e200, 0), 2L, 2L, byrow = TRUE)),
    "finite numeric vector"
  )

  extreme <- tryCatch(
    moment_metamorphic_fit(
      fixture,
      controls = list(
        epsilon = 1e-14,
        reg_marginals = c(1e-12, 1e12),
        max_iter = 2L,
        max_iter_ot = 2L,
        tol = 1e-14,
        tol_ot = 1e-14
      )
    ),
    error = identity
  )
  if (inherits(extreme, "error")) {
    expect_match(conditionMessage(extreme), "finite|numerical|represent|inner")
  } else {
    expect_false(extreme$converged)
    expect_false(extreme$certificate$outer_stationarity$certified)
    expect_false(extreme$final_block_optimality_certified)
    expect_identical(extreme$global_optimality, "not_claimed")
    expect_false(extreme$status %in% c(
      "converged_stationary", "converged_approximate_stationary"
    ))
  }
})
