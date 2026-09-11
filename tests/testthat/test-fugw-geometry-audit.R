test_that("exact geometry constructors verify metric invariants", {
  coordinates <- matrix(
    c(0, 0, 1, 0, 0, 1, 1, 1, .3, .7),
    ncol = 2,
    byrow = TRUE
  )
  audit <- geometry_audit(sqeuclidean_cost(coordinates))

  expect_true(audit$certified)
  expect_true(audit$exact)
  expect_identical(audit$status, "exact_verified")
  expect_identical(audit$pair_design$mode, "exact_all_unordered_pairs")
  expect_equal(audit$symmetry_residual, 0, tolerance = 1e-14)
  expect_equal(audit$nonnegativity_violation, 0, tolerance = 1e-14)
  expect_equal(audit$diagonal_residual, 0, tolerance = 1e-14)
  expect_true(audit$invariants_certified)
  expect_true(audit$constructor_verified_exact)
})

test_that("computed audits reject metadata-only exactness and bad self-costs", {
  coordinates <- cbind(seq_len(5), c(0, 1, 0, 1, .5))
  exact_cost <- sqeuclidean_cost(coordinates)
  terms <- rfugw:::.cost_native_terms(exact_cost)
  metadata_only <- factorized_cost(
    terms$left,
    terms$right,
    row = terms$row,
    column = terms$column,
    exact = TRUE,
    provenance = list(metric = "claimed_exact_without_constructor_receipt")
  )
  invalid_diagonal <- factorized_cost(
    terms$left,
    terms$right,
    row = terms$row + .125,
    column = terms$column + .125,
    exact = TRUE,
    provenance = list(metric = "falsified_self_distance")
  )
  metadata_audit <- geometry_audit(metadata_only)
  invalid_audit <- geometry_audit(invalid_diagonal)

  expect_true(metadata_audit$invariants_certified)
  expect_false(metadata_audit$certified)
  expect_identical(
    metadata_audit$status,
    "unverified_or_falsified_exact_claim"
  )
  expect_false(metadata_audit$constructor_verified_exact)
  expect_false(invalid_audit$certified)
  expect_identical(invalid_audit$status, "invalid_geometry_invariants")
  expect_gte(invalid_audit$diagonal_residual, .25 - 1e-12)
})

test_that("embedded geometry retains reproducible train and held-out evidence", {
  coordinates <- cbind(seq_len(12) / 12, (seq_len(12) %% 4) / 4)
  embedding <- 1.01 * coordinates
  reference <- sqeuclidean_cost(coordinates)
  cost <- embedded_geometry_cost(
    embedding,
    reference,
    method = "scaled_test_embedding",
    reference_metric = "sqeuclidean_original",
    seed = 7727L,
    train_pairs = 10L,
    holdout_pairs = 9L,
    relative_error_tolerance = .03
  )
  set.seed(991)
  rng_before <- .Random.seed
  first <- geometry_audit(cost, exact_pair_limit = 1L)
  rng_after <- .Random.seed
  second <- geometry_audit(cost, exact_pair_limit = 1L)

  expect_identical(rng_after, rng_before)
  expect_true(first$certified)
  expect_true(first$approximate)
  expect_false(first$exact)
  expect_identical(first$status, "approximate_verified")
  expect_equal(first$weighted_stress, .0201, tolerance = 1e-12)
  expect_equal(first$heldout_error, .0201, tolerance = 1e-12)
  expect_identical(first$embedding_method, "scaled_test_embedding")
  expect_identical(first$embedding_rank, 2L)
  expect_identical(first$pair_design$seed, 7727L)
  expect_identical(
    first$pair_design$mode,
    "reproducible_disjoint_train_holdout"
  )
  expect_length(
    intersect(
      apply(first$pair_design$train_pairs, 1, paste, collapse = ":"),
      apply(first$pair_design$holdout_pairs, 1, paste, collapse = ":")
    ),
    0L
  )
  expect_identical(first$pair_design, second$pair_design)
  expect_identical(first$reference_metric, "sqeuclidean_original")
  expect_identical(
    first$reference_provenance$metric,
    "sqeuclidean"
  )
})

test_that("held-out pairs catch a train-only geometry fit", {
  n <- 20L
  seed <- 181L
  train_count <- 6L
  holdout_count <- 6L
  coordinates <- cbind(seq_len(n) / n, (seq_len(n) / n)^2)
  represented <- cost_block(sqeuclidean_cost(coordinates))
  design <- rfugw:::.geometry_sample_pairs(
    n, train_count + holdout_count, seed
  )
  heldout <- design[seq.int(train_count + 1L, nrow(design)), , drop = FALSE]
  reference <- represented
  for (index in seq_len(nrow(heldout))) {
    i <- heldout[index, 1L]
    j <- heldout[index, 2L]
    reference[i, j] <- 2 * reference[i, j]
    reference[j, i] <- reference[i, j]
  }
  cost <- embedded_geometry_cost(
    coordinates,
    reference,
    method = "train_only_mutant",
    seed = seed,
    train_pairs = train_count,
    holdout_pairs = holdout_count,
    relative_error_tolerance = .1
  )
  audit <- geometry_audit(cost, exact_pair_limit = 1L)

  expect_equal(audit$weighted_stress, 0, tolerance = 1e-14)
  expect_gte(audit$heldout_error, .49)
  expect_false(audit$certified)
  expect_identical(audit$status, "heldout_approximation_failure")
  expect_identical(audit$pair_design$train_pairs, design[seq_len(6L), ])
  expect_identical(audit$pair_design$holdout_pairs, heldout)
})

test_that("verified approximate geometry changes status without solver error", {
  source_coordinates <- matrix(
    c(0, 0, 1, 0, 0, 1, 1, 1), 4L, 2L, byrow = TRUE
  )
  target_coordinates <- matrix(
    c(0, 0, 1.1, 0, 0, .9, 1, 1), 4L, 2L, byrow = TRUE
  )
  source <- embedded_geometry_cost(
    1.01 * source_coordinates,
    sqeuclidean_cost(source_coordinates),
    method = "small_verified_embedding",
    relative_error_tolerance = .03
  )
  fit <- fugw_factorized(
    source,
    sqeuclidean_cost(target_coordinates),
    epsilon = .3,
    reg_marginals = c(4, 6),
    max_iter = 40L,
    tol = 1e-8,
    max_iter_ot = 5000L,
    tol_ot = 1e-10,
    block_size = 2L
  )

  expect_true(fit$converged)
  expect_true(fit$geometry_certified)
  expect_false(fit$geometry_exact)
  expect_identical(fit$status, "converged_approximate_stationary")
  expect_identical(fit$certificate$geometry$status, "approximate_verified")
  expect_true(fit$certificate$geometry$source$relative_error_verified)
  expect_equal(fit$geometry_heldout_error, .0201, tolerance = 1e-12)
  expect_true(fit$certificate$outer_stationarity$certified)
})

geometry_hierarchy_fixture <- function() {
  coarse_coordinates <- cbind(c(0, 1), c(0, 0))
  fine_coordinates <- coarse_coordinates[c(1, 1, 2, 2), , drop = FALSE]
  coarse_weights <- c(.4, .6)
  fine_weights <- c(.1, .3, .2, .4)
  coarse_features <- matrix(c(1, 2, 3, 4), 2L, 2L, byrow = TRUE)
  fine_features <- coarse_features[c(1, 1, 2, 2), , drop = FALSE]
  list(
    coarse = fugw_domain(
      sqeuclidean_cost(coarse_coordinates),
      coarse_features,
      coarse_weights,
      name = "coarse"
    ),
    fine = fugw_domain(
      sqeuclidean_cost(fine_coordinates),
      fine_features,
      fine_weights,
      parent = c(1, 1, 2, 2),
      name = "fine",
      feature_aggregation_policy = "weighted_mean",
      feature_aggregation_tolerance = 1e-12,
      geometry_aggregation_tolerance = 1e-12
    )
  )
}

test_that("hierarchy audits cover aggregation, geometry, and round trips", {
  fixture <- geometry_hierarchy_fixture()
  fit <- fugw_multiscale(
    list(fixture$coarse, fixture$fine),
    list(fixture$coarse, fixture$fine),
    feature_metric = "sqeuclidean",
    epsilon = c(.3, .2),
    max_iter = c(25L, 35L),
    tol = 1e-8,
    max_iter_ot = 5000L,
    tol_ot = 1e-10,
    block_size = 2L
  )
  transition <- fit$certificate$multiscale$hierarchy$source[[1L]]

  expect_true(fit$multiscale_certified)
  expect_true(fit$hierarchy_transfer_certified)
  expect_true(transition$parent_coverage$all_coarse_points_covered)
  expect_equal(transition$weight_residual, 0, tolerance = 1e-14)
  expect_true(transition$feature_aggregation$certified)
  expect_identical(
    transition$feature_aggregation$policy,
    "weighted_mean"
  )
  expect_equal(
    transition$feature_aggregation$relative_error, 0, tolerance = 1e-14
  )
  expect_true(transition$geometry_agreement$certified)
  expect_equal(
    transition$geometry_agreement$weighted_stress, 0, tolerance = 1e-14
  )
  expect_true(transition$restrict_prolong$certified)
  expect_equal(
    transition$restrict_prolong$value_round_trip_residual,
    0,
    tolerance = 1e-14
  )
  expect_equal(fit$hierarchy_transfer_error, 0, tolerance = 1e-14)
  expect_true(is.finite(fit$residual))
  expect_true(is.finite(fit$geometry_relative_error))
  expect_identical(fit$support_tail_bound, 0)
})

test_that("inconsistent coarse geometry and feature policy fail closed", {
  fixture <- geometry_hierarchy_fixture()
  bad_geometry <- fixture$coarse
  bad_geometry$structure <- scaled_cost(bad_geometry$structure, 100)
  bad_features <- fixture$fine
  bad_features$features[[1L]] <- bad_features$features[[1L]] + 1

  expect_error(
    fugw_multiscale(
      list(bad_geometry, fixture$fine),
      list(fixture$coarse, fixture$fine),
      max_iter = 1L
    ),
    "coarse geometry disagrees"
  )
  expect_error(
    fugw_multiscale(
      list(fixture$coarse, bad_features),
      list(fixture$coarse, fixture$fine),
      max_iter = 1L
    ),
    "weighted-mean aggregation policy"
  )
})
