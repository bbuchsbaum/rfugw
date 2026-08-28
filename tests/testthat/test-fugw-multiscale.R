source(bench_test_resource("benchmark_moment_fugw_multiscale.R"))

multiscale_domain_fixture <- function() {
  source_coarse_coordinates <- matrix(
    c(0, 0, 1, 0, 2, 0), 3, 2, byrow = TRUE
  )
  target_coarse_coordinates <- matrix(
    c(0, 0, 1.2, 0), 2, 2, byrow = TRUE
  )
  source_fine_coordinates <- rbind(
    source_coarse_coordinates[1, ], source_coarse_coordinates[1, ] + c(.1, .05),
    source_coarse_coordinates[2, ], source_coarse_coordinates[2, ] + c(.1, -.05),
    source_coarse_coordinates[3, ], source_coarse_coordinates[3, ] + c(.1, .02)
  )
  target_fine_coordinates <- rbind(
    target_coarse_coordinates[1, ], target_coarse_coordinates[1, ] + c(.1, .03),
    target_coarse_coordinates[2, ], target_coarse_coordinates[2, ] + c(.1, -.03)
  )
  source_coarse_weights <- c(.2, .3, .5)
  target_coarse_weights <- c(.4, .6)
  source_fine_weights <- c(.05, .15, .1, .2, .2, .3)
  target_fine_weights <- c(.1, .3, .2, .4)
  source_coarse_features <- cbind(
    source_coarse_coordinates[, 1], c(.1, .4, .8), c(.7, .2, .5)
  )
  target_coarse_features <- cbind(
    target_coarse_coordinates[, 1], c(.15, .75), c(.65, .45)
  )
  source_fine_features <- cbind(
    source_fine_coordinates[, 1],
    c(.08, .12, .35, .45, .76, .84),
    c(.72, .68, .24, .16, .48, .52)
  )
  target_fine_features <- cbind(
    target_fine_coordinates[, 1],
    c(.12, .18, .7, .8),
    c(.68, .62, .48, .42)
  )
  list(
    source = list(
      fugw_domain(
        sqeuclidean_cost(source_coarse_coordinates),
        source_coarse_features,
        source_coarse_weights,
        name = "source-coarse"
      ),
      fugw_domain(
        sqeuclidean_cost(source_fine_coordinates),
        source_fine_features,
        source_fine_weights,
        parent = rep(seq_len(3), each = 2),
        name = "source-fine"
      )
    ),
    target = list(
      fugw_domain(
        sqeuclidean_cost(target_coarse_coordinates),
        target_coarse_features,
        target_coarse_weights,
        name = "target-coarse"
      ),
      fugw_domain(
        sqeuclidean_cost(target_fine_coordinates),
        target_fine_features,
        target_fine_weights,
        parent = rep(seq_len(2), each = 2),
        name = "target-fine"
      )
    )
  )
}

test_that("fugw domains validate factorized geometry and hierarchy metadata", {
  coordinates <- cbind(seq_len(4), c(0, 1, 0, 1))
  domain <- fugw_domain(
    sqeuclidean_cost(coordinates),
    features = cbind(coordinates, seq_len(4)),
    weights = c(1, 2, 3, 4),
    parent = c(1, 1, 2, 2),
    name = "fine"
  )

  expect_s3_class(domain, "rfugw_fugw_domain")
  expect_equal(sum(domain$weights), 1)
  expect_identical(domain$parent, c(1L, 1L, 2L, 2L))
  expect_match(paste(capture.output(print(domain)), collapse = "\n"), "points:.*4")
  expect_error(
    fugw_domain(sqeuclidean_cost(coordinates, coordinates[1:3, ])),
    "square"
  )
  expect_error(
    fugw_domain(block_cost(4, 4, function(i, j) outer(i, j, "-"))),
    "no affine-bilinear factors"
  )
  expect_error(
    fugw_domain(sqeuclidean_cost(coordinates), matrix(1, 3, 2)),
    "one row"
  )
  expect_error(
    fugw_domain(sqeuclidean_cost(coordinates), parent = c(1, 1, 0, 2)),
    "positive integer"
  )
})

test_that("conditional hierarchy lift has an independent dense oracle", {
  fixture <- multiscale_domain_fixture()
  source_coarse <- fixture$source[[1]]
  target_coarse <- fixture$target[[1]]
  source_fine <- fixture$source[[2]]
  target_fine <- fixture$target[[2]]
  coarse <- fugw_factorized(
    source_coarse$structure,
    target_coarse$structure,
    wx = source_coarse$weights,
    wy = target_coarse$weights,
    M = sqeuclidean_cost(
      source_coarse$features, target_coarse$features
    ),
    epsilon = .2,
    max_iter = 25L,
    tol = 1e-8,
    max_iter_ot = 3000L,
    tol_ot = 1e-10,
    feature_weight = .3,
    structure_weight = 1,
    block_size = 2L
  )
  lifted <- rfugw:::.lift_plan_summary(
    coarse$implicit_plans$sample,
    coarse$potentials$sample,
    source_coarse,
    target_coarse,
    source_fine,
    target_fine
  )

  coarse_plan <- transport_plan_materialize(coarse$implicit_plans$sample)
  source_conditional <- source_fine$weights /
    source_coarse$weights[source_fine$parent]
  target_conditional <- target_fine$weights /
    target_coarse$weights[target_fine$parent]
  explicit_lift <- coarse_plan[
    source_fine$parent, target_fine$parent, drop = FALSE
  ] * outer(source_conditional, target_conditional)
  source_right <- rfugw:::.cost_factor_pair(source_fine$structure)$right
  target_right <- rfugw:::.cost_factor_pair(target_fine$structure)$right
  positive <- explicit_lift > 0
  entropy_sum <- sum(explicit_lift[positive] * log(explicit_lift[positive]))
  product_reference <- outer(source_fine$weights, target_fine$weights)
  product_kl <- sum(
    explicit_lift[positive] *
      log(explicit_lift[positive] / product_reference[positive])
  ) - sum(explicit_lift) + 1

  expect_equal(lifted$stats$source_marginal, rowSums(explicit_lift),
               tolerance = 2e-12)
  expect_equal(lifted$stats$target_marginal, colSums(explicit_lift),
               tolerance = 2e-12)
  expect_equal(lifted$stats$mass, sum(explicit_lift), tolerance = 2e-12)
  expect_equal(lifted$stats$entropy_sum, entropy_sum, tolerance = 2e-12)
  expect_equal(lifted$stats$plan_product_kl, product_kl, tolerance = 2e-12)
  expect_equal(
    lifted$moments$H,
    crossprod(source_right, explicit_lift %*% target_right),
    tolerance = 2e-12
  )
  expect_equal(
    lifted$moments$G_source,
    crossprod(source_right, source_right * rowSums(explicit_lift)),
    tolerance = 2e-12
  )
  expect_equal(
    lifted$moments$G_target,
    crossprod(target_right, target_right * colSums(explicit_lift)),
    tolerance = 2e-12
  )
  raw_source <- coarse$potentials$sample$source_bar[source_fine$parent]
  raw_target <- coarse$potentials$sample$target_bar[target_fine$parent]
  expected_shift <- 0.5 * (
    sum(target_fine$weights * raw_target) -
      sum(source_fine$weights * raw_source)
  )
  expect_equal(lifted$potentials$source_bar, raw_source + expected_shift)
  expect_equal(lifted$potentials$target_bar, raw_target - expected_shift)
  expect_true(lifted$gauge$certified)
  expect_lte(lifted$gauge$centered_mean_residual, 1e-14)
  expect_lte(lifted$gauge$log_kernel_invariance_residual, 1e-14)

  shifted_potentials <- coarse$potentials$sample
  shifted_potentials$source_bar <- shifted_potentials$source_bar + 7
  shifted_potentials$target_bar <- shifted_potentials$target_bar - 7
  shifted <- rfugw:::.lift_plan_summary(
    coarse$implicit_plans$sample,
    shifted_potentials,
    source_coarse,
    target_coarse,
    source_fine,
    target_fine
  )
  expect_equal(shifted$potentials, lifted$potentials, tolerance = 2e-14)
  expect_equal(shifted$stats, lifted$stats, tolerance = 2e-14)
  expect_equal(shifted$moments, lifted$moments, tolerance = 2e-14)
})

test_that("one-level multiscale is exactly the factorized solver", {
  fixture <- multiscale_domain_fixture()
  source <- fixture$source[[1]]
  target <- fixture$target[[1]]
  feature <- sqeuclidean_cost(source$features, target$features)
  direct <- fugw_factorized(
    source$structure,
    target$structure,
    wx = source$weights,
    wy = target$weights,
    M = feature,
    reg_marginals = c(4, 6),
    epsilon = .2,
    max_iter = 30L,
    tol = 1e-8,
    max_iter_ot = 3000L,
    tol_ot = 1e-10,
    feature_weight = .3,
    structure_weight = 1,
    block_size = 2L
  )
  multiscale <- fugw_multiscale(
    source,
    target,
    feature_metric = "sqeuclidean",
    feature_weight = .3,
    structure_weight = 1,
    reg_marginals = c(4, 6),
    epsilon = .2,
    max_iter = 30L,
    tol = 1e-8,
    max_iter_ot = 3000L,
    tol_ot = 1e-10,
    block_size = 2L
  )

  expect_equal(multiscale$fugw_cost, direct$fugw_cost, tolerance = 2e-12)
  expect_equal(
    transport_plan_materialize(multiscale$plans$sample),
    transport_plan_materialize(direct$plans$sample),
    tolerance = 2e-12
  )
  expect_identical(multiscale$hierarchy$depth, 1L)
  expect_identical(multiscale$formulation, "fugw_multiscale_joint_kl")
  expect_identical(multiscale$level_trace$initialization, "independent_product")
})

test_that("multiscale transfers potentials, masses, and moments without dense plans", {
  fixture <- multiscale_domain_fixture()
  fit <- fugw_multiscale(
    fixture$source,
    fixture$target,
    feature_metric = "sqeuclidean",
    feature_weight = .1,
    structure_weight = 1,
    reg_marginals = c(5, 5),
    epsilon = c(.3, .2),
    max_iter = c(30L, 40L),
    tol = 1e-7,
    max_iter_ot = 5000L,
    tol_ot = 1e-9,
    block_size = 2L
  )

  expect_length(fit$level_results, 2L)
  expect_equal(fit$level_trace$epsilon, c(.3, .2))
  expect_identical(
    fit$level_trace$initialization,
    c("independent_product", "transferred_summary")
  )
  expect_true(fit$level_results[[2]]$initialization$transferred)
  expect_true(fit$transfers[[1]]$potentials_prolonged)
  expect_true(fit$transfers[[1]]$potentials_gauge_recentered)
  expect_true(fit$transfers[[1]]$sample_gauge$certified)
  expect_true(fit$transfers[[1]]$feature_gauge$certified)
  expect_true(fit$transfers[[1]]$moments_recomputed_in_fine_basis)
  expect_true(fit$transfers[[1]]$mass_preserved)
  expect_lte(fit$transfers[[1]]$source_weight_aggregation_residual, 1e-14)
  expect_lte(fit$transfers[[1]]$target_weight_aggregation_residual, 1e-14)
  expect_identical(
    transport_plan_representation(fit$plans$sample)$representation,
    "implicit_operator"
  )
  expect_false(fit$runtime_provenance$multiscale$intermediate_dense_plans)
  expect_identical(fit$certificate$multiscale$support$mode, "full_implicit")
  expect_equal(
    fit$certificate$multiscale$support$omitted_kernel_mass_bound, 0
  )
  expect_identical(fit$certificate$multiscale$global_optimality, "not_claimed")
  expect_equal(fit$transfers[[1]]$epsilon_schedule$from, .3)
  expect_equal(fit$transfers[[1]]$epsilon_schedule$to, .2)
  expect_true(fit$transfers[[1]]$epsilon_schedule$changed)
  expect_identical(
    fit$transfers[[1]]$epsilon_schedule$potential_units, "cost"
  )
  expect_identical(
    fit$transfers[[1]]$epsilon_schedule$potential_rescaling, "none"
  )
  expect_equal(fit$hierarchy$controls$epsilon, c(.3, .2))
  expect_equal(fit$hierarchy$controls$feature_weight, c(.1, .1))
  expect_equal(fit$hierarchy$controls$structure_weight, c(1, 1))
  expect_identical(
    fit$hierarchy$controls$level_contract,
    rep("certified_endpoint", 2L)
  )
  expect_true(all(vapply(
    fit$level_results,
    function(x) inherits(x$implicit_plans$sample, "rfugw_transport_plan"),
    logical(1)
  )))
})

test_that("multiscale hierarchy and unsupported adaptive support fail closed", {
  fixture <- multiscale_domain_fixture()
  missing_parent <- fixture$source
  missing_parent[[2]]$parent <- NULL
  expect_error(
    fugw_multiscale(missing_parent, fixture$target, max_iter = 1),
    "must define a parent map"
  )

  bad_weights <- fixture$source
  bad_weights[[2]]$weights <- rep(1 / bad_weights[[2]]$n, bad_weights[[2]]$n)
  expect_error(
    fugw_multiscale(bad_weights, fixture$target, max_iter = 1),
    "do not aggregate"
  )

  expect_error(
    fugw_multiscale(
      fixture$source, fixture$target,
      support = "adaptive", max_iter = 1
    ),
    "not implemented"
  )
  expect_error(
    fugw_multiscale(fixture$source, fixture$target[[1]], max_iter = 1),
    "equal depth"
  )
  mismatched <- fixture$target
  mismatched[[2]]$features <- mismatched[[2]]$features[, 1:2, drop = FALSE]
  expect_error(
    fugw_multiscale(fixture$source, mismatched, max_iter = 1),
    "feature dimensions differ"
  )
})

test_that("multiscale status distinguishes final-level-only certification", {
  fixture <- multiscale_domain_fixture()
  common <- list(
    source = fixture$source,
    target = fixture$target,
    feature_metric = "sqeuclidean",
    feature_weight = .1,
    structure_weight = 1,
    reg_marginals = c(5, 5),
    epsilon = c(.3, .2),
    tol = 1e-7,
    max_iter_ot = 5000L,
    tol_ot = 1e-9,
    block_size = 2L
  )
  final_only <- do.call(
    fugw_multiscale,
    c(common, list(max_iter = c(1L, 40L)))
  )
  final_starved <- do.call(
    fugw_multiscale,
    c(common, list(max_iter = c(30L, 1L)))
  )

  expect_true(final_only$converged)
  expect_false(final_only$multiscale_certified)
  expect_identical(final_only$status, "converged_final_level_only")
  expect_identical(
    final_only$certificate$multiscale$classification,
    "final_level_only"
  )
  expect_true(final_only$certificate$multiscale$final_level_stationary)
  expect_false(final_only$certificate$multiscale$all_levels_stationary)
  expect_identical(
    final_only$certificate$outer_stationarity$scope,
    "final_level_only"
  )
  expect_identical(final_only$warning_payload$code, "final_level_only")
  expect_identical(final_only$global_optimality, "not_claimed")

  expect_false(final_starved$converged)
  expect_false(final_starved$multiscale_certified)
  expect_identical(
    final_starved$certificate$multiscale$classification,
    "uncertified"
  )
  expect_false(final_starved$certificate$multiscale$final_level_stationary)
})

test_that("identity-parent repeated levels preserve the certified endpoint", {
  set.seed(417L)
  coordinates <- matrix(rnorm(12), 6L, 2L)
  features <- matrix(rnorm(18), 6L, 3L)
  weights <- runif(6L)
  weights <- weights / sum(weights)
  coarse <- fugw_domain(
    sqeuclidean_cost(coordinates), features, weights, name = "identity-coarse"
  )
  fine <- fugw_domain(
    sqeuclidean_cost(coordinates), features, weights,
    parent = seq_len(6L), name = "identity-fine"
  )
  controls <- list(
    reg_marginals = c(4, 7),
    epsilon = .25,
    max_iter = 50L,
    tol = 1e-8,
    max_iter_ot = 5000L,
    tol_ot = 1e-10,
    feature_weight = .4,
    structure_weight = 1,
    block_size = 3L
  )
  direct <- do.call(fugw_factorized, c(
    list(
      Cx = coarse$structure,
      Cy = coarse$structure,
      M = sqeuclidean_cost(features, features),
      wx = weights,
      wy = weights
    ),
    controls
  ))
  repeated <- do.call(fugw_multiscale, c(
    list(
      source = list(coarse, fine),
      target = list(coarse, fine),
      feature_metric = "sqeuclidean"
    ),
    controls
  ))
  target_values <- matrix(rnorm(24), 6L, 4L)
  source_values <- matrix(rnorm(18), 6L, 3L)

  expect_true(direct$converged)
  expect_true(repeated$multiscale_certified)
  expect_identical(repeated$status, "converged_stationary")
  expect_equal(repeated$fugw_cost, direct$fugw_cost, tolerance = 1e-7)
  expect_equal(
    transport_plan_apply(repeated$plans$sample, target_values),
    transport_plan_apply(direct$plans$sample, target_values),
    tolerance = 1e-7
  )
  expect_equal(
    transport_plan_adjoint(repeated$plans$sample, source_values),
    transport_plan_adjoint(direct$plans$sample, source_values),
    tolerance = 1e-7
  )
  expect_true(all(repeated$level_trace$contract_satisfied))
})

test_that("budgeted warm starts have a distinct honest final-level status", {
  fixture <- multiscale_domain_fixture()
  fit <- fugw_multiscale(
    fixture$source,
    fixture$target,
    feature_metric = "sqeuclidean",
    feature_weight = .1,
    structure_weight = 1,
    reg_marginals = c(5, 5),
    epsilon = c(.3, .2),
    max_iter = c(1L, 40L),
    tol = 1e-7,
    max_iter_ot = 5000L,
    tol_ot = 1e-9,
    block_size = 2L,
    level_contract = c("budgeted_warm_start", "certified_endpoint")
  )

  expect_true(fit$converged)
  expect_false(fit$multiscale_certified)
  expect_true(fit$certificate$multiscale$execution_contract_certified)
  expect_identical(
    fit$status, "converged_final_level_with_budgeted_warm_starts"
  )
  expect_identical(
    fit$certificate$multiscale$classification,
    "final_level_with_budgeted_warm_starts"
  )
  expect_identical(
    fit$certificate$outer_stationarity$scope,
    "final_level_with_budgeted_warm_starts"
  )
  expect_identical(
    fit$level_trace$level_contract,
    c("budgeted_warm_start", "certified_endpoint")
  )
  expect_true(all(fit$level_trace$contract_satisfied))
  expect_false(fit$level_trace$outer_stationarity_certified[[1L]])
  expect_true(fit$level_trace$outer_stationarity_certified[[2L]])
  expect_identical(fit$warning_payload$code, "budgeted_warm_starts")
  expect_error(
    fugw_multiscale(
      fixture$source,
      fixture$target,
      max_iter = 1L,
      level_contract = "budgeted_warm_start"
    ),
    "finest level"
  )
})

test_that("planted transfer efficacy meets frozen work and accuracy gates", {
  evidence <- utils::read.csv(
    bench_test_resource("moment-fugw-multiscale-evidence.csv"),
    stringsAsFactors = FALSE
  )

  expect_equal(evidence$seed, 1701:1710)
  expect_equal(nrow(evidence), 10L)
  expect_true(all(evidence$same_final_contract))
  expect_true(all(evidence$cold_certified))
  expect_true(all(evidence$transfer_final_certified))
  expect_true(all(evidence$transfer_multiscale_certified))
  expect_gte(mean(evidence$inner_iteration_reduction > 0), .8)
  expect_gte(stats::median(evidence$inner_iteration_reduction), .2)
  expect_lte(max(evidence$objective_relative_difference), 1e-4)
  expect_gte(min(evidence$heldout_score_ratio), .99)
  expect_true(all(evidence$efficacy_gate))

  live <- bench_run_moment_fugw_multiscale_case(1701L)
  expect_true(live$same_final_contract)
  expect_true(live$cold_certified)
  expect_true(live$transfer_multiscale_certified)
  expect_gt(live$inner_iteration_reduction, .2)
  expect_lte(live$objective_relative_difference, 1e-4)
  expect_gte(live$heldout_score_ratio, .99)
})
