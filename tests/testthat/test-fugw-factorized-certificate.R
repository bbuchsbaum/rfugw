moment_certificate_fixture <- function(seed = 731L) {
  set.seed(seed)
  source_coordinates <- matrix(rnorm(8), 4L, 2L)
  target_coordinates <- matrix(rnorm(10), 5L, 2L)
  source_features <- matrix(rnorm(12), 4L, 3L)
  target_features <- matrix(rnorm(15), 5L, 3L)
  wx <- runif(4L)
  wy <- runif(5L)
  list(
    Cx = sqeuclidean_cost(source_coordinates),
    Cy = sqeuclidean_cost(target_coordinates),
    M = sqeuclidean_cost(source_features, target_features),
    wx = wx / sum(wx),
    wy = wy / sum(wy)
  )
}

moment_certificate_controls <- function(max_iter_ot = 10000L) {
  list(
    reg_marginals = c(3, 6),
    epsilon = .3,
    feature_weight = .4,
    structure_weight = 1,
    max_iter = 80L,
    tol = 1e-8,
    max_iter_ot = max_iter_ot,
    tol_ot = 1e-10,
    rescale_plan = TRUE,
    block_size = 2L
  )
}

moment_certificate_fit <- function(
    fixture = moment_certificate_fixture(), max_iter_ot = 10000L) {
  do.call(fugw_factorized, c(
    list(
      Cx = fixture$Cx,
      Cy = fixture$Cy,
      M = fixture$M,
      wx = fixture$wx,
      wy = fixture$wy
    ),
    moment_certificate_controls(max_iter_ot)
  ))
}

moment_certificate_final_costs <- function(fit, fixture) {
  controls <- moment_certificate_controls()
  source_left <- rfugw:::.cost_factor_pair(fixture$Cx)$left
  target_left <- rfugw:::.cost_factor_pair(fixture$Cy)$left
  feature_terms <- rfugw:::.cost_native_terms(fixture$M)
  sample_state <- rfugw:::.factorized_plan_state(fit$implicit_plans$sample)
  feature_state <- rfugw:::.factorized_plan_state(fit$implicit_plans$feature)
  build <- function(moments, stats, role) {
    rfugw:::.fugw_dynamic_cost(
      moments,
      stats,
      source_left,
      target_left,
      feature_terms,
      fixture$wx,
      fixture$wy,
      controls$reg_marginals,
      controls$epsilon,
      controls$feature_weight,
      controls$structure_weight,
      list(exact = TRUE, role = role)
    )
  }
  list(
    sample = build(fit$moments$feature, feature_state$stats, "sample_test"),
    feature = build(fit$moments$sample, sample_state$stats, "feature_test")
  )
}

moment_certificate_gkl <- function(value, reference) {
  if (any(value < 0) || any(reference < 0) ||
      any(value > 0 & reference == 0)) return(Inf)
  active <- value > 0
  sum(value[active] * log(value[active] / reference[active])) -
    sum(value) + sum(reference)
}

test_that("final post-rescale KKT audits match an explicit derivative oracle", {
  fixture <- moment_certificate_fixture()
  controls <- moment_certificate_controls()
  fit <- moment_certificate_fit(fixture)
  final_costs <- moment_certificate_final_costs(fit, fixture)
  audit <- fit$certificate$outer_stationarity$final_blocks$sample
  plan <- transport_plan_materialize(fit$implicit_plans$sample)
  opponent_mass <- transport_plan_mass(fit$implicit_plans$feature)
  epsilon_block <- controls$epsilon * opponent_mass
  rho_block <- controls$reg_marginals * opponent_mass
  source_active <- fixture$wx > 0
  target_active <- fixture$wy > 0
  active_plan <- plan[source_active, target_active, drop = FALSE]
  reference <- outer(
    fixture$wx[source_active], fixture$wy[target_active]
  )
  derivative <- cost_block(final_costs$sample)[
    source_active, target_active, drop = FALSE
  ] +
    outer(
      rho_block[[1L]] * log(
        rowSums(plan)[source_active] / fixture$wx[source_active]
      ),
      rep(1, sum(target_active))
    ) +
    outer(
      rep(1, sum(source_active)),
      rho_block[[2L]] * log(
        colSums(plan)[target_active] / fixture$wy[target_active]
      )
    ) +
    epsilon_block * log(active_plan / reference)
  direct_objective <- sum(cost_block(final_costs$sample) * plan) +
    rho_block[[1L]] * moment_certificate_gkl(
      rowSums(plan), fixture$wx
    ) +
    rho_block[[2L]] * moment_certificate_gkl(
      colSums(plan), fixture$wy
    ) +
    epsilon_block * moment_certificate_gkl(
      as.numeric(plan), as.numeric(outer(fixture$wx, fixture$wy))
    )

  expect_true(fit$converged)
  expect_identical(fit$status, "converged_stationary")
  expect_identical(
    fit$certificate$outer_stationarity$classification, "stationary"
  )
  expect_true(fit$certificate$outer_stationarity$certified)
  expect_true(fit$final_block_optimality_certified)
  expect_true(fit$equal_mass_certified)
  expect_true(audit$certified)
  expect_identical(audit$method, "post_rescale_full_support_primal_kkt")
  expect_lt(abs(audit$kkt_residual - max(abs(derivative))), 2e-12)
  expect_equal(audit$block_objective, direct_objective, tolerance = 2e-12)
  expect_lte(audit$kkt_residual, audit$effective_tolerance)
  expect_lte(audit$max_tile_elements, controls$block_size^2)
  expect_lt(audit$max_tile_elements, audit$full_matrix_elements)
  expect_equal(
    audit$workspace_elements,
    2 * controls$block_size^2 + 2 * (length(fixture$wx) + length(fixture$wy))
  )
  expect_match(audit$workspace_scope, "Peak live native doubles")
  expect_identical(fit$global_optimality, "not_claimed")
})

test_that("cheap iterate screening never replaces the final full audit", {
  fixture <- moment_certificate_fixture(seed = 901L)
  controls <- moment_certificate_controls(max_iter_ot = 10L)
  controls$max_iter <- 1L
  controls$tol <- 1e-12
  fit <- do.call(fugw_factorized, c(
    list(
      Cx = fixture$Cx,
      Cy = fixture$Cy,
      M = fixture$M,
      wx = fixture$wx,
      wy = fixture$wy
    ),
    controls
  ))
  blocks <- fit$certificate$outer_stationarity$final_blocks

  expect_false(fit$converged)
  expect_identical(fit$status, "max_iter")
  expect_false(fit$plan_update_audited_trace[[1L]])
  expect_false(fit$outer_plan_residual_exact)
  expect_true(is.infinite(fit$plan_screening_trace[[1L]]))
  expect_gt(blocks$sample$tile_count, 0L)
  expect_gt(blocks$feature$tile_count, 0L)
  expect_true(all(is.finite(c(
    blocks$sample$workspace_elements,
    blocks$feature$workspace_elements
  ))))
  expect_true(
    fit$runtime_provenance$screening_contract$
      full_final_two_sided_kkt_mandatory
  )
  expect_equal(
    fit$runtime_provenance$screening_contract$exact_plan_update_audits,
    0L
  )
  expect_equal(
    fit$runtime_provenance$screening_contract$screened_plan_updates,
    1L
  )
})

test_that("cost perturbation and post-solve rescaling invalidate block KKT", {
  fixture <- moment_certificate_fixture()
  controls <- moment_certificate_controls()
  fit <- moment_certificate_fit(fixture)
  final_cost <- moment_certificate_final_costs(fit, fixture)$sample
  plan <- fit$implicit_plans$sample
  opponent_mass <- transport_plan_mass(fit$implicit_plans$feature)
  terms <- rfugw:::.cost_native_terms(final_cost)
  perturbed_cost <- factorized_cost(
    terms$left,
    terms$right,
    row = terms$row + .01,
    column = terms$column,
    exact = TRUE,
    provenance = list(role = "deliberate_certificate_mutant")
  )
  audit_args <- list(
    epsilon = controls$epsilon * opponent_mass,
    rho = controls$reg_marginals * opponent_mass,
    tolerance = controls$tol_ot,
    block_size = controls$block_size
  )
  perturbed <- do.call(
    rfugw:::.factorized_block_kkt_audit,
    c(list(plan = plan, final_cost = perturbed_cost), audit_args)
  )
  rescaled <- rfugw:::.factorized_rescale_plan(plan, 1.02)
  post_rescale <- do.call(
    rfugw:::.factorized_block_kkt_audit,
    c(
      list(
        plan = rescaled,
        final_cost = final_cost,
        post_solve_scale = 1.02
      ),
      audit_args
    )
  )

  expect_false(perturbed$certified)
  expect_gte(perturbed$maximum_cost_change, .01 - 1e-12)
  expect_gt(perturbed$kkt_residual, perturbed$effective_tolerance)
  expect_false(post_rescale$certified)
  expect_identical(post_rescale$post_solve_scale, 1.02)
  expect_gt(post_rescale$kkt_residual, post_rescale$effective_tolerance)
})

test_that("nonfinite and fully underflowed implicit states fail closed", {
  fixture <- moment_certificate_fixture()
  controls <- moment_certificate_controls()
  fit <- moment_certificate_fit(fixture)
  final_cost <- moment_certificate_final_costs(fit, fixture)$sample
  plan <- fit$implicit_plans$sample
  opponent_mass <- transport_plan_mass(fit$implicit_plans$feature)
  audit <- function(candidate) {
    rfugw:::.factorized_block_kkt_audit(
      candidate,
      final_cost,
      controls$epsilon * opponent_mass,
      controls$reg_marginals * opponent_mass,
      controls$tol_ot,
      controls$block_size
    )
  }

  underflow_state <- rfugw:::.factorized_plan_state(plan)
  underflow_state$source_bar <- underflow_state$source_bar -
    1000 * underflow_state$epsilon
  underflow_state$target_bar <- underflow_state$target_bar -
    1000 * underflow_state$epsilon
  underflow_state$stats <- rfugw:::.factorized_plan_stats(underflow_state)
  underflow <- audit(rfugw:::.factorized_plan_from_state(underflow_state))

  nonfinite_state <- rfugw:::.factorized_plan_state(plan)
  nonfinite_state$source_bar[[1L]] <- Inf
  nonfinite <- audit(rfugw:::.factorized_plan_from_state(nonfinite_state))

  expect_false(underflow$certified)
  expect_false(underflow$numerical_finite)
  expect_false(underflow$support_consistent)
  expect_false(underflow$active_marginals_representable)
  expect_false(nonfinite$certified)
  expect_false(nonfinite$numerical_finite)
})

test_that("inner failure history is retained while a later retry can recover", {
  fixture <- moment_certificate_fixture()
  recovered <- moment_certificate_fit(fixture, max_iter_ot = 3L)
  starved <- moment_certificate_fit(fixture, max_iter_ot = 1L)

  expect_true(recovered$converged)
  expect_true(recovered$inner_uot_certified)
  expect_true(recovered$final_inner_solve_certified)
  expect_gt(recovered$historical_inner_failure_count, 0L)
  expect_false(recovered$historical_all_inner_certified)
  expect_gt(nrow(recovered$certificate$inner_uot$historical_failures), 0L)
  expect_lte(recovered$inner_residual, recovered$feasibility_tolerance)
  expect_gt(
    recovered$historical_max_inner_residual,
    recovered$inner_residual
  )

  expect_false(starved$converged)
  expect_false(starved$inner_uot_certified)
  expect_false(starved$final_inner_solve_certified)
  expect_identical(
    starved$certificate$outer_stationarity$classification, "uncertified"
  )
  expect_false(starved$certificate$outer_stationarity$certified)
  expect_identical(starved$global_optimality, "not_claimed")
})

test_that("the final certificate uses factors and bounded tiles only", {
  fixture <- moment_certificate_fixture(seed = 1702L)
  fixture$Cx$block <- function(...) stop("source structure block evaluated")
  fixture$Cy$block <- function(...) stop("target structure block evaluated")
  fixture$M$block <- function(...) stop("feature block evaluated")
  fit <- moment_certificate_fit(fixture)
  blocks <- fit$certificate$outer_stationarity$final_blocks

  expect_true(fit$converged)
  expect_true(blocks$certified)
  expect_lte(blocks$sample$max_tile_elements, 4L)
  expect_lte(blocks$feature$max_tile_elements, 4L)
  expect_equal(blocks$sample$full_matrix_elements, 20)
  expect_equal(blocks$feature$full_matrix_elements, 20)
  expect_identical(
    transport_plan_representation(fit$plans$sample)$representation,
    "implicit_operator"
  )
})
