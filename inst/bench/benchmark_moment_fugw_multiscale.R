bench_multiscale_group_means <- function(values, parent, weights, n_groups) {
  out <- matrix(0, n_groups, ncol(values))
  totals <- as.numeric(rowsum(weights, parent, reorder = FALSE))
  for (column in seq_len(ncol(values))) {
    out[, column] <- as.numeric(rowsum(
      weights * values[, column], parent, reorder = FALSE
    )) / totals
  }
  out
}

bench_make_moment_fugw_multiscale_case <- function(
    seed, n_groups = 8L, children = 3L) {
  set.seed(seed)
  n_groups <- as.integer(n_groups)
  children <- as.integer(children)
  n_fine <- n_groups * children
  source_parent <- rep(seq_len(n_groups), each = children)
  theta <- 2 * pi * (seq_len(n_groups) - 1) / n_groups
  coarse_latent <- cbind(cos(theta), sin(theta))
  child_angle <- rep(
    2 * pi * (seq_len(children) - 1) / children,
    times = n_groups
  ) + rep(theta, each = children)
  child_radius <- 0.025 * (1 + 0.1 * sin(seq_len(n_fine) + seed))
  source_coordinates <- coarse_latent[source_parent, , drop = FALSE] +
    child_radius * cbind(cos(child_angle), sin(child_angle))
  child_code <- rep(seq(-1, 1, length.out = children), times = n_groups)
  latent_features <- cbind(
    cos(theta[source_parent]),
    sin(theta[source_parent]),
    cos(2 * theta[source_parent]),
    child_code
  )
  source_features <- latent_features + matrix(
    rnorm(n_fine * ncol(latent_features), sd = 0.002),
    n_fine,
    ncol(latent_features)
  )
  heldout_source <- cbind(
    sin(3 * theta[source_parent]) + 0.2 * child_code,
    cos(theta[source_parent] + 0.4 * child_code),
    child_code^2 + 0.1 * sin(theta[source_parent])
  )

  group_latent <- sample(seq_len(n_groups))
  target_latent_index <- unlist(lapply(group_latent, function(group) {
    members <- which(source_parent == group)
    sample(members)
  }), use.names = FALSE)
  target_parent <- rep(seq_len(n_groups), each = children)
  rotation_angle <- runif(1L, -pi, pi)
  rotation <- matrix(
    c(cos(rotation_angle), -sin(rotation_angle),
      sin(rotation_angle), cos(rotation_angle)),
    2L,
    2L,
    byrow = TRUE
  )
  target_coordinates <- source_coordinates[target_latent_index, , drop = FALSE] %*%
    t(rotation)
  target_features <- latent_features[target_latent_index, , drop = FALSE] +
    matrix(
      rnorm(n_fine * ncol(latent_features), sd = 0.002),
      n_fine,
      ncol(latent_features)
    )
  heldout_target <- heldout_source[target_latent_index, , drop = FALSE]
  source_weights <- rep(1 / n_fine, n_fine)
  target_weights <- rep(1 / n_fine, n_fine)
  source_coarse_weights <- rep(1 / n_groups, n_groups)
  target_coarse_weights <- rep(1 / n_groups, n_groups)
  source_coarse_coordinates <- bench_multiscale_group_means(
    source_coordinates, source_parent, source_weights, n_groups
  )
  target_coarse_coordinates <- bench_multiscale_group_means(
    target_coordinates, target_parent, target_weights, n_groups
  )
  source_coarse_features <- bench_multiscale_group_means(
    source_features, source_parent, source_weights, n_groups
  )
  target_coarse_features <- bench_multiscale_group_means(
    target_features, target_parent, target_weights, n_groups
  )

  make_domain <- function(structure, features, weights, parent, name) {
    fugw_domain(
      structure = structure,
      features = features,
      weights = weights,
      parent = parent,
      name = name,
      feature_aggregation_policy = "weighted_mean",
      feature_aggregation_tolerance = 1e-10,
      geometry_aggregation_tolerance = 0.05,
      provenance = list(fixture = "planted_rigid_permutation", seed = seed)
    )
  }
  list(
    seed = seed,
    source = list(
      make_domain(
        sqeuclidean_cost(source_coarse_coordinates),
        source_coarse_features,
        source_coarse_weights,
        NULL,
        "source-coarse"
      ),
      make_domain(
        sqeuclidean_cost(source_coordinates),
        source_features,
        source_weights,
        source_parent,
        "source-fine"
      )
    ),
    target = list(
      make_domain(
        sqeuclidean_cost(target_coarse_coordinates),
        target_coarse_features,
        target_coarse_weights,
        NULL,
        "target-coarse"
      ),
      make_domain(
        sqeuclidean_cost(target_coordinates),
        target_features,
        target_weights,
        target_parent,
        "target-fine"
      )
    ),
    heldout_source = heldout_source,
    heldout_target = heldout_target,
    truth_source_to_target = match(seq_len(n_fine), target_latent_index),
    target_latent_index = target_latent_index
  )
}

bench_moment_fugw_multiscale_controls <- function() {
  list(
    feature_metric = "sqeuclidean",
    feature_weight = 0.7,
    structure_weight = 1,
    reg_marginals = c(8, 8),
    epsilon = c(0.12, 0.06),
    max_iter = c(60L, 80L),
    tol = c(1e-8, 1e-8),
    max_iter_ot = c(5000L, 5000L),
    tol_ot = c(1e-10, 1e-10),
    check_every = 1L,
    block_size = 16L,
    level_contract = c("certified_endpoint", "certified_endpoint")
  )
}

bench_moment_fugw_heldout_score <- function(plan, target_values) {
  source_mass <- transport_plan_mass(plan, "source")
  prediction <- transport_plan_apply(plan, target_values)
  active <- source_mass > sqrt(.Machine$double.eps)
  prediction[active, ] <- prediction[active, , drop = FALSE] /
    source_mass[active]
  list(prediction = prediction, active = active)
}

bench_run_moment_fugw_multiscale_case <- function(seed) {
  problem <- bench_make_moment_fugw_multiscale_case(seed)
  controls <- bench_moment_fugw_multiscale_controls()
  fine_source <- problem$source[[2L]]
  fine_target <- problem$target[[2L]]
  cold_controls <- controls
  cold_controls$feature_metric <- NULL
  cold_controls$epsilon <- tail(controls$epsilon, 1L)
  cold_controls$max_iter <- tail(controls$max_iter, 1L)
  cold_controls$tol <- tail(controls$tol, 1L)
  cold_controls$max_iter_ot <- tail(controls$max_iter_ot, 1L)
  cold_controls$tol_ot <- tail(controls$tol_ot, 1L)
  cold_controls$check_every <- tail(controls$check_every, 1L)
  cold_controls$block_size <- tail(controls$block_size, 1L)
  cold_controls$level_contract <- NULL
  cold <- do.call(fugw_factorized, c(
    list(
      Cx = fine_source$structure,
      Cy = fine_target$structure,
      M = sqeuclidean_cost(fine_source$features, fine_target$features),
      wx = fine_source$weights,
      wy = fine_target$weights
    ),
    cold_controls
  ))
  transferred <- do.call(fugw_multiscale, c(
    list(source = problem$source, target = problem$target),
    controls
  ))
  fine <- transferred$level_results[[2L]]
  cold_heldout <- bench_moment_fugw_heldout_score(
    cold$plans$sample, problem$heldout_target
  )
  transfer_heldout <- bench_moment_fugw_heldout_score(
    fine$plans$sample, problem$heldout_target
  )
  active <- cold_heldout$active & transfer_heldout$active
  cold_rmse <- sqrt(mean(
    (cold_heldout$prediction[active, , drop = FALSE] -
       problem$heldout_source[active, , drop = FALSE])^2
  ))
  transfer_rmse <- sqrt(mean(
    (transfer_heldout$prediction[active, , drop = FALSE] -
       problem$heldout_source[active, , drop = FALSE])^2
  ))
  cold_score <- 1 / (1 + cold_rmse)
  transfer_score <- 1 / (1 + transfer_rmse)
  objective_relative_difference <-
    (fine$fugw_cost - cold$fugw_cost) / max(1, abs(cold$fugw_cost))
  plan_action_relative_difference <- sqrt(sum(
    (transfer_heldout$prediction - cold_heldout$prediction)^2
  )) / max(1, sqrt(sum(cold_heldout$prediction^2)))
  final_controls <- transferred$hierarchy$controls[2L, , drop = FALSE]
  same_final_contract <- isTRUE(all.equal(
    unname(unlist(final_controls[c(
      "epsilon", "max_iter", "tol", "max_iter_ot", "tol_ot",
      "check_every", "block_size", "feature_weight", "structure_weight",
      "rho_source", "rho_target", "rescale_plan"
    )])),
    unname(c(
      cold$regularization,
      cold_controls$max_iter,
      cold_controls$tol,
      cold_controls$max_iter_ot,
      cold_controls$tol_ot,
      cold_controls$check_every,
      cold_controls$block_size,
      cold$feature_weight,
      cold$structure_weight,
      cold$reg_marginals,
      cold$rescale_plan
    )),
    tolerance = 0
  ))
  data.frame(
    seed = seed,
    n_source = fine_source$n,
    n_target = fine_target$n,
    same_final_contract = same_final_contract,
    final_epsilon = cold$regularization,
    final_feature_weight = cold$feature_weight,
    final_structure_weight = cold$structure_weight,
    final_rho_source = cold$reg_marginals[[1L]],
    final_rho_target = cold$reg_marginals[[2L]],
    final_outer_tolerance = cold_controls$tol,
    final_inner_tolerance = cold_controls$tol_ot,
    final_outer_iteration_limit = cold_controls$max_iter,
    final_inner_iteration_limit = cold_controls$max_iter_ot,
    final_block_size = cold_controls$block_size,
    final_support = "full_implicit",
    cold_status = cold$status,
    transfer_status = transferred$status,
    cold_certified = isTRUE(cold$final_stationarity_certified),
    transfer_final_certified = isTRUE(fine$final_stationarity_certified),
    transfer_multiscale_certified = isTRUE(transferred$multiscale_certified),
    cold_outer_iterations = cold$iterations,
    transfer_fine_outer_iterations = fine$iterations,
    cold_inner_iterations = cold$inner_iterations,
    transfer_fine_inner_iterations = fine$inner_iterations,
    inner_iteration_reduction = 1 -
      fine$inner_iterations / cold$inner_iterations,
    cold_objective = cold$fugw_cost,
    transfer_objective = fine$fugw_cost,
    objective_relative_difference = objective_relative_difference,
    cold_heldout_score = cold_score,
    transfer_heldout_score = transfer_score,
    heldout_score_ratio = transfer_score / cold_score,
    plan_action_relative_difference = plan_action_relative_difference,
    cold_plan_residual = cold$outer_plan_residual,
    transfer_plan_residual = fine$outer_plan_residual,
    cold_moment_residual = cold$moment_residual,
    transfer_moment_residual = fine$moment_residual,
    cold_inner_residual = cold$inner_residual,
    transfer_inner_residual = fine$inner_residual,
    stringsAsFactors = FALSE
  )
}

bench_run_moment_fugw_multiscale_efficacy <- function(
    seeds = 1701:1710,
    output = file.path(
      "inst", "bench", "moment-fugw-multiscale-evidence.csv"
    )) {
  rows <- do.call(rbind, lapply(seeds, function(seed) {
    message("Moment-FUGW multiscale efficacy seed ", seed)
    bench_run_moment_fugw_multiscale_case(seed)
  }))
  improved_fraction <- mean(rows$inner_iteration_reduction > 0)
  median_reduction <- stats::median(rows$inner_iteration_reduction)
  rows$improved_fraction <- improved_fraction
  rows$median_inner_iteration_reduction <- median_reduction
  rows$efficacy_gate <-
    all(rows$same_final_contract) &&
    all(rows$cold_certified) &&
    all(rows$transfer_final_certified) &&
    all(rows$transfer_multiscale_certified) &&
    improved_fraction >= 0.8 &&
    median_reduction >= 0.2 &&
    all(rows$objective_relative_difference <= 1e-4) &&
    all(rows$heldout_score_ratio >= 0.99)
  dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(rows, output, row.names = FALSE)
  rows
}
