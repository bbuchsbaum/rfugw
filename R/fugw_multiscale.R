.validate_domain_parent <- function(parent, n) {
  if (is.null(parent)) return(NULL)
  if (!is.numeric(parent) || length(parent) != n ||
      any(!is.finite(parent)) || any(parent != floor(parent)) ||
      any(parent < 1) || any(parent > .Machine$integer.max)) {
    stop(
      "`parent` must contain one positive integer per domain point.",
      call. = FALSE
    )
  }
  as.integer(parent)
}

#' Construct a Matrix-Free FUGW Domain
#'
#' Packages one resolution of a domain for [fugw_multiscale()]. The structure
#' cost must have affine-bilinear factors; it is never materialized. A
#' hierarchy is represented as a coarse-to-fine list of domains. For every
#' level after the first, `parent[i]` identifies the parent of fine point `i`
#' in the preceding level.
#'
#' @param structure Square affine-bilinear [factorized_cost()] operator.
#' @param features Optional finite matrix with one feature profile per row.
#' @param weights Optional nonnegative point weights, normalized to sum to one.
#' @param parent Optional positive integer parent index for each point. Parent
#'   bounds and weight aggregation are checked when a hierarchy is solved.
#' @param name Optional scalar label for the level.
#' @param provenance Optional serializable construction metadata.
#' @param feature_aggregation_policy Whether coarse features are independently
#'   level-specific or must equal conditional weighted means of fine features.
#' @param feature_aggregation_tolerance Relative tolerance when weighted-mean
#'   feature aggregation is required.
#' @param geometry_aggregation_tolerance Maximum weighted stress between each
#'   coarse structure and its conditional fine-level restriction.
#' @param hierarchy_seed Integer seed for reproducible sampled hierarchy audits.
#' @param hierarchy_pair_limit Maximum number of coarse pairs used by a
#'   hierarchy geometry audit; smaller pair sets are enumerated exactly.
#' @return An `rfugw_fugw_domain`.
#' @examples
#' coordinates <- cbind(seq_len(6), rep(0, 6))
#' domain <- fugw_domain(
#'   sqeuclidean_cost(coordinates),
#'   features = coordinates,
#'   weights = rep(1, 6)
#' )
#' @export
fugw_domain <- function(
    structure,
    features = NULL,
    weights = NULL,
    parent = NULL,
    name = NULL,
    provenance = list(),
    feature_aggregation_policy = c("level_specific", "weighted_mean"),
    feature_aggregation_tolerance = 1e-8,
    geometry_aggregation_tolerance = 0.5,
    hierarchy_seed = 130363L,
    hierarchy_pair_limit = 100000L) {
  feature_aggregation_policy <- match.arg(feature_aggregation_policy)
  if (!inherits(structure, "rfugw_cost_operator")) {
    stop("`structure` must be an rfugw cost operator.", call. = FALSE)
  }
  shape <- cost_shape(structure)
  if (shape[["source"]] != shape[["target"]]) {
    stop("`structure` must be square.", call. = FALSE)
  }
  .cost_native_terms(structure, "structure")
  n <- unname(shape[["source"]])
  if (!is.null(features)) {
    features <- .feature_matrix(features, "features")
    if (nrow(features) != n) {
      stop("`features` must have one row per structure point.", call. = FALSE)
    }
  }
  if (is.null(weights)) weights <- rep(1 / n, n)
  weights <- .assert_prob(weights, n, "weights")
  parent <- .validate_domain_parent(parent, n)
  if (!is.null(name) &&
      (!is.character(name) || length(name) != 1L || is.na(name) || !nzchar(name))) {
    stop("`name` must be NULL or a nonempty character scalar.", call. = FALSE)
  }
  if (!is.list(provenance)) {
    stop("`provenance` must be a list.", call. = FALSE)
  }
  feature_aggregation_tolerance <- .validate_positive_scalar(
    feature_aggregation_tolerance, "feature_aggregation_tolerance"
  )
  geometry_aggregation_tolerance <- .validate_positive_scalar(
    geometry_aggregation_tolerance, "geometry_aggregation_tolerance"
  )
  hierarchy_seed <- .validate_count(hierarchy_seed, "hierarchy_seed")
  hierarchy_pair_limit <- .validate_count(
    hierarchy_pair_limit, "hierarchy_pair_limit"
  )
  hierarchy_policy <- list(
    feature_aggregation = feature_aggregation_policy,
    feature_tolerance = feature_aggregation_tolerance,
    geometry_aggregation = "conditional_pair_mean",
    geometry_tolerance = geometry_aggregation_tolerance,
    seed = hierarchy_seed,
    pair_limit = hierarchy_pair_limit
  )
  domain_provenance <- utils::modifyList(
    list(
      representation = "factorized_fugw_domain",
      structure = cost_provenance(structure),
      feature_dimension = if (is.null(features)) 0L else ncol(features),
      has_parent = !is.null(parent),
      hierarchy_policy = hierarchy_policy
    ),
    provenance
  )
  structure(
    list(
      structure = structure,
      features = features,
      weights = weights,
      parent = parent,
      name = name,
      n = n,
      provenance = domain_provenance,
      hierarchy_policy = hierarchy_policy
    ),
    class = "rfugw_fugw_domain"
  )
}

#' @export
print.rfugw_fugw_domain <- function(x, ...) {
  cat("<rfugw_fugw_domain>\n")
  cat(sprintf("  name:             %s\n", x$name %||% "<unnamed>"))
  cat(sprintf("  points:           %d\n", x$n))
  cat(sprintf(
    "  feature columns:  %d\n",
    if (is.null(x$features)) 0L else ncol(x$features)
  ))
  cat(sprintf("  structure rank:   %d\n", ncol(.cost_native_terms(x$structure)$left)))
  cat(sprintf("  parent map:       %s\n", if (is.null(x$parent)) "no" else "yes"))
  invisible(x)
}

.as_fugw_hierarchy <- function(x, name) {
  if (inherits(x, "rfugw_fugw_domain")) return(list(x))
  if (!is.list(x) || !length(x) ||
      any(!vapply(x, inherits, logical(1), "rfugw_fugw_domain"))) {
    stop(
      sprintf("`%s` must be a fugw domain or a coarse-to-fine list of them.", name),
      call. = FALSE
    )
  }
  unname(x)
}

.multiscale_group_sum <- function(values, group, n_groups) {
  out <- numeric(n_groups)
  if (length(values)) {
    summed <- rowsum(values, group, reorder = FALSE)
    out[as.integer(rownames(summed))] <- as.numeric(summed)
  }
  out
}

.hierarchy_conditional_weights <- function(weights, parent, n_coarse) {
  totals <- .multiscale_group_sum(weights, parent, n_coarse)
  counts <- tabulate(parent, nbins = n_coarse)
  conditional <- numeric(length(weights))
  positive <- totals[parent] > 0
  conditional[positive] <- weights[positive] / totals[parent[positive]]
  zero_mass <- !positive
  conditional[zero_mass] <- 1 / counts[parent[zero_mass]]
  list(weights = conditional, totals = totals)
}

.hierarchy_group_means <- function(values, conditional, parent, n_coarse) {
  if (!ncol(values)) return(matrix(0, nrow = n_coarse, ncol = 0L))
  out <- matrix(0, nrow = n_coarse, ncol = ncol(values))
  for (column in seq_len(ncol(values))) {
    out[, column] <- .multiscale_group_sum(
      conditional * values[, column], parent, n_coarse
    )
  }
  out
}

.hierarchy_restricted_cost <- function(
    fine_cost, conditional, parent, n_coarse) {
  terms <- .cost_native_terms(fine_cost, "fine$structure")
  factorized_cost(
    .hierarchy_group_means(terms$left, conditional, parent, n_coarse),
    .hierarchy_group_means(terms$right, conditional, parent, n_coarse),
    row = .multiscale_group_sum(
      conditional * terms$row, parent, n_coarse
    ),
    column = .multiscale_group_sum(
      conditional * terms$column, parent, n_coarse
    ),
    exact = isTRUE(cost_provenance(fine_cost)$exact),
    provenance = list(
      role = "conditional_parent_pair_restriction",
      fine_provenance = cost_provenance(fine_cost)
    )
  )
}

.hierarchy_geometry_agreement <- function(
    coarse, fine, conditional, parent, policy) {
  restricted <- .hierarchy_restricted_cost(
    fine$structure, conditional, parent, coarse$n
  )
  total_pairs <- as.double(coarse$n) * as.double(coarse$n - 1L) / 2
  exact <- total_pairs <= policy$pair_limit
  pairs <- if (exact) {
    .geometry_all_pairs(coarse$n)
  } else {
    .geometry_sample_pairs(
      coarse$n, min(4096L, policy$pair_limit), policy$seed
    )
  }
  coarse_value <- .cost_pair_values(
    coarse$structure, pairs[, 1L], pairs[, 2L]
  )
  restricted_value <- .cost_pair_values(
    restricted, pairs[, 1L], pairs[, 2L]
  )
  error <- .geometry_error_summary(
    coarse_value,
    restricted_value,
    coarse$weights[pairs[, 1L]] * coarse$weights[pairs[, 2L]]
  )
  coarse_diagonal <- .cost_pair_values(
    coarse$structure, seq_len(coarse$n), seq_len(coarse$n)
  )
  restricted_diagonal <- .cost_pair_values(
    restricted, seq_len(coarse$n), seq_len(coarse$n)
  )
  list(
    certified = is.finite(error$weighted_stress) &&
      error$weighted_stress <= policy$geometry_tolerance,
    weighted_stress = error$weighted_stress,
    relative_distance_error = error$relative_distance_error,
    maximum_relative_error = error$maximum_relative_error,
    diagonal_difference = max(abs(coarse_diagonal - restricted_diagonal)),
    tolerance = policy$geometry_tolerance,
    pair_design = list(
      mode = if (exact) "exact_all_coarse_pairs" else {
        "reproducible_sampled_coarse_pairs"
      },
      seed = policy$seed,
      pair_count = nrow(pairs),
      total_available_pairs = total_pairs,
      pairs = if (exact) NULL else pairs
    ),
    restriction = "conditional_parent_pair_mean"
  )
}

.hierarchy_feature_agreement <- function(
    coarse, fine, conditional, parent, policy) {
  both_absent <- is.null(coarse$features) && is.null(fine$features)
  if (both_absent) {
    return(list(
      certified = TRUE,
      policy = policy$feature_aggregation,
      relative_error = 0,
      maximum_error = 0,
      tolerance = policy$feature_tolerance,
      applicable = FALSE
    ))
  }
  compatible <- !is.null(coarse$features) && !is.null(fine$features) &&
    ncol(coarse$features) == ncol(fine$features)
  if (!compatible) {
    equality_required <- identical(
      policy$feature_aggregation, "weighted_mean"
    )
    return(list(
      certified = !equality_required,
      policy = policy$feature_aggregation,
      relative_error = if (equality_required) Inf else NA_real_,
      maximum_error = if (equality_required) Inf else NA_real_,
      tolerance = policy$feature_tolerance,
      applicable = equality_required,
      compatible = FALSE,
      equality_required = equality_required
    ))
  }
  restricted <- .hierarchy_group_means(
    fine$features, conditional, parent, coarse$n
  )
  difference <- coarse$features - restricted
  relative_error <- sqrt(sum(difference^2)) /
    max(1, sqrt(sum(restricted^2)))
  maximum_error <- max(abs(difference))
  equality_required <- identical(
    policy$feature_aggregation, "weighted_mean"
  )
  list(
    certified = !equality_required ||
      relative_error <= policy$feature_tolerance,
    policy = policy$feature_aggregation,
    relative_error = relative_error,
    maximum_error = maximum_error,
    tolerance = policy$feature_tolerance,
    applicable = TRUE,
    compatible = TRUE,
    equality_required = equality_required
  )
}

.hierarchy_transition <- function(coarse, fine, name, level, tolerance) {
  parent <- fine$parent
  if (is.null(parent)) {
    stop(
      sprintf("`%s[[%d]]` must define a parent map to level %d.",
              name, level, level - 1L),
      call. = FALSE
    )
  }
  if (any(parent > coarse$n)) {
    stop(
      sprintf("`%s[[%d]]$parent` exceeds the preceding level size.", name, level),
      call. = FALSE
    )
  }
  child_counts <- tabulate(parent, nbins = coarse$n)
  if (any(child_counts == 0L)) {
    stop(
      sprintf("Every `%s` point at level %d must have a child at level %d.",
              name, level - 1L, level),
      call. = FALSE
    )
  }
  aggregated <- .multiscale_group_sum(fine$weights, parent, coarse$n)
  residual <- max(abs(aggregated - coarse$weights))
  if (!is.finite(residual) || residual > tolerance) {
    stop(
      sprintf(
        "`%s` weights do not aggregate through the level-%d parent map (residual %.3g > %.3g).",
        name, level, residual, tolerance
      ),
      call. = FALSE
    )
  }
  policy <- fine$hierarchy_policy
  conditional <- .hierarchy_conditional_weights(
    fine$weights, parent, coarse$n
  )$weights
  conditional_sums <- .multiscale_group_sum(
    conditional, parent, coarse$n
  )
  round_trip_residual <- max(abs(conditional_sums - 1))
  mass_prolongation_residual <- max(abs(
    coarse$weights[parent] * conditional - fine$weights
  ))
  round_trip_certified <- is.finite(round_trip_residual) &&
    is.finite(mass_prolongation_residual) &&
    max(round_trip_residual, mass_prolongation_residual) <= tolerance
  feature_agreement <- .hierarchy_feature_agreement(
    coarse, fine, conditional, parent, policy
  )
  geometry_agreement <- .hierarchy_geometry_agreement(
    coarse, fine, conditional, parent, policy
  )
  if (!round_trip_certified) {
    stop(
      sprintf(
        "`%s` level-%d restrict-prolong round trip is inconsistent.",
        name, level
      ),
      call. = FALSE
    )
  }
  if (!isTRUE(feature_agreement$certified)) {
    stop(
      sprintf(
        "`%s` level-%d features violate the declared weighted-mean aggregation policy (relative error %.3g > %.3g).",
        name, level, feature_agreement$relative_error,
        feature_agreement$tolerance
      ),
      call. = FALSE
    )
  }
  if (!isTRUE(geometry_agreement$certified)) {
    stop(
      sprintf(
        "`%s` level-%d coarse geometry disagrees with conditional fine geometry (weighted stress %.3g > %.3g).",
        name, level, geometry_agreement$weighted_stress,
        geometry_agreement$tolerance
      ),
      call. = FALSE
    )
  }
  list(
    certified = TRUE,
    parent = parent,
    parent_coverage = list(
      all_coarse_points_covered = all(child_counts > 0L),
      child_count_min = min(child_counts),
      child_count_max = max(child_counts)
    ),
    aggregated_weights = aggregated,
    weight_residual = residual,
    weight_tolerance = tolerance,
    child_count_min = min(child_counts),
    child_count_max = max(child_counts),
    restrict_prolong = list(
      certified = round_trip_certified,
      value_round_trip_residual = round_trip_residual,
      mass_prolongation_residual = mass_prolongation_residual,
      tolerance = tolerance
    ),
    feature_aggregation = feature_agreement,
    geometry_agreement = geometry_agreement,
    transfer_error = max(
      residual,
      round_trip_residual,
      mass_prolongation_residual,
      if (is.finite(feature_agreement$relative_error)) {
        feature_agreement$relative_error
      } else 0,
      geometry_agreement$weighted_stress
    )
  )
}

.validate_fugw_hierarchy <- function(levels, name, tolerance) {
  if (!is.null(levels[[1L]]$parent)) {
    stop(sprintf("The coarsest `%s` domain must not define `parent`.", name),
         call. = FALSE)
  }
  if (length(levels) == 1L) return(list())
  lapply(seq.int(2L, length(levels)), function(level) {
    .hierarchy_transition(
      levels[[level - 1L]], levels[[level]], name, level, tolerance
    )
  })
}

.multiscale_schedule <- function(x, n_levels, name, integer = FALSE) {
  if (!is.numeric(x) || !length(x) ||
      !(length(x) %in% c(1L, n_levels))) {
    stop(
      sprintf("`%s` must have length one or the number of hierarchy levels.", name),
      call. = FALSE
    )
  }
  if (length(x) == 1L) x <- rep(x, n_levels)
  if (integer) {
    return(vapply(x, .validate_count, integer(1), name = name))
  }
  vapply(x, .validate_positive_scalar, numeric(1), name = name)
}

.multiscale_feature_cost <- function(source, target, metric, level) {
  if (is.null(source$features) && is.null(target$features)) {
    return(.zero_factorized_cost(
      source$n, target$n,
      list(role = "multiscale_feature_cost", level = level, metric = "none")
    ))
  }
  if (xor(is.null(source$features), is.null(target$features))) {
    stop(
      sprintf("Both domains must provide features at hierarchy level %d.", level),
      call. = FALSE
    )
  }
  if (ncol(source$features) != ncol(target$features)) {
    stop(
      sprintf("Source and target feature dimensions differ at level %d.", level),
      call. = FALSE
    )
  }
  provenance <- list(role = "multiscale_feature_cost", level = level)
  switch(
    metric,
    sqeuclidean = sqeuclidean_cost(
      source$features, target$features, provenance = provenance
    ),
    correlation = correlation_cost(
      source$features, target$features, provenance = provenance
    ),
    cosine = cosine_cost(
      source$features, target$features, provenance = provenance
    )
  )
}

.conditional_child_weights <- function(fine_weights, parent, n_coarse) {
  totals <- .multiscale_group_sum(fine_weights, parent, n_coarse)
  out <- numeric(length(fine_weights))
  active <- totals[parent] > 0
  out[active] <- fine_weights[active] / totals[parent[active]]
  list(weights = out, totals = totals)
}

.conditional_group_means <- function(values, conditional, parent, n_coarse) {
  if (!ncol(values)) return(matrix(0, nrow = n_coarse, ncol = 0L))
  out <- matrix(0, nrow = n_coarse, ncol = ncol(values))
  for (column in seq_len(ncol(values))) {
    out[, column] <- .multiscale_group_sum(
      conditional * values[, column], parent, n_coarse
    )
  }
  out
}

.multiscale_recenter_potentials <- function(
    source_bar, target_bar, source_weights, target_weights) {
  source_mean <- sum(source_weights * source_bar)
  target_mean <- sum(target_weights * target_bar)
  shift <- 0.5 * (target_mean - source_mean)
  source_centered <- as.numeric(source_bar + shift)
  target_centered <- as.numeric(target_bar - shift)
  centered_source_mean <- sum(source_weights * source_centered)
  centered_target_mean <- sum(target_weights * target_centered)
  list(
    source_bar = source_centered,
    target_bar = target_centered,
    audit = list(
      method = "equal_weighted_mean_gauge",
      shift = shift,
      source_mean_before = source_mean,
      target_mean_before = target_mean,
      source_mean_after = centered_source_mean,
      target_mean_after = centered_target_mean,
      centered_mean_residual = abs(
        centered_source_mean - centered_target_mean
      ),
      log_kernel_invariance_residual = abs(
        (source_centered[[1L]] - source_bar[[1L]]) +
          (target_centered[[1L]] - target_bar[[1L]])
      ),
      certified = all(is.finite(c(
        source_centered,
        target_centered,
        centered_source_mean,
        centered_target_mean
      ))) && abs(centered_source_mean - centered_target_mean) <=
        100 * .Machine$double.eps * max(
          1, abs(centered_source_mean), abs(centered_target_mean)
        )
    )
  )
}

.lifted_product_kl <- function(
    entropy_sum, source_marginal, target_marginal, source, target, mass) {
  if (any(source_marginal > 0 & source == 0) ||
      any(target_marginal > 0 & target == 0)) {
    return(Inf)
  }
  source_active <- source_marginal > 0
  target_active <- target_marginal > 0
  value <- entropy_sum -
    sum(source_marginal[source_active] * log(source[source_active])) -
    sum(target_marginal[target_active] * log(target[target_active])) -
    mass + sum(source) * sum(target)
  tolerance <- 1e-10 * max(1, abs(entropy_sum), mass)
  if (value < -tolerance) {
    stop("Lifted plan has a numerically invalid product KL divergence.",
         call. = FALSE)
  }
  max(0, value)
}

.lift_plan_summary <- function(
    coarse_plan,
    coarse_potentials,
    source_coarse,
    target_coarse,
    source_fine,
    target_fine) {
  state <- .factorized_plan_state(coarse_plan, "coarse_plan")
  source_parent <- source_fine$parent
  target_parent <- target_fine$parent
  source_conditional <- .conditional_child_weights(
    source_fine$weights, source_parent, source_coarse$n
  )$weights
  target_conditional <- .conditional_child_weights(
    target_fine$weights, target_parent, target_coarse$n
  )$weights
  source_marginal <-
    state$stats$source_marginal[source_parent] * source_conditional
  target_marginal <-
    state$stats$target_marginal[target_parent] * target_conditional
  source_active <- source_conditional > 0 & source_marginal > 0
  target_active <- target_conditional > 0 & target_marginal > 0
  entropy_sum <- state$stats$entropy_sum +
    sum(source_marginal[source_active] * log(source_conditional[source_active])) +
    sum(target_marginal[target_active] * log(target_conditional[target_active]))

  source_right <- .cost_factor_pair(
    source_fine$structure, "source_fine$structure"
  )$right
  target_right <- .cost_factor_pair(
    target_fine$structure, "target_fine$structure"
  )$right
  source_means <- .conditional_group_means(
    source_right, source_conditional, source_parent, source_coarse$n
  )
  target_means <- .conditional_group_means(
    target_right, target_conditional, target_parent, target_coarse$n
  )
  cross <- if (!ncol(source_right) || !ncol(target_right)) {
    matrix(0, nrow = ncol(source_right), ncol = ncol(target_right))
  } else {
    crossprod(
      source_means,
      transport_plan_apply(coarse_plan, target_means)
    )
  }
  source_gram <- if (!ncol(source_right)) {
    matrix(0, 0L, 0L)
  } else {
    crossprod(source_right, source_right * source_marginal)
  }
  target_gram <- if (!ncol(target_right)) {
    matrix(0, 0L, 0L)
  } else {
    crossprod(target_right, target_right * target_marginal)
  }
  mass <- state$stats$mass
  raw_source_bar <- as.numeric(
    coarse_potentials$source_bar[source_parent]
  )
  raw_target_bar <- as.numeric(
    coarse_potentials$target_bar[target_parent]
  )
  recentered <- .multiscale_recenter_potentials(
    raw_source_bar,
    raw_target_bar,
    source_fine$weights,
    target_fine$weights
  )
  stats <- list(
    source_marginal = as.numeric(source_marginal),
    target_marginal = as.numeric(target_marginal),
    mass = mass,
    entropy_sum = entropy_sum,
    plan_product_kl = .lifted_product_kl(
      entropy_sum,
      source_marginal,
      target_marginal,
      source_fine$weights,
      target_fine$weights,
      mass
    )
  )
  moments <- list(
    source_mass = stats$source_marginal,
    target_mass = stats$target_marginal,
    H = unname(cross),
    G_source = unname(source_gram),
    G_target = unname(target_gram),
    transported_mass = mass
  )
  list(
    stats = stats,
    moments = moments,
    potentials = list(
      source_bar = recentered$source_bar,
      target_bar = recentered$target_bar
    ),
    gauge = recentered$audit
  )
}

.lift_fugw_initialization <- function(
    coarse_fit,
    source_coarse,
    target_coarse,
    source_fine,
    target_fine,
    source_transition,
    target_transition,
    from_level,
    to_level,
    from_epsilon = NULL,
    to_epsilon = NULL) {
  sample <- .lift_plan_summary(
    coarse_fit$implicit_plans$sample,
    coarse_fit$potentials$sample,
    source_coarse,
    target_coarse,
    source_fine,
    target_fine
  )
  feature <- .lift_plan_summary(
    coarse_fit$implicit_plans$feature,
    coarse_fit$potentials$feature,
    source_coarse,
    target_coarse,
    source_fine,
    target_fine
  )
  provenance <- list(
    origin = "multiscale_conditional_lift",
    from_level = from_level,
    to_level = to_level,
    potentials_prolonged = TRUE,
    potentials_gauge_recentered = TRUE,
    sample_gauge = sample$gauge,
    feature_gauge = feature$gauge,
    moments_recomputed_in_fine_basis = TRUE,
    mass_preserved = TRUE,
    source_weight_aggregation_residual = source_transition$weight_residual,
    target_weight_aggregation_residual = target_transition$weight_residual,
    source_hierarchy_transfer_error = source_transition$transfer_error,
    target_hierarchy_transfer_error = target_transition$transfer_error,
    source_geometry_agreement = source_transition$geometry_agreement,
    target_geometry_agreement = target_transition$geometry_agreement,
    source_restrict_prolong = source_transition$restrict_prolong,
    target_restrict_prolong = target_transition$restrict_prolong,
    lift = "conditional_product_within_parent_pairs",
    epsilon_schedule = list(
      from = from_epsilon,
      to = to_epsilon,
      changed = !is.null(from_epsilon) && !is.null(to_epsilon) &&
        !isTRUE(all.equal(from_epsilon, to_epsilon)),
      potential_units = "cost",
      potential_rescaling = "none"
    )
  )
  list(sample = sample, feature = feature, provenance = provenance)
}

.multiscale_level_contract <- function(level_contract, n_levels) {
  choices <- c("certified_endpoint", "budgeted_warm_start")
  if (is.null(level_contract)) {
    return(rep("certified_endpoint", n_levels))
  }
  if (!is.character(level_contract) ||
      !(length(level_contract) %in% c(1L, n_levels)) ||
      any(is.na(level_contract)) || any(!level_contract %in% choices)) {
    stop(
      paste0(
        "`level_contract` must contain `certified_endpoint` or ",
        "`budgeted_warm_start`, with length one or the hierarchy depth."
      ),
      call. = FALSE
    )
  }
  if (length(level_contract) == 1L) {
    level_contract <- rep(level_contract, n_levels)
  }
  if (!identical(level_contract[[n_levels]], "certified_endpoint")) {
    stop("The finest level must be a `certified_endpoint`.", call. = FALSE)
  }
  unname(level_contract)
}

#' Matrix-Free Multiscale Fused Unbalanced Gromov-Wasserstein
#'
#' Solves [fugw_factorized()] over an arbitrary-depth coarse-to-fine hierarchy.
#' Between levels, both translation-invariant potential pairs are prolonged
#' through explicit parent maps. The coarse implicit couplings are lifted by
#' within-parent conditional weights, which preserves transported mass and
#' exactly recomputes their marginal, entropy, cross-moment, and Gram-moment
#' state in the fine structure basis. No coupling or cost matrix is assembled.
#'
#' The current implementation uses complete implicit support. Adaptive sparse
#' support is deliberately rejected until omitted-mass or reduced-cost bounds
#' are implemented.
#'
#' @section Computational contract:
#' Every level inherits the matrix-free memory contract of
#' [fugw_factorized()]: factors, potentials, marginals, small moments, and
#' reusable tiles are retained, but no full cost or coupling is stored. Fine
#' levels still perform full-support rank-dependent pair arithmetic, typically
#' `O(n_source * n_target * r)` per sweep. Multiscale transfer is a warm-start
#' mechanism and does not change that asymptotic arithmetic claim.
#'
#' @section Certificate layers:
#' The returned object keeps inner UOT and final two-sided stationarity under
#' `certificate`, geometry and parent-map/aggregation evidence under
#' `certificate$multiscale`, and a complete-support omitted-mass bound of zero.
#' Scientific validation is external benchmark evidence, not a per-fit
#' certificate. The method is experimental and never claims global optimality.
#'
#' @param source,target A [fugw_domain()] or an equal-depth coarse-to-fine list
#'   of domains.
#' @param feature_metric Feature-profile cost constructed independently at each
#'   level.
#' @param feature_weight,structure_weight Nonnegative FUGW coefficients, not
#'   normalized shares. They cannot both be zero.
#' @param reg_marginals One or two positive marginal KL penalties.
#' @param epsilon Positive scalar or one value per hierarchy level.
#' @param max_iter,tol,max_iter_ot,tol_ot,check_every,block_size Scalar or
#'   per-level schedules forwarded to [fugw_factorized()].
#' @param rescale_plan Whether to equalize the two alternating coupling masses.
#' @param plan Return the final couplings as implicit operators or explicitly
#'   materialized dense matrices. Intermediate levels always stay implicit.
#' @param support Only `"full"` is implemented. `"adaptive"` fails closed.
#' @param certify Retain complete per-level and layered certificates.
#' @param level_contract Optional per-level character schedule. A
#'   `"certified_endpoint"` must reach its full numerical and stationarity
#'   certificate. A `"budgeted_warm_start"` is retained as finite provenance
#'   but cannot support an all-level stationary status. The finest level must
#'   always be a certified endpoint. The default requires every level to be
#'   certified.
#' @param hierarchy_tolerance Maximum absolute weight-aggregation residual for
#'   every parent transition.
#' @return An `rfugw_result` whose primary plan is the finest sample coupling.
#'   `level_results`, `level_trace`, and `transfers` retain multiscale evidence.
#' @examples
#' coarse_x <- matrix(c(0, 0, 1, 0), ncol = 2, byrow = TRUE)
#' fine_x <- matrix(c(0, 0, 0.4, 0, 0.8, 0, 1.2, 0), ncol = 2, byrow = TRUE)
#' source <- list(
#'   fugw_domain(sqeuclidean_cost(coarse_x)),
#'   fugw_domain(
#'     sqeuclidean_cost(fine_x), parent = c(1, 1, 2, 2)
#'   )
#' )
#' fit <- fugw_multiscale(
#'   source, source, epsilon = c(0.2, 0.1),
#'   max_iter = 3, max_iter_ot = 200
#' )
#' certificate_layers <- list(
#'   inner_uot = fit$certificate$inner_uot$certified,
#'   final_two_sided_stationarity =
#'     fit$certificate$outer_stationarity$certified,
#'   geometry = fit$certificate$multiscale$all_geometry_certified,
#'   hierarchy = fit$certificate$multiscale$hierarchy_certified,
#'   support = fit$certificate$multiscale$support,
#'   scientific_validation = "external_frozen_benchmark_not_per_fit"
#' )
#' certificate_layers
#' @export
fugw_multiscale <- function(
    source,
    target,
    feature_metric = c("sqeuclidean", "correlation", "cosine"),
    feature_weight = 1,
    structure_weight = 1,
    reg_marginals = c(10, 10),
    epsilon = 1e-2,
    max_iter = 100L,
    tol = 1e-7,
    max_iter_ot = 2000L,
    tol_ot = 1e-8,
    rescale_plan = TRUE,
    check_every = 1L,
    plan = c("operator", "dense"),
    support = c("full", "adaptive"),
    block_size = 256L,
    certify = TRUE,
    level_contract = NULL,
    hierarchy_tolerance = 1e-10) {
  feature_metric <- match.arg(feature_metric)
  plan <- match.arg(plan)
  support <- match.arg(support)
  if (identical(support, "adaptive")) {
    stop(
      paste0(
        "`support = \"adaptive\"` is not implemented: use complete implicit ",
        "support until omitted-mass or reduced-cost certification is available."
      ),
      call. = FALSE
    )
  }
  if (!is.logical(certify) || length(certify) != 1L || is.na(certify)) {
    stop("`certify` must be TRUE or FALSE.", call. = FALSE)
  }
  hierarchy_tolerance <- .validate_positive_scalar(
    hierarchy_tolerance, "hierarchy_tolerance"
  )
  source_levels <- .as_fugw_hierarchy(source, "source")
  target_levels <- .as_fugw_hierarchy(target, "target")
  n_levels <- length(source_levels)
  if (length(target_levels) != n_levels) {
    stop("Source and target hierarchies must have equal depth.", call. = FALSE)
  }
  source_transitions <- .validate_fugw_hierarchy(
    source_levels, "source", hierarchy_tolerance
  )
  target_transitions <- .validate_fugw_hierarchy(
    target_levels, "target", hierarchy_tolerance
  )
  epsilon <- .multiscale_schedule(epsilon, n_levels, "epsilon")
  max_iter <- .multiscale_schedule(
    max_iter, n_levels, "max_iter", integer = TRUE
  )
  tol <- .multiscale_schedule(tol, n_levels, "tol")
  max_iter_ot <- .multiscale_schedule(
    max_iter_ot, n_levels, "max_iter_ot", integer = TRUE
  )
  tol_ot <- .multiscale_schedule(tol_ot, n_levels, "tol_ot")
  check_every <- .multiscale_schedule(
    check_every, n_levels, "check_every", integer = TRUE
  )
  block_size <- .multiscale_schedule(
    block_size, n_levels, "block_size", integer = TRUE
  )
  level_contract <- .multiscale_level_contract(level_contract, n_levels)

  level_results <- vector("list", n_levels)
  transfers <- vector("list", max(0L, n_levels - 1L))
  init <- NULL
  for (level in seq_len(n_levels)) {
    source_level <- source_levels[[level]]
    target_level <- target_levels[[level]]
    feature_cost <- .multiscale_feature_cost(
      source_level, target_level, feature_metric, level
    )
    level_plan <- if (level == n_levels) plan else "operator"
    fit <- fugw_factorized(
      source_level$structure,
      target_level$structure,
      wx = source_level$weights,
      wy = target_level$weights,
      reg_marginals = reg_marginals,
      epsilon = epsilon[[level]],
      M = feature_cost,
      max_iter = max_iter[[level]],
      tol = tol[[level]],
      max_iter_ot = max_iter_ot[[level]],
      tol_ot = tol_ot[[level]],
      rescale_plan = rescale_plan,
      check_every = check_every[[level]],
      plan = level_plan,
      block_size = block_size[[level]],
      certify = certify,
      init = init,
      feature_weight = feature_weight,
      structure_weight = structure_weight
    )
    level_results[[level]] <- fit
    if (level < n_levels) {
      init <- .lift_fugw_initialization(
        fit,
        source_level,
        target_level,
        source_levels[[level + 1L]],
        target_levels[[level + 1L]],
        source_transitions[[level]],
        target_transitions[[level]],
        level,
        level + 1L,
        epsilon[[level]],
        epsilon[[level + 1L]]
      )
      transfers[[level]] <- init$provenance
    }
  }

  level_trace <- do.call(rbind, lapply(seq_len(n_levels), function(level) {
    result <- level_results[[level]]
    data.frame(
      level = level,
      source_size = source_levels[[level]]$n,
      target_size = target_levels[[level]]$n,
      epsilon = epsilon[[level]],
      objective = result$fugw_cost,
      iterations = result$iterations,
      inner_iterations = result$inner_iterations,
      residual = result$residual,
      plan_update_residual = result$outer_plan_residual,
      plan_update_residual_exact = isTRUE(result$outer_plan_residual_exact),
      moment_residual = result$moment_residual,
      max_inner_residual = result$max_inner_residual,
      converged = isTRUE(result$converged),
      inner_uot_certified = isTRUE(result$inner_uot_certified),
      outer_stationarity_certified = isTRUE(
        result$certificate$outer_stationarity$certified
      ),
      stationarity_classification =
        result$certificate$outer_stationarity$classification,
      geometry_certified = isTRUE(result$geometry_certified),
      geometry_status = result$geometry_status,
      geometry_exact = isTRUE(result$geometry_exact),
      geometry_relative_error = result$geometry_relative_error,
      geometry_heldout_error = result$geometry_heldout_error,
      initialization = result$initialization$mode,
      level_contract = level_contract[[level]],
      budgeted_admissible = is.finite(result$fugw_cost) &&
        isTRUE(result$certificate$outer_stationarity$numerical_finite) &&
        !(result$status %in% c("numerical_failure", "objective_mismatch")),
      status = result$status,
      stringsAsFactors = FALSE
    )
  }))
  level_trace$contract_satisfied <- ifelse(
    level_trace$level_contract == "certified_endpoint",
    level_trace$inner_uot_certified &
      level_trace$outer_stationarity_certified,
    level_trace$budgeted_admissible
  )
  final <- level_results[[n_levels]]
  all_inner <- all(level_trace$inner_uot_certified)
  all_stationary <- all(level_trace$outer_stationarity_certified)
  all_geometry_certified <- all(level_trace$geometry_certified)
  all_geometry_exact <- all(level_trace$geometry_exact)
  transition_audits <- c(source_transitions, target_transitions)
  hierarchy_certified <- if (length(transition_audits)) {
    all(vapply(transition_audits, function(x) isTRUE(x$certified), logical(1)))
  } else TRUE
  hierarchy_errors <- if (length(transition_audits)) {
    vapply(transition_audits, `[[`, numeric(1), "transfer_error")
  } else 0
  hierarchy_transfer_error <- max(hierarchy_errors)
  geometry_errors <- level_trace$geometry_relative_error
  maximum_geometry_error <- if (any(is.finite(geometry_errors))) {
    max(geometry_errors[is.finite(geometry_errors)])
  } else {
    NA_real_
  }
  final$formulation <- "fugw_multiscale_joint_kl"
  final$backend <- "cpp_blocked_affine_bilinear_ti_multiscale"
  final$level_results <- level_results
  final$level_trace <- level_trace
  final$transfers <- transfers
  final$hierarchy <- list(
    depth = n_levels,
    source_sizes = vapply(source_levels, `[[`, integer(1), "n"),
    target_sizes = vapply(target_levels, `[[`, integer(1), "n"),
    feature_metric = feature_metric,
    epsilon = epsilon,
    level_contract = level_contract,
    controls = data.frame(
      level = seq_len(n_levels),
      epsilon = epsilon,
      max_iter = max_iter,
      tol = tol,
      max_iter_ot = max_iter_ot,
      tol_ot = tol_ot,
      check_every = check_every,
      block_size = block_size,
      level_contract = level_contract,
      feature_weight = rep(feature_weight, n_levels),
      structure_weight = rep(structure_weight, n_levels),
      rho_source = rep(reg_marginals[[1L]], n_levels),
      rho_target = rep(reg_marginals[[min(2L, length(reg_marginals))]], n_levels),
      rescale_plan = rep(isTRUE(rescale_plan), n_levels),
      feature_metric = rep(feature_metric, n_levels),
      support = rep("full_implicit", n_levels),
      stringsAsFactors = FALSE
    ),
    support = "full_implicit",
    hierarchy_tolerance = hierarchy_tolerance
  )
  final_level_certified <- isTRUE(
    final$certificate$outer_stationarity$certified
  )
  final$hierarchy_transfer_error <- hierarchy_transfer_error
  final$hierarchy_transfer_certified <- hierarchy_certified
  final$multiscale_certified <- all_inner && all_stationary &&
    all_geometry_certified && hierarchy_certified
  execution_contract_certified <- all(level_trace$contract_satisfied) &&
    all_geometry_certified && hierarchy_certified
  has_budgeted_warm_starts <- any(
    level_contract[-n_levels] == "budgeted_warm_start"
  )
  multiscale_classification <- if (!final_level_certified) {
    "uncertified"
  } else if (!isTRUE(final$geometry_certified)) {
    "stationary_geometry_uncertified"
  } else if (!all_geometry_certified) {
    "multiscale_geometry_uncertified"
  } else if (!hierarchy_certified) {
    "hierarchy_transfer_uncertified"
  } else if (has_budgeted_warm_starts && !execution_contract_certified) {
    "multiscale_contract_failure"
  } else if (has_budgeted_warm_starts) {
    "final_level_with_budgeted_warm_starts"
  } else if (!all_inner || !all_stationary) {
    "final_level_only"
  } else if (!isTRUE(final$multiscale_certified)) {
    "final_level_only"
  } else if (!all_geometry_exact) {
    "approximate_stationary"
  } else {
    "stationary"
  }
  final$certificate$outer_stationarity$classification <-
    multiscale_classification
  final$certificate$outer_stationarity$scope <- if (
      identical(multiscale_classification, "final_level_only")) {
    "final_level_only"
  } else if (identical(
      multiscale_classification,
      "final_level_with_budgeted_warm_starts")) {
    "final_level_with_budgeted_warm_starts"
  } else {
    "multiscale"
  }
  if (identical(multiscale_classification, "stationary")) {
    final$status <- "converged_stationary"
    final$warning_payload <- NULL
  } else if (identical(
      multiscale_classification, "approximate_stationary")) {
    final$status <- "converged_approximate_stationary"
    final$warning_payload <- NULL
  } else if (identical(multiscale_classification, "final_level_only")) {
    final$status <- "converged_final_level_only"
    final$termination_reason <- "final_level_tolerance_only"
    final$warning_payload <- list(
      code = "final_level_only",
      message = paste0(
        "The finest Moment-FUGW level is stationary, but at least one ",
        "coarser hierarchy level is not certified."
      )
    )
  } else if (identical(
      multiscale_classification,
      "final_level_with_budgeted_warm_starts")) {
    final$status <- "converged_final_level_with_budgeted_warm_starts"
    final$termination_reason <- "final_level_certified_budgeted_warm_starts"
    final$warning_payload <- list(
      code = "budgeted_warm_starts",
      message = paste0(
        "The finest Moment-FUGW level is stationary and every declared ",
        "execution contract passed; budgeted intermediate levels are not ",
        "claimed as stationary endpoints."
      )
    )
  } else if (identical(
      multiscale_classification, "multiscale_contract_failure")) {
    final$status <- "converged_uncertified_multiscale_contract"
    final$termination_reason <- "multiscale_contract_failure"
    final$warning_payload <- list(
      code = "multiscale_contract_failure",
      message = paste0(
        "The finest Moment-FUGW level is stationary, but at least one ",
        "declared intermediate-level execution contract failed."
      )
    )
  } else if (multiscale_classification %in% c(
      "stationary_geometry_uncertified",
      "multiscale_geometry_uncertified")) {
    final$status <- "converged_uncertified_geometry"
    final$warning_payload <- list(
      code = "geometry_uncertified",
      message = paste0(
        "The represented finest-level objective is stationary, but the ",
        "requested geometry is not certified at every required level."
      )
    )
  } else if (identical(
      multiscale_classification, "hierarchy_transfer_uncertified")) {
    final$status <- "converged_uncertified_hierarchy"
    final$warning_payload <- list(
      code = "hierarchy_transfer_uncertified",
      message = "The finest-level objective is stationary but hierarchy transfer is not certified."
    )
  }
  final$certificate$multiscale <- list(
    certified = isTRUE(final$multiscale_certified),
    execution_contract_certified = execution_contract_certified,
    classification = multiscale_classification,
    final_level_stationary = final_level_certified,
    all_inner_uot_certified = all_inner,
    all_levels_stationary = all_stationary,
    has_budgeted_warm_starts = has_budgeted_warm_starts,
    level_contract = level_contract,
    level_contract_satisfied = level_trace$contract_satisfied,
    all_geometry_certified = all_geometry_certified,
    all_geometry_exact = all_geometry_exact,
    maximum_geometry_relative_error = maximum_geometry_error,
    hierarchy_certified = hierarchy_certified,
    hierarchy_transfer_error = hierarchy_transfer_error,
    hierarchy = list(
      source = source_transitions,
      target = target_transitions
    ),
    level_trace = level_trace,
    transfers = transfers,
    support = list(
      mode = "full_implicit",
      complete = TRUE,
      omitted_kernel_mass_bound = 0,
      adaptive_refinement = FALSE
    ),
    global_optimality = "not_claimed"
  )
  final$runtime_provenance$multiscale <- list(
    depth = n_levels,
    intermediate_dense_plans = FALSE,
    final_dense_plan_requested = identical(plan, "dense"),
    transferred_state = c("potentials", "mass", "marginals", "H", "G"),
    potential_gauge = "equal_weighted_mean",
    epsilon_transfer = "potentials_preserved_in_cost_units",
    level_contract = level_contract
  )
  final
}
