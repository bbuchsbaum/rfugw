# Independent tiny-problem oracle for the joint-KL two-coupling FUGW
# objective. This file deliberately contains no calls to rfugw objective,
# moment, or solver helpers. Its O(ns^2 * nt^2) loops are intended only for
# small certification fixtures.

moment_fugw_oracle_scalar_gkl <- function(value, reference) {
  if (!is.numeric(value) || length(value) != 1L || !is.finite(value) ||
      value < 0 || !is.numeric(reference) || length(reference) != 1L ||
      !is.finite(reference) || reference < 0) {
    stop("Generalized-KL scalar inputs must be finite and nonnegative.",
         call. = FALSE)
  }
  if (value == 0) return(reference)
  if (reference == 0) return(Inf)
  value * log(value / reference) - value + reference
}

moment_fugw_oracle_product_gkl <- function(
    left, right, reference_left, reference_right, mass_terms = TRUE) {
  left <- as.numeric(left)
  right <- as.numeric(right)
  reference_left <- as.numeric(reference_left)
  reference_right <- as.numeric(reference_right)
  if (length(left) != length(reference_left) ||
      length(right) != length(reference_right) ||
      any(!is.finite(c(left, right, reference_left, reference_right))) ||
      any(c(left, right, reference_left, reference_right) < 0)) {
    stop("Product-KL inputs must be compatible, finite, and nonnegative.",
         call. = FALSE)
  }

  total <- 0
  for (a in seq_along(left)) {
    for (b in seq_along(right)) {
      value <- left[[a]] * right[[b]]
      reference <- reference_left[[a]] * reference_right[[b]]
      contribution <- if (isTRUE(mass_terms)) {
        moment_fugw_oracle_scalar_gkl(value, reference)
      } else if (value == 0) {
        0
      } else if (reference == 0) {
        Inf
      } else {
        value * log(value / reference)
      }
      total <- total + contribution
    }
  }
  total
}

moment_fugw_oracle_structure <- function(Cx, Cy, sample_plan, feature_plan) {
  total <- 0
  for (i in seq_len(nrow(Cx))) {
    for (j in seq_len(nrow(Cy))) {
      sample_weight <- sample_plan[i, j]
      if (sample_weight == 0) next
      for (k in seq_len(ncol(Cx))) {
        for (l in seq_len(ncol(Cy))) {
          feature_weight <- feature_plan[k, l]
          if (feature_weight == 0) next
          discrepancy <- Cx[i, k] - Cy[j, l]
          total <- total + discrepancy * discrepancy *
            sample_weight * feature_weight
        }
      }
    }
  }
  total
}

moment_fugw_oracle_wrong_cross_structure <- function(
    Cx, Cy, sample_plan, feature_plan) {
  total <- 0
  for (i in seq_len(nrow(Cx))) {
    for (j in seq_len(nrow(Cy))) {
      for (k in seq_len(ncol(Cx))) {
        for (l in seq_len(ncol(Cy))) {
          # Deliberate mutant: +2 Cx Cy instead of -2 Cx Cy.
          discrepancy_mutant <- Cx[i, k] + Cy[j, l]
          total <- total + discrepancy_mutant * discrepancy_mutant *
            sample_plan[i, j] * feature_plan[k, l]
        }
      }
    }
  }
  total
}

moment_fugw_oracle_validate <- function(
    Cx, Cy, M, wx, wy, sample_plan, feature_plan) {
  matrices <- list(
    Cx = Cx, Cy = Cy, M = M,
    sample_plan = sample_plan, feature_plan = feature_plan
  )
  if (any(!vapply(matrices, is.matrix, logical(1))) ||
      any(!vapply(matrices, is.numeric, logical(1))) ||
      any(!vapply(matrices, function(x) all(is.finite(x)), logical(1)))) {
    stop("Oracle matrices must be finite numeric matrices.", call. = FALSE)
  }
  ns <- nrow(Cx)
  nt <- nrow(Cy)
  if (ncol(Cx) != ns || ncol(Cy) != nt ||
      !identical(dim(M), c(ns, nt)) ||
      !identical(dim(sample_plan), c(ns, nt)) ||
      !identical(dim(feature_plan), c(ns, nt))) {
    stop("Oracle costs and plans have incompatible shapes.", call. = FALSE)
  }
  if (any(sample_plan < 0) || any(feature_plan < 0)) {
    stop("Oracle plans must be nonnegative.", call. = FALSE)
  }
  wx <- as.numeric(wx)
  wy <- as.numeric(wy)
  if (length(wx) != ns || length(wy) != nt ||
      any(!is.finite(c(wx, wy))) || any(c(wx, wy) < 0) ||
      sum(wx) <= 0 || sum(wy) <= 0) {
    stop("Oracle reference weights are invalid.", call. = FALSE)
  }
  list(
    Cx = unname(Cx), Cy = unname(Cy), M = unname(M),
    wx = wx, wy = wy,
    sample_plan = unname(sample_plan), feature_plan = unname(feature_plan)
  )
}

moment_fugw_oracle <- function(
    Cx,
    Cy,
    M,
    wx,
    wy,
    sample_plan,
    feature_plan,
    reg_marginals,
    epsilon,
    feature_weight = 1,
    structure_weight = 1) {
  data <- moment_fugw_oracle_validate(
    Cx, Cy, M, wx, wy, sample_plan, feature_plan
  )
  reg_marginals <- as.numeric(reg_marginals)
  if (length(reg_marginals) == 1L) reg_marginals <- rep(reg_marginals, 2L)
  scalars <- c(epsilon, feature_weight, structure_weight, reg_marginals)
  if (length(reg_marginals) != 2L || any(!is.finite(scalars)) ||
      epsilon < 0 || feature_weight < 0 || structure_weight < 0 ||
      any(reg_marginals < 0)) {
    stop("Oracle coefficients must be finite and nonnegative.", call. = FALSE)
  }

  sample_source <- rowSums(data$sample_plan)
  feature_source <- rowSums(data$feature_plan)
  sample_target <- colSums(data$sample_plan)
  feature_target <- colSums(data$feature_plan)

  structure_unweighted <- moment_fugw_oracle_structure(
    data$Cx, data$Cy, data$sample_plan, data$feature_plan
  )
  feature_sample <- sum(data$M * data$sample_plan)
  feature_feature <- sum(data$M * data$feature_plan)
  feature_unweighted <- 0.5 * (feature_sample + feature_feature)
  source_divergence <- moment_fugw_oracle_product_gkl(
    sample_source, feature_source, data$wx, data$wx
  )
  target_divergence <- moment_fugw_oracle_product_gkl(
    sample_target, feature_target, data$wy, data$wy
  )
  reference_plan <- as.numeric(outer(data$wx, data$wy))
  plan_divergence <- moment_fugw_oracle_product_gkl(
    as.numeric(data$sample_plan), as.numeric(data$feature_plan),
    reference_plan, reference_plan
  )
  regularization <- reg_marginals[[1L]] * source_divergence +
    reg_marginals[[2L]] * target_divergence +
    epsilon * plan_divergence

  list(
    structure_unweighted = structure_unweighted,
    feature_sample = feature_sample,
    feature_feature = feature_feature,
    feature_unweighted = feature_unweighted,
    source_marginal_divergence = source_divergence,
    target_marginal_divergence = target_divergence,
    plan_divergence = plan_divergence,
    regularization = regularization,
    structure_weighted = structure_weight * structure_unweighted,
    feature_weighted = feature_weight * feature_unweighted,
    fugw_cost = structure_weight * structure_unweighted +
      feature_weight * feature_unweighted + regularization
  )
}

moment_fugw_oracle_mutants <- function(
    Cx,
    Cy,
    M,
    wx,
    wy,
    sample_plan,
    feature_plan,
    reg_marginals,
    epsilon,
    feature_weight = 1,
    structure_weight = 1) {
  truth <- moment_fugw_oracle(
    Cx, Cy, M, wx, wy, sample_plan, feature_plan,
    reg_marginals, epsilon, feature_weight, structure_weight
  )
  reference_plan <- as.numeric(outer(wx, wy))
  plan_log_only <- moment_fugw_oracle_product_gkl(
    as.numeric(sample_plan), as.numeric(feature_plan),
    reference_plan, reference_plan, mass_terms = FALSE
  )
  wrong_cross <- moment_fugw_oracle_wrong_cross_structure(
    Cx, Cy, sample_plan, feature_plan
  )
  c(
    wrong_structure_cross_sign =
      truth$fugw_cost + structure_weight *
        (wrong_cross - truth$structure_unweighted),
    omitted_source_product_kl =
      truth$fugw_cost - reg_marginals[[1L]] *
        truth$source_marginal_divergence,
    omitted_plan_kl_mass_terms =
      truth$fugw_cost + epsilon * (plan_log_only - truth$plan_divergence),
    normalized_coefficient_convention =
      (1 - feature_weight) * truth$structure_unweighted +
        feature_weight * truth$feature_unweighted + truth$regularization
  )
}
