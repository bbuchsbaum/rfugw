.ti_sparse_cost_edges <- function(cost, n_source = NULL, n_target = NULL) {
  representation <- "edge_list"
  if (is.matrix(cost)) {
    cost <- .validate_finite_matrix(cost, "cost")
    n_source <- nrow(cost)
    n_target <- ncol(cost)
    edges <- data.frame(
      source = rep(seq_len(n_source), each = n_target),
      target = rep.int(seq_len(n_target), times = n_source),
      cost = as.numeric(t(cost))
    )
    representation <- "dense_complete"
  } else if (inherits(cost, "sparseMatrix")) {
    if (!requireNamespace("Matrix", quietly = TRUE)) {
      stop("Package `Matrix` is required for a sparse Matrix cost.", call. = FALSE)
    }
    n_source <- nrow(cost)
    n_target <- ncol(cost)
    entries <- Matrix::summary(methods::as(cost, "CsparseMatrix"))
    edge_cost <- if ("x" %in% names(entries)) {
      entries$x
    } else {
      rep(1, nrow(entries))
    }
    edges <- data.frame(
      source = entries$i, target = entries$j, cost = edge_cost
    )
    representation <- "matrix_structural_support"
  } else if (is.data.frame(cost) || is.list(cost)) {
    if (all(c("row_ptr", "col_idx", "cost", "n_rows", "n_cols") %in%
            names(cost))) {
      n_source <- cost$n_rows
      n_target <- cost$n_cols
      row_ptr <- as.numeric(cost$row_ptr)
      offset <- if (length(row_ptr) && identical(row_ptr[[1L]], 0)) 0 else 1
      if (length(row_ptr) != n_source + 1L ||
          any(!is.finite(row_ptr)) || any(row_ptr != floor(row_ptr)) ||
          any(diff(row_ptr) < 0) || row_ptr[[1L]] - offset != 0 ||
          tail(row_ptr, 1L) - offset != length(cost$col_idx)) {
        stop("Sparse `row_ptr` must consistently span every CSR edge.", call. = FALSE)
      }
      edges <- data.frame(
        source = rep.int(seq_len(n_source), diff(row_ptr)),
        target = cost$col_idx,
        cost = cost$cost
      )
      representation <- "csr"
    } else {
      fields <- if (all(c("source", "target", "cost") %in% names(cost))) {
        c("source", "target", "cost")
      } else if (all(c("i", "j", "x") %in% names(cost))) {
        c("i", "j", "x")
      } else {
        stop(
          paste0(
            "Sparse costs need source/target/cost, i/j/x, or canonical ",
            "CSR fields."
          ),
          call. = FALSE
        )
      }
      edges <- data.frame(
        source = cost[[fields[[1L]]]],
        target = cost[[fields[[2L]]]],
        cost = cost[[fields[[3L]]]]
      )
      n_source <- n_source %||% cost$n_source %||% cost$n_rows
      n_target <- n_target %||% cost$n_target %||% cost$n_cols
    }
  } else {
    stop("`cost` must be a dense matrix, sparse Matrix, edge list, or CSR list.",
         call. = FALSE)
  }

  shape <- c(n_source, n_target)
  if (length(shape) != 2L || any(!is.finite(shape)) || any(shape < 1) ||
      any(shape != floor(shape))) {
    stop(
      "Sparse edge costs require positive integer `n_source` and `n_target`.",
      call. = FALSE
    )
  }
  n_source <- as.integer(n_source)
  n_target <- as.integer(n_target)
  if (!is.numeric(edges$source) || !is.numeric(edges$target) ||
      !is.numeric(edges$cost) ||
      length(edges$source) != length(edges$target) ||
      length(edges$source) != length(edges$cost)) {
    stop("Sparse cost edge fields must be equal-length numeric vectors.",
         call. = FALSE)
  }
  if (any(!is.finite(edges$source)) ||
      any(edges$source != floor(edges$source)) ||
      any(edges$source < 1 | edges$source > n_source)) {
    stop("Sparse cost source indices are out of range.", call. = FALSE)
  }
  if (any(!is.finite(edges$target)) ||
      any(edges$target != floor(edges$target)) ||
      any(edges$target < 1 | edges$target > n_target)) {
    stop("Sparse cost target indices are out of range.", call. = FALSE)
  }
  if (any(!is.finite(edges$cost))) {
    stop("Sparse edge costs must be finite.", call. = FALSE)
  }
  edges$source <- as.integer(edges$source)
  edges$target <- as.integer(edges$target)
  edges$cost <- as.numeric(edges$cost)
  ordering <- order(edges$source, edges$target)
  edges <- edges[ordering, , drop = FALSE]
  duplicate <- duplicated(edges[c("source", "target")])
  if (any(duplicate)) {
    stop("Sparse cost support must not contain duplicate edges.", call. = FALSE)
  }

  row_count <- tabulate(edges$source, nbins = n_source)
  col_order <- order(edges$target, edges$source)
  col_edges <- edges[col_order, , drop = FALSE]
  col_count <- tabulate(col_edges$target, nbins = n_target)
  list(
    edges = edges,
    col_edges = col_edges,
    col_order = col_order,
    row_ptr = as.integer(c(0, cumsum(row_count))),
    col_idx = edges$target,
    row_cost = edges$cost,
    col_ptr = as.integer(c(0, cumsum(col_count))),
    row_idx = col_edges$source,
    col_cost = col_edges$cost,
    n_source = n_source,
    n_target = n_target,
    representation = representation,
    complete = nrow(edges) == n_source * n_target
  )
}

.ti_active_support_certificate <- function(prepared, source, target, tolerance) {
  edges <- prepared$edges
  active_edge <- source[edges$source] > 0 & target[edges$target] > 0
  active_source_degree <- tabulate(
    edges$source[active_edge], nbins = prepared$n_source
  )
  active_target_degree <- tabulate(
    edges$target[active_edge], nbins = prepared$n_target
  )
  uncovered_source <- which(source > 0 & active_source_degree == 0L)
  uncovered_target <- which(target > 0 & active_target_degree == 0L)
  finite_potential_support <- !length(uncovered_source) && !length(uncovered_target)
  flow <- cpp_bipartite_transport_max_flow(
    prepared$n_source,
    prepared$n_target,
    edges$source,
    edges$target,
    source,
    target,
    tolerance
  )
  list(
    uot_primal_feasible = TRUE,
    finite_potential_support = finite_potential_support,
    uncovered_source = uncovered_source,
    uncovered_target = uncovered_target,
    active_support_size = sum(active_edge),
    balanced_marginal_feasible = isTRUE(flow$feasible),
    balanced_max_flow = flow$max_flow,
    balanced_flow_deficit = flow$deficit,
    balanced_equal_total_mass = isTRUE(flow$equal_total),
    balanced_check_method = flow$method,
    interpretation = paste0(
      "KL-UOT remains primal-feasible through the zero coupling; finite TI ",
      "potentials additionally require active support coverage. Exact balanced ",
      "marginal feasibility is reported separately by global max flow."
    )
  )
}

.ti_logsumexp_slice <- function(values) {
  if (!length(values)) return(-Inf)
  maximum <- max(values)
  if (!is.finite(maximum)) return(maximum)
  maximum + log(sum(exp(values - maximum)))
}

.ti_group_logsumexp <- function(values, pointer) {
  vapply(seq_len(length(pointer) - 1L), function(index) {
    start <- pointer[[index]] + 1L
    end <- pointer[[index + 1L]]
    if (start > end) -Inf else .ti_logsumexp_slice(values[start:end])
  }, numeric(1))
}

.ti_softmin_vector <- function(log_weight, value, temperature) {
  -temperature * .ti_logsumexp_slice(log_weight - value / temperature)
}

.ti_fixed_point_residual <- function(
    prepared, source, target, epsilon, rho, source_bar, target_bar) {
  log_source <- ifelse(source > 0, log(source), -Inf)
  log_target <- ifelse(target > 0, log(target), -Inf)
  denominator <- epsilon + sum(rho)
  xi12 <- epsilon * rho[[2L]] / (rho[[1L]] * denominator)
  xi21 <- epsilon * rho[[1L]] / (rho[[2L]] * denominator)
  k11 <- epsilon / (epsilon + rho[[1L]]) * rho[[1L]] / sum(rho)
  k22 <- epsilon / (epsilon + rho[[2L]]) * rho[[2L]] / sum(rho)
  a1 <- rho[[1L]] / (rho[[1L]] + epsilon)
  a2 <- rho[[2L]] / (rho[[2L]] + epsilon)
  target_scalar <- .ti_softmin_vector(log_target, target_bar, rho[[2L]])
  source_tmp <- vapply(seq_len(prepared$n_source), function(i) {
    start <- prepared$row_ptr[[i]] + 1L
    end <- prepared$row_ptr[[i + 1L]]
    if (start > end) return(if (source[[i]] == 0) 0 else Inf)
    index <- start:end
    active <- target[prepared$col_idx[index]] > 0
    if (!any(active)) return(if (source[[i]] == 0) 0 else Inf)
    index <- index[active]
    soft <- -epsilon * .ti_logsumexp_slice(
      log_target[prepared$col_idx[index]] +
        (target_bar[prepared$col_idx[index]] - prepared$row_cost[index]) /
          epsilon
    )
    a1 * soft - k11 * target_scalar
  }, numeric(1))
  source_tmp[source == 0] <- 0
  source_next <- source_tmp + xi12 *
    .ti_softmin_vector(log_source, source_tmp, rho[[1L]])
  source_next[source == 0] <- 0
  source_scalar <- .ti_softmin_vector(log_source, source_next, rho[[1L]])
  target_tmp <- vapply(seq_len(prepared$n_target), function(j) {
    start <- prepared$col_ptr[[j]] + 1L
    end <- prepared$col_ptr[[j + 1L]]
    if (start > end) return(if (target[[j]] == 0) 0 else Inf)
    index <- start:end
    active <- source[prepared$row_idx[index]] > 0
    if (!any(active)) return(if (target[[j]] == 0) 0 else Inf)
    index <- index[active]
    soft <- -epsilon * .ti_logsumexp_slice(
      log_source[prepared$row_idx[index]] +
        (source_next[prepared$row_idx[index]] - prepared$col_cost[index]) /
          epsilon
    )
    a2 * soft - k22 * source_scalar
  }, numeric(1))
  target_tmp[target == 0] <- 0
  target_next <- target_tmp + xi21 *
    .ti_softmin_vector(log_target, target_tmp, rho[[2L]])
  target_next[target == 0] <- 0
  max(abs(c(source_next - source_bar, target_next - target_bar)))
}

.ti_edge_log_weights <- function(
    prepared, source, target, source_bar, target_bar, epsilon) {
  edges <- prepared$edges
  log_source <- ifelse(source > 0, log(source), -Inf)
  log_target <- ifelse(target > 0, log(target), -Inf)
  log_source[edges$source] + log_target[edges$target] +
    (source_bar[edges$source] + target_bar[edges$target] - edges$cost) /
      epsilon
}

.ti_edge_weights <- function(
    prepared, source, target, source_bar, target_bar, epsilon) {
  exp(.ti_edge_log_weights(
    prepared, source, target, source_bar, target_bar, epsilon
  ))
}

.ti_operator_apply <- function(
    prepared, weight, values, orientation = c("forward", "adjoint")) {
  orientation <- match.arg(orientation)
  edges <- prepared$edges
  if (identical(orientation, "forward")) {
    .transport_edge_apply(
      transform(edges, weight = weight), prepared$n_source, values,
      index = "target", group = "source"
    )
  } else {
    .transport_edge_apply(
      transform(edges, weight = weight), prepared$n_target, values,
      index = "source", group = "target"
    )
  }
}

.ti_plan_from_potentials <- function(
    prepared, source, target, source_bar, target_bar, epsilon,
    representation = c("operator", "sparse", "dense")) {
  representation <- match.arg(representation)
  weight_function <- function() {
    .ti_edge_weights(
      prepared, source, target, source_bar, target_bar, epsilon
    )
  }
  weights <- weight_function()
  source_mass <- .weighted_transport_tabulate(
    prepared$edges$source, weights, prepared$n_source
  )
  target_mass <- .weighted_transport_tabulate(
    prepared$edges$target, weights, prepared$n_target
  )
  edge_plan <- function() {
    as_transport_plan(
      data.frame(
        source = prepared$edges$source,
        target = prepared$edges$target,
        weight = weight_function()
      ),
      prepared$n_source,
      prepared$n_target,
      duplicates = "error"
    )
  }
  if (identical(representation, "sparse")) return(edge_plan())
  if (identical(representation, "dense")) {
    return(transport_plan_materialize(edge_plan()))
  }
  transport_operator(
    prepared$n_source,
    prepared$n_target,
    apply = function(x) {
      .ti_operator_apply(prepared, weight_function(), x, "forward")
    },
    adjoint = function(x) {
      .ti_operator_apply(prepared, weight_function(), x, "adjoint")
    },
    source_mass = source_mass,
    target_mass = target_mass,
    materialize = edge_plan,
    metadata = list(
      formulation = "translation_invariant_kl_uot",
      support_size = nrow(prepared$edges),
      cost_representation = prepared$representation,
      coupling_stored = FALSE
    )
  )
}

.ti_certificate <- function(
    prepared, source, target, epsilon, rho, native, tolerance) {
  log_weight <- .ti_edge_log_weights(
    prepared, source, target, native$source_bar, native$target_bar, epsilon
  )
  weight <- exp(log_weight)
  finite_plan <- all(is.finite(weight)) && all(weight >= 0)
  source_marginal <- .weighted_transport_tabulate(
    prepared$edges$source, weight, prepared$n_source
  )
  target_marginal <- .weighted_transport_tabulate(
    prepared$edges$target, weight, prepared$n_target
  )
  log_source_marginal <- .ti_group_logsumexp(log_weight, prepared$row_ptr)
  log_target_marginal <- .ti_group_logsumexp(
    log_weight[prepared$col_order], prepared$col_ptr
  )
  active_source <- source > 0
  active_target <- target > 0
  source_kkt <- max(abs(
    log_source_marginal[active_source] -
      (log(source[active_source]) - native$source_potential[active_source] /
         rho[[1L]])
  ))
  target_kkt <- max(abs(
    log_target_marginal[active_target] -
      (log(target[active_target]) - native$target_potential[active_target] /
         rho[[2L]])
  ))
  kkt_residual <- max(source_kkt, target_kkt)
  fixed_point_residual <- .ti_fixed_point_residual(
    prepared, source, target, epsilon, rho,
    native$source_bar, native$target_bar
  )
  source_kl <- .generalized_kl_vector(source_marginal, source)
  target_kl <- .generalized_kl_vector(target_marginal, target)
  positive_weight <- weight > 0
  product_kl <- sum(weight[positive_weight] *
    (log_weight[positive_weight] -
       log(source[prepared$edges$source[positive_weight]]) -
       log(target[prepared$edges$target[positive_weight]]))) -
    sum(weight) + sum(source) * sum(target)
  transport <- sum(weight * prepared$edges$cost)
  primal <- transport + rho[[1L]] * source_kl +
    rho[[2L]] * target_kl + epsilon * product_kl
  dual_source_mass <- sum(source * exp(-native$source_potential / rho[[1L]]))
  dual_target_mass <- sum(target * exp(-native$target_potential / rho[[2L]]))
  dual <- rho[[1L]] * (sum(source) - dual_source_mass) +
    rho[[2L]] * (sum(target) - dual_target_mass) +
    epsilon * (sum(source) * sum(target) - sum(weight))
  gap <- primal - dual
  objective_tolerance <- max(1e-8, 100 * tolerance) * max(1, abs(primal))
  kkt_tolerance <- max(1e-7, 100 * tolerance)
  fixed_point_tolerance <- max(1e-8, 20 * tolerance)
  gauge_shift <- 1.23456789
  shifted_log_weight <- .ti_edge_log_weights(
    prepared, source, target,
    native$source_bar + gauge_shift,
    native$target_bar - gauge_shift,
    epsilon
  )
  finite_log_weight <- is.finite(log_weight)
  gauge_residual <- if (any(finite_log_weight)) {
    max(abs(shifted_log_weight[finite_log_weight] -
      log_weight[finite_log_weight]))
  } else {
    0
  }
  list(
    weight = weight,
    source_marginal = source_marginal,
    target_marginal = target_marginal,
    finite_plan = finite_plan,
    fixed_point_residual = fixed_point_residual,
    fixed_point_tolerance = fixed_point_tolerance,
    fixed_point_consistent = is.finite(fixed_point_residual) &&
      fixed_point_residual <= fixed_point_tolerance,
    source_kkt_residual = source_kkt,
    target_kkt_residual = target_kkt,
    kkt_residual = kkt_residual,
    kkt_tolerance = kkt_tolerance,
    kkt_consistent = is.finite(kkt_residual) && kkt_residual <= kkt_tolerance,
    transport = transport,
    source_kl = source_kl,
    target_kl = target_kl,
    product_kl = product_kl,
    primal = primal,
    dual = dual,
    gap = gap,
    objective_tolerance = objective_tolerance,
    objective_consistent = is.finite(gap) && gap >= -objective_tolerance &&
      abs(gap) <= objective_tolerance,
    transported_mass = sum(weight),
    dual_source_mass = dual_source_mass,
    dual_target_mass = dual_target_mass,
    gauge_residual = gauge_residual,
    gauge_invariant = is.finite(gauge_residual) && gauge_residual <= 1e-12
  )
}

#' Translation-Invariant Sparse KL-UOT
#'
#' Solves finite-measure entropic unbalanced optimal transport with generalized
#' KL penalties and translation-invariant Sinkhorn updates. Dense costs use
#' complete support. Sparse costs may be a `Matrix::sparseMatrix` (structural
#' nonzeros are the support), a 1-based `source,target,cost` edge list, or the
#' CSR list used by manifoldalign. Edge-list shapes are supplied with
#' `n_source` and `n_target` or embedded in the list.
#'
#' The optimized objective is
#' \deqn{\langle C, \pi\rangle + \rho_1 KL(\pi_1|p) +
#' \rho_2 KL(\pi_2|q) + \epsilon KL(\pi|p\otimes q).}
#' KL terms use finite-measure generalized KL, with the product-reference
#' constant retained even outside sparse support. Potentials are returned in a
#' canonical translated gauge. The plan is implicit by default, so apply and
#' adjoint operations do not require dense materialization.
#'
#' @param cost Dense matrix, sparse Matrix, edge list, CSR sparse-cost list, or
#'   an affine-bilinear `rfugw_cost_operator`. Cost operators use complete
#'   implicit support and native blocked reductions.
#' @param p,q Finite nonnegative source and target measures. Defaults are
#'   probability measures; positive total mass is required.
#' @param epsilon Positive entropic regularization.
#' @param rho One or two positive marginal KL penalties.
#' @param max_iter Maximum TI-Sinkhorn iterations.
#' @param tol Positive iterate tolerance. Independent KKT, fixed-point, and
#'   primal-dual tolerances are derived and reported.
#' @param n_source,n_target Sparse edge-list shape; ignored for matrix costs.
#' @param plan Returned plan representation: coupling-free `"operator"`,
#'   explicit `"sparse"` edge plan, or explicitly allocated `"dense"` matrix.
#' @param init_potentials Optional list with `source_bar` and `target_bar` for a
#'   factorized-cost warm start. Explicit edge and matrix costs do not yet
#'   accept this state.
#' @param block_size Positive tile size for factorized complete-support costs.
#' @return An `rfugw_result` with potentials, sparse-support diagnostics,
#'   objective components, KKT and primal-dual certificates, mass diagnostics,
#'   plan/operator representation, and runtime provenance.
#' @examples
#' cost <- matrix(c(0, 1, 1, 0), 2)
#' fit <- ot_sinkhorn_unbalanced_ti(
#'   cost, p = c(2, 1), q = c(1, 3), epsilon = 0.2, rho = c(2, 3)
#' )
#' transport_plan_adjoint(rfugw_plan(fit), c(10, 20))
#' @export
ot_sinkhorn_unbalanced_ti <- function(
    cost,
    p = NULL,
    q = NULL,
    epsilon = 0.05,
    rho = 10,
    max_iter = 2000L,
    tol = 1e-8,
    n_source = NULL,
    n_target = NULL,
    plan = c("operator", "sparse", "dense"),
    init_potentials = NULL,
    block_size = 256L) {
  plan <- match.arg(plan)
  epsilon <- .validate_positive_scalar(epsilon, "epsilon")
  max_iter <- .validate_count(max_iter, "max_iter")
  tol <- .validate_positive_scalar(tol, "tol")
  if (length(rho) == 1L) rho <- rep(rho, 2L)
  if (!is.numeric(rho) || length(rho) != 2L ||
      any(!is.finite(rho)) || any(rho <= 0)) {
    stop("`rho` must be one or two finite positive numbers.", call. = FALSE)
  }

  if (inherits(cost, "rfugw_cost_operator")) {
    shape <- cost_shape(cost)
    if (!is.null(n_source)) {
      n_source <- .validate_cost_shape_scalar(n_source, "n_source")
      if (!identical(n_source, shape[[1L]])) {
        stop("`n_source` conflicts with the cost-operator shape.", call. = FALSE)
      }
    }
    if (!is.null(n_target)) {
      n_target <- .validate_cost_shape_scalar(n_target, "n_target")
      if (!identical(n_target, shape[[2L]])) {
        stop("`n_target` conflicts with the cost-operator shape.", call. = FALSE)
      }
    }
    return(.ot_sinkhorn_unbalanced_ti_factorized(
      cost = cost,
      p = p,
      q = q,
      epsilon = epsilon,
      rho = rho,
      max_iter = max_iter,
      tol = tol,
      plan = plan,
      init_potentials = init_potentials,
      block_size = block_size
    ))
  }
  if (!is.null(init_potentials)) {
    stop(
      "`init_potentials` is currently supported only for factorized costs.",
      call. = FALSE
    )
  }

  setup_started <- proc.time()[["elapsed"]]
  prepared <- .ti_sparse_cost_edges(cost, n_source, n_target)
  p <- .validate_finite_measure(
    p, prepared$n_source, "p", rep(1 / prepared$n_source, prepared$n_source)
  )
  q <- .validate_finite_measure(
    q, prepared$n_target, "q", rep(1 / prepared$n_target, prepared$n_target)
  )
  if (sum(p) <= 0 || sum(q) <= 0) {
    stop(
      "TI-UOT requires positive total mass on both measures; use `ot_sinkhorn_unbalanced()` for a zero-measure closed form.",
      call. = FALSE
    )
  }
  flow_tolerance <- max(1e-12, tol * max(1, sum(p), sum(q)))
  support <- .ti_active_support_certificate(prepared, p, q, flow_tolerance)
  setup_seconds <- proc.time()[["elapsed"]] - setup_started
  if (!isTRUE(support$finite_potential_support)) {
    stop(
      paste0(
        "Sparse support cannot yield finite TI potentials: active source or ",
        "target nodes are uncovered."
      ),
      call. = FALSE
    )
  }

  native <- cpp_ot_sinkhorn_unbalanced_ti_sparse(
    prepared$row_ptr,
    prepared$col_idx,
    prepared$row_cost,
    prepared$col_ptr,
    prepared$row_idx,
    prepared$col_cost,
    prepared$n_source,
    prepared$n_target,
    p,
    q,
    epsilon,
    rho[[1L]],
    rho[[2L]],
    max_iter,
    tol
  )
  certificate_started <- proc.time()[["elapsed"]]
  certificate <- .ti_certificate(
    prepared, p, q, epsilon, rho, native, tol
  )
  result_plan <- .ti_plan_from_potentials(
    prepared, p, q, native$source_bar, native$target_bar, epsilon, plan
  )
  certificate_seconds <- proc.time()[["elapsed"]] - certificate_started

  inner_residual <- max(
    certificate$fixed_point_residual,
    certificate$kkt_residual
  )
  inner_tolerance <- max(
    certificate$fixed_point_tolerance,
    certificate$kkt_tolerance
  )
  inner_converged <- isTRUE(certificate$fixed_point_consistent) &&
    isTRUE(certificate$kkt_consistent) &&
    isTRUE(certificate$objective_consistent) &&
    isTRUE(certificate$gauge_invariant)
  out <- list(
    plan = result_plan,
    formulation = "ot_sinkhorn_unbalanced_ti",
    backend = native$backend,
    objective = certificate$primal,
    regularized_objective = certificate$primal,
    dual_objective = certificate$dual,
    primal_dual_gap = certificate$gap,
    transport_cost = certificate$transport,
    source_marginal_kl = certificate$source_kl,
    target_marginal_kl = certificate$target_kl,
    plan_product_kl = certificate$product_kl,
    transported_mass = certificate$transported_mass,
    source_marginal = certificate$source_marginal,
    target_marginal = certificate$target_marginal,
    source_measure_effective = p,
    target_measure_effective = q,
    original_source_mass = sum(p),
    original_target_mass = sum(q),
    effective_source_mass = sum(p),
    effective_target_mass = sum(q),
    measure_normalization = "none",
    regularization = epsilon,
    rho = rho,
    source_bar = as.numeric(native$source_bar),
    target_bar = as.numeric(native$target_bar),
    fbar = as.numeric(native$source_bar),
    gbar = as.numeric(native$target_bar),
    translation = as.numeric(native$translation),
    source_potential = as.numeric(native$source_potential),
    target_potential = as.numeric(native$target_potential),
    f = as.numeric(native$source_potential),
    g = as.numeric(native$target_potential),
    fixed_point_residual = certificate$fixed_point_residual,
    fixed_point_tolerance = certificate$fixed_point_tolerance,
    fixed_point_consistent = certificate$fixed_point_consistent,
    source_kkt_residual = certificate$source_kkt_residual,
    target_kkt_residual = certificate$target_kkt_residual,
    kkt_residual = certificate$kkt_residual,
    kkt_tolerance = certificate$kkt_tolerance,
    kkt_consistent = certificate$kkt_consistent,
    gauge_residual = certificate$gauge_residual,
    gauge_invariant = certificate$gauge_invariant,
    support_certificate = support,
    cost_representation = prepared$representation,
    support_size = nrow(prepared$edges),
    support_density = nrow(prepared$edges) /
      (prepared$n_source * prepared$n_target),
    plan_representation = transport_plan_representation(result_plan),
    iterations = native$iterations,
    error = max(native$residual, inner_residual),
    native_residual = native$residual,
    inner_residual = inner_residual,
    max_inner_residual = inner_residual,
    inner_iterations = native$iterations,
    inner_converged = inner_converged,
    inner_status = if (inner_converged) "converged" else "certificate_failure"
  )
  ans <- .attach_solver_diagnostics(
    out,
    residual = max(native$residual, inner_residual),
    converged = isTRUE(native$converged),
    iterations = native$iterations,
    max_iter = max_iter,
    plan = result_plan,
    inner_residual = inner_residual,
    max_inner_residual = inner_residual,
    inner_iterations = native$iterations,
    inner_converged = inner_converged,
    inner_status = out$inner_status,
    feasibility = "unbalanced",
    feasibility_tol = inner_tolerance,
    objective_recomputed = certificate$primal,
    objective_tolerance = certificate$objective_tolerance,
    objective_components = list(
      transport_cost = certificate$transport,
      source_marginal_kl = certificate$source_kl,
      target_marginal_kl = certificate$target_kl,
      plan_product_kl = certificate$product_kl
    ),
    numerical_ok = isTRUE(native$numerical_ok) && certificate$finite_plan
  )
  if (!isTRUE(native$numerical_ok) || !certificate$finite_plan) {
    ans$status <- "numerical_failure"
    ans$converged <- FALSE
  } else if (!isTRUE(support$finite_potential_support)) {
    ans$status <- "infeasible"
    ans$converged <- FALSE
  } else if (!isTRUE(certificate$objective_consistent)) {
    ans$status <- "objective_mismatch"
    ans$converged <- FALSE
  }
  ans$termination_reason <- .termination_reason_from_result(ans, max_iter)
  ans$uot_certificate <- certificate[names(certificate) != "weight"]
  ans$runtime_provenance$timing <- list(
    setup_seconds = unname(setup_seconds),
    solve_seconds = unname(native$solve_seconds),
    certificate_seconds = unname(certificate_seconds)
  )
  ans$runtime_provenance$support <- list(
    size = nrow(prepared$edges),
    density = ans$support_density,
    representation = prepared$representation,
    coupling_stored = !identical(plan, "operator")
  )
  ans
}
