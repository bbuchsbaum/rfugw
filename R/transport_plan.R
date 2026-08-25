.is_transport_plan <- function(x) inherits(x, "rfugw_transport_plan")

.canonical_transport_edges <- function(
    source, target, weight, n_source, n_target,
    duplicates = c("sum", "error")) {
  duplicates <- match.arg(duplicates)
  if (!is.numeric(source) || !is.numeric(target) || !is.numeric(weight)) {
    stop("Edge indices and weights must be numeric.", call. = FALSE)
  }
  if (length(source) != length(target) || length(source) != length(weight)) {
    stop("Edge source, target, and weight vectors must have equal length.", call. = FALSE)
  }
  if (any(!is.finite(source)) || any(source != floor(source)) ||
      any(source < 1) || any(source > n_source)) {
    stop("Edge source indices must be finite 1-based integers in range.", call. = FALSE)
  }
  if (any(!is.finite(target)) || any(target != floor(target)) ||
      any(target < 1) || any(target > n_target)) {
    stop("Edge target indices must be finite 1-based integers in range.", call. = FALSE)
  }
  if (any(!is.finite(weight))) {
    stop("Edge weights must be finite.", call. = FALSE)
  }
  if (any(weight < 0)) {
    stop("Edge weights must be nonnegative.", call. = FALSE)
  }
  if (!length(weight)) {
    return(data.frame(
      source = integer(), target = integer(), weight = numeric()
    ))
  }
  ordering <- order(source, target)
  source <- as.integer(source[ordering])
  target <- as.integer(target[ordering])
  weight <- as.numeric(weight[ordering])
  duplicated_edge <- c(FALSE, source[-1L] == source[-length(source)] &
    target[-1L] == target[-length(target)])
  if (any(duplicated_edge) && identical(duplicates, "error")) {
    stop("Duplicate transport edges are not allowed.", call. = FALSE)
  }
  if (any(duplicated_edge)) {
    group <- cumsum(!duplicated_edge)
    first <- !duplicated_edge
    weight <- as.numeric(rowsum(weight, group, reorder = FALSE))
    source <- source[first]
    target <- target[first]
  }
  keep <- weight > 0
  data.frame(
    source = source[keep], target = target[keep], weight = weight[keep]
  )
}

.new_transport_plan <- function(
    data, n_source, n_target, representation, metadata = list()) {
  structure(
    list(
      data = data,
      n_source = as.integer(n_source),
      n_target = as.integer(n_target),
      representation = representation,
      materialized = identical(representation, "dense_materialized"),
      sparse = representation %in% c("sparse_csc", "edge_list"),
      implicit = identical(representation, "implicit_operator"),
      pruned = isTRUE(metadata$pruned),
      lost_mass = metadata$lost_mass %||% 0,
      metadata = metadata
    ),
    class = "rfugw_transport_plan"
  )
}

#' Construct a Transport Plan Representation
#'
#' Wraps a dense numeric matrix, a `Matrix::sparseMatrix`, or a 1-based edge
#' list without changing its mathematical coupling. Edge lists are sorted by
#' source then target. Duplicate edges are summed by default (or rejected), and
#' resulting zero-weight entries are removed. Missing rows or columns are
#' valid and have explicit zero transported mass.
#'
#' @param x Dense matrix, `Matrix::sparseMatrix`, edge data frame/list, existing
#'   transport plan, or result containing a plan. Edge columns may be named
#'   `source,target,weight` or `i,j,x`.
#' @param n_source,n_target Required positive shape for an edge list; inferred
#'   from matrix inputs.
#' @param duplicates Edge-list duplicate policy: canonicalize by summing or
#'   reject.
#' @param ... Ignored by the matrix-conversion and print methods.
#' @return A serializable `rfugw_transport_plan`.
#' @examples
#' edges <- data.frame(source = c(2, 1, 1), target = c(1, 2, 2),
#'                     weight = c(0.4, 0.2, 0.3))
#' plan <- as_transport_plan(edges, 2, 2)
#' transport_plan_mass(plan, "source")
#' @export
as_transport_plan <- function(
    x, n_source = NULL, n_target = NULL,
    duplicates = c("sum", "error")) {
  duplicates <- match.arg(duplicates)
  if (.is_transport_plan(x)) return(x)
  if (inherits(x, "rfugw_result") ||
      (is.list(x) && !is.null(x$plan) && !is.data.frame(x))) {
    return(as_transport_plan(
      rfugw_plan(x), n_source, n_target, duplicates = duplicates
    ))
  }
  if (is.matrix(x)) {
    if (!is.numeric(x) || any(!is.finite(x)) || any(x < 0)) {
      stop("Dense transport plans must be finite and nonnegative.", call. = FALSE)
    }
    if (nrow(x) < 1L || ncol(x) < 1L) {
      stop("Transport-plan shape must be nonempty.", call. = FALSE)
    }
    return(.new_transport_plan(
      unname(x), nrow(x), ncol(x), "dense_materialized"
    ))
  }
  if (inherits(x, "sparseMatrix")) {
    if (!requireNamespace("Matrix", quietly = TRUE)) {
      stop("Package `Matrix` is required for sparse plans.", call. = FALSE)
    }
    if (nrow(x) < 1L || ncol(x) < 1L) {
      stop("Transport-plan shape must be nonempty.", call. = FALSE)
    }
    x <- methods::as(x, "CsparseMatrix")
    entries <- Matrix::summary(x)
    values <- if ("x" %in% names(entries)) entries$x else rep(1, nrow(entries))
    if (any(!is.finite(values)) || any(values < 0)) {
      stop("Sparse transport plans must be finite and nonnegative.", call. = FALSE)
    }
    x <- Matrix::drop0(x)
    return(.new_transport_plan(
      x, nrow(x), ncol(x), "sparse_csc",
      metadata = list(storage_class = class(x)[[1L]])
    ))
  }
  if (is.data.frame(x) || is.list(x)) {
    if (is.null(n_source) || is.null(n_target) ||
        length(n_source) != 1L || length(n_target) != 1L ||
        !is.finite(n_source) || !is.finite(n_target) ||
        n_source < 1 || n_target < 1 ||
        n_source != floor(n_source) || n_target != floor(n_target)) {
      stop(
        "Edge-list plans require positive integer `n_source` and `n_target`.",
        call. = FALSE
      )
    }
    names_x <- names(x)
    fields <- if (all(c("source", "target", "weight") %in% names_x)) {
      c("source", "target", "weight")
    } else if (all(c("i", "j", "x") %in% names_x)) {
      c("i", "j", "x")
    } else {
      stop(
        "Edge lists need source/target/weight or i/j/x fields.",
        call. = FALSE
      )
    }
    edges <- .canonical_transport_edges(
      x[[fields[[1L]]]], x[[fields[[2L]]]], x[[fields[[3L]]]],
      as.integer(n_source), as.integer(n_target), duplicates
    )
    return(.new_transport_plan(
      edges, n_source, n_target, "edge_list",
      metadata = list(
        canonical_order = "source_then_target",
        duplicate_policy = duplicates,
        zero_edges_removed = TRUE
      )
    ))
  }
  stop(
    "`x` must be a dense matrix, sparseMatrix, edge list, or transport plan.",
    call. = FALSE
  )
}

#' Construct an Implicit Transport Operator
#'
#' Creates a coupling-like linear operator when entries need not be stored.
#' Both forward and adjoint callbacks are mandatory; supplied row and column
#' masses define barycentric abstention. `adjoint` is the mathematical
#' transpose action, not an inverse or reverse conditional. A materializer is
#' optional and is called only by explicit [transport_plan_materialize()].
#'
#' @param n_source,n_target Positive operator shape.
#' @param apply,adjoint Functions accepting a numeric destination/source matrix
#'   and returning the forward/adjoint product.
#' @param source_mass,target_mass Finite nonnegative row and column masses with
#'   equal totals.
#' @param materialize Optional zero-argument function returning a dense or
#'   sparse nonnegative plan of the declared shape.
#' @param metadata Serializable provenance list.
#' @return An implicit `rfugw_transport_plan`.
#' @export
transport_operator <- function(
    n_source, n_target, apply, adjoint, source_mass, target_mass,
    materialize = NULL, metadata = list()) {
  if (length(n_source) != 1L || length(n_target) != 1L ||
      !is.finite(n_source) || !is.finite(n_target) ||
      n_source < 1 || n_target < 1 ||
      n_source != floor(n_source) || n_target != floor(n_target)) {
    stop("Operator shape must contain positive integers.", call. = FALSE)
  }
  if (!is.function(apply) || !is.function(adjoint)) {
    stop("`apply` and `adjoint` must be functions.", call. = FALSE)
  }
  if (!is.null(materialize) && !is.function(materialize)) {
    stop("`materialize` must be NULL or a function.", call. = FALSE)
  }
  if (!is.numeric(source_mass) || length(source_mass) != n_source ||
      any(!is.finite(source_mass)) || any(source_mass < 0)) {
    stop("`source_mass` must be finite, nonnegative, and match shape.", call. = FALSE)
  }
  if (!is.numeric(target_mass) || length(target_mass) != n_target ||
      any(!is.finite(target_mass)) || any(target_mass < 0)) {
    stop("`target_mass` must be finite, nonnegative, and match shape.", call. = FALSE)
  }
  mass_tolerance <- 1e-10 * max(1, sum(source_mass), sum(target_mass))
  if (abs(sum(source_mass) - sum(target_mass)) > mass_tolerance) {
    stop("Source and target operator masses must have equal totals.", call. = FALSE)
  }
  data <- list(
    apply = apply,
    adjoint = adjoint,
    source_mass = as.numeric(source_mass),
    target_mass = as.numeric(target_mass),
    materialize = materialize
  )
  .new_transport_plan(
    data, n_source, n_target, "implicit_operator", metadata = metadata
  )
}

#' Inspect Transport Plan Shape
#'
#' @param plan Transport plan, operator, result, dense matrix, or sparse matrix.
#' @return Integer vector `c(source, target)`.
#' @export
transport_plan_shape <- function(plan) {
  plan <- as_transport_plan(plan)
  c(source = plan$n_source, target = plan$n_target)
}

.weighted_transport_tabulate <- function(index, weight, nbins) {
  out <- numeric(nbins)
  if (length(index)) {
    summed <- rowsum(weight, index, reorder = FALSE)
    out[as.integer(rownames(summed))] <- as.numeric(summed)
  }
  out
}

#' Inspect Transported Mass
#'
#' @inheritParams transport_plan_shape
#' @param margin `"total"`, `"source"` (row masses), or `"target"` (column
#'   masses).
#' @return Numeric scalar or vector.
#' @export
transport_plan_mass <- function(
    plan, margin = c("total", "source", "target")) {
  margin <- match.arg(margin)
  plan <- as_transport_plan(plan)
  if (identical(plan$representation, "implicit_operator")) {
    source <- plan$data$source_mass
    target <- plan$data$target_mass
  } else if (identical(plan$representation, "dense_materialized")) {
    source <- rowSums(plan$data)
    target <- colSums(plan$data)
  } else if (identical(plan$representation, "sparse_csc")) {
    source <- as.numeric(Matrix::rowSums(plan$data))
    target <- as.numeric(Matrix::colSums(plan$data))
  } else {
    source <- .weighted_transport_tabulate(
      plan$data$source, plan$data$weight, plan$n_source
    )
    target <- .weighted_transport_tabulate(
      plan$data$target, plan$data$weight, plan$n_target
    )
  }
  switch(margin, total = sum(source), source = source, target = target)
}

.transport_matrix_input <- function(x, expected_rows, name) {
  vector_input <- is.numeric(x) && is.null(dim(x))
  if (vector_input) x <- matrix(x, ncol = 1L)
  x <- .validate_finite_matrix(x, name)
  if (nrow(x) != expected_rows) {
    stop(
      sprintf("`%s` must have %d rows.", name, expected_rows),
      call. = FALSE
    )
  }
  list(value = x, vector = vector_input)
}

.transport_edge_apply <- function(edges, n_output, values, index, group) {
  out <- matrix(0, n_output, ncol(values))
  if (nrow(edges)) {
    for (column in seq_len(ncol(values))) {
      out[, column] <- .weighted_transport_tabulate(
        edges[[group]],
        edges$weight * values[edges[[index]], column],
        n_output
      )
    }
  }
  out
}

#' Apply a Transport Plan as a Linear Operator
#'
#' Computes `plan %*% x` without materializing edge-list, sparse, or implicit
#' plans. This is the forward action from target-indexed values to source rows;
#' it is not row-normalized.
#'
#' @inheritParams transport_plan_shape
#' @param x Numeric vector or matrix with one row per target.
#' @return Numeric vector or matrix with one row per source.
#' @export
transport_plan_apply <- function(plan, x) {
  plan <- as_transport_plan(plan)
  input <- .transport_matrix_input(x, plan$n_target, "x")
  out <- if (identical(plan$representation, "implicit_operator")) {
    plan$data$apply(input$value)
  } else if (identical(plan$representation, "dense_materialized")) {
    plan$data %*% input$value
  } else if (identical(plan$representation, "sparse_csc")) {
    as.matrix(plan$data %*% input$value)
  } else {
    .transport_edge_apply(
      plan$data, plan$n_source, input$value,
      index = "target", group = "source"
    )
  }
  if (is.null(dim(out))) out <- matrix(out, ncol = 1L)
  if (!is.numeric(out) || !identical(dim(out), c(plan$n_source, ncol(input$value))) ||
      any(!is.finite(out))) {
    stop("Transport `apply` returned an incompatible or nonfinite value.", call. = FALSE)
  }
  if (input$vector) as.numeric(out) else unname(out)
}

#' Apply the Transport Adjoint
#'
#' Computes `t(plan) %*% x` without materialization. The adjoint reverses the
#' linear action; it is not an inverse, a reverse conditional distribution, or
#' a barycentric reverse map.
#'
#' @inheritParams transport_plan_shape
#' @param x Numeric vector or matrix with one row per source.
#' @return Numeric vector or matrix with one row per target.
#' @export
transport_plan_adjoint <- function(plan, x) {
  plan <- as_transport_plan(plan)
  input <- .transport_matrix_input(x, plan$n_source, "x")
  out <- if (identical(plan$representation, "implicit_operator")) {
    plan$data$adjoint(input$value)
  } else if (identical(plan$representation, "dense_materialized")) {
    t(plan$data) %*% input$value
  } else if (identical(plan$representation, "sparse_csc")) {
    as.matrix(Matrix::t(plan$data) %*% input$value)
  } else {
    .transport_edge_apply(
      plan$data, plan$n_target, input$value,
      index = "source", group = "target"
    )
  }
  if (is.null(dim(out))) out <- matrix(out, ncol = 1L)
  if (!is.numeric(out) || !identical(dim(out), c(plan$n_target, ncol(input$value))) ||
      any(!is.finite(out))) {
    stop("Transport `adjoint` returned an incompatible or nonfinite value.", call. = FALSE)
  }
  if (input$vector) as.numeric(out) else unname(out)
}

#' Barycentric Projection Through a Transport Plan
#'
#' Normalizes the forward or adjoint action by transported row/column mass.
#' Empty conditionals abstain as `NaN`, zero, or an error; they are never
#' replaced by a uniform match. Reverse barycentric projection is not an
#' inverse coupling.
#'
#' @inheritParams transport_plan_shape
#' @param points Destination point matrix.
#' @param orientation `"source_to_target"` or `"target_to_source"`.
#' @param zero_mass Empty-conditional abstention policy.
#' @return Projected point matrix.
#' @export
transport_plan_barycentric <- function(
    plan, points,
    orientation = c("source_to_target", "target_to_source"),
    zero_mass = c("nan", "zero", "error")) {
  orientation <- match.arg(orientation)
  zero_mass <- match.arg(zero_mass)
  plan <- as_transport_plan(plan)
  if (identical(orientation, "source_to_target")) {
    input <- .transport_matrix_input(points, plan$n_target, "points")
    mass <- transport_plan_mass(plan, "source")
    projected <- transport_plan_apply(plan, input$value)
  } else {
    input <- .transport_matrix_input(points, plan$n_source, "points")
    mass <- transport_plan_mass(plan, "target")
    projected <- transport_plan_adjoint(plan, input$value)
  }
  empty <- mass == 0
  if (any(empty) && identical(zero_mass, "error")) {
    stop("Barycentric projection has zero transported mass.", call. = FALSE)
  }
  fill <- if (identical(zero_mass, "nan")) NaN else 0
  out <- matrix(fill, length(mass), ncol(input$value))
  keep <- !empty
  if (any(keep)) {
    out[keep, ] <- sweep(
      projected[keep, , drop = FALSE], 1L, mass[keep], "/"
    )
  }
  out
}

.transport_plan_edges <- function(plan) {
  plan <- as_transport_plan(plan)
  if (identical(plan$representation, "implicit_operator")) {
    stop(
      "Implicit operators have no enumerable support; materialize explicitly first.",
      call. = FALSE
    )
  }
  if (identical(plan$representation, "edge_list")) return(plan$data)
  if (identical(plan$representation, "sparse_csc")) {
    entries <- Matrix::summary(plan$data)
    weight <- if ("x" %in% names(entries)) entries$x else rep(1, nrow(entries))
    return(.canonical_transport_edges(
      entries$i, entries$j, weight, plan$n_source, plan$n_target, "sum"
    ))
  }
  index <- which(plan$data > 0, arr.ind = TRUE)
  if (!nrow(index)) {
    return(data.frame(
      source = integer(), target = integer(), weight = numeric()
    ))
  }
  .canonical_transport_edges(
    index[, 1L], index[, 2L], plan$data[index],
    plan$n_source, plan$n_target, "sum"
  )
}

#' Explicitly Materialize a Transport Plan
#'
#' This is the only generic operation that intentionally allocates a dense
#' source-by-target matrix. Implicit operators must provide a materializer.
#'
#' @inheritParams transport_plan_shape
#' @return Dense numeric matrix.
#' @export
transport_plan_materialize <- function(plan) {
  plan <- as_transport_plan(plan)
  if (identical(plan$representation, "dense_materialized")) {
    return(plan$data)
  }
  if (identical(plan$representation, "sparse_csc")) {
    return(as.matrix(plan$data))
  }
  if (identical(plan$representation, "implicit_operator")) {
    if (is.null(plan$data$materialize)) {
      stop("This implicit operator has no materializer.", call. = FALSE)
    }
    materialized <- as_transport_plan(plan$data$materialize())
    if (!identical(
      transport_plan_shape(materialized), transport_plan_shape(plan)
    )) {
      stop("Operator materializer returned the wrong shape.", call. = FALSE)
    }
    return(transport_plan_materialize(materialized))
  }
  out <- matrix(0, plan$n_source, plan$n_target)
  if (nrow(plan$data)) {
    out[cbind(plan$data$source, plan$data$target)] <- plan$data$weight
  }
  out
}

#' Inspect Transport Representation
#'
#' @inheritParams transport_plan_shape
#' @return A list naming representation, shape, materialization, sparsity,
#'   implicit/pruned state, and lost mass.
#' @export
transport_plan_representation <- function(plan) {
  plan <- as_transport_plan(plan)
  list(
    representation = plan$representation,
    shape = transport_plan_shape(plan),
    materialized = plan$materialized,
    sparse = plan$sparse,
    implicit = plan$implicit,
    pruned = plan$pruned,
    lost_mass = plan$lost_mass
  )
}

#' Prune a Transport Plan and Invalidate Stale Certificates
#'
#' Removes explicitly stored entries below `min_weight`, reports lost mass, and
#' returns a canonical edge-list plan. When pruning an `rfugw_result`, any
#' positive lost mass clears convergence/feasibility/objective certificates and
#' sets status to `"requires_revalidation"`.
#'
#' @param plan Explicit transport plan or `rfugw_result`.
#' @param min_weight Finite nonnegative retention threshold. Entries strictly
#'   below it are removed.
#' @return Pruned plan or result.
#' @export
transport_plan_prune <- function(plan, min_weight) {
  if (!is.numeric(min_weight) || length(min_weight) != 1L ||
      !is.finite(min_weight) || min_weight < 0) {
    stop("`min_weight` must be one finite nonnegative number.", call. = FALSE)
  }
  is_result <- inherits(plan, "rfugw_result")
  source_plan <- if (is_result) rfugw_plan(plan) else plan
  wrapped <- as_transport_plan(source_plan)
  edges <- .transport_plan_edges(wrapped)
  original_mass <- sum(edges$weight)
  kept <- edges$weight >= min_weight
  pruned_edges <- edges[kept, , drop = FALSE]
  lost_mass <- original_mass - sum(pruned_edges$weight)
  pruned <- as_transport_plan(
    pruned_edges, wrapped$n_source, wrapped$n_target, duplicates = "error"
  )
  pruned$pruned <- lost_mass > 0
  pruned$lost_mass <- lost_mass
  pruned$metadata$pruned <- lost_mass > 0
  pruned$metadata$lost_mass <- lost_mass
  pruned$metadata$original_mass <- original_mass
  pruned$metadata$retained_mass <- original_mass - lost_mass
  pruned$metadata$min_weight <- min_weight
  if (!is_result) return(pruned)

  out <- plan
  out$plan <- pruned
  out$plan_representation <- transport_plan_representation(pruned)
  out$pruning_lost_mass <- lost_mass
  out$certificate_invalidated_by_pruning <- lost_mass > 0
  if (lost_mass > 0) {
    out$converged <- FALSE
    out$feasible <- FALSE
    out$objective_consistent <- FALSE
    out$objective_components_consistent <- FALSE
    out$status <- "requires_revalidation"
    out$termination_reason <- "pruning_changed_plan"
    out$residual <- Inf
    out$error <- Inf
    out$warning_payload <- list(
      code = "requires_revalidation",
      message = sprintf(
        "Pruning removed mass %.6g; all prior plan certificates are stale.",
        lost_mass
      )
    )
  }
  out
}

#' @rdname as_transport_plan
#' @export
dim.rfugw_transport_plan <- function(x) {
  c(x$n_source, x$n_target)
}

#' @rdname as_transport_plan
#' @export
as.matrix.rfugw_transport_plan <- function(x, ...) {
  transport_plan_materialize(x)
}

#' @rdname as_transport_plan
#' @export
print.rfugw_transport_plan <- function(x, ...) {
  cat("<rfugw_transport_plan>\n")
  cat(sprintf("  representation: %s\n", x$representation))
  cat(sprintf("  shape:          %d x %d\n", x$n_source, x$n_target))
  cat(sprintf("  mass:           %s\n", format(transport_plan_mass(x), digits = 6)))
  if (isTRUE(x$pruned)) {
    cat(sprintf("  pruned loss:    %s\n", format(x$lost_mass, digits = 6)))
  }
  invisible(x)
}
