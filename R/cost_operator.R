.validate_cost_shape_scalar <- function(x, name) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) ||
      x < 1 || x != floor(x) || x > .Machine$integer.max) {
    stop(sprintf("`%s` must be a positive integer.", name), call. = FALSE)
  }
  as.integer(x)
}

.validate_cost_index <- function(index, size, name) {
  if (is.null(index)) return(seq_len(size))
  if (!is.numeric(index) || any(!is.finite(index)) ||
      any(index != floor(index)) || any(index < 1 | index > size)) {
    stop(sprintf("`%s` contains invalid indices.", name), call. = FALSE)
  }
  as.integer(index)
}

.validate_cost_factor_matrix <- function(x, name) {
  if (!is.matrix(x) || !is.numeric(x) || nrow(x) < 1L ||
      any(!is.finite(x))) {
    stop(
      sprintf("`%s` must be a finite numeric matrix with at least one row.", name),
      call. = FALSE
    )
  }
  unname(x)
}

.validate_cost_term <- function(x, size, name) {
  if (is.null(x)) return(numeric(size))
  if (!is.numeric(x) || length(x) != size || any(!is.finite(x))) {
    stop(sprintf("`%s` must be a finite numeric vector of length %d.", name, size),
         call. = FALSE)
  }
  as.numeric(x)
}

.new_cost_operator <- function(
    n_source,
    n_target,
    block,
    pair = NULL,
    factors = NULL,
    provenance = list(),
    geometry_reference = NULL,
    geometry_audit_spec = NULL) {
  n_source <- .validate_cost_shape_scalar(n_source, "n_source")
  n_target <- .validate_cost_shape_scalar(n_target, "n_target")
  if (!is.function(block)) {
    stop("`block` must be a function.", call. = FALSE)
  }
  if (!is.null(pair) && !is.function(pair)) {
    stop("`pair` must be NULL or a function.", call. = FALSE)
  }
  if (!is.list(provenance)) {
    stop("`provenance` must be a list.", call. = FALSE)
  }
  structure(
    list(
      n_source = n_source,
      n_target = n_target,
      block = block,
      pair = pair,
      factors = factors,
      provenance = provenance,
      geometry_reference = geometry_reference,
      geometry_audit_spec = geometry_audit_spec
    ),
    class = "rfugw_cost_operator"
  )
}

#' Matrix-Free Cost Operators
#'
#' Cost operators represent a source-by-target cost without requiring a dense
#' matrix. `factorized_cost()` stores an affine-bilinear representation
#'
#' \deqn{C_{ij} = r_i + c_j + L_i R_j^T,}
#'
#' while `block_cost()` accepts a callback for generic block evaluation.
#' Affine-bilinear operators support native matrix-free solvers and exact
#' weighted reductions. Generic block operators remain matrix-free for block
#' evaluation and reductions, but are not accepted by solvers that require
#' algebraic factors.
#'
#' `sqeuclidean_cost()`, `cosine_cost()`, and `correlation_cost()` construct
#' exact affine-bilinear feature costs. `sum_costs()` and `scaled_cost()`
#' preserve factors whenever their inputs are factorized.
#'
#' @param left,right Finite matrices with the same number of columns. Rows
#'   index source and target observations, respectively.
#' @param row,column Optional additive source and target terms.
#' @param constant Finite scalar added to every entry.
#' @param exact Whether the represented cost is exact for the requested
#'   geometry or feature metric.
#' @param provenance Serializable metadata describing construction and any
#'   approximation.
#' @param geometry_reference Optional square requested-geometry matrix or cost
#'   operator retained for an independent approximation audit.
#' @param geometry_audit Optional list describing how a factorized geometry was
#'   constructed and should be audited.
#' @param n_source,n_target Positive operator dimensions.
#' @param evaluate Function of `(rows, columns)` returning the requested
#'   numeric cost block.
#' @param evaluate_pairs Optional function of paired row and column index
#'   vectors returning one finite cost per pair.
#' @param source,target Finite matrices whose rows are observations.
#' @param ... Cost operators, or a single list of cost operators.
#' @param cost A cost operator.
#' @param scale Finite scalar multiplier.
#' @param rows,columns Optional 1-based indices. `NULL` selects a full margin.
#' @param weights Optional finite weights for the opposite margin. The default
#'   is a vector of ones.
#' @param block_size Positive block size for a generic callback reduction.
#' @return An `rfugw_cost_operator`, a numeric matrix block, a weighted
#'   reduction vector, a named shape vector, or a provenance list.
#' @examples
#' x <- matrix(c(0, 0, 1, 0, 0, 2), ncol = 2, byrow = TRUE)
#' y <- matrix(c(0, 1, 2, 0), ncol = 2, byrow = TRUE)
#' cost <- sqeuclidean_cost(x, y)
#' cost_block(cost)
#' cost_row_reduce(cost, c(0.25, 0.75))
#' @name rfugw_cost_operator
NULL

#' @rdname rfugw_cost_operator
#' @export
factorized_cost <- function(
    left,
    right,
    row = NULL,
    column = NULL,
    constant = 0,
    exact = TRUE,
    provenance = list(),
    geometry_reference = NULL,
    geometry_audit = list()) {
  left <- .validate_cost_factor_matrix(left, "left")
  right <- .validate_cost_factor_matrix(right, "right")
  if (ncol(left) != ncol(right)) {
    stop("`left` and `right` must have the same number of columns.",
         call. = FALSE)
  }
  row <- .validate_cost_term(row, nrow(left), "row")
  column <- .validate_cost_term(column, nrow(right), "column")
  if (!is.numeric(constant) || length(constant) != 1L || !is.finite(constant)) {
    stop("`constant` must be a finite scalar.", call. = FALSE)
  }
  if (!is.logical(exact) || length(exact) != 1L || is.na(exact)) {
    stop("`exact` must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.list(provenance)) {
    stop("`provenance` must be a list.", call. = FALSE)
  }
  if (!is.list(geometry_audit)) {
    stop("`geometry_audit` must be a list.", call. = FALSE)
  }
  row <- row + constant
  factors <- list(row = row, column = column, left = left, right = right)
  block <- function(rows, columns) {
    out <- outer(row[rows], column[columns], "+")
    if (ncol(left)) {
      out <- out + tcrossprod(left[rows, , drop = FALSE],
                              right[columns, , drop = FALSE])
    }
    unname(out)
  }
  pair <- function(rows, columns) {
    out <- row[rows] + column[columns]
    if (ncol(left)) {
      out <- out + rowSums(
        left[rows, , drop = FALSE] * right[columns, , drop = FALSE]
      )
    }
    as.numeric(out)
  }
  metadata <- utils::modifyList(
    list(
      representation = "affine_bilinear",
      rank = ncol(left),
      exact = exact,
      constant_folded_into_row = constant
    ),
    provenance
  )
  .new_cost_operator(
    nrow(left), nrow(right), block, pair = pair, factors = factors,
    provenance = metadata,
    geometry_reference = geometry_reference,
    geometry_audit_spec = utils::modifyList(
      list(
        constructor = "generic_factorized",
        constructor_verified_exact = FALSE
      ),
      geometry_audit
    )
  )
}

#' @rdname rfugw_cost_operator
#' @export
block_cost <- function(
    n_source,
    n_target,
    evaluate,
    exact = NA,
    provenance = list(),
    evaluate_pairs = NULL) {
  n_source <- .validate_cost_shape_scalar(n_source, "n_source")
  n_target <- .validate_cost_shape_scalar(n_target, "n_target")
  if (!is.function(evaluate)) {
    stop("`evaluate` must be a function.", call. = FALSE)
  }
  if (!is.null(evaluate_pairs) && !is.function(evaluate_pairs)) {
    stop("`evaluate_pairs` must be NULL or a function.", call. = FALSE)
  }
  if (!is.logical(exact) || length(exact) != 1L) {
    stop("`exact` must be TRUE, FALSE, or NA.", call. = FALSE)
  }
  block <- function(rows, columns) {
    value <- evaluate(rows, columns)
    if (!is.matrix(value) || !is.numeric(value) ||
        nrow(value) != length(rows) || ncol(value) != length(columns) ||
        any(!is.finite(value))) {
      stop(
        "The block callback must return a finite numeric matrix of the requested shape.",
        call. = FALSE
      )
    }
    unname(value)
  }
  metadata <- utils::modifyList(
    list(representation = "block_callback", rank = NA_integer_, exact = exact),
    provenance
  )
  pair <- if (is.null(evaluate_pairs)) NULL else function(rows, columns) {
    value <- evaluate_pairs(rows, columns)
    if (!is.numeric(value) || length(value) != length(rows) ||
        any(!is.finite(value))) {
      stop(
        "The pair callback must return one finite numeric value per pair.",
        call. = FALSE
      )
    }
    as.numeric(value)
  }
  .new_cost_operator(
    n_source, n_target, block, pair = pair, factors = NULL,
    provenance = metadata
  )
}

.feature_matrix <- function(x, name) {
  if (!is.matrix(x) || !is.numeric(x) || nrow(x) < 1L || ncol(x) < 1L ||
      any(!is.finite(x))) {
    stop(sprintf("`%s` must be a nonempty finite numeric matrix.", name),
         call. = FALSE)
  }
  unname(x)
}

#' @rdname rfugw_cost_operator
#' @export
sqeuclidean_cost <- function(source, target = source, provenance = list()) {
  source <- .feature_matrix(source, "source")
  target <- .feature_matrix(target, "target")
  if (ncol(source) != ncol(target)) {
    stop("`source` and `target` must have the same number of columns.",
         call. = FALSE)
  }
  factorized_cost(
    left = -2 * source,
    right = target,
    row = rowSums(source^2),
    column = rowSums(target^2),
    exact = TRUE,
    provenance = utils::modifyList(
      list(metric = "sqeuclidean", input_dimension = ncol(source)),
      provenance
    ),
    geometry_audit = list(
      constructor = "sqeuclidean_cost",
      constructor_verified_exact = TRUE
    )
  )
}

.normalize_feature_rows <- function(x, center, name) {
  if (center && ncol(x) < 2L) {
    stop(sprintf("`%s` needs at least two columns for correlation cost.", name),
         call. = FALSE)
  }
  if (center) x <- x - rowMeans(x)
  norm <- sqrt(rowSums(x^2))
  if (any(!is.finite(norm)) || any(norm <= sqrt(.Machine$double.eps))) {
    stop(sprintf("`%s` contains a zero-norm feature profile.", name),
         call. = FALSE)
  }
  x / norm
}

#' @rdname rfugw_cost_operator
#' @export
cosine_cost <- function(source, target = source, provenance = list()) {
  source <- .feature_matrix(source, "source")
  target <- .feature_matrix(target, "target")
  if (ncol(source) != ncol(target)) {
    stop("`source` and `target` must have the same number of columns.",
         call. = FALSE)
  }
  source <- .normalize_feature_rows(source, FALSE, "source")
  target <- .normalize_feature_rows(target, FALSE, "target")
  factorized_cost(
    -source, target, row = rep(1, nrow(source)), exact = TRUE,
    provenance = utils::modifyList(
      list(metric = "cosine", input_dimension = ncol(source)), provenance
    ),
    geometry_audit = list(
      constructor = "cosine_cost",
      constructor_verified_exact = TRUE
    )
  )
}

#' @rdname rfugw_cost_operator
#' @export
correlation_cost <- function(source, target = source, provenance = list()) {
  source <- .feature_matrix(source, "source")
  target <- .feature_matrix(target, "target")
  if (ncol(source) != ncol(target)) {
    stop("`source` and `target` must have the same number of columns.",
         call. = FALSE)
  }
  source <- .normalize_feature_rows(source, TRUE, "source")
  target <- .normalize_feature_rows(target, TRUE, "target")
  factorized_cost(
    -source, target, row = rep(1, nrow(source)), exact = TRUE,
    provenance = utils::modifyList(
      list(metric = "correlation", input_dimension = ncol(source)), provenance
    ),
    geometry_audit = list(
      constructor = "correlation_cost",
      constructor_verified_exact = TRUE
    )
  )
}

#' @rdname rfugw_cost_operator
#' @export
sum_costs <- function(...) {
  costs <- list(...)
  if (length(costs) == 1L && is.list(costs[[1L]]) &&
      !inherits(costs[[1L]], "rfugw_cost_operator")) {
    costs <- costs[[1L]]
  }
  if (!length(costs) ||
      any(!vapply(costs, inherits, logical(1), "rfugw_cost_operator"))) {
    stop("`...` must contain one or more cost operators.", call. = FALSE)
  }
  shapes <- lapply(costs, cost_shape)
  if (any(!vapply(shapes, identical, logical(1), shapes[[1L]]))) {
    stop("All cost operators must have the same shape.", call. = FALSE)
  }
  factorized <- all(vapply(costs, function(x) !is.null(x$factors), logical(1)))
  provenance <- list(
    representation = if (factorized) "affine_bilinear_sum" else "block_sum",
    exact = all(vapply(costs, function(x) isTRUE(x$provenance$exact), logical(1))),
    components = lapply(costs, cost_provenance)
  )
  if (factorized) {
    factors <- lapply(costs, `[[`, "factors")
    return(factorized_cost(
      left = do.call(cbind, lapply(factors, `[[`, "left")),
      right = do.call(cbind, lapply(factors, `[[`, "right")),
      row = Reduce(`+`, lapply(factors, `[[`, "row")),
      column = Reduce(`+`, lapply(factors, `[[`, "column")),
      exact = provenance$exact,
      provenance = provenance,
      geometry_audit = list(
        constructor = "sum_costs",
        constructor_verified_exact = all(vapply(
          costs,
          function(x) isTRUE(
            x$geometry_audit_spec$constructor_verified_exact
          ),
          logical(1)
        ))
      )
    ))
  }
  block_cost(
    shapes[[1L]][["source"]], shapes[[1L]][["target"]],
    function(rows, columns) {
      Reduce(`+`, lapply(costs, cost_block, rows = rows, columns = columns))
    },
    exact = provenance$exact,
    provenance = provenance
  )
}

#' @rdname rfugw_cost_operator
#' @export
scaled_cost <- function(cost, scale) {
  if (!inherits(cost, "rfugw_cost_operator")) {
    stop("`cost` must be an rfugw cost operator.", call. = FALSE)
  }
  if (!is.numeric(scale) || length(scale) != 1L || !is.finite(scale)) {
    stop("`scale` must be a finite scalar.", call. = FALSE)
  }
  provenance <- list(
    representation = paste0(cost$provenance$representation, "_scaled"),
    exact = cost$provenance$exact,
    scale = scale,
    input = cost_provenance(cost)
  )
  if (!is.null(cost$factors)) {
    return(factorized_cost(
      left = scale * cost$factors$left,
      right = cost$factors$right,
      row = scale * cost$factors$row,
      column = scale * cost$factors$column,
      exact = isTRUE(cost$provenance$exact),
      provenance = provenance,
      geometry_audit = list(
        constructor = "scaled_cost",
        constructor_verified_exact = isTRUE(
          cost$geometry_audit_spec$constructor_verified_exact
        )
      )
    ))
  }
  shape <- cost_shape(cost)
  block_cost(
    shape[["source"]], shape[["target"]],
    function(rows, columns) scale * cost_block(cost, rows, columns),
    exact = cost$provenance$exact,
    provenance = provenance
  )
}

#' @rdname rfugw_cost_operator
#' @export
cost_shape <- function(cost) {
  if (!inherits(cost, "rfugw_cost_operator")) {
    stop("`cost` must be an rfugw cost operator.", call. = FALSE)
  }
  c(source = cost$n_source, target = cost$n_target)
}

#' @rdname rfugw_cost_operator
#' @export
cost_block <- function(cost, rows = NULL, columns = NULL) {
  shape <- cost_shape(cost)
  rows <- .validate_cost_index(rows, shape[["source"]], "rows")
  columns <- .validate_cost_index(columns, shape[["target"]], "columns")
  cost$block(rows, columns)
}

.cost_pair_values <- function(cost, rows, columns, block_size = 32L) {
  shape <- cost_shape(cost)
  rows <- .validate_cost_index(rows, shape[["source"]], "rows")
  columns <- .validate_cost_index(columns, shape[["target"]], "columns")
  if (length(rows) != length(columns)) {
    stop("`rows` and `columns` must identify the same number of pairs.",
         call. = FALSE)
  }
  if (!length(rows)) return(numeric())
  if (is.function(cost$pair)) {
    return(cost$pair(rows, columns))
  }
  block_size <- .validate_count(block_size, "block_size")
  out <- numeric(length(rows))
  starts <- seq.int(1L, length(rows), by = block_size)
  for (start in starts) {
    index <- seq.int(start, min(length(rows), start + block_size - 1L))
    value <- cost_block(cost, rows[index], columns[index])
    out[index] <- diag(value)
  }
  out
}

.geometry_reference_operator <- function(reference, n) {
  if (inherits(reference, "rfugw_cost_operator")) {
    shape <- cost_shape(reference)
    if (!identical(unname(shape), c(n, n))) {
      stop("The geometry reference must have the same square shape.",
           call. = FALSE)
    }
    return(reference)
  }
  if (is.matrix(reference) && is.numeric(reference) &&
      identical(dim(reference), c(n, n)) && all(is.finite(reference))) {
    reference <- unname(reference)
    return(.new_cost_operator(
      n,
      n,
      block = function(rows, columns) {
        reference[rows, columns, drop = FALSE]
      },
      pair = function(rows, columns) {
        as.numeric(reference[cbind(rows, columns)])
      },
      provenance = list(
        representation = "explicit_geometry_reference",
        exact = TRUE
      )
    ))
  }
  if (is.function(reference)) {
    return(block_cost(
      n,
      n,
      reference,
      exact = NA,
      provenance = list(representation = "geometry_reference_callback")
    ))
  }
  stop(
    "`reference` must be a square numeric matrix, cost operator, or block callback.",
    call. = FALSE
  )
}

#' Embedded Geometry with a Computed Approximation Audit
#'
#' Constructs squared-Euclidean geometry in an embedding while retaining a
#' separate requested-geometry operator for reproducible train and held-out
#' error audits. The reference can be an explicit matrix, a cost operator, or
#' a block callback. It is never copied into the factorized solver state.
#'
#' @param embedding Finite matrix with one embedded coordinate vector per row.
#' @param reference Requested square geometry as a numeric matrix,
#'   [rfugw_cost_operator], or block callback.
#' @param method Nonempty embedding-method label.
#' @param reference_metric Nonempty requested-metric label.
#' @param seed Integer seed for the reproducible pair design.
#' @param train_pairs,holdout_pairs Positive pair-sample counts.
#' @param relative_error_tolerance Maximum held-out weighted stress accepted by
#'   the approximation certificate.
#' @param provenance Additional serializable provenance.
#' @return A factorized squared-Euclidean cost carrying a geometry audit
#'   contract and reference operator.
#' @export
embedded_geometry_cost <- function(
    embedding,
    reference,
    method = "user_embedding",
    reference_metric = "requested_geometry",
    seed = 104729L,
    train_pairs = 2048L,
    holdout_pairs = 2048L,
    relative_error_tolerance = 0.1,
    provenance = list()) {
  embedding <- .feature_matrix(embedding, "embedding")
  if (!is.character(method) || length(method) != 1L || is.na(method) ||
      !nzchar(method) || !is.character(reference_metric) ||
      length(reference_metric) != 1L || is.na(reference_metric) ||
      !nzchar(reference_metric)) {
    stop("Embedding method and reference metric must be nonempty labels.",
         call. = FALSE)
  }
  seed <- .validate_count(seed, "seed")
  train_pairs <- .validate_count(train_pairs, "train_pairs")
  holdout_pairs <- .validate_count(holdout_pairs, "holdout_pairs")
  relative_error_tolerance <- .validate_positive_scalar(
    relative_error_tolerance, "relative_error_tolerance"
  )
  reference <- .geometry_reference_operator(reference, nrow(embedding))
  cost <- sqeuclidean_cost(
    embedding,
    provenance = utils::modifyList(
      list(
        metric = "embedded_sqeuclidean",
        embedding_method = method,
        embedding_rank = ncol(embedding),
        reference_metric = reference_metric
      ),
      provenance
    )
  )
  cost$provenance$exact <- FALSE
  cost$geometry_reference <- reference
  cost$geometry_audit_spec <- list(
    constructor = "embedded_geometry_cost",
    constructor_verified_exact = FALSE,
    embedding_method = method,
    embedding_rank = ncol(embedding),
    seed = seed,
    train_pairs = train_pairs,
    holdout_pairs = holdout_pairs,
    relative_error_tolerance = relative_error_tolerance,
    reference_metric = reference_metric,
    reference_provenance = cost_provenance(reference)
  )
  cost
}

.with_geometry_seed <- function(seed, expression) {
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv)
  on.exit({
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  set.seed(seed)
  force(expression)
}

.geometry_all_pairs <- function(n) {
  if (n < 2L) return(matrix(integer(), nrow = 0L, ncol = 2L))
  t(utils::combn(n, 2L))
}

.geometry_sample_pairs <- function(n, count, seed) {
  total <- as.double(n) * as.double(n - 1L) / 2
  count <- as.integer(min(as.double(count), total))
  if (!count) return(matrix(integer(), nrow = 0L, ncol = 2L))
  if (total <= count) return(.geometry_all_pairs(n))
  .with_geometry_seed(seed, {
    pairs <- matrix(integer(), nrow = 0L, ncol = 2L)
    while (nrow(pairs) < count) {
      needed <- count - nrow(pairs)
      batch <- max(64L, min(.Machine$integer.max, 3L * needed))
      first <- sample.int(n, batch, replace = TRUE)
      second <- sample.int(n, batch, replace = TRUE)
      keep <- first != second
      candidate <- cbind(
        pmin.int(first[keep], second[keep]),
        pmax.int(first[keep], second[keep])
      )
      pairs <- unique(rbind(pairs, candidate))
    }
    pairs[seq_len(count), , drop = FALSE]
  })
}

.geometry_error_summary <- function(value, reference, pair_weights) {
  if (!length(value)) {
    return(list(weighted_stress = 0, relative_distance_error = 0,
                maximum_relative_error = 0))
  }
  finite <- all(is.finite(value)) && all(is.finite(reference)) &&
    all(is.finite(pair_weights)) && all(pair_weights >= 0)
  if (!finite) {
    return(list(weighted_stress = Inf, relative_distance_error = Inf,
                maximum_relative_error = Inf))
  }
  if (sum(pair_weights) <= 0) pair_weights <- rep(1, length(value))
  difference <- value - reference
  squared_denominator <- sum(pair_weights * reference^2)
  absolute_denominator <- sum(pair_weights * abs(reference))
  scale <- max(1, max(abs(reference)))
  list(
    weighted_stress = sqrt(
      sum(pair_weights * difference^2) /
        max(squared_denominator, .Machine$double.eps)
    ),
    relative_distance_error =
      sum(pair_weights * abs(difference)) /
        max(absolute_denominator, .Machine$double.eps),
    maximum_relative_error = max(abs(difference)) / scale
  )
}

#' Audit a Matrix-Free Geometry Operator
#'
#' Computes structural invariants and, when a requested-geometry reference is
#' attached, weighted train and held-out approximation errors. Small domains
#' enumerate all unordered pairs without assembling a matrix; large domains
#' use a disjoint reproducible pair design.
#'
#' @param cost Square geometry cost operator.
#' @param weights Optional nonnegative node weights.
#' @param seed Optional audit seed overriding constructor provenance.
#' @param train_pairs,holdout_pairs Optional pair counts.
#' @param exact_pair_limit Maximum number of unordered pairs enumerated exactly.
#' @param block_size Pair-callback fallback tile size.
#' @param relative_error_tolerance Optional held-out stress threshold.
#' @return A layered geometry audit certificate.
#' @export
geometry_audit <- function(
    cost,
    weights = NULL,
    seed = NULL,
    train_pairs = NULL,
    holdout_pairs = NULL,
    exact_pair_limit = 100000L,
    block_size = 32L,
    relative_error_tolerance = NULL) {
  shape <- cost_shape(cost)
  if (shape[["source"]] != shape[["target"]]) {
    stop("Geometry audits require a square cost operator.", call. = FALSE)
  }
  n <- unname(shape[["source"]])
  if (is.null(weights)) weights <- rep(1 / n, n)
  if (!is.numeric(weights) || length(weights) != n ||
      any(!is.finite(weights)) || any(weights < 0) || sum(weights) <= 0) {
    stop("`weights` must be finite, nonnegative, and have positive mass.",
         call. = FALSE)
  }
  weights <- as.numeric(weights / sum(weights))
  exact_pair_limit <- .validate_count(exact_pair_limit, "exact_pair_limit")
  block_size <- .validate_count(block_size, "block_size")
  spec <- cost$geometry_audit_spec %||% list()
  seed <- seed %||% spec$seed %||% 104729L
  train_pairs <- train_pairs %||% spec$train_pairs %||% 2048L
  holdout_pairs <- holdout_pairs %||% spec$holdout_pairs %||% 2048L
  relative_error_tolerance <- relative_error_tolerance %||%
    spec$relative_error_tolerance %||% 0.1
  seed <- .validate_count(seed, "seed")
  train_pairs <- .validate_count(train_pairs, "train_pairs")
  holdout_pairs <- .validate_count(holdout_pairs, "holdout_pairs")
  relative_error_tolerance <- .validate_positive_scalar(
    relative_error_tolerance, "relative_error_tolerance"
  )

  total_pairs <- as.double(n) * as.double(n - 1L) / 2
  exact_mode <- total_pairs <= exact_pair_limit
  if (exact_mode) {
    all_pairs <- .geometry_all_pairs(n)
    train <- all_pairs
    heldout <- all_pairs
  } else {
    combined_count <- train_pairs + holdout_pairs
    combined <- .geometry_sample_pairs(n, combined_count, seed)
    train_count <- min(train_pairs, nrow(combined))
    train <- combined[seq_len(train_count), , drop = FALSE]
    heldout_index <- if (nrow(combined) > train_count) {
      seq.int(train_count + 1L, nrow(combined))
    } else integer()
    heldout <- combined[heldout_index, , drop = FALSE]
  }
  audited_pairs <- if (exact_mode) train else rbind(train, heldout)
  forward <- .cost_pair_values(
    cost, audited_pairs[, 1L], audited_pairs[, 2L], block_size
  )
  reverse <- .cost_pair_values(
    cost, audited_pairs[, 2L], audited_pairs[, 1L], block_size
  )
  diagonal <- .cost_pair_values(cost, seq_len(n), seq_len(n), block_size)
  value_scale <- max(1, abs(forward), abs(reverse), abs(diagonal))
  invariant_tolerance <- max(
    1e-10,
    1000 * .Machine$double.eps * value_scale
  )
  symmetry_residual <- if (length(forward)) {
    max(abs(forward - reverse))
  } else 0
  nonnegativity_violation <- max(0, -min(c(forward, reverse, diagonal, 0)))
  diagonal_residual <- max(abs(diagonal))
  invariants_certified <- all(is.finite(c(
    symmetry_residual, nonnegativity_violation, diagonal_residual
  ))) && symmetry_residual <= invariant_tolerance &&
    nonnegativity_violation <= invariant_tolerance &&
    diagonal_residual <= invariant_tolerance

  reference <- cost$geometry_reference
  reference_available <- inherits(reference, "rfugw_cost_operator")
  train_position <- seq_len(nrow(train))
  heldout_position <- if (exact_mode) {
    train_position
  } else if (nrow(heldout)) {
    seq.int(nrow(train) + 1L, nrow(train) + nrow(heldout))
  } else integer()
  if (reference_available) {
    reference_forward <- .cost_pair_values(
      reference,
      audited_pairs[, 1L],
      audited_pairs[, 2L],
      block_size
    )
    train_error <- .geometry_error_summary(
      forward[train_position],
      reference_forward[train_position],
      weights[train[, 1L]] * weights[train[, 2L]]
    )
    heldout_error <- .geometry_error_summary(
      forward[heldout_position],
      reference_forward[heldout_position],
      weights[heldout[, 1L]] * weights[heldout[, 2L]]
    )
  } else {
    train_error <- heldout_error <- list(
      weighted_stress = NA_real_,
      relative_distance_error = NA_real_,
      maximum_relative_error = NA_real_
    )
  }

  provenance <- cost_provenance(cost)
  reported_exact <- isTRUE(provenance$exact)
  constructor_verified <- isTRUE(spec$constructor_verified_exact)
  exact_reference_consistent <- reference_available &&
    is.finite(train_error$weighted_stress) &&
    is.finite(heldout_error$weighted_stress) &&
    max(train_error$weighted_stress, heldout_error$weighted_stress) <=
      invariant_tolerance
  exact_verified <- reported_exact && invariants_certified &&
    if (reference_available) exact_reference_consistent else constructor_verified
  approximation_certified <- !reported_exact && reference_available &&
    invariants_certified && is.finite(heldout_error$weighted_stress) &&
    heldout_error$weighted_stress <= relative_error_tolerance
  certified <- exact_verified || approximation_certified
  status <- if (!invariants_certified) {
    "invalid_geometry_invariants"
  } else if (reported_exact && !exact_verified) {
    "unverified_or_falsified_exact_claim"
  } else if (!reported_exact && !reference_available) {
    "approximation_unverified_missing_reference"
  } else if (!reported_exact && !approximation_certified) {
    "heldout_approximation_failure"
  } else if (exact_verified) {
    "exact_verified"
  } else {
    "approximate_verified"
  }
  pair_design <- list(
    mode = if (exact_mode) "exact_all_unordered_pairs" else {
      "reproducible_disjoint_train_holdout"
    },
    seed = seed,
    total_available_pairs = total_pairs,
    train_pair_count = nrow(train),
    holdout_pair_count = nrow(heldout),
    train_pairs = if (exact_mode) NULL else train,
    holdout_pairs = if (exact_mode) NULL else heldout
  )
  reported_error <- provenance$relative_error %||% NA_real_
  relative_error <- if (reference_available) {
    heldout_error$weighted_stress
  } else if (is.numeric(reported_error) && length(reported_error) == 1L) {
    reported_error
  } else {
    NA_real_
  }

  list(
    certified = certified,
    status = status,
    exact = exact_verified,
    approximate = approximation_certified,
    relative_error = relative_error,
    relative_error_verified = reference_available,
    relative_error_tolerance = relative_error_tolerance,
    weighted_stress = if (reference_available) {
      train_error$weighted_stress
    } else NA_real_,
    heldout_error = if (reference_available) {
      heldout_error$weighted_stress
    } else NA_real_,
    train_error = train_error,
    heldout_error_components = heldout_error,
    symmetry_residual = symmetry_residual,
    nonnegativity_violation = nonnegativity_violation,
    diagonal_residual = diagonal_residual,
    invariant_tolerance = invariant_tolerance,
    invariants_certified = invariants_certified,
    reported_exact = reported_exact,
    constructor_verified_exact = constructor_verified,
    reference_available = reference_available,
    reference_metric = spec$reference_metric %||%
      provenance$metric %||% "unspecified",
    embedding_method = spec$embedding_method %||% NA_character_,
    embedding_rank = spec$embedding_rank %||% NA_integer_,
    pair_design = pair_design,
    provenance = provenance,
    reference_provenance = if (reference_available) {
      cost_provenance(reference)
    } else NULL,
    interpretation = paste0(
      "Geometry error is independent of solver residual and support error. ",
      "Exactness requires computed invariants plus a verified constructor or ",
      "an agreeing requested-geometry reference."
    )
  )
}

.cost_reduction_block_size <- function(block_size) {
  .validate_count(block_size, "block_size")
}

#' @rdname rfugw_cost_operator
#' @export
cost_row_reduce <- function(cost, weights = NULL, block_size = 256L) {
  shape <- cost_shape(cost)
  if (is.null(weights)) weights <- rep(1, shape[["target"]])
  weights <- .validate_cost_term(weights, shape[["target"]], "weights")
  if (!is.null(cost$factors)) {
    factors <- cost$factors
    return(as.numeric(
      factors$row * sum(weights) + sum(factors$column * weights) +
        if (ncol(factors$left)) {
          factors$left %*% crossprod(factors$right, weights)
        } else {
          0
        }
    ))
  }
  block_size <- .cost_reduction_block_size(block_size)
  out <- numeric(shape[["source"]])
  row_starts <- seq.int(1L, shape[["source"]], by = block_size)
  col_starts <- seq.int(1L, shape[["target"]], by = block_size)
  for (row_start in row_starts) {
    rows <- seq.int(row_start, min(shape[["source"]], row_start + block_size - 1L))
    value <- numeric(length(rows))
    for (col_start in col_starts) {
      columns <- seq.int(
        col_start, min(shape[["target"]], col_start + block_size - 1L)
      )
      value <- value + as.numeric(
        cost_block(cost, rows, columns) %*% weights[columns]
      )
    }
    out[rows] <- value
  }
  out
}

#' @rdname rfugw_cost_operator
#' @export
cost_col_reduce <- function(cost, weights = NULL, block_size = 256L) {
  shape <- cost_shape(cost)
  if (is.null(weights)) weights <- rep(1, shape[["source"]])
  weights <- .validate_cost_term(weights, shape[["source"]], "weights")
  if (!is.null(cost$factors)) {
    factors <- cost$factors
    return(as.numeric(
      factors$column * sum(weights) + sum(factors$row * weights) +
        if (ncol(factors$left)) {
          factors$right %*% crossprod(factors$left, weights)
        } else {
          0
        }
    ))
  }
  block_size <- .cost_reduction_block_size(block_size)
  out <- numeric(shape[["target"]])
  row_starts <- seq.int(1L, shape[["source"]], by = block_size)
  col_starts <- seq.int(1L, shape[["target"]], by = block_size)
  for (col_start in col_starts) {
    columns <- seq.int(
      col_start, min(shape[["target"]], col_start + block_size - 1L)
    )
    value <- numeric(length(columns))
    for (row_start in row_starts) {
      rows <- seq.int(
        row_start, min(shape[["source"]], row_start + block_size - 1L)
      )
      value <- value + as.numeric(
        crossprod(cost_block(cost, rows, columns), weights[rows])
      )
    }
    out[columns] <- value
  }
  out
}

#' @rdname rfugw_cost_operator
#' @export
cost_provenance <- function(cost) {
  cost_shape(cost)
  cost$provenance
}

.cost_native_terms <- function(cost, name = "cost") {
  if (!inherits(cost, "rfugw_cost_operator")) {
    stop(sprintf("`%s` must be an rfugw cost operator.", name), call. = FALSE)
  }
  if (is.null(cost$factors)) {
    stop(
      sprintf(
        "`%s` is block-evaluable but has no affine-bilinear factors required by this solver.",
        name
      ),
      call. = FALSE
    )
  }
  cost$factors
}

.cost_factor_pair <- function(cost, name = "cost") {
  factors <- .cost_native_terms(cost, name)
  n_source <- length(factors$row)
  n_target <- length(factors$column)
  left <- cbind(factors$row, rep(1, n_source), factors$left)
  right <- cbind(rep(1, n_target), factors$column, factors$right)
  keep <- vapply(seq_len(ncol(left)), function(index) {
    any(left[, index] != 0) && any(right[, index] != 0)
  }, logical(1))
  list(
    left = left[, keep, drop = FALSE],
    right = right[, keep, drop = FALSE]
  )
}

.zero_factorized_cost <- function(n_source, n_target, provenance = list()) {
  factorized_cost(
    matrix(numeric(), nrow = n_source, ncol = 0L),
    matrix(numeric(), nrow = n_target, ncol = 0L),
    provenance = utils::modifyList(list(metric = "zero"), provenance)
  )
}

#' @export
print.rfugw_cost_operator <- function(x, ...) {
  shape <- cost_shape(x)
  cat("<rfugw_cost_operator>\n")
  cat(sprintf("  shape:          %d x %d\n", shape[[1L]], shape[[2L]]))
  cat(sprintf("  representation: %s\n", x$provenance$representation %||% "unknown"))
  if (!is.null(x$provenance$rank) && !is.na(x$provenance$rank)) {
    cat(sprintf("  bilinear rank:  %d\n", x$provenance$rank))
  }
  cat(sprintf("  exact:          %s\n", format(x$provenance$exact %||% NA)))
  invisible(x)
}
