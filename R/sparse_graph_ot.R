.is_sparse_graph_matrix <- function(x) {
  inherits(x, "sparseMatrix")
}

.validate_similarity_graph <- function(W, name = "W") {
  if (is.matrix(W)) {
    if (!is.numeric(W)) {
      stop(sprintf("`%s` must be numeric.", name), call. = FALSE)
    }
    if (any(!is.finite(W)) || any(W < 0)) {
      stop(sprintf("`%s` must be finite and nonnegative.", name), call. = FALSE)
    }
    return(W)
  }

  if (.is_sparse_graph_matrix(W)) {
    xvals <- methods::slot(W, "x")
    if (length(xvals) > 0L && (any(!is.finite(xvals)) || any(xvals < 0))) {
      stop(sprintf("`%s` must be finite and nonnegative.", name), call. = FALSE)
    }
    return(W)
  }

  stop(sprintf("`%s` must be a numeric matrix or sparseMatrix.", name), call. = FALSE)
}

.canonicalize_eigvec_sign <- function(V) {
  if (!is.matrix(V) || ncol(V) == 0L) {
    return(V)
  }
  out <- V
  for (j in seq_len(ncol(out))) {
    i <- which.max(abs(out[, j]))
    if (is.finite(out[i, j]) && out[i, j] < 0) {
      out[, j] <- -out[, j]
    }
  }
  out
}

.row_distance_from_coords <- function(X, norms2, i, metric = c("euclidean", "sqeuclidean")) {
  metric <- match.arg(metric)
  xi <- X[i, , drop = FALSE]
  d2 <- norms2 + norms2[i] - 2 * as.vector(X %*% t(xi))
  d2[!is.finite(d2)] <- 0
  d2[d2 < 0] <- 0
  if (identical(metric, "euclidean")) {
    return(sqrt(d2))
  }
  d2
}

.deterministic_topk_indices <- function(prob, k) {
  prob <- as.numeric(prob)
  n <- length(prob)
  if (n == 0L) return(integer())
  kk <- as.integer(min(max(0L, k), n))
  if (kk <= 0L) return(integer())
  prob[!is.finite(prob) | prob < 0] <- 0
  ord <- order(-prob, seq_len(n))
  ord[seq_len(kk)]
}

.estimate_gw_square_from_coords <- function(X1, X2, T, metric = c("euclidean", "sqeuclidean"), n_samples = 2000L) {
  metric <- match.arg(metric)
  ns <- nrow(X1)
  nt <- nrow(X2)
  n_samples <- as.integer(n_samples)[1]
  if (!is.finite(n_samples) || n_samples < 1L) {
    return(NA_real_)
  }

  prob <- as.vector(T)
  sprob <- sum(prob)
  if (!is.finite(sprob) || sprob <= 0) {
    return(NA_real_)
  }
  prob <- prob / sprob

  idx_a <- sample.int(length(prob), size = n_samples, replace = TRUE, prob = prob)
  idx_b <- sample.int(length(prob), size = n_samples, replace = TRUE, prob = prob)

  ia <- ((idx_a - 1L) %% ns) + 1L
  ja <- ((idx_a - 1L) %/% ns) + 1L
  ib <- ((idx_b - 1L) %% ns) + 1L
  jb <- ((idx_b - 1L) %/% ns) + 1L

  d1_sq <- rowSums((X1[ia, , drop = FALSE] - X1[ib, , drop = FALSE])^2)
  d2_sq <- rowSums((X2[ja, , drop = FALSE] - X2[jb, , drop = FALSE])^2)

  d1_sq[d1_sq < 0] <- 0
  d2_sq[d2_sq < 0] <- 0
  if (identical(metric, "euclidean")) {
    d1 <- sqrt(d1_sq)
    d2 <- sqrt(d2_sq)
  } else {
    d1 <- d1_sq
    d2 <- d2_sq
  }
  mean((d1 - d2)^2)
}

#' Diffusion Coordinates from a Similarity Graph
#'
#' Build diffusion-map style coordinates from a dense or sparse nonnegative
#' similarity graph. This provides a memory-aware route from sparse kNN graphs
#' to geometry vectors that can be aligned with sampled GW without materializing
#' dense `n x n` structure matrices.
#'
#' @param W Similarity graph (`n x n`) as a dense numeric matrix or a
#'   sparse `Matrix::sparseMatrix`.
#' @param n_components Number of diffusion coordinates to return.
#' @param diffusion_time Diffusion-time exponent applied to eigenvalues.
#' @param symmetrize If `TRUE`, symmetrize graph as `(W + t(W))/2`.
#' @param self_loop Optional diagonal value added to all nodes.
#' @param tol Degree floor used for validity checks.
#' @param log If `TRUE`, return diagnostics (`eigenvalues`, `degree`, `solver`).
#' @return If `log = FALSE`, returns a dense matrix of coordinates (`n x k`).
#'   If `log = TRUE`, returns a list with `coords`, `eigenvalues`, `degree`,
#'   and `solver`.
#' @export
graph_diffusion_coordinates <- function(
    W,
    n_components = 16L,
    diffusion_time = 1,
    symmetrize = TRUE,
    self_loop = 0,
    tol = 1e-12,
    log = FALSE) {
  W <- .validate_similarity_graph(W, "W")
  sparse_input <- .is_sparse_graph_matrix(W)
  if (sparse_input && !requireNamespace("Matrix", quietly = TRUE)) {
    stop("Sparse graph input requires package `Matrix`.", call. = FALSE)
  }

  n <- nrow(W)
  if (ncol(W) != n) {
    stop("`W` must be square.", call. = FALSE)
  }
  if (n < 2L) {
    stop("`W` must contain at least 2 nodes.", call. = FALSE)
  }

  n_components <- as.integer(n_components)[1]
  if (!is.finite(n_components) || n_components < 1L) {
    stop("`n_components` must be an integer >= 1.", call. = FALSE)
  }
  diffusion_time <- as.numeric(diffusion_time)[1]
  if (!is.finite(diffusion_time) || diffusion_time < 0) {
    stop("`diffusion_time` must be finite and >= 0.", call. = FALSE)
  }
  tol <- as.numeric(tol)[1]
  if (!is.finite(tol) || tol <= 0) {
    stop("`tol` must be finite and > 0.", call. = FALSE)
  }
  self_loop <- as.numeric(self_loop)[1]
  if (!is.finite(self_loop) || self_loop < 0) {
    stop("`self_loop` must be finite and >= 0.", call. = FALSE)
  }

  if (isTRUE(symmetrize)) {
    if (sparse_input) {
      W <- (W + Matrix::t(W)) * 0.5
    } else {
      W <- (W + t(W)) * 0.5
    }
  }

  if (self_loop > 0) {
    if (sparse_input) {
      W <- W + Matrix::Diagonal(n = n, x = rep(self_loop, n))
    } else {
      diag(W) <- diag(W) + self_loop
    }
  }

  degree <- if (sparse_input) as.numeric(Matrix::rowSums(W)) else rowSums(W)
  if (any(!is.finite(degree)) || any(degree <= tol)) {
    stop("`W` has isolated or near-isolated nodes; add `self_loop` or drop isolates.", call. = FALSE)
  }
  d_inv_sqrt <- 1 / sqrt(pmax(degree, tol))

  S <- if (sparse_input) {
    D <- Matrix::Diagonal(x = d_inv_sqrt)
    D %*% W %*% D
  } else {
    sweep(sweep(W, 1, d_inv_sqrt, `*`), 2, d_inv_sqrt, `*`)
  }

  k_eff <- min(as.integer(n - 1L), n_components)
  if (k_eff < 1L) {
    stop("Need at least one non-trivial diffusion component.", call. = FALSE)
  }
  k_solver <- min(n, k_eff + 1L)

  eig <- NULL
  eig_solver <- "base::eigen"
  if (requireNamespace("RSpectra", quietly = TRUE) && k_solver < n) {
    eig <- tryCatch(
      RSpectra::eigs_sym(S, k = k_solver, which = "LM"),
      error = function(e) NULL
    )
    if (!is.null(eig)) {
      eig_solver <- "RSpectra::eigs_sym"
    }
  }
  if (is.null(eig)) {
    eig <- eigen(as.matrix(S), symmetric = TRUE)
  }

  vals <- Re(eig$values)
  vecs <- Re(eig$vectors)
  ord <- order(vals, decreasing = TRUE)
  vals <- vals[ord]
  vecs <- vecs[, ord, drop = FALSE]

  keep <- seq.int(2L, min(length(vals), k_eff + 1L))
  if (length(keep) < 1L) {
    stop("Could not compute non-trivial diffusion components.", call. = FALSE)
  }
  vals_keep <- pmax(vals[keep], 0)
  vecs_keep <- .canonicalize_eigvec_sign(vecs[, keep, drop = FALSE])
  coords <- sweep(vecs_keep, 2, vals_keep^diffusion_time, `*`)
  coords <- as.matrix(coords)

  if (!isTRUE(log)) {
    return(coords)
  }
  list(
    coords = coords,
    eigenvalues = vals_keep,
    degree = degree,
    solver = eig_solver
  )
}

#' Sampled Gromov-Wasserstein Directly from Coordinates
#'
#' POT-style sampled GW solver that operates on point coordinates and computes
#' required structure-distance rows on the fly. This avoids storing dense
#' `n x n` distance matrices and is suitable after sparse-graph diffusion
#' embeddings.
#'
#' @param X1 Source coordinates (`ns x d`).
#' @param X2 Target coordinates (`nt x d`).
#' @param p Source weights (default uniform).
#' @param q Target weights (default uniform).
#' @param metric Ground metric between coordinates (`"euclidean"` or
#'   `"sqeuclidean"`).
#' @param loss_fun Currently only `"square_loss"` is supported.
#' @param nb_samples_grad Number of sampled gradient points, or length-2 vector
#'   `(n_source_samples, n_target_samples)`. Values below 1 error. Source or
#'   target counts above `ns` / `nt` warn and clamp.
#' @param epsilon Entropic regularization for projection step.
#' @param max_iter Maximum stochastic iterations.
#' @param log If `TRUE`, return diagnostics.
#' @param verbose If `TRUE`, print iterative diagnostics.
#' @param random_state Optional seed.
#' @param sinkhorn_max_iter Sinkhorn iterations when `epsilon > 0`.
#' @param sinkhorn_tol Sinkhorn tolerance when `epsilon > 0`.
#' @param lp_solver LP backend when `epsilon <= 0`.
#' @param lp_scale Integer scaling for LP marginals.
#' @param use_cpp If `TRUE` and C++ kernel is available, use the optimized
#'   entropic coordinate-native kernel when `epsilon > 0`.
#' @param objective_samples Number of Monte Carlo samples for
#'   `gw_dist_estimated` when `log = TRUE`.
#' @param sampling Sampling policy for stochastic updates. `"stochastic"` uses
#'   weighted random sampling (POT-style). `"deterministic"` uses top-k
#'   probability selections for reproducible parity checks.
#' @return If `log = FALSE`, returns coupling matrix `T`. If `log = TRUE`,
#'   returns a list with `plan`, `gw_dist_estimated`, and `iterations`.
#'
#' @section Experimental:
#' Coordinate-native sampled GW is experimental. The certified 0.1
#' envelope matches [sampled_gromov_wasserstein()]: a full budget
#' `(ns, nt)` is closer to dense entropic GW than a tiny budget.
#' Intermediate budgets are not certified as monotone. Inputs scale as
#' `O(n d)` rather than `O(n^2)` structure costs. See
#' `inst/bench/sampled-budget-curves.md`.
#' @export
sampled_gromov_wasserstein_coords <- function(
    X1,
    X2,
    p = NULL,
    q = NULL,
    metric = c("euclidean", "sqeuclidean"),
    loss_fun = "square_loss",
    nb_samples_grad = 100L,
    epsilon = 1,
    max_iter = 500L,
    log = FALSE,
    verbose = FALSE,
    random_state = NULL,
    sinkhorn_max_iter = 200L,
    sinkhorn_tol = 1e-9,
    lp_solver = c("lp_matrix", "lp_transport"),
    lp_scale = 1e6,
    use_cpp = TRUE,
    objective_samples = 2000L,
    sampling = c("stochastic", "deterministic")) {
  .check_square_loss(loss_fun)
  .assert_matrix(X1, "X1")
  .assert_matrix(X2, "X2")
  metric <- match.arg(metric)
  sampling <- match.arg(sampling)
  deterministic_sampling <- identical(sampling, "deterministic")

  if (ncol(X1) != ncol(X2)) {
    stop("`X1` and `X2` must have the same number of columns.", call. = FALSE)
  }
  ns <- nrow(X1)
  nt <- nrow(X2)
  if (ns < 2L || nt < 2L) {
    stop("`X1` and `X2` must each have at least 2 rows.", call. = FALSE)
  }

  if (is.null(p)) p <- rep(1 / ns, ns)
  if (is.null(q)) q <- rep(1 / nt, nt)
  p <- .assert_prob(p, ns, "p")
  q <- .assert_prob(q, nt, "q")

  budget <- .parse_sampled_budget(nb_samples_grad, ns, nt)
  nb_p <- budget$nb_p
  nb_q <- budget$nb_q

  max_iter <- .validate_count(max_iter, "max_iter")
  sinkhorn_max_iter <- .validate_count(sinkhorn_max_iter, "sinkhorn_max_iter")

  if (!is.null(random_state)) set.seed(as.integer(random_state))

  T <- p %o% q
  it_last <- 0L
  continue_small <- 0L
  lp_solver <- match.arg(lp_solver)
  lp_direction <- if (epsilon <= 0) .make_transport_lp_solver(p, q, scale = lp_scale, solver = lp_solver) else NULL

  metric_code <- if (identical(metric, "euclidean")) 0L else 1L
  use_cpp <- isTRUE(use_cpp)

  if (epsilon > 0 && use_cpp &&
      exists("cpp_sampled_gromov_wasserstein_coords_entropic_square", mode = "function")) {
    use_mixed_precision <- .runtime_env_bool("RFUGW_SAMPLED_MIXED", FALSE)
    out_cpp <- cpp_sampled_gromov_wasserstein_coords_entropic_square(
      X1 = X1,
      X2 = X2,
      p = p,
      q = q,
      nb_p = as.integer(nb_p),
      nb_q = as.integer(nb_q),
      epsilon = epsilon,
      max_iter = as.integer(max_iter),
      sinkhorn_max_iter = as.integer(sinkhorn_max_iter),
      sinkhorn_tol = sinkhorn_tol,
      metric_code = as.integer(metric_code),
      init_plan = T,
      verbose = isTRUE(verbose),
      use_mixed_precision = use_mixed_precision,
      deterministic_sampling = deterministic_sampling
    )
    T <- out_cpp$plan
    it_last <- as.integer(out_cpp$iterations)
  } else {
    X1_norm2 <- rowSums(X1 * X1)
    X2_norm2 <- rowSums(X2 * X2)
    ones_ns <- rep(1, ns)
    ones_nt <- rep(1, nt)

    for (it in seq_len(max_iter)) {
      it_last <- it
      idx0 <- if (deterministic_sampling) {
        .deterministic_topk_indices(p, min(nb_p, ns))
      } else {
        sample.int(ns, size = min(nb_p, ns), prob = p, replace = FALSE)
      }
      Lik <- matrix(0, nrow = ns, ncol = nt)

      for (i in idx0) {
        row_prob <- T[i, ]
        srow <- sum(row_prob)
        if (!is.finite(srow) || srow <= 0) {
          row_prob <- q
        } else {
          row_prob <- row_prob / srow
        }
        kq <- min(nb_q, nt)
        idx1 <- if (deterministic_sampling) {
          .deterministic_topk_indices(row_prob, kq)
        } else {
          replace_q <- sum(row_prob > 0) < kq
          sample.int(nt, size = kq, prob = row_prob, replace = replace_q)
        }

        d1_i <- .row_distance_from_coords(X1, X1_norm2, i, metric = metric)
        mu2 <- numeric(nt)
        mu2_sq <- numeric(nt)
        for (j in idx1) {
          d2_j <- .row_distance_from_coords(X2, X2_norm2, j, metric = metric)
          mu2 <- mu2 + d2_j
          mu2_sq <- mu2_sq + d2_j^2
        }
        denom <- max(1L, length(idx1))
        mu2 <- mu2 / denom
        mu2_sq <- mu2_sq / denom

        block <- tcrossprod(d1_i^2, ones_nt) +
          tcrossprod(ones_ns, mu2_sq) -
          2 * tcrossprod(d1_i, mu2)
        Lik <- Lik + block
      }

      max_lik <- suppressWarnings(max(Lik, na.rm = TRUE))
      if (is.finite(max_lik) && max_lik > 0) {
        Lik <- Lik / max_lik
      }

      if (epsilon > 0) {
        logT <- log(pmax(T, exp(-200)))
        logT[logT <= -200] <- -Inf
        Lik_eff <- Lik - epsilon * logT
        new_T <- .sinkhorn_balanced(
          a = p,
          b = q,
          M = Lik_eff,
          reg = epsilon,
          max_iter = sinkhorn_max_iter,
          tol = sinkhorn_tol
        )
      } else {
        new_T <- lp_direction(Lik)
      }

      change_T <- mean((T - new_T)^2)
      if (!is.finite(change_T)) {
        break
      }

      if (change_T <= 1e-19) {
        continue_small <- continue_small + 1L
        if (continue_small > 100L) {
          T <- new_T
          break
        }
      } else {
        continue_small <- 0L
      }

      if (isTRUE(verbose) && (it %% 10L == 0L || it == 1L)) {
        cat(sprintf("iter=%d change=%.8e\n", it, change_T))
      }
      T <- new_T
    }
  }

  if (!isTRUE(log)) {
    return(T)
  }

  gw_est <- .estimate_gw_square_from_coords(
    X1 = X1,
    X2 = X2,
    T = T,
    metric = metric,
    n_samples = objective_samples
  )
  list(
    plan = T,
    gw_dist_estimated = gw_est,
    iterations = as.integer(it_last)
  )
}

#' Sampled GW from Sparse or Dense Similarity Graphs
#'
#' Convenience wrapper that converts similarity graphs into diffusion
#' coordinates and solves sampled GW in coordinate space.
#'
#' @param W1 Source similarity graph (`n1 x n1`) as dense matrix or sparse
#'   `Matrix::sparseMatrix`.
#' @param W2 Target similarity graph (`n2 x n2`) as dense matrix or sparse
#'   `Matrix::sparseMatrix`.
#' @param n_components Number of diffusion coordinates.
#' @param diffusion_time Diffusion-time exponent.
#' @param symmetrize If `TRUE`, symmetrize each graph.
#' @param self_loop Optional diagonal regularization added before diffusion map.
#' @param tol Degree floor for graph validity checks.
#' @param p Source weights (default uniform).
#' @param q Target weights (default uniform).
#' @param metric Ground metric for coordinate GW.
#' @param loss_fun Currently only `"square_loss"` is supported.
#' @param nb_samples_grad Number of sampled gradient points, or length-2 vector
#'   `(n_source_samples, n_target_samples)`. Values below 1 error. Source or
#'   target counts above `ns` / `nt` warn and clamp.
#' @param epsilon Entropic regularization.
#' @param max_iter Maximum stochastic iterations.
#' @param log If `TRUE`, return diagnostics.
#' @param verbose If `TRUE`, print iterative diagnostics.
#' @param random_state Optional seed.
#' @param sinkhorn_max_iter Sinkhorn iterations when `epsilon > 0`.
#' @param sinkhorn_tol Sinkhorn tolerance when `epsilon > 0`.
#' @param lp_solver LP backend when `epsilon <= 0`.
#' @param lp_scale Integer scaling for LP marginals.
#' @param use_cpp Use C++ fast path for entropic coordinate solver when
#'   available.
#' @param objective_samples Monte Carlo samples for objective estimation.
#' @param return_embeddings If `TRUE`, include diffusion embeddings in output.
#' @param sampling Sampling policy passed to
#'   [sampled_gromov_wasserstein_coords()].
#' @return Coupling matrix or diagnostics list, mirroring
#'   `sampled_gromov_wasserstein_coords`.
#'
#' @section Experimental:
#' Graph-then-sampled GW is experimental. The certified 0.1 envelope is
#' the same quality-versus-budget claim as
#' [sampled_gromov_wasserstein_coords()], plus that a sparse similarity
#' graph plus `k` diffusion coordinates stores less than two dense
#' `n x n` structure costs. See `inst/bench/sampled-budget-curves.md`.
#' @export
sampled_gw_from_graphs <- function(
    W1,
    W2,
    n_components = 16L,
    diffusion_time = 1,
    symmetrize = TRUE,
    self_loop = 1e-6,
    tol = 1e-12,
    p = NULL,
    q = NULL,
    metric = c("euclidean", "sqeuclidean"),
    loss_fun = "square_loss",
    nb_samples_grad = 100L,
    epsilon = 1,
    max_iter = 500L,
    log = FALSE,
    verbose = FALSE,
    random_state = NULL,
    sinkhorn_max_iter = 200L,
    sinkhorn_tol = 1e-9,
    lp_solver = c("lp_matrix", "lp_transport"),
    lp_scale = 1e6,
    use_cpp = TRUE,
    objective_samples = 2000L,
    return_embeddings = FALSE,
    sampling = c("stochastic", "deterministic")) {
  metric <- match.arg(metric)
  sampling <- match.arg(sampling)
  E1 <- graph_diffusion_coordinates(
    W = W1,
    n_components = as.integer(n_components),
    diffusion_time = diffusion_time,
    symmetrize = symmetrize,
    self_loop = self_loop,
    tol = tol,
    log = FALSE
  )
  E2 <- graph_diffusion_coordinates(
    W = W2,
    n_components = as.integer(n_components),
    diffusion_time = diffusion_time,
    symmetrize = symmetrize,
    self_loop = self_loop,
    tol = tol,
    log = FALSE
  )

  out <- sampled_gromov_wasserstein_coords(
    X1 = E1,
    X2 = E2,
    p = p,
    q = q,
    metric = metric,
    loss_fun = loss_fun,
    nb_samples_grad = nb_samples_grad,
    epsilon = epsilon,
    max_iter = max_iter,
    log = log,
    verbose = verbose,
    random_state = random_state,
    sinkhorn_max_iter = sinkhorn_max_iter,
    sinkhorn_tol = sinkhorn_tol,
    lp_solver = lp_solver,
    lp_scale = lp_scale,
    use_cpp = use_cpp,
    objective_samples = objective_samples,
    sampling = sampling
  )

  if (!isTRUE(return_embeddings)) {
    return(out)
  }
  if (!isTRUE(log)) {
    return(list(plan = out, X1_embed = E1, X2_embed = E2))
  }
  out$X1_embed <- E1
  out$X2_embed <- E2
  out
}
