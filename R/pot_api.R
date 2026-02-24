.check_square_loss <- function(loss_fun) {
  loss_fun <- as.character(loss_fun)[1]
  if (!identical(loss_fun, "square_loss")) {
    stop("Only `loss_fun = \"square_loss\"` is currently supported.", call. = FALSE)
  }
  invisible(TRUE)
}

.is_symmetric_cost <- function(C, tol = 1e-10) {
  max(abs(C - t(C))) <= tol
}

.init_semirelaxed_square <- function(C1, C2, p) {
  ns <- nrow(C1)
  nt <- nrow(C2)
  constC <- tcrossprod(as.vector((C1^2) %*% p), rep(1, nt))
  list(
    constC = constC,
    hC1 = C1,
    hC2 = 2 * C2,
    fC2t = t(C2^2),
    ones_p = rep(1, ns)
  )
}

.validate_semirelaxed_init <- function(G0, p, ns, nt) {
  if (is.null(G0)) {
    return(tcrossprod(p, rep(1 / nt, nt)))
  }
  .assert_matrix(G0, "G0")
  if (nrow(G0) != ns || ncol(G0) != nt) {
    stop("`G0` must have shape nrow(C1) x nrow(C2).", call. = FALSE)
  }
  if (any(!is.finite(G0)) || any(G0 < 0)) {
    stop("`G0` must be finite and nonnegative.", call. = FALSE)
  }
  if (max(abs(rowSums(G0) - p)) > 1e-8) {
    stop("`G0` must satisfy row marginals equal to `p`.", call. = FALSE)
  }
  G0
}

.rowwise_kl_mirror_update <- function(G, grad, p, epsilon, tiny = 1e-300) {
  Z <- log(pmax(G, tiny)) - grad / epsilon
  row_max <- apply(Z, 1L, max)
  W <- exp(Z - row_max)
  rs <- rowSums(W)
  rs[!is.finite(rs) | rs <= 0] <- 1
  sweep(W, 1L, p / rs, `*`)
}

.entropic_semirelaxed_fgw_core <- function(
    M,
    C1,
    C2,
    p = NULL,
    epsilon = 0.1,
    alpha = 0.5,
    symmetric = NULL,
    G0 = NULL,
    max_iter = 10000L,
    tol = 1e-9,
    check_every = 10L,
    verbose = FALSE) {
  .assert_matrix(M, "M")
  .assert_matrix(C1, "C1")
  .assert_matrix(C2, "C2")

  ns <- nrow(C1)
  nt <- nrow(C2)
  if (ncol(C1) != ns) stop("`C1` must be square.", call. = FALSE)
  if (ncol(C2) != nt) stop("`C2` must be square.", call. = FALSE)
  if (nrow(M) != ns || ncol(M) != nt) {
    stop("`M` must have shape nrow(C1) x nrow(C2).", call. = FALSE)
  }
  if (!is.finite(alpha) || alpha < 0 || alpha > 1) {
    stop("`alpha` must be in [0, 1].", call. = FALSE)
  }
  if (!is.finite(epsilon) || epsilon <= 0) {
    stop("`epsilon` must be positive.", call. = FALSE)
  }
  max_iter <- as.integer(max_iter)[1]
  if (!is.finite(max_iter) || max_iter < 1L) {
    stop("`max_iter` must be an integer >= 1.", call. = FALSE)
  }
  check_every <- as.integer(check_every)[1]
  if (!is.finite(check_every) || check_every < 1L) {
    stop("`check_every` must be an integer >= 1.", call. = FALSE)
  }

  if (is.null(p)) p <- rep(1 / ns, ns)
  p <- .assert_prob(p, ns, "p")
  if (is.null(symmetric)) {
    symmetric <- .is_symmetric_cost(C1) && .is_symmetric_cost(C2)
  } else {
    symmetric <- isTRUE(symmetric)
  }
  G <- .validate_semirelaxed_init(G0, p, ns, nt)

  core <- .init_semirelaxed_square(C1, C2, p)
  constC <- core$constC
  hC1 <- core$hC1
  hC2 <- core$hC2
  fC2t <- core$fC2t
  ones_p <- core$ones_p

  if (!symmetric) {
    core_t <- .init_semirelaxed_square(t(C1), t(C2), p)
    constCt <- core_t$constC
    hC1t <- core_t$hC1
    hC2t <- core_t$hC2
    fC2 <- core_t$fC2t
  }

  err <- Inf
  it <- 0L
  for (k in seq_len(max_iter)) {
    Gprev <- G
    qG <- colSums(G)

    if (symmetric) {
      marginal <- tcrossprod(ones_p, as.vector(qG %*% fC2t))
      grad_quad <- 2 * (constC + marginal - hC1 %*% G %*% t(hC2))
    } else {
      marginal_1 <- tcrossprod(ones_p, as.vector(qG %*% fC2t))
      grad_1 <- 2 * (constC + marginal_1 - hC1 %*% G %*% t(hC2))
      marginal_2 <- tcrossprod(ones_p, as.vector(qG %*% fC2))
      grad_2 <- 2 * (constCt + marginal_2 - hC1t %*% G %*% t(hC2t))
      grad_quad <- 0.5 * (grad_1 + grad_2)
    }

    grad <- alpha * grad_quad + (1 - alpha) * M
    G <- .rowwise_kl_mirror_update(G = G, grad = grad, p = p, epsilon = epsilon)

    do_check <- (k %% check_every) == 0L || k == 1L
    if (do_check) {
      err <- sqrt(sum((G - Gprev)^2))
      if (isTRUE(verbose)) {
        if ((k - 1L) %% 200L == 0L) {
          cat(sprintf("%5s|%12s\n", "It.", "Err"), "-------------------\n", sep = "")
        }
        cat(sprintf("%5d|%8e|\n", k - 1L, err))
      }
      if (is.finite(err) && err <= tol) {
        it <- k
        break
      }
    }
    it <- k
  }

  qG <- colSums(G)
  marginal <- tcrossprod(ones_p, as.vector(qG %*% fC2t))
  quad_raw <- sum((constC + marginal - hC1 %*% G %*% t(hC2)) * G)
  lin_loss <- (1 - alpha) * sum(M * G)
  quad_loss <- alpha * quad_raw
  srfgw_dist <- lin_loss + quad_loss

  list(
    plan = G,
    q = qG,
    lin_loss = lin_loss,
    quad_loss = quad_loss,
    srfgw_dist = srfgw_dist,
    srgw_dist = quad_raw,
    iterations = as.integer(it),
    error = as.numeric(err),
    symmetric = symmetric
  )
}

#' Entropic Gromov-Wasserstein (square loss)
#'
#' POT-compatible GW wrapper using the optimized FGW entropic core with a zero
#' feature cost and `alpha = 1`.
#'
#' @param C1 Source structure cost matrix (`ns x ns`).
#' @param C2 Target structure cost matrix (`nt x nt`).
#' @param p Source weights (default uniform).
#' @param q Target weights (default uniform).
#' @param loss_fun Inner loss, currently only `"square_loss"` is supported.
#' @param epsilon Entropic regularization.
#' @param symmetric Assume `C1` and `C2` are symmetric if `TRUE`.
#' @param G0 Optional initial coupling (`ns x nt`) warm start.
#' @param max_iter Max outer GW iterations.
#' @param tol Outer stopping tolerance on plan updates.
#' @param solver Either `"PGD"` or `"PPA"`.
#' @param sinkhorn_max_iter Max Sinkhorn iterations per outer step.
#' @param sinkhorn_tol Sinkhorn stopping tolerance.
#' @param sinkhorn_method Sinkhorn variant (`"scaling"` or `"log"`).
#' @param precision Numeric precision mode (`"mixed"` or `"double"`).
#' @param check_every Evaluate stopping criterion every `check_every` iterations.
#' @param structure_rank Optional low-rank rank for structure tensor product.
#' @return A list with `plan`, `gw_dist`, `iterations`, `error`.
#' @export
entropic_gromov_wasserstein <- function(
    C1,
    C2,
    p = NULL,
    q = NULL,
    loss_fun = "square_loss",
    epsilon = 0.1,
    symmetric = TRUE,
    G0 = NULL,
    max_iter = 1000L,
    tol = 1e-9,
    solver = c("PGD", "PPA"),
    sinkhorn_max_iter = 1000L,
    sinkhorn_tol = 1e-9,
    sinkhorn_method = c("scaling", "log"),
    precision = c("mixed", "double"),
    check_every = 1L,
    structure_rank = 0L) {
  .check_square_loss(loss_fun)
  .assert_matrix(C1, "C1")
  .assert_matrix(C2, "C2")
  ns <- nrow(C1)
  nt <- nrow(C2)
  M0 <- matrix(0, nrow = ns, ncol = nt)
  out <- fgw_entropic(
    M = M0,
    C1 = C1,
    C2 = C2,
    p = p,
    q = q,
    alpha = 1,
    epsilon = epsilon,
    max_iter = max_iter,
    tol = tol,
    sinkhorn_max_iter = sinkhorn_max_iter,
    sinkhorn_tol = sinkhorn_tol,
    init_plan = G0,
    structure_rank = structure_rank,
    sinkhorn_method = sinkhorn_method,
    precision = precision,
    symmetric = symmetric,
    solver = solver,
    check_every = check_every
  )
  out$gw_dist <- out$fgw_dist
  out
}

#' Entropic Gromov-Wasserstein Objective Value
#'
#' @inheritParams entropic_gromov_wasserstein
#' @return Numeric scalar GW value.
#' @export
entropic_gromov_wasserstein2 <- function(...) {
  out <- entropic_gromov_wasserstein(...)
  out$gw_dist
}

#' Exact Gromov-Wasserstein via Conditional Gradient (square loss)
#'
#' POT-compatible exact GW wrapper using the optimized exact FGW core with
#' zero feature cost and `alpha = 1`.
#'
#' @inheritParams entropic_gromov_wasserstein
#' @param tol_rel Relative stopping tolerance for exact CG.
#' @param tol_abs Absolute stopping tolerance for exact CG.
#' @param lp_solver LP backend (`"cpp_transport"`, `"lp_transport"`, `"lp_matrix"`).
#' @param lp_max_iter Maximum iterations for C++ transport simplex LP backend.
#' @param lp_tol Optimality tolerance for C++ transport simplex LP backend.
#' @return A list with `plan`, `gw_dist`, `iterations`, `error`, `loss_trace`.
#' @export
gromov_wasserstein <- function(
    C1,
    C2,
    p = NULL,
    q = NULL,
    loss_fun = "square_loss",
    symmetric = TRUE,
    G0 = NULL,
    max_iter = 500L,
    tol_rel = 1e-9,
    tol_abs = 1e-9,
    lp_solver = c("cpp_transport", "lp_transport", "lp_matrix"),
    lp_max_iter = 20000L,
    lp_tol = 1e-12) {
  .check_square_loss(loss_fun)
  .assert_matrix(C1, "C1")
  .assert_matrix(C2, "C2")
  ns <- nrow(C1)
  nt <- nrow(C2)
  M0 <- matrix(0, nrow = ns, ncol = nt)
  out <- fgw_exact_cg(
    M = M0,
    C1 = C1,
    C2 = C2,
    p = p,
    q = q,
    alpha = 1,
    symmetric = symmetric,
    G0 = G0,
    max_iter = max_iter,
    tol_rel = tol_rel,
    tol_abs = tol_abs,
    lp_solver = lp_solver,
    lp_max_iter = lp_max_iter,
    lp_tol = lp_tol
  )
  out$gw_dist <- out$fgw_dist
  out
}

#' Exact Gromov-Wasserstein Objective Value
#'
#' @inheritParams gromov_wasserstein
#' @return Numeric scalar GW value.
#' @export
gromov_wasserstein2 <- function(...) {
  out <- gromov_wasserstein(...)
  out$gw_dist
}

#' POT Alias for Entropic FGW
#'
#' @inheritParams fgw_entropic
#' @return A list with `plan`, `fgw_dist`, `iterations`, `error`.
#' @export
entropic_fused_gromov_wasserstein <- function(...) {
  fgw_entropic(...)
}

#' POT Alias for Entropic FGW Value
#'
#' @inheritParams fgw_entropic
#' @return Numeric scalar FGW value.
#' @export
entropic_fused_gromov_wasserstein2 <- function(...) {
  fgw_entropic2(...)
}

#' POT Alias for Exact FGW
#'
#' @inheritParams fgw_exact_cg
#' @return A list with `plan`, `fgw_dist`, `iterations`, `error`, `loss_trace`.
#' @export
fused_gromov_wasserstein <- function(...) {
  fgw_exact_cg(...)
}

#' POT Alias for Exact FGW Value
#'
#' @inheritParams fgw_exact_cg
#' @return Numeric scalar FGW value.
#' @export
fused_gromov_wasserstein2 <- function(...) {
  out <- fused_gromov_wasserstein(...)
  out$fgw_dist
}

#' POT Alias for Fused Unbalanced GW
#'
#' @inheritParams fugw_kl
#' @return A list with `pi_samp`, `pi_feat`, `fugw_cost`, `linear_cost`,
#'   `iterations`, and `error`.
#' @export
fused_unbalanced_gromov_wasserstein <- function(...) {
  fugw_kl(...)
}

#' POT Alias for Fused Unbalanced GW Value
#'
#' @inheritParams fugw_kl
#' @return Numeric scalar FUGW value.
#' @export
fused_unbalanced_gromov_wasserstein2 <- function(...) {
  fugw_kl2(...)
}

#' Entropic Semirelaxed Gromov-Wasserstein (square loss)
#'
#' Mirror-descent solver in KL geometry with a fixed source marginal and free
#' target marginal.
#'
#' @param C1 Source structure cost matrix (`ns x ns`).
#' @param C2 Target structure cost matrix (`nt x nt`).
#' @param p Source weights (default uniform).
#' @param loss_fun Inner loss, currently only `"square_loss"` is supported.
#' @param epsilon Entropic regularization.
#' @param symmetric If `NULL`, inferred from `C1` and `C2`; otherwise logical.
#' @param G0 Optional initial semirelaxed plan (`ns x nt`) satisfying
#'   `rowSums(G0) == p`.
#' @param max_iter Max mirror-descent iterations.
#' @param tol Stopping tolerance on Frobenius norm of plan updates.
#' @param check_every Evaluate stopping criterion every `check_every` iterations.
#' @param verbose If `TRUE`, print iterative error diagnostics.
#' @return A list with `plan`, `q`, `srgw_dist`, `iterations`, and `error`.
#' @export
entropic_semirelaxed_gromov_wasserstein <- function(
    C1,
    C2,
    p = NULL,
    loss_fun = "square_loss",
    epsilon = 0.1,
    symmetric = NULL,
    G0 = NULL,
    max_iter = 10000L,
    tol = 1e-9,
    check_every = 10L,
    verbose = FALSE,
    precision = c("mixed", "double"),
    backend = c("cpp", "r")) {
  .check_square_loss(loss_fun)
  precision <- match.arg(precision)
  backend <- match.arg(backend)
  .assert_matrix(C1, "C1")
  .assert_matrix(C2, "C2")
  if (is.null(symmetric)) {
    symmetric <- .is_symmetric_cost(C1) && .is_symmetric_cost(C2)
  } else {
    symmetric <- isTRUE(symmetric)
  }
  ns <- nrow(C1)
  nt <- nrow(C2)
  M0 <- matrix(0, nrow = ns, ncol = nt)
  out <- if (backend == "cpp") {
    if (is.null(p)) p <- rep(1 / ns, ns)
    p <- .assert_prob(p, ns, "p")
    if (is.null(G0)) {
      G0 <- matrix(numeric(0), nrow = 0, ncol = 0)
    } else {
      .assert_matrix(G0, "G0")
    }
    cpp_entropic_semirelaxed_fgw_square(
      M = M0,
      C1 = C1,
      C2 = C2,
      p = p,
      epsilon = epsilon,
      alpha = 1,
      symmetric = symmetric,
      use_mixed_precision = identical(precision, "mixed"),
      max_iter = as.integer(max_iter),
      tol = tol,
      check_every = as.integer(check_every),
      init_plan = G0
    )
  } else {
    .entropic_semirelaxed_fgw_core(
      M = M0,
      C1 = C1,
      C2 = C2,
      p = p,
      epsilon = epsilon,
      alpha = 1,
      symmetric = symmetric,
      G0 = G0,
      max_iter = max_iter,
      tol = tol,
      check_every = check_every,
      verbose = verbose
    )
  }
  list(
    plan = out$plan,
    q = out$q,
    srgw_dist = out$srgw_dist,
    iterations = out$iterations,
    error = out$error,
    symmetric = out$symmetric
  )
}

#' Entropic Semirelaxed Gromov-Wasserstein Objective Value
#'
#' @inheritParams entropic_semirelaxed_gromov_wasserstein
#' @return Numeric scalar semirelaxed GW value.
#' @export
entropic_semirelaxed_gromov_wasserstein2 <- function(...) {
  out <- entropic_semirelaxed_gromov_wasserstein(...)
  out$srgw_dist
}

#' Entropic Semirelaxed Fused Gromov-Wasserstein (square loss)
#'
#' Mirror-descent solver in KL geometry with fixed source marginal and free
#' target marginal.
#'
#' @param M Cross-domain feature cost matrix (`ns x nt`).
#' @param C1 Source structure cost matrix (`ns x ns`).
#' @param C2 Target structure cost matrix (`nt x nt`).
#' @param p Source weights (default uniform).
#' @param loss_fun Inner loss, currently only `"square_loss"` is supported.
#' @param symmetric If `NULL`, inferred from `C1` and `C2`; otherwise logical.
#' @param epsilon Entropic regularization.
#' @param alpha Feature/structure trade-off in `[0, 1]`.
#' @param G0 Optional initial semirelaxed plan (`ns x nt`) with row marginals `p`.
#' @param max_iter Max mirror-descent iterations.
#' @param tol Stopping tolerance on Frobenius norm of plan updates.
#' @param check_every Evaluate stopping criterion every `check_every` iterations.
#' @param verbose If `TRUE`, print iterative error diagnostics.
#' @return A list with `plan`, `q`, `srfgw_dist`, `lin_loss`, `quad_loss`,
#'   `iterations`, and `error`.
#' @export
entropic_semirelaxed_fused_gromov_wasserstein <- function(
    M,
    C1,
    C2,
    p = NULL,
    loss_fun = "square_loss",
    symmetric = NULL,
    epsilon = 0.1,
    alpha = 0.5,
    G0 = NULL,
    max_iter = 10000L,
    tol = 1e-9,
    check_every = 10L,
    verbose = FALSE,
    precision = c("mixed", "double"),
    backend = c("cpp", "r")) {
  .check_square_loss(loss_fun)
  precision <- match.arg(precision)
  backend <- match.arg(backend)
  .assert_matrix(M, "M")
  .assert_matrix(C1, "C1")
  .assert_matrix(C2, "C2")
  if (is.null(symmetric)) {
    symmetric <- .is_symmetric_cost(C1) && .is_symmetric_cost(C2)
  } else {
    symmetric <- isTRUE(symmetric)
  }
  ns <- nrow(C1)
  out <- if (backend == "cpp") {
    if (is.null(p)) p <- rep(1 / ns, ns)
    p <- .assert_prob(p, ns, "p")
    if (is.null(G0)) {
      G0 <- matrix(numeric(0), nrow = 0, ncol = 0)
    } else {
      .assert_matrix(G0, "G0")
    }
    cpp_entropic_semirelaxed_fgw_square(
      M = M,
      C1 = C1,
      C2 = C2,
      p = p,
      epsilon = epsilon,
      alpha = alpha,
      symmetric = symmetric,
      use_mixed_precision = identical(precision, "mixed"),
      max_iter = as.integer(max_iter),
      tol = tol,
      check_every = as.integer(check_every),
      init_plan = G0
    )
  } else {
    .entropic_semirelaxed_fgw_core(
      M = M,
      C1 = C1,
      C2 = C2,
      p = p,
      epsilon = epsilon,
      alpha = alpha,
      symmetric = symmetric,
      G0 = G0,
      max_iter = max_iter,
      tol = tol,
      check_every = check_every,
      verbose = verbose
    )
  }
  out
}

#' Entropic Semirelaxed Fused Gromov-Wasserstein Objective Value
#'
#' @inheritParams entropic_semirelaxed_fused_gromov_wasserstein
#' @return Numeric scalar semirelaxed FGW value.
#' @export
entropic_semirelaxed_fused_gromov_wasserstein2 <- function(...) {
  out <- entropic_semirelaxed_fused_gromov_wasserstein(...)
  out$srfgw_dist
}

#' Fixed-Support FGW Barycenters
#'
#' Learns a fixed-size barycenter support by alternating between (entropic) FGW
#' couplings and closed-form barycenter updates of structure and features.
#'
#' @param N Number of barycenter nodes.
#' @param Ys List of feature matrices, one per sample (`ns x d`).
#' @param Cs List of structure matrices, one per sample (`ns x ns`).
#' @param ps Optional list of sample weights for each sample.
#' @param lambdas Optional sample weights across sets (default uniform).
#' @param alpha FGW feature/structure trade-off.
#' @param fixed_structure If `TRUE`, keep barycenter structure fixed at `init_C`.
#' @param fixed_features If `TRUE`, keep barycenter features fixed at `init_X`.
#' @param p Optional barycenter weights (`N`) (default uniform).
#' @param loss_fun Inner loss, currently only `"square_loss"` is supported.
#' @param epsilon Entropic regularization used for inner FGW solves.
#' @param symmetric Assume symmetric structure matrices in inner FGW solves.
#' @param max_iter Maximum outer barycenter iterations.
#' @param tol Stopping tolerance on barycenter updates.
#' @param warmstartT If `TRUE`, warm-start inner FGW with previous couplings.
#' @param solver Inner FGW solver (`"PGD"` or `"PPA"`).
#' @param sinkhorn_max_iter Inner Sinkhorn iterations for FGW solves.
#' @param sinkhorn_tol Inner Sinkhorn tolerance for FGW solves.
#' @param sinkhorn_method Sinkhorn variant (`"scaling"` or `"log"`).
#' @param precision Numeric precision mode for inner FGW solves.
#' @param check_every Outer FGW check interval for inner solves.
#' @param init_C Optional initial barycenter structure (`N x N`).
#' @param init_X Optional initial barycenter feature matrix (`N x d`).
#' @param random_state Optional seed used only when `init_C` is `NULL`.
#' @param feature_cost_metric Feature cost metric (`"euclidean"` or
#'   `"sqeuclidean"`).
#' @param feature_cost_normalize If `TRUE`, normalize feature costs to `[0, 1]`.
#' @param barycenter_eps Stabilizer for divisions in barycenter updates.
#' @param verbose If `TRUE`, print outer-iteration diagnostics.
#' @return A list with `X`, `C`, `p`, `couplings`, `history`, `iterations`,
#'   and `objective`.
#' @export
fgw_barycenters <- function(
    N,
    Ys,
    Cs,
    ps = NULL,
    lambdas = NULL,
    alpha = 0.5,
    fixed_structure = FALSE,
    fixed_features = FALSE,
    p = NULL,
    loss_fun = "square_loss",
    epsilon = 0.05,
    symmetric = TRUE,
    max_iter = 100L,
    tol = 1e-9,
    warmstartT = TRUE,
    solver = c("PGD", "PPA"),
    sinkhorn_max_iter = 500L,
    sinkhorn_tol = 1e-9,
    sinkhorn_method = c("scaling", "log"),
    precision = c("mixed", "double"),
    check_every = 10L,
    init_C = NULL,
    init_X = NULL,
    random_state = 1L,
    feature_cost_metric = c("euclidean", "sqeuclidean"),
    feature_cost_normalize = TRUE,
    barycenter_eps = 1e-12,
    verbose = FALSE) {
  .check_square_loss(loss_fun)
  solver <- match.arg(solver)
  sinkhorn_method <- match.arg(sinkhorn_method)
  precision <- match.arg(precision)
  feature_cost_metric <- match.arg(feature_cost_metric)

  if (!is.list(Cs) || length(Cs) == 0L) {
    stop("`Cs` must be a non-empty list of square matrices.", call. = FALSE)
  }
  if (!is.list(Ys) || length(Ys) != length(Cs)) {
    stop("`Ys` must be a list with the same length as `Cs`.", call. = FALSE)
  }

  S <- length(Cs)
  N <- as.integer(N)[1]
  if (!is.finite(N) || N < 2L) {
    stop("`N` must be an integer >= 2.", call. = FALSE)
  }
  if (!is.finite(alpha) || alpha < 0 || alpha > 1) {
    stop("`alpha` must be in [0, 1].", call. = FALSE)
  }
  if (!is.finite(epsilon) || epsilon <= 0) {
    stop("`epsilon` must be positive.", call. = FALSE)
  }

  Cs <- lapply(seq_along(Cs), function(i) {
    Ci <- Cs[[i]]
    .assert_matrix(Ci, sprintf("Cs[[%d]]", i))
    if (nrow(Ci) != ncol(Ci)) {
      stop(sprintf("`Cs[[%d]]` must be square.", i), call. = FALSE)
    }
    Ci
  })
  ns <- vapply(Cs, nrow, integer(1))

  Ys <- lapply(seq_along(Ys), function(i) {
    Yi <- Ys[[i]]
    .assert_matrix(Yi, sprintf("Ys[[%d]]", i))
    if (nrow(Yi) != ns[[i]]) {
      stop(sprintf("`nrow(Ys[[%d]])` must equal `nrow(Cs[[%d]])`.", i, i), call. = FALSE)
    }
    Yi
  })
  dvals <- unique(vapply(Ys, ncol, integer(1)))
  if (length(dvals) != 1L) {
    stop("All feature matrices in `Ys` must have the same number of columns.", call. = FALSE)
  }
  d <- dvals[[1]]

  if (is.null(ps)) {
    ps <- lapply(ns, function(n) rep(1 / n, n))
  } else {
    if (!is.list(ps) || length(ps) != S) {
      stop("`ps` must be `NULL` or a list with one weight vector per sample.", call. = FALSE)
    }
    ps <- lapply(seq_along(ps), function(i) {
      .assert_prob(ps[[i]], ns[[i]], sprintf("ps[[%d]]", i))
    })
  }

  if (is.null(lambdas)) {
    lambdas <- rep(1 / S, S)
  }
  lambdas <- .assert_prob(lambdas, S, "lambdas")

  if (is.null(p)) {
    p <- rep(1 / N, N)
  }
  p <- .assert_prob(p, N, "p")

  if (isTRUE(fixed_structure)) {
    if (is.null(init_C)) {
      stop("If `fixed_structure = TRUE`, `init_C` must be provided.", call. = FALSE)
    }
    .assert_matrix(init_C, "init_C")
    if (nrow(init_C) != N || ncol(init_C) != N) {
      stop("`init_C` must have shape N x N.", call. = FALSE)
    }
    C <- init_C
  } else if (is.null(init_C)) {
    if (!is.null(random_state)) set.seed(as.integer(random_state))
    x0 <- matrix(stats::rnorm(N * 2L), nrow = N, ncol = 2L)
    C <- as.matrix(stats::dist(x0))
  } else {
    .assert_matrix(init_C, "init_C")
    if (nrow(init_C) != N || ncol(init_C) != N) {
      stop("`init_C` must have shape N x N.", call. = FALSE)
    }
    C <- init_C
  }

  if (isTRUE(fixed_features)) {
    if (is.null(init_X)) {
      stop("If `fixed_features = TRUE`, `init_X` must be provided.", call. = FALSE)
    }
    .assert_matrix(init_X, "init_X")
    if (nrow(init_X) != N || ncol(init_X) != d) {
      stop("`init_X` must have shape N x d.", call. = FALSE)
    }
    X <- init_X
  } else if (is.null(init_X)) {
    X <- matrix(0, nrow = N, ncol = d)
  } else {
    .assert_matrix(init_X, "init_X")
    if (nrow(init_X) != N || ncol(init_X) != d) {
      stop("`init_X` must have shape N x d.", call. = FALSE)
    }
    X <- init_X
  }

  Ms <- lapply(Ys, function(Y) {
    .cross_feature_cost(X, Y, metric = feature_cost_metric, rescale01 = feature_cost_normalize)
  })
  T <- vector("list", S)
  objectives <- rep(NA_real_, S)
  history <- vector("list", as.integer(max_iter))
  outer_done <- 0L

  for (it in seq_len(as.integer(max_iter))) {
    Cprev <- C
    Xprev <- X

    for (s in seq_len(S)) {
      init_plan <- if (isTRUE(warmstartT)) T[[s]] else NULL
      out <- fgw_entropic(
        M = Ms[[s]],
        C1 = C,
        C2 = Cs[[s]],
        p = p,
        q = ps[[s]],
        alpha = alpha,
        epsilon = epsilon,
        max_iter = as.integer(max_iter),
        tol = tol,
        sinkhorn_max_iter = as.integer(sinkhorn_max_iter),
        sinkhorn_tol = sinkhorn_tol,
        init_plan = init_plan,
        structure_rank = 0L,
        sinkhorn_method = sinkhorn_method,
        precision = precision,
        symmetric = symmetric,
        solver = solver,
        check_every = as.integer(check_every)
      )
      T[[s]] <- out$plan
      objectives[[s]] <- out$fgw_dist
    }

    inv_p <- 1 / pmax(p, barycenter_eps)
    if (!isTRUE(fixed_features)) {
      numX <- Reduce(
        `+`,
        lapply(seq_len(S), function(s) lambdas[[s]] * (T[[s]] %*% Ys[[s]]))
      )
      X <- sweep(numX, 1L, inv_p, `*`)
      Ms <- lapply(Ys, function(Y) {
        .cross_feature_cost(X, Y, metric = feature_cost_metric, rescale01 = feature_cost_normalize)
      })
    }

    if (!isTRUE(fixed_structure)) {
      numC <- Reduce(
        `+`,
        lapply(seq_len(S), function(s) lambdas[[s]] * (T[[s]] %*% Cs[[s]] %*% t(T[[s]])))
      )
      C <- numC * tcrossprod(inv_p, inv_p)
      C <- 0.5 * (C + t(C))
      diag(C) <- 0
    }

    err_feature <- if (isTRUE(fixed_features)) 0 else sqrt(sum((X - Xprev)^2))
    err_structure <- if (isTRUE(fixed_structure)) 0 else sqrt(sum((C - Cprev)^2))
    obj_total <- sum(lambdas * objectives)
    history[[it]] <- data.frame(
      iter = it,
      objective = obj_total,
      err_feature = err_feature,
      err_structure = err_structure,
      stringsAsFactors = FALSE
    )
    outer_done <- it

    if (isTRUE(verbose)) {
      cat(sprintf(
        "iter=%d objective=%.8e err_feature=%.3e err_structure=%.3e\n",
        it, obj_total, err_feature, err_structure
      ))
    }
    if (max(err_feature, err_structure) <= tol) {
      break
    }
  }

  hist_df <- do.call(rbind, history[seq_len(outer_done)])
  list(
    X = X,
    C = C,
    p = p,
    couplings = T,
    objectives = objectives,
    objective = sum(lambdas * objectives),
    history = hist_df,
    iterations = as.integer(outer_done),
    error = if (outer_done > 0L) {
      max(hist_df$err_feature[[outer_done]], hist_df$err_structure[[outer_done]])
    } else {
      NA_real_
    }
  )
}
