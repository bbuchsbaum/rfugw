.pairwise_sqeuclidean <- function(X, Y = NULL) {
  .assert_matrix(X, "X")
  if (is.null(Y)) {
    Y <- X
  } else {
    .assert_matrix(Y, "Y")
    if (ncol(X) != ncol(Y)) {
      stop("`X` and `Y` must have the same number of columns.", call. = FALSE)
    }
  }
  x2 <- rowSums(X * X)
  y2 <- rowSums(Y * Y)
  D2 <- outer(x2, y2, "+") - 2 * tcrossprod(X, Y)
  D2[D2 < 0] <- 0
  D2
}

.maybe_rescale01 <- function(x, enabled = TRUE) {
  if (!enabled) return(x)
  mx <- max(x)
  if (!is.finite(mx) || mx <= 0) return(x)
  x / mx
}

.normalize_feature_bundle <- function(feature_list, mode = c("none", "zscore")) {
  mode <- match.arg(mode)
  if (mode == "none") return(feature_list)
  keep <- vapply(feature_list, function(x) !is.null(x), logical(1))
  if (!any(keep)) return(feature_list)

  stacked <- do.call(rbind, feature_list[keep])
  mu <- colMeans(stacked)
  sdv <- apply(stacked, 2, stats::sd)
  sdv[!is.finite(sdv) | sdv == 0] <- 1

  out <- feature_list
  for (i in seq_along(out)) {
    if (!is.null(out[[i]])) {
      out[[i]] <- sweep(out[[i]], 2, mu, "-")
      out[[i]] <- sweep(out[[i]], 2, sdv, "/")
    }
  }
  out
}

.with_single_blas_threads <- function(expr, enable = TRUE) {
  if (!isTRUE(enable)) {
    return(eval.parent(substitute(expr)))
  }
  vars <- c(
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
    "BLIS_NUM_THREADS"
  )
  old <- Sys.getenv(vars, unset = NA_character_)
  on.exit({
    restore <- old[!is.na(old)]
    unset <- names(old)[is.na(old)]
    if (length(restore) > 0L) {
      do.call(Sys.setenv, as.list(restore))
    }
    if (length(unset) > 0L) {
      Sys.unsetenv(unset)
    }
  }, add = TRUE)
  do.call(Sys.setenv, as.list(setNames(rep("1", length(vars)), vars)))
  eval.parent(substitute(expr))
}

.validate_set <- function(x, name = "set") {
  if (!is.list(x)) {
    stop(sprintf("`%s` must be a list.", name), call. = FALSE)
  }
  if (is.null(x$C)) {
    stop(sprintf("`%s$C` is required.", name), call. = FALSE)
  }
  .assert_matrix(x$C, sprintf("%s$C", name))
  n <- nrow(x$C)
  if (ncol(x$C) != n) {
    stop(sprintf("`%s$C` must be square.", name), call. = FALSE)
  }

  F <- x$F
  if (!is.null(F)) {
    .assert_matrix(F, sprintf("%s$F", name))
    if (nrow(F) != n) {
      stop(sprintf("`%s$F` must have nrow = nrow(%s$C).", name, name), call. = FALSE)
    }
  }

  w <- x$w
  if (is.null(w)) {
    w <- rep(1 / n, n)
  }
  w <- .assert_prob(w, n, sprintf("%s$w", name))

  id <- x$id
  if (is.null(id)) {
    id <- name
  } else {
    id <- as.character(id)[1]
  }

  list(
    C = x$C,
    F = F,
    w = w,
    id = id
  )
}

.cross_feature_cost <- function(F1, F2, metric = c("euclidean", "sqeuclidean"), rescale01 = TRUE) {
  metric <- match.arg(metric)
  if (is.null(F1) || is.null(F2)) {
    return(matrix(0, nrow = if (is.null(F1)) 0L else nrow(F1), ncol = if (is.null(F2)) 0L else nrow(F2)))
  }
  D2 <- .pairwise_sqeuclidean(F1, F2)
  out <- if (metric == "euclidean") sqrt(D2) else D2
  .maybe_rescale01(out, enabled = rescale01)
}

.knn_sparsify_structure <- function(C, k, fill = NULL) {
  .assert_matrix(C, "C")
  if (nrow(C) != ncol(C)) {
    stop("`C` must be square for kNN sparsification.", call. = FALSE)
  }
  k <- as.integer(k)[1]
  if (!is.finite(k) || k < 1L) {
    return(C)
  }
  n <- nrow(C)
  keep <- min(n, k + 1L)
  if (keep >= n) {
    return(C)
  }
  if (is.null(fill)) {
    fill <- max(C[is.finite(C)], na.rm = TRUE)
    if (!is.finite(fill)) {
      fill <- 1
    }
  }
  out <- matrix(fill, nrow = n, ncol = n, dimnames = dimnames(C))
  for (i in seq_len(n)) {
    idx <- order(C[i, ])[seq_len(keep)]
    out[i, idx] <- C[i, idx]
  }
  out <- pmin(out, t(out))
  dimnames(out) <- dimnames(C)
  diag(out) <- diag(C)
  out
}

.structure_cost_to_similarity <- function(C, sigma = NULL) {
  .assert_matrix(C, "C")
  if (nrow(C) != ncol(C)) {
    stop("`C` must be square.", call. = FALSE)
  }
  if (is.null(sigma)) {
    vals <- C[is.finite(C) & C > 0]
    sigma <- if (length(vals) > 0L) stats::median(vals) else 1
  }
  sigma <- as.numeric(sigma)[1]
  if (!is.finite(sigma) || sigma <= 0) {
    sigma <- 1
  }
  W <- exp(-(C^2) / (2 * sigma^2))
  W <- 0.5 * (W + t(W))
  diag(W) <- 0
  W
}

.subject_ids_unique <- function(ids) {
  if (length(ids) == 0L) return(ids)
  if (length(unique(ids)) == length(ids)) return(ids)
  make.unique(ids)
}

.autotune_multialign_controls <- function(
    ns_source,
    n_template,
    precision,
    max_iter,
    sinkhorn_max_iter,
    check_every,
    level = c("moderate", "aggressive")) {
  level <- match.arg(level)
  n_ref <- as.integer(stats::quantile(c(ns_source, n_template), probs = 0.75, type = 1))
  n_ref <- max(1L, n_ref)
  check_mode <- "max"

  if (identical(level, "aggressive")) {
    check_target <- if (identical(precision, "mixed")) {
      check_mode <- "min"
      if (n_ref >= 2000L) 6L else if (n_ref >= 1200L) 5L else if (n_ref >= 800L) 4L else if (n_ref >= 400L) 4L else 3L
    } else {
      if (n_ref >= 2000L) 12L else if (n_ref >= 1200L) 10L else if (n_ref >= 800L) 8L else if (n_ref >= 400L) 6L else 4L
    }
    sink_cap <- if (identical(precision, "mixed")) {
      if (n_ref >= 2000L) 120L else if (n_ref >= 1200L) 160L else if (n_ref >= 800L) 200L else if (n_ref >= 400L) 240L else 320L
    } else {
      if (n_ref >= 2000L) 160L else if (n_ref >= 1200L) 200L else if (n_ref >= 800L) 240L else if (n_ref >= 400L) 300L else 380L
    }
    outer_cap <- if (n_ref >= 2000L) 80L else if (n_ref >= 1200L) 100L else if (n_ref >= 800L) 120L else as.integer(max_iter)
  } else {
    check_target <- if (identical(precision, "mixed")) {
      check_mode <- "min"
      if (n_ref >= 2000L) 8L else if (n_ref >= 1200L) 6L else if (n_ref >= 800L) 5L else if (n_ref >= 400L) 4L else 3L
    } else {
      if (n_ref >= 2000L) 10L else if (n_ref >= 1200L) 8L else if (n_ref >= 800L) 6L else if (n_ref >= 400L) 4L else 2L
    }
    sink_cap <- if (identical(precision, "mixed")) {
      if (n_ref >= 2000L) 160L else if (n_ref >= 1200L) 200L else if (n_ref >= 800L) 240L else if (n_ref >= 400L) 280L else 360L
    } else {
      if (n_ref >= 2000L) 200L else if (n_ref >= 1200L) 240L else if (n_ref >= 800L) 300L else if (n_ref >= 400L) 360L else 450L
    }
    outer_cap <- if (n_ref >= 2000L) 100L else if (n_ref >= 1200L) 140L else if (n_ref >= 800L) 180L else as.integer(max_iter)
  }
  check_every_eff <- if (identical(check_mode, "min")) {
    min(as.integer(check_every), as.integer(check_target))
  } else {
    max(as.integer(check_every), as.integer(check_target))
  }
  check_every_eff <- max(1L, as.integer(check_every_eff))

  list(
    n_ref = n_ref,
    check_every = check_every_eff,
    sinkhorn_max_iter = min(as.integer(sinkhorn_max_iter), as.integer(sink_cap)),
    max_iter = min(as.integer(max_iter), as.integer(outer_cap))
  )
}

.resolve_structure_rank_request <- function(structure_rank) {
  if (is.character(structure_rank)) {
    tag <- tolower(trimws(structure_rank[[1]]))
    if (identical(tag, "auto")) {
      return(list(rank = 0L, auto = TRUE, requested = "auto"))
    }
    if (tag %in% c("off", "none", "disabled", "disable", "0")) {
      return(list(rank = 0L, auto = FALSE, requested = structure_rank[[1]]))
    }
    stop("`structure_rank` must be a nonnegative integer or one of {'auto', 'off'}.", call. = FALSE)
  }
  rank <- as.integer(structure_rank)[1]
  if (!is.finite(rank) || rank < 0) {
    stop("`structure_rank` must be a nonnegative integer or one of {'auto', 'off'}.", call. = FALSE)
  }
  list(rank = rank, auto = FALSE, requested = rank)
}

.autotune_structure_rank <- function(
    ns_source,
    n_template,
    precision,
    level = c("moderate", "aggressive")) {
  level <- match.arg(level)
  n_ref <- as.integer(stats::quantile(c(ns_source, n_template), probs = 0.75, type = 1))
  n_ref <- max(1L, n_ref)

  if (identical(level, "aggressive")) {
    rank <- if (n_ref >= 2400L) 160L else if (n_ref >= 1800L) 128L else if (n_ref >= 1200L) 96L else if (n_ref >= 800L) 64L else if (n_ref >= 500L) 48L else if (n_ref >= 350L) 32L else 0L
  } else {
    rank <- if (n_ref >= 2600L) 128L else if (n_ref >= 1800L) 96L else if (n_ref >= 1200L) 80L else if (n_ref >= 900L) 64L else if (n_ref >= 600L) 48L else if (n_ref >= 400L) 32L else 0L
  }
  if (identical(precision, "double") && rank > 0L && n_ref >= 1600L) {
    rank <- rank + 16L
  }
  if (identical(precision, "mixed") && rank > 0L && n_ref < 800L) {
    rank <- 0L
  }
  rank_cap <- min(as.integer(n_template), as.integer(max(ns_source)))
  rank <- as.integer(min(rank, rank_cap))
  if (rank > 0L) {
    dense_work <- as.numeric(n_ref) * as.numeric(n_template) * as.numeric(n_ref + n_template)
    low_work <- 4.0 * as.numeric(rank) * (
      as.numeric(n_ref) * as.numeric(n_template) +
        as.numeric(n_ref) * as.numeric(rank) +
        as.numeric(n_template) * as.numeric(rank)
    )
    if (!is.finite(dense_work) || !is.finite(low_work) || low_work >= 0.90 * dense_work) {
      rank <- 0L
    }
  }
  list(
    n_ref = n_ref,
    structure_rank = max(0L, rank),
    rank_cap = max(1L, rank_cap)
  )
}

.prepare_multialign_subject_cache <- function(sets, structure_normalize) {
  C1_list <- lapply(sets, function(s) .maybe_rescale01(s$C, enabled = structure_normalize))
  p_list <- lapply(sets, function(s) s$w)
  ns <- vapply(C1_list, nrow, integer(1))
  list(
    C1_list = C1_list,
    p_list = p_list,
    ns = ns
  )
}

.precompute_c1_lowrank_cache <- function(C1_list, structure_rank) {
  if (!is.finite(structure_rank) || structure_rank <= 0L) {
    return(NULL)
  }
  A_scaled_list <- vector("list", length(C1_list))
  Bt_list <- vector("list", length(C1_list))
  for (i in seq_along(C1_list)) {
    C1 <- C1_list[[i]]
    r <- min(as.integer(structure_rank), nrow(C1), ncol(C1))
    if (!is.finite(r) || r <= 0L) {
      next
    }
    sv <- tryCatch(
      svd(C1, nu = r, nv = r),
      error = function(e) NULL
    )
    if (is.null(sv) || length(sv$d) < r || is.null(sv$u) || is.null(sv$v)) {
      next
    }
    d <- sv$d[seq_len(r)]
    A_scaled_list[[i]] <- sweep(sv$u[, seq_len(r), drop = FALSE], 2L, d, `*`)
    Bt_list[[i]] <- t(sv$v[, seq_len(r), drop = FALSE])
  }
  list(
    A_scaled_list = A_scaled_list,
    Bt_list = Bt_list
  )
}

.multialign_align_fixed_once <- function(
    sets,
    ids,
    template,
    method,
    alpha,
    epsilon,
    feature_cost_metric,
    structure_normalize,
    feature_cost_normalize,
    solver,
    precision,
    sinkhorn_method,
    max_iter,
    tol,
    sinkhorn_max_iter,
    sinkhorn_tol,
    reg_marginals,
    max_iter_ot,
    tol_ot,
    rescale_plan,
    check_every,
    structure_rank,
    use_cpp_batch,
    use_cpp_feature_fused,
    n_threads,
    batch_min_subjects,
    coarse_init_mode = "none",
    coarse_init_settings = NULL,
    warm_plans = NULL,
    subject_cache = NULL,
    c1_lowrank_cache = NULL) {
  if (is.null(subject_cache)) {
    ns <- vapply(sets, function(s) nrow(s$C), integer(1))
  } else {
    ns <- as.integer(subject_cache$ns)
  }
  C2 <- .maybe_rescale01(template$C, enabled = structure_normalize)
  q <- .assert_prob(template$w, nrow(C2), "template$w")
  use_coarse_init <- identical(method, "fgw_entropic") &&
    identical(coarse_init_mode, "sampled_graph") &&
    exists("sampled_gw_from_graphs", mode = "function")
  W2_coarse <- if (use_coarse_init) {
    .structure_cost_to_similarity(C2)
  } else {
    NULL
  }

  use_cpp_batch_effective <- isTRUE(use_cpp_batch) && length(sets) >= as.integer(batch_min_subjects)
  pin_blas_threads <- use_cpp_batch_effective &&
    is.finite(n_threads) &&
    (n_threads > 1L) &&
    !identical(Sys.getenv("RFUGW_PIN_BLAS_THREADS", unset = "1"), "0")
  C1_list <- if (is.null(subject_cache)) vector("list", length(sets)) else subject_cache$C1_list
  need_M_list <- TRUE
  p_list <- if (is.null(subject_cache)) vector("list", length(sets)) else subject_cache$p_list
  has_template_features <- !is.null(template$F)
  has_all_subject_features <- has_template_features && all(vapply(sets, function(s) !is.null(s$F), logical(1)))
  use_cpp_feature_fused <- isTRUE(use_cpp_feature_fused) &&
    use_cpp_batch_effective &&
    identical(method, "fgw_entropic") &&
    has_all_subject_features
  use_cpp_feature_batch <- use_cpp_batch_effective &&
    !use_cpp_feature_fused &&
    has_all_subject_features
  if (use_cpp_feature_fused && identical(method, "fgw_entropic")) {
    need_M_list <- FALSE
  }
  M_list <- if (need_M_list) vector("list", length(sets)) else NULL
  if (use_cpp_feature_batch) {
    M_list <- cpp_feature_cost_batch(
      F1_list = lapply(sets, function(s) s$F),
      F2 = template$F,
      use_euclidean = identical(feature_cost_metric, "euclidean"),
      rescale01 = isTRUE(feature_cost_normalize),
      n_threads = as.integer(n_threads)
    )
  }
  init_plan_list <- vector("list", length(sets))
  for (i in seq_along(sets)) {
    subj <- sets[[i]]
    C1 <- if (is.null(subject_cache)) .maybe_rescale01(subj$C, enabled = structure_normalize) else C1_list[[i]]
    if (use_cpp_feature_fused) {
      M <- NULL
    } else if (use_cpp_feature_batch) {
      M <- M_list[[i]]
    } else if (is.null(subj$F) || is.null(template$F)) {
      M <- matrix(0, nrow = nrow(C1), ncol = nrow(C2))
    } else {
      M <- .cross_feature_cost(
        subj$F,
        template$F,
        metric = feature_cost_metric,
        rescale01 = feature_cost_normalize
      )
    }
    if (is.null(subject_cache)) {
      C1_list[[i]] <- C1
      p_list[[i]] <- subj$w
    }
    if (need_M_list) {
      M_list[[i]] <- M
    }
    warm_i <- NULL
    if (!is.null(warm_plans)) {
      if (!is.null(names(warm_plans)) && ids[[i]] %in% names(warm_plans)) {
        warm_i <- warm_plans[[ids[[i]]]]
      } else if (length(warm_plans) >= i) {
        warm_i <- warm_plans[[i]]
      }
      if (!is.null(warm_i)) {
        if (!is.matrix(warm_i) || !is.numeric(warm_i) || nrow(warm_i) != nrow(C1) || ncol(warm_i) != nrow(C2)) {
          warm_i <- NULL
        }
      }
    }
    if (is.null(warm_i) && use_coarse_init) {
      W1_coarse <- .structure_cost_to_similarity(C1)
      coarse_args <- c(
        list(
          W1 = W1_coarse,
          W2 = W2_coarse,
          p = p_list[[i]],
          q = q,
          log = FALSE,
          return_embeddings = FALSE
        ),
        coarse_init_settings
      )
      warm_try <- tryCatch(
        do.call(sampled_gw_from_graphs, coarse_args),
        error = function(e) NULL
      )
      if (is.matrix(warm_try) &&
          is.numeric(warm_try) &&
          nrow(warm_try) == nrow(C1) &&
          ncol(warm_try) == nrow(C2) &&
          all(is.finite(warm_try)) &&
          all(warm_try >= 0)) {
        warm_i <- warm_try
      }
    }
    init_plan_list[i] <- list(warm_i)
  }

  objectives <- numeric(length(sets))
  iterations <- integer(length(sets))
  errors <- numeric(length(sets))
  couplings <- vector("list", length(sets))
  names(couplings) <- ids
  used_threads <- 1L
  used_cpp_batch <- FALSE
  batch_summary <- NULL
  batch_kernel_ms <- rep(NA_real_, length(sets))
  batch_feature_ms <- rep(NA_real_, length(sets))
  batch_solve_ms <- rep(NA_real_, length(sets))
  batch_lowrank_used <- rep(NA_integer_, length(sets))
  batch_c1_cache_used <- rep(NA_integer_, length(sets))
  batch_c2_cache_used <- rep(NA_integer_, length(sets))
  batch_square_cache_used <- rep(NA_integer_, length(sets))

  if (identical(method, "fgw_entropic") && use_cpp_batch_effective) {
    use_ppa <- identical(solver, "PPA")
    use_log_sinkhorn <- identical(sinkhorn_method, "log")
    use_mixed_precision <- identical(precision, "mixed")
    if (use_cpp_feature_fused) {
      batch <- .with_single_blas_threads(
        cpp_fgw_entropic_square_batch_features(
          F1_list = lapply(sets, function(s) s$F),
          C1_list = C1_list,
          p_list = p_list,
          F2 = template$F,
          C2 = C2,
          q = q,
          use_euclidean = identical(feature_cost_metric, "euclidean"),
          rescale01 = isTRUE(feature_cost_normalize),
          alpha = alpha,
          epsilon = epsilon,
          max_iter = as.integer(max_iter),
          tol = tol,
          sinkhorn_max_iter = as.integer(sinkhorn_max_iter),
          sinkhorn_tol = sinkhorn_tol,
          symmetric = TRUE,
          use_ppa = use_ppa,
          use_log_sinkhorn = use_log_sinkhorn,
          use_mixed_precision = use_mixed_precision,
          check_every = as.integer(check_every),
          init_plan_list = init_plan_list,
          c1_A_scaled_list = if (is.null(c1_lowrank_cache)) list() else c1_lowrank_cache$A_scaled_list,
          c1_Bt_list = if (is.null(c1_lowrank_cache)) list() else c1_lowrank_cache$Bt_list,
          approx_rank = as.integer(structure_rank),
          n_threads = as.integer(n_threads)
        ),
        enable = pin_blas_threads
      )
    } else {
      batch <- .with_single_blas_threads(
        cpp_fgw_entropic_square_batch(
          M_list = M_list,
          C1_list = C1_list,
          p_list = p_list,
          C2 = C2,
          q = q,
          alpha = alpha,
          epsilon = epsilon,
          max_iter = as.integer(max_iter),
          tol = tol,
          sinkhorn_max_iter = as.integer(sinkhorn_max_iter),
          sinkhorn_tol = sinkhorn_tol,
          symmetric = TRUE,
          use_ppa = use_ppa,
          use_log_sinkhorn = use_log_sinkhorn,
          use_mixed_precision = use_mixed_precision,
          check_every = as.integer(check_every),
          init_plan_list = init_plan_list,
          c1_A_scaled_list = if (is.null(c1_lowrank_cache)) list() else c1_lowrank_cache$A_scaled_list,
          c1_Bt_list = if (is.null(c1_lowrank_cache)) list() else c1_lowrank_cache$Bt_list,
          approx_rank = as.integer(structure_rank),
          n_threads = as.integer(n_threads)
        ),
        enable = pin_blas_threads
      )
    }
    for (i in seq_along(sets)) {
      couplings[[i]] <- batch$plans[[i]]
    }
    objectives[] <- as.numeric(batch$fgw_dist)
    iterations[] <- as.integer(batch$iterations)
    errors[] <- as.numeric(batch$error)
    used_threads <- as.integer(batch$used_threads)
    if (!is.null(batch$kernel_ms)) {
      batch_kernel_ms[] <- as.numeric(batch$kernel_ms)
    }
    if (!is.null(batch$feature_ms)) {
      batch_feature_ms[] <- as.numeric(batch$feature_ms)
    }
    if (!is.null(batch$solve_ms)) {
      batch_solve_ms[] <- as.numeric(batch$solve_ms)
    }
    if (!is.null(batch$lowrank_used)) {
      batch_lowrank_used[] <- as.integer(batch$lowrank_used)
    }
    if (!is.null(batch$c1_cache_used)) {
      batch_c1_cache_used[] <- as.integer(batch$c1_cache_used)
    }
    if (!is.null(batch$c2_cache_used)) {
      batch_c2_cache_used[] <- as.integer(batch$c2_cache_used)
    }
    if (!is.null(batch$square_cache_used)) {
      batch_square_cache_used[] <- as.integer(batch$square_cache_used)
    }
    batch_summary <- list(
      n_lowrank = if (!is.null(batch$n_lowrank)) as.integer(batch$n_lowrank) else sum(batch_lowrank_used, na.rm = TRUE),
      n_c1_cache = if (!is.null(batch$n_c1_cache)) as.integer(batch$n_c1_cache) else sum(batch_c1_cache_used, na.rm = TRUE),
      n_c2_cache = if (!is.null(batch$n_c2_cache)) as.integer(batch$n_c2_cache) else sum(batch_c2_cache_used, na.rm = TRUE),
      n_square_cache = if (!is.null(batch$n_square_cache)) as.integer(batch$n_square_cache) else sum(batch_square_cache_used, na.rm = TRUE)
    )
    used_cpp_batch <- TRUE
  } else if (identical(method, "fugw_kl") && use_cpp_batch_effective) {
    batch <- .with_single_blas_threads(
      cpp_fugw_kl_square_batch(
        Cx_list = C1_list,
        wx_list = p_list,
        M_list = M_list,
        Cy = C2,
        wy = q,
        reg_marginals = as.numeric(reg_marginals),
        epsilon = epsilon,
        alpha = alpha,
        max_iter = as.integer(max_iter),
        tol = tol,
        max_iter_ot = as.integer(max_iter_ot),
        tol_ot = tol_ot,
        rescale_plan = isTRUE(rescale_plan),
        check_every = as.integer(check_every),
        init_pi_list = init_plan_list,
        n_threads = as.integer(n_threads)
      ),
      enable = pin_blas_threads
    )
    for (i in seq_along(sets)) {
      couplings[[i]] <- batch$pi_samp[[i]]
    }
    objectives[] <- as.numeric(batch$fugw_cost)
    iterations[] <- as.integer(batch$iterations)
    errors[] <- as.numeric(batch$error)
    used_threads <- as.integer(batch$used_threads)
    if (!is.null(batch$kernel_ms)) {
      batch_kernel_ms[] <- as.numeric(batch$kernel_ms)
    }
    batch_summary <- list(
      n_lowrank = 0L,
      n_c1_cache = 0L,
      n_c2_cache = 0L,
      n_square_cache = 0L
    )
    used_cpp_batch <- TRUE
  } else {
    for (i in seq_along(sets)) {
      if (method == "fgw_entropic") {
        out <- fgw_entropic(
          M = M_list[[i]],
          C1 = C1_list[[i]],
          C2 = C2,
          p = p_list[[i]],
          q = q,
          alpha = alpha,
          epsilon = epsilon,
          max_iter = max_iter,
          tol = tol,
          sinkhorn_max_iter = sinkhorn_max_iter,
          sinkhorn_tol = sinkhorn_tol,
          init_plan = init_plan_list[[i]],
          sinkhorn_method = sinkhorn_method,
          precision = precision,
          solver = solver,
          symmetric = TRUE,
          check_every = check_every,
          structure_rank = structure_rank
        )
        couplings[[i]] <- out$plan
        objectives[[i]] <- out$fgw_dist
        iterations[[i]] <- out$iterations
        errors[[i]] <- out$error
      } else {
        out <- fugw_kl(
          Cx = C1_list[[i]],
          Cy = C2,
          wx = p_list[[i]],
          wy = q,
          reg_marginals = reg_marginals,
          epsilon = epsilon,
          alpha = alpha,
          M = M_list[[i]],
          init_pi = init_plan_list[[i]],
          max_iter = max_iter,
          tol = tol,
          max_iter_ot = max_iter_ot,
          tol_ot = tol_ot,
          rescale_plan = rescale_plan,
          check_every = check_every
        )
        couplings[[i]] <- out$pi_samp
        objectives[[i]] <- out$fugw_cost
        iterations[[i]] <- out$iterations
        errors[[i]] <- out$error
      }
    }
  }

  diagnostics <- data.frame(
    id = ids,
    n_nodes = ns,
    objective = objectives,
    iterations = iterations,
    error = errors,
    plan_mass = vapply(couplings, sum, numeric(1)),
    kernel_ms = batch_kernel_ms,
    feature_ms = batch_feature_ms,
    solve_ms = batch_solve_ms,
    lowrank_used = batch_lowrank_used,
    c1_cache_used = batch_c1_cache_used,
    c2_cache_used = batch_c2_cache_used,
    square_cache_used = batch_square_cache_used,
    stringsAsFactors = FALSE
  )

  list(
    couplings = couplings,
    objectives = objectives,
    objective_total = sum(objectives),
    diagnostics = diagnostics,
    batch_summary = batch_summary,
    used_threads = used_threads,
    used_cpp_batch = used_cpp_batch,
    used_cpp_feature_fused = use_cpp_feature_fused && used_cpp_batch
  )
}

.multialign_update_template <- function(
    sets,
    ids,
    couplings,
    template,
    structure_normalize,
    barycenter_eps = 1e-12) {
  K <- nrow(template$C)
  S <- length(sets)
  numC <- matrix(0, nrow = K, ncol = K)
  denC <- matrix(0, nrow = K, ncol = K)
  masses <- matrix(0, nrow = K, ncol = S)

  update_features <- !is.null(template$F) && all(vapply(sets, function(s) !is.null(s$F), logical(1)))
  if (update_features) {
    d <- ncol(template$F)
    numF <- matrix(0, nrow = K, ncol = d)
    denF <- rep(0, K)
  } else {
    numF <- NULL
    denF <- NULL
  }

  for (i in seq_along(sets)) {
    P <- couplings[[ids[[i]]]]
    if (is.null(P)) {
      P <- couplings[[i]]
    }
    C1 <- .maybe_rescale01(sets[[i]]$C, enabled = structure_normalize)
    ms <- colSums(P)
    numC <- numC + t(P) %*% C1 %*% P
    denC <- denC + tcrossprod(ms, ms)
    masses[, i] <- ms
    if (update_features) {
      numF <- numF + t(P) %*% sets[[i]]$F
      denF <- denF + ms
    }
  }

  C_new <- numC / pmax(denC, barycenter_eps)
  C_new <- 0.5 * (C_new + t(C_new))
  diag(C_new) <- 0
  C_new <- .maybe_rescale01(C_new, enabled = structure_normalize)

  w_new <- rowMeans(masses)
  if (!all(is.finite(w_new)) || sum(w_new) <= 0) {
    w_new <- template$w
  }
  w_new <- w_new / sum(w_new)

  F_new <- template$F
  if (update_features) {
    F_new <- numF / pmax(denF, barycenter_eps)
  }

  list(C = C_new, F = F_new, w = w_new, id = template$id)
}

.multialign_template_delta <- function(template_old, template_new, eps = 1e-12) {
  dC <- sqrt(sum((template_old$C - template_new$C)^2)) / max(eps, sqrt(sum(template_old$C^2)))
  dW <- sqrt(sum((template_old$w - template_new$w)^2)) / max(eps, sqrt(sum(template_old$w^2)))
  dF <- 0
  if (!is.null(template_old$F) && !is.null(template_new$F)) {
    dF <- sqrt(sum((template_old$F - template_new$F)^2)) / max(eps, sqrt(sum(template_old$F^2)))
  }
  max(dC, dW, dF)
}

#' Build a Fixed Template Support for Multiset Alignment
#'
#' Builds a fixed template support from pooled subject node features using
#' k-means, then defines template geometry as pairwise distances between
#' template feature centroids.
#'
#' @param subjects A list of subject sets, each with at least `C` and optional
#'   `F`, `w`, `id`.
#' @param k Number of template nodes. Defaults to rounded median subject size.
#' @param feature_normalization Feature normalization mode.
#' @param feature_cost_metric Feature distance metric.
#' @param seed Optional random seed for k-means reproducibility.
#' @return A template list with `C`, `F`, `w`, `id`.
#' @export
multialign_make_template <- function(
    subjects,
    k = NULL,
    feature_normalization = c("none", "zscore"),
    feature_cost_metric = c("euclidean", "sqeuclidean"),
    seed = 1L) {
  if (!is.list(subjects) || length(subjects) == 0L) {
    stop("`subjects` must be a non-empty list.", call. = FALSE)
  }
  feature_normalization <- match.arg(feature_normalization)
  feature_cost_metric <- match.arg(feature_cost_metric)

  sets <- lapply(seq_along(subjects), function(i) {
    .validate_set(subjects[[i]], sprintf("subjects[[%d]]", i))
  })
  features <- lapply(sets, function(s) s$F)
  if (!all(vapply(features, function(F) !is.null(F), logical(1)))) {
    stop("All subjects must include `F` when `template` is not provided.", call. = FALSE)
  }
  dvals <- unique(vapply(features, ncol, integer(1)))
  if (length(dvals) != 1L) {
    stop("All subject feature matrices must have the same number of columns.", call. = FALSE)
  }

  features <- .normalize_feature_bundle(features, mode = feature_normalization)
  pooled <- do.call(rbind, features)
  if (is.null(k)) {
    k <- as.integer(round(stats::median(vapply(sets, function(s) nrow(s$C), integer(1)))))
  }
  k <- as.integer(k)
  if (!is.finite(k) || k < 2L) {
    stop("`k` must be an integer >= 2.", call. = FALSE)
  }
  if (nrow(pooled) < k) {
    stop("`k` cannot exceed the total number of pooled nodes.", call. = FALSE)
  }

  if (!is.null(seed)) set.seed(as.integer(seed))
  km <- stats::kmeans(pooled, centers = k, nstart = 5, iter.max = 100)
  F <- km$centers
  C <- .cross_feature_cost(F, F, metric = feature_cost_metric, rescale01 = TRUE)
  list(
    C = C,
    F = F,
    w = rep(1 / nrow(F), nrow(F)),
    id = "template"
  )
}

#' Fit Multiset Alignment to a Fixed Template
#'
#' Aligns a list of attributed metric-measure sets to a fixed template using
#' FGW (`method = "fgw_entropic"`) or unbalanced FUGW
#' (`method = "fugw_kl"`). This is the first stage of general multiset
#' alignment and returns subject-to-template couplings.
#'
#' @param subjects A non-empty list of subject sets. Each subject is a list with
#'   `C` (required), optional `F`, optional `w`, optional `id`.
#' @param template_mode Template strategy: `"fixed"` (fit couplings to a fixed
#'   template) or `"learned"` (alternating align-and-update barycenter-style
#'   template refinement).
#' @param template Fixed template set with the same fields as a subject list.
#'   If `NULL`, a template is built from pooled subject features via
#'   [multialign_make_template()].
#' @param k_template Number of template nodes used only when `template = NULL`.
#' @param method Alignment backend: `"fgw_entropic"` or `"fugw_kl"`.
#' @param alpha Feature/structure trade-off parameter.
#' @param epsilon Entropic regularization for the selected backend.
#' @param feature_normalization Feature normalization mode before cost assembly.
#' @param feature_cost_metric Feature cost metric used to form cross-set `M`.
#' @param structure_normalize If `TRUE`, each structure matrix is scaled to have
#'   max value 1 before solving.
#' @param structure_knn Optional kNN sparsification level for structure matrices.
#'   If positive, rows keep `k+1` smallest entries (including diagonal) and
#'   other entries are set to a large fill value before solving.
#' @param structure_knn_min_n Minimum matrix size for applying `structure_knn`.
#' @param feature_cost_normalize If `TRUE`, each cross-feature cost matrix is
#'   scaled to have max value 1.
#' @param solver `fgw_entropic` solver (`"PGD"` or `"PPA"`).
#' @param precision Numeric precision for `fgw_entropic` (defaults to `"mixed"`
#'   for speed in multiset alignment).
#' @param sinkhorn_method Sinkhorn variant for `fgw_entropic`.
#' @param max_iter Outer iterations (`fgw_entropic`) or BCD iterations (`fugw_kl`).
#' @param tol Outer stopping tolerance for selected backend.
#' @param sinkhorn_max_iter Inner Sinkhorn iterations for `fgw_entropic`.
#' @param sinkhorn_tol Inner Sinkhorn tolerance for `fgw_entropic`.
#' @param reg_marginals Marginal relaxation for `fugw_kl`.
#' @param max_iter_ot Inner Sinkhorn iterations for `fugw_kl`.
#' @param tol_ot Inner Sinkhorn tolerance for `fugw_kl`.
#' @param rescale_plan Plan rescaling option for `fugw_kl`.
#' @param check_every Iteration check interval for selected backend.
#' @param coarse_init Optional coarse coupling initialization mode. `"none"`
#'   disables coarse initialization. `"sampled_graph"` computes a fast
#'   sampled-GW plan on graph-derived diffusion coordinates to initialize FGW.
#' @param coarse_init_components Diffusion components for `"sampled_graph"`.
#' @param coarse_init_diffusion_time Diffusion-time exponent for coarse init.
#' @param coarse_init_self_loop Diagonal stabilizer for coarse graph init.
#' @param coarse_init_nb_samples Sample counts for coarse sampled-GW gradients.
#' @param coarse_init_epsilon Entropic regularization for coarse init.
#' @param coarse_init_max_iter Iterations for coarse sampled-GW init.
#' @param coarse_init_use_cpp Use C++ fast path for coarse sampled-GW init.
#' @param coarse_init_sampling Sampling mode for coarse init.
#' @param structure_rank Optional low-rank rank used in the FGW tensor-product
#'   update. Use `0`/`"off"` to disable or `"auto"` for size-aware rank
#'   selection.
#' @param autotune If `TRUE`, apply size-aware caps to FGW iterations for faster
#'   large-`n` runs.
#' @param autotune_level Autotune aggressiveness.
#' @param seed Optional seed used only when `template = NULL`.
#' @param use_cpp_batch If `TRUE`, use the C++ batched FGW path for fixed-template
#'   `method = "fgw_entropic"`.
#' @param use_cpp_feature_fused If `TRUE` and conditions allow, use a fused C++
#'   path that computes feature costs and solves FGW in one batched call without
#'   materializing `M_list` in R.
#' @param use_warm_start If `TRUE`, reuse previous couplings as initialization
#'   across learned-template outer iterations (and for `fugw_kl` fixed solves).
#' @param n_threads Requested subject-parallel threads for the C++ batched path.
#' @param batch_min_subjects Minimum number of subjects before enabling batch path.
#' @param template_max_iter Maximum alternating updates for `template_mode = "learned"`.
#' @param template_tol Relative stopping tolerance for learned-template updates.
#' @param template_relax Relaxation in `(0, 1]` when updating learned templates.
#' @param barycenter_eps Stabilizer used in learned-template denominator updates.
#' @param final_refit If `TRUE`, run one final alignment against the last updated
#'   learned template. If `FALSE` (default), skip this extra pass for speed.
#' @return A list with `couplings`, `objectives`, `template`, `diagnostics`,
#'   and run settings.
#' @export
multialign_fit <- function(
    subjects,
    template_mode = c("fixed", "learned"),
    template = NULL,
    k_template = NULL,
    method = c("fgw_entropic", "fugw_kl"),
    alpha = 0.5,
    epsilon = 0.05,
    feature_normalization = c("none", "zscore"),
    feature_cost_metric = c("euclidean", "sqeuclidean"),
    structure_normalize = TRUE,
    structure_knn = 0L,
    structure_knn_min_n = 600L,
    feature_cost_normalize = TRUE,
    solver = c("PGD", "PPA"),
    precision = c("mixed", "double"),
    sinkhorn_method = c("scaling", "log"),
    max_iter = 200L,
    tol = 1e-9,
    sinkhorn_max_iter = 500L,
    sinkhorn_tol = 1e-9,
    reg_marginals = c(10, 10),
    max_iter_ot = 500L,
    tol_ot = 1e-7,
    rescale_plan = TRUE,
    check_every = 10L,
    coarse_init = c("none", "sampled_graph"),
    coarse_init_components = 10L,
    coarse_init_diffusion_time = 1,
    coarse_init_self_loop = 1e-6,
    coarse_init_nb_samples = c(8L, 2L),
    coarse_init_epsilon = 0.1,
    coarse_init_max_iter = 40L,
    coarse_init_use_cpp = TRUE,
    coarse_init_sampling = c("stochastic", "deterministic"),
    structure_rank = 0L,
    autotune = FALSE,
    autotune_level = c("moderate", "aggressive"),
    seed = 1L,
    use_cpp_batch = TRUE,
    use_cpp_feature_fused = TRUE,
    use_warm_start = TRUE,
    n_threads = 1L,
    batch_min_subjects = 2L,
    template_max_iter = 8L,
    template_tol = 1e-6,
    template_relax = 1.0,
    barycenter_eps = 1e-12,
    final_refit = FALSE) {
  if (!is.list(subjects) || length(subjects) == 0L) {
    stop("`subjects` must be a non-empty list.", call. = FALSE)
  }
  template_mode <- match.arg(template_mode)
  if (!is.finite(template_relax) || template_relax <= 0 || template_relax > 1) {
    stop("`template_relax` must be in (0, 1].", call. = FALSE)
  }
  if (!is.finite(template_tol) || template_tol <= 0) {
    stop("`template_tol` must be positive.", call. = FALSE)
  }
  template_max_iter <- as.integer(template_max_iter)
  if (template_max_iter < 1L) {
    stop("`template_max_iter` must be >= 1.", call. = FALSE)
  }
  method <- match.arg(method)
  solver <- match.arg(solver)
  precision <- match.arg(precision)
  sinkhorn_method <- match.arg(sinkhorn_method)
  coarse_init <- match.arg(coarse_init)
  coarse_init_sampling <- match.arg(coarse_init_sampling)
  autotune_level <- match.arg(autotune_level)
  feature_normalization <- match.arg(feature_normalization)
  feature_cost_metric <- match.arg(feature_cost_metric)
  structure_knn <- as.integer(structure_knn)[1]
  structure_knn_min_n <- as.integer(structure_knn_min_n)[1]
  if (!is.finite(structure_knn) || structure_knn < 0L) {
    stop("`structure_knn` must be a nonnegative integer.", call. = FALSE)
  }
  if (!is.finite(structure_knn_min_n) || structure_knn_min_n < 2L) {
    stop("`structure_knn_min_n` must be an integer >= 2.", call. = FALSE)
  }
  structure_rank_spec <- .resolve_structure_rank_request(structure_rank)
  structure_rank <- as.integer(structure_rank_spec$rank)
  structure_rank_auto <- isTRUE(structure_rank_spec$auto)
  structure_rank_requested <- structure_rank_spec$requested

  sets <- lapply(seq_along(subjects), function(i) {
    .validate_set(subjects[[i]], sprintf("subjects[[%d]]", i))
  })
  ids <- .subject_ids_unique(vapply(sets, function(s) s$id, character(1)))
  for (i in seq_along(sets)) sets[[i]]$id <- ids[[i]]

  if (is.null(template)) {
    template <- multialign_make_template(
      subjects = sets,
      k = k_template,
      feature_normalization = feature_normalization,
      feature_cost_metric = feature_cost_metric,
      seed = seed
    )
  }
  template <- .validate_set(template, "template")

  if (structure_knn > 0L) {
    sets <- lapply(sets, function(s) {
      if (nrow(s$C) >= structure_knn_min_n) {
        s$C <- .knn_sparsify_structure(s$C, k = structure_knn)
      }
      s
    })
    if (nrow(template$C) >= structure_knn_min_n) {
      template$C <- .knn_sparsify_structure(template$C, k = structure_knn)
    }
  }

  dims <- vapply(
    c(lapply(sets, function(s) s$F), list(template$F)),
    function(F) if (is.null(F)) NA_integer_ else ncol(F),
    integer(1)
  )
  dims <- unique(dims[!is.na(dims)])
  if (length(dims) > 1L) {
    stop("All non-null feature matrices must have the same number of columns.", call. = FALSE)
  }

  if (feature_normalization == "zscore" && length(dims) == 1L) {
    bundle <- c(lapply(sets, function(s) s$F), list(template$F))
    bundle <- .normalize_feature_bundle(bundle, mode = "zscore")
    for (i in seq_along(sets)) sets[[i]]$F <- bundle[[i]]
    template$F <- bundle[[length(bundle)]]
  }

  template$C <- .maybe_rescale01(template$C, enabled = structure_normalize)
  template$w <- .assert_prob(template$w, nrow(template$C), "template$w")
  subject_cache <- .prepare_multialign_subject_cache(
    sets = sets,
    structure_normalize = structure_normalize
  )
  ns_source <- subject_cache$ns

  max_iter_eff <- as.integer(max_iter)
  sinkhorn_max_iter_eff <- as.integer(sinkhorn_max_iter)
  check_every_eff <- as.integer(check_every)
  coarse_init_settings <- NULL
  if (identical(coarse_init, "sampled_graph")) {
    coarse_init_settings <- list(
      n_components = as.integer(coarse_init_components),
      diffusion_time = as.numeric(coarse_init_diffusion_time),
      self_loop = as.numeric(coarse_init_self_loop),
      nb_samples_grad = as.integer(coarse_init_nb_samples),
      epsilon = as.numeric(coarse_init_epsilon),
      max_iter = as.integer(coarse_init_max_iter),
      use_cpp = isTRUE(coarse_init_use_cpp),
      sampling = coarse_init_sampling
    )
  }
  structure_rank_eff <- structure_rank
  tuned_controls <- NULL
  if (isTRUE(autotune) && identical(method, "fgw_entropic")) {
    tuned <- .autotune_multialign_controls(
      ns_source = ns_source,
      n_template = nrow(template$C),
      precision = precision,
      max_iter = max_iter_eff,
      sinkhorn_max_iter = sinkhorn_max_iter_eff,
      check_every = check_every_eff,
      level = autotune_level
    )
    max_iter_eff <- tuned$max_iter
    sinkhorn_max_iter_eff <- tuned$sinkhorn_max_iter
    check_every_eff <- tuned$check_every
    tuned_controls <- tuned
  }
  structure_rank_adaptive <- identical(method, "fgw_entropic") &&
    (isTRUE(structure_rank_auto) || (isTRUE(autotune) && structure_rank == 0L))
  if (structure_rank_adaptive) {
    tuned_rank <- .autotune_structure_rank(
      ns_source = ns_source,
      n_template = nrow(template$C),
      precision = precision,
      level = autotune_level
    )
    structure_rank_eff <- tuned_rank$structure_rank
    expected_template_passes <- if (identical(template_mode, "learned")) {
      as.integer(template_max_iter) + if (isTRUE(final_refit)) 1L else 0L
    } else {
      1L
    }
    expected_template_passes <- max(1L, as.integer(expected_template_passes))
    rank_gate_reason <- NULL
    if (identical(precision, "mixed") &&
        structure_rank_eff > 0L &&
        expected_template_passes <= 2L &&
        tuned_rank$n_ref < 1600L) {
      structure_rank_eff <- 0L
      rank_gate_reason <- sprintf(
        "disabled mixed low-rank (n_ref=%d, expected_passes=%d)",
        tuned_rank$n_ref,
        expected_template_passes
      )
    }
    if (is.null(tuned_controls)) {
      tuned_controls <- tuned_rank
    } else {
      tuned_controls$structure_rank <- tuned_rank$structure_rank
      tuned_controls$structure_rank_n_ref <- tuned_rank$n_ref
      tuned_controls$structure_rank_cap <- tuned_rank$rank_cap
    }
    tuned_controls$structure_rank_effective <- structure_rank_eff
    tuned_controls$structure_rank_expected_passes <- expected_template_passes
    tuned_controls$structure_rank_gated <- isTRUE(structure_rank_eff < tuned_rank$structure_rank)
    tuned_controls$structure_rank_gate_reason <- rank_gate_reason
  }
  if (!isTRUE(structure_rank_adaptive) && !is.null(tuned_controls)) {
    tuned_controls$structure_rank <- structure_rank_eff
  }
  c1_lowrank_cache <- NULL
  if (identical(method, "fgw_entropic") &&
      isTRUE(use_cpp_batch) &&
      structure_rank_eff > 0L) {
    c1_lowrank_cache <- .precompute_c1_lowrank_cache(
      C1_list = subject_cache$C1_list,
      structure_rank = structure_rank_eff
    )
  }

  fit_once <- function(template_arg, warm_plans = NULL) {
    template_solve <- template_arg
    if (structure_knn > 0L && nrow(template_solve$C) >= structure_knn_min_n) {
      template_solve$C <- .knn_sparsify_structure(template_solve$C, k = structure_knn)
    }
    .multialign_align_fixed_once(
      sets = sets,
      ids = ids,
      template = template_solve,
      method = method,
      alpha = alpha,
      epsilon = epsilon,
      feature_cost_metric = feature_cost_metric,
      structure_normalize = structure_normalize,
      feature_cost_normalize = feature_cost_normalize,
      solver = solver,
      precision = precision,
      sinkhorn_method = sinkhorn_method,
      max_iter = max_iter_eff,
      tol = tol,
      sinkhorn_max_iter = sinkhorn_max_iter_eff,
      sinkhorn_tol = sinkhorn_tol,
      reg_marginals = reg_marginals,
      max_iter_ot = max_iter_ot,
      tol_ot = tol_ot,
      rescale_plan = rescale_plan,
      check_every = check_every_eff,
      structure_rank = structure_rank_eff,
      use_cpp_batch = use_cpp_batch,
      use_cpp_feature_fused = use_cpp_feature_fused,
      n_threads = n_threads,
      batch_min_subjects = batch_min_subjects,
      coarse_init_mode = coarse_init,
      coarse_init_settings = coarse_init_settings,
      warm_plans = warm_plans,
      subject_cache = subject_cache,
      c1_lowrank_cache = c1_lowrank_cache
    )
  }

  if (identical(template_mode, "fixed")) {
    fit <- fit_once(template)
    return(structure(list(
      couplings = fit$couplings,
      objectives = fit$objectives,
      objective_total = fit$objective_total,
      template = template,
      diagnostics = fit$diagnostics,
      template_mode = template_mode,
      method = method,
      alpha = alpha,
      epsilon = epsilon,
      feature_normalization = feature_normalization,
      feature_cost_metric = feature_cost_metric,
      used_threads = fit$used_threads,
      used_cpp_batch = fit$used_cpp_batch,
      used_cpp_feature_fused = fit$used_cpp_feature_fused,
      batch_summary = fit$batch_summary,
      structure_rank = structure_rank_eff,
      structure_rank_requested = structure_rank_requested,
      structure_rank_adaptive = structure_rank_adaptive,
      structure_knn = structure_knn,
      structure_knn_min_n = structure_knn_min_n,
      autotune = isTRUE(autotune),
      autotuned_controls = tuned_controls,
      coarse_init = coarse_init,
      template_history = NULL,
      outer_iterations = 1L
    ), class = "rfugw_multialign"))
  }

  history <- vector("list", template_max_iter)
  tpl <- template
  prev_obj <- NA_real_
  outer_done <- 0L
  warm_plans <- NULL
  fit_last <- NULL
  for (outer in seq_len(template_max_iter)) {
    fit <- fit_once(tpl, warm_plans = if (isTRUE(use_warm_start)) warm_plans else NULL)
    fit_last <- fit
    tpl_new <- .multialign_update_template(
      sets = sets,
      ids = ids,
      couplings = fit$couplings,
      template = tpl,
      structure_normalize = structure_normalize,
      barycenter_eps = barycenter_eps
    )
    if (template_relax < 1) {
      tpl_new$C <- (1 - template_relax) * tpl$C + template_relax * tpl_new$C
      tpl_new$w <- (1 - template_relax) * tpl$w + template_relax * tpl_new$w
      tpl_new$w <- tpl_new$w / sum(tpl_new$w)
      if (!is.null(tpl$F) && !is.null(tpl_new$F)) {
        tpl_new$F <- (1 - template_relax) * tpl$F + template_relax * tpl_new$F
      }
    }
    delta_template <- .multialign_template_delta(tpl, tpl_new, eps = barycenter_eps)
    delta_objective <- if (is.finite(prev_obj)) {
      abs(fit$objective_total - prev_obj) / max(1, abs(prev_obj))
    } else {
      Inf
    }

    history[[outer]] <- data.frame(
      outer_iter = outer,
      objective_total = fit$objective_total,
      delta_objective = delta_objective,
      delta_template = delta_template,
      stringsAsFactors = FALSE
    )
    outer_done <- outer
    if (isTRUE(use_warm_start)) {
      warm_plans <- fit$couplings
    }
    converged <- (outer > 1L && delta_template < template_tol && delta_objective < template_tol)
    reached_max <- (outer == template_max_iter)
    if (converged || reached_max) {
      if (isTRUE(final_refit)) {
        tpl <- tpl_new
      }
      prev_obj <- fit$objective_total
      break
    }
    tpl <- tpl_new
    prev_obj <- fit$objective_total
  }
  history <- do.call(rbind, history[seq_len(outer_done)])

  fit_final <- if (isTRUE(final_refit)) {
    fit_once(
      tpl,
      warm_plans = if (isTRUE(use_warm_start)) warm_plans else NULL
    )
  } else {
    fit_last
  }
  structure(list(
    couplings = fit_final$couplings,
    objectives = fit_final$objectives,
    objective_total = fit_final$objective_total,
    template = tpl,
    diagnostics = fit_final$diagnostics,
    template_mode = template_mode,
    method = method,
    alpha = alpha,
    epsilon = epsilon,
    feature_normalization = feature_normalization,
    feature_cost_metric = feature_cost_metric,
    used_threads = fit_final$used_threads,
    used_cpp_batch = fit_final$used_cpp_batch,
    used_cpp_feature_fused = fit_final$used_cpp_feature_fused,
    batch_summary = fit_final$batch_summary,
    structure_rank = structure_rank_eff,
    structure_rank_requested = structure_rank_requested,
    structure_rank_adaptive = structure_rank_adaptive,
    structure_knn = structure_knn,
    structure_knn_min_n = structure_knn_min_n,
    autotune = isTRUE(autotune),
    autotuned_controls = tuned_controls,
    coarse_init = coarse_init,
    template_history = history,
    outer_iterations = outer_done
  ), class = "rfugw_multialign")
}

#' Normalize a Coupling Matrix
#'
#' @param P Coupling/transport matrix.
#' @param mode Normalization mode.
#' @param eps Stabilizer for division.
#' @return Normalized coupling matrix.
#' @export
multialign_normalize_plan <- function(P, mode = c("none", "row", "col", "sum"), eps = 1e-15) {
  .assert_matrix(P, "P")
  mode <- match.arg(mode)
  if (mode == "none") return(P)
  if (!is.finite(eps) || eps <= 0) {
    stop("`eps` must be positive.", call. = FALSE)
  }

  if (mode == "row") {
    rs <- rowSums(P)
    return(P / pmax(rs, eps))
  }
  if (mode == "col") {
    cs <- colSums(P)
    return(sweep(P, 2, pmax(cs, eps), "/"))
  }
  P / max(sum(P), eps)
}

#' Project a Subject Matrix into Template Space
#'
#' Computes `t(Pn) %*% R %*% Pn` where `Pn` is the normalized coupling.
#'
#' @param P Subject-to-template coupling (`n_subject x n_template`).
#' @param R Subject square matrix (`n_subject x n_subject`).
#' @param normalize Plan normalization mode passed to
#'   [multialign_normalize_plan()].
#' @param eps Stabilizer for normalization.
#' @return Projected template-space matrix.
#' @export
multialign_project_matrix <- function(P, R, normalize = c("none", "row", "col", "sum"), eps = 1e-15) {
  .assert_matrix(P, "P")
  .assert_matrix(R, "R")
  if (nrow(R) != ncol(R)) {
    stop("`R` must be square.", call. = FALSE)
  }
  if (nrow(P) != nrow(R)) {
    stop("`nrow(P)` must equal `nrow(R)`.", call. = FALSE)
  }
  Pn <- multialign_normalize_plan(P, mode = normalize, eps = eps)
  t(Pn) %*% R %*% Pn
}

#' Project Template Features to Subject Nodes
#'
#' Computes `Pn %*% E_template`, typically with row-normalized `Pn`.
#'
#' @param P Subject-to-template coupling (`n_subject x n_template`).
#' @param E_template Template feature/embedding matrix (`n_template x d`).
#' @param normalize Plan normalization mode.
#' @param eps Stabilizer for normalization.
#' @return Subject-space projected feature matrix (`n_subject x d`).
#' @export
multialign_project_features <- function(P, E_template, normalize = c("row", "none", "col", "sum"), eps = 1e-15) {
  .assert_matrix(P, "P")
  .assert_matrix(E_template, "E_template")
  normalize <- match.arg(normalize)
  if (ncol(P) != nrow(E_template)) {
    stop("`ncol(P)` must equal `nrow(E_template)`.", call. = FALSE)
  }
  Pn <- multialign_normalize_plan(P, mode = normalize, eps = eps)
  Pn %*% E_template
}
