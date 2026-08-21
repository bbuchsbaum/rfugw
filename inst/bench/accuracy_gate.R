suppressPackageStartupMessages({
  rlib <- Sys.getenv("RFUGW_RLIB", unset = "")
  if (nzchar(rlib)) {
    .libPaths(c(rlib, .libPaths()))
  }
  library(rfugw)
  library(jsonlite)
})

if (!"graph_diffusion_coordinates" %in% getNamespaceExports("rfugw") &&
    requireNamespace("devtools", quietly = TRUE)) {
  devtools::load_all(".", quiet = TRUE)
}

.resolve_fixture_dir <- function(fixture_dir) {
  candidates <- unique(c(
    fixture_dir,
    "inst/extdata/fixtures",
    "rfugw/inst/extdata/fixtures"
  ))
  hit <- candidates[file.exists(candidates)]
  if (length(hit) == 0L) {
    stop(
      "Could not locate fixture directory. Tried: ",
      paste(candidates, collapse = ", "),
      call. = FALSE
    )
  }
  hit[[1L]]
}

.fixture_path <- function(name, fixture_dir = "inst/extdata/fixtures") {
  file.path(.resolve_fixture_dir(fixture_dir), name)
}

.read_fixture <- function(name, fixture_dir = "inst/extdata/fixtures") {
  fromJSON(.fixture_path(name, fixture_dir), simplifyVector = TRUE)
}

.check_close <- function(x, y, tol, name) {
  d <- max(abs(x - y))
  ok <- is.finite(d) && d <= tol
  list(ok = ok, max_abs = d, name = name, tol = tol)
}

run_accuracy_gate <- function(
    fixture_dir = "inst/extdata/fixtures",
    stop_on_fail = TRUE,
    verbose = TRUE) {
  results <- list()
  idx <- 1L
  push <- function(r) {
    results[[idx]] <<- r
    idx <<- idx + 1L
  }

  # 1) FGW entropic vs POT fixture. The fixture's original regularization is
  # outside the certified scaling-domain envelope, so use the log backend for
  # the reference comparison.
  fx <- .read_fixture("fgw_entropic_square_fixture.json", fixture_dir = fixture_dir)
  out_reference <- rfugw::fgw_entropic(
    M = fx$inputs$M,
    C1 = fx$inputs$C1,
    C2 = fx$inputs$C2,
    p = fx$inputs$p,
    q = fx$inputs$q,
    alpha = fx$params$alpha,
    epsilon = fx$params$epsilon,
    max_iter = fx$params$max_iter,
    tol = fx$params$tol,
    sinkhorn_max_iter = fx$params$sinkhorn_numItermax,
    sinkhorn_tol = fx$params$sinkhorn_stopThr,
    sinkhorn_method = "log",
    precision = "double",
    symmetric = TRUE,
    solver = fx$params$solver
  )
  push(cbind(
    check = "fgw_log_vs_pot_obj",
    as.data.frame(.check_close(out_reference$fgw_dist, fx$outputs$fgw_dist, 1e-6, "fgw_dist"), stringsAsFactors = FALSE)
  ))
  push(cbind(
    check = "fgw_log_vs_pot_plan",
    as.data.frame(.check_close(out_reference$plan, fx$outputs$plan, 1e-5, "plan"), stringsAsFactors = FALSE)
  ))
  push(cbind(
    check = "fgw_log_marginals_row",
    as.data.frame(.check_close(rowSums(out_reference$plan), fx$inputs$p, 1e-7, "row_marg"), stringsAsFactors = FALSE)
  ))
  push(cbind(
    check = "fgw_log_marginals_col",
    as.data.frame(.check_close(colSums(out_reference$plan), fx$inputs$q, 1e-7, "col_marg"), stringsAsFactors = FALSE)
  ))

  # 2) FGW log Sinkhorn should remain close to scaling inside scaling's
  # certified envelope. Relax the stopping tolerances as well so section 3
  # exercises the actual mixed-precision path instead of automatic promotion.
  safe_epsilon <- 4 * fx$params$epsilon
  safe_tol <- 1e-6
  out_scaling <- rfugw::fgw_entropic(
    M = fx$inputs$M,
    C1 = fx$inputs$C1,
    C2 = fx$inputs$C2,
    p = fx$inputs$p,
    q = fx$inputs$q,
    alpha = fx$params$alpha,
    epsilon = safe_epsilon,
    max_iter = fx$params$max_iter,
    tol = safe_tol,
    sinkhorn_max_iter = fx$params$sinkhorn_numItermax,
    sinkhorn_tol = safe_tol,
    sinkhorn_method = "scaling",
    precision = "double",
    symmetric = TRUE,
    solver = fx$params$solver
  )
  push(cbind(
    check = "fgw_scaling_marginals_row",
    as.data.frame(.check_close(rowSums(out_scaling$plan), fx$inputs$p, 2e-6, "row_marg"), stringsAsFactors = FALSE)
  ))
  push(cbind(
    check = "fgw_scaling_marginals_col",
    as.data.frame(.check_close(colSums(out_scaling$plan), fx$inputs$q, 2e-6, "col_marg"), stringsAsFactors = FALSE)
  ))
  out_log <- rfugw::fgw_entropic(
    M = fx$inputs$M,
    C1 = fx$inputs$C1,
    C2 = fx$inputs$C2,
    p = fx$inputs$p,
    q = fx$inputs$q,
    alpha = fx$params$alpha,
    epsilon = safe_epsilon,
    max_iter = fx$params$max_iter,
    tol = safe_tol,
    sinkhorn_max_iter = fx$params$sinkhorn_numItermax,
    sinkhorn_tol = safe_tol,
    sinkhorn_method = "log",
    precision = "double",
    symmetric = TRUE,
    solver = fx$params$solver
  )
  push(cbind(
    check = "fgw_log_vs_scaling_obj",
    as.data.frame(.check_close(out_log$fgw_dist, out_scaling$fgw_dist, 1e-6, "fgw_dist"), stringsAsFactors = FALSE)
  ))
  push(cbind(
    check = "fgw_log_vs_scaling_plan",
    as.data.frame(.check_close(out_log$plan, out_scaling$plan, 1e-5, "plan"), stringsAsFactors = FALSE)
  ))

  # 3) FGW mixed precision should remain close to double precision.
  out_mixed <- rfugw::fgw_entropic(
    M = fx$inputs$M,
    C1 = fx$inputs$C1,
    C2 = fx$inputs$C2,
    p = fx$inputs$p,
    q = fx$inputs$q,
    alpha = fx$params$alpha,
    epsilon = safe_epsilon,
    max_iter = fx$params$max_iter,
    tol = safe_tol,
    sinkhorn_max_iter = fx$params$sinkhorn_numItermax,
    sinkhorn_tol = safe_tol,
    sinkhorn_method = "scaling",
    precision = "mixed",
    symmetric = TRUE,
    solver = fx$params$solver
  )
  push(cbind(
    check = "fgw_mixed_vs_double_obj",
    as.data.frame(.check_close(out_mixed$fgw_dist, out_scaling$fgw_dist, 5e-5, "fgw_dist"), stringsAsFactors = FALSE)
  ))
  push(cbind(
    check = "fgw_mixed_marginals_row",
    as.data.frame(.check_close(rowSums(out_mixed$plan), fx$inputs$p, 2e-6, "row_marg"), stringsAsFactors = FALSE)
  ))
  push(cbind(
    check = "fgw_mixed_marginals_col",
    as.data.frame(.check_close(colSums(out_mixed$plan), fx$inputs$q, 2e-6, "col_marg"), stringsAsFactors = FALSE)
  ))

  # 4) FGW low-rank full-rank fallback should match dense solve.
  out_lowrank_full <- rfugw::fgw_entropic(
    M = fx$inputs$M,
    C1 = fx$inputs$C1,
    C2 = fx$inputs$C2,
    p = fx$inputs$p,
    q = fx$inputs$q,
    alpha = fx$params$alpha,
    epsilon = fx$params$epsilon,
    max_iter = fx$params$max_iter,
    tol = fx$params$tol,
    sinkhorn_max_iter = fx$params$sinkhorn_numItermax,
    sinkhorn_tol = fx$params$sinkhorn_stopThr,
    sinkhorn_method = "log",
    precision = "double",
    symmetric = TRUE,
    solver = fx$params$solver,
    structure_rank = max(nrow(fx$inputs$C1), nrow(fx$inputs$C2))
  )
  push(cbind(
    check = "fgw_lowrank_fullrank_obj",
    as.data.frame(.check_close(out_lowrank_full$fgw_dist, out_reference$fgw_dist, 1e-6, "fgw_dist"), stringsAsFactors = FALSE)
  ))
  push(cbind(
    check = "fgw_lowrank_fullrank_plan",
    as.data.frame(.check_close(out_lowrank_full$plan, out_reference$plan, 1e-5, "plan"), stringsAsFactors = FALSE)
  ))

  # 5) Multiset FGW mixed precision should remain close to double.
  set.seed(20260221L)
  make_subject <- function(n, d_struct = 3L, d_feat = 3L) {
    X <- matrix(rnorm(n * d_struct), n, d_struct)
    F <- matrix(rnorm(n * d_feat), n, d_feat)
    C <- as.matrix(dist(X))
    C <- C / max(C)
    list(C = C, F = F, w = rep(1 / n, n))
  }
  ms_subjects <- list(
    c(make_subject(18L), list(id = "s1")),
    c(make_subject(16L), list(id = "s2")),
    c(make_subject(15L), list(id = "s3"))
  )
  ms_template <- rfugw::multialign_make_template(
    subjects = ms_subjects,
    k = 14L,
    feature_normalization = "zscore",
    seed = 999L
  )
  ms_double <- rfugw::multialign_fit(
    subjects = ms_subjects,
    template_mode = "fixed",
    template = ms_template,
    method = "fgw_entropic",
    alpha = 0.5,
    epsilon = 0.05,
    precision = "double",
    max_iter = 80L,
    sinkhorn_max_iter = 250L,
    sinkhorn_tol = 1e-9,
    use_cpp_batch = TRUE,
    n_threads = 1L
  )
  ms_mixed <- rfugw::multialign_fit(
    subjects = ms_subjects,
    template_mode = "fixed",
    template = ms_template,
    method = "fgw_entropic",
    alpha = 0.5,
    epsilon = 0.05,
    precision = "mixed",
    max_iter = 80L,
    sinkhorn_max_iter = 250L,
    sinkhorn_tol = 1e-9,
    use_cpp_batch = TRUE,
    n_threads = 1L
  )
  push(cbind(
    check = "multiset_mixed_vs_double_obj",
    as.data.frame(.check_close(ms_mixed$objective_total, ms_double$objective_total, 5e-4, "objective_total"), stringsAsFactors = FALSE)
  ))
  mixed_plan <- do.call(c, lapply(ms_mixed$couplings, as.numeric))
  double_plan <- do.call(c, lapply(ms_double$couplings, as.numeric))
  push(cbind(
    check = "multiset_mixed_vs_double_plan",
    as.data.frame(.check_close(mixed_plan, double_plan, 2e-4, "couplings"), stringsAsFactors = FALSE)
  ))

  # 5b) Larger synthetic multiset auto-rank should remain close to dense and
  # respect permutation equivariance of template node order.
  set.seed(20260222L)
  ms_auto_subjects <- list(
    c(make_subject(420L), list(id = "a1")),
    c(make_subject(390L), list(id = "a2")),
    c(make_subject(360L), list(id = "a3"))
  )
  ms_auto_template <- rfugw::multialign_make_template(
    subjects = ms_auto_subjects,
    k = 400L,
    feature_normalization = "zscore",
    seed = 1234L
  )
  ms_auto_dense <- rfugw::multialign_fit(
    subjects = ms_auto_subjects,
    template_mode = "fixed",
    template = ms_auto_template,
    method = "fgw_entropic",
    alpha = 0.5,
    epsilon = 0.05,
    precision = "mixed",
    structure_rank = 0L,
    max_iter = 40L,
    sinkhorn_max_iter = 180L,
    sinkhorn_tol = 1e-9,
    use_cpp_batch = TRUE,
    n_threads = 1L
  )
  ms_auto <- rfugw::multialign_fit(
    subjects = ms_auto_subjects,
    template_mode = "fixed",
    template = ms_auto_template,
    method = "fgw_entropic",
    alpha = 0.5,
    epsilon = 0.05,
    precision = "mixed",
    structure_rank = "auto",
    max_iter = 40L,
    sinkhorn_max_iter = 180L,
    sinkhorn_tol = 1e-9,
    use_cpp_batch = TRUE,
    n_threads = 1L
  )
  push(data.frame(
    check = "multiset_auto_rank_adaptive",
    ok = isTRUE(ms_auto$structure_rank_adaptive),
    max_abs = if (isTRUE(ms_auto$structure_rank_adaptive)) 0 else 1,
    name = "structure_rank_adaptive",
    tol = 0,
    stringsAsFactors = FALSE
  ))
  push(cbind(
    check = "multiset_auto_vs_dense_obj",
    as.data.frame(.check_close(ms_auto$objective_total, ms_auto_dense$objective_total, 2e-3, "objective_total"), stringsAsFactors = FALSE)
  ))
  auto_plan <- do.call(c, lapply(ms_auto$couplings, as.numeric))
  dense_plan <- do.call(c, lapply(ms_auto_dense$couplings, as.numeric))
  push(cbind(
    check = "multiset_auto_vs_dense_plan",
    as.data.frame(.check_close(auto_plan, dense_plan, 3e-3, "couplings"), stringsAsFactors = FALSE)
  ))

  perm <- rev(seq_len(nrow(ms_auto_template$C)))
  inv_perm <- order(perm)
  tpl_perm <- list(
    C = ms_auto_template$C[perm, perm, drop = FALSE],
    F = ms_auto_template$F[perm, , drop = FALSE],
    w = ms_auto_template$w[perm],
    id = "template_perm"
  )
  ms_auto_perm <- rfugw::multialign_fit(
    subjects = ms_auto_subjects,
    template_mode = "fixed",
    template = tpl_perm,
    method = "fgw_entropic",
    alpha = 0.5,
    epsilon = 0.05,
    precision = "mixed",
    structure_rank = "auto",
    max_iter = 40L,
    sinkhorn_max_iter = 180L,
    sinkhorn_tol = 1e-9,
    use_cpp_batch = TRUE,
    n_threads = 1L
  )
  push(cbind(
    check = "multiset_auto_template_permute_obj",
    as.data.frame(.check_close(ms_auto_perm$objective_total, ms_auto$objective_total, 2e-6, "objective_total"), stringsAsFactors = FALSE)
  ))
  ids_auto <- names(ms_auto$couplings)
  auto_plan_perm_aligned <- do.call(c, lapply(ids_auto, function(id) {
    as.numeric(ms_auto_perm$couplings[[id]][, inv_perm, drop = FALSE])
  }))
  push(cbind(
    check = "multiset_auto_template_permute_plan",
    as.data.frame(.check_close(auto_plan_perm_aligned, auto_plan, 5e-4, "couplings"), stringsAsFactors = FALSE)
  ))

  # 6) FUGW KL vs POT fixture.
  fu <- .read_fixture("fugw_kl_sinkhorn_fixture.json", fixture_dir = fixture_dir)
  out_fugw <- rfugw::fugw_kl(
    Cx = fu$inputs$Cx,
    Cy = fu$inputs$Cy,
    wx = fu$inputs$wx,
    wy = fu$inputs$wy,
    reg_marginals = unlist(fu$params$reg_marginals),
    epsilon = fu$params$epsilon,
    alpha = fu$params$alpha,
    M = fu$inputs$M,
    max_iter = fu$params$max_iter,
    tol = fu$params$tol,
    max_iter_ot = fu$params$max_iter_ot,
    tol_ot = fu$params$tol_ot
  )
  push(cbind(
    check = "fugw_vs_pot_obj",
    as.data.frame(.check_close(out_fugw$fugw_cost, fu$outputs$fugw_cost, 5e-4, "fugw_cost"), stringsAsFactors = FALSE)
  ))
  push(cbind(
    check = "fugw_vs_pot_pi_samp",
    as.data.frame(.check_close(out_fugw$pi_samp, fu$outputs$pi_samp, 2e-4, "pi_samp"), stringsAsFactors = FALSE)
  ))

  # 7) Exact FGW CG vs POT exact fixture.
  fe <- .read_fixture("fgw_exact_square_fixture.json", fixture_dir = fixture_dir)
  out_exact <- rfugw::fgw_exact_cg(
    M = fe$inputs$M,
    C1 = fe$inputs$C1,
    C2 = fe$inputs$C2,
    p = fe$inputs$p,
    q = fe$inputs$q,
    alpha = fe$params$alpha,
    max_iter = fe$params$max_iter,
    tol_rel = fe$params$tol_rel,
    tol_abs = fe$params$tol_abs,
    lp_solver = "cpp_transport"
  )
  push(cbind(
    check = "fgw_exact_vs_pot_obj",
    as.data.frame(.check_close(out_exact$fgw_dist, fe$outputs$fgw_dist, 1e-5, "fgw_exact_obj"), stringsAsFactors = FALSE)
  ))
  push(cbind(
    check = "fgw_exact_vs_pot_plan",
    as.data.frame(.check_close(out_exact$plan, fe$outputs$plan, 1e-4, "fgw_exact_plan"), stringsAsFactors = FALSE)
  ))

  # 8) Sparse graph diffusion + coordinate sampled GW path.
  make_knn_similarity <- function(X, k = 5L) {
    D <- as.matrix(stats::dist(X))
    sigma <- stats::median(D[D > 0])
    if (!is.finite(sigma) || sigma <= 0) sigma <- 1
    W <- exp(-(D^2) / (2 * sigma^2))
    diag(W) <- 0
    n <- nrow(W)
    out <- matrix(0, nrow = n, ncol = n)
    for (i in seq_len(n)) {
      ord <- order(W[i, ], decreasing = TRUE)
      ord <- ord[ord != i]
      keep <- ord[seq_len(min(k, length(ord)))]
      out[i, keep] <- W[i, keep]
    }
    out <- pmax(out, t(out))
    diag(out) <- 0
    out
  }

  set.seed(20260224L)
  Xs <- matrix(rnorm(18L * 3L), 18L, 3L)
  Xt <- matrix(rnorm(17L * 3L), 17L, 3L)
  p_sparse <- rep(1 / nrow(Xs), nrow(Xs))
  q_sparse <- rep(1 / nrow(Xt), nrow(Xt))
  C1_sparse <- as.matrix(stats::dist(Xs))
  C2_sparse <- as.matrix(stats::dist(Xt))
  old_mixed <- Sys.getenv("RFUGW_SAMPLED_MIXED", unset = NA_character_)
  Sys.setenv(RFUGW_SAMPLED_MIXED = "0")
  on.exit({
    if (is.na(old_mixed)) {
      Sys.unsetenv("RFUGW_SAMPLED_MIXED")
    } else {
      Sys.setenv(RFUGW_SAMPLED_MIXED = old_mixed)
    }
  }, add = TRUE)

  Wd <- make_knn_similarity(Xs, k = 5L)
  if (requireNamespace("Matrix", quietly = TRUE)) {
    Ws <- Matrix::Matrix(Wd, sparse = TRUE)
    E_dense <- rfugw::graph_diffusion_coordinates(Wd, n_components = 4L, self_loop = 1e-6)
    E_sparse <- rfugw::graph_diffusion_coordinates(Ws, n_components = 4L, self_loop = 1e-6)
    push(cbind(
      check = "sparse_diffusion_dense_vs_sparse",
      as.data.frame(.check_close(as.matrix(stats::dist(E_dense)), as.matrix(stats::dist(E_sparse)), 1e-5, "diffusion_dist"), stringsAsFactors = FALSE)
    ))
  }

  out_coord_cpp <- rfugw::sampled_gromov_wasserstein_coords(
    X1 = Xs,
    X2 = Xt,
    p = p_sparse,
    q = q_sparse,
    nb_samples_grad = c(8L, 2L),
    epsilon = 0.05,
    max_iter = 45L,
    random_state = 19L,
    use_cpp = TRUE,
    sampling = "deterministic",
    log = TRUE
  )
  out_coord_r <- rfugw::sampled_gromov_wasserstein_coords(
    X1 = Xs,
    X2 = Xt,
    p = p_sparse,
    q = q_sparse,
    nb_samples_grad = c(8L, 2L),
    epsilon = 0.05,
    max_iter = 45L,
    random_state = 19L,
    use_cpp = FALSE,
    sampling = "deterministic",
    log = TRUE
  )
  val_cpp <- rfugw:::.gw_square_value(C1_sparse, C2_sparse, out_coord_cpp$plan, p_sparse, q_sparse, symmetric = TRUE)
  val_r <- rfugw:::.gw_square_value(C1_sparse, C2_sparse, out_coord_r$plan, p_sparse, q_sparse, symmetric = TRUE)
  push(cbind(
    check = "sampled_coords_cpp_vs_r_obj",
    as.data.frame(.check_close(val_cpp, val_r, 0.03, "gw_obj"), stringsAsFactors = FALSE)
  ))
  push(cbind(
    check = "sampled_coords_cpp_row_marginals",
    as.data.frame(.check_close(rowSums(out_coord_cpp$plan), p_sparse, 8e-4, "row_marg"), stringsAsFactors = FALSE)
  ))
  push(cbind(
    check = "sampled_coords_cpp_col_marginals",
    as.data.frame(.check_close(colSums(out_coord_cpp$plan), q_sparse, 8e-4, "col_marg"), stringsAsFactors = FALSE)
  ))

  if (requireNamespace("Matrix", quietly = TRUE)) {
    W1 <- Matrix::Matrix(make_knn_similarity(Xs, k = 5L), sparse = TRUE)
    W2 <- Matrix::Matrix(make_knn_similarity(Xt, k = 5L), sparse = TRUE)
    out_wrap <- rfugw::sampled_gw_from_graphs(
      W1 = W1,
      W2 = W2,
      n_components = 4L,
      self_loop = 1e-6,
      nb_samples_grad = c(8L, 2L),
      epsilon = 0.05,
      max_iter = 40L,
      random_state = 29L,
      use_cpp = FALSE,
      log = TRUE,
      return_embeddings = TRUE
    )
    out_manual <- rfugw::sampled_gromov_wasserstein_coords(
      X1 = out_wrap$X1_embed,
      X2 = out_wrap$X2_embed,
      nb_samples_grad = c(8L, 2L),
      epsilon = 0.05,
      max_iter = 40L,
      random_state = 29L,
      use_cpp = FALSE,
      log = TRUE
    )
    push(cbind(
      check = "sampled_graph_wrapper_vs_manual_plan",
      as.data.frame(.check_close(out_wrap$plan, out_manual$plan, 1e-10, "plan"), stringsAsFactors = FALSE)
    ))
  }

  tbl <- do.call(rbind, results)
  tbl$ok <- as.logical(tbl$ok)
  tbl$max_abs <- as.numeric(tbl$max_abs)
  tbl$tol <- as.numeric(tbl$tol)

  if (verbose) {
    cat("== Accuracy Gate ==\n")
    print(tbl[, c("check", "ok", "max_abs", "tol")], row.names = FALSE)
  }

  if (stop_on_fail && any(!tbl$ok)) {
    bad <- tbl[!tbl$ok, c("check", "max_abs", "tol"), drop = FALSE]
    stop(
      paste0(
        "Accuracy gate failed:\n",
        paste(
          sprintf(
            "- %s (max_abs=%g, tol=%g)",
            bad$check, bad$max_abs, bad$tol
          ),
          collapse = "\n"
        )
      ),
      call. = FALSE
    )
  }

  invisible(tbl)
}

if (sys.nframe() == 0L) {
  run_accuracy_gate(
    fixture_dir = Sys.getenv("RFUGW_FIXTURE_DIR", unset = "inst/extdata/fixtures"),
    stop_on_fail = TRUE,
    verbose = TRUE
  )
}
