pairwise_euclidean <- function(X, Y) {
  x2 <- rowSums(X * X)
  y2 <- rowSums(Y * Y)
  D2 <- outer(x2, y2, "+") - 2 * tcrossprod(X, Y)
  D2[D2 < 0] <- 0
  sqrt(D2)
}

make_set <- function(n, d_struct = 3L, d_feat = 2L, seed = 1L) {
  set.seed(seed)
  X <- matrix(rnorm(n * d_struct), n, d_struct)
  F <- matrix(rnorm(n * d_feat), n, d_feat)
  C <- as.matrix(dist(X))
  C <- C / max(C)
  list(C = C, F = F, w = rep(1 / n, n))
}

test_that("multialign_fit (FGW) matches direct pairwise solve", {
  s <- make_set(12L, seed = 11L)
  t <- make_set(9L, seed = 22L)
  M <- pairwise_euclidean(s$F, t$F)

  direct <- rfugw::fgw_entropic(
    M = M,
    C1 = s$C,
    C2 = t$C,
    p = s$w,
    q = t$w,
    alpha = 0.4,
    epsilon = 0.05,
    max_iter = 80L,
    tol = 1e-9,
    sinkhorn_max_iter = 300L,
    sinkhorn_tol = 1e-9,
    sinkhorn_method = "scaling",
    solver = "PGD",
    precision = "double",
    check_every = 1L
  )

  fit <- rfugw::multialign_fit(
    subjects = list(c(s, list(id = "s1"))),
    template = c(t, list(id = "tpl")),
    method = "fgw_entropic",
    alpha = 0.4,
    epsilon = 0.05,
    feature_normalization = "none",
    feature_cost_metric = "euclidean",
    structure_normalize = FALSE,
    feature_cost_normalize = FALSE,
    solver = "PGD",
    precision = "double",
    sinkhorn_method = "scaling",
    max_iter = 80L,
    tol = 1e-9,
    sinkhorn_max_iter = 300L,
    sinkhorn_tol = 1e-9,
    check_every = 1L
  )

  expect_equal(fit$objectives[[1]], direct$fgw_dist, tolerance = 1e-12)
  expect_equal(fit$couplings[[1]], direct$plan, tolerance = 1e-12)
  legacy_diagnostic_names <- c(
    "id", "n_nodes", "objective", "iterations", "error", "plan_mass",
    "kernel_ms", "feature_ms", "solve_ms", "lowrank_used", "c1_cache_used",
    "c2_cache_used", "square_cache_used"
  )
  expect_identical(
    head(names(fit$diagnostics), length(legacy_diagnostic_names)),
    legacy_diagnostic_names
  )
  diag_fields <- c(
    "status", "converged", "residual", "inner_status", "inner_converged",
    "inner_residual", "max_inner_residual", "inner_iterations"
  )
  expect_true(all(diag_fields %in% names(fit$diagnostics)))
  expect_identical(fit$diagnostics$status[[1]], direct$status)
  expect_identical(fit$diagnostics$converged[[1]], direct$converged)
  expect_identical(fit$diagnostics$inner_status[[1]], direct$inner_status)
  expect_identical(fit$diagnostics$inner_converged[[1]], direct$inner_converged)
  expect_equal(fit$diagnostics$inner_residual[[1]], direct$inner_residual)
  expect_equal(fit$diagnostics$max_inner_residual[[1]], direct$max_inner_residual)
})

test_that("multialign sampled-graph coarse init matches explicit warm-started FGW", {
  s <- make_set(15L, seed = 811L)
  t <- make_set(11L, seed = 812L)
  M <- pairwise_euclidean(s$F, t$F)
  coarse_ctrl <- list(
    n_components = 6L,
    diffusion_time = 1,
    self_loop = 1e-6,
    nb_samples_grad = c(6L, 2L),
    epsilon = 0.08,
    max_iter = 25L,
    use_cpp = FALSE,
    sampling = "deterministic"
  )

  warm <- rfugw::sampled_gw_from_graphs(
    W1 = rfugw:::.structure_cost_to_similarity(s$C),
    W2 = rfugw:::.structure_cost_to_similarity(t$C),
    p = s$w,
    q = t$w,
    n_components = coarse_ctrl$n_components,
    diffusion_time = coarse_ctrl$diffusion_time,
    self_loop = coarse_ctrl$self_loop,
    nb_samples_grad = coarse_ctrl$nb_samples_grad,
    epsilon = coarse_ctrl$epsilon,
    max_iter = coarse_ctrl$max_iter,
    use_cpp = coarse_ctrl$use_cpp,
    sampling = coarse_ctrl$sampling,
    log = FALSE,
    return_embeddings = FALSE
  )

  direct <- rfugw::fgw_entropic(
    M = M,
    C1 = s$C,
    C2 = t$C,
    p = s$w,
    q = t$w,
    alpha = 0.5,
    epsilon = 0.05,
    max_iter = 80L,
    tol = 1e-9,
    sinkhorn_max_iter = 300L,
    sinkhorn_tol = 1e-9,
    sinkhorn_method = "scaling",
    solver = "PGD",
    precision = "double",
    check_every = 1L,
    init_plan = warm
  )

  fit <- rfugw::multialign_fit(
    subjects = list(c(s, list(id = "s1"))),
    template = c(t, list(id = "tpl")),
    method = "fgw_entropic",
    alpha = 0.5,
    epsilon = 0.05,
    feature_normalization = "none",
    feature_cost_metric = "euclidean",
    structure_normalize = FALSE,
    feature_cost_normalize = FALSE,
    solver = "PGD",
    precision = "double",
    sinkhorn_method = "scaling",
    max_iter = 80L,
    tol = 1e-9,
    sinkhorn_max_iter = 300L,
    sinkhorn_tol = 1e-9,
    check_every = 1L,
    use_cpp_batch = FALSE,
    coarse_init = "sampled_graph",
    coarse_init_components = coarse_ctrl$n_components,
    coarse_init_diffusion_time = coarse_ctrl$diffusion_time,
    coarse_init_self_loop = coarse_ctrl$self_loop,
    coarse_init_nb_samples = coarse_ctrl$nb_samples_grad,
    coarse_init_epsilon = coarse_ctrl$epsilon,
    coarse_init_max_iter = coarse_ctrl$max_iter,
    coarse_init_use_cpp = coarse_ctrl$use_cpp,
    coarse_init_sampling = coarse_ctrl$sampling
  )

  expect_equal(fit$coarse_init, "sampled_graph")
  expect_equal(fit$objectives[[1]], direct$fgw_dist, tolerance = 1e-10)
  expect_equal(fit$couplings[[1]], direct$plan, tolerance = 1e-10)
})

test_that("multialign sampled-graph coarse init produces valid batch couplings", {
  s1 <- c(make_set(14L, seed = 821L), list(id = "s1"))
  s2 <- c(make_set(13L, seed = 822L), list(id = "s2"))
  s3 <- c(make_set(12L, seed = 823L), list(id = "s3"))
  tpl <- c(make_set(10L, seed = 824L), list(id = "tpl"))

  fit <- rfugw::multialign_fit(
    subjects = list(s1, s2, s3),
    template_mode = "fixed",
    template = tpl,
    method = "fgw_entropic",
    alpha = 0.5,
    epsilon = 0.05,
    precision = "mixed",
    max_iter = 60L,
    sinkhorn_max_iter = 220L,
    sinkhorn_tol = 1e-9,
    use_cpp_batch = TRUE,
    n_threads = 2L,
    coarse_init = "sampled_graph",
    coarse_init_components = 8L,
    coarse_init_nb_samples = c(6L, 2L),
    coarse_init_epsilon = 0.08,
    coarse_init_max_iter = 25L,
    coarse_init_use_cpp = TRUE,
    coarse_init_sampling = "deterministic"
  )

  expect_true(isTRUE(fit$used_cpp_batch))
  expect_equal(fit$coarse_init, "sampled_graph")
  expect_equal(sort(names(fit$couplings)), c("s1", "s2", "s3"))
  ref_weights <- list(s1 = s1$w, s2 = s2$w, s3 = s3$w)

  for (id in names(fit$couplings)) {
    G <- fit$couplings[[id]]
    expect_true(all(is.finite(G)))
    expect_true(all(G >= 0))
    expect_equal(rowSums(G), ref_weights[[id]], tolerance = 3e-5)
    expect_equal(colSums(G), tpl$w, tolerance = 3e-5)
  }
})

test_that("cpp feature-cost batch matches R feature-cost assembly", {
  set.seed(901L)
  F1_list <- list(
    matrix(rnorm(12L * 3L), 12L, 3L),
    matrix(rnorm(9L * 3L), 9L, 3L),
    matrix(rnorm(7L * 3L), 7L, 3L)
  )
  F2 <- matrix(rnorm(8L * 3L), 8L, 3L)

  out_eu <- rfugw:::cpp_feature_cost_batch(
    F1_list = F1_list,
    F2 = F2,
    use_euclidean = TRUE,
    rescale01 = TRUE,
    n_threads = 2L
  )
  ref_eu <- lapply(F1_list, function(F1) rfugw:::.cross_feature_cost(F1, F2, metric = "euclidean", rescale01 = TRUE))
  expect_gte(attr(out_eu, "used_threads"), 1L)
  expect_equal(out_eu[[1]], ref_eu[[1]], tolerance = 1e-12)
  expect_equal(out_eu[[2]], ref_eu[[2]], tolerance = 1e-12)
  expect_equal(out_eu[[3]], ref_eu[[3]], tolerance = 1e-12)

  out_sq <- rfugw:::cpp_feature_cost_batch(
    F1_list = F1_list,
    F2 = F2,
    use_euclidean = FALSE,
    rescale01 = FALSE,
    n_threads = 1L
  )
  ref_sq <- lapply(F1_list, function(F1) rfugw:::.cross_feature_cost(F1, F2, metric = "sqeuclidean", rescale01 = FALSE))
  expect_equal(out_sq[[1]], ref_sq[[1]], tolerance = 1e-12)
  expect_equal(out_sq[[2]], ref_sq[[2]], tolerance = 1e-12)
  expect_equal(out_sq[[3]], ref_sq[[3]], tolerance = 1e-12)
})

test_that("multialign_fit (FUGW) matches direct pairwise solve", {
  s <- make_set(10L, seed = 31L)
  t <- make_set(8L, seed = 41L)
  M <- pairwise_euclidean(s$F, t$F)

  direct <- rfugw::fugw_kl(
    Cx = s$C,
    Cy = t$C,
    wx = s$w,
    wy = t$w,
    reg_marginals = c(25, 10),
    epsilon = 1e-2,
    alpha = 0.6,
    M = M,
    max_iter = 40L,
    tol = 1e-8,
    max_iter_ot = 200L,
    tol_ot = 1e-8,
    check_every = 1L
  )

  fit <- rfugw::multialign_fit(
    subjects = list(c(s, list(id = "s1"))),
    template = c(t, list(id = "tpl")),
    method = "fugw_kl",
    alpha = 0.6,
    epsilon = 1e-2,
    feature_normalization = "none",
    feature_cost_metric = "euclidean",
    structure_normalize = FALSE,
    feature_cost_normalize = FALSE,
    max_iter = 40L,
    tol = 1e-8,
    reg_marginals = c(25, 10),
    max_iter_ot = 200L,
    tol_ot = 1e-8,
    check_every = 1L
  )

  expect_equal(fit$objectives[[1]], direct$fugw_cost, tolerance = 1e-12)
  expect_equal(fit$couplings[[1]], direct$pi_samp, tolerance = 1e-12)
})

test_that("multialign_fit returns one coupling per subject with valid marginals", {
  s1 <- c(make_set(11L, seed = 101L), list(id = "a"))
  s2 <- c(make_set(13L, seed = 102L), list(id = "b"))
  tpl <- c(make_set(9L, seed = 103L), list(id = "tpl"))

  fit <- rfugw::multialign_fit(
    subjects = list(s1, s2),
    template = tpl,
    method = "fgw_entropic",
    alpha = 0.5,
    epsilon = 0.05,
    max_iter = 60L,
    sinkhorn_max_iter = 250L,
    sinkhorn_tol = 1e-9,
    check_every = 1L
  )

  expect_equal(sort(names(fit$couplings)), c("a", "b"))
  expect_equal(length(fit$couplings), 2L)
  expect_equal(dim(fit$couplings$a), c(11L, 9L))
  expect_equal(dim(fit$couplings$b), c(13L, 9L))
  expect_equal(fit$diagnostics$id, c("a", "b"))

  expect_equal(rowSums(fit$couplings$a), s1$w, tolerance = 2e-6)
  expect_equal(rowSums(fit$couplings$b), s2$w, tolerance = 2e-6)
  expect_equal(colSums(fit$couplings$a), tpl$w, tolerance = 2e-6)
  expect_equal(colSums(fit$couplings$b), tpl$w, tolerance = 2e-6)
})

test_that("multialign_make_template and template=NULL flow work", {
  s1 <- c(make_set(10L, seed = 201L), list(id = "s1"))
  s2 <- c(make_set(12L, seed = 202L), list(id = "s2"))

  tpl <- rfugw::multialign_make_template(
    subjects = list(s1, s2),
    k = 7L,
    feature_normalization = "zscore",
    seed = 7L
  )
  expect_equal(dim(tpl$C), c(7L, 7L))
  expect_equal(dim(tpl$F), c(7L, 2L))
  expect_equal(sum(tpl$w), 1, tolerance = 1e-12)

  fit <- rfugw::multialign_fit(
    subjects = list(s1, s2),
    template = NULL,
    k_template = 7L,
    method = "fgw_entropic",
    alpha = 0.5,
    epsilon = 0.05,
    feature_normalization = "zscore",
    max_iter = 30L,
    sinkhorn_max_iter = 150L
  )
  expect_equal(dim(fit$template$C), c(7L, 7L))
  expect_equal(dim(fit$couplings$s1), c(10L, 7L))
  expect_equal(dim(fit$couplings$s2), c(12L, 7L))
})

test_that("projection helpers preserve expected dimensions", {
  P <- matrix(
    c(0.5, 0.2, 0.3,
      0.2, 0.7, 0.1,
      0.0, 0.1, 0.9,
      0.6, 0.0, 0.4),
    nrow = 4, byrow = TRUE
  )
  R <- matrix(
    c(1.0, 0.2, 0.3, 0.1,
      0.2, 1.1, 0.4, 0.0,
      0.3, 0.4, 0.9, 0.5,
      0.1, 0.0, 0.5, 1.3),
    nrow = 4, byrow = TRUE
  )
  Et <- matrix(
    c(0.0, 1.0,
      0.5, 0.5,
      1.0, 0.0),
    nrow = 3, byrow = TRUE
  )

  Pn <- rfugw::multialign_normalize_plan(P, mode = "row")
  expect_equal(rowSums(Pn), rep(1, 4), tolerance = 1e-12)

  Rt <- rfugw::multialign_project_matrix(P, R, normalize = "row")
  expect_equal(dim(Rt), c(3L, 3L))
  expect_equal(Rt, t(Rt), tolerance = 1e-12)

  Es <- rfugw::multialign_project_features(P, Et, normalize = "row")
  expect_equal(dim(Es), c(4L, 2L))
  expect_equal(Es[1, ], c(0.4, 0.6), tolerance = 1e-12)
})

test_that("multialign fixed C++ batch matches R-loop baseline", {
  s1 <- c(make_set(14L, seed = 301L), list(id = "s1"))
  s2 <- c(make_set(13L, seed = 302L), list(id = "s2"))
  s3 <- c(make_set(11L, seed = 303L), list(id = "s3"))
  tpl <- c(make_set(9L, seed = 304L), list(id = "tpl"))

  out_r <- rfugw::multialign_fit(
    subjects = list(s1, s2, s3),
    template_mode = "fixed",
    template = tpl,
    method = "fgw_entropic",
    alpha = 0.5,
    epsilon = 0.05,
    max_iter = 80L,
    sinkhorn_max_iter = 250L,
    sinkhorn_tol = 1e-9,
    use_cpp_batch = FALSE
  )
  out_cpp <- rfugw::multialign_fit(
    subjects = list(s1, s2, s3),
    template_mode = "fixed",
    template = tpl,
    method = "fgw_entropic",
    alpha = 0.5,
    epsilon = 0.05,
    max_iter = 80L,
    sinkhorn_max_iter = 250L,
    sinkhorn_tol = 1e-9,
    use_cpp_batch = TRUE,
    n_threads = 2L
  )

  expect_true(isTRUE(out_cpp$used_cpp_batch))
  expect_gte(out_cpp$used_threads, 1L)
  expect_equal(out_cpp$objective_total, out_r$objective_total, tolerance = 1e-10)
  expect_equal(out_cpp$couplings$s1, out_r$couplings$s1, tolerance = 1e-10)
  expect_equal(out_cpp$couplings$s2, out_r$couplings$s2, tolerance = 1e-10)
  expect_equal(out_cpp$couplings$s3, out_r$couplings$s3, tolerance = 1e-10)
  diag_fields <- c(
    "status", "converged", "residual", "inner_status", "inner_converged",
    "inner_residual", "max_inner_residual", "inner_iterations"
  )
  expect_true(all(diag_fields %in% names(out_r$diagnostics)))
  expect_true(all(diag_fields %in% names(out_cpp$diagnostics)))
  expect_true(all(!is.na(out_r$diagnostics$status)))
  expect_true(all(out_cpp$diagnostics$inner_status == "unavailable"))
  expect_true(all(is.na(out_cpp$diagnostics$inner_converged)))
})

test_that("multialign fused feature+solve batch matches non-fused batch (double)", {
  s1 <- c(make_set(14L, seed = 341L), list(id = "s1"))
  s2 <- c(make_set(13L, seed = 342L), list(id = "s2"))
  s3 <- c(make_set(12L, seed = 343L), list(id = "s3"))
  tpl <- c(make_set(10L, seed = 344L), list(id = "tpl"))

  out_fused <- rfugw::multialign_fit(
    subjects = list(s1, s2, s3),
    template_mode = "fixed",
    template = tpl,
    method = "fgw_entropic",
    alpha = 0.5,
    epsilon = 0.05,
    precision = "double",
    max_iter = 80L,
    sinkhorn_max_iter = 250L,
    sinkhorn_tol = 1e-9,
    use_cpp_batch = TRUE,
    use_cpp_feature_fused = TRUE,
    n_threads = 2L
  )
  out_split <- rfugw::multialign_fit(
    subjects = list(s1, s2, s3),
    template_mode = "fixed",
    template = tpl,
    method = "fgw_entropic",
    alpha = 0.5,
    epsilon = 0.05,
    precision = "double",
    max_iter = 80L,
    sinkhorn_max_iter = 250L,
    sinkhorn_tol = 1e-9,
    use_cpp_batch = TRUE,
    use_cpp_feature_fused = FALSE,
    n_threads = 2L
  )

  expect_true(isTRUE(out_fused$used_cpp_feature_fused))
  expect_equal(out_fused$objective_total, out_split$objective_total, tolerance = 1e-10)
  expect_equal(out_fused$couplings$s1, out_split$couplings$s1, tolerance = 1e-10)
  expect_equal(out_fused$couplings$s2, out_split$couplings$s2, tolerance = 1e-10)
  expect_equal(out_fused$couplings$s3, out_split$couplings$s3, tolerance = 1e-10)
})

test_that("multialign fused feature+solve batch matches non-fused batch (mixed)", {
  s1 <- c(make_set(14L, seed = 351L), list(id = "s1"))
  s2 <- c(make_set(13L, seed = 352L), list(id = "s2"))
  s3 <- c(make_set(12L, seed = 353L), list(id = "s3"))
  tpl <- c(make_set(10L, seed = 354L), list(id = "tpl"))

  out_fused <- rfugw::multialign_fit(
    subjects = list(s1, s2, s3),
    template_mode = "fixed",
    template = tpl,
    method = "fgw_entropic",
    alpha = 0.5,
    epsilon = 0.05,
    precision = "mixed",
    max_iter = 80L,
    sinkhorn_max_iter = 250L,
    sinkhorn_tol = 1e-9,
    use_cpp_batch = TRUE,
    use_cpp_feature_fused = TRUE,
    n_threads = 2L
  )
  out_split <- rfugw::multialign_fit(
    subjects = list(s1, s2, s3),
    template_mode = "fixed",
    template = tpl,
    method = "fgw_entropic",
    alpha = 0.5,
    epsilon = 0.05,
    precision = "mixed",
    max_iter = 80L,
    sinkhorn_max_iter = 250L,
    sinkhorn_tol = 1e-9,
    use_cpp_batch = TRUE,
    use_cpp_feature_fused = FALSE,
    n_threads = 2L
  )

  expect_true(isTRUE(out_fused$used_cpp_feature_fused))
  expect_equal(out_fused$objective_total, out_split$objective_total, tolerance = 5e-4)
  expect_equal(out_fused$couplings$s1, out_split$couplings$s1, tolerance = 2e-4)
  expect_equal(out_fused$couplings$s2, out_split$couplings$s2, tolerance = 2e-4)
  expect_equal(out_fused$couplings$s3, out_split$couplings$s3, tolerance = 2e-4)
})

test_that("multialign batch diagnostics expose kernel/profile counters", {
  s1 <- c(make_set(12L, seed = 331L), list(id = "s1"))
  s2 <- c(make_set(11L, seed = 332L), list(id = "s2"))
  s3 <- c(make_set(10L, seed = 333L), list(id = "s3"))
  tpl <- c(make_set(9L, seed = 334L), list(id = "tpl"))

  out <- rfugw::multialign_fit(
    subjects = list(s1, s2, s3),
    template_mode = "fixed",
    template = tpl,
    method = "fgw_entropic",
    alpha = 0.5,
    epsilon = 0.05,
    max_iter = 40L,
    sinkhorn_max_iter = 180L,
    sinkhorn_tol = 1e-9,
    use_cpp_batch = TRUE,
    n_threads = 2L
  )

  expect_true(isTRUE(out$used_cpp_batch))
  expect_true(is.data.frame(out$diagnostics))
  expect_true(all(c("kernel_ms", "lowrank_used", "c1_cache_used", "c2_cache_used", "square_cache_used") %in% names(out$diagnostics)))
  expect_true(is.list(out$batch_summary))
  expect_true(all(c("n_lowrank", "n_c1_cache", "n_c2_cache", "n_square_cache") %in% names(out$batch_summary)))
  expect_equal(out$batch_summary$n_lowrank, 0L)
  expect_equal(out$batch_summary$n_square_cache, 3L)
})

test_that("multialign fixed C++ batch low-rank cache matches R-loop baseline", {
  s1 <- c(make_set(14L, seed = 601L), list(id = "s1"))
  s2 <- c(make_set(13L, seed = 602L), list(id = "s2"))
  s3 <- c(make_set(11L, seed = 603L), list(id = "s3"))
  tpl <- c(make_set(9L, seed = 604L), list(id = "tpl"))

  out_r <- rfugw::multialign_fit(
    subjects = list(s1, s2, s3),
    template_mode = "fixed",
    template = tpl,
    method = "fgw_entropic",
    alpha = 0.5,
    epsilon = 0.05,
    precision = "double",
    structure_rank = 6L,
    max_iter = 80L,
    sinkhorn_max_iter = 250L,
    sinkhorn_tol = 1e-9,
    use_cpp_batch = FALSE
  )
  out_cpp <- rfugw::multialign_fit(
    subjects = list(s1, s2, s3),
    template_mode = "fixed",
    template = tpl,
    method = "fgw_entropic",
    alpha = 0.5,
    epsilon = 0.05,
    precision = "double",
    structure_rank = 6L,
    max_iter = 80L,
    sinkhorn_max_iter = 250L,
    sinkhorn_tol = 1e-9,
    use_cpp_batch = TRUE,
    n_threads = 2L
  )

  expect_true(isTRUE(out_cpp$used_cpp_batch))
  expect_equal(out_cpp$objective_total, out_r$objective_total, tolerance = 1e-10)
  expect_equal(out_cpp$couplings$s1, out_r$couplings$s1, tolerance = 1e-10)
  expect_equal(out_cpp$couplings$s2, out_r$couplings$s2, tolerance = 1e-10)
  expect_equal(out_cpp$couplings$s3, out_r$couplings$s3, tolerance = 1e-10)
})

test_that("multialign fixed C++ batch (FUGW) matches R-loop baseline", {
  s1 <- c(make_set(12L, seed = 611L), list(id = "s1"))
  s2 <- c(make_set(11L, seed = 612L), list(id = "s2"))
  s3 <- c(make_set(10L, seed = 613L), list(id = "s3"))
  tpl <- c(make_set(9L, seed = 614L), list(id = "tpl"))

  out_r <- rfugw::multialign_fit(
    subjects = list(s1, s2, s3),
    template_mode = "fixed",
    template = tpl,
    method = "fugw_kl",
    alpha = 0.5,
    epsilon = 1e-2,
    reg_marginals = c(20, 10),
    max_iter = 30L,
    max_iter_ot = 180L,
    tol = 1e-8,
    tol_ot = 1e-8,
    use_cpp_batch = FALSE
  )
  out_cpp <- rfugw::multialign_fit(
    subjects = list(s1, s2, s3),
    template_mode = "fixed",
    template = tpl,
    method = "fugw_kl",
    alpha = 0.5,
    epsilon = 1e-2,
    reg_marginals = c(20, 10),
    max_iter = 30L,
    max_iter_ot = 180L,
    tol = 1e-8,
    tol_ot = 1e-8,
    use_cpp_batch = TRUE,
    n_threads = 2L
  )

  expect_true(isTRUE(out_cpp$used_cpp_batch))
  expect_equal(out_cpp$objective_total, out_r$objective_total, tolerance = 1e-10)
  expect_equal(out_cpp$couplings$s1, out_r$couplings$s1, tolerance = 1e-10)
  expect_equal(out_cpp$couplings$s2, out_r$couplings$s2, tolerance = 1e-10)
  expect_equal(out_cpp$couplings$s3, out_r$couplings$s3, tolerance = 1e-10)
})

test_that("multialign mixed precision remains close to double", {
  s1 <- c(make_set(14L, seed = 701L), list(id = "s1"))
  s2 <- c(make_set(13L, seed = 702L), list(id = "s2"))
  s3 <- c(make_set(12L, seed = 703L), list(id = "s3"))
  tpl <- c(make_set(10L, seed = 704L), list(id = "tpl"))

  out_double <- rfugw::multialign_fit(
    subjects = list(s1, s2, s3),
    template_mode = "fixed",
    template = tpl,
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
  out_mixed <- rfugw::multialign_fit(
    subjects = list(s1, s2, s3),
    template_mode = "fixed",
    template = tpl,
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

  expect_equal(out_mixed$objective_total, out_double$objective_total, tolerance = 5e-4)
  expect_equal(out_mixed$couplings$s1, out_double$couplings$s1, tolerance = 2e-4)
  expect_equal(out_mixed$couplings$s2, out_double$couplings$s2, tolerance = 2e-4)
  expect_equal(out_mixed$couplings$s3, out_double$couplings$s3, tolerance = 2e-4)
})

test_that("multialign fixed C++ batch low-rank mixed cache matches non-batch mixed", {
  s1 <- c(make_set(14L, seed = 705L), list(id = "s1"))
  s2 <- c(make_set(13L, seed = 706L), list(id = "s2"))
  s3 <- c(make_set(12L, seed = 707L), list(id = "s3"))
  tpl <- c(make_set(10L, seed = 708L), list(id = "tpl"))

  out_r <- rfugw::multialign_fit(
    subjects = list(s1, s2, s3),
    template_mode = "fixed",
    template = tpl,
    method = "fgw_entropic",
    alpha = 0.5,
    epsilon = 0.05,
    precision = "mixed",
    structure_rank = 6L,
    max_iter = 80L,
    sinkhorn_max_iter = 250L,
    sinkhorn_tol = 1e-9,
    use_cpp_batch = FALSE
  )
  out_cpp <- rfugw::multialign_fit(
    subjects = list(s1, s2, s3),
    template_mode = "fixed",
    template = tpl,
    method = "fgw_entropic",
    alpha = 0.5,
    epsilon = 0.05,
    precision = "mixed",
    structure_rank = 6L,
    max_iter = 80L,
    sinkhorn_max_iter = 250L,
    sinkhorn_tol = 1e-9,
    use_cpp_batch = TRUE,
    n_threads = 2L
  )

  expect_true(isTRUE(out_cpp$used_cpp_batch))
  expect_equal(out_cpp$objective_total, out_r$objective_total, tolerance = 5e-4)
  expect_equal(out_cpp$couplings$s1, out_r$couplings$s1, tolerance = 2e-4)
  expect_equal(out_cpp$couplings$s2, out_r$couplings$s2, tolerance = 2e-4)
  expect_equal(out_cpp$couplings$s3, out_r$couplings$s3, tolerance = 2e-4)
})

test_that("multialign low-rank full-rank fallback matches dense FGW", {
  s1 <- c(make_set(12L, seed = 711L), list(id = "s1"))
  s2 <- c(make_set(11L, seed = 712L), list(id = "s2"))
  tpl <- c(make_set(9L, seed = 713L), list(id = "tpl"))

  out_dense <- rfugw::multialign_fit(
    subjects = list(s1, s2),
    template_mode = "fixed",
    template = tpl,
    method = "fgw_entropic",
    alpha = 0.5,
    epsilon = 0.05,
    precision = "double",
    structure_rank = 0L,
    max_iter = 80L,
    sinkhorn_max_iter = 250L,
    sinkhorn_tol = 1e-9,
    use_cpp_batch = TRUE,
    n_threads = 1L
  )
  out_full <- rfugw::multialign_fit(
    subjects = list(s1, s2),
    template_mode = "fixed",
    template = tpl,
    method = "fgw_entropic",
    alpha = 0.5,
    epsilon = 0.05,
    precision = "double",
    structure_rank = max(nrow(s1$C), nrow(tpl$C)),
    max_iter = 80L,
    sinkhorn_max_iter = 250L,
    sinkhorn_tol = 1e-9,
    use_cpp_batch = TRUE,
    n_threads = 1L
  )

  expect_equal(out_full$objective_total, out_dense$objective_total, tolerance = 1e-6)
  expect_equal(out_full$couplings$s1, out_dense$couplings$s1, tolerance = 1e-5)
  expect_equal(out_full$couplings$s2, out_dense$couplings$s2, tolerance = 1e-5)
})

test_that("multialign autotune exposes tuned controls", {
  s1 <- c(make_set(22L, seed = 721L), list(id = "s1"))
  s2 <- c(make_set(20L, seed = 722L), list(id = "s2"))
  s3 <- c(make_set(19L, seed = 723L), list(id = "s3"))
  tpl <- c(make_set(18L, seed = 724L), list(id = "tpl"))

  out <- rfugw::multialign_fit(
    subjects = list(s1, s2, s3),
    template_mode = "fixed",
    template = tpl,
    method = "fgw_entropic",
    alpha = 0.5,
    epsilon = 0.05,
    precision = "mixed",
    autotune = TRUE,
    autotune_level = "aggressive",
    max_iter = 200L,
    sinkhorn_max_iter = 600L,
    check_every = 1L,
    use_cpp_batch = TRUE,
    n_threads = 1L
  )

  expect_true(isTRUE(out$autotune))
  expect_true(is.list(out$autotuned_controls))
  expect_gte(out$autotuned_controls$check_every, 1L)
  expect_lte(out$autotuned_controls$max_iter, 200L)
  expect_lte(out$autotuned_controls$sinkhorn_max_iter, 600L)
  expect_true("structure_rank" %in% names(out$autotuned_controls))
})

test_that("multialign structure_rank auto mode resolves rank and metadata", {
  s1 <- c(make_set(22L, seed = 731L), list(id = "s1"))
  s2 <- c(make_set(20L, seed = 732L), list(id = "s2"))
  s3 <- c(make_set(19L, seed = 733L), list(id = "s3"))
  tpl <- c(make_set(18L, seed = 734L), list(id = "tpl"))

  out <- rfugw::multialign_fit(
    subjects = list(s1, s2, s3),
    template_mode = "fixed",
    template = tpl,
    method = "fgw_entropic",
    alpha = 0.5,
    epsilon = 0.05,
    precision = "mixed",
    structure_rank = "auto",
    max_iter = 80L,
    sinkhorn_max_iter = 300L,
    use_cpp_batch = TRUE,
    n_threads = 1L
  )

  expect_true(isTRUE(out$structure_rank_adaptive))
  expect_equal(out$structure_rank_requested, "auto")
  expect_equal(out$structure_rank, 0L)
})

test_that("structure-rank autotune helper is size-aware", {
  small <- rfugw:::.autotune_structure_rank(
    ns_source = c(180L, 220L, 200L),
    n_template = 210L,
    precision = "mixed",
    level = "aggressive"
  )
  large <- rfugw:::.autotune_structure_rank(
    ns_source = c(1800L, 2200L, 2000L),
    n_template = 2100L,
    precision = "mixed",
    level = "aggressive"
  )

  expect_equal(small$structure_rank, 0L)
  expect_gt(large$structure_rank, 0L)
  expect_lte(large$structure_rank, large$rank_cap)
})

test_that("kNN sparsification helper preserves symmetry and diagonal", {
  set.seed(811L)
  X <- matrix(rnorm(12L * 3L), 12L, 3L)
  C <- as.matrix(dist(X))
  out <- rfugw:::.knn_sparsify_structure(C, k = 3L)
  expect_equal(dim(out), dim(C))
  expect_equal(out, t(out), tolerance = 0)
  expect_equal(diag(out), diag(C), tolerance = 0)
})

test_that("multialign learned-template mode returns history and updates template", {
  s1 <- c(make_set(12L, seed = 401L), list(id = "s1"))
  s2 <- c(make_set(10L, seed = 402L), list(id = "s2"))
  s3 <- c(make_set(11L, seed = 403L), list(id = "s3"))
  tpl0 <- rfugw::multialign_make_template(
    subjects = list(s1, s2, s3),
    k = 8L,
    feature_normalization = "zscore",
    seed = 123L
  )

  out <- rfugw::multialign_fit(
    subjects = list(s1, s2, s3),
    template_mode = "learned",
    template = tpl0,
    method = "fgw_entropic",
    alpha = 0.5,
    epsilon = 0.05,
    feature_normalization = "zscore",
    max_iter = 50L,
    sinkhorn_max_iter = 200L,
    template_max_iter = 3L,
    template_tol = 1e-6,
    use_cpp_batch = TRUE,
    n_threads = 2L
  )

  expect_equal(out$template_mode, "learned")
  expect_true(is.data.frame(out$template_history))
  expect_equal(nrow(out$template_history), out$outer_iterations)
  expect_true(all(c("outer_iter", "objective_total", "delta_objective", "delta_template") %in% names(out$template_history)))
  expect_equal(dim(out$template$C), c(8L, 8L))
  expect_equal(sum(out$template$w), 1, tolerance = 1e-10)

  delta_c <- max(abs(out$template$C - tpl0$C))
  expect_gt(delta_c, 1e-12)
})

test_that("multialign learned-template mode works for fugw_kl backend", {
  s1 <- c(make_set(9L, seed = 501L), list(id = "s1"))
  s2 <- c(make_set(8L, seed = 502L), list(id = "s2"))
  tpl0 <- rfugw::multialign_make_template(
    subjects = list(s1, s2),
    k = 6L,
    feature_normalization = "zscore",
    seed = 77L
  )

  out <- rfugw::multialign_fit(
    subjects = list(s1, s2),
    template_mode = "learned",
    template = tpl0,
    method = "fugw_kl",
    alpha = 0.5,
    epsilon = 1e-2,
    reg_marginals = c(20, 10),
    feature_normalization = "zscore",
    max_iter = 20L,
    max_iter_ot = 120L,
    template_max_iter = 2L,
    template_tol = 1e-6
  )

  expect_equal(out$template_mode, "learned")
  expect_equal(out$method, "fugw_kl")
  expect_equal(length(out$couplings), 2L)
  expect_equal(nrow(out$template_history), out$outer_iterations)
})
