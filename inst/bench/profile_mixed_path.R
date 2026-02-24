suppressPackageStartupMessages({
  rlib <- Sys.getenv("RFUGW_RLIB", unset = "")
  if (nzchar(rlib)) {
    .libPaths(c(rlib, .libPaths()))
  }
  library(rfugw)
})

args <- commandArgs(trailingOnly = TRUE)
iters <- if (length(args) >= 1) as.integer(args[[1]]) else 2L
out_csv <- if (length(args) >= 2) args[[2]] else "inst/bench/results/profile_mixed_path_latest.csv"
seed <- if (length(args) >= 3) as.integer(args[[3]]) else 20260222L
S <- if (length(args) >= 4) as.integer(args[[4]]) else 8L
n <- if (length(args) >= 5) as.integer(args[[5]]) else 600L
threads <- if (length(args) >= 6) as.integer(args[[6]]) else as.integer(Sys.getenv("OMP_NUM_THREADS", unset = "2"))
if (!is.finite(threads) || threads < 1L) threads <- 1L

Sys.setenv(
  OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1",
  OMP_NUM_THREADS = as.character(threads)
)

set.seed(seed)

make_subject_set <- function(n, d_struct = 3L, d_feat = 3L) {
  X <- matrix(rnorm(n * d_struct), n, d_struct)
  F <- matrix(rnorm(n * d_feat), n, d_feat)
  C <- as.matrix(dist(X))
  C <- C / max(C)
  list(C = C, F = F, w = rep(1 / n, n))
}

make_subjects <- function(S, n, jitter = 0.04) {
  out <- vector("list", S)
  for (s in seq_len(S)) {
    ns <- max(20L, as.integer(round(n * (1 + runif(1, -jitter, jitter)))))
    out[[s]] <- c(make_subject_set(ns), list(id = sprintf("s%02d", s)))
  }
  out
}

pairwise_euclidean <- function(X, Y) {
  x2 <- rowSums(X * X)
  y2 <- rowSums(Y * Y)
  D2 <- outer(x2, y2, "+") - 2 * tcrossprod(X, Y)
  D2[D2 < 0] <- 0
  sqrt(D2)
}

rescale01 <- function(M) {
  m <- max(M)
  if (is.finite(m) && m > 0) M / m else M
}

subjects <- make_subjects(S, n)
template <- rfugw::multialign_make_template(
  subjects = subjects,
  k = n,
  feature_normalization = "zscore",
  seed = seed + 1L
)

mk_batch_inputs <- function(subjects, template) {
  C2 <- rescale01(template$C)
  q <- template$w / sum(template$w)
  C1_list <- vector("list", length(subjects))
  p_list <- vector("list", length(subjects))
  M_list <- vector("list", length(subjects))
  for (i in seq_along(subjects)) {
    si <- subjects[[i]]
    C1 <- rescale01(si$C)
    M <- rescale01(pairwise_euclidean(si$F, template$F))
    C1_list[[i]] <- C1
    p_list[[i]] <- si$w / sum(si$w)
    M_list[[i]] <- M
  }
  list(M_list = M_list, C1_list = C1_list, p_list = p_list, C2 = C2, q = q)
}

batch_inputs <- mk_batch_inputs(subjects, template)

run_case <- function(name, fn) {
  times <- numeric(iters)
  objs <- numeric(iters)
  lowrank_counts <- rep(NA_real_, iters)
  kernel_means <- rep(NA_real_, iters)
  feature_means <- rep(NA_real_, iters)
  solve_means <- rep(NA_real_, iters)
  for (i in seq_len(iters)) {
    gc()
    t0 <- proc.time()[["elapsed"]]
    out <- fn()
    t1 <- proc.time()[["elapsed"]]
    times[[i]] <- (t1 - t0) * 1000
    objs[[i]] <- out$objective
    if (!is.null(out$n_lowrank)) {
      lowrank_counts[[i]] <- as.numeric(out$n_lowrank)
    }
    if (!is.null(out$mean_kernel_ms)) {
      kernel_means[[i]] <- as.numeric(out$mean_kernel_ms)
    }
    if (!is.null(out$mean_feature_ms)) {
      feature_means[[i]] <- as.numeric(out$mean_feature_ms)
    }
    if (!is.null(out$mean_solve_ms)) {
      solve_means[[i]] <- as.numeric(out$mean_solve_ms)
    }
  }
  km <- if (all(!is.finite(kernel_means))) NA_real_ else stats::median(kernel_means[is.finite(kernel_means)])
  fm <- if (all(!is.finite(feature_means))) NA_real_ else stats::median(feature_means[is.finite(feature_means)])
  sm <- if (all(!is.finite(solve_means))) NA_real_ else stats::median(solve_means[is.finite(solve_means)])
  fs <- if (is.finite(fm) && is.finite(sm) && (fm + sm) > 0) 100.0 * fm / (fm + sm) else NA_real_
  ss <- if (is.finite(fm) && is.finite(sm) && (fm + sm) > 0) 100.0 * sm / (fm + sm) else NA_real_
  data.frame(
    method = name,
    min_ms = min(times),
    median_ms = stats::median(times),
    objective = stats::median(objs),
    n_lowrank = if (all(!is.finite(lowrank_counts))) NA_real_ else stats::median(lowrank_counts[is.finite(lowrank_counts)]),
    kernel_ms = km,
    feature_ms = fm,
    solve_ms = sm,
    feature_share_pct = fs,
    solve_share_pct = ss,
    stringsAsFactors = FALSE
  )
}

fit_wrap <- function(...) {
  out <- rfugw::multialign_fit(
    subjects = subjects,
    template_mode = "fixed",
    template = template,
    method = "fgw_entropic",
    alpha = 0.5,
    epsilon = 0.05,
    feature_normalization = "zscore",
    max_iter = 60L,
    sinkhorn_max_iter = 220L,
    sinkhorn_tol = 1e-9,
    use_cpp_batch = TRUE,
    n_threads = threads,
    ...
  )
  list(
    objective = out$objective_total,
    n_lowrank = if (!is.null(out$batch_summary)) out$batch_summary$n_lowrank else NA_real_,
    mean_kernel_ms = if (!is.null(out$diagnostics$kernel_ms)) mean(out$diagnostics$kernel_ms, na.rm = TRUE) else NA_real_,
    mean_feature_ms = if (!is.null(out$diagnostics$feature_ms)) mean(out$diagnostics$feature_ms, na.rm = TRUE) else NA_real_,
    mean_solve_ms = if (!is.null(out$diagnostics$solve_ms)) mean(out$diagnostics$solve_ms, na.rm = TRUE) else NA_real_
  )
}

cpp_wrap <- function(approx_rank = 0L) {
  out <- rfugw:::cpp_fgw_entropic_square_batch(
    M_list = batch_inputs$M_list,
    C1_list = batch_inputs$C1_list,
    p_list = batch_inputs$p_list,
    C2 = batch_inputs$C2,
    q = batch_inputs$q,
    alpha = 0.5,
    epsilon = 0.05,
    max_iter = 60L,
    tol = 1e-9,
    sinkhorn_max_iter = 220L,
    sinkhorn_tol = 1e-9,
    symmetric = TRUE,
    use_ppa = FALSE,
    use_log_sinkhorn = FALSE,
    use_mixed_precision = TRUE,
    check_every = 10L,
    init_plan_list = vector("list", length(batch_inputs$M_list)),
    c1_A_scaled_list = list(),
    c1_Bt_list = list(),
    approx_rank = as.integer(approx_rank),
    n_threads = as.integer(threads)
  )
  list(
    objective = sum(unlist(out$fgw_dist), na.rm = TRUE),
    n_lowrank = if (!is.null(out$n_lowrank)) as.numeric(out$n_lowrank) else NA_real_,
    mean_kernel_ms = if (!is.null(out$kernel_ms)) mean(as.numeric(out$kernel_ms), na.rm = TRUE) else NA_real_
  )
}

rows <- list(
  run_case("multialign_batch_mixed_dense_split", function() fit_wrap(precision = "mixed", structure_rank = 0L, use_cpp_feature_fused = FALSE)),
  run_case("multialign_batch_mixed_dense", function() fit_wrap(precision = "mixed", structure_rank = 0L)),
  run_case("multialign_batch_mixed_lr64", function() fit_wrap(precision = "mixed", structure_rank = 64L)),
  run_case("multialign_batch_mixed_auto", function() fit_wrap(precision = "mixed", structure_rank = "auto", autotune = TRUE, autotune_level = "aggressive")),
  run_case("cpp_batch_mixed_dense", function() cpp_wrap(approx_rank = 0L)),
  run_case("cpp_batch_mixed_lr64", function() cpp_wrap(approx_rank = 64L))
)

out <- do.call(rbind, rows)
ref <- out$median_ms[out$method == "multialign_batch_mixed_dense"]
if (length(ref) == 1L && is.finite(ref) && ref > 0) {
  out$speedup_vs_mixed_dense <- ref / out$median_ms
} else {
  out$speedup_vs_mixed_dense <- NA_real_
}
out$S <- S
out$n <- n
out$threads <- threads

cat("== Mixed Path Profile ==\n")
print(out)

dir.create(dirname(out_csv), recursive = TRUE, showWarnings = FALSE)
write.csv(out, out_csv, row.names = FALSE)
cat("Wrote mixed-path profile CSV:", out_csv, "\n")
