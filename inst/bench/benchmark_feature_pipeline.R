suppressPackageStartupMessages({
  rlib <- Sys.getenv("RFUGW_RLIB", unset = "")
  if (nzchar(rlib)) {
    .libPaths(c(rlib, .libPaths()))
  }
  library(rfugw)
})

args <- commandArgs(trailingOnly = TRUE)
iters <- if (length(args) >= 1) as.integer(args[[1]]) else 3L
out_csv <- if (length(args) >= 2) args[[2]] else "inst/bench/results/feature_pipeline_latest.csv"
seed <- if (length(args) >= 3) as.integer(args[[3]]) else 20260222L
S <- if (length(args) >= 4) as.integer(args[[4]]) else 8L
n <- if (length(args) >= 5) as.integer(args[[5]]) else 600L
threads <- if (length(args) >= 6) as.integer(args[[6]]) else as.integer(Sys.getenv("OMP_NUM_THREADS", unset = "2"))
precision <- if (length(args) >= 7) tolower(args[[7]]) else "mixed"
if (!precision %in% c("mixed", "double")) {
  stop("`precision` must be 'mixed' or 'double'.", call. = FALSE)
}
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

time_many <- function(fn, iters = 3L) {
  times <- numeric(iters)
  for (i in seq_len(iters)) {
    gc()
    t0 <- proc.time()[["elapsed"]]
    invisible(fn())
    t1 <- proc.time()[["elapsed"]]
    times[[i]] <- (t1 - t0) * 1000
  }
  stats::median(times)
}

subjects <- make_subjects(S, n)
template <- rfugw::multialign_make_template(
  subjects = subjects,
  k = n,
  feature_normalization = "zscore",
  seed = seed + 1L
)

bundle <- rfugw:::.normalize_feature_bundle(c(lapply(subjects, function(s) s$F), list(template$F)), mode = "zscore")
for (i in seq_along(subjects)) {
  subjects[[i]]$F <- bundle[[i]]
}
template$F <- bundle[[length(bundle)]]

C1_list <- lapply(subjects, function(s) rfugw:::.maybe_rescale01(s$C, enabled = TRUE))
p_list <- lapply(subjects, function(s) s$w / sum(s$w))
F1_list <- lapply(subjects, function(s) s$F)
F2 <- template$F
C2 <- rfugw:::.maybe_rescale01(template$C, enabled = TRUE)
q <- template$w / sum(template$w)

solve_batch <- function(M_list) {
  rfugw:::cpp_fgw_entropic_square_batch(
    M_list = M_list,
    C1_list = C1_list,
    p_list = p_list,
    C2 = C2,
    q = q,
    alpha = 0.5,
    epsilon = 0.05,
    max_iter = 60L,
    tol = 1e-9,
    sinkhorn_max_iter = 220L,
    sinkhorn_tol = 1e-9,
    symmetric = TRUE,
    use_ppa = FALSE,
    use_log_sinkhorn = FALSE,
    use_mixed_precision = identical(precision, "mixed"),
    check_every = 8L,
    init_plan_list = vector("list", length(M_list)),
    c1_A_scaled_list = list(),
    c1_Bt_list = list(),
    approx_rank = 0L,
    n_threads = as.integer(threads)
  )
}

solve_fused <- function() {
  rfugw:::cpp_fgw_entropic_square_batch_features(
    F1_list = F1_list,
    C1_list = C1_list,
    p_list = p_list,
    F2 = F2,
    C2 = C2,
    q = q,
    use_euclidean = TRUE,
    rescale01 = TRUE,
    alpha = 0.5,
    epsilon = 0.05,
    max_iter = 60L,
    tol = 1e-9,
    sinkhorn_max_iter = 220L,
    sinkhorn_tol = 1e-9,
    symmetric = TRUE,
    use_ppa = FALSE,
    use_log_sinkhorn = FALSE,
    use_mixed_precision = identical(precision, "mixed"),
    check_every = 8L,
    init_plan_list = vector("list", length(F1_list)),
    c1_A_scaled_list = list(),
    c1_Bt_list = list(),
    approx_rank = 0L,
    n_threads = as.integer(threads)
  )
}

feature_r_ms <- time_many(function() {
  lapply(F1_list, function(F1) rfugw:::.cross_feature_cost(F1, F2, metric = "euclidean", rescale01 = TRUE))
}, iters = iters)

feature_cpp_ms <- time_many(function() {
  rfugw:::cpp_feature_cost_batch(
    F1_list = F1_list,
    F2 = F2,
    use_euclidean = TRUE,
    rescale01 = TRUE,
    n_threads = as.integer(threads)
  )
}, iters = iters)

M_cpp <- rfugw:::cpp_feature_cost_batch(
  F1_list = F1_list,
  F2 = F2,
  use_euclidean = TRUE,
  rescale01 = TRUE,
  n_threads = as.integer(threads)
)

solve_only_ms <- time_many(function() {
  solve_batch(M_cpp)
}, iters = iters)

nonfused_total_ms <- time_many(function() {
  M <- rfugw:::cpp_feature_cost_batch(
    F1_list = F1_list,
    F2 = F2,
    use_euclidean = TRUE,
    rescale01 = TRUE,
    n_threads = as.integer(threads)
  )
  solve_batch(M)
}, iters = iters)

fused_total_ms <- time_many(function() {
  solve_fused()
}, iters = iters)

out_nonfused <- solve_batch(M_cpp)
out_fused <- solve_fused()

obj_diff <- max(abs(as.numeric(out_nonfused$fgw_dist) - as.numeric(out_fused$fgw_dist)))
plan_diff <- max(vapply(seq_along(out_nonfused$plans), function(i) {
  max(abs(out_nonfused$plans[[i]] - out_fused$plans[[i]]))
}, numeric(1)))

tol_obj <- if (identical(precision, "mixed")) 5e-4 else 1e-10
tol_plan <- if (identical(precision, "mixed")) 2e-4 else 1e-10
if (!is.finite(obj_diff) || obj_diff > tol_obj) {
  stop(sprintf("Objective drift too large for fused path: %.3e > %.3e", obj_diff, tol_obj), call. = FALSE)
}
if (!is.finite(plan_diff) || plan_diff > tol_plan) {
  stop(sprintf("Plan drift too large for fused path: %.3e > %.3e", plan_diff, tol_plan), call. = FALSE)
}

nonfused_components_ms <- feature_cpp_ms + solve_only_ms
feature_share_nonfused_pct <- 100 * feature_cpp_ms / nonfused_components_ms
solve_share_nonfused_pct <- 100 * solve_only_ms / nonfused_components_ms

summary_row <- data.frame(
  S = S,
  n = n,
  threads = threads,
  precision = precision,
  feature_r_ms = feature_r_ms,
  feature_cpp_ms = feature_cpp_ms,
  solve_only_ms = solve_only_ms,
  nonfused_total_ms = nonfused_total_ms,
  fused_total_ms = fused_total_ms,
  feature_cpp_speedup_vs_r = feature_r_ms / feature_cpp_ms,
  fused_speedup_vs_nonfused = nonfused_total_ms / fused_total_ms,
  feature_share_nonfused_pct = feature_share_nonfused_pct,
  solve_share_nonfused_pct = solve_share_nonfused_pct,
  fused_mean_feature_ms = mean(as.numeric(out_fused$feature_ms), na.rm = TRUE),
  fused_mean_solve_ms = mean(as.numeric(out_fused$solve_ms), na.rm = TRUE),
  objective_diff = obj_diff,
  plan_diff = plan_diff,
  feature_threads_used = as.integer(attr(M_cpp, "used_threads")),
  fused_threads_used = as.integer(out_fused$used_threads),
  stringsAsFactors = FALSE
)

cat("== Feature Pipeline Benchmark ==\n")
print(summary_row)

dir.create(dirname(out_csv), recursive = TRUE, showWarnings = FALSE)
write.csv(summary_row, out_csv, row.names = FALSE)
cat("Wrote feature pipeline benchmark CSV:", out_csv, "\n")
