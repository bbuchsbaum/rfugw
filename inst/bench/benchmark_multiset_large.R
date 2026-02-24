suppressPackageStartupMessages({
  rlib <- Sys.getenv("RFUGW_RLIB", unset = "")
  if (nzchar(rlib)) {
    .libPaths(c(rlib, .libPaths()))
  }
  library(rfugw)
})

if (Sys.getenv("RFUGW_SKIP_ACCURACY_GATE", unset = "0") != "1") {
  gate_script <- if (file.exists("inst/bench/accuracy_gate.R")) {
    "inst/bench/accuracy_gate.R"
  } else {
    "rfugw/inst/bench/accuracy_gate.R"
  }
  source(gate_script)
  run_accuracy_gate(
    fixture_dir = Sys.getenv("RFUGW_FIXTURE_DIR", unset = "inst/extdata/fixtures"),
    stop_on_fail = TRUE,
    verbose = TRUE
  )
}

args <- commandArgs(trailingOnly = TRUE)
iters <- if (length(args) >= 1) as.integer(args[[1]]) else 1L
out_csv <- if (length(args) >= 2) args[[2]] else "inst/bench/results/benchmark_multiset_large_latest.csv"
seed <- if (length(args) >= 3) as.integer(args[[3]]) else 42L
threads_arg <- if (length(args) >= 4) args[[4]] else "1 2 4 8"
thread_grid <- as.integer(strsplit(threads_arg, "[,[:space:]]+")[[1]])
thread_grid <- thread_grid[is.finite(thread_grid) & thread_grid >= 1L]
if (length(thread_grid) == 0L) thread_grid <- 1L
n_arg <- if (length(args) >= 5) args[[5]] else "500 1000 2000"
n_grid <- as.integer(strsplit(n_arg, "[,[:space:]]+")[[1]])
n_grid <- n_grid[is.finite(n_grid) & n_grid >= 50L]
if (length(n_grid) == 0L) n_grid <- c(500L, 1000L, 2000L)

Sys.setenv(
  OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1"
)

set.seed(seed)

make_subject_set <- function(n, d_struct = 3L, d_feat = 3L) {
  X <- matrix(rnorm(n * d_struct), n, d_struct)
  F <- matrix(rnorm(n * d_feat), n, d_feat)
  C <- as.matrix(dist(X))
  C <- C / max(C)
  list(C = C, F = F, w = rep(1 / n, n))
}

make_subjects <- function(S, n, jitter = 0.05) {
  out <- vector("list", S)
  for (s in seq_len(S)) {
    ns <- max(8L, as.integer(round(n * (1 + runif(1, -jitter, jitter)))))
    out[[s]] <- c(make_subject_set(ns), list(id = sprintf("s%02d", s)))
  }
  out
}

rank_for_n <- function(n) {
  as.integer(max(32L, min(128L, round(sqrt(n) * 3))))
}

time_many <- function(fn, iters = 1L) {
  vals <- numeric(iters)
  objs <- numeric(iters)
  for (i in seq_len(iters)) {
    gc()
    t0 <- proc.time()[["elapsed"]]
    out <- fn()
    t1 <- proc.time()[["elapsed"]]
    vals[[i]] <- (t1 - t0) * 1000
    objs[[i]] <- out$objective_total
  }
  list(
    min_ms = min(vals),
    median_ms = stats::median(vals),
    objective = stats::median(objs)
  )
}

add_row <- function(rows, suite, S, n, method, thread_count, res) {
  rows[[length(rows) + 1L]] <- data.frame(
    suite = suite,
    S = S,
    n = n,
    method = method,
    thread_count = thread_count,
    min_ms = res$min_ms,
    median_ms = res$median_ms,
    iter_per_sec = if (res$median_ms > 0) 1000 / res$median_ms else NA_real_,
    objective = res$objective,
    runs = iters,
    stringsAsFactors = FALSE
  )
  rows
}

S <- 10L
rows <- list()

for (n in n_grid) {
  message("Preparing data for S=", S, ", n=", n)
  subjects <- make_subjects(S, n)
  tpl <- rfugw::multialign_make_template(
    subjects = subjects,
    k = n,
    feature_normalization = "zscore",
    seed = seed + n
  )

  if (n <= 500L) {
    res <- time_many(function() rfugw::multialign_fit(
      subjects = subjects,
      template_mode = "fixed",
      template = tpl,
      method = "fgw_entropic",
      alpha = 0.5,
      epsilon = 0.05,
      precision = "double",
      feature_normalization = "zscore",
      max_iter = 80L,
      sinkhorn_max_iter = 300L,
      use_cpp_batch = FALSE
    ), iters = iters)
    rows <- add_row(rows, "multiset_large", S, n, "multiset_fixed_rloop", NA_integer_, res)
  }

  for (th in thread_grid) {
    Sys.setenv(OMP_NUM_THREADS = as.character(th))
    lr_rank <- rank_for_n(n)
    res_fixed <- time_many(function() rfugw::multialign_fit(
      subjects = subjects,
      template_mode = "fixed",
      template = tpl,
      method = "fgw_entropic",
      alpha = 0.5,
      epsilon = 0.05,
      precision = "double",
      feature_normalization = "zscore",
      max_iter = 80L,
      sinkhorn_max_iter = 300L,
      use_cpp_batch = TRUE,
      n_threads = th
    ), iters = iters)
    rows <- add_row(rows, "multiset_large", S, n, "multiset_fixed_cpp_batch", th, res_fixed)

    res_fixed_mixed <- time_many(function() rfugw::multialign_fit(
      subjects = subjects,
      template_mode = "fixed",
      template = tpl,
      method = "fgw_entropic",
      alpha = 0.5,
      epsilon = 0.05,
      precision = "mixed",
      feature_normalization = "zscore",
      max_iter = 80L,
      sinkhorn_max_iter = 300L,
      use_cpp_batch = TRUE,
      n_threads = th
    ), iters = iters)
    rows <- add_row(rows, "multiset_large", S, n, "multiset_fixed_cpp_batch_mixed", th, res_fixed_mixed)

    res_fixed_mixed_auto <- time_many(function() rfugw::multialign_fit(
      subjects = subjects,
      template_mode = "fixed",
      template = tpl,
      method = "fgw_entropic",
      alpha = 0.5,
      epsilon = 0.05,
      precision = "mixed",
      feature_normalization = "zscore",
      max_iter = 80L,
      sinkhorn_max_iter = 300L,
      autotune = TRUE,
      autotune_level = "aggressive",
      use_cpp_batch = TRUE,
      n_threads = th
    ), iters = iters)
    rows <- add_row(rows, "multiset_large", S, n, "multiset_fixed_cpp_batch_mixed_auto", th, res_fixed_mixed_auto)

    if (n >= 500L) {
      res_fixed_mixed_auto_lr <- time_many(function() rfugw::multialign_fit(
        subjects = subjects,
        template_mode = "fixed",
        template = tpl,
        method = "fgw_entropic",
        alpha = 0.5,
        epsilon = 0.05,
        precision = "mixed",
        feature_normalization = "zscore",
        max_iter = 80L,
        sinkhorn_max_iter = 300L,
        autotune = TRUE,
        autotune_level = "aggressive",
        structure_rank = lr_rank,
        use_cpp_batch = TRUE,
        n_threads = th
      ), iters = iters)
      rows <- add_row(rows, "multiset_large", S, n, sprintf("multiset_fixed_cpp_batch_mixed_auto_lr%d", lr_rank), th, res_fixed_mixed_auto_lr)
    }

    if (n <= 1000L) {
      res_learned <- time_many(function() rfugw::multialign_fit(
        subjects = subjects,
        template_mode = "learned",
        template = tpl,
        method = "fgw_entropic",
        alpha = 0.5,
        epsilon = 0.05,
        precision = "double",
        feature_normalization = "zscore",
        max_iter = 60L,
        sinkhorn_max_iter = 250L,
        use_cpp_batch = TRUE,
        n_threads = th,
        template_max_iter = 2L,
        template_tol = 1e-6,
        final_refit = FALSE
      ), iters = iters)
      rows <- add_row(rows, "multiset_large", S, n, "multiset_learned_cpp_batch", th, res_learned)

      res_learned_mixed <- time_many(function() rfugw::multialign_fit(
        subjects = subjects,
        template_mode = "learned",
        template = tpl,
        method = "fgw_entropic",
        alpha = 0.5,
        epsilon = 0.05,
        precision = "mixed",
        feature_normalization = "zscore",
        max_iter = 60L,
        sinkhorn_max_iter = 250L,
        use_cpp_batch = TRUE,
        n_threads = th,
        template_max_iter = 2L,
        template_tol = 1e-6,
        final_refit = FALSE
      ), iters = iters)
      rows <- add_row(rows, "multiset_large", S, n, "multiset_learned_cpp_batch_mixed", th, res_learned_mixed)

      res_learned_mixed_auto <- time_many(function() rfugw::multialign_fit(
        subjects = subjects,
        template_mode = "learned",
        template = tpl,
        method = "fgw_entropic",
        alpha = 0.5,
        epsilon = 0.05,
        precision = "mixed",
        feature_normalization = "zscore",
        max_iter = 60L,
        sinkhorn_max_iter = 250L,
        autotune = TRUE,
        autotune_level = "aggressive",
        use_cpp_batch = TRUE,
        n_threads = th,
        template_max_iter = 2L,
        template_tol = 1e-6,
        final_refit = FALSE
      ), iters = iters)
      rows <- add_row(rows, "multiset_large", S, n, "multiset_learned_cpp_batch_mixed_auto", th, res_learned_mixed_auto)

      if (n >= 500L) {
        res_learned_mixed_auto_lr <- time_many(function() rfugw::multialign_fit(
          subjects = subjects,
          template_mode = "learned",
          template = tpl,
          method = "fgw_entropic",
          alpha = 0.5,
          epsilon = 0.05,
          precision = "mixed",
          feature_normalization = "zscore",
          max_iter = 60L,
          sinkhorn_max_iter = 250L,
          autotune = TRUE,
          autotune_level = "aggressive",
          structure_rank = lr_rank,
          use_cpp_batch = TRUE,
          n_threads = th,
          template_max_iter = 2L,
          template_tol = 1e-6,
          final_refit = FALSE
        ), iters = iters)
        rows <- add_row(rows, "multiset_large", S, n, sprintf("multiset_learned_cpp_batch_mixed_auto_lr%d", lr_rank), th, res_learned_mixed_auto_lr)
      }
    }
  }
}

out <- do.call(rbind, rows)
dir.create(dirname(out_csv), recursive = TRUE, showWarnings = FALSE)
write.csv(out, out_csv, row.names = FALSE)
cat("Wrote large multiset benchmark CSV:", out_csv, "\n")
print(out)
