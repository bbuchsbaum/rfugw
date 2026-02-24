suppressPackageStartupMessages({
  rlib <- Sys.getenv("RFUGW_RLIB", unset = "")
  if (nzchar(rlib)) {
    .libPaths(c(rlib, .libPaths()))
  }
  library(rfugw)
  library(bench)
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
iters <- if (length(args) >= 1) as.integer(args[[1]]) else 2L
out_csv <- if (length(args) >= 2) args[[2]] else "inst/bench/results/benchmark_multiset_latest.csv"
seed <- if (length(args) >= 3) as.integer(args[[3]]) else 42L

set.seed(seed)

make_subject_set <- function(n, d_struct = 3L, d_feat = 3L) {
  X <- matrix(rnorm(n * d_struct), n, d_struct)
  F <- matrix(rnorm(n * d_feat), n, d_feat)
  C <- as.matrix(dist(X))
  C <- C / max(C)
  list(C = C, F = F, w = rep(1 / n, n))
}

make_subjects <- function(S, n, jitter = 0.2) {
  out <- vector("list", S)
  for (s in seq_len(S)) {
    ns <- max(6L, as.integer(round(n * (1 + runif(1, -jitter, jitter)))))
    out[[s]] <- c(make_subject_set(ns), list(id = sprintf("s%02d", s)))
  }
  out
}

bench_one <- function(label, fn) {
  b <- bench::mark(
    run = fn(),
    iterations = iters,
    check = FALSE,
    memory = TRUE,
    filter_gc = TRUE
  )
  tb <- tibble::as_tibble(b)[1, c("min", "median", "itr/sec", "mem_alloc")]
  data.frame(
    method = label,
    min_ms = as.numeric(tb$min[[1]]) * 1000,
    median_ms = as.numeric(tb$median[[1]]) * 1000,
    iter_per_sec = as.numeric(tb$`itr/sec`[[1]]),
    mem_bytes = as.numeric(tb$mem_alloc[[1]])
  )
}

grid <- expand.grid(
  S = c(3L, 5L, 10L),
  n = c(80L, 120L),
  stringsAsFactors = FALSE
)
grid <- subset(grid, !(S == 10L & n == 120L))

rows <- list()
k <- 1L
for (g in seq_len(nrow(grid))) {
  S <- grid$S[[g]]
  n <- grid$n[[g]]
  subjects <- make_subjects(S, n)
  tpl <- rfugw::multialign_make_template(
    subjects = subjects,
    k = n,
    feature_normalization = "zscore",
    seed = seed + g
  )

  warm <- rfugw::multialign_fit(
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
    n_threads = as.integer(Sys.getenv("OMP_NUM_THREADS", unset = "1"))
  )

  rows[[k]] <- cbind(
    data.frame(
      suite = "multiset",
      S = S,
      n = n,
      warm_obj = warm$objective_total,
      warm_outer = warm$outer_iterations
    ),
    bench_one("multiset_fixed_rloop", function() rfugw::multialign_fit(
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
    ))
  )
  k <- k + 1L

  rows[[k]] <- cbind(
    data.frame(
      suite = "multiset",
      S = S,
      n = n,
      warm_obj = warm$objective_total,
      warm_outer = warm$outer_iterations
    ),
    bench_one("multiset_fixed_cpp_batch_mixed_auto", function() rfugw::multialign_fit(
      subjects = subjects,
      template_mode = "fixed",
      template = tpl,
      method = "fgw_entropic",
      alpha = 0.5,
      epsilon = 0.05,
      precision = "mixed",
      autotune = TRUE,
      autotune_level = "aggressive",
      feature_normalization = "zscore",
      max_iter = 80L,
      sinkhorn_max_iter = 300L,
      use_cpp_batch = TRUE,
      n_threads = as.integer(Sys.getenv("OMP_NUM_THREADS", unset = "1"))
    ))
  )
  k <- k + 1L

  rows[[k]] <- cbind(
    data.frame(
      suite = "multiset",
      S = S,
      n = n,
      warm_obj = warm$objective_total,
      warm_outer = warm$outer_iterations
    ),
    bench_one("multiset_fixed_cpp_batch", function() rfugw::multialign_fit(
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
      n_threads = as.integer(Sys.getenv("OMP_NUM_THREADS", unset = "1"))
    ))
  )
  k <- k + 1L

  rows[[k]] <- cbind(
    data.frame(
      suite = "multiset",
      S = S,
      n = n,
      warm_obj = warm$objective_total,
      warm_outer = warm$outer_iterations
    ),
    bench_one("multiset_fixed_cpp_batch_mixed", function() rfugw::multialign_fit(
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
      n_threads = as.integer(Sys.getenv("OMP_NUM_THREADS", unset = "1"))
    ))
  )
  k <- k + 1L

  rows[[k]] <- cbind(
    data.frame(
      suite = "multiset",
      S = S,
      n = n,
      warm_obj = warm$objective_total,
      warm_outer = warm$outer_iterations
    ),
    bench_one("multiset_learned_cpp_batch", function() rfugw::multialign_fit(
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
      n_threads = as.integer(Sys.getenv("OMP_NUM_THREADS", unset = "1")),
      template_max_iter = 3L,
      template_tol = 1e-6
    ))
  )
  k <- k + 1L

  rows[[k]] <- cbind(
    data.frame(
      suite = "multiset",
      S = S,
      n = n,
      warm_obj = warm$objective_total,
      warm_outer = warm$outer_iterations
    ),
    bench_one("multiset_learned_cpp_batch_mixed", function() rfugw::multialign_fit(
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
      n_threads = as.integer(Sys.getenv("OMP_NUM_THREADS", unset = "1")),
      template_max_iter = 3L,
      template_tol = 1e-6
    ))
  )
  k <- k + 1L

  rows[[k]] <- cbind(
    data.frame(
      suite = "multiset",
      S = S,
      n = n,
      warm_obj = warm$objective_total,
      warm_outer = warm$outer_iterations
    ),
    bench_one("multiset_learned_cpp_batch_mixed_auto", function() rfugw::multialign_fit(
      subjects = subjects,
      template_mode = "learned",
      template = tpl,
      method = "fgw_entropic",
      alpha = 0.5,
      epsilon = 0.05,
      precision = "mixed",
      autotune = TRUE,
      autotune_level = "aggressive",
      feature_normalization = "zscore",
      max_iter = 60L,
      sinkhorn_max_iter = 250L,
      use_cpp_batch = TRUE,
      n_threads = as.integer(Sys.getenv("OMP_NUM_THREADS", unset = "1")),
      template_max_iter = 3L,
      template_tol = 1e-6
    ))
  )
  k <- k + 1L
}

out <- do.call(rbind, rows)
dir.create(dirname(out_csv), recursive = TRUE, showWarnings = FALSE)
write.csv(out, out_csv, row.names = FALSE)
cat("Wrote multiset benchmark CSV:", out_csv, "\n")
print(out)
