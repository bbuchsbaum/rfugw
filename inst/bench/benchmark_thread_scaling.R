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
default_out_csv <- if (file.exists("inst/bench/results")) {
  "inst/bench/results/thread_scaling_latest.csv"
} else {
  "rfugw/inst/bench/results/thread_scaling_latest.csv"
}
out_csv <- if (length(args) >= 1) args[[1]] else default_out_csv
iters <- if (length(args) >= 2) as.integer(args[[2]]) else 3L
seed <- if (length(args) >= 3) as.integer(args[[3]]) else 42L
thread_count <- as.integer(Sys.getenv("OMP_NUM_THREADS", unset = "1"))

set.seed(seed)

make_fgw_problem <- function(n, d_struct = 3L, d_feat = 5L) {
  X1 <- matrix(rnorm(n * d_struct), n, d_struct)
  X2 <- matrix(rnorm(n * d_struct), n, d_struct)
  Y1 <- matrix(rnorm(n * d_feat), n, d_feat)
  Y2 <- matrix(rnorm(n * d_feat), n, d_feat)
  C1 <- as.matrix(dist(X1)); C1 <- C1 / max(C1)
  C2 <- as.matrix(dist(X2)); C2 <- C2 / max(C2)
  M <- as.matrix(dist(rbind(Y1, Y2)))[1:n, (n + 1):(2 * n)]
  M <- M / max(M)
  list(C1 = C1, C2 = C2, M = M)
}

make_fugw_problem <- function(n) {
  X1 <- matrix(rnorm(n * 3), n, 3)
  X2 <- matrix(rnorm(n * 3), n, 3)
  Cx <- as.matrix(dist(X1)); Cx <- Cx / max(Cx)
  Cy <- as.matrix(dist(X2)); Cy <- Cy / max(Cy)
  M <- matrix(1, nrow = n, ncol = n)
  idx <- seq_len(n)
  M[cbind(idx, rev(idx))] <- 0
  list(Cx = Cx, Cy = Cy, M = M)
}

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

fgw <- make_fgw_problem(160L)
fugw <- make_fugw_problem(90L)
subjects <- make_subjects(S = 10L, n = 200L)
lr_rank <- rank_for_n(200L)
tpl <- rfugw::multialign_make_template(
  subjects = subjects,
  k = 200L,
  feature_normalization = "zscore",
  seed = seed + 200L
)

b <- bench::mark(
  fgw_pgd = rfugw::fgw_entropic(
    fgw$M, fgw$C1, fgw$C2,
    alpha = 0.5, epsilon = 0.05, solver = "PGD",
    max_iter = 200, sinkhorn_max_iter = 500, sinkhorn_tol = 1e-9
  ),
  fgw_pgd_log = rfugw::fgw_entropic(
    fgw$M, fgw$C1, fgw$C2,
    alpha = 0.5, epsilon = 0.05, solver = "PGD", sinkhorn_method = "log",
    max_iter = 200, sinkhorn_max_iter = 500, sinkhorn_tol = 1e-9
  ),
  fgw_pgd_mixed = rfugw::fgw_entropic(
    fgw$M, fgw$C1, fgw$C2,
    alpha = 0.5, epsilon = 0.05, solver = "PGD", precision = "mixed",
    max_iter = 200, sinkhorn_max_iter = 500, sinkhorn_tol = 1e-9
  ),
  fugw_kl = rfugw::fugw_kl(
    fugw$Cx, fugw$Cy,
    reg_marginals = c(100, 50), epsilon = 1e-2, alpha = 0.5, M = fugw$M,
    max_iter = 80, max_iter_ot = 300, tol = 1e-8, tol_ot = 1e-8
  ),
  multiset_fixed_cpp_batch = rfugw::multialign_fit(
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
    n_threads = thread_count
  ),
  multiset_fixed_cpp_batch_mixed = rfugw::multialign_fit(
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
    n_threads = thread_count
  ),
  multiset_fixed_cpp_batch_mixed_auto = rfugw::multialign_fit(
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
    n_threads = thread_count
  ),
  multiset_fixed_cpp_batch_mixed_auto_lr = rfugw::multialign_fit(
    subjects = subjects,
    template_mode = "fixed",
    template = tpl,
    method = "fgw_entropic",
    alpha = 0.5,
    epsilon = 0.05,
    precision = "mixed",
    autotune = TRUE,
    autotune_level = "aggressive",
    structure_rank = lr_rank,
    feature_normalization = "zscore",
    max_iter = 80L,
    sinkhorn_max_iter = 300L,
    use_cpp_batch = TRUE,
    n_threads = thread_count
  ),
  iterations = iters,
  check = FALSE,
  memory = TRUE,
  filter_gc = TRUE
)

tb <- tibble::as_tibble(b)[, c("expression", "min", "median", "itr/sec", "mem_alloc")]
expr_desc <- attr(tb$expression, "description")
if (is.null(expr_desc)) {
  expr_desc <- vapply(tb$expression, function(e) paste(deparse(e), collapse = ""), character(1))
}

rows <- data.frame(
  thread_count = thread_count,
  suite = ifelse(
    grepl("^fgw_", expr_desc), "fgw",
    ifelse(grepl("^multiset_", expr_desc), "multiset", "fugw")
  ),
  method = expr_desc,
  min_ms = as.numeric(tb$min) * 1000,
  median_ms = as.numeric(tb$median) * 1000,
  iter_per_sec = as.numeric(tb$`itr/sec`),
  mem_bytes = as.numeric(tb$mem_alloc),
  runs = iters
)

dir.create(dirname(out_csv), recursive = TRUE, showWarnings = FALSE)
write.csv(rows, out_csv, row.names = FALSE)
cat("Wrote thread benchmark CSV:", out_csv, "\n")
print(rows)
