suppressPackageStartupMessages({
  rlib <- Sys.getenv("RFUGW_RLIB", unset = "")
  if (nzchar(rlib)) {
    .libPaths(c(rlib, .libPaths()))
  }
  library(rfugw)
  library(bench)
  library(jsonlite)
})

if (!"sampled_gw_from_graphs" %in% getNamespaceExports("rfugw") &&
    requireNamespace("devtools", quietly = TRUE)) {
  devtools::load_all(".", quiet = TRUE)
}

args <- commandArgs(trailingOnly = TRUE)
iters <- if (length(args) >= 1) as.integer(args[[1]]) else 3L
out_csv <- if (length(args) >= 2) args[[2]] else "inst/bench/results/benchmark_sparse_sampled_latest.csv"
seed <- if (length(args) >= 3) as.integer(args[[3]]) else 123L
n_values <- if (length(args) >= 4) {
  as.integer(strsplit(args[[4]], ",", fixed = TRUE)[[1]])
} else {
  c(120L, 220L, 400L, 800L)
}
rfugw_sinkhorn_max_iter <- if (length(args) >= 5) as.integer(args[[5]]) else 80L
rfugw_sinkhorn_tol <- if (length(args) >= 6) as.numeric(args[[6]]) else 1e-6
n_values <- unique(n_values[is.finite(n_values) & n_values >= 20L])
if (length(n_values) == 0L) {
  stop("`n_values` must contain at least one integer >= 20.", call. = FALSE)
}
if (!is.finite(rfugw_sinkhorn_max_iter) || rfugw_sinkhorn_max_iter < 1L) {
  stop("`rfugw_sinkhorn_max_iter` must be an integer >= 1.", call. = FALSE)
}
if (!is.finite(rfugw_sinkhorn_tol) || rfugw_sinkhorn_tol <= 0) {
  stop("`rfugw_sinkhorn_tol` must be positive.", call. = FALSE)
}
set.seed(seed)

deterministic_points <- function(n, d = 3L, offset = 0L) {
  i <- seq_len(n) + as.integer(offset)
  out <- matrix(0, nrow = n, ncol = d)
  out[, 1] <- sin(i * 0.37) + cos(i * 0.11)
  if (d >= 2L) out[, 2] <- cos(i * 0.23) + sin(i * 0.17)
  if (d >= 3L) out[, 3] <- sin(i * 0.07) * cos(i * 0.13)
  if (d > 3L) {
    for (j in seq.int(4L, d)) {
      out[, j] <- sin(i * (0.03 * j + 0.01)) + cos(i * (0.02 * j + 0.05))
    }
  }
  out
}

make_knn_similarity <- function(X, k = 10L) {
  D <- as.matrix(dist(X))
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

bench_one <- function(label, fn, suite, n, gw_obj = NA_real_, obj_gap = NA_real_, runs = iters) {
  b <- bench::mark(
    run = fn(),
    iterations = runs,
    check = FALSE,
    memory = TRUE,
    filter_gc = TRUE
  )
  bt <- tibble::as_tibble(b)
  data.frame(
    suite = suite,
    n = n,
    method = label,
    min_ms = as.numeric(bt$min[[1]]) * 1000,
    median_ms = as.numeric(bt$median[[1]]) * 1000,
    iter_per_sec = as.numeric(bt$`itr/sec`[[1]]),
    mem_bytes = as.numeric(bt$mem_alloc[[1]]),
    gw_obj = gw_obj,
    objective_gap_dense = obj_gap,
    runs = runs,
    stringsAsFactors = FALSE
  )
}

python_bin <- Sys.getenv("RFUGW_PYTHON", unset = "python3")
pot_ready <- TRUE
probe <- suppressWarnings(system2(
  command = python_bin,
  args = c("-c", shQuote("import numpy,ot")),
  stdout = TRUE,
  stderr = TRUE
))
probe_status <- attr(probe, "status")
if (!is.null(probe_status) && probe_status != 0L) {
  pot_ready <- FALSE
  warning(
    "Skipping POT sparse sampled benchmark because Python dependencies are unavailable: ",
    paste(probe, collapse = "\n"),
    call. = FALSE
  )
}

run_pot_sampled <- function(mode, A, B, iters, seed, python_bin) {
  script <- file.path("inst", "bench", "benchmark_sparse_sampled_pot.py")
  if (!file.exists(script)) {
    stop("POT benchmark helper not found: ", script, call. = FALSE)
  }
  a_file <- tempfile(pattern = "rfugw_a_", fileext = ".csv")
  b_file <- tempfile(pattern = "rfugw_b_", fileext = ".csv")
  on.exit(unlink(c(a_file, b_file), force = TRUE), add = TRUE)
  write.table(A, file = a_file, sep = ",", row.names = FALSE, col.names = FALSE)
  write.table(B, file = b_file, sep = ",", row.names = FALSE, col.names = FALSE)

  out <- system2(
    command = python_bin,
    args = c(script, mode, a_file, b_file, as.integer(iters), as.integer(seed)),
    stdout = TRUE,
    stderr = TRUE
  )
  status <- attr(out, "status")
  if (!is.null(status) && status != 0L) {
    stop("POT sparse sampled benchmark failed: ", paste(out, collapse = "\n"), call. = FALSE)
  }
  fromJSON(out[[length(out)]], simplifyVector = TRUE)
}

rows <- list()
k <- 1L

if (!requireNamespace("Matrix", quietly = TRUE)) {
  stop("Package `Matrix` is required for sparse benchmark.", call. = FALSE)
}

for (n in n_values) {
  X1 <- deterministic_points(n, d = 3L, offset = 11L)
  X2 <- deterministic_points(n, d = 3L, offset = 37L)
  W1 <- Matrix::Matrix(make_knn_similarity(X1, k = 10L), sparse = TRUE)
  W2 <- Matrix::Matrix(make_knn_similarity(X2, k = 10L), sparse = TRUE)

  E1 <- rfugw::graph_diffusion_coordinates(W1, n_components = 12L, diffusion_time = 1, self_loop = 1e-6)
  E2 <- rfugw::graph_diffusion_coordinates(W2, n_components = 12L, diffusion_time = 1, self_loop = 1e-6)
  p <- rep(1 / nrow(E1), nrow(E1))
  q <- rep(1 / nrow(E2), nrow(E2))
  C1 <- as.matrix(dist(E1))
  C2 <- as.matrix(dist(E2))

  # Suite 1: solver-only (precomputed costs/embeddings)
  ctrl_dense_solver <- list(
    p = p,
    q = q,
    nb_samples_grad = c(16L, 2L),
    epsilon = 0.1,
    max_iter = 120L,
    sinkhorn_max_iter = rfugw_sinkhorn_max_iter,
    sinkhorn_tol = rfugw_sinkhorn_tol,
    random_state = seed
  )
  ctrl_coords_solver <- list(
    X1 = E1,
    X2 = E2,
    p = p,
    q = q,
    nb_samples_grad = c(16L, 2L),
    epsilon = 0.1,
    max_iter = 120L,
    sinkhorn_max_iter = rfugw_sinkhorn_max_iter,
    sinkhorn_tol = rfugw_sinkhorn_tol,
    random_state = seed,
    use_cpp = TRUE
  )
  warm_dense_solver <- do.call(
    rfugw::sampled_gromov_wasserstein,
    c(ctrl_dense_solver, list(C1 = C1, C2 = C2, log = TRUE))
  )
  warm_coords_solver <- do.call(rfugw::sampled_gromov_wasserstein_coords, c(ctrl_coords_solver, list(log = TRUE)))
  dense_ref_solver <- warm_dense_solver$gw_dist_estimated

  rows[[k]] <- bench_one(
    "rfugw_sampled_coords_solver",
    function() do.call(rfugw::sampled_gromov_wasserstein_coords, c(ctrl_coords_solver, list(log = FALSE))),
    "sampled_sparse_solver_only",
    n,
    gw_obj = warm_coords_solver$gw_dist_estimated,
    obj_gap = abs(warm_coords_solver$gw_dist_estimated - dense_ref_solver) / max(1, abs(dense_ref_solver))
  ); k <- k + 1L

  rows[[k]] <- bench_one(
    "rfugw_sampled_dense_solver",
    function() do.call(
      rfugw::sampled_gromov_wasserstein,
      c(ctrl_dense_solver, list(C1 = C1, C2 = C2, log = FALSE))
    ),
    "sampled_sparse_solver_only",
    n,
    gw_obj = warm_dense_solver$gw_dist_estimated,
    obj_gap = 0
  ); k <- k + 1L

  if (isTRUE(pot_ready)) {
    pot_solver <- tryCatch(
      run_pot_sampled("cost", C1, C2, iters, seed, python_bin),
      error = function(e) {
        warning(conditionMessage(e), call. = FALSE)
        NULL
      }
    )
    if (!is.null(pot_solver)) {
      rows[[k]] <- data.frame(
        suite = "sampled_sparse_solver_only",
        n = n,
        method = as.character(pot_solver$method),
        min_ms = as.numeric(pot_solver$min_ms),
        median_ms = as.numeric(pot_solver$median_ms),
        iter_per_sec = 1000 / as.numeric(pot_solver$median_ms),
        mem_bytes = as.numeric(pot_solver$mem_bytes),
        gw_obj = as.numeric(pot_solver$gw_obj),
        objective_gap_dense = abs(as.numeric(pot_solver$gw_obj) - dense_ref_solver) / max(1, abs(dense_ref_solver)),
        runs = as.integer(pot_solver$runs),
        stringsAsFactors = FALSE
      )
      k <- k + 1L
    }
  }

  # Suite 2: end-to-end from embeddings (dense includes dist build in timed block)
  run_coords_e2e <- function(log = FALSE) {
    do.call(rfugw::sampled_gromov_wasserstein_coords, c(ctrl_coords_solver, list(log = log)))
  }
  run_dense_e2e <- function(log = FALSE) {
    C1e <- as.matrix(dist(E1))
    C2e <- as.matrix(dist(E2))
    do.call(
      rfugw::sampled_gromov_wasserstein,
      c(ctrl_dense_solver, list(C1 = C1e, C2 = C2e, log = log))
    )
  }
  ctrl_pipe <- list(
    W1 = W1,
    W2 = W2,
    n_components = 12L,
    diffusion_time = 1,
    self_loop = 1e-6,
    nb_samples_grad = c(16L, 2L),
    epsilon = 0.1,
    max_iter = 120L,
    sinkhorn_max_iter = rfugw_sinkhorn_max_iter,
    sinkhorn_tol = rfugw_sinkhorn_tol,
    random_state = seed,
    use_cpp = TRUE
  )

  warm_dense_e2e <- run_dense_e2e(log = TRUE)
  warm_coords_e2e <- run_coords_e2e(log = TRUE)
  warm_pipe <- do.call(rfugw::sampled_gw_from_graphs, c(ctrl_pipe, list(log = TRUE, return_embeddings = FALSE)))
  dense_ref_e2e <- warm_dense_e2e$gw_dist_estimated

  rows[[k]] <- bench_one(
    "rfugw_sampled_coords_end_to_end",
    function() run_coords_e2e(log = FALSE),
    "sampled_sparse_end_to_end",
    n,
    gw_obj = warm_coords_e2e$gw_dist_estimated,
    obj_gap = abs(warm_coords_e2e$gw_dist_estimated - dense_ref_e2e) / max(1, abs(dense_ref_e2e))
  ); k <- k + 1L

  rows[[k]] <- bench_one(
    "rfugw_sampled_dense_end_to_end",
    function() run_dense_e2e(log = FALSE),
    "sampled_sparse_end_to_end",
    n,
    gw_obj = warm_dense_e2e$gw_dist_estimated,
    obj_gap = 0
  ); k <- k + 1L

  rows[[k]] <- bench_one(
    "rfugw_sparse_graph_pipeline",
    function() do.call(rfugw::sampled_gw_from_graphs, c(ctrl_pipe, list(log = FALSE, return_embeddings = FALSE))),
    "sampled_sparse_end_to_end",
    n,
    gw_obj = warm_pipe$gw_dist_estimated,
    obj_gap = abs(warm_pipe$gw_dist_estimated - dense_ref_e2e) / max(1, abs(dense_ref_e2e))
  ); k <- k + 1L

  if (isTRUE(pot_ready)) {
    pot_e2e <- tryCatch(
      run_pot_sampled("coords", E1, E2, iters, seed, python_bin),
      error = function(e) {
        warning(conditionMessage(e), call. = FALSE)
        NULL
      }
    )
    if (!is.null(pot_e2e)) {
      rows[[k]] <- data.frame(
        suite = "sampled_sparse_end_to_end",
        n = n,
        method = as.character(pot_e2e$method),
        min_ms = as.numeric(pot_e2e$min_ms),
        median_ms = as.numeric(pot_e2e$median_ms),
        iter_per_sec = 1000 / as.numeric(pot_e2e$median_ms),
        mem_bytes = as.numeric(pot_e2e$mem_bytes),
        gw_obj = as.numeric(pot_e2e$gw_obj),
        objective_gap_dense = abs(as.numeric(pot_e2e$gw_obj) - dense_ref_e2e) / max(1, abs(dense_ref_e2e)),
        runs = as.integer(pot_e2e$runs),
        stringsAsFactors = FALSE
      )
      k <- k + 1L
    }
  }
}

out <- do.call(rbind, rows)
dir.create(dirname(out_csv), recursive = TRUE, showWarnings = FALSE)
write.csv(out, out_csv, row.names = FALSE)

cat("== sparse sampled benchmark ==\n")
print(out[, c("suite", "n", "method", "median_ms", "mem_bytes", "gw_obj", "objective_gap_dense")])
cat(sprintf("\nWrote benchmark CSV: %s\n", out_csv))
