#!/usr/bin/env Rscript
# Quality-versus-budget, rank, and memory curves for sampled / approximate GW.

args <- commandArgs(trailingOnly = TRUE)
out_dir <- if (length(args) >= 1L) args[[1]] else "inst/bench/results"
seed <- if (length(args) >= 2L) as.integer(args[[2]]) else 20260816L
n_seeds <- if (length(args) >= 3L) as.integer(args[[3]]) else 5L

if (!file.exists("DESCRIPTION")) {
  stop("Run inst/bench/sampled_budget_curves.R from the package root.", call. = FALSE)
}

suppressPackageStartupMessages({
  rlib <- Sys.getenv("RFUGW_RLIB", unset = "")
  if (nzchar(rlib)) .libPaths(c(rlib, .libPaths()))
  library(rfugw)
})

source("inst/bench/protocol.R")
threads <- bench_pin_threads(1L)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

make_pair <- function(n, seed, d = 2L) {
  set.seed(as.integer(seed) + as.integer(n))
  X1 <- matrix(rnorm(n * d), n, d)
  X2 <- matrix(rnorm(n * d), n, d)
  C1 <- as.matrix(dist(X1)); C1 <- C1 / max(C1)
  C2 <- as.matrix(dist(X2)); C2 <- C2 / max(C2)
  list(X1 = X1, X2 = X2, C1 = C1, C2 = C2,
       p = rep(1 / n, n), q = rep(1 / n, n))
}

make_knn <- function(X, k = 5L) {
  D <- as.matrix(dist(X))
  W <- matrix(0, nrow(X), nrow(X))
  for (i in seq_len(nrow(X))) {
    keep <- order(D[i, ])[seq_len(k + 1L)]
    keep <- keep[keep != i]
    W[i, keep] <- 1
  }
  pmax(W, t(W))
}

lowrank_rel_error <- function(T, rank) {
  r <- as.integer(rank)
  sv <- svd(T, nu = r, nv = r)
  s <- pmax(sv$d[seq_len(r)], 0)
  S <- diag(sqrt(s), nrow = r, ncol = r)
  Q <- sv$u[, seq_len(r), drop = FALSE] %*% S
  R <- sv$v[, seq_len(r), drop = FALSE] %*% S
  Tr <- Q %*% t(R)
  sqrt(sum((T - Tr)^2)) / sqrt(sum(T^2))
}

elapsed_ms <- function(fn) {
  t0 <- proc.time()[["elapsed"]]
  fn()
  (proc.time()[["elapsed"]] - t0) * 1000
}

n <- 16L
d <- make_pair(n, seed)
dense <- entropic_gromov_wasserstein(
  d$C1, d$C2, d$p, d$q,
  epsilon = 0.1, max_iter = 80L, tol = 1e-7
)
gw_ref <- ot_gw_square(dense$plan, d$C1, d$C2)
budgets <- list(c(2L, 1L), c(4L, 1L), c(8L, 2L), c(n, n))
seeds <- as.integer(seed) + seq_len(n_seeds)

budget_rows <- list()
k <- 1L
for (b in budgets) {
  gaps <- numeric(n_seeds)
  gws <- numeric(n_seeds)
  frobs <- numeric(n_seeds)
  times <- numeric(n_seeds)
  for (i in seq_len(n_seeds)) {
    t0 <- proc.time()[["elapsed"]]
    out <- sampled_gromov_wasserstein(
      d$C1, d$C2, d$p, d$q,
      nb_samples_grad = b,
      epsilon = 0.1,
      max_iter = 40L,
      random_state = seeds[[i]],
      log = TRUE
    )
    times[[i]] <- (proc.time()[["elapsed"]] - t0) * 1000
    gws[[i]] <- ot_gw_square(out$plan, d$C1, d$C2)
    gaps[[i]] <- abs(gws[[i]] - gw_ref)
    frobs[[i]] <- sqrt(sum((out$plan - dense$plan)^2))
  }
  budget_rows[[k]] <- data.frame(
    n = n,
    nb_p = b[[1]],
    nb_q = b[[2]],
    median_gap = stats::median(gaps),
    mean_gap = mean(gaps),
    median_frob = stats::median(frobs),
    median_gw = stats::median(gws),
    median_ms = stats::median(times),
    n_seeds = n_seeds,
    stringsAsFactors = FALSE
  )
  k <- k + 1L
}
budget_df <- do.call(rbind, budget_rows)

rank_vals <- c(1L, 2L, 4L, 8L, n)
rank_rows <- lapply(rank_vals, function(r) {
  t_ms <- elapsed_ms(function() {
    invisible(lowrank_rel_error(dense$plan, r))
  })
  data.frame(
    n = n,
    rank = r,
    rel_error = lowrank_rel_error(dense$plan, r),
    svd_ms = t_ms,
    stringsAsFactors = FALSE
  )
})
rank_df <- do.call(rbind, rank_rows)

mem_sizes <- c(32L, 64L, 96L)
mem_rows <- lapply(mem_sizes, function(nn) {
  dd <- make_pair(nn, seed, d = 3L)
  W1 <- make_knn(dd$X1, k = 6L)
  W2 <- make_knn(dd$X2, k = 6L)
  if (requireNamespace("Matrix", quietly = TRUE)) {
    Ws1 <- Matrix::Matrix(W1, sparse = TRUE)
    Ws2 <- Matrix::Matrix(W2, sparse = TRUE)
    graph_bytes <- as.numeric(object.size(Ws1)) + as.numeric(object.size(Ws2))
  } else {
    graph_bytes <- as.numeric(object.size(W1)) + as.numeric(object.size(W2))
  }
  coords <- graph_diffusion_coordinates(W1, n_components = 8L, diffusion_time = 1)
  data.frame(
    n = nn,
    dense_cost_bytes = as.numeric(object.size(dd$C1)) + as.numeric(object.size(dd$C2)),
    coord_bytes = as.numeric(object.size(dd$X1)) + as.numeric(object.size(dd$X2)),
    graph_bytes = graph_bytes,
    embed_bytes = as.numeric(object.size(coords)),
    plan_bytes = as.numeric(object.size(matrix(0, nn, nn))),
    stringsAsFactors = FALSE
  )
})
mem_df <- do.call(rbind, mem_rows)

budget_csv <- file.path(out_dir, "sampled_budget_curves.csv")
rank_csv <- file.path(out_dir, "sampled_rank_curves.csv")
mem_csv <- file.path(out_dir, "sampled_memory_scaling.csv")
utils::write.csv(budget_df, budget_csv, row.names = FALSE)
utils::write.csv(rank_df, rank_csv, row.names = FALSE)
utils::write.csv(mem_df, mem_csv, row.names = FALSE)

cat("threads=", threads, "\n", sep = "")
cat("dense_gw=", gw_ref, "\n", sep = "")
print(budget_df)
print(rank_df)
print(mem_df)
cat("wrote ", budget_csv, "\n", sep = "")
