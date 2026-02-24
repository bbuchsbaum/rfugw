suppressPackageStartupMessages({
  rlib <- Sys.getenv("RFUGW_RLIB", unset = "")
  if (nzchar(rlib)) {
    .libPaths(c(rlib, .libPaths()))
  }
  library(rfugw)
  library(bench)
})

args <- commandArgs(trailingOnly = TRUE)
iters <- if (length(args) >= 1) as.integer(args[[1]]) else 3L
out_csv <- if (length(args) >= 2) args[[2]] else "inst/bench/results/benchmark_ports_latest.csv"
seed <- if (length(args) >= 3) as.integer(args[[3]]) else 123L
set.seed(seed)

bench_one <- function(label, fn, suite, n, warm_obj = NA_real_, warm_iters = NA_integer_, runs = iters) {
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
    warm_obj = warm_obj,
    warm_iters = warm_iters,
    runs = runs,
    stringsAsFactors = FALSE
  )
}

make_pair_cost <- function(n, d = 3L) {
  X1 <- matrix(rnorm(n * d), n, d)
  X2 <- matrix(rnorm(n * d), n, d)
  C1 <- as.matrix(dist(X1)); C1 <- C1 / max(C1)
  C2 <- as.matrix(dist(X2)); C2 <- C2 / max(C2)
  list(C1 = C1, C2 = C2)
}

make_feature_cost <- function(ns, nt, d = 2L) {
  F1 <- matrix(rnorm(ns * d), ns, d)
  F2 <- matrix(rnorm(nt * d), nt, d)
  M <- as.matrix(dist(rbind(F1, F2)))[seq_len(ns), ns + seq_len(nt)]
  M <- M / max(M)
  M
}

make_ucoot_problem <- function(ns, nf, nt, mf) {
  X <- matrix(rnorm(ns * nf), ns, nf)
  Y <- matrix(rnorm(nt * mf), nt, mf)
  list(X = X, Y = Y)
}

rows <- list()
k <- 1L

# Partial exact/entropic
for (n in c(25L, 35L)) {
  ns <- n
  nt <- n + 4L
  costs <- make_pair_cost(max(ns, nt))
  C1 <- costs$C1[seq_len(ns), seq_len(ns), drop = FALSE]
  C2 <- costs$C2[seq_len(nt), seq_len(nt), drop = FALSE]
  M <- make_feature_cost(ns, nt)
  p <- rep(1 / ns, ns)
  q <- rep(1 / nt, nt)

  ctrl_epgw <- list(C1 = C1, C2 = C2, p = p, q = q, reg = 0.2, m = 0.7, numItermax = 80L, tol = 1e-7)
  ctrl_epfgw <- c(ctrl_epgw, list(M = M, alpha = 0.6))

  warm_epgw <- do.call(rfugw::entropic_partial_gromov_wasserstein, c(ctrl_epgw, list(log = TRUE)))
  warm_epfgw <- do.call(rfugw::entropic_partial_fused_gromov_wasserstein, c(ctrl_epfgw, list(log = TRUE)))

  rows[[k]] <- bench_one(
    "rfugw_entropic_partial_gw",
    function() do.call(rfugw::entropic_partial_gromov_wasserstein, c(ctrl_epgw, list(log = FALSE))),
    "partial",
    n,
    warm_epgw$partial_gw_dist,
    warm_epgw$iterations
  ); k <- k + 1L

  rows[[k]] <- bench_one(
    "rfugw_entropic_partial_fgw",
    function() do.call(rfugw::entropic_partial_fused_gromov_wasserstein, c(ctrl_epfgw, list(log = FALSE))),
    "partial",
    n,
    warm_epfgw$partial_fgw_dist,
    warm_epfgw$iterations
  ); k <- k + 1L

  if (requireNamespace("lpSolve", quietly = TRUE)) {
    ctrl_pgw <- list(C1 = C1, C2 = C2, p = p, q = q, m = 0.7, nb_dummies = 5L, numItermax = 120L, tol = 1e-8, lp_solver = "lp_matrix")
    ctrl_pfgw <- c(ctrl_pgw, list(M = M, alpha = 0.5))

    warm_pgw <- do.call(rfugw::partial_gromov_wasserstein, c(ctrl_pgw, list(log = TRUE)))
    warm_pfgw <- do.call(rfugw::partial_fused_gromov_wasserstein, c(ctrl_pfgw, list(log = TRUE)))

    rows[[k]] <- bench_one(
      "rfugw_partial_gw_exact",
      function() do.call(rfugw::partial_gromov_wasserstein, c(ctrl_pgw, list(log = FALSE))),
      "partial",
      n,
      warm_pgw$partial_gw_dist,
      warm_pgw$iterations
    ); k <- k + 1L

    rows[[k]] <- bench_one(
      "rfugw_partial_fgw_exact",
      function() do.call(rfugw::partial_fused_gromov_wasserstein, c(ctrl_pfgw, list(log = FALSE))),
      "partial",
      n,
      warm_pfgw$partial_fgw_dist,
      warm_pfgw$iterations
    ); k <- k + 1L
  }
}

# Non-entropic semirelaxed
for (n in c(80L, 140L, 220L)) {
  ns <- n
  nt <- as.integer(round(n * 0.8))
  costs <- make_pair_cost(max(ns, nt))
  C1 <- costs$C1[seq_len(ns), seq_len(ns), drop = FALSE]
  C2 <- costs$C2[seq_len(nt), seq_len(nt), drop = FALSE]
  M <- make_feature_cost(ns, nt)
  p <- rep(1 / ns, ns)

  ctrl_srgw <- list(C1 = C1, C2 = C2, p = p, max_iter = 200L, tol_rel = 1e-9, tol_abs = 1e-9)
  ctrl_srfgw <- c(ctrl_srgw, list(M = M, alpha = 0.55))

  warm_srgw <- do.call(rfugw::semirelaxed_gromov_wasserstein, ctrl_srgw)
  warm_srfgw <- do.call(rfugw::semirelaxed_fused_gromov_wasserstein, ctrl_srfgw)

  rows[[k]] <- bench_one(
    "rfugw_semirelaxed_gw",
    function() do.call(rfugw::semirelaxed_gromov_wasserstein, ctrl_srgw),
    "semirelaxed",
    n,
    warm_srgw$srgw_dist,
    warm_srgw$iterations
  ); k <- k + 1L

  rows[[k]] <- bench_one(
    "rfugw_semirelaxed_fgw",
    function() do.call(rfugw::semirelaxed_fused_gromov_wasserstein, ctrl_srfgw),
    "semirelaxed",
    n,
    warm_srfgw$srfgw_dist,
    warm_srfgw$iterations
  ); k <- k + 1L
}

# Sampled GW
for (n in c(90L, 160L, 260L)) {
  costs <- make_pair_cost(n)
  p <- rep(1 / n, n)
  q <- rep(1 / n, n)

  ctrl <- list(
    C1 = costs$C1,
    C2 = costs$C2,
    p = p,
    q = q,
    nb_samples_grad = c(16L, 2L),
    epsilon = 0.1,
    max_iter = 120L,
    random_state = seed
  )

  warm <- do.call(rfugw::sampled_gromov_wasserstein, c(ctrl, list(log = TRUE)))

  rows[[k]] <- bench_one(
    "rfugw_sampled_gw",
    function() do.call(rfugw::sampled_gromov_wasserstein, c(ctrl, list(log = FALSE))),
    "sampled",
    n,
    warm$gw_dist_estimated,
    warm$iterations
  ); k <- k + 1L
}

# UCOOT/across-space
for (n in c(40L, 70L, 110L)) {
  prob <- make_ucoot_problem(ns = n, nf = max(10L, as.integer(round(n * 0.6))), nt = n - 5L, mf = max(9L, as.integer(round(n * 0.5))))

  ctrl <- list(
    X = prob$X,
    Y = prob$Y,
    reg_marginals = c(10, 8),
    epsilon = c(0.05, 0.03),
    divergence = "kl",
    unbalanced_solver = "sinkhorn",
    max_iter = 40L,
    tol = 1e-7,
    max_iter_ot = 200L,
    tol_ot = 1e-7
  )

  warm <- do.call(rfugw::unbalanced_co_optimal_transport, c(ctrl, list(log = TRUE)))

  rows[[k]] <- bench_one(
    "rfugw_ucoot_kl",
    function() do.call(rfugw::unbalanced_co_optimal_transport, c(ctrl, list(log = FALSE))),
    "ucoot",
    n,
    warm$ucoot_cost,
    length(warm$error)
  ); k <- k + 1L
}

out <- do.call(rbind, rows)
dir.create(dirname(out_csv), recursive = TRUE, showWarnings = FALSE)
write.csv(out, out_csv, row.names = FALSE)

cat("== rfugw benchmark ports ==\n")
print(out[, c("suite", "n", "method", "median_ms", "iter_per_sec", "warm_obj", "warm_iters")])
cat(sprintf("\nWrote benchmark CSV: %s\n", out_csv))
