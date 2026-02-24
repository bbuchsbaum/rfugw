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
out_csv <- if (length(args) >= 2) args[[2]] else "inst/bench/results/benchmark_semirelaxed_latest.csv"
seed <- if (length(args) >= 3) as.integer(args[[3]]) else 42L
set.seed(seed)

make_problem <- function(n, d_struct = 3L, d_feat = 5L) {
  X1 <- matrix(rnorm(n * d_struct), n, d_struct)
  X2 <- matrix(rnorm(n * d_struct), n, d_struct)
  Y1 <- matrix(rnorm(n * d_feat), n, d_feat)
  Y2 <- matrix(rnorm(n * d_feat), n, d_feat)

  C1 <- as.matrix(dist(X1))
  C2 <- as.matrix(dist(X2))
  C1 <- C1 / max(C1)
  C2 <- C2 / max(C2)

  M <- as.matrix(dist(rbind(Y1, Y2)))[1:n, (n + 1):(2 * n)]
  M <- M / max(M)
  p <- rep(1 / n, n)

  list(M = M, C1 = C1, C2 = C2, p = p)
}

bench_one <- function(label, fn, suite, n, warm_obj, warm_iters, warm_error, runs) {
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
    warm_error = warm_error,
    runs = runs,
    stringsAsFactors = FALSE
  )
}

sizes <- c(80L, 120L, 180L, 240L)
rows <- list()
k <- 1L

for (n in sizes) {
  p <- make_problem(n)
  M0 <- matrix(0, nrow = n, ncol = n)
  ctrl_srfgw <- list(
    M = p$M,
    C1 = p$C1,
    C2 = p$C2,
    p = p$p,
    epsilon = 0.05,
    alpha = 0.5,
    symmetric = TRUE,
    max_iter = 200L,
    tol = 1e-9,
    check_every = 10L
  )
  ctrl_srgw <- list(
    C1 = p$C1,
    C2 = p$C2,
    p = p$p,
    epsilon = 0.05,
    symmetric = TRUE,
    max_iter = 200L,
    tol = 1e-9,
    check_every = 10L
  )

  warm_srfgw <- do.call(
    rfugw::entropic_semirelaxed_fused_gromov_wasserstein,
    c(ctrl_srfgw, list(backend = "cpp", precision = "mixed"))
  )
  warm_srgw <- do.call(
    rfugw::entropic_semirelaxed_gromov_wasserstein,
    c(ctrl_srgw, list(backend = "cpp", precision = "mixed"))
  )

  rows[[k]] <- bench_one(
    "rfugw_srfgw_cpp_mixed",
    function() do.call(
      rfugw::entropic_semirelaxed_fused_gromov_wasserstein,
      c(ctrl_srfgw, list(backend = "cpp", precision = "mixed"))
    ),
    "srfgw",
    n,
    warm_srfgw$srfgw_dist,
    warm_srfgw$iterations,
    warm_srfgw$error,
    iters
  )
  k <- k + 1L

  rows[[k]] <- bench_one(
    "rfugw_srfgw_cpp_double",
    function() do.call(
      rfugw::entropic_semirelaxed_fused_gromov_wasserstein,
      c(ctrl_srfgw, list(backend = "cpp", precision = "double"))
    ),
    "srfgw",
    n,
    warm_srfgw$srfgw_dist,
    warm_srfgw$iterations,
    warm_srfgw$error,
    iters
  )
  k <- k + 1L

  if (n <= 180L) {
    rows[[k]] <- bench_one(
      "rfugw_srfgw_r",
      function() do.call(
        rfugw::entropic_semirelaxed_fused_gromov_wasserstein,
        c(ctrl_srfgw, list(backend = "r", precision = "double"))
      ),
      "srfgw",
      n,
      warm_srfgw$srfgw_dist,
      warm_srfgw$iterations,
      warm_srfgw$error,
      iters
    )
    k <- k + 1L
  }

  rows[[k]] <- bench_one(
    "rfugw_srgw_cpp_mixed",
    function() do.call(
      rfugw::entropic_semirelaxed_gromov_wasserstein,
      c(ctrl_srgw, list(backend = "cpp", precision = "mixed"))
    ),
    "srgw",
    n,
    warm_srgw$srgw_dist,
    warm_srgw$iterations,
    warm_srgw$error,
    iters
  )
  k <- k + 1L

  rows[[k]] <- bench_one(
    "rfugw_srgw_cpp_double",
    function() do.call(
      rfugw::entropic_semirelaxed_gromov_wasserstein,
      c(ctrl_srgw, list(backend = "cpp", precision = "double"))
    ),
    "srgw",
    n,
    warm_srgw$srgw_dist,
    warm_srgw$iterations,
    warm_srgw$error,
    iters
  )
  k <- k + 1L
}

out <- do.call(rbind, rows)
dir.create(dirname(out_csv), recursive = TRUE, showWarnings = FALSE)
write.csv(out, out_csv, row.names = FALSE)

cat("== rfugw semirelaxed benchmark ==\n")
print(out[, c("suite", "n", "method", "median_ms", "iter_per_sec", "warm_obj", "warm_iters")])
cat(sprintf("\nWrote benchmark CSV: %s\n", out_csv))
