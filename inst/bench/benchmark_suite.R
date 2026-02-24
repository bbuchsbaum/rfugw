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
iters <- if (length(args) >= 1) as.integer(args[[1]]) else 3L
out_csv <- if (length(args) >= 2) args[[2]] else ""
set.seed(if (length(args) >= 3) as.integer(args[[3]]) else 42L)

make_fgw_problem <- function(ns, nt, d_struct = 3L, d_feat = 5L) {
  X1 <- matrix(rnorm(ns * d_struct), ns, d_struct)
  X2 <- matrix(rnorm(nt * d_struct), nt, d_struct)
  Y1 <- matrix(rnorm(ns * d_feat), ns, d_feat)
  Y2 <- matrix(rnorm(nt * d_feat), nt, d_feat)

  C1 <- as.matrix(dist(X1))
  C2 <- as.matrix(dist(X2))
  C1 <- C1 / max(C1)
  C2 <- C2 / max(C2)

  M <- as.matrix(dist(rbind(Y1, Y2)))[1:ns, (ns + 1):(ns + nt)]
  M <- M / max(M)

  list(C1 = C1, C2 = C2, M = M)
}

make_fugw_problem <- function(n) {
  X1 <- matrix(rnorm(n * 3), n, 3)
  X2 <- matrix(rnorm(n * 3), n, 3)
  Cx <- as.matrix(dist(X1))
  Cy <- as.matrix(dist(X2))
  Cx <- Cx / max(Cx)
  Cy <- Cy / max(Cy)

  M <- matrix(1, nrow = n, ncol = n)
  idx <- seq_len(n)
  M[cbind(idx, rev(idx))] <- 0
  list(Cx = Cx, Cy = Cy, M = M)
}

run_fgw_scaling <- function(iters) {
  sizes <- c(40L, 80L, 120L, 160L)
  rows <- list()
  k <- 1L

  bench_one <- function(label, fn) {
    b <- bench::mark(
      run = fn(),
      iterations = iters,
      check = FALSE,
      memory = TRUE,
      filter_gc = TRUE
    )
    bd <- tibble::as_tibble(b)[1, c("min", "median", "itr/sec", "mem_alloc")]
    data.frame(
      method = label,
      min_ms = as.numeric(bd$min[[1]]) * 1000,
      median_ms = as.numeric(bd$median[[1]]) * 1000,
      iter_per_sec = as.numeric(bd$`itr/sec`[[1]]),
      mem_bytes = as.numeric(bd$mem_alloc[[1]])
    )
  }

  for (n in sizes) {
    p <- make_fgw_problem(n, n)
    warm <- rfugw::fgw_entropic(
      p$M, p$C1, p$C2,
      alpha = 0.5,
      epsilon = 0.05,
      solver = "PGD",
      max_iter = 200,
      sinkhorn_max_iter = 500,
      sinkhorn_tol = 1e-9
    )

    rows[[k]] <- cbind(
      data.frame(
        suite = "fgw", n = n,
        warm_iters = warm$iterations,
        warm_error = warm$error,
        warm_obj = warm$fgw_dist
      ),
      bench_one("fgw_pgd", function() rfugw::fgw_entropic(
        p$M, p$C1, p$C2,
        alpha = 0.5,
        epsilon = 0.05,
        solver = "PGD",
        max_iter = 200,
        sinkhorn_max_iter = 500,
        sinkhorn_tol = 1e-9
      ))
    )
    k <- k + 1L

    rows[[k]] <- cbind(
      data.frame(
        suite = "fgw", n = n,
        warm_iters = warm$iterations,
        warm_error = warm$error,
        warm_obj = warm$fgw_dist
      ),
      bench_one("fgw_pgd_double", function() rfugw::fgw_entropic(
        p$M, p$C1, p$C2,
        alpha = 0.5,
        epsilon = 0.05,
        solver = "PGD",
        precision = "double",
        max_iter = 200,
        sinkhorn_max_iter = 500,
        sinkhorn_tol = 1e-9
      ))
    )
    k <- k + 1L

    rows[[k]] <- cbind(
      data.frame(
        suite = "fgw", n = n,
        warm_iters = warm$iterations,
        warm_error = warm$error,
        warm_obj = warm$fgw_dist
      ),
      bench_one("fgw_pgd_mixed", function() rfugw::fgw_entropic(
        p$M, p$C1, p$C2,
        alpha = 0.5,
        epsilon = 0.05,
        solver = "PGD",
        precision = "mixed",
        max_iter = 200,
        sinkhorn_max_iter = 500,
        sinkhorn_tol = 1e-9
      ))
    )
    k <- k + 1L

    if (n >= 80L) {
      rows[[k]] <- cbind(
        data.frame(
          suite = "fgw", n = n,
          warm_iters = warm$iterations,
          warm_error = warm$error,
          warm_obj = warm$fgw_dist
        ),
        bench_one("fgw_pgd_mixed_lr64", function() rfugw::fgw_entropic(
          p$M, p$C1, p$C2,
          alpha = 0.5,
          epsilon = 0.05,
          solver = "PGD",
          precision = "mixed",
          structure_rank = min(64L, n),
          max_iter = 200,
          sinkhorn_max_iter = 500,
          sinkhorn_tol = 1e-9
        ))
      )
      k <- k + 1L
    }

    rows[[k]] <- cbind(
      data.frame(
        suite = "fgw", n = n,
        warm_iters = warm$iterations,
        warm_error = warm$error,
        warm_obj = warm$fgw_dist
      ),
      bench_one("fgw_pgd_log", function() rfugw::fgw_entropic(
        p$M, p$C1, p$C2,
        alpha = 0.5,
        epsilon = 0.05,
        solver = "PGD",
        sinkhorn_method = "log",
        max_iter = 200,
        sinkhorn_max_iter = 500,
        sinkhorn_tol = 1e-9
      ))
    )
    k <- k + 1L

    rows[[k]] <- cbind(
      data.frame(
        suite = "fgw", n = n,
        warm_iters = warm$iterations,
        warm_error = warm$error,
        warm_obj = warm$fgw_dist
      ),
      bench_one("fgw_ppa", function() rfugw::fgw_entropic(
        p$M, p$C1, p$C2,
        alpha = 0.5,
        epsilon = 0.05,
        solver = "PPA",
        max_iter = 200,
        sinkhorn_max_iter = 500,
        sinkhorn_tol = 1e-9
      ))
    )
    k <- k + 1L

    if (n <= 80L) {
      rows[[k]] <- cbind(
        data.frame(
          suite = "fgw", n = n,
          warm_iters = warm$iterations,
          warm_error = warm$error,
          warm_obj = warm$fgw_dist
        ),
        bench_one("fgw_exact_cg", function() rfugw::fgw_exact_cg(
          p$M, p$C1, p$C2,
          alpha = 0.5,
          symmetric = TRUE,
          max_iter = 120,
          tol_rel = 1e-9,
          tol_abs = 1e-9,
          lp_solver = "cpp_transport"
        ))
      )
      k <- k + 1L
    }
  }

  do.call(rbind, rows)
}

run_fugw_scaling <- function(iters) {
  sizes <- c(30L, 50L, 70L, 90L)
  rows <- list()
  k <- 1L
  include_mixed <- identical(Sys.getenv("RFUGW_BENCH_FUGW_MIXED", unset = "0"), "1")

  bench_one <- function(label, fn) {
    b <- bench::mark(
      run = fn(),
      iterations = iters,
      check = FALSE,
      memory = TRUE,
      filter_gc = TRUE
    )
    bd <- tibble::as_tibble(b)[1, c("min", "median", "itr/sec", "mem_alloc")]
    data.frame(
      method = label,
      min_ms = as.numeric(bd$min[[1]]) * 1000,
      median_ms = as.numeric(bd$median[[1]]) * 1000,
      iter_per_sec = as.numeric(bd$`itr/sec`[[1]]),
      mem_bytes = as.numeric(bd$mem_alloc[[1]])
    )
  }

  for (i in seq_along(sizes)) {
    n <- sizes[[i]]
    p <- make_fugw_problem(n)
    warm_d <- rfugw::fugw_kl(
      p$Cx, p$Cy,
      reg_marginals = c(100, 50),
      epsilon = 1e-2,
      alpha = 0.5,
      M = p$M,
      precision = "double",
      max_iter = 80,
      max_iter_ot = 300,
      tol = 1e-8,
      tol_ot = 1e-8
    )

    rows[[k]] <- cbind(
      data.frame(
        suite = "fugw",
        n = n,
        warm_iters = warm_d$iterations,
        warm_error = warm_d$error,
        warm_obj = warm_d$fugw_cost,
        warm_inner_total = if (!is.null(warm_d$inner_iters_total)) warm_d$inner_iters_total else NA_integer_,
        warm_inner_feat_mean = if (!is.null(warm_d$inner_iters_feat)) mean(warm_d$inner_iters_feat) else NA_real_,
        warm_inner_samp_mean = if (!is.null(warm_d$inner_iters_samp)) mean(warm_d$inner_iters_samp) else NA_real_,
        warm_inner_warm_frac_feat = if (!is.null(warm_d$inner_warm_feat)) mean(as.numeric(warm_d$inner_warm_feat)) else NA_real_,
        warm_inner_warm_frac_samp = if (!is.null(warm_d$inner_warm_samp)) mean(as.numeric(warm_d$inner_warm_samp)) else NA_real_,
        warm_inner_fallback_frac_feat = if (!is.null(warm_d$inner_warm_fallback_feat)) mean(as.numeric(warm_d$inner_warm_fallback_feat)) else NA_real_,
        warm_inner_fallback_frac_samp = if (!is.null(warm_d$inner_warm_fallback_samp)) mean(as.numeric(warm_d$inner_warm_fallback_samp)) else NA_real_
      ),
      bench_one("fugw_kl", function() rfugw::fugw_kl(
        p$Cx, p$Cy,
        reg_marginals = c(100, 50),
        epsilon = 1e-2,
        alpha = 0.5,
        M = p$M,
        precision = "double",
        max_iter = 80,
        max_iter_ot = 300,
        tol = 1e-8,
        tol_ot = 1e-8
      ))
    )
    k <- k + 1L

    if (include_mixed) {
      warm_m <- rfugw::fugw_kl(
        p$Cx, p$Cy,
        reg_marginals = c(100, 50),
        epsilon = 1e-2,
        alpha = 0.5,
        M = p$M,
        precision = "mixed",
        max_iter = 80,
        max_iter_ot = 300,
        tol = 1e-8,
        tol_ot = 1e-8
      )

      rows[[k]] <- cbind(
        data.frame(
          suite = "fugw",
          n = n,
          warm_iters = warm_m$iterations,
          warm_error = warm_m$error,
          warm_obj = warm_m$fugw_cost,
          warm_inner_total = if (!is.null(warm_m$inner_iters_total)) warm_m$inner_iters_total else NA_integer_,
          warm_inner_feat_mean = if (!is.null(warm_m$inner_iters_feat)) mean(warm_m$inner_iters_feat) else NA_real_,
          warm_inner_samp_mean = if (!is.null(warm_m$inner_iters_samp)) mean(warm_m$inner_iters_samp) else NA_real_,
          warm_inner_warm_frac_feat = if (!is.null(warm_m$inner_warm_feat)) mean(as.numeric(warm_m$inner_warm_feat)) else NA_real_,
          warm_inner_warm_frac_samp = if (!is.null(warm_m$inner_warm_samp)) mean(as.numeric(warm_m$inner_warm_samp)) else NA_real_,
          warm_inner_fallback_frac_feat = if (!is.null(warm_m$inner_warm_fallback_feat)) mean(as.numeric(warm_m$inner_warm_fallback_feat)) else NA_real_,
          warm_inner_fallback_frac_samp = if (!is.null(warm_m$inner_warm_fallback_samp)) mean(as.numeric(warm_m$inner_warm_fallback_samp)) else NA_real_
        ),
        bench_one("fugw_kl_mixed", function() rfugw::fugw_kl(
          p$Cx, p$Cy,
          reg_marginals = c(100, 50),
          epsilon = 1e-2,
          alpha = 0.5,
          M = p$M,
          precision = "mixed",
          max_iter = 80,
          max_iter_ot = 300,
          tol = 1e-8,
          tol_ot = 1e-8
        ))
      )
      k <- k + 1L
    }
  }

  do.call(rbind, rows)
}

fgw_res <- run_fgw_scaling(iters)
fugw_res <- run_fugw_scaling(iters)
all_cols <- union(names(fgw_res), names(fugw_res))
for (nm in setdiff(all_cols, names(fgw_res))) fgw_res[[nm]] <- NA
for (nm in setdiff(all_cols, names(fugw_res))) fugw_res[[nm]] <- NA
all_res <- rbind(fgw_res[all_cols], fugw_res[all_cols])

if (nzchar(out_csv)) {
  dir.create(dirname(out_csv), recursive = TRUE, showWarnings = FALSE)
  write.csv(all_res, out_csv, row.names = FALSE)
}

cat("== FGW scaling ==\n")
print(fgw_res)
cat("\n== FUGW scaling ==\n")
print(fugw_res)
if (nzchar(out_csv)) {
  cat("\nWrote benchmark CSV:", out_csv, "\n")
}
