suppressPackageStartupMessages({
  rlib <- Sys.getenv("RFUGW_RLIB", unset = "")
  if (nzchar(rlib)) {
    .libPaths(c(rlib, .libPaths()))
  }
})

args <- commandArgs(trailingOnly = TRUE)
out_dir <- if (length(args) >= 1) args[[1]] else "inst/bench/results/nightly"
seed <- if (length(args) >= 2) as.integer(args[[2]]) else 20260222L
threads <- if (length(args) >= 3) as.integer(args[[3]]) else 2L
if (!is.finite(threads) || threads < 1L) threads <- 2L

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

Sys.setenv(
  OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1",
  OMP_NUM_THREADS = as.character(threads)
)

checks <- list()
add_check <- function(name, ok, value = NA_real_, threshold = NA_real_) {
  checks[[length(checks) + 1L]] <<- data.frame(
    check = name,
    ok = isTRUE(ok),
    value = as.numeric(value),
    threshold = as.numeric(threshold),
    stringsAsFactors = FALSE
  )
}

run_cmd <- function(cmd, args, env = character()) {
  status <- system2(cmd, args = args, env = env)
  if (!identical(status, 0L)) {
    stop(sprintf("Command failed (%s): %s", status, paste(c(cmd, args), collapse = " ")), call. = FALSE)
  }
}

bench_root <- if (file.exists("inst/bench")) "inst/bench" else "rfugw/inst/bench"

# 1) Strict accuracy gate.
source(file.path(bench_root, "accuracy_gate.R"))
run_accuracy_gate(
  fixture_dir = Sys.getenv("RFUGW_FIXTURE_DIR", unset = "inst/extdata/fixtures"),
  stop_on_fail = TRUE,
  verbose = TRUE
)
add_check("accuracy_gate", TRUE, 1, 1)

# 2) Mixed-path micro profile.
profile_csv <- file.path(out_dir, "profile_mixed_path_nightly.csv")
run_cmd(
  file.path(R.home("bin"), "Rscript"),
  c(file.path(bench_root, "profile_mixed_path.R"), "1", profile_csv, as.character(seed), "8", "600", as.character(threads))
)
prof <- read.csv(profile_csv, stringsAsFactors = FALSE)
need_methods <- c("multialign_batch_mixed_dense", "multialign_batch_mixed_auto", "cpp_batch_mixed_dense")
add_check("profile_methods_present", all(need_methods %in% prof$method), sum(need_methods %in% prof$method), length(need_methods))

if (all(need_methods %in% prof$method)) {
  dense_ms <- prof$median_ms[prof$method == "multialign_batch_mixed_dense"][1]
  auto_ms <- prof$median_ms[prof$method == "multialign_batch_mixed_auto"][1]
  cpp_ms <- prof$median_ms[prof$method == "cpp_batch_mixed_dense"][1]
  add_check("auto_not_slower_than_dense_x1.35", auto_ms <= 1.35 * dense_ms, auto_ms / dense_ms, 1.35)
  add_check("cpp_kernel_within_dense_x1.20", cpp_ms <= 1.20 * dense_ms, cpp_ms / dense_ms, 1.20)
  if ("n_lowrank" %in% names(prof)) {
    auto_lr <- prof$n_lowrank[prof$method == "multialign_batch_mixed_auto"][1]
    add_check("mixed_auto_lowrank_disabled", isTRUE(is.na(auto_lr) || auto_lr <= 0), auto_lr, 0)
  }
}

# 3) Larger multiset sanity run.
large_csv <- file.path(out_dir, "benchmark_multiset_large_nightly.csv")
run_cmd(
  file.path(R.home("bin"), "Rscript"),
  c(file.path(bench_root, "benchmark_multiset_large.R"), "1", large_csv, as.character(seed), "1,2", "500,1000"),
  env = c("RFUGW_SKIP_ACCURACY_GATE=1")
)
large <- read.csv(large_csv, stringsAsFactors = FALSE)
add_check("large_objective_finite", all(is.finite(large$objective)), sum(is.finite(large$objective)), nrow(large))

for (nn in unique(large$n)) {
  for (tt in unique(large$thread_count)) {
    dd <- subset(large, n == nn & thread_count == tt & method == "multiset_fixed_cpp_batch")
    dm <- subset(large, n == nn & thread_count == tt & method == "multiset_fixed_cpp_batch_mixed")
    if (nrow(dd) == 1L && nrow(dm) == 1L) {
      ratio <- dm$median_ms[[1]] / dd$median_ms[[1]]
      add_check(sprintf("large_mixed_faster_n%d_t%d", nn, tt), ratio <= 0.75, ratio, 0.75)
    }
  }
}

# 4) POT parity suite.
suite_csv <- file.path(out_dir, "benchmark_suite_nightly.csv")
pot_csv <- file.path(out_dir, "pot_benchmark_nightly.csv")
thread_csv <- file.path(out_dir, "thread_scaling_nightly.csv")
run_cmd(
  file.path(R.home("bin"), "Rscript"),
  c(file.path(bench_root, "benchmark_suite.R"), "1", suite_csv, as.character(seed)),
  env = c("RFUGW_SKIP_ACCURACY_GATE=1")
)
run_cmd(
  "bash",
  c(file.path(bench_root, "run_thread_scaling.sh"), thread_csv, "1", as.character(seed), "1,2"),
  env = c("RFUGW_SKIP_ACCURACY_GATE=1")
)

python_bin <- Sys.which("python")
if (nzchar(python_bin)) {
  run_cmd(
    python_bin,
    c(file.path(bench_root, "benchmark_pot_reference.py"), "1", pot_csv, as.character(seed))
  )
  run_cmd(
    python_bin,
    c(file.path(bench_root, "make_benchmark_report.py"), suite_csv, pot_csv, thread_csv, file.path(out_dir, "benchmark_report_nightly.md"))
  )

  suite <- read.csv(suite_csv, stringsAsFactors = FALSE)
  pot <- read.csv(pot_csv, stringsAsFactors = FALSE)
  compare_one <- function(suite_name, rf_method, pot_method, max_ratio) {
    rf <- subset(suite, suite == suite_name & method == rf_method)
    pt <- subset(pot, suite == suite_name & method == pot_method)
    ok_all <- TRUE
    vals <- c()
    for (n in intersect(rf$n, pt$n)) {
      rf_ms <- rf$median_ms[rf$n == n][1]
      pt_ms <- pt$median_ms[pt$n == n][1]
      ratio <- rf_ms / pt_ms
      vals <- c(vals, ratio)
      if (!is.finite(ratio) || ratio > max_ratio) {
        ok_all <- FALSE
      }
    }
    add_check(sprintf("pot_parity_%s_%s", suite_name, rf_method), ok_all, if (length(vals)) max(vals) else NA_real_, max_ratio)
  }
  compare_one("fgw", "fgw_pgd_mixed", "pot_fgw_pgd", 1.25)
  compare_one("fugw", "fugw_kl", "pot_fugw_kl", 1.25)
} else {
  add_check("pot_parity_skipped_no_python", TRUE, 1, 1)
}

check_tbl <- do.call(rbind, checks)
check_csv <- file.path(out_dir, "nightly_checks.csv")
write.csv(check_tbl, check_csv, row.names = FALSE)
cat("Wrote nightly checks:", check_csv, "\n")
print(check_tbl)

if (any(!check_tbl$ok)) {
  bad <- check_tbl[!check_tbl$ok, , drop = FALSE]
  stop(
    paste0(
      "Nightly guard failed:\n",
      paste(sprintf("- %s (value=%g threshold=%g)", bad$check, bad$value, bad$threshold), collapse = "\n")
    ),
    call. = FALSE
  )
}

cat("Nightly guard passed.\n")
