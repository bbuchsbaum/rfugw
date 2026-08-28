#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
profile <- if (length(args) >= 1L) args[[1L]] else "smoke"
results_root <- if (length(args) >= 2L) {
  args[[2L]]
} else file.path("inst", "bench", "results", "moment_fugw")
package_mode <- if (length(args) >= 3L) args[[3L]] else "source"

if (identical(package_mode, "source")) {
  if (!requireNamespace("pkgload", quietly = TRUE)) {
    stop("Source mode requires pkgload.", call. = FALSE)
  }
  pkgload::load_all(".", quiet = TRUE, export_all = FALSE)
} else if (identical(package_mode, "installed")) {
  suppressPackageStartupMessages(library(rfugw))
} else {
  stop("package_mode must be source or installed.", call. = FALSE)
}

protocol <- if (file.exists(file.path("inst", "bench", "protocol.R"))) {
  file.path("inst", "bench", "protocol.R")
} else {
  system.file("bench", "protocol.R", package = "rfugw")
}
benchmark <- if (file.exists(file.path(
    "inst", "bench", "benchmark_moment_fugw.R"))) {
  file.path("inst", "bench", "benchmark_moment_fugw.R")
} else {
  system.file("bench", "benchmark_moment_fugw.R", package = "rfugw")
}
source(protocol, local = .GlobalEnv)
source(benchmark, local = .GlobalEnv)

result <- bench_run_moment_fugw_protocol(
  profile = profile,
  results_root = results_root,
  package_mode = package_mode,
  package_root = "."
)
print(result$summary)
cat("Moment-FUGW artifacts:", normalizePath(result$out_dir), "\n")
