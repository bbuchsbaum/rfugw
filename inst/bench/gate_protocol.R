#!/usr/bin/env Rscript
# Quality-first protocol gate. Exit 0 ok, 1 solver/regression, 2 infrastructure.

args <- commandArgs(trailingOnly = TRUE)
out_dir <- if (length(args) >= 1L) args[[1]] else "inst/bench/results/current"
caps_path <- if (length(args) >= 2L) args[[2]] else "inst/bench/ci_time_caps.json"
scale <- if (length(args) >= 3L) args[[3]] else Sys.getenv("RFUGW_PROTOCOL_SCALE", unset = "pr")
`%||%` <- function(x, y) if (is.null(x)) y else x

infra_fail <- function(...) {
  message("INFRA: ", paste(..., sep = ""))
  quit(status = 2L)
}
solver_fail <- function(...) {
  message("SOLVER: ", paste(..., sep = ""))
  quit(status = 1L)
}

if (!dir.exists(out_dir)) {
  infra_fail("missing output directory: ", out_dir)
}
need <- c("meta.json", "runs.csv", "quality.csv")
for (fn in need) {
  if (!file.exists(file.path(out_dir, fn))) {
    infra_fail("protocol did not write ", fn)
  }
}

meta <- tryCatch(
  jsonlite::fromJSON(file.path(out_dir, "meta.json"), simplifyVector = TRUE),
  error = function(e) NULL
)
runs <- tryCatch(utils::read.csv(file.path(out_dir, "runs.csv"), stringsAsFactors = FALSE), error = function(e) NULL)
qual <- tryCatch(utils::read.csv(file.path(out_dir, "quality.csv"), stringsAsFactors = FALSE), error = function(e) NULL)
if (is.null(meta) || is.null(runs) || is.null(qual)) {
  infra_fail("unreadable protocol artifacts")
}
if (!nrow(runs)) {
  infra_fail("runs.csv has no rows")
}
required_run_columns <- c(
  "evidence_class", "certified", "comparison_eligible",
  "performance_regression_eligible", "prepare_ms", "setup_ms", "solve_ms", "e2e_ms",
  "mem_solve_bytes", "mem_e2e_bytes", "requested_precision",
  "effective_precision", "requested_threads", "used_threads"
)
missing_run_columns <- setdiff(required_run_columns, names(runs))
if (length(missing_run_columns)) {
  infra_fail("runs.csv missing evidence columns: ", paste(missing_run_columns, collapse = ", "))
}
if (nrow(qual) != nrow(runs) || !"status" %in% names(qual)) {
  infra_fail("quality.csv must align one-for-one with runs.csv and include status")
}
thresholds_path <- file.path(dirname(caps_path), "thresholds.json")
if (!file.exists(thresholds_path)) {
  infra_fail("missing thresholds.json beside time caps")
}
thresholds_md5 <- unname(tools::md5sum(thresholds_path)[[1]])
if (!identical(as.character(meta$thresholds_md5 %||% ""), thresholds_md5)) {
  infra_fail("benchmark artifact threshold digest does not match this checkout")
}

caps <- list()
if (file.exists(caps_path)) {
  all_caps <- tryCatch(
    jsonlite::fromJSON(caps_path, simplifyVector = FALSE),
    error = function(e) NULL
  )
  if (is.null(all_caps)) {
    infra_fail("unreadable time-cap file: ", caps_path)
  }
  caps <- all_caps[[scale]] %||% list()
}

find_cap <- function(suite, method, n) {
  for (row in caps) {
    if (identical(row$suite, suite) &&
        identical(row$method, method) &&
        identical(as.integer(row$n), as.integer(n))) {
      return(as.numeric(row$solve_ms_max))
    }
  }
  NA_real_
}

lines <- c(
  "# Flagship protocol gate",
  "",
  sprintf("- scale: `%s`", scale),
  sprintf("- commit: `%s`", meta$commit %||% NA),
  sprintf("- seed: %s", meta$seed %||% NA),
  sprintf("- threads: %s", meta$threads %||% NA),
  sprintf("- profile: `%s`", meta$profile %||% NA),
  sprintf("- suites: `%s`", meta$suites %||% NA),
  sprintf("- certified comparison rows: %d", sum(runs$evidence_class == "certified_comparison")),
  sprintf("- fixed-budget performance rows: %d (not certification)", sum(runs$evidence_class == "fixed_budget_performance")),
  "",
  "| suite | method | n | evidence | certified | comparison eligible | performance eligible | solve_ms | cap_ms |",
  "|---|---|---|---|---|---|---|---|---|"
)

solver_reasons <- character()
for (i in seq_len(nrow(runs))) {
  row <- runs[i, ]
  qrow <- qual[i, ]
  cap <- find_cap(row$suite, row$method, row$n)
  is_certified_class <- identical(as.character(row$evidence_class), "certified_comparison")
  is_performance_class <- identical(as.character(row$evidence_class), "fixed_budget_performance")
  eligible <- if (is_certified_class) {
    isTRUE(as.logical(row$comparison_eligible))
  } else if (is_performance_class) {
    isTRUE(as.logical(row$performance_regression_eligible))
  } else {
    FALSE
  }
  cap_txt <- if (is.finite(cap)) sprintf("%.0f", cap) else ""
  lines <- c(lines, sprintf(
    "| %s | %s | %s | %s | %s | %s | %s | %s | %s |",
    row$suite, row$method, row$n, row$evidence_class, row$certified,
    row$comparison_eligible, row$performance_regression_eligible,
    if (is.finite(row$solve_ms)) sprintf("%.1f", row$solve_ms) else "",
    cap_txt
  ))
  if (is_certified_class && !identical(as.character(qrow$status), "converged")) {
    solver_reasons <- c(
      solver_reasons,
      sprintf("%s/%s n=%s: certified row has status=%s",
              row$suite, row$method, row$n, qrow$status)
    )
  }
  if (!eligible) {
    solver_reasons <- c(
      solver_reasons,
      sprintf("%s/%s n=%s [%s]: %s", row$suite, row$method, row$n,
              row$evidence_class, row$reject_reason)
    )
  } else if (is.finite(cap) && is.finite(row$solve_ms) && row$solve_ms > cap) {
    solver_reasons <- c(
      solver_reasons,
      sprintf("%s/%s n=%s: solve_ms %.1f exceeded cap %.0f",
              row$suite, row$method, row$n, row$solve_ms, cap)
    )
  }
}

report <- file.path(out_dir, "gate_report.md")
writeLines(lines, report)
cat(paste(lines, collapse = "\n"), "\n")
cat("Wrote ", report, "\n", sep = "")

if (length(solver_reasons)) {
  solver_fail(paste(solver_reasons, collapse = "\n"))
}
cat(sprintf(
  "Protocol gate passed: %d certified comparison row(s), %d non-certified fixed-budget row(s).\n",
  sum(runs$evidence_class == "certified_comparison"),
  sum(runs$evidence_class == "fixed_budget_performance")
))
