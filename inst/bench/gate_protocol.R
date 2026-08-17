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
  "",
  "| suite | method | n | valid | solve_ms | cap_ms |",
  "|---|---|---|---|---|---|"
)

solver_reasons <- character()
for (i in seq_len(nrow(runs))) {
  row <- runs[i, ]
  cap <- find_cap(row$suite, row$method, row$n)
  cap_txt <- if (is.finite(cap)) sprintf("%.0f", cap) else ""
  lines <- c(lines, sprintf(
    "| %s | %s | %s | %s | %s | %s |",
    row$suite, row$method, row$n, row$valid,
    if (is.finite(row$solve_ms)) sprintf("%.1f", row$solve_ms) else "",
    cap_txt
  ))
  if (!isTRUE(as.logical(row$valid))) {
    solver_reasons <- c(
      solver_reasons,
      sprintf("%s/%s n=%s: %s", row$suite, row$method, row$n, row$reject_reason)
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
cat("Flagship protocol gate passed.\n")
