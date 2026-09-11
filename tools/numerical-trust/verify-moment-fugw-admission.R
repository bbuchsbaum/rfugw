#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
values <- function(prefix) {
  hits <- args[startsWith(args, prefix)]
  sub(prefix, "", hits, fixed = TRUE)
}
value <- function(prefix, default = NULL) {
  hits <- values(prefix)
  if (!length(hits)) default else hits[[1L]]
}
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- normalizePath(
  sub("^--file=", "", script_arg[[1L]]), mustWork = TRUE
)
repository_root <- normalizePath(
  file.path(dirname(script_path), "..", ".."), mustWork = TRUE
)
source(file.path(
  repository_root, "inst", "bench", "moment_fugw_admission.R"
))

mode <- value("--mode=", "candidate")
manifest <- value("--manifest=")
evidence_root <- value("--evidence-root=")
review <- value("--review=")
output <- value(
  "--output=",
  file.path(".gate", "moment-fugw-admission", paste0(mode, "-receipt.json"))
)
expected_commit <- value("--commit=", Sys.getenv("GITHUB_SHA", ""))
if (is.null(manifest) || is.null(evidence_root)) {
  stop("Required arguments: --manifest=PATH --evidence-root=DIR",
       call. = FALSE)
}
if (!nzchar(expected_commit)) expected_commit <- NULL
receipt <- moment_fugw_verify_admission(
  manifest_path = manifest,
  evidence_root = evidence_root,
  review_path = review,
  mode = mode,
  expected_commit = expected_commit,
  output_path = output
)
cat(sprintf(
  "Moment-FUGW admission: %s; receipt: %s\n",
  receipt$status, normalizePath(output, mustWork = TRUE)
))
if (!isTRUE(receipt$passed)) quit(status = 1L)
