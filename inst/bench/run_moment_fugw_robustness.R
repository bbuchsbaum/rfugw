args <- commandArgs(trailingOnly = TRUE)
`%||%` <- function(x, y) {
  if (is.null(x) || !length(x)) y else x
}
profile <- if (length(args) >= 1L) args[[1L]] else "local"
output_dir <- if (length(args) >= 2L) {
  args[[2L]]
} else file.path("inst", "bench", "results", "moment_fugw_robustness")
package_mode <- if (length(args) >= 3L) args[[3L]] else "source"

if (!profile %in% c("local", "release")) {
  stop("profile must be `local` or `release`.", call. = FALSE)
}
if (!package_mode %in% c("source", "installed")) {
  stop("package_mode must be `source` or `installed`.", call. = FALSE)
}
required_packages <- c("digest", "jsonlite", "testthat")
missing_packages <- required_packages[!vapply(
  required_packages, requireNamespace, logical(1), quietly = TRUE
)]
if (length(missing_packages)) {
  stop(
    "Moment-FUGW robustness replay requires: ",
    paste(missing_packages, collapse = ", "),
    call. = FALSE
  )
}

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (!length(script_arg)) {
  stop("Could not resolve the robustness runner path.", call. = FALSE)
}
script_path <- normalizePath(
  sub("^--file=", "", script_arg[[1L]]), mustWork = TRUE
)
repository_root <- normalizePath(
  file.path(dirname(script_path), "..", ".."), mustWork = TRUE
)
contract_path <- file.path(
  repository_root, "inst", "bench", "moment-fugw-robustness-contract.json"
)
provenance_path <- file.path(
  repository_root, "inst", "numerical-trust", "provenance-lib.R"
)
source(provenance_path, local = TRUE)
contract <- jsonlite::fromJSON(contract_path, simplifyVector = FALSE)

thread_environment <- unlist(contract$thread_environment, use.names = TRUE)
do.call(Sys.setenv, as.list(thread_environment))

if (identical(package_mode, "source")) {
  if (!requireNamespace("devtools", quietly = TRUE)) {
    stop("Source-mode robustness replay requires devtools.", call. = FALSE)
  }
  devtools::load_all(repository_root, quiet = TRUE)
} else {
  suppressPackageStartupMessages(library(rfugw))
}
suppressPackageStartupMessages(library(testthat))

result_rows <- list()
file_receipts <- list()
for (entry in contract$test_files) {
  relative_path <- as.character(entry$path[[1L]])
  test_path <- file.path(repository_root, relative_path)
  if (!file.exists(test_path)) {
    stop("Missing robustness test file: ", relative_path, call. = FALSE)
  }
  result <- testthat::test_file(test_path, reporter = "silent")
  file_receipts[[length(file_receipts) + 1L]] <- list(
    path = relative_path,
    sha256 = unname(digest::digest(test_path, algo = "sha256", file = TRUE))
  )
  for (case in result) {
    expectations <- case$results %||% list()
    class_count <- function(class_name) {
      sum(vapply(expectations, inherits, logical(1), class_name))
    }
    result_rows[[length(result_rows) + 1L]] <- data.frame(
      file = relative_path,
      case = as.character(case$test),
      passed_expectations = class_count("expectation_success"),
      failures = class_count("expectation_failure"),
      errors = class_count("expectation_error"),
      warnings = class_count("expectation_warning"),
      skips = class_count("expectation_skip"),
      elapsed_seconds = as.numeric(case$real),
      stringsAsFactors = FALSE
    )
  }
}
cases <- do.call(rbind, result_rows)
required_cases <- unique(unlist(lapply(
  contract$test_files,
  function(entry) unlist(entry$required_cases, use.names = FALSE)
), use.names = FALSE))
property_cases <- unlist(contract$property_cases, use.names = TRUE)
case_passed <- function(case_name) {
  rows <- cases[cases$case == case_name, , drop = FALSE]
  nrow(rows) == 1L &&
    sum(rows$failures + rows$errors + rows$warnings + rows$skips) == 0L &&
    rows$passed_expectations[[1L]] > 0L
}
properties <- vapply(property_cases, case_passed, logical(1))
all_required_cases_present <- setequal(
  required_cases,
  intersect(required_cases, cases$case)
)
test_passed <- nrow(cases) >= length(required_cases) &&
  all_required_cases_present &&
  all(properties) &&
  sum(cases$failures + cases$errors + cases$warnings + cases$skips) == 0L &&
  sum(cases$passed_expectations) > 0L

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
cases_path <- file.path(output_dir, "cases.csv")
utils::write.csv(cases, cases_path, row.names = FALSE)

git <- rfugw_git_provenance(repository_root)
expected_commit <- Sys.getenv("GITHUB_SHA", unset = "")
commit_matches_host <- rfugw_valid_commit(expected_commit) &&
  rfugw_valid_commit(git$commit) &&
  identical(tolower(expected_commit), tolower(git$commit))
manifest_path <- file.path(repository_root, "release-artifact.json")
manifest <- if (file.exists(manifest_path)) {
  tryCatch(
    jsonlite::fromJSON(manifest_path, simplifyVector = FALSE),
    error = function(error) list()
  )
} else {
  list()
}
artifact_name <- as.character(manifest$artifact %||% "")[[1L]]
artifact_path <- if (nzchar(artifact_name)) {
  file.path(repository_root, artifact_name)
} else {
  ""
}
artifact_sha256 <- if (nzchar(artifact_path) && file.exists(artifact_path)) {
  unname(digest::digest(artifact_path, algo = "sha256", file = TRUE))
} else {
  NA_character_
}
manifest_sha256 <- if (file.exists(manifest_path)) {
  unname(digest::digest(manifest_path, algo = "sha256", file = TRUE))
} else {
  NA_character_
}
manifest_commit <- as.character(manifest$source_commit %||% "")[[1L]]
manifest_bound <- rfugw_valid_commit(manifest_commit) &&
  rfugw_valid_commit(git$commit) &&
  identical(tolower(manifest_commit), tolower(git$commit)) &&
  nzchar(as.character(manifest$sha256 %||% "")[[1L]]) &&
  identical(
    artifact_sha256,
    as.character(manifest$sha256 %||% "")[[1L]]
  )

installed_version <- as.character(utils::packageVersion("rfugw"))
manifest_version <- as.character(manifest$version %||% "")[[1L]]
admission_evaluable <- test_passed &&
  identical(profile, "release") &&
  identical(package_mode, "installed") &&
  isTRUE(git$exact_commit_evidence) &&
  commit_matches_host && manifest_bound &&
  identical(installed_version, manifest_version)

summary <- list(
  schema_version = 1L,
  method_family = "moment_fugw",
  status = if (test_passed) "passed" else "failed",
  admission_evaluable = admission_evaluable,
  profile = profile,
  package_mode = package_mode,
  candidate = list(
    git_commit = git$commit,
    git_dirty = git$working_tree_dirty,
    git_status_entry_count = git$git_status_entry_count,
    git_provenance_complete = git$git_provenance_complete,
    git_commit_command_status = git$git_commit_command_status,
    git_status_command_status = git$git_status_command_status,
    commit_matches_host = commit_matches_host,
    package_version = installed_version
  ),
  source_artifact = list(
    manifest = if (file.exists(manifest_path)) basename(manifest_path) else NA,
    manifest_sha256 = manifest_sha256,
    artifact = if (nzchar(artifact_name)) artifact_name else NA,
    artifact_sha256 = artifact_sha256,
    manifest_commit = if (nzchar(manifest_commit)) manifest_commit else NA,
    bound = manifest_bound
  ),
  contract = list(
    artifact = basename(contract_path),
    sha256 = unname(digest::digest(contract_path, algo = "sha256", file = TRUE)),
    schema_version = as.integer(contract$schema_version)
  ),
  tests = list(
    artifact = basename(cases_path),
    observed_cases = nrow(cases),
    required_cases = length(required_cases),
    all_required_cases_present = all_required_cases_present,
    passed_expectations = sum(cases$passed_expectations),
    failures = sum(cases$failures),
    errors = sum(cases$errors),
    warnings = sum(cases$warnings),
    skips = sum(cases$skips),
    files = file_receipts
  ),
  properties = as.list(properties),
  environment = list(
    r_version = R.version.string,
    platform = R.version$platform,
    sys_info = as.list(Sys.info()),
    threads = as.list(Sys.getenv(
      names(thread_environment), unset = NA_character_, names = TRUE
    ))
  ),
  global_optimality = "not_claimed"
)
summary_path <- file.path(output_dir, "summary.json")
jsonlite::write_json(
  summary, summary_path, auto_unbox = TRUE, pretty = TRUE,
  null = "null", na = "null", dataframe = "rows"
)
cat("Wrote Moment-FUGW robustness evidence to", output_dir, "\n")
cat("Summary:", summary_path, "\n")

if (!test_passed) {
  stop("Moment-FUGW robustness replay failed.", call. = FALSE)
}
if (identical(profile, "release") && !admission_evaluable) {
  stop("Release robustness evidence is not exact-artifact admissible.",
       call. = FALSE)
}
