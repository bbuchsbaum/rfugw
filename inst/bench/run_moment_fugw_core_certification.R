args <- commandArgs(trailingOnly = TRUE)
`%||%` <- function(x, y) {
  if (is.null(x) || !length(x)) y else x
}
profile <- if (length(args) >= 1L) args[[1L]] else "local"
output_dir <- if (length(args) >= 2L) {
  args[[2L]]
} else file.path("inst", "bench", "results", "moment_fugw_core")
package_mode <- if (length(args) >= 3L) args[[3L]] else "source"

if (!profile %in% c("local", "nightly", "release")) {
  stop("profile must be `local`, `nightly`, or `release`.", call. = FALSE)
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
    "Moment-FUGW core certification requires: ",
    paste(missing_packages, collapse = ", "), call. = FALSE
  )
}

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (!length(script_arg)) {
  stop("Could not resolve the core-certification runner path.", call. = FALSE)
}
script_path <- normalizePath(
  sub("^--file=", "", script_arg[[1L]]), mustWork = TRUE
)
repository_root <- normalizePath(
  file.path(dirname(script_path), "..", ".."), mustWork = TRUE
)
old_working_directory <- setwd(repository_root)
on.exit(setwd(old_working_directory), add = TRUE)
output_dir <- normalizePath(
  output_dir, mustWork = FALSE
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
output_dir <- normalizePath(output_dir, mustWork = TRUE)

contract_path <- file.path(
  repository_root, "inst", "bench", "moment-fugw-core-contract.json"
)
provenance_path <- file.path(
  repository_root, "inst", "numerical-trust", "provenance-lib.R"
)
source(provenance_path, local = TRUE)
contract <- jsonlite::fromJSON(contract_path, simplifyVector = FALSE)
accuracy_design_path <- file.path(
  repository_root, as.character(contract$accuracy$design[[1L]])
)
accuracy_design <- utils::read.csv(
  accuracy_design_path, stringsAsFactors = FALSE
)
thread_environment <- unlist(contract$thread_environment, use.names = TRUE)
do.call(Sys.setenv, as.list(thread_environment))

if (identical(package_mode, "source")) {
  if (!requireNamespace("devtools", quietly = TRUE)) {
    stop("Source-mode core certification requires devtools.", call. = FALSE)
  }
  devtools::load_all(repository_root, quiet = TRUE)
} else {
  suppressPackageStartupMessages(library(rfugw))
}
suppressPackageStartupMessages(library(testthat))

accuracy_path <- file.path(output_dir, "accuracy.csv")
accuracy_script <- file.path(
  repository_root, as.character(contract$accuracy$script[[1L]])
)
accuracy_environment <- if (identical(package_mode, "source")) {
  "RFUGW_DEV_SOURCE=1"
} else {
  "RFUGW_DEV_SOURCE=0"
}
accuracy_status <- system2(
  file.path(R.home("bin"), "Rscript"),
  c(
    accuracy_script,
    accuracy_path,
    as.character(contract$accuracy$profile[[1L]]),
    as.character(contract$accuracy$seed[[1L]])
  ),
  env = accuracy_environment
)
if (!identical(accuracy_status, 0L) || !file.exists(accuracy_path)) {
  stop("The full Moment-FUGW differential matrix failed.", call. = FALSE)
}
accuracy <- utils::read.csv(accuracy_path, stringsAsFactors = FALSE)
accuracy_thresholds <- unlist(
  contract$accuracy$thresholds, use.names = TRUE
)
accuracy_columns <- c(
  "case_id", "ns", "nt", "zero_weight", "near_zero_weight",
  "dense_converged", "implicit_converged", "inner_certified",
  "outer_certified", names(accuracy_thresholds), "pass", "reject_reason"
)
accuracy_pairs <- paste0(accuracy$ns, "x", accuracy$nt)
expected_pairs <- unlist(contract$accuracy$expected_pairs, use.names = FALSE)
direct_rows <- pmax(accuracy$ns, accuracy$nt) <= 8L
design_columns <- names(accuracy_design)
accuracy_design_gate <-
  nrow(accuracy) == nrow(accuracy_design) &&
  all(design_columns %in% names(accuracy)) &&
  all(vapply(design_columns, function(column) {
    isTRUE(all.equal(
      accuracy[[column]], accuracy_design[[column]],
      tolerance = 0, check.attributes = FALSE
    ))
  }, logical(1))) &&
  all(accuracy$seed ==
    as.integer(contract$accuracy$seed) + accuracy_design$seed_offset)
accuracy_gate <-
  accuracy_design_gate &&
  all(accuracy_columns %in% names(accuracy)) &&
  nrow(accuracy) == as.integer(contract$accuracy$rows) &&
  setequal(accuracy_pairs, expected_pairs) &&
  !anyDuplicated(accuracy$case_id) &&
  all(accuracy$pass) && all(accuracy$dense_converged) &&
  all(accuracy$implicit_converged) && all(accuracy$inner_certified) &&
  all(accuracy$outer_certified) &&
  any(accuracy$zero_weight) && any(accuracy$near_zero_weight) &&
  all(!accuracy$zero_weight | !accuracy$near_zero_weight) &&
  all(is.finite(accuracy$direct_oracle_relative_error[direct_rows])) &&
  all(
    accuracy$direct_oracle_relative_error[direct_rows] <=
      accuracy_thresholds[["direct_oracle_relative_error"]]
  ) &&
  all(is.na(accuracy$reject_reason) | !nzchar(accuracy$reject_reason))
if (accuracy_gate) {
  for (metric in names(accuracy_thresholds)) {
    values <- accuracy[[metric]]
    if (identical(metric, "direct_oracle_relative_error")) {
      values <- values[is.finite(values)]
    }
    accuracy_gate <- accuracy_gate && length(values) > 0L &&
      all(is.finite(values)) && all(values <= accuracy_thresholds[[metric]])
  }
}

multiscale_path <- file.path(output_dir, "multiscale.csv")
multiscale_script <- file.path(
  repository_root, as.character(contract$multiscale$script[[1L]])
)
multiscale_environment <- new.env(parent = asNamespace("rfugw"))
sys.source(multiscale_script, envir = multiscale_environment)
multiscale <- multiscale_environment$
  bench_run_moment_fugw_multiscale_efficacy(
    seeds = as.integer(unlist(contract$multiscale$seeds)),
    output = multiscale_path
  )
improved_fraction <- mean(multiscale$inner_iteration_reduction > 0)
median_reduction <- stats::median(multiscale$inner_iteration_reduction)
multiscale_gate <-
  nrow(multiscale) == as.integer(contract$multiscale$rows) &&
  setequal(multiscale$seed, as.integer(unlist(contract$multiscale$seeds))) &&
  !anyDuplicated(multiscale$seed) &&
  all(multiscale$same_final_contract) &&
  all(multiscale$cold_certified) &&
  all(multiscale$transfer_final_certified) &&
  all(multiscale$transfer_multiscale_certified) &&
  all(multiscale$efficacy_gate) &&
  improved_fraction >=
    as.numeric(contract$multiscale$minimum_improved_fraction) &&
  median_reduction >= as.numeric(
    contract$multiscale$minimum_median_inner_iteration_reduction
  ) &&
  all(abs(multiscale$objective_relative_difference) <= as.numeric(
    contract$multiscale$maximum_objective_relative_difference
  )) &&
  all(multiscale$heldout_score_ratio >= as.numeric(
    contract$multiscale$minimum_heldout_score_ratio
  )) &&
  all(multiscale$plan_action_relative_difference <= as.numeric(
    contract$multiscale$maximum_identity_action_difference
  ))

test_environment <- new.env(parent = globalenv())
helper_receipts <- lapply(contract$helper_files, function(relative_path) {
  relative_path <- as.character(relative_path[[1L]])
  path <- file.path(repository_root, relative_path)
  sys.source(path, envir = test_environment)
  list(
    path = relative_path,
    sha256 = unname(digest::digest(path, algo = "sha256", file = TRUE))
  )
})
case_rows <- list()
test_receipts <- list()
for (entry in contract$test_files) {
  relative_path <- as.character(entry$path[[1L]])
  test_path <- file.path(repository_root, relative_path)
  result <- testthat::test_file(
    test_path, reporter = "silent", env = test_environment
  )
  test_receipts[[length(test_receipts) + 1L]] <- list(
    path = relative_path,
    sha256 = unname(digest::digest(test_path, algo = "sha256", file = TRUE))
  )
  for (case in result) {
    expectations <- case$results %||% list()
    class_count <- function(class_name) {
      sum(vapply(expectations, inherits, logical(1), class_name))
    }
    case_rows[[length(case_rows) + 1L]] <- data.frame(
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
cases <- do.call(rbind, case_rows)
cases_path <- file.path(output_dir, "cases.csv")
utils::write.csv(cases, cases_path, row.names = FALSE)
required_cases <- do.call(rbind, lapply(contract$test_files, function(entry) {
  values <- unlist(entry$required_cases, use.names = FALSE)
  data.frame(
    file = rep(as.character(entry$path[[1L]]), length(values)),
    case = as.character(values), stringsAsFactors = FALSE
  )
}))
required_case_pass <- vapply(seq_len(nrow(required_cases)), function(index) {
  hit <- cases$file == required_cases$file[[index]] &
    cases$case == required_cases$case[[index]]
  sum(hit) == 1L && cases$passed_expectations[hit] > 0L &&
    sum(unlist(cases[hit, c(
      "failures", "errors", "warnings", "skips"
    )], use.names = FALSE)) == 0L
}, logical(1))
test_gate <- all(required_case_pass) &&
  sum(cases$failures + cases$errors + cases$warnings + cases$skips) == 0L &&
  sum(cases$passed_expectations) > 0L

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
core_passed <- accuracy_gate && multiscale_gate && test_gate
admission_evaluable <- core_passed && identical(profile, "release") &&
  identical(package_mode, "installed") &&
  isTRUE(git$exact_commit_evidence) && commit_matches_host &&
  manifest_bound && identical(installed_version, manifest_version)

executable_paths <- c(
  as.character(contract$accuracy$script[[1L]]),
  as.character(contract$multiscale$script[[1L]])
)
executables <- lapply(executable_paths, function(relative_path) {
  path <- file.path(repository_root, relative_path)
  list(
    path = relative_path,
    sha256 = unname(digest::digest(path, algo = "sha256", file = TRUE))
  )
})
input_paths <- as.character(contract$accuracy$design[[1L]])
inputs <- lapply(input_paths, function(relative_path) {
  path <- file.path(repository_root, relative_path)
  list(
    path = relative_path,
    sha256 = unname(digest::digest(path, algo = "sha256", file = TRUE))
  )
})
metric_maximum <- function(metric) {
  values <- accuracy[[metric]]
  values <- values[is.finite(values)]
  if (length(values)) max(values) else NA_real_
}
summary <- list(
  schema_version = 1L,
  method_family = "moment_fugw",
  status = if (core_passed) "passed" else "failed",
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
  accuracy = list(
    gate = accuracy_gate,
    artifact = basename(accuracy_path),
    rows = nrow(accuracy),
    zero_weight_rows = sum(accuracy$zero_weight),
    near_zero_weight_rows = sum(accuracy$near_zero_weight),
    maximum_objective_relative_error = metric_maximum(
      "objective_relative_error"
    ),
    maximum_component_relative_error = metric_maximum(
      "component_relative_error"
    ),
    maximum_action_relative_error = metric_maximum("action_relative_error"),
    maximum_adjoint_relative_error = metric_maximum("adjoint_relative_error"),
    maximum_mass_relative_error = metric_maximum("mass_relative_error"),
    maximum_moment_relative_error = metric_maximum("moment_relative_error")
  ),
  multiscale = list(
    gate = multiscale_gate,
    artifact = basename(multiscale_path),
    rows = nrow(multiscale),
    improved_fraction = improved_fraction,
    median_inner_iteration_reduction = median_reduction,
    maximum_objective_relative_difference = max(
      abs(multiscale$objective_relative_difference)
    ),
    minimum_heldout_score_ratio = min(multiscale$heldout_score_ratio),
    maximum_identity_action_difference = max(
      multiscale$plan_action_relative_difference
    )
  ),
  tests = list(
    gate = test_gate,
    artifact = basename(cases_path),
    observed_cases = nrow(cases),
    required_cases = nrow(required_cases),
    passed_expectations = sum(cases$passed_expectations),
    failures = sum(cases$failures),
    errors = sum(cases$errors),
    warnings = sum(cases$warnings),
    skips = sum(cases$skips),
    files = test_receipts,
    helpers = helper_receipts
  ),
  executables = executables,
  inputs = inputs,
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
cat("Wrote Moment-FUGW core certification to", output_dir, "\n")
cat("Summary:", summary_path, "\n")

if (!core_passed) {
  stop("Moment-FUGW core certification failed.", call. = FALSE)
}
if (identical(profile, "release") && !admission_evaluable) {
  stop("Release core certification is not exact-artifact admissible.",
       call. = FALSE)
}
