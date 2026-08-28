args <- commandArgs(trailingOnly = TRUE)
split <- if (length(args) >= 1L) args[[1L]] else "evaluation"
output_dir <- if (length(args) >= 2L) {
  args[[2L]]
} else file.path("inst", "bench", "results", "moment_fugw_scientific", split)
package_mode <- if (length(args) >= 3L) args[[3L]] else "source"

if (!split %in% c("tuning", "evaluation")) {
  stop("split must be `tuning` or `evaluation`.", call. = FALSE)
}
if (!package_mode %in% c("source", "installed")) {
  stop("package_mode must be `source` or `installed`.", call. = FALSE)
}
if (!requireNamespace("jsonlite", quietly = TRUE)) {
  stop("Scientific validation requires jsonlite.", call. = FALSE)
}

Sys.setenv(
  OMP_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1"
)

if (identical(package_mode, "source")) {
  if (!requireNamespace("devtools", quietly = TRUE)) {
    stop("Source-mode validation requires devtools.", call. = FALSE)
  }
  devtools::load_all(quiet = TRUE)
} else {
  suppressPackageStartupMessages(library(rfugw))
}

bench_root <- if (identical(package_mode, "installed")) {
  system.file("bench", package = "rfugw")
} else if (file.exists(file.path("inst", "bench"))) {
  file.path("inst", "bench")
} else {
  system.file("bench", package = "rfugw")
}
if (!nzchar(bench_root)) {
  stop("Could not locate installed Moment-FUGW benchmark resources.", call. = FALSE)
}
validation_path <- file.path(bench_root, "moment_fugw_scientific_validation.R")
validation <- new.env(parent = asNamespace("rfugw"))
sys.source(validation_path, envir = validation)

protocol <- validation$moment_validation_protocol()
if (identical(split, "evaluation") &&
    !identical(protocol$status, "frozen_before_evaluation")) {
  stop("Evaluation requires a frozen_before_evaluation protocol.", call. = FALSE)
}
image_seeds <- protocol$image[[paste0(split, "_seeds")]]
surface_seeds <- protocol$surface[[paste0(split, "_seeds")]]
if (!length(image_seeds) || !length(surface_seeds)) {
  stop("The requested split has no frozen seeds.", call. = FALSE)
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
image_path <- file.path(output_dir, paste0("moment-fugw-image-", split, ".csv"))
surface_path <- file.path(
  output_dir, paste0("moment-fugw-surface-", split, ".csv")
)
rank_path <- file.path(
  output_dir, paste0("moment-fugw-geometry-rank-", split, ".csv")
)

started_at <- format(Sys.time(), tz = "UTC", usetz = TRUE)
image <- validation$moment_run_image_suite(image_seeds, image_path)
surface <- validation$moment_run_surface_suite(
  surface_seeds, surface_path, rank_path
)
image_summary <- attr(image, "summary")
surface_parity <- attr(surface, "dense_parity")
rank_curve <- attr(surface, "rank_curve")

image_gate <- nrow(image) ==
  length(image_seeds) * length(protocol$image$overlap_percent) * 2L &&
  all(image$solver_success) && all(image$suite_gate) &&
  all(image_summary$paired_complete) &&
  all(image_summary$false_mass_gate) &&
  all(image_summary$endpoint_gate) &&
  all(image_summary$recovered_mass_gate)
surface_gate <- nrow(surface) == length(surface_seeds) &&
  all(surface$solver_success) && all(surface$suite_gate) &&
  isTRUE(surface_parity$passed)
rank_gate <- nrow(rank_curve) == length(protocol$surface$embedding_curve$ranks) &&
  all(rank_curve$solver_success) &&
  all(is.finite(rank_curve$geometry_relative_error)) &&
  all(is.finite(rank_curve$heldout_correlation)) &&
  all(rank_curve$numerical_and_geometry_error_separate)
fixture_provenance_path <- validation$moment_validation_resource(
  protocol$surface$fixture_provenance_resource
)

git <- validation$moment_git_provenance(".")
scientific_passed <- image_gate && surface_gate && rank_gate
admission_evaluable <- scientific_passed &&
  identical(split, "evaluation") &&
  identical(package_mode, "installed") &&
  isTRUE(git$git_provenance_complete) &&
  identical(git$git_dirty, FALSE)

summary <- list(
  schema_version = 1L,
  status = if (scientific_passed) "passed" else "failed",
  admission_evaluable = admission_evaluable,
  split = split,
  package_mode = package_mode,
  started_at_utc = started_at,
  completed_at_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  protocol = list(
    schema_version = protocol$schema_version,
    status = protocol$status,
    md5 = unname(tools::md5sum(
      validation$moment_validation_resource(
        "moment-fugw-validation-protocol.json"
      )
    )),
    evaluation_seed_accessed_before_freeze =
      protocol$pre_evaluation_amendment$evaluation_seed_accessed_before_freeze
  ),
  candidate = list(
    git_commit = git$git_commit,
    git_dirty = git$git_dirty,
    git_status_entry_count = git$git_status_entry_count,
    git_provenance_complete = git$git_provenance_complete,
    git_commit_command_status = git$git_commit_command_status,
    git_status_command_status = git$git_status_command_status,
    package_version = as.character(utils::packageVersion("rfugw"))
  ),
  environment = list(
    r_version = R.version.string,
    platform = R.version$platform,
    sys_info = as.list(Sys.info()),
    omp_threads = Sys.getenv("OMP_NUM_THREADS", unset = NA_character_),
    blas_threads = list(
      openblas = Sys.getenv("OPENBLAS_NUM_THREADS", unset = NA_character_),
      mkl = Sys.getenv("MKL_NUM_THREADS", unset = NA_character_),
      veclib = Sys.getenv("VECLIB_MAXIMUM_THREADS", unset = NA_character_)
    )
  ),
  image = list(
    gate = image_gate,
    requested_rows = length(image_seeds) *
      length(protocol$image$overlap_percent) * 2L,
    completed_rows = nrow(image),
    failed_rows = sum(!image$solver_success),
    summary = image_summary,
    artifact = basename(image_path)
  ),
  surface = list(
    gate = surface_gate,
    requested_rows = length(surface_seeds),
    completed_rows = nrow(surface),
    failed_rows = sum(!surface$solver_success),
    fixture = list(
      mode = protocol$surface$fixture_resource_mode,
      evaluation_fingerprint_md5 =
        protocol$surface$fixture_fingerprint_md5,
      dense_parity_fingerprint_md5 =
        protocol$surface$dense_parity_fixture_fingerprint_md5,
      provenance_resource = basename(fixture_provenance_path),
      provenance_resource_md5 = unname(tools::md5sum(
        fixture_provenance_path
      ))
    ),
    dense_parity = surface_parity,
    artifact = basename(surface_path)
  ),
  geometry_rank_curve = list(
    gate = rank_gate,
    completed_rows = nrow(rank_curve),
    artifact = basename(rank_path)
  ),
  global_optimality = "not_claimed"
)
summary_path <- file.path(
  output_dir, paste0("moment-fugw-scientific-summary-", split, ".json")
)
jsonlite::write_json(summary, summary_path, auto_unbox = TRUE, pretty = TRUE,
                     na = "null")
cat("Wrote scientific validation evidence to", output_dir, "\n")
cat("Summary:", summary_path, "\n")

if (!identical(summary$status, "passed")) {
  stop("Moment-FUGW scientific validation failed.", call. = FALSE)
}
