#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- if (length(script_arg)) {
  normalizePath(sub("^--file=", "", script_arg[[1L]]), mustWork = TRUE)
} else {
  normalizePath("tools/numerical-trust/collect-evidence.R", mustWork = TRUE)
}
repository_root <- normalizePath(
  file.path(dirname(script_path), "..", ".."), mustWork = TRUE
)
source(file.path(
  repository_root, "inst", "numerical-trust", "provenance-lib.R"
))
value_arg <- function(prefix, default) {
  hit <- args[startsWith(args, prefix)]
  if (!length(hit)) return(default)
  sub(prefix, "", hit[[1L]], fixed = TRUE)
}
output_dir <- value_arg("--output=", file.path(".gate", "numerical-trust"))
scope <- value_arg("--scope=", Sys.getenv("RFUGW_TRUST_SCOPE", "pr"))
channel <- value_arg("--channel=", "local")
installed <- "--installed" %in% args
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

if (installed) {
  library(rfugw)
} else if (requireNamespace("devtools", quietly = TRUE)) {
  devtools::load_all(".", quiet = TRUE)
} else {
  library(rfugw)
}

digest_file <- function(path) {
  if (!file.exists(path)) return(NULL)
  sha <- if (requireNamespace("digest", quietly = TRUE)) {
    digest::digest(path, algo = "sha256", file = TRUE)
  } else if (nzchar(Sys.which("shasum"))) {
    sub("\\s.*$", "", system2("shasum", c("-a", "256", path), stdout = TRUE)[1])
  } else if (nzchar(Sys.which("sha256sum"))) {
    sub("\\s.*$", "", system2("sha256sum", path, stdout = TRUE)[1])
  } else {
    stop("SHA-256 evidence requires digest, shasum, or sha256sum.", call. = FALSE)
  }
  list(path = path, sha256 = sha)
}

M <- matrix(c(0, 1, 1, 0), 2, 2)
emd <- ot_emd(M)
fgw <- fgw_entropic(
  M, M, M, epsilon = 0.1, max_iter = 20L, tol = 1e-7,
  sinkhorn_tol = 1e-8, sinkhorn_method = "auto",
  precision = "strict_double"
)
certificate_fields <- function(x) {
  keep <- c(
    "formulation", "backend", "status", "converged", "termination_reason",
    "residual", "row_residual", "col_residual", "feasible",
    "feasibility_residual", "feasibility_tolerance", "inner_converged",
    "inner_residual", "max_inner_residual", "objective_recomputed",
    "objective_residual", "objective_consistent", "requested_tol",
    "effective_tol", "requested_inner_tol", "effective_inner_tol",
    "requested_precision", "effective_precision", "compute_precision",
    "requested_sinkhorn_method", "effective_sinkhorn_method",
    "backend_transition", "sinkhorn_backend_transition"
  )
  x[intersect(keep, names(x))]
}

tarballs <- list.files(pattern = "^rfugw_[^/]+\\.tar\\.gz$", full.names = TRUE)
path_matrix_file <- if (file.exists("inst/numerical-path-matrix.csv")) {
  "inst/numerical-path-matrix.csv"
} else {
  system.file("numerical-path-matrix.csv", package = "rfugw")
}
if (!nzchar(path_matrix_file)) {
  stop("Cannot locate the installed numerical-path matrix.", call. = FALSE)
}
path_matrix <- utils::read.csv(path_matrix_file, stringsAsFactors = FALSE)
capabilities <- transport_capabilities()
support_evidence <- getFromNamespace(
  ".release_support_evidence", "rfugw"
)(capabilities, path_matrix)
selected_env <- c(
  "RFUGW_TRUST_SCOPE", "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS",
  "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "RFUGW_OPENMP_FLAGS",
  "BLIS_NUM_THREADS", "RFUGW_FAST_FLAGS", "RFUGW_EXTRA_CXXFLAGS",
  "RFUGW_EXTRA_LIBS", "RFUGW_OPENMP_LIBS"
)
verified_support <- support_evidence$verified_support
experimental_boundaries <- support_evidence$experimental_boundaries
git <- rfugw_git_provenance(".")
working_tree_dirty <- git$working_tree_dirty
expected_host_commit <- Sys.getenv("GITHUB_SHA", "")
host_commit_valid <- rfugw_valid_commit(expected_host_commit)
commit_matches_host <- host_commit_valid &&
  isTRUE(git$git_provenance_complete) &&
  identical(git$commit, tolower(expected_host_commit))
exact_commit_evidence <- isTRUE(git$exact_commit_evidence) &&
  (!identical(channel, "hosted") || commit_matches_host)
limitations <- c(
  if (identical(channel, "local"))
    "Local evidence only: the hosted Linux/macOS/Windows matrix and sanitizer jobs are not evaluated by this bundle.",
  if (!isTRUE(git$git_provenance_complete))
    "Git commit or status provenance is unavailable, so exact-commit evidence is false.",
  if (isTRUE(working_tree_dirty))
    "The source identity is a dirty working tree, so this is not exact-release-commit evidence.",
  if (identical(channel, "hosted") && !commit_matches_host)
    "The checked-out commit does not match a valid 40-hex GITHUB_SHA.",
  "Publication status is not evaluated by the numerical-trust workflow.",
  if (!isTRUE(fgw$runtime_provenance$native_effective$openmp_available))
    "This machine has no active OpenMP runtime; multi-thread execution requires hosted or toolchain-specific evidence."
)
report <- list(
  schema_version = 2L,
  package = "rfugw",
  scope = scope,
  evidence_channel = channel,
  timestamp_utc = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
  commit = git$commit,
  working_tree_dirty = working_tree_dirty,
  git_provenance_complete = git$git_provenance_complete,
  git_status_entry_count = git$git_status_entry_count,
  git_commit_command_status = git$git_commit_command_status,
  git_status_command_status = git$git_status_command_status,
  expected_host_commit = if (host_commit_valid) {
    tolower(expected_host_commit)
  } else NA_character_,
  commit_matches_host = commit_matches_host,
  exact_commit_evidence = exact_commit_evidence,
  hosted_run = if (identical(channel, "hosted")) Sys.getenv("GITHUB_RUN_ID", NA_character_) else NA_character_,
  publication_status = "not_evaluated_by_numerical_trust_job",
  sys_info = as.list(Sys.info()),
  seed_contract = list(
    base = 20260820L,
    transport = "410000 + case_id",
    replay = "RFUGW_REPLAY_SEED=<seed> Rscript tools/numerical-trust/run-laws.R --family=transport"
  ),
  environment = as.list(Sys.getenv(selected_env, unset = NA_character_)),
  runtime_provenance = fgw$runtime_provenance,
  session = capture.output(sessionInfo()),
  representative_certificates = list(exact_transport = certificate_fields(emd),
                                     entropic_fgw = certificate_fields(fgw)),
  capabilities = capabilities,
  verified_support = verified_support,
  experimental_boundaries = experimental_boundaries,
  limitations = limitations,
  path_matrix = path_matrix,
  artifacts = lapply(tarballs, digest_file)
)

json_path <- file.path(output_dir, sprintf("numerical-trust-%s.json", scope))
txt_path <- file.path(output_dir, sprintf("numerical-trust-%s.txt", scope))
dossier_path <- file.path(output_dir, sprintf("release-dossier-%s.md", scope))
if (requireNamespace("jsonlite", quietly = TRUE)) {
  jsonlite::write_json(report, json_path, auto_unbox = TRUE, pretty = TRUE,
                       na = "null", dataframe = "rows")
} else {
  dput(report, file = json_path)
}
writeLines(c(
  sprintf("rfugw numerical trust evidence: %s", scope),
  sprintf("channel: %s", channel),
  sprintf("commit: %s", report$commit),
  sprintf("publication: %s", report$publication_status),
  sprintf("exact transport: %s", emd$status),
  sprintf("entropic FGW: %s", fgw$status),
  sprintf("replay: %s", report$seed_contract$replay)
), txt_path)
artifact_lines <- if (length(report$artifacts)) {
  vapply(report$artifacts, function(x) {
    sprintf("- `%s`: SHA-256 `%s`", x$path, x$sha256)
  }, character(1))
} else {
  "- No source artifact was present in the evidence workspace."
}
writeLines(c(
  "# rfugw release dossier",
  "",
  sprintf("- Evidence scope: `%s`", scope),
  sprintf("- Evidence channel: `%s`", channel),
  sprintf("- Commit: `%s`", report$commit),
  sprintf("- Working tree dirty: `%s`", report$working_tree_dirty),
  sprintf("- Exact-commit evidence: `%s`", report$exact_commit_evidence),
  sprintf("- Publication status: `%s`", report$publication_status),
  sprintf("- Generated: `%s`", report$timestamp_utc),
  "",
  "## Verified support in this bundle",
  "",
  paste0("- ", verified_support),
  "",
  "Representative exact-transport and entropic-FGW certificates are recorded in the adjacent JSON evidence ledger.",
  "",
  "## Experimental and deferred boundaries",
  "",
  paste0("- ", experimental_boundaries),
  "",
  "## Limitations",
  "",
  paste0("- ", limitations),
  "",
  "## Artifacts",
  "",
  artifact_lines,
  "",
  sprintf("Replay seeded trust laws with `%s`.", report$seed_contract$replay)
), dossier_path)
cat(sprintf("wrote %s, %s, and %s\n", json_path, txt_path, dossier_path))
