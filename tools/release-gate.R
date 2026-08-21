#!/usr/bin/env Rscript
# Local release-quality gate for rfugw.
#
# Builds the clean tarball, runs source consistency checks, R CMD check,
# and a skip inventory. Writes machine-readable evidence under .gate/.
# This is local evidence only. Hosted evidence is GitHub Actions.
# Publication evidence is a separate CRAN / R-universe record.

args <- commandArgs(trailingOnly = TRUE)
as_cran <- !"--no-as-cran" %in% args
channel <- Sys.getenv("RFUGW_EVIDENCE_CHANNEL", "local")

if (!file.exists("DESCRIPTION")) {
  stop("Run tools/release-gate.R from the package root.", call. = FALSE)
}

desc <- read.dcf("DESCRIPTION")
pkg <- unname(desc[1, "Package"])
version <- unname(desc[1, "Version"])
commit <- tryCatch(
  system2("git", c("rev-parse", "HEAD"), stdout = TRUE, stderr = FALSE),
  error = function(e) NA_character_
)
if (length(commit) != 1L) commit <- NA_character_
dirty_lines <- tryCatch(
  system2("git", c("status", "--porcelain"), stdout = TRUE, stderr = FALSE),
  error = function(e) character()
)
working_tree_dirty <- length(dirty_lines) > 0L
source_identity <- if (working_tree_dirty) paste0(commit, "+dirty") else commit

gate_dir <- file.path(".gate")
dir.create(gate_dir, showWarnings = FALSE)
stamp <- format(Sys.time(), "%Y%m%dT%H%M%SZ", tz = "UTC")

reviewed_skip_patterns <- c(
  "\\{lpSolve\\} is not installed",
  "requireNamespace\\(\"lpSolve\""
)

suggests <- trimws(unlist(strsplit(gsub("\\n", " ", desc[1, "Suggests"]), ",")))
suggests <- sub("\\s*\\(.*$", "", suggests)
suggests <- suggests[nzchar(suggests)]
# Complete --as-cran evidence needs every Suggests package. Bench-only
# extras can be missing only when --allow-missing-suggests is set.
bench_only <- c("bench")
allow_missing <- "--allow-missing-suggests" %in% args
required_suggests <- if (allow_missing) setdiff(suggests, bench_only) else suggests
missing_required <- required_suggests[!vapply(required_suggests, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_required)) {
  repos <- getOption("repos")
  if (is.null(repos) || identical(unname(repos)[1], "@CRAN@")) {
    options(repos = c(CRAN = "https://cloud.r-project.org"))
  }
  try(install.packages(missing_required, quiet = TRUE), silent = TRUE)
  missing_required <- required_suggests[!vapply(required_suggests, requireNamespace, logical(1), quietly = TRUE)]
}
missing_optional <- setdiff(suggests, required_suggests)
missing_optional <- missing_optional[!vapply(missing_optional, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_required) && !allow_missing) {
  stop(
    "Release gate needs these Suggests packages installed: ",
    paste(missing_required, collapse = ", "),
    call. = FALSE
  )
}

gate_fail <- function(...) {
  list(ok = FALSE, detail = paste0(...))
}
gate_ok <- function(detail = "") {
  list(ok = TRUE, detail = detail)
}

run_script <- function(rel) {
  status <- system2("Rscript", rel)
  if (!identical(status, 0L)) {
    return(gate_fail(rel, " exited with status ", status))
  }
  gate_ok(rel)
}

generated <- run_script("tools/check-generated.R")
docs <- run_script("tools/check-docs.R")
old_scope <- Sys.getenv("RFUGW_TRUST_SCOPE", unset = NA_character_)
Sys.setenv(RFUGW_TRUST_SCOPE = "release")
mutation <- run_script("tools/numerical-trust/run-mutation-proof.R")

tarball <- sprintf("%s_%s.tar.gz", pkg, version)
if (file.exists(tarball)) {
  unlink(tarball)
}
build_status <- system2("R", c("CMD", "build", ".", "--md5"))
if (!identical(build_status, 0L) || !file.exists(tarball)) {
  build <- gate_fail("R CMD build failed")
  digest <- NA_character_
} else {
  digest <- unname(tools::md5sum(tarball))
  sha256 <- if (nzchar(Sys.which("shasum"))) {
    sub("\\s.*$", "", system2("shasum", c("-a", "256", tarball), stdout = TRUE)[1])
  } else if (nzchar(Sys.which("sha256sum"))) {
    sub("\\s.*$", "", system2("sha256sum", tarball, stdout = TRUE)[1])
  } else {
    digest
  }
  build <- gate_ok(sprintf("%s md5=%s sha256=%s", tarball, digest, sha256))
}

check_args <- c("CMD", "check")
if (as_cran) check_args <- c(check_args, "--as-cran")
check_args <- c(check_args, "--no-manual", tarball)
check_status <- if (isTRUE(build$ok)) system2("R", check_args) else 1L
check_dir <- paste0(pkg, ".Rcheck")
log_00 <- file.path(check_dir, "00check.log")
errors <- character()
warnings <- character()
notes <- character()
if (file.exists(log_00)) {
  lines <- readLines(log_00, warn = FALSE)
  errors <- grep("^ERROR", lines, value = TRUE)
  warnings <- grep("^WARNING", lines, value = TRUE)
  notes <- grep("^NOTE", lines, value = TRUE)
  # rcmdcheck-style: also capture "* checking ... ERROR" blocks
  err_idx <- grep("ERROR$", lines)
  warn_idx <- grep("WARNING$", lines)
  note_idx <- grep("NOTE$", lines)
  take_block <- function(idx) {
    vapply(idx, function(i) {
      end <- min(i + 8L, length(lines))
      paste(lines[i:end], collapse = "\n")
    }, character(1))
  }
  if (!length(errors) && length(err_idx)) errors <- take_block(err_idx)
  if (!length(warnings) && length(warn_idx)) warnings <- take_block(warn_idx)
  if (!length(notes) && length(note_idx)) notes <- take_block(note_idx)
}
check <- if (identical(check_status, 0L) && !length(errors) && !length(warnings)) {
  gate_ok(sprintf("R CMD check status=%s", check_status))
} else {
  gate_fail(sprintf(
    "R CMD check status=%s errors=%s warnings=%s",
    check_status, length(errors), length(warnings)
  ))
}

skips <- tryCatch(
  {
    res <- testthat::test_local(".", reporter = "silent", stop_on_failure = FALSE)
    msgs <- unlist(lapply(res, function(file_res) {
      vapply(file_res$results, function(r) {
        if (inherits(r, "expectation_skip")) as.character(r$message) else NA_character_
      }, character(1))
    }), use.names = FALSE)
    unique(msgs[!is.na(msgs)])
  },
  error = function(e) {
    test_log <- file.path(check_dir, "tests", "testthat.Rout")
    if (file.exists(test_log)) {
      unique(trimws(grep("Reason:", readLines(test_log, warn = FALSE), value = TRUE)))
    } else {
      character()
    }
  }
)
if (is.na(old_scope)) Sys.unsetenv("RFUGW_TRUST_SCOPE") else
  Sys.setenv(RFUGW_TRUST_SCOPE = old_scope)
unreviewed <- skips
for (pat in reviewed_skip_patterns) {
  unreviewed <- unreviewed[!grepl(pat, unreviewed)]
}
skips_ok <- !length(unreviewed)

status <- if (
  isTRUE(generated$ok) &&
    isTRUE(docs$ok) &&
    isTRUE(mutation$ok) &&
    isTRUE(build$ok) &&
    isTRUE(check$ok) &&
    skips_ok &&
    !length(missing_required)
) {
  "pass"
} else {
  "fail"
}

report <- list(
  package = pkg,
  version = version,
  commit = commit,
  source_identity = source_identity,
  working_tree_dirty = working_tree_dirty,
  build_environment = as.list(Sys.getenv(c(
    "RFUGW_FAST_FLAGS", "RFUGW_EXTRA_CXXFLAGS", "RFUGW_EXTRA_LIBS",
    "RFUGW_OPENMP_FLAGS", "RFUGW_OPENMP_LIBS", "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS",
    "BLIS_NUM_THREADS"
  ), unset = NA_character_)),
  channel = channel,
  timestamp = stamp,
  tarball = if (file.exists(tarball)) tarball else NA_character_,
  md5 = if (exists("digest")) digest else NA_character_,
  sha256 = if (exists("sha256")) sha256 else NA_character_,
  errors = length(errors),
  warnings = length(warnings),
  notes = length(notes),
  error_text = errors,
  warning_text = warnings,
  note_text = notes,
  skips = skips,
  reviewed_skip_patterns = reviewed_skip_patterns,
  unreviewed_skips = unreviewed,
  missing_required_suggests = missing_required,
  missing_optional_suggests = missing_optional,
  generated = generated,
  docs = docs,
  mutation = mutation,
  runtime_provenance = tryCatch({
    M_provenance <- matrix(c(0, 1, 1, 0), 2, 2)
    fgw_entropic(
      M_provenance, M_provenance, M_provenance,
      epsilon = 0.1, max_iter = 20L, tol = 1e-7,
      sinkhorn_tol = 1e-8, sinkhorn_method = "auto",
      precision = "strict_double"
    )$runtime_provenance
  }, error = function(e) list(error = conditionMessage(e))),
  build = build,
  check = check,
  status = status
)

json_path <- file.path(gate_dir, "latest.json")
txt_path <- file.path(gate_dir, "latest.txt")
hist_json <- file.path(gate_dir, paste0(stamp, ".json"))
if (requireNamespace("jsonlite", quietly = TRUE)) {
  jsonlite::write_json(report, json_path, auto_unbox = TRUE, pretty = TRUE)
  file.copy(json_path, hist_json, overwrite = TRUE)
} else {
  dput(report, file = json_path)
}

txt <- c(
  sprintf("rfugw release-quality gate: %s", status),
  sprintf("package: %s %s", pkg, version),
  sprintf("commit:  %s", commit),
  sprintf("source:  %s", source_identity),
  sprintf("dirty:   %s", working_tree_dirty),
  sprintf("channel: %s", channel),
  sprintf("tarball: %s", report$tarball),
  sprintf("md5:     %s", report$md5),
  sprintf("sha256:  %s", report$sha256),
  sprintf("errors/warnings/notes: %s / %s / %s", length(errors), length(warnings), length(notes)),
  sprintf("skips: %s", if (length(skips)) paste(skips, collapse = " | ") else "none"),
  sprintf("unreviewed skips: %s", if (length(unreviewed)) paste(unreviewed, collapse = " | ") else "none"),
  sprintf("generated: %s", generated$detail),
  sprintf("docs: %s", docs$detail),
  sprintf("mutation: %s", mutation$detail),
  sprintf("build: %s", build$detail),
  sprintf("check: %s", check$detail)
)
writeLines(txt, txt_path)
writeLines(txt, file.path(gate_dir, paste0(stamp, ".txt")))
cat(paste(txt, collapse = "\n"), "\n", sep = "")

if (!identical(status, "pass")) {
  quit(status = 1L)
}
