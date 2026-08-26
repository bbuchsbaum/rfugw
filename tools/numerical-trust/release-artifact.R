#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
value_arg <- function(prefix, default = NULL) {
  hit <- args[startsWith(args, prefix)]
  if (!length(hit)) return(default)
  sub(prefix, "", hit[[1L]], fixed = TRUE)
}

mode <- value_arg("--mode=", "verify")
manifest_path <- value_arg("--manifest=", "release-artifact.json")
artifact_path <- value_arg("--artifact=", "auto")
expected_commit <- value_arg("--commit=", Sys.getenv("GITHUB_SHA", ""))

if (!requireNamespace("digest", quietly = TRUE)) {
  stop("Release-artifact verification requires the `digest` package.", call. = FALSE)
}
if (!requireNamespace("jsonlite", quietly = TRUE)) {
  stop("Release-artifact verification requires the `jsonlite` package.", call. = FALSE)
}

find_artifact <- function(path, directory = ".") {
  if (!identical(path, "auto")) return(path)
  candidates <- list.files(
    directory, pattern = "^rfugw_[^/]+\\.tar\\.gz$", full.names = TRUE
  )
  if (length(candidates) != 1L) {
    stop(
      sprintf(
        "Expected exactly one rfugw source artifact in `%s`; found %d.",
        directory, length(candidates)
      ),
      call. = FALSE
    )
  }
  candidates[[1L]]
}

sha256_file <- function(path) {
  if (!file.exists(path)) {
    stop("Release artifact does not exist: ", path, call. = FALSE)
  }
  unname(digest::digest(path, algo = "sha256", file = TRUE))
}

if (identical(mode, "manifest")) {
  artifact_path <- find_artifact(artifact_path)
  if (!nzchar(expected_commit)) {
    stop("Manifest creation requires an exact source commit.", call. = FALSE)
  }
  filename <- basename(artifact_path)
  version <- sub("^rfugw_(.+)\\.tar\\.gz$", "\\1", filename)
  if (identical(version, filename)) {
    stop("Artifact name does not encode the rfugw package version.", call. = FALSE)
  }
  manifest <- list(
    schema_version = 1L,
    package = "rfugw",
    version = version,
    artifact = filename,
    sha256 = sha256_file(artifact_path),
    bytes = unname(file.info(artifact_path)$size),
    source_commit = expected_commit,
    created_utc = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
  )
  jsonlite::write_json(
    manifest, manifest_path, auto_unbox = TRUE, pretty = TRUE
  )
  cat(sprintf(
    "manifested %s at commit %s with SHA-256 %s\n",
    filename, expected_commit, manifest$sha256
  ))
} else if (identical(mode, "verify")) {
  if (!file.exists(manifest_path)) {
    stop("Release-artifact manifest does not exist: ", manifest_path,
         call. = FALSE)
  }
  manifest <- jsonlite::read_json(manifest_path, simplifyVector = TRUE)
  required <- c(
    "schema_version", "package", "version", "artifact", "sha256", "bytes",
    "source_commit"
  )
  missing <- setdiff(required, names(manifest))
  if (length(missing)) {
    stop("Release-artifact manifest is missing: ", paste(missing, collapse = ", "),
         call. = FALSE)
  }
  if (!identical(as.integer(manifest$schema_version), 1L) ||
      !identical(manifest$package, "rfugw")) {
    stop("Unsupported release-artifact manifest identity.", call. = FALSE)
  }
  manifest_dir <- dirname(normalizePath(manifest_path, mustWork = TRUE))
  artifact_path <- if (identical(artifact_path, "auto")) {
    file.path(manifest_dir, manifest$artifact)
  } else {
    artifact_path
  }
  if (!identical(basename(artifact_path), manifest$artifact)) {
    stop("Release artifact name does not match its manifest.", call. = FALSE)
  }
  observed_sha <- sha256_file(artifact_path)
  if (!identical(observed_sha, manifest$sha256)) {
    stop(
      sprintf(
        "Release artifact SHA-256 mismatch: expected %s, observed %s.",
        manifest$sha256, observed_sha
      ),
      call. = FALSE
    )
  }
  observed_bytes <- unname(file.info(artifact_path)$size)
  if (!identical(as.numeric(observed_bytes), as.numeric(manifest$bytes))) {
    stop("Release artifact byte size does not match its manifest.", call. = FALSE)
  }
  if (nzchar(expected_commit) &&
      !identical(manifest$source_commit, expected_commit)) {
    stop(
      sprintf(
        "Release artifact commit mismatch: expected %s, observed %s.",
        expected_commit, manifest$source_commit
      ),
      call. = FALSE
    )
  }
  cat(sprintf(
    "verified %s at commit %s with SHA-256 %s\n",
    manifest$artifact, manifest$source_commit, observed_sha
  ))
} else {
  stop("`--mode` must be `manifest` or `verify`.", call. = FALSE)
}
