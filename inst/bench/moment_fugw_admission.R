.mfa_or <- function(value, default) if (is.null(value)) default else value
`%||%` <- .mfa_or

.mfa_source_dir <- local({
  source_files <- unlist(lapply(sys.frames(), function(frame) {
    value <- frame$ofile
    if (is.null(value) || !length(value)) character() else as.character(value)
  }), use.names = FALSE)
  source_file <- if (length(source_files)) tail(source_files, 1L) else NULL
  if (is.null(source_file) || !nzchar(source_file)) return(NULL)
  dirname(normalizePath(source_file, mustWork = TRUE))
})

moment_admission_resource <- function(name) {
  if (!is.null(.mfa_source_dir)) {
    sibling <- file.path(.mfa_source_dir, name)
    if (file.exists(sibling)) return(sibling)
  }
  search_root <- normalizePath(".", mustWork = TRUE)
  for (depth in 0:4) {
    source_path <- file.path(search_root, "inst", "bench", name)
    description <- file.path(search_root, "DESCRIPTION")
    if (file.exists(source_path) && file.exists(description) &&
        any(grepl("^Package:[[:space:]]+rfugw$", readLines(
          description, warn = FALSE
        )))) {
      return(source_path)
    }
    parent <- dirname(search_root)
    if (identical(parent, search_root)) break
    search_root <- parent
  }
  installed <- system.file("bench", name, package = "rfugw")
  if (nzchar(installed) && file.exists(installed)) return(installed)
  stop("Moment-FUGW admission resource not found: ", name, call. = FALSE)
}

moment_admission_repository_file <- function(relative_path) {
  if (!is.character(relative_path) || length(relative_path) != 1L ||
      is.na(relative_path) || !nzchar(relative_path) ||
      grepl("^([A-Za-z]:)?[/\\\\]", relative_path) ||
      any(strsplit(gsub("\\\\", "/", relative_path), "/", fixed = TRUE)[[1L]] ==
          "..")) {
    stop("Invalid repository-relative admission path.", call. = FALSE)
  }
  roots <- character()
  if (!is.null(.mfa_source_dir)) {
    package_parent <- file.path(.mfa_source_dir, "..", "..")
    roots <- c(
      roots,
      package_parent,
      file.path(package_parent, "00_pkg_src", "rfugw")
    )
  }
  search_root <- normalizePath(".", mustWork = TRUE)
  for (depth in 0:4) {
    roots <- c(
      roots,
      search_root,
      file.path(search_root, "00_pkg_src", "rfugw")
    )
    parent <- dirname(search_root)
    if (identical(parent, search_root)) break
    search_root <- parent
  }
  roots <- unique(roots)
  for (root in roots) {
    candidate <- file.path(root, relative_path)
    description <- file.path(root, "DESCRIPTION")
    if (file.exists(candidate) && file.exists(description) &&
        any(grepl("^Package:[[:space:]]+rfugw$", readLines(
          description, warn = FALSE
        )))) {
      return(normalizePath(candidate, mustWork = TRUE))
    }
  }
  stop("Moment-FUGW repository file not found: ", relative_path,
       call. = FALSE)
}

.mfa_require_packages <- function() {
  for (package in c("digest", "jsonlite")) {
    if (!requireNamespace(package, quietly = TRUE)) {
      stop("Moment-FUGW admission verification requires `", package, "`.",
           call. = FALSE)
    }
  }
}

.mfa_valid_commit <- function(value) {
  is.character(value) && length(value) == 1L && !is.na(value) &&
    grepl("^[[:xdigit:]]{40}$", value)
}

.mfa_character <- function(value, default = NA_character_) {
  if (is.null(value) || !length(value)) return(default)
  as.character(value[[1L]])
}

.mfa_numeric <- function(value, default = NA_real_) {
  if (is.null(value) || !length(value)) return(default)
  suppressWarnings(as.numeric(value[[1L]]))
}

.mfa_integer <- function(value, default = NA_integer_) {
  if (is.null(value) || !length(value)) return(default)
  suppressWarnings(as.integer(value[[1L]]))
}

.mfa_same_numeric <- function(left, right, tolerance = 1e-12) {
  left <- .mfa_numeric(left)
  right <- .mfa_numeric(right)
  is.finite(left) && is.finite(right) &&
    abs(left - right) <= tolerance * max(1, abs(right))
}

.mfa_sha256 <- function(path) {
  if (!file.exists(path) || dir.exists(path)) return(NA_character_)
  unname(digest::digest(path, algo = "sha256", file = TRUE))
}

.mfa_read_json <- function(path) {
  if (!file.exists(path)) stop("Missing JSON artifact: ", path, call. = FALSE)
  jsonlite::fromJSON(path, simplifyVector = TRUE)
}

.mfa_tree_fingerprint <- function(root, files = NULL) {
  root <- normalizePath(root, mustWork = TRUE)
  if (is.null(files)) {
    files <- list.files(root, recursive = TRUE, full.names = TRUE,
                        all.files = FALSE, no.. = TRUE)
  } else {
    files <- file.path(root, files)
  }
  files <- files[file.exists(files) & !dir.exists(files)]
  files <- sort(unique(normalizePath(files, mustWork = TRUE)))
  relative <- substring(files, nchar(root) + 2L)
  hashes <- vapply(files, .mfa_sha256, character(1))
  bytes <- as.numeric(file.info(files)$size)
  canonical <- paste(relative, hashes, format(bytes, scientific = FALSE),
                     sep = "\t", collapse = "\n")
  list(
    sha256 = digest::digest(canonical, algo = "sha256", serialize = FALSE),
    files = data.frame(
      path = relative, sha256 = hashes, bytes = bytes,
      stringsAsFactors = FALSE
    )
  )
}

.mfa_platform <- function(sysname) {
  value <- .mfa_character(sysname)
  if (value %in% c("Darwin", "Linux", "Windows")) value else NA_character_
}

.mfa_expected_counts <- function(contract) {
  profile <- contract$profiles$full
  scaling <- contract$scaling
  fixed <- length(profile$domain_sizes) *
    (as.integer(profile$warmup) + as.integer(profile$repetitions))
  allocations <- length(profile$domain_sizes)
  memory_actions <- length(profile$domain_sizes) * 3L
  ranks <- length(scaling$structure_ranks)
  converged <- length(profile$converged_sizes) *
    (as.integer(profile$warmup) + as.integer(profile$repetitions))
  operators <- length(scaling$heldout_batches) * 2L
  workers <- 1L + fixed + allocations + memory_actions + ranks +
    converged + operators
  list(workers = workers, rows = workers - 1L, raw_files = workers * 4L)
}

.mfa_loglog_exponent <- function(size, value) {
  keep <- is.finite(size) & size > 0 & is.finite(value) & value > 0
  size <- as.numeric(size[keep])
  value <- as.numeric(value[keep])
  if (length(unique(size)) < 2L) return(NA_real_)
  unname(stats::coef(stats::lm(log2(value) ~ log2(size)))[[2L]])
}

.mfa_curve_medians <- function(runs, curve, action, value) {
  keep <- runs$curve == curve & runs$action == action &
    runs$phase == "measured" & runs$valid %in% TRUE &
    is.finite(runs[[value]])
  selected <- runs[keep, , drop = FALSE]
  if (!nrow(selected)) {
    return(data.frame(n_source = integer(), value = numeric()))
  }
  output <- stats::aggregate(
    selected[[value]], list(n_source = selected$n_source), stats::median,
    na.rm = TRUE
  )
  names(output) <- c("n_source", "value")
  output[order(output$n_source), , drop = FALSE]
}

.mfa_design_keys <- function(contract) {
  profile <- contract$profiles$full
  fixed_sizes <- as.integer(profile$domain_sizes)
  converged_sizes <- as.integer(profile$converged_sizes)
  ranks <- as.integer(contract$scaling$structure_ranks)
  batches <- as.integer(contract$scaling$heldout_batches)
  key <- function(curve, action, n, phase, repetition, rank, batch) {
    paste(curve, action, n, phase, repetition, rank, batch, sep = "|")
  }
  repeated_solve <- function(curve, sizes) {
    unlist(lapply(sizes, function(n) c(
      if (profile$warmup > 0L) {
        vapply(seq_len(profile$warmup), function(repetition) {
          key(curve, "solve", n, "warmup", repetition,
              contract$problem$structure_rank, NA_integer_)
        }, character(1))
      },
      vapply(seq_len(profile$repetitions), function(repetition) {
        key(curve, "solve", n, "measured", repetition,
            contract$problem$structure_rank, NA_integer_)
      }, character(1))
    )), use.names = FALSE)
  }
  c(
    repeated_solve("fixed_work_scaling", fixed_sizes),
    vapply(fixed_sizes, function(n) {
      key("r_visible_allocation", "solve", n, "measured", 1L,
          contract$problem$structure_rank, NA_integer_)
    }, character(1)),
    unlist(lapply(fixed_sizes, function(n) {
      vapply(c("objective", "apply", "adjoint"), function(action) {
        key("matrix_free_memory", action, n, "measured", 1L,
            contract$problem$structure_rank, 1L)
      }, character(1))
    }), use.names = FALSE),
    vapply(ranks, function(rank) {
      key("rank_scaling", "solve", profile$rank_size, "measured", 1L,
          rank, NA_integer_)
    }, character(1)),
    repeated_solve("end_to_end_convergence", converged_sizes),
    unlist(lapply(batches, function(batch) {
      vapply(c("apply", "adjoint"), function(action) {
        key("operator_throughput", action, profile$operator_size,
            "measured", 1L, contract$problem$structure_rank, batch)
      }, character(1))
    }), use.names = FALSE)
  )
}

.mfa_validate_source <- function(manifest_path, expected_commit = NULL) {
  reasons <- character()
  manifest <- tryCatch(.mfa_read_json(manifest_path), error = function(error) {
    reasons <<- c(reasons, "source_manifest_missing_or_invalid")
    list()
  })
  if (!identical(.mfa_integer(manifest$schema_version), 1L) ||
      !identical(.mfa_character(manifest$package, ""), "rfugw")) {
    reasons <- c(reasons, "source_manifest_identity_invalid")
  }
  if (!isTRUE(manifest$source_tree_clean) ||
      !isTRUE(manifest$git_provenance_complete) ||
      !isTRUE(manifest$exact_commit_evidence) ||
      !identical(.mfa_integer(manifest$git_status_entry_count), 0L) ||
      !identical(.mfa_integer(manifest$git_commit_command_status), 0L) ||
      !identical(.mfa_integer(manifest$git_status_command_status), 0L)) {
    reasons <- c(reasons, "source_manifest_not_exact_clean")
  }
  version <- .mfa_character(manifest$version, "")
  if (!grepl("^[0-9]+(?:[.-][0-9A-Za-z]+)+$", version)) {
    reasons <- c(reasons, "source_manifest_version_invalid")
  }
  commit <- .mfa_character(manifest$source_commit)
  if (!.mfa_valid_commit(commit)) {
    reasons <- c(reasons, "source_manifest_commit_invalid")
  } else {
    commit <- tolower(commit)
  }
  if (!is.null(expected_commit) &&
      (!.mfa_valid_commit(expected_commit) ||
       !identical(commit, tolower(expected_commit)))) {
    reasons <- c(reasons, "source_manifest_commit_mismatch")
  }
  artifact_name <- .mfa_character(manifest$artifact, "")
  if (!nzchar(artifact_name) ||
      !identical(basename(artifact_name), artifact_name) ||
      !identical(
        artifact_name, paste0("rfugw_", version, ".tar.gz")
      )) {
    reasons <- c(reasons, "source_artifact_name_invalid")
  }
  artifact_path <- file.path(dirname(manifest_path), artifact_name)
  observed_sha <- .mfa_sha256(artifact_path)
  expected_sha <- tolower(.mfa_character(manifest$sha256, ""))
  if (!grepl("^[[:xdigit:]]{64}$", expected_sha) ||
      is.na(observed_sha) || !identical(observed_sha, expected_sha)) {
    reasons <- c(reasons, "source_artifact_sha256_mismatch")
  }
  observed_bytes <- if (file.exists(artifact_path) && !dir.exists(artifact_path)) {
    as.numeric(file.info(artifact_path)$size)
  } else NA_real_
  if (!is.finite(observed_bytes) ||
      !identical(observed_bytes, .mfa_numeric(manifest$bytes))) {
    reasons <- c(reasons, "source_artifact_size_mismatch")
  }
  list(
    passed = !length(unique(reasons)), reasons = unique(reasons),
    commit = commit, version = version, artifact = artifact_path,
    sha256 = observed_sha, bytes = observed_bytes,
    manifest_sha256 = .mfa_sha256(manifest_path)
  )
}

.mfa_validate_performance <- function(
    summary_path, source_commit, source_version, contract) {
  reasons <- character()
  root <- dirname(normalizePath(summary_path, mustWork = TRUE))
  meta_path <- file.path(root, "meta.json")
  runs_path <- file.path(root, "runs.csv")
  raw_dir <- file.path(root, "raw")
  summary <- tryCatch(.mfa_read_json(summary_path), error = function(error) {
    reasons <<- c(reasons, "performance_summary_invalid")
    list()
  })
  meta <- tryCatch(.mfa_read_json(meta_path), error = function(error) {
    reasons <<- c(reasons, "performance_meta_invalid")
    list()
  })
  runs <- tryCatch(
    utils::read.csv(runs_path, stringsAsFactors = FALSE),
    error = function(error) {
      reasons <<- c(reasons, "performance_runs_invalid")
      data.frame()
    }
  )
  platform <- .mfa_platform(meta$sysname)
  if (is.na(platform)) reasons <- c(reasons, "performance_platform_invalid")
  backends <- c(Darwin = "macos_time_l_bytes", Linux = "linux_time_v_kib")
  expected_backend <- if (!is.na(platform) && platform %in% names(backends)) {
    unname(backends[[platform]])
  } else {
    NULL
  }
  if (is.null(expected_backend) ||
      !identical(.mfa_character(meta$rss_backend, ""),
                 expected_backend)) {
    reasons <- c(reasons, "performance_rss_backend_invalid")
  }
  if (!identical(.mfa_integer(summary$schema_version), 1L) ||
      !identical(.mfa_character(summary$profile, ""), "full") ||
      !identical(.mfa_character(summary$package_mode, ""),
                 "installed")) {
    reasons <- c(reasons, "performance_not_full_installed")
  }
  commit <- .mfa_character(meta$commit)
  if (!.mfa_valid_commit(commit) ||
      !identical(tolower(commit), source_commit)) {
    reasons <- c(reasons, "performance_commit_mismatch")
  }
  if (!isTRUE(meta$git_provenance_complete) ||
      !identical(meta$git_dirty, FALSE) ||
      !identical(.mfa_integer(meta$git_status_entry_count), 0L) ||
      !identical(.mfa_integer(meta$git_commit_command_status), 0L) ||
      !identical(.mfa_integer(meta$git_status_command_status), 0L)) {
    reasons <- c(reasons, "performance_source_not_exact_clean")
  }
  if (!identical(.mfa_integer(meta$benchmark_schema_version), 1L) ||
      !identical(.mfa_character(meta$package, ""), "rfugw") ||
      !identical(.mfa_character(meta$version, ""), source_version) ||
      !identical(.mfa_character(meta$profile, ""),
                 "moment_fugw_full_installed") ||
      !identical(.mfa_character(meta$package_mode, ""), "installed")) {
    reasons <- c(reasons, "performance_meta_not_installed")
  }
  if (!identical(.mfa_character(summary$rss_backend, ""), expected_backend)) {
    reasons <- c(reasons, "performance_summary_rss_backend_invalid")
  }
  threshold_md5 <- unname(tools::md5sum(
    moment_admission_resource("moment-fugw-thresholds.json")
  )[[1L]])
  if (!identical(.mfa_character(meta$thresholds_md5, ""),
                 threshold_md5)) {
    reasons <- c(reasons, "performance_threshold_contract_mismatch")
  }
  required_gates <- c(
    "design", "full_profile_complete", "exact_commit_provenance",
    "clean_exact_source_tree", "release_hardware",
    "fixed_work_runtime_envelope", "near_linear_peak_memory",
    "prohibited_dense_allocations_absent", "r_visible_allocation_complete",
    "rank_curve_complete", "operator_batch_curve_complete",
    "larger_certified_endpoint"
  )
  gates <- summary$gates %||% list()
  gate_pass <- vapply(required_gates, function(name) isTRUE(gates[[name]]),
                      logical(1))
  if (!all(gate_pass) || !isTRUE(summary$admission_evaluable)) {
    reasons <- c(reasons, paste0(
      "performance_gate_", required_gates[!gate_pass]
    ))
    if (!isTRUE(summary$admission_evaluable)) {
      reasons <- c(reasons, "performance_not_admission_evaluable")
    }
  }
  full_sizes <- as.integer(contract$scaling$domain_sizes)
  if (!identical(as.integer(summary$completed_fixed_work_sizes), full_sizes) ||
      !identical(as.integer(summary$configured_full_sizes), full_sizes)) {
    reasons <- c(reasons, "performance_size_curve_incomplete")
  }
  runtime <- .mfa_numeric(summary$fixed_work_runtime_exponent)
  memory <- .mfa_numeric(summary$incremental_memory_exponent)
  endpoint <- .mfa_integer(summary$largest_completed_certified_endpoint)
  if (!is.finite(runtime) ||
      runtime < contract$scaling$fixed_work_runtime_exponent_min ||
      runtime > contract$scaling$fixed_work_runtime_exponent_max) {
    reasons <- c(reasons, "performance_runtime_exponent_outside_contract")
  }
  if (!is.finite(memory) ||
      memory >= contract$scaling$memory_doubling_exponent_max) {
    reasons <- c(reasons, "performance_memory_exponent_outside_contract")
  }
  if (!is.finite(endpoint) || endpoint < contract$problem$minimum_certified_endpoint) {
    reasons <- c(reasons, "performance_certified_endpoint_too_small")
  }
  counts <- .mfa_expected_counts(contract)
  if (nrow(runs) != counts$rows) {
    reasons <- c(reasons, "performance_run_count_mismatch")
  }
  required_run_columns <- c(
    "curve", "action", "n_source", "n_target", "phase", "repetition",
    "structure_rank", "batch", "valid", "certified", "timing_eligible",
    "infrastructure_error", "process_status", "reportable_seconds",
    "peak_rss_bytes", "incremental_peak_rss_bytes",
    "dense_sentinel_passed", "full_block_calls", "max_tile_elements",
    "full_matrix_elements", "r_visible_allocation_bytes", "spec_artifact",
    "result_artifact", "time_artifact"
  )
  missing_run_columns <- setdiff(required_run_columns, names(runs))
  if (length(missing_run_columns)) {
    reasons <- c(reasons, paste0(
      "performance_runs_column_missing_", missing_run_columns
    ))
  }
  if (nrow(runs) && !length(missing_run_columns)) {
    observed_design <- paste(
      runs$curve, runs$action, runs$n_source, runs$phase, runs$repetition,
      runs$structure_rank, runs$batch, sep = "|"
    )
    expected_design <- .mfa_design_keys(contract)
    if (!identical(sort(observed_design), sort(expected_design)) ||
        any(runs$n_source != runs$n_target)) {
      reasons <- c(reasons, "performance_design_rows_mismatch")
    }
    if (!all(runs$valid %in% TRUE) ||
        any(runs$infrastructure_error %in% TRUE) ||
        !all(runs$process_status %in% 0L)) {
      reasons <- c(reasons, "performance_invalid_or_infrastructure_row")
    }
    if (!all(runs$dense_sentinel_passed %in% TRUE) ||
        !all(is.finite(runs$full_block_calls)) ||
        !all(runs$full_block_calls == 0L)) {
      reasons <- c(reasons, "performance_dense_allocation_sentinel_failed")
    }
    allocation <- runs[
      runs$curve %in% c("r_visible_allocation", "matrix_free_memory"),
      , drop = FALSE
    ]
    expected_pairs <- as.vector(outer(
      full_sizes, as.character(contract$scaling$r_visible_allocation_actions),
      function(n, action) paste(n, action, sep = ":")
    ))
    observed_pairs <- paste(allocation$n_source, allocation$action, sep = ":")
    if (!identical(sort(expected_pairs), sort(observed_pairs)) ||
        any(!is.finite(allocation$r_visible_allocation_bytes)) ||
        any(allocation$r_visible_allocation_bytes <= 0)) {
      reasons <- c(reasons, "performance_r_visible_allocation_incomplete")
    }
    fixed_time <- .mfa_curve_medians(
      runs, "fixed_work_scaling", "solve", "reportable_seconds"
    )
    fixed_memory <- .mfa_curve_medians(
      runs, "fixed_work_scaling", "solve", "incremental_peak_rss_bytes"
    )
    recomputed_runtime_all <- .mfa_loglog_exponent(
      fixed_time$n_source, fixed_time$value
    )
    fit_count <- min(3L, nrow(fixed_time))
    recomputed_runtime <- if (fit_count >= 2L) {
      .mfa_loglog_exponent(
        tail(fixed_time$n_source, fit_count),
        tail(fixed_time$value, fit_count)
      )
    } else {
      NA_real_
    }
    recomputed_memory <- .mfa_loglog_exponent(
      fixed_memory$n_source, fixed_memory$value
    )
    runtime_all <- .mfa_numeric(
      summary$fixed_work_runtime_exponent_all_points
    )
    slope_tolerance <- 5e-4
    if (!is.finite(recomputed_runtime) ||
        abs(runtime - recomputed_runtime) > slope_tolerance ||
        !is.finite(recomputed_runtime_all) ||
        abs(runtime_all - recomputed_runtime_all) > slope_tolerance ||
        !is.finite(recomputed_memory) ||
        abs(memory - recomputed_memory) > slope_tolerance ||
        !identical(.mfa_character(
          summary$fixed_work_runtime_exponent_fit, ""
        ), contract$scaling$runtime_exponent_primary_fit)) {
      reasons <- c(reasons, "performance_reported_slopes_not_reproducible")
    }
    if (isTRUE(contract$scaling$runtime_exponent_all_points_required) &&
        !is.finite(runtime_all)) {
      reasons <- c(reasons, "performance_all_points_runtime_missing")
    }
    endpoint_rows <- runs[
      runs$curve == "end_to_end_convergence" &
        runs$phase == "measured" & runs$certified %in% TRUE &
        runs$timing_eligible %in% TRUE,
      , drop = FALSE
    ]
    recomputed_endpoint <- if (nrow(endpoint_rows)) {
      max(endpoint_rows$n_source)
    } else {
      NA_integer_
    }
    if (!identical(endpoint, as.integer(recomputed_endpoint)) ||
        !identical(.mfa_integer(summary$claim_ceiling_n), endpoint)) {
      reasons <- c(reasons, "performance_certified_endpoint_not_reproducible")
    }
    if (any(!is.finite(runs$max_tile_elements)) ||
        any(!is.finite(runs$full_matrix_elements)) ||
        any(runs$max_tile_elements > runs$full_matrix_elements) ||
        any(runs$max_tile_elements >
              pmin(runs$n_source, contract$profiles$full$block_size) *
                runs$n_target)) {
      reasons <- c(reasons, "performance_tile_sentinel_invalid")
    }
  }
  actual_raw <- if (dir.exists(raw_dir)) {
    sort(list.files(raw_dir, recursive = TRUE))
  } else character()
  declared_raw <- sort(as.character(meta$raw_artifacts %||% character()))
  if (!identical(actual_raw, declared_raw) ||
      length(actual_raw) != counts$raw_files) {
    reasons <- c(reasons, "performance_raw_artifact_manifest_mismatch")
  }
  suffix_counts <- vapply(
    c(".spec.json", ".result.json", ".stdout.log", ".time.log"),
    function(suffix) sum(endsWith(actual_raw, suffix)), integer(1)
  )
  if (any(suffix_counts != counts$workers)) {
    reasons <- c(reasons, "performance_raw_worker_receipts_incomplete")
  }
  raw_stems <- lapply(
    c(".spec.json", ".result.json", ".stdout.log", ".time.log"),
    function(suffix) {
      selected <- actual_raw[endsWith(actual_raw, suffix)]
      substr(selected, 1L, nchar(selected) - nchar(suffix))
    }
  )
  if (length(raw_stems) &&
      !all(vapply(raw_stems[-1L], identical, logical(1), raw_stems[[1L]]))) {
    reasons <- c(reasons, "performance_raw_worker_stems_mismatch")
  }
  json_receipts <- file.path(
    raw_dir,
    actual_raw[endsWith(actual_raw, ".spec.json") |
                 endsWith(actual_raw, ".result.json")]
  )
  json_valid <- vapply(json_receipts, function(path) {
    !inherits(try(jsonlite::fromJSON(path), silent = TRUE), "try-error")
  }, logical(1))
  if (length(json_valid) != counts$workers * 2L || !all(json_valid)) {
    reasons <- c(reasons, "performance_raw_json_receipt_invalid")
  }
  if (nrow(runs) && !length(missing_run_columns)) {
    run_stems <- sub(
      "[.]spec[.]json$", "", basename(as.character(runs$spec_artifact))
    )
    result_stems <- sub(
      "[.]result[.]json$", "", basename(as.character(runs$result_artifact))
    )
    time_stems <- sub(
      "[.]time[.]log$", "", basename(as.character(runs$time_artifact))
    )
    if (length(unique(run_stems)) != nrow(runs) ||
        !identical(run_stems, result_stems) ||
        !identical(run_stems, time_stems) ||
        !all(run_stems %in% raw_stems[[1L]]) ||
        length(setdiff(raw_stems[[1L]], run_stems)) != 1L) {
      reasons <- c(reasons, "performance_run_raw_receipt_binding_invalid")
    }
  }
  fingerprint <- .mfa_tree_fingerprint(
    root, c("summary.json", "meta.json", "runs.csv",
            file.path("raw", actual_raw))
  )
  list(
    passed = !length(unique(reasons)), reasons = unique(reasons),
    platform = platform, root = root, commit = commit,
    fingerprint = fingerprint$sha256, file_count = nrow(fingerprint$files),
    rows = nrow(runs), raw_files = length(actual_raw),
    runtime_exponent = runtime, memory_exponent = memory,
    certified_endpoint = endpoint
  )
}

.mfa_validate_science <- function(summary_path, source_commit, source_version) {
  reasons <- character()
  root <- dirname(normalizePath(summary_path, mustWork = TRUE))
  summary <- tryCatch(.mfa_read_json(summary_path), error = function(error) {
    reasons <<- c(reasons, "scientific_summary_invalid")
    list()
  })
  platform <- .mfa_platform(summary$environment$sys_info$sysname)
  if (is.na(platform)) reasons <- c(reasons, "scientific_platform_invalid")
  if (!identical(.mfa_integer(summary$schema_version), 1L) ||
      !identical(.mfa_character(summary$status, ""), "passed") ||
      !isTRUE(summary$admission_evaluable) ||
      !identical(.mfa_character(summary$split, ""), "evaluation") ||
      !identical(.mfa_character(summary$package_mode, ""),
                 "installed")) {
    reasons <- c(reasons, "scientific_not_admissible_installed_evaluation")
  }
  candidate <- summary$candidate %||% list()
  commit <- .mfa_character(candidate$git_commit)
  if (!.mfa_valid_commit(commit) ||
      !identical(tolower(commit), source_commit)) {
    reasons <- c(reasons, "scientific_commit_mismatch")
  }
  if (!isTRUE(candidate$git_provenance_complete) ||
      !identical(candidate$git_dirty, FALSE) ||
      !identical(.mfa_integer(candidate$git_status_entry_count), 0L) ||
      !identical(.mfa_integer(candidate$git_commit_command_status), 0L) ||
      !identical(.mfa_integer(candidate$git_status_command_status), 0L)) {
    reasons <- c(reasons, "scientific_source_not_exact_clean")
  }
  if (!identical(.mfa_character(candidate$package_version, ""),
                 source_version)) {
    reasons <- c(reasons, "scientific_package_version_mismatch")
  }
  gate_values <- c(
    image = isTRUE(summary$image$gate),
    surface = isTRUE(summary$surface$gate),
    geometry_rank = isTRUE(summary$geometry_rank_curve$gate)
  )
  if (!all(gate_values)) {
    reasons <- c(reasons, paste0("scientific_gate_", names(gate_values)[
      !gate_values
    ]))
  }
  protocol_path <- moment_admission_resource(
    "moment-fugw-validation-protocol.json"
  )
  protocol <- .mfa_read_json(protocol_path)
  expected_image_rows <- length(protocol$image$evaluation_seeds) *
    length(protocol$image$overlap_percent) * 2L
  expected_surface_rows <- length(protocol$surface$evaluation_seeds)
  expected_rank_rows <- length(protocol$surface$embedding_curve$ranks)
  requested_image <- .mfa_integer(summary$image$requested_rows)
  completed_image <- .mfa_integer(summary$image$completed_rows)
  failed_image <- .mfa_integer(summary$image$failed_rows)
  requested_surface <- .mfa_integer(summary$surface$requested_rows)
  completed_surface <- .mfa_integer(summary$surface$completed_rows)
  failed_surface <- .mfa_integer(summary$surface$failed_rows)
  completed_rank <- .mfa_integer(summary$geometry_rank_curve$completed_rows)
  if (!identical(requested_image, expected_image_rows) ||
      !identical(completed_image, expected_image_rows) ||
      !identical(failed_image, 0L) ||
      !identical(requested_surface, expected_surface_rows) ||
      !identical(completed_surface, expected_surface_rows) ||
      !identical(failed_surface, 0L) ||
      !identical(completed_rank, expected_rank_rows)) {
    reasons <- c(reasons, "scientific_requested_rows_incomplete")
  }
  if (!identical(.mfa_integer(summary$protocol$schema_version),
                 .mfa_integer(protocol$schema_version)) ||
      !identical(.mfa_character(summary$protocol$status, ""),
                 "frozen_before_evaluation") ||
      !identical(.mfa_character(summary$protocol$md5, ""),
                 unname(tools::md5sum(protocol_path)[[1L]]))) {
    reasons <- c(reasons, "scientific_protocol_mismatch")
  }
  artifacts <- c(
    .mfa_character(summary$image$artifact, ""),
    .mfa_character(summary$surface$artifact, ""),
    .mfa_character(summary$geometry_rank_curve$artifact, "")
  )
  expected_artifacts <- c(
    "moment-fugw-image-evaluation.csv",
    "moment-fugw-surface-evaluation.csv",
    "moment-fugw-geometry-rank-evaluation.csv"
  )
  artifacts_available <- !any(!nzchar(artifacts)) &&
    !any(basename(artifacts) != artifacts) &&
    all(file.exists(file.path(root, artifacts)))
  if (!identical(artifacts, expected_artifacts) || !artifacts_available) {
    reasons <- c(reasons, "scientific_artifact_missing")
  }
  tables <- lapply(file.path(root, artifacts), function(path) {
    tryCatch(
      utils::read.csv(path, stringsAsFactors = FALSE),
      error = function(error) data.frame()
    )
  })
  image <- tables[[1L]]
  surface <- tables[[2L]]
  rank_curve <- tables[[3L]]
  image_metrics <- c(
    "transported_mass", "discarded_source_mass", "discarded_target_mass",
    "false_occluded_source_mass", "weighted_endpoint_error",
    "heldout_rmse", "heldout_correlation", "local_distortion",
    "coupling_l1_discrepancy", "geometry_relative_error"
  )
  image_columns <- c(
    "method", "solver_success", "suite_gate", "seed", "overlap_percent",
    "evaluation_split", "protocol_schema_version", "protocol_status",
    "recovered_overlap_fraction", image_metrics
  )
  surface_metrics <- c(
    "geometry_relative_error", "geometry_heldout_error",
    "heldout_correlation", "heldout_r2",
    "expected_geodesic_endpoint_error", "transported_mass",
    "discarded_source_mass", "discarded_target_mass",
    "coupling_l1_discrepancy"
  )
  surface_columns <- c(
    "solver_success", "suite_gate", "evaluation_split",
    "dense_parity_passed", "dense_parity_maximum_error",
    "dense_parity_tolerance", "protocol_schema_version", "protocol_status",
    surface_metrics
  )
  rank_columns <- c(
    "solver_success", "numerical_and_geometry_error_separate",
    "embedding_rank", "protocol_schema_version", "protocol_status",
    surface_metrics
  )
  if (!nrow(image) || !all(image_columns %in% names(image)) ||
      nrow(image) != expected_image_rows ||
      !all(image$solver_success %in% TRUE) ||
      !all(image$suite_gate %in% TRUE) ||
      !all(image$evaluation_split == "frozen_evaluation") ||
      !all(image$protocol_schema_version == protocol$schema_version) ||
      !all(image$protocol_status == "frozen_before_evaluation") ||
      any(!is.finite(as.matrix(image[intersect(image_metrics, names(image))])))) {
    reasons <- c(reasons, "scientific_image_csv_invalid")
  } else {
    expected_method_rows <- length(protocol$image$evaluation_seeds) *
      length(protocol$image$overlap_percent)
    method_counts <- table(image$method)
    if (!identical(
        as.integer(method_counts[c("balanced_fgw", "unbalanced_fugw")]),
        rep(expected_method_rows, 2L)
      )) {
      reasons <- c(reasons, "scientific_image_method_matrix_incomplete")
    }
    for (overlap in setdiff(protocol$image$overlap_percent, 100L)) {
      unbalanced <- image[
        image$overlap_percent == overlap &
          image$method == "unbalanced_fugw", , drop = FALSE
      ]
      balanced <- image[
        image$overlap_percent == overlap & image$method == "balanced_fgw",
        , drop = FALSE
      ]
      endpoint_denominator <- stats::median(
        balanced$weighted_endpoint_error
      )
      endpoint_ratio <- if (is.finite(endpoint_denominator) &&
          endpoint_denominator > 0) {
        stats::median(unbalanced$weighted_endpoint_error) /
          endpoint_denominator
      } else {
        Inf
      }
      if (nrow(unbalanced) != length(protocol$image$evaluation_seeds) ||
          nrow(balanced) != length(protocol$image$evaluation_seeds) ||
          stats::median(unbalanced$false_occluded_source_mass) >=
            stats::median(balanced$false_occluded_source_mass) ||
          endpoint_ratio >
            protocol$image$partial_overlap_gate$median_endpoint_error_ratio_max ||
          stats::median(unbalanced$recovered_overlap_fraction) <
            protocol$image$partial_overlap_gate$median_recovered_overlap_fraction_min) {
        reasons <- c(
          reasons, paste0("scientific_image_partial_overlap_", overlap)
        )
      }
    }
  }
  if (!nrow(surface) || !all(surface_columns %in% names(surface)) ||
      nrow(surface) != expected_surface_rows ||
      !all(surface$solver_success %in% TRUE) ||
      !all(surface$suite_gate %in% TRUE) ||
      !all(surface$evaluation_split == "frozen_evaluation") ||
      !all(surface$dense_parity_passed %in% TRUE) ||
      any(surface$dense_parity_maximum_error > surface$dense_parity_tolerance) ||
      !all(surface$protocol_schema_version == protocol$schema_version) ||
      !all(surface$protocol_status == "frozen_before_evaluation") ||
      any(!is.finite(as.matrix(surface[
        intersect(surface_metrics, names(surface))
      ])))) {
    reasons <- c(reasons, "scientific_surface_csv_invalid")
  }
  if (!nrow(rank_curve) || !all(rank_columns %in% names(rank_curve)) ||
      nrow(rank_curve) != expected_rank_rows ||
      !all(rank_curve$solver_success %in% TRUE) ||
      !all(rank_curve$numerical_and_geometry_error_separate %in% TRUE) ||
      !setequal(rank_curve$embedding_rank,
                protocol$surface$embedding_curve$ranks) ||
      !all(rank_curve$protocol_schema_version == protocol$schema_version) ||
      !all(rank_curve$protocol_status == "frozen_before_evaluation") ||
      any(!is.finite(as.matrix(rank_curve[
        intersect(surface_metrics, names(rank_curve))
      ])))) {
    reasons <- c(reasons, "scientific_geometry_rank_csv_invalid")
  }
  parity <- summary$surface$dense_parity
  if (!isTRUE(parity$passed) ||
      .mfa_numeric(parity$maximum_parity_error) >
        .mfa_numeric(parity$tolerance) ||
      .mfa_numeric(parity$tolerance) > 1e-6) {
    reasons <- c(reasons, "scientific_dense_parity_invalid")
  }
  files <- c(basename(summary_path), artifacts)
  files <- files[file.exists(file.path(root, files))]
  fingerprint <- .mfa_tree_fingerprint(root, files)
  list(
    passed = !length(unique(reasons)), reasons = unique(reasons),
    platform = platform, root = root, commit = commit,
    fingerprint = fingerprint$sha256, file_count = nrow(fingerprint$files),
    image_rows = completed_image, surface_rows = completed_surface,
    rank_rows = completed_rank
  )
}

.mfa_validate_robustness <- function(summary_path, source) {
  reasons <- character()
  root <- dirname(normalizePath(summary_path, mustWork = TRUE))
  summary <- tryCatch(.mfa_read_json(summary_path), error = function(error) {
    reasons <<- c(reasons, "robustness_summary_invalid")
    list()
  })
  platform <- .mfa_platform(summary$environment$sys_info$sysname)
  if (is.na(platform)) reasons <- c(reasons, "robustness_platform_invalid")
  if (!identical(.mfa_integer(summary$schema_version), 1L) ||
      !identical(.mfa_character(summary$method_family, ""), "moment_fugw") ||
      !identical(.mfa_character(summary$status, ""), "passed") ||
      !isTRUE(summary$admission_evaluable) ||
      !identical(.mfa_character(summary$profile, ""), "release") ||
      !identical(.mfa_character(summary$package_mode, ""), "installed") ||
      !identical(.mfa_character(summary$global_optimality, ""),
                 "not_claimed")) {
    reasons <- c(reasons, "robustness_not_admissible_installed_release")
  }

  candidate <- summary$candidate %||% list()
  commit <- tolower(.mfa_character(candidate$git_commit, ""))
  if (!.mfa_valid_commit(commit) ||
      !identical(commit, source$commit) ||
      !identical(.mfa_character(candidate$package_version, ""),
                 source$version)) {
    reasons <- c(reasons, "robustness_commit_or_version_mismatch")
  }
  if (!isTRUE(candidate$git_provenance_complete) ||
      !identical(candidate$git_dirty, FALSE) ||
      !identical(.mfa_integer(candidate$git_status_entry_count), 0L) ||
      !identical(.mfa_integer(candidate$git_commit_command_status), 0L) ||
      !identical(.mfa_integer(candidate$git_status_command_status), 0L) ||
      !isTRUE(candidate$commit_matches_host)) {
    reasons <- c(reasons, "robustness_source_not_exact_clean")
  }
  artifact <- summary$source_artifact %||% list()
  if (!isTRUE(artifact$bound) ||
      !identical(tolower(.mfa_character(artifact$manifest_commit, "")),
                 source$commit) ||
      !identical(.mfa_character(artifact$artifact_sha256, ""),
                 source$sha256) ||
      !identical(.mfa_character(artifact$manifest_sha256, ""),
                 source$manifest_sha256) ||
      !identical(
        .mfa_character(artifact$artifact, ""),
        paste0("rfugw_", source$version, ".tar.gz")
      )) {
    reasons <- c(reasons, "robustness_source_artifact_mismatch")
  }

  contract_path <- moment_admission_resource(
    "moment-fugw-robustness-contract.json"
  )
  contract <- tryCatch(
    jsonlite::fromJSON(contract_path, simplifyVector = FALSE),
    error = function(error) {
      reasons <<- c(reasons, "robustness_contract_invalid")
      list()
    }
  )
  if (!identical(.mfa_integer(contract$schema_version), 1L) ||
      !identical(.mfa_character(contract$method_family, ""), "moment_fugw") ||
      !identical(.mfa_character(contract$profile, ""), "release") ||
      !identical(.mfa_character(contract$required_package_mode, ""),
                 "installed") ||
      !identical(.mfa_character(summary$contract$artifact, ""),
                 basename(contract_path)) ||
      !identical(.mfa_character(summary$contract$sha256, ""),
                 .mfa_sha256(contract_path)) ||
      !identical(.mfa_integer(summary$contract$schema_version), 1L)) {
    reasons <- c(reasons, "robustness_contract_mismatch")
  }

  contract_files <- contract$test_files %||% list()
  expected_paths <- tryCatch(vapply(
    contract_files,
    function(entry) .mfa_character(entry$path, ""),
    character(1)
  ), error = function(error) character())
  expected_hashes <- tryCatch(vapply(
    expected_paths,
    function(path) .mfa_sha256(moment_admission_repository_file(path)),
    character(1)
  ), error = function(error) {
    reasons <<- c(reasons, "robustness_test_source_missing")
    character()
  })
  observed_files <- summary$tests$files
  if (!is.data.frame(observed_files) ||
      !all(c("path", "sha256") %in% names(observed_files)) ||
      nrow(observed_files) != length(expected_paths) ||
      anyDuplicated(as.character(observed_files$path))) {
    reasons <- c(reasons, "robustness_test_file_receipts_invalid")
  } else {
    observed_hashes <- setNames(
      as.character(observed_files$sha256), as.character(observed_files$path)
    )
    if (!setequal(names(observed_hashes), expected_paths) ||
        length(expected_hashes) != length(expected_paths) ||
        !identical(
          unname(observed_hashes[sort(names(observed_hashes))]),
          unname(expected_hashes[sort(names(expected_hashes))])
        )) {
      reasons <- c(reasons, "robustness_test_source_hash_mismatch")
    }
  }

  cases_name <- .mfa_character(summary$tests$artifact, "")
  cases_path <- file.path(root, cases_name)
  cases <- tryCatch(
    utils::read.csv(cases_path, stringsAsFactors = FALSE),
    error = function(error) data.frame()
  )
  required_columns <- c(
    "file", "case", "passed_expectations", "failures", "errors",
    "warnings", "skips", "elapsed_seconds"
  )
  cases_valid <- identical(cases_name, "cases.csv") &&
    identical(basename(cases_name), cases_name) &&
    nrow(cases) > 0L && all(required_columns %in% names(cases)) &&
    !anyDuplicated(cases[c("file", "case")]) &&
    all(as.character(cases$file) %in% expected_paths)
  count_columns <- c(
    "passed_expectations", "failures", "errors", "warnings", "skips"
  )
  if (cases_valid) {
    counts <- lapply(cases[count_columns], function(value) {
      suppressWarnings(as.numeric(value))
    })
    cases_valid <- all(vapply(counts, function(value) {
      all(is.finite(value)) && all(value >= 0) && all(value == floor(value))
    }, logical(1))) &&
      all(is.finite(suppressWarnings(as.numeric(cases$elapsed_seconds)))) &&
      all(suppressWarnings(as.numeric(cases$elapsed_seconds)) >= 0)
    if (cases_valid) cases[count_columns] <- counts
  }
  if (!cases_valid) reasons <- c(reasons, "robustness_cases_invalid")

  required <- do.call(rbind, lapply(contract_files, function(entry) {
    values <- unlist(entry$required_cases %||% list(), use.names = FALSE)
    data.frame(
      file = rep(.mfa_character(entry$path, ""), length(values)),
      case = as.character(values), stringsAsFactors = FALSE
    )
  }))
  required_pass <- logical(nrow(required))
  if (cases_valid && nrow(required)) {
    required_pass <- vapply(seq_len(nrow(required)), function(index) {
      hit <- cases$file == required$file[[index]] &
        cases$case == required$case[[index]]
      sum(hit) == 1L &&
        cases$passed_expectations[hit] > 0L &&
        sum(unlist(cases[hit, c(
          "failures", "errors", "warnings", "skips"
        )], use.names = FALSE)) == 0L
    }, logical(1))
  }
  if (!nrow(required) || !all(required_pass)) {
    reasons <- c(reasons, "robustness_required_cases_failed_or_missing")
  }
  if (cases_valid &&
      (sum(cases$failures + cases$errors + cases$warnings + cases$skips) != 0L ||
       sum(cases$passed_expectations) <= 0L)) {
    reasons <- c(reasons, "robustness_unexpected_test_outcome")
  }

  property_cases <- unlist(contract$property_cases %||% list(), use.names = TRUE)
  observed_properties <- unlist(summary$properties %||% list(), use.names = TRUE)
  properties_valid <- length(property_cases) > 0L &&
    setequal(names(observed_properties), names(property_cases)) &&
    all(observed_properties %in% TRUE)
  if (properties_valid && cases_valid) {
    properties_valid <- all(vapply(property_cases, function(case_name) {
      hit <- cases$case == case_name
      sum(hit) == 1L && cases$passed_expectations[hit] > 0L &&
        sum(unlist(cases[hit, c(
          "failures", "errors", "warnings", "skips"
        )], use.names = FALSE)) == 0L
    }, logical(1)))
  }
  if (!properties_valid) reasons <- c(reasons, "robustness_properties_invalid")

  if (cases_valid &&
      (!identical(.mfa_integer(summary$tests$observed_cases), nrow(cases)) ||
       !identical(.mfa_integer(summary$tests$required_cases), nrow(required)) ||
       !isTRUE(summary$tests$all_required_cases_present) ||
       !identical(.mfa_integer(summary$tests$passed_expectations),
                  as.integer(sum(cases$passed_expectations))) ||
       !identical(.mfa_integer(summary$tests$failures),
                  as.integer(sum(cases$failures))) ||
       !identical(.mfa_integer(summary$tests$errors),
                  as.integer(sum(cases$errors))) ||
       !identical(.mfa_integer(summary$tests$warnings),
                  as.integer(sum(cases$warnings))) ||
       !identical(.mfa_integer(summary$tests$skips),
                  as.integer(sum(cases$skips)))) ) {
    reasons <- c(reasons, "robustness_summary_counts_not_reproducible")
  }

  expected_threads <- unlist(contract$thread_environment %||% list(),
                             use.names = TRUE)
  observed_threads <- unlist(summary$environment$threads %||% list(),
                             use.names = TRUE)
  if (!setequal(names(observed_threads), names(expected_threads)) ||
      !identical(
        as.character(observed_threads[sort(names(observed_threads))]),
        as.character(expected_threads[sort(names(expected_threads))])
      )) {
    reasons <- c(reasons, "robustness_thread_contract_mismatch")
  }

  retained <- c(basename(summary_path), cases_name)
  retained <- retained[nzchar(retained) & file.exists(file.path(root, retained))]
  fingerprint <- .mfa_tree_fingerprint(root, retained)
  list(
    passed = !length(unique(reasons)), reasons = unique(reasons),
    platform = platform, root = root, commit = commit,
    fingerprint = fingerprint$sha256, file_count = nrow(fingerprint$files),
    observed_cases = if (cases_valid) nrow(cases) else NA_integer_,
    required_properties = length(property_cases)
  )
}

.mfa_validate_core <- function(summary_path, source) {
  reasons <- character()
  root <- dirname(normalizePath(summary_path, mustWork = TRUE))
  summary <- tryCatch(.mfa_read_json(summary_path), error = function(error) {
    reasons <<- c(reasons, "core_summary_invalid")
    list()
  })
  platform <- .mfa_platform(summary$environment$sys_info$sysname)
  if (is.na(platform)) reasons <- c(reasons, "core_platform_invalid")
  if (!identical(.mfa_integer(summary$schema_version), 1L) ||
      !identical(.mfa_character(summary$method_family, ""), "moment_fugw") ||
      !identical(.mfa_character(summary$status, ""), "passed") ||
      !isTRUE(summary$admission_evaluable) ||
      !identical(.mfa_character(summary$profile, ""), "release") ||
      !identical(.mfa_character(summary$package_mode, ""), "installed") ||
      !identical(.mfa_character(summary$global_optimality, ""),
                 "not_claimed")) {
    reasons <- c(reasons, "core_not_admissible_installed_release")
  }

  candidate <- summary$candidate %||% list()
  commit <- tolower(.mfa_character(candidate$git_commit, ""))
  if (!.mfa_valid_commit(commit) ||
      !identical(commit, source$commit) ||
      !identical(.mfa_character(candidate$package_version, ""),
                 source$version)) {
    reasons <- c(reasons, "core_commit_or_version_mismatch")
  }
  if (!isTRUE(candidate$git_provenance_complete) ||
      !identical(candidate$git_dirty, FALSE) ||
      !identical(.mfa_integer(candidate$git_status_entry_count), 0L) ||
      !identical(.mfa_integer(candidate$git_commit_command_status), 0L) ||
      !identical(.mfa_integer(candidate$git_status_command_status), 0L) ||
      !isTRUE(candidate$commit_matches_host)) {
    reasons <- c(reasons, "core_source_not_exact_clean")
  }
  artifact <- summary$source_artifact %||% list()
  if (!isTRUE(artifact$bound) ||
      !identical(tolower(.mfa_character(artifact$manifest_commit, "")),
                 source$commit) ||
      !identical(.mfa_character(artifact$artifact_sha256, ""),
                 source$sha256) ||
      !identical(.mfa_character(artifact$manifest_sha256, ""),
                 source$manifest_sha256) ||
      !identical(
        .mfa_character(artifact$artifact, ""),
        paste0("rfugw_", source$version, ".tar.gz")
      )) {
    reasons <- c(reasons, "core_source_artifact_mismatch")
  }

  contract_path <- moment_admission_resource("moment-fugw-core-contract.json")
  contract <- tryCatch(
    jsonlite::fromJSON(contract_path, simplifyVector = FALSE),
    error = function(error) {
      reasons <<- c(reasons, "core_contract_invalid")
      list()
    }
  )
  required_platforms <- unlist(
    contract$required_platforms %||% list(), use.names = FALSE
  )
  if (!identical(.mfa_integer(contract$schema_version), 1L) ||
      !identical(.mfa_character(contract$method_family, ""), "moment_fugw") ||
      !identical(.mfa_character(contract$profile, ""), "release") ||
      !identical(.mfa_character(contract$required_package_mode, ""),
                 "installed") ||
      !setequal(required_platforms, c("Darwin", "Linux", "Windows")) ||
      !identical(.mfa_character(summary$contract$artifact, ""),
                 basename(contract_path)) ||
      !identical(.mfa_character(summary$contract$sha256, ""),
                 .mfa_sha256(contract_path)) ||
      !identical(.mfa_integer(summary$contract$schema_version), 1L)) {
    reasons <- c(reasons, "core_contract_mismatch")
  }

  validate_hash_receipts <- function(observed, expected_paths, label) {
    expected_paths <- as.character(expected_paths)
    expected_hashes <- tryCatch(vapply(
      expected_paths,
      function(path) .mfa_sha256(moment_admission_repository_file(path)),
      character(1)
    ), error = function(error) character())
    valid <- is.data.frame(observed) &&
      all(c("path", "sha256") %in% names(observed)) &&
      nrow(observed) == length(expected_paths) &&
      !anyDuplicated(as.character(observed$path)) &&
      length(expected_hashes) == length(expected_paths)
    if (valid) {
      observed_hashes <- setNames(
        as.character(observed$sha256), as.character(observed$path)
      )
      valid <- setequal(names(observed_hashes), expected_paths) &&
        identical(
          unname(observed_hashes[sort(names(observed_hashes))]),
          unname(expected_hashes[sort(names(expected_hashes))])
        )
    }
    if (!valid) reasons <<- c(reasons, paste0("core_", label, "_hash_mismatch"))
    valid
  }
  test_files <- contract$test_files %||% list()
  test_paths <- tryCatch(vapply(
    test_files, function(entry) .mfa_character(entry$path, ""), character(1)
  ), error = function(error) character())
  helper_paths <- as.character(unlist(
    contract$helper_files %||% list(), use.names = FALSE
  ))
  executable_paths <- c(
    .mfa_character(contract$accuracy$script, ""),
    .mfa_character(contract$multiscale$script, "")
  )
  input_paths <- .mfa_character(contract$accuracy$design, "")
  validate_hash_receipts(summary$tests$files, test_paths, "test_source")
  validate_hash_receipts(summary$tests$helpers, helper_paths, "helper_source")
  validate_hash_receipts(summary$executables, executable_paths, "executable")
  validate_hash_receipts(summary$inputs, input_paths, "input")

  accuracy_name <- .mfa_character(summary$accuracy$artifact, "")
  accuracy_path <- file.path(root, accuracy_name)
  accuracy <- tryCatch(
    utils::read.csv(accuracy_path, stringsAsFactors = FALSE),
    error = function(error) data.frame()
  )
  accuracy_design <- tryCatch(
    utils::read.csv(
      moment_admission_repository_file(input_paths),
      stringsAsFactors = FALSE
    ),
    error = function(error) data.frame()
  )
  thresholds <- unlist(
    contract$accuracy$thresholds %||% list(), use.names = TRUE
  )
  accuracy_columns <- c(
    "case_id", "seed", "ns", "nt", "zero_weight", "near_zero_weight", "profile",
    "dense_converged", "implicit_converged", "inner_certified",
    "outer_certified", names(thresholds), "pass", "reject_reason"
  )
  accuracy_valid <- identical(accuracy_name, "accuracy.csv") &&
    nrow(accuracy) == .mfa_integer(contract$accuracy$rows) &&
    nrow(accuracy_design) == .mfa_integer(contract$accuracy$rows) &&
    all(accuracy_columns %in% names(accuracy)) &&
    all(names(accuracy_design) %in% names(accuracy)) &&
    !anyDuplicated(accuracy$case_id)
  if (accuracy_valid) {
    expected_pairs <- as.character(unlist(
      contract$accuracy$expected_pairs, use.names = FALSE
    ))
    observed_pairs <- paste0(accuracy$ns, "x", accuracy$nt)
    logical_columns <- c(
      "zero_weight", "near_zero_weight", "dense_converged",
      "implicit_converged", "inner_certified", "outer_certified", "pass"
    )
    accuracy_valid <- setequal(observed_pairs, expected_pairs) &&
      all(vapply(accuracy[logical_columns], is.logical, logical(1))) &&
      all(accuracy$dense_converged) && all(accuracy$implicit_converged) &&
      all(accuracy$inner_certified) && all(accuracy$outer_certified) &&
      all(accuracy$pass) && all(accuracy$profile == "full") &&
      any(accuracy$zero_weight) && any(accuracy$near_zero_weight) &&
      all(!accuracy$zero_weight | !accuracy$near_zero_weight) &&
      all(is.na(accuracy$reject_reason) | !nzchar(accuracy$reject_reason)) &&
      all(vapply(names(accuracy_design), function(column) {
        isTRUE(all.equal(
          accuracy[[column]], accuracy_design[[column]],
          tolerance = 0, check.attributes = FALSE
        ))
      }, logical(1))) &&
      all(accuracy$seed ==
        .mfa_integer(contract$accuracy$seed) + accuracy_design$seed_offset)
  }
  if (accuracy_valid) {
    for (metric in names(thresholds)) {
      values <- suppressWarnings(as.numeric(accuracy[[metric]]))
      if (identical(metric, "direct_oracle_relative_error")) {
        direct <- pmax(accuracy$ns, accuracy$nt) <= 8L
        accuracy_valid <- accuracy_valid && all(is.finite(values[direct])) &&
          all(values[direct] <= thresholds[[metric]]) &&
          all(is.na(values[!direct]) | is.finite(values[!direct])) &&
          all(values[!direct][is.finite(values[!direct])] <= thresholds[[metric]])
      } else {
        accuracy_valid <- accuracy_valid && all(is.finite(values)) &&
          all(values >= 0) && all(values <= thresholds[[metric]])
      }
    }
  }
  if (!accuracy_valid) reasons <- c(reasons, "core_accuracy_csv_invalid")

  metric_maximum <- function(metric) {
    values <- suppressWarnings(as.numeric(accuracy[[metric]]))
    values <- values[is.finite(values)]
    if (length(values)) max(values) else NA_real_
  }
  if (accuracy_valid) {
    accuracy_summary_valid <- isTRUE(summary$accuracy$gate) &&
      identical(.mfa_integer(summary$accuracy$rows), nrow(accuracy)) &&
      identical(.mfa_integer(summary$accuracy$zero_weight_rows),
                as.integer(sum(accuracy$zero_weight))) &&
      identical(.mfa_integer(summary$accuracy$near_zero_weight_rows),
                as.integer(sum(accuracy$near_zero_weight))) &&
      .mfa_same_numeric(
        summary$accuracy$maximum_objective_relative_error,
        metric_maximum("objective_relative_error")
      ) &&
      .mfa_same_numeric(
        summary$accuracy$maximum_component_relative_error,
        metric_maximum("component_relative_error")
      ) &&
      .mfa_same_numeric(
        summary$accuracy$maximum_action_relative_error,
        metric_maximum("action_relative_error")
      ) &&
      .mfa_same_numeric(
        summary$accuracy$maximum_adjoint_relative_error,
        metric_maximum("adjoint_relative_error")
      ) &&
      .mfa_same_numeric(
        summary$accuracy$maximum_mass_relative_error,
        metric_maximum("mass_relative_error")
      ) &&
      .mfa_same_numeric(
        summary$accuracy$maximum_moment_relative_error,
        metric_maximum("moment_relative_error")
      )
    if (!accuracy_summary_valid) {
      reasons <- c(reasons, "core_accuracy_summary_not_reproducible")
    }
  }

  multiscale_name <- .mfa_character(summary$multiscale$artifact, "")
  multiscale_path <- file.path(root, multiscale_name)
  multiscale <- tryCatch(
    utils::read.csv(multiscale_path, stringsAsFactors = FALSE),
    error = function(error) data.frame()
  )
  multiscale_columns <- c(
    "seed", "same_final_contract", "cold_certified",
    "transfer_final_certified", "transfer_multiscale_certified",
    "inner_iteration_reduction", "objective_relative_difference",
    "heldout_score_ratio", "plan_action_relative_difference", "efficacy_gate"
  )
  multiscale_valid <- identical(multiscale_name, "multiscale.csv") &&
    nrow(multiscale) == .mfa_integer(contract$multiscale$rows) &&
    all(multiscale_columns %in% names(multiscale)) &&
    !anyDuplicated(multiscale$seed) &&
    setequal(
      as.integer(multiscale$seed),
      as.integer(unlist(contract$multiscale$seeds, use.names = FALSE))
    )
  if (multiscale_valid) {
    logical_columns <- c(
      "same_final_contract", "cold_certified", "transfer_final_certified",
      "transfer_multiscale_certified", "efficacy_gate"
    )
    numeric_columns <- c(
      "inner_iteration_reduction", "objective_relative_difference",
      "heldout_score_ratio", "plan_action_relative_difference"
    )
    multiscale_valid <-
      all(vapply(multiscale[logical_columns], is.logical, logical(1))) &&
      all(vapply(multiscale[numeric_columns], function(value) {
        all(is.finite(suppressWarnings(as.numeric(value))))
      }, logical(1))) &&
      all(multiscale$same_final_contract) && all(multiscale$cold_certified) &&
      all(multiscale$transfer_final_certified) &&
      all(multiscale$transfer_multiscale_certified) &&
      all(multiscale$efficacy_gate)
  }
  if (multiscale_valid) {
    improved_fraction <- mean(multiscale$inner_iteration_reduction > 0)
    median_reduction <- stats::median(multiscale$inner_iteration_reduction)
    maximum_objective <- max(abs(multiscale$objective_relative_difference))
    minimum_heldout <- min(multiscale$heldout_score_ratio)
    maximum_action <- max(multiscale$plan_action_relative_difference)
    multiscale_valid <-
      improved_fraction >= .mfa_numeric(
        contract$multiscale$minimum_improved_fraction
      ) &&
      median_reduction >= .mfa_numeric(
        contract$multiscale$minimum_median_inner_iteration_reduction
      ) &&
      maximum_objective <= .mfa_numeric(
        contract$multiscale$maximum_objective_relative_difference
      ) &&
      minimum_heldout >= .mfa_numeric(
        contract$multiscale$minimum_heldout_score_ratio
      ) &&
      maximum_action <= .mfa_numeric(
        contract$multiscale$maximum_identity_action_difference
      )
  }
  if (!multiscale_valid) reasons <- c(reasons, "core_multiscale_csv_invalid")
  if (multiscale_valid) {
    multiscale_summary_valid <- isTRUE(summary$multiscale$gate) &&
      identical(.mfa_integer(summary$multiscale$rows), nrow(multiscale)) &&
      .mfa_same_numeric(
        summary$multiscale$improved_fraction, improved_fraction
      ) &&
      .mfa_same_numeric(
        summary$multiscale$median_inner_iteration_reduction, median_reduction
      ) &&
      .mfa_same_numeric(
        summary$multiscale$maximum_objective_relative_difference,
        maximum_objective
      ) &&
      .mfa_same_numeric(
        summary$multiscale$minimum_heldout_score_ratio, minimum_heldout
      ) &&
      .mfa_same_numeric(
        summary$multiscale$maximum_identity_action_difference, maximum_action
      )
    if (!multiscale_summary_valid) {
      reasons <- c(reasons, "core_multiscale_summary_not_reproducible")
    }
  }

  cases_name <- .mfa_character(summary$tests$artifact, "")
  cases_path <- file.path(root, cases_name)
  cases <- tryCatch(
    utils::read.csv(cases_path, stringsAsFactors = FALSE),
    error = function(error) data.frame()
  )
  case_columns <- c(
    "file", "case", "passed_expectations", "failures", "errors",
    "warnings", "skips", "elapsed_seconds"
  )
  cases_valid <- identical(cases_name, "cases.csv") && nrow(cases) > 0L &&
    all(case_columns %in% names(cases)) &&
    !anyDuplicated(cases[c("file", "case")]) &&
    all(as.character(cases$file) %in% test_paths)
  count_columns <- c(
    "passed_expectations", "failures", "errors", "warnings", "skips"
  )
  if (cases_valid) {
    counts <- lapply(cases[count_columns], function(value) {
      suppressWarnings(as.numeric(value))
    })
    cases_valid <- all(vapply(counts, function(value) {
      all(is.finite(value)) && all(value >= 0) && all(value == floor(value))
    }, logical(1))) &&
      all(is.finite(suppressWarnings(as.numeric(cases$elapsed_seconds)))) &&
      all(suppressWarnings(as.numeric(cases$elapsed_seconds)) >= 0)
    if (cases_valid) cases[count_columns] <- counts
  }
  required <- tryCatch(do.call(rbind, lapply(test_files, function(entry) {
    values <- unlist(entry$required_cases %||% list(), use.names = FALSE)
    data.frame(
      file = rep(.mfa_character(entry$path, ""), length(values)),
      case = as.character(values), stringsAsFactors = FALSE
    )
  })), error = function(error) data.frame())
  if (cases_valid && nrow(required)) {
    required_pass <- vapply(seq_len(nrow(required)), function(index) {
      hit <- cases$file == required$file[[index]] &
        cases$case == required$case[[index]]
      sum(hit) == 1L && cases$passed_expectations[hit] > 0L &&
        sum(unlist(cases[hit, c(
          "failures", "errors", "warnings", "skips"
        )], use.names = FALSE)) == 0L
    }, logical(1))
    cases_valid <- all(required_pass) &&
      sum(cases$failures + cases$errors + cases$warnings + cases$skips) == 0L &&
      sum(cases$passed_expectations) > 0L
  } else {
    cases_valid <- FALSE
  }
  if (!cases_valid) reasons <- c(reasons, "core_required_cases_failed_or_missing")
  if (cases_valid) {
    case_summary_valid <- isTRUE(summary$tests$gate) &&
      identical(.mfa_integer(summary$tests$observed_cases), nrow(cases)) &&
      identical(.mfa_integer(summary$tests$required_cases), nrow(required)) &&
      identical(.mfa_integer(summary$tests$passed_expectations),
                as.integer(sum(cases$passed_expectations))) &&
      identical(.mfa_integer(summary$tests$failures), 0L) &&
      identical(.mfa_integer(summary$tests$errors), 0L) &&
      identical(.mfa_integer(summary$tests$warnings), 0L) &&
      identical(.mfa_integer(summary$tests$skips), 0L)
    if (!case_summary_valid) {
      reasons <- c(reasons, "core_test_summary_not_reproducible")
    }
  }

  expected_threads <- unlist(
    contract$thread_environment %||% list(), use.names = TRUE
  )
  observed_threads <- unlist(
    summary$environment$threads %||% list(), use.names = TRUE
  )
  if (!setequal(names(observed_threads), names(expected_threads)) ||
      !identical(
        as.character(observed_threads[sort(names(observed_threads))]),
        as.character(expected_threads[sort(names(expected_threads))])
      )) {
    reasons <- c(reasons, "core_thread_contract_mismatch")
  }

  retained <- c(
    basename(summary_path), accuracy_name, multiscale_name, cases_name
  )
  retained <- retained[nzchar(retained) & file.exists(file.path(root, retained))]
  fingerprint <- .mfa_tree_fingerprint(root, retained)
  list(
    passed = !length(unique(reasons)), reasons = unique(reasons),
    platform = platform, root = root, commit = commit,
    fingerprint = fingerprint$sha256, file_count = nrow(fingerprint$files),
    accuracy_rows = if (accuracy_valid) nrow(accuracy) else NA_integer_,
    multiscale_rows = if (multiscale_valid) nrow(multiscale) else NA_integer_,
    required_cases = if (nrow(required)) nrow(required) else NA_integer_
  )
}

.mfa_validate_release_trust <- function(report_path, source) {
  reasons <- character()
  report_root <- dirname(normalizePath(report_path, mustWork = TRUE))
  bundle_root <- dirname(dirname(report_root))
  report <- tryCatch(.mfa_read_json(report_path), error = function(error) {
    reasons <<- c(reasons, "release_trust_report_invalid")
    list()
  })
  platform <- .mfa_platform(report$sys_info$sysname)
  if (is.na(platform)) reasons <- c(reasons, "release_trust_platform_invalid")
  if (!identical(.mfa_integer(report$schema_version), 2L) ||
      !identical(.mfa_character(report$package, ""), "rfugw") ||
      !identical(.mfa_character(report$scope, ""), "release") ||
      !identical(.mfa_character(report$evidence_channel, ""), "hosted") ||
      !identical(
        .mfa_character(report$publication_status, ""),
        "not_evaluated_by_numerical_trust_job"
      )) {
    reasons <- c(reasons, "release_trust_identity_invalid")
  }
  commit <- .mfa_character(report$commit)
  if (!.mfa_valid_commit(commit) ||
      !identical(tolower(commit), source$commit) ||
      !identical(.mfa_character(report$expected_host_commit), source$commit)) {
    reasons <- c(reasons, "release_trust_commit_mismatch")
  }
  if (!isTRUE(report$git_provenance_complete) ||
      !identical(report$working_tree_dirty, FALSE) ||
      !identical(.mfa_integer(report$git_status_entry_count), 0L) ||
      !identical(.mfa_integer(report$git_commit_command_status), 0L) ||
      !identical(.mfa_integer(report$git_status_command_status), 0L) ||
      !isTRUE(report$commit_matches_host) ||
      !isTRUE(report$exact_commit_evidence) ||
      !nzchar(.mfa_character(report$hosted_run, ""))) {
    reasons <- c(reasons, "release_trust_not_exact_hosted")
  }
  source_artifact <- file.path(
    bundle_root, paste0("rfugw_", source$version, ".tar.gz")
  )
  copied_manifest <- file.path(bundle_root, "release-artifact.json")
  check_log <- file.path(bundle_root, "rfugw.Rcheck", "00check.log")
  trust_txt <- file.path(report_root, "numerical-trust-release.txt")
  dossier <- file.path(report_root, "release-dossier-release.md")
  if (!identical(.mfa_sha256(source_artifact), source$sha256) ||
      !identical(.mfa_sha256(copied_manifest), source$manifest_sha256)) {
    reasons <- c(reasons, "release_trust_source_artifact_mismatch")
  }
  artifact_hashes <- if (is.data.frame(report$artifacts)) {
    as.character(report$artifacts$sha256)
  } else if (is.list(report$artifacts)) {
    unname(unlist(lapply(report$artifacts, function(artifact) {
      artifact$sha256 %||% character()
    }), use.names = FALSE))
  } else {
    character()
  }
  if (!source$sha256 %in% artifact_hashes) {
    reasons <- c(reasons, "release_trust_artifact_ledger_mismatch")
  }
  check_lines <- if (file.exists(check_log)) {
    readLines(check_log, warn = FALSE)
  } else {
    character()
  }
  status_lines <- grep("^Status:", check_lines, value = TRUE)
  if (!length(check_lines) ||
      !any(grepl("^[*] checking tests .* OK$", check_lines)) ||
      !any(grepl("^[*] DONE$", check_lines)) ||
      length(status_lines) > 1L ||
      any(grepl("ERROR|WARNING", status_lines))) {
    reasons <- c(reasons, "release_trust_r_cmd_check_failed")
  }
  required_files <- c(
    report_path, trust_txt, dossier, source_artifact, copied_manifest, check_log
  )
  if (!all(file.exists(required_files))) {
    reasons <- c(reasons, "release_trust_bundle_incomplete")
  }
  fingerprint <- .mfa_tree_fingerprint(
    bundle_root,
    substring(required_files[file.exists(required_files)],
              nchar(bundle_root) + 2L)
  )
  list(
    passed = !length(unique(reasons)), reasons = unique(reasons),
    platform = platform, root = bundle_root, commit = commit,
    fingerprint = fingerprint$sha256,
    file_count = nrow(fingerprint$files),
    check_status = if (length(status_lines)) status_lines[[1L]] else "clean"
  )
}

.mfa_validate_review <- function(path, source_commit, evidence_fingerprint) {
  if (is.null(path) || !nzchar(path) || !file.exists(path)) {
    return(list(
      passed = FALSE, reasons = "independent_review_missing",
      reviewer = NA_character_, sha256 = NA_character_
    ))
  }
  reasons <- character()
  review <- tryCatch(.mfa_read_json(path), error = function(error) {
    reasons <<- c(reasons, "independent_review_invalid_json")
    list()
  })
  if (!identical(as.integer(review$schema_version %||% NA_integer_), 1L)) {
    reasons <- c(reasons, "independent_review_schema_invalid")
  }
  if (!identical(tolower(as.character(
      review$source_commit %||% ""
    )[[1L]]), source_commit) ||
      !identical(as.character(
        review$evidence_fingerprint_sha256 %||% ""
      )[[1L]], evidence_fingerprint)) {
    reasons <- c(reasons, "independent_review_evidence_binding_mismatch")
  }
  reviewer <- review$reviewer %||% list()
  reviewer_name <- .mfa_character(reviewer$name, "")
  if (!nzchar(reviewer_name) ||
      !nzchar(.mfa_character(reviewer$affiliation, "")) ||
      !nzchar(.mfa_character(reviewer$contact_or_profile, "")) ||
      !nzchar(.mfa_character(reviewer$conflicts_disclosed, "")) ||
      !isTRUE(reviewer$independent_of_implementation) ||
      !identical(reviewer$implementation_contributor, FALSE)) {
    reasons <- c(reasons, "independent_reviewer_identity_or_conflict_invalid")
  }
  scopes <- review$scopes %||% list()
  for (scope in c("numerical", "performance", "scientific")) {
    decision <- scopes[[scope]] %||% list()
    if (!identical(as.character(decision$decision %||% "")[[1L]],
                   "approve") ||
        !nzchar(as.character(decision$notes %||% "")[[1L]])) {
      reasons <- c(reasons, paste0("independent_review_scope_", scope))
    }
  }
  signed_at <- .mfa_character(review$signed_at_utc, "")
  if (!grepl(
      "^[0-9]{4}-[0-9]{2}-[0-9]{2}T[0-9]{2}:[0-9]{2}:[0-9]{2}Z$",
      signed_at
    )) {
    reasons <- c(reasons, "independent_review_signature_time_missing")
  }
  required_limitations <- c(
    "global_optimality_not_claimed", "adaptive_support_not_admitted",
    "gpu_and_100k_scale_not_claimed"
  )
  reproduction <- review$reproduction %||% list()
  if (!nzchar(.mfa_character(reproduction$commands_or_run_url, "")) ||
      !setequal(as.character(
        reproduction$limitations_confirmed %||% character()
      ), required_limitations)) {
    reasons <- c(reasons, "independent_review_reproduction_incomplete")
  }
  list(
    passed = !length(unique(reasons)), reasons = unique(reasons),
    reviewer = reviewer_name, sha256 = .mfa_sha256(path)
  )
}

moment_fugw_verify_admission <- function(
    manifest_path,
    evidence_root,
    review_path = NULL,
    mode = c("candidate", "promotion"),
    expected_commit = NULL,
    output_path = NULL) {
  .mfa_require_packages()
  mode <- match.arg(mode)
  contract <- jsonlite::fromJSON(
    moment_admission_resource("moment-fugw-thresholds.json"),
    simplifyVector = TRUE
  )
  source <- .mfa_validate_source(manifest_path, expected_commit)
  source_commit <- source$commit
  evidence_root <- normalizePath(evidence_root, mustWork = TRUE)
  paths <- list.files(
    evidence_root, recursive = TRUE, full.names = TRUE,
    all.files = TRUE, no.. = TRUE
  )
  normalized <- gsub("\\\\", "/", paths)
  performance_paths <- paths[grepl(
    "moment-fugw-performance/current/summary\\.json$", normalized
  )]
  scientific_paths <- paths[grepl(
    paste0(
      "moment-fugw-scientific/",
      "moment-fugw-scientific-summary-evaluation\\.json$"
    ), normalized
  )]
  robustness_paths <- paths[grepl(
    "moment-fugw-robustness/summary\\.json$", normalized
  )]
  core_paths <- paths[grepl(
    "moment-fugw-core/summary\\.json$", normalized
  )]
  trust_paths <- paths[grepl(
    "numerical-trust/numerical-trust-release\\.json$", normalized
  )]
  performance <- lapply(
    performance_paths, .mfa_validate_performance,
    source_commit = source_commit, source_version = source$version,
    contract = contract
  )
  scientific <- lapply(
    scientific_paths, .mfa_validate_science,
    source_commit = source_commit, source_version = source$version
  )
  robustness <- lapply(
    robustness_paths, .mfa_validate_robustness, source = source
  )
  core <- lapply(
    core_paths, .mfa_validate_core, source = source
  )
  release_trust <- lapply(
    trust_paths, .mfa_validate_release_trust, source = source
  )
  performance_platforms <- vapply(
    performance, function(value) value$platform, character(1)
  )
  scientific_platforms <- vapply(
    scientific, function(value) value$platform, character(1)
  )
  robustness_platforms <- vapply(
    robustness, function(value) value$platform, character(1)
  )
  core_platforms <- vapply(
    core, function(value) value$platform, character(1)
  )
  trust_platforms <- vapply(
    release_trust, function(value) value$platform, character(1)
  )
  platform_reasons <- character()
  if (!setequal(performance_platforms, c("Darwin", "Linux")) ||
      anyDuplicated(performance_platforms)) {
    platform_reasons <- c(
      platform_reasons, "performance_platform_matrix_incomplete"
    )
  }
  if (!setequal(scientific_platforms, c("Darwin", "Linux", "Windows")) ||
      anyDuplicated(scientific_platforms)) {
    platform_reasons <- c(
      platform_reasons, "scientific_platform_matrix_incomplete"
    )
  }
  if (!setequal(robustness_platforms, c("Darwin", "Linux", "Windows")) ||
      anyDuplicated(robustness_platforms)) {
    platform_reasons <- c(
      platform_reasons, "robustness_platform_matrix_incomplete"
    )
  }
  if (!setequal(core_platforms, c("Darwin", "Linux", "Windows")) ||
      anyDuplicated(core_platforms)) {
    platform_reasons <- c(
      platform_reasons, "core_platform_matrix_incomplete"
    )
  }
  if (!setequal(trust_platforms, c("Darwin", "Linux", "Windows")) ||
      anyDuplicated(trust_platforms)) {
    platform_reasons <- c(
      platform_reasons, "release_trust_platform_matrix_incomplete"
    )
  }
  performance <- performance[order(performance_platforms)]
  scientific <- scientific[order(scientific_platforms)]
  robustness <- robustness[order(robustness_platforms)]
  core <- core[order(core_platforms)]
  release_trust <- release_trust[order(trust_platforms)]
  evidence_lines <- c(
    paste0("source_commit=", source_commit),
    paste0("source_sha256=", source$sha256),
    paste0("thresholds_md5=", unname(tools::md5sum(
      moment_admission_resource("moment-fugw-thresholds.json")
    )[[1L]])),
    vapply(performance, function(value) paste0(
      "performance_", value$platform, "=", value$fingerprint
    ), character(1)),
    vapply(scientific, function(value) paste0(
      "scientific_", value$platform, "=", value$fingerprint
    ), character(1)),
    vapply(robustness, function(value) paste0(
      "robustness_", value$platform, "=", value$fingerprint
    ), character(1)),
    vapply(core, function(value) paste0(
      "core_", value$platform, "=", value$fingerprint
    ), character(1)),
    vapply(release_trust, function(value) paste0(
      "release_trust_", value$platform, "=", value$fingerprint
    ), character(1))
  )
  evidence_fingerprint <- digest::digest(
    paste(evidence_lines, collapse = "\n"),
    algo = "sha256", serialize = FALSE
  )
  review <- .mfa_validate_review(
    review_path, source_commit, evidence_fingerprint
  )
  technical_reasons <- unique(c(
    source$reasons,
    unlist(lapply(performance, `[[`, "reasons"), use.names = FALSE),
    unlist(lapply(scientific, `[[`, "reasons"), use.names = FALSE),
    unlist(lapply(robustness, `[[`, "reasons"), use.names = FALSE),
    unlist(lapply(core, `[[`, "reasons"), use.names = FALSE),
    unlist(lapply(release_trust, `[[`, "reasons"), use.names = FALSE),
    platform_reasons
  ))
  technical_passed <- !length(technical_reasons)
  promotion_eligible <- technical_passed && isTRUE(review$passed)
  passed <- if (identical(mode, "candidate")) {
    technical_passed
  } else {
    promotion_eligible
  }
  receipt <- list(
    schema_version = 1L,
    method_family = "moment_fugw",
    mode = mode,
    status = if (promotion_eligible) {
      "promotion_eligible_independently_reviewed"
    } else if (technical_passed) {
      "technical_candidate_passed_review_pending"
    } else {
      "rejected"
    },
    passed = passed,
    technical_gates_passed = technical_passed,
    independent_review_passed = isTRUE(review$passed),
    promotion_eligible = promotion_eligible,
    supported = FALSE,
    public_status = "experimental",
    global_optimality = "not_claimed",
    source = source,
    performance = performance,
    scientific = scientific,
    robustness = robustness,
    core = core,
    release_trust = release_trust,
    evidence_fingerprint_sha256 = evidence_fingerprint,
    review = review,
    technical_failure_reasons = technical_reasons,
    unmet_promotion_gates = unique(c(
      if (!technical_passed) technical_reasons,
      if (!isTRUE(review$passed)) review$reasons
    )),
    generated_at_utc = format(
      Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"
    )
  )
  if (!is.null(output_path)) {
    dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
    jsonlite::write_json(
      receipt, output_path, auto_unbox = TRUE, pretty = TRUE,
      null = "null", na = "null", dataframe = "rows"
    )
  }
  receipt
}
