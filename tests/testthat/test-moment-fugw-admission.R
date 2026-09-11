.mfa_test_write_json <- function(value, path) {
  jsonlite::write_json(
    value, path, auto_unbox = TRUE, pretty = TRUE, null = "null",
    na = "null", dataframe = "rows"
  )
}

.mfa_test_runs <- function(contract, raw_dir) {
  rows <- list()
  add <- function(curve, action, n, phase = "measured", repetition = 1L,
                  rank = contract$problem$structure_rank,
                  batch = NA_integer_, seconds = 0.1,
                  allocation = NA_real_) {
    increment <- 8e7 * (as.numeric(n) / 625)^1.2
    rows[[length(rows) + 1L]] <<- data.frame(
      profile = "full",
      curve = curve,
      phase = phase,
      repetition = as.integer(repetition),
      action = action,
      method = "fugw_factorized",
      n_source = as.integer(n),
      n_target = as.integer(n),
      structure_rank = as.integer(rank),
      batch = as.integer(batch),
      process_status = 0L,
      infrastructure_error = FALSE,
      valid = TRUE,
      certified = identical(curve, "end_to_end_convergence"),
      timing_eligible = identical(curve, "end_to_end_convergence") &&
        identical(phase, "measured"),
      reportable_seconds = as.numeric(seconds),
      peak_rss_bytes = 1.2e8 + increment,
      incremental_peak_rss_bytes = increment,
      r_visible_allocation_bytes = as.numeric(allocation),
      dense_sentinel_passed = TRUE,
      full_block_calls = 0L,
      max_tile_elements = as.numeric(min(n, 256L) * n),
      full_matrix_elements = as.numeric(n) * as.numeric(n),
      stringsAsFactors = FALSE
    )
  }
  profile <- contract$profiles$full
  for (n in as.integer(profile$domain_sizes)) {
    add(
      "fixed_work_scaling", "solve", n, "warmup", 1L,
      seconds = 0.98 * (n / 625)^2
    )
    for (repetition in seq_len(profile$repetitions)) {
      add(
        "fixed_work_scaling", "solve", n, "measured", repetition,
        seconds = c(0.98, 1, 1.02)[[repetition]] * (n / 625)^2
      )
    }
  }
  for (n in as.integer(profile$domain_sizes)) {
    add(
      "r_visible_allocation", "solve", n,
      allocation = 1e7 + n * n
    )
    for (action in c("objective", "apply", "adjoint")) {
      add(
        "matrix_free_memory", action, n, batch = 1L,
        allocation = 1e6 + n * match(action, c("objective", "apply", "adjoint"))
      )
    }
  }
  for (rank in as.integer(contract$scaling$structure_ranks)) {
    add("rank_scaling", "solve", profile$rank_size, rank = rank)
  }
  for (n in as.integer(profile$converged_sizes)) {
    add("end_to_end_convergence", "solve", n, "warmup", 1L)
    for (repetition in seq_len(profile$repetitions)) {
      add(
        "end_to_end_convergence", "solve", n, "measured", repetition
      )
    }
  }
  for (batch in as.integer(contract$scaling$heldout_batches)) {
    for (action in c("apply", "adjoint")) {
      add(
        "operator_throughput", action, profile$operator_size,
        batch = batch
      )
    }
  }
  output <- do.call(rbind, rows)
  stopifnot(nrow(output) == 60L)
  stems <- sprintf("worker-%03d", seq_len(nrow(output)) + 1L)
  output$spec_artifact <- file.path(raw_dir, paste0(stems, ".spec.json"))
  output$result_artifact <- file.path(raw_dir, paste0(stems, ".result.json"))
  output$time_artifact <- file.path(raw_dir, paste0(stems, ".time.log"))
  output
}

.mfa_test_performance <- function(root, platform, commit, version, contract) {
  current <- file.path(
    root, paste0("performance-", tolower(platform)),
    "moment-fugw-performance", "current"
  )
  raw_dir <- file.path(current, "raw")
  dir.create(raw_dir, recursive = TRUE)
  workers <- sprintf("worker-%03d", seq_len(61L))
  for (stem in workers) {
    writeLines("{}", file.path(raw_dir, paste0(stem, ".spec.json")))
    writeLines("{}", file.path(raw_dir, paste0(stem, ".result.json")))
    writeLines(character(), file.path(raw_dir, paste0(stem, ".stdout.log")))
    writeLines("rss=1", file.path(raw_dir, paste0(stem, ".time.log")))
  }
  runs <- .mfa_test_runs(contract, raw_dir)
  utils::write.csv(runs, file.path(current, "runs.csv"), row.names = FALSE)
  backend <- if (identical(platform, "Darwin")) {
    "macos_time_l_bytes"
  } else {
    "linux_time_v_kib"
  }
  required_gates <- c(
    "design", "full_profile_complete", "exact_commit_provenance",
    "clean_exact_source_tree", "release_hardware",
    "fixed_work_runtime_envelope", "near_linear_peak_memory",
    "prohibited_dense_allocations_absent", "r_visible_allocation_complete",
    "rank_curve_complete", "operator_batch_curve_complete",
    "larger_certified_endpoint"
  )
  gates <- as.list(rep(TRUE, length(required_gates)))
  names(gates) <- required_gates
  summary <- list(
    schema_version = 1L,
    profile = "full",
    package_mode = "installed",
    admission_evaluable = TRUE,
    rss_backend = backend,
    fixed_work_runtime_exponent = 2,
    fixed_work_runtime_exponent_fit = "largest_three_points",
    fixed_work_runtime_exponent_all_points = 2,
    incremental_memory_exponent = 1.2,
    absolute_peak_memory_exponent = 1,
    completed_fixed_work_sizes = as.integer(contract$scaling$domain_sizes),
    configured_full_sizes = as.integer(contract$scaling$domain_sizes),
    largest_completed_certified_endpoint = 512L,
    claim_ceiling_n = 512L,
    gates = gates
  )
  meta <- list(
    benchmark_schema_version = 1L,
    package = "rfugw",
    version = version,
    profile = "moment_fugw_full_installed",
    package_mode = "installed",
    commit = commit,
    git_dirty = FALSE,
    git_status_entry_count = 0L,
    git_provenance_complete = TRUE,
    git_commit_command_status = 0L,
    git_status_command_status = 0L,
    sysname = platform,
    release = "fixture-release",
    machine = "fixture-machine",
    nodename = "fixture-node",
    r_version = R.version.string,
    r_platform = R.version$platform,
    rss_backend = backend,
    thresholds_md5 = unname(tools::md5sum(
      bench_test_resource("moment-fugw-thresholds.json")
    )[[1L]]),
    raw_artifacts = sort(list.files(raw_dir, recursive = TRUE))
  )
  .mfa_test_write_json(summary, file.path(current, "summary.json"))
  .mfa_test_write_json(meta, file.path(current, "meta.json"))
  current
}

.mfa_test_science_tables <- function(protocol) {
  image <- expand.grid(
    seed = as.integer(protocol$image$evaluation_seeds),
    overlap_percent = as.integer(protocol$image$overlap_percent),
    method = c("balanced_fgw", "unbalanced_fugw"),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  image$solver_success <- TRUE
  image$suite_gate <- TRUE
  image$evaluation_split <- "frozen_evaluation"
  image$protocol_schema_version <- protocol$schema_version
  image$protocol_status <- "frozen_before_evaluation"
  partial <- image$overlap_percent < 100L
  image$transported_mass <- 0.9
  image$discarded_source_mass <- 0.1
  image$discarded_target_mass <- 0.1
  image$false_occluded_source_mass <- ifelse(
    partial & image$method == "unbalanced_fugw", 0.1,
    ifelse(partial, 0.2, 0)
  )
  image$weighted_endpoint_error <- ifelse(
    image$method == "unbalanced_fugw", 0.9, 1
  )
  image$heldout_rmse <- 0.1
  image$heldout_correlation <- 0.95
  image$local_distortion <- 0.05
  image$coupling_l1_discrepancy <- 1e-8
  image$geometry_relative_error <- 0
  image$recovered_overlap_fraction <- 0.9

  surface <- data.frame(
    solver_success = TRUE,
    suite_gate = TRUE,
    evaluation_split = "frozen_evaluation",
    dense_parity_passed = TRUE,
    dense_parity_maximum_error = 1e-8,
    dense_parity_tolerance = 1e-6,
    protocol_schema_version = protocol$schema_version,
    protocol_status = "frozen_before_evaluation",
    geometry_relative_error = 0,
    geometry_heldout_error = 0,
    heldout_correlation = 0.99,
    heldout_r2 = 0.98,
    expected_geodesic_endpoint_error = 0.02,
    transported_mass = 0.98,
    discarded_source_mass = 0.02,
    discarded_target_mass = 0.02,
    coupling_l1_discrepancy = 1e-8,
    seed = as.integer(protocol$surface$evaluation_seeds),
    stringsAsFactors = FALSE
  )
  ranks <- as.integer(protocol$surface$embedding_curve$ranks)
  rank_curve <- data.frame(
    solver_success = TRUE,
    numerical_and_geometry_error_separate = TRUE,
    embedding_rank = ranks,
    protocol_schema_version = protocol$schema_version,
    protocol_status = "frozen_before_evaluation",
    geometry_relative_error = seq(0.3, 0.05, length.out = length(ranks)),
    geometry_heldout_error = seq(0.3, 0.05, length.out = length(ranks)),
    heldout_correlation = seq(0.9, 0.99, length.out = length(ranks)),
    heldout_r2 = seq(0.8, 0.98, length.out = length(ranks)),
    expected_geodesic_endpoint_error = seq(0.1, 0.02, length.out = length(ranks)),
    transported_mass = 0.98,
    discarded_source_mass = 0.02,
    discarded_target_mass = 0.02,
    coupling_l1_discrepancy = 1e-8,
    stringsAsFactors = FALSE
  )
  list(image = image, surface = surface, rank = rank_curve)
}

.mfa_test_science <- function(root, platform, commit, version, protocol) {
  output <- file.path(
    root, paste0("science-", tolower(platform)), "moment-fugw-scientific"
  )
  dir.create(output, recursive = TRUE)
  tables <- .mfa_test_science_tables(protocol)
  paths <- c(
    image = "moment-fugw-image-evaluation.csv",
    surface = "moment-fugw-surface-evaluation.csv",
    rank = "moment-fugw-geometry-rank-evaluation.csv"
  )
  utils::write.csv(
    tables$image, file.path(output, paths[["image"]]), row.names = FALSE
  )
  utils::write.csv(
    tables$surface, file.path(output, paths[["surface"]]), row.names = FALSE
  )
  utils::write.csv(
    tables$rank, file.path(output, paths[["rank"]]), row.names = FALSE
  )
  summary <- list(
    schema_version = 1L,
    status = "passed",
    admission_evaluable = TRUE,
    split = "evaluation",
    package_mode = "installed",
    protocol = list(
      schema_version = protocol$schema_version,
      status = "frozen_before_evaluation",
      md5 = unname(tools::md5sum(
        bench_test_resource("moment-fugw-validation-protocol.json")
      )[[1L]])
    ),
    candidate = list(
      git_commit = commit,
      git_dirty = FALSE,
      git_status_entry_count = 0L,
      git_provenance_complete = TRUE,
      git_commit_command_status = 0L,
      git_status_command_status = 0L,
      package_version = version
    ),
    environment = list(sys_info = list(sysname = platform)),
    image = list(
      gate = TRUE, requested_rows = nrow(tables$image),
      completed_rows = nrow(tables$image), failed_rows = 0L,
      artifact = paths[["image"]]
    ),
    surface = list(
      gate = TRUE, requested_rows = nrow(tables$surface),
      completed_rows = nrow(tables$surface), failed_rows = 0L,
      dense_parity = list(
        passed = TRUE, maximum_parity_error = 1e-8, tolerance = 1e-6
      ),
      artifact = paths[["surface"]]
    ),
    geometry_rank_curve = list(
      gate = TRUE, completed_rows = nrow(tables$rank),
      artifact = paths[["rank"]]
    )
  )
  .mfa_test_write_json(
    summary,
    file.path(output, "moment-fugw-scientific-summary-evaluation.json")
  )
  output
}

.mfa_test_repository_file <- function(relative_path) {
  bench_root <- dirname(bench_test_resource(
    "moment-fugw-robustness-contract.json"
  ))
  package_parent <- file.path(bench_root, "..", "..")
  candidates <- c(
    file.path(package_parent, relative_path),
    file.path(package_parent, "00_pkg_src", "rfugw", relative_path)
  )
  candidates <- candidates[file.exists(candidates)]
  if (!length(candidates)) {
    stop("Admission fixture source not found: ", relative_path, call. = FALSE)
  }
  normalizePath(candidates[[1L]], mustWork = TRUE)
}

.mfa_test_robustness <- function(
    root, platform, commit, version, artifact, manifest) {
  output <- file.path(
    root, paste0("robustness-", tolower(platform)),
    "moment-fugw-robustness"
  )
  dir.create(output, recursive = TRUE)
  contract_path <- bench_test_resource(
    "moment-fugw-robustness-contract.json"
  )
  contract <- jsonlite::fromJSON(contract_path, simplifyVector = FALSE)
  required <- do.call(rbind, lapply(contract$test_files, function(entry) {
    values <- unlist(entry$required_cases, use.names = FALSE)
    data.frame(
      file = rep(as.character(entry$path[[1L]]), length(values)),
      case = as.character(values), stringsAsFactors = FALSE
    )
  }))
  cases <- data.frame(
    required,
    passed_expectations = rep(1L, nrow(required)),
    failures = rep(0L, nrow(required)),
    errors = rep(0L, nrow(required)),
    warnings = rep(0L, nrow(required)),
    skips = rep(0L, nrow(required)),
    elapsed_seconds = rep(0.01, nrow(required)),
    stringsAsFactors = FALSE
  )
  cases_path <- file.path(output, "cases.csv")
  utils::write.csv(cases, cases_path, row.names = FALSE)
  test_files <- lapply(contract$test_files, function(entry) {
    relative <- as.character(entry$path[[1L]])
    path <- .mfa_test_repository_file(relative)
    list(
      path = relative,
      sha256 = unname(digest::digest(path, algo = "sha256", file = TRUE))
    )
  })
  property_names <- names(contract$property_cases)
  properties <- stats::setNames(
    as.list(rep(TRUE, length(property_names))), property_names
  )
  summary <- list(
    schema_version = 1L,
    method_family = "moment_fugw",
    status = "passed",
    admission_evaluable = TRUE,
    profile = "release",
    package_mode = "installed",
    candidate = list(
      git_commit = commit,
      git_dirty = FALSE,
      git_status_entry_count = 0L,
      git_provenance_complete = TRUE,
      git_commit_command_status = 0L,
      git_status_command_status = 0L,
      commit_matches_host = TRUE,
      package_version = version
    ),
    source_artifact = list(
      manifest = basename(manifest),
      manifest_sha256 = unname(digest::digest(
        manifest, algo = "sha256", file = TRUE
      )),
      artifact = basename(artifact),
      artifact_sha256 = unname(digest::digest(
        artifact, algo = "sha256", file = TRUE
      )),
      manifest_commit = commit,
      bound = TRUE
    ),
    contract = list(
      artifact = basename(contract_path),
      sha256 = unname(digest::digest(
        contract_path, algo = "sha256", file = TRUE
      )),
      schema_version = 1L
    ),
    tests = list(
      artifact = basename(cases_path),
      observed_cases = nrow(cases),
      required_cases = nrow(required),
      all_required_cases_present = TRUE,
      passed_expectations = sum(cases$passed_expectations),
      failures = 0L,
      errors = 0L,
      warnings = 0L,
      skips = 0L,
      files = test_files
    ),
    properties = properties,
    environment = list(
      r_version = R.version.string,
      platform = R.version$platform,
      sys_info = list(sysname = platform),
      threads = contract$thread_environment
    ),
    global_optimality = "not_claimed"
  )
  .mfa_test_write_json(summary, file.path(output, "summary.json"))
  output
}

.mfa_test_core <- function(
    root, platform, commit, version, artifact, manifest) {
  output <- file.path(
    root, paste0("core-", tolower(platform)), "moment-fugw-core"
  )
  dir.create(output, recursive = TRUE)
  contract_path <- bench_test_resource("moment-fugw-core-contract.json")
  contract <- jsonlite::fromJSON(contract_path, simplifyVector = FALSE)

  accuracy_design_path <- .mfa_test_repository_file(
    as.character(contract$accuracy$design[[1L]])
  )
  accuracy_design <- utils::read.csv(
    accuracy_design_path, stringsAsFactors = FALSE
  )
  pairs <- unlist(contract$accuracy$expected_pairs, use.names = FALSE)
  accuracy <- data.frame(
    accuracy_design,
    seed = as.integer(contract$accuracy$seed) +
      accuracy_design$seed_offset,
    profile = rep("full", length(pairs)),
    dense_converged = TRUE,
    implicit_converged = TRUE,
    inner_certified = TRUE,
    outer_certified = TRUE,
    objective_relative_error = 1e-10,
    component_relative_error = 1e-10,
    direct_oracle_relative_error = ifelse(
      pmax(accuracy_design$ns, accuracy_design$nt) <= 8L,
      1e-10, NA_real_
    ),
    action_relative_error = 1e-10,
    adjoint_relative_error = 1e-10,
    mass_relative_error = 1e-10,
    moment_relative_error = 1e-10,
    pass = TRUE,
    reject_reason = NA_character_,
    stringsAsFactors = FALSE
  )
  accuracy_path <- file.path(output, "accuracy.csv")
  utils::write.csv(accuracy, accuracy_path, row.names = FALSE)

  seeds <- as.integer(unlist(contract$multiscale$seeds, use.names = FALSE))
  multiscale <- data.frame(
    seed = seeds,
    same_final_contract = TRUE,
    cold_certified = TRUE,
    transfer_final_certified = TRUE,
    transfer_multiscale_certified = TRUE,
    inner_iteration_reduction = 0.3,
    objective_relative_difference = rep(c(-1e-10, 1e-10),
                                         length.out = length(seeds)),
    heldout_score_ratio = 1,
    plan_action_relative_difference = 1e-10,
    efficacy_gate = TRUE,
    stringsAsFactors = FALSE
  )
  multiscale_path <- file.path(output, "multiscale.csv")
  utils::write.csv(multiscale, multiscale_path, row.names = FALSE)

  required <- do.call(rbind, lapply(contract$test_files, function(entry) {
    values <- unlist(entry$required_cases, use.names = FALSE)
    data.frame(
      file = rep(as.character(entry$path[[1L]]), length(values)),
      case = as.character(values), stringsAsFactors = FALSE
    )
  }))
  cases <- data.frame(
    required,
    passed_expectations = rep(1L, nrow(required)),
    failures = rep(0L, nrow(required)),
    errors = rep(0L, nrow(required)),
    warnings = rep(0L, nrow(required)),
    skips = rep(0L, nrow(required)),
    elapsed_seconds = rep(0.01, nrow(required)),
    stringsAsFactors = FALSE
  )
  cases_path <- file.path(output, "cases.csv")
  utils::write.csv(cases, cases_path, row.names = FALSE)

  hash_receipts <- function(paths) lapply(paths, function(relative) {
    relative <- as.character(relative)
    path <- .mfa_test_repository_file(relative)
    list(
      path = relative,
      sha256 = unname(digest::digest(path, algo = "sha256", file = TRUE))
    )
  })
  test_paths <- vapply(
    contract$test_files,
    function(entry) as.character(entry$path[[1L]]), character(1)
  )
  helper_paths <- unlist(contract$helper_files, use.names = FALSE)
  executable_paths <- c(
    as.character(contract$accuracy$script[[1L]]),
    as.character(contract$multiscale$script[[1L]])
  )
  input_paths <- as.character(contract$accuracy$design[[1L]])
  summary <- list(
    schema_version = 1L,
    method_family = "moment_fugw",
    status = "passed",
    admission_evaluable = TRUE,
    profile = "release",
    package_mode = "installed",
    candidate = list(
      git_commit = commit,
      git_dirty = FALSE,
      git_status_entry_count = 0L,
      git_provenance_complete = TRUE,
      git_commit_command_status = 0L,
      git_status_command_status = 0L,
      commit_matches_host = TRUE,
      package_version = version
    ),
    source_artifact = list(
      manifest = basename(manifest),
      manifest_sha256 = unname(digest::digest(
        manifest, algo = "sha256", file = TRUE
      )),
      artifact = basename(artifact),
      artifact_sha256 = unname(digest::digest(
        artifact, algo = "sha256", file = TRUE
      )),
      manifest_commit = commit,
      bound = TRUE
    ),
    contract = list(
      artifact = basename(contract_path),
      sha256 = unname(digest::digest(
        contract_path, algo = "sha256", file = TRUE
      )),
      schema_version = 1L
    ),
    accuracy = list(
      gate = TRUE,
      artifact = basename(accuracy_path),
      rows = nrow(accuracy),
      zero_weight_rows = sum(accuracy$zero_weight),
      near_zero_weight_rows = sum(accuracy$near_zero_weight),
      maximum_objective_relative_error = 1e-10,
      maximum_component_relative_error = 1e-10,
      maximum_action_relative_error = 1e-10,
      maximum_adjoint_relative_error = 1e-10,
      maximum_mass_relative_error = 1e-10,
      maximum_moment_relative_error = 1e-10
    ),
    multiscale = list(
      gate = TRUE,
      artifact = basename(multiscale_path),
      rows = nrow(multiscale),
      improved_fraction = 1,
      median_inner_iteration_reduction = 0.3,
      maximum_objective_relative_difference = 1e-10,
      minimum_heldout_score_ratio = 1,
      maximum_identity_action_difference = 1e-10
    ),
    tests = list(
      gate = TRUE,
      artifact = basename(cases_path),
      observed_cases = nrow(cases),
      required_cases = nrow(required),
      passed_expectations = sum(cases$passed_expectations),
      failures = 0L,
      errors = 0L,
      warnings = 0L,
      skips = 0L,
      files = hash_receipts(test_paths),
      helpers = hash_receipts(helper_paths)
    ),
    executables = hash_receipts(executable_paths),
    inputs = hash_receipts(input_paths),
    environment = list(
      r_version = R.version.string,
      platform = R.version$platform,
      sys_info = list(sysname = platform),
      threads = contract$thread_environment
    ),
    global_optimality = "not_claimed"
  )
  .mfa_test_write_json(summary, file.path(output, "summary.json"))
  output
}

.mfa_test_release_trust <- function(
    root, platform, commit, version, artifact, manifest) {
  bundle <- file.path(root, paste0("release-", tolower(platform)))
  trust <- file.path(bundle, ".gate", "numerical-trust")
  check <- file.path(bundle, "rfugw.Rcheck")
  dir.create(trust, recursive = TRUE)
  dir.create(check, recursive = TRUE)
  artifact_copy <- file.path(bundle, paste0("rfugw_", version, ".tar.gz"))
  manifest_copy <- file.path(bundle, "release-artifact.json")
  stopifnot(file.copy(artifact, artifact_copy))
  stopifnot(file.copy(manifest, manifest_copy))
  writeLines(
    c("* checking tests ... OK", "* DONE", "Status: 1 NOTE"),
    file.path(check, "00check.log")
  )
  writeLines("hosted release evidence", file.path(
    trust, "numerical-trust-release.txt"
  ))
  writeLines("# Hosted release dossier", file.path(
    trust, "release-dossier-release.md"
  ))
  report <- list(
    schema_version = 2L,
    package = "rfugw",
    scope = "release",
    evidence_channel = "hosted",
    commit = commit,
    working_tree_dirty = FALSE,
    git_provenance_complete = TRUE,
    git_status_entry_count = 0L,
    git_commit_command_status = 0L,
    git_status_command_status = 0L,
    expected_host_commit = commit,
    commit_matches_host = TRUE,
    exact_commit_evidence = TRUE,
    hosted_run = paste0("fixture-", tolower(platform)),
    publication_status = "not_evaluated_by_numerical_trust_job",
    sys_info = list(sysname = platform),
    artifacts = list(list(
      path = basename(artifact_copy),
      sha256 = unname(digest::digest(
        artifact_copy, algo = "sha256", file = TRUE
      ))
    ))
  )
  .mfa_test_write_json(
    report, file.path(trust, "numerical-trust-release.json")
  )
  bundle
}

.mfa_test_fixture <- function() {
  root <- tempfile("moment-fugw-admission-")
  dir.create(root)
  source_root <- file.path(root, "source")
  evidence_root <- file.path(root, "evidence")
  dir.create(source_root)
  dir.create(evidence_root)
  commit <- paste(rep("a", 40L), collapse = "")
  version <- "0.1.0"
  artifact <- file.path(source_root, paste0("rfugw_", version, ".tar.gz"))
  writeBin(charToRaw("exact canonical source artifact"), artifact)
  manifest <- list(
    schema_version = 1L,
    package = "rfugw",
    version = version,
    artifact = basename(artifact),
    sha256 = unname(digest::digest(artifact, algo = "sha256", file = TRUE)),
    bytes = as.numeric(file.info(artifact)$size),
    source_commit = commit,
    source_tree_clean = TRUE,
    git_provenance_complete = TRUE,
    git_status_entry_count = 0L,
    git_commit_command_status = 0L,
    git_status_command_status = 0L,
    exact_commit_evidence = TRUE
  )
  manifest_path <- file.path(source_root, "release-artifact.json")
  .mfa_test_write_json(manifest, manifest_path)
  contract <- jsonlite::fromJSON(
    bench_test_resource("moment-fugw-thresholds.json"), simplifyVector = TRUE
  )
  protocol <- jsonlite::fromJSON(
    bench_test_resource("moment-fugw-validation-protocol.json"),
    simplifyVector = TRUE
  )
  performance <- vapply(c("Darwin", "Linux"), function(platform) {
    .mfa_test_performance(
      evidence_root, platform, commit, version, contract
    )
  }, character(1))
  science <- vapply(c("Darwin", "Linux", "Windows"), function(platform) {
    .mfa_test_science(evidence_root, platform, commit, version, protocol)
  }, character(1))
  robustness <- vapply(c("Darwin", "Linux", "Windows"), function(platform) {
    .mfa_test_robustness(
      evidence_root, platform, commit, version, artifact, manifest_path
    )
  }, character(1))
  core <- vapply(c("Darwin", "Linux", "Windows"), function(platform) {
    .mfa_test_core(
      evidence_root, platform, commit, version, artifact, manifest_path
    )
  }, character(1))
  release_trust <- vapply(
    c("Darwin", "Linux", "Windows"), function(platform) {
      .mfa_test_release_trust(
        evidence_root, platform, commit, version, artifact, manifest_path
      )
    }, character(1)
  )
  list(
    root = root, manifest = manifest_path, artifact = artifact,
    evidence = evidence_root, commit = commit, version = version,
    performance = performance, science = science, robustness = robustness,
    core = core, release_trust = release_trust
  )
}

test_that("Moment-FUGW admission is exact-artifact bound and fail closed", {
  skip_if_not_installed("digest")
  skip_if_not_installed("jsonlite")
  admission <- new.env(parent = globalenv())
  sys.source(
    bench_test_resource("moment_fugw_admission.R"), envir = admission
  )
  fixture <- .mfa_test_fixture()
  on.exit(unlink(fixture$root, recursive = TRUE, force = TRUE), add = TRUE)
  receipt_path <- file.path(fixture$root, "candidate-receipt.json")

  candidate <- admission$moment_fugw_verify_admission(
    fixture$manifest, fixture$evidence, mode = "candidate",
    expected_commit = fixture$commit, output_path = receipt_path
  )
  expect_true(
    candidate$passed,
    info = paste(candidate$technical_failure_reasons, collapse = ", ")
  )
  expect_true(candidate$technical_gates_passed)
  expect_false(candidate$independent_review_passed)
  expect_false(candidate$promotion_eligible)
  expect_identical(candidate$status, "technical_candidate_passed_review_pending")
  expect_false(candidate$supported)
  expect_identical(candidate$public_status, "experimental")
  expect_identical(candidate$global_optimality, "not_claimed")
  expect_length(candidate$robustness, 3L)
  expect_true(all(vapply(
    candidate$robustness,
    function(value) value$passed && value$required_properties == 15L,
    logical(1)
  )))
  expect_length(candidate$core, 3L)
  expect_true(all(vapply(
    candidate$core,
    function(value) {
      value$passed && value$accuracy_rows == 12L &&
        value$multiscale_rows == 10L && value$required_cases == 27L
    },
    logical(1)
  )))
  expect_true(file.exists(receipt_path))

  wrapper <- testthat::test_path(
    "..", "..", "tools", "numerical-trust",
    "verify-moment-fugw-admission.R"
  )
  if (file.exists(wrapper)) {
    cli_receipt <- file.path(fixture$root, "cli-candidate-receipt.json")
    cli_output <- suppressWarnings(system2(
      file.path(R.home("bin"), "Rscript"),
      c(
        wrapper, "--mode=candidate",
        paste0("--manifest=", fixture$manifest),
        paste0("--evidence-root=", fixture$evidence),
        paste0("--commit=", fixture$commit),
        paste0("--output=", cli_receipt)
      ),
      stdout = TRUE, stderr = TRUE
    ))
    cli_status <- attr(cli_output, "status", exact = TRUE)
    if (is.null(cli_status)) cli_status <- 0L
    expect_identical(as.integer(cli_status), 0L, info = paste(cli_output))
    expect_true(file.exists(cli_receipt))
  }

  review <- list(
    schema_version = 1L,
    source_commit = fixture$commit,
    evidence_fingerprint_sha256 = candidate$evidence_fingerprint_sha256,
    reviewer = list(
      name = "Independent Reviewer",
      affiliation = "External Numerical Methods Lab",
      contact_or_profile = "https://example.org/reviewer",
      independent_of_implementation = TRUE,
      implementation_contributor = FALSE,
      conflicts_disclosed = "none"
    ),
    scopes = list(
      numerical = list(decision = "approve", notes = "Oracle replay passed."),
      performance = list(decision = "approve", notes = "RSS replay passed."),
      scientific = list(decision = "approve", notes = "Held-out replay passed.")
    ),
    reproduction = list(
      commands_or_run_url = "https://example.org/immutable-run",
      limitations_confirmed = c(
        "global_optimality_not_claimed", "adaptive_support_not_admitted",
        "gpu_and_100k_scale_not_claimed"
      )
    ),
    signed_at_utc = "2026-08-27T12:00:00Z"
  )
  review_path <- file.path(fixture$root, "independent-review.json")
  .mfa_test_write_json(review, review_path)
  promotion <- admission$moment_fugw_verify_admission(
    fixture$manifest, fixture$evidence, review_path = review_path,
    mode = "promotion", expected_commit = fixture$commit
  )
  expect_true(promotion$passed)
  expect_true(promotion$promotion_eligible)
  expect_true(promotion$independent_review_passed)
  expect_false(promotion$supported)

  darwin_runs_path <- file.path(fixture$performance[["Darwin"]], "runs.csv")
  original_runs <- readLines(darwin_runs_path, warn = FALSE)
  runs <- utils::read.csv(darwin_runs_path, stringsAsFactors = FALSE)
  runs$reportable_seconds[
    runs$curve == "fixed_work_scaling" & runs$phase == "measured" &
      runs$n_source == 5000L
  ] <- 1000
  utils::write.csv(runs, darwin_runs_path, row.names = FALSE)
  drift <- admission$moment_fugw_verify_admission(
    fixture$manifest, fixture$evidence, mode = "candidate",
    expected_commit = fixture$commit
  )
  expect_false(drift$passed)
  expect_true(
    "performance_reported_slopes_not_reproducible" %in%
      drift$technical_failure_reasons
  )
  writeLines(original_runs, darwin_runs_path)

  raw_path <- file.path(
    fixture$performance[["Linux"]], "raw", "worker-061.result.json"
  )
  unlink(raw_path)
  missing_raw <- admission$moment_fugw_verify_admission(
    fixture$manifest, fixture$evidence, mode = "candidate",
    expected_commit = fixture$commit
  )
  expect_false(missing_raw$passed)
  expect_true(
    "performance_raw_artifact_manifest_mismatch" %in%
      missing_raw$technical_failure_reasons
  )
  writeLines("{}", raw_path)
  writeLines("not-json", raw_path)
  invalid_raw <- admission$moment_fugw_verify_admission(
    fixture$manifest, fixture$evidence, mode = "candidate",
    expected_commit = fixture$commit
  )
  expect_false(invalid_raw$passed)
  expect_true(
    "performance_raw_json_receipt_invalid" %in%
      invalid_raw$technical_failure_reasons
  )
  writeLines("{}", raw_path)

  linux_meta_path <- file.path(fixture$performance[["Linux"]], "meta.json")
  linux_meta <- jsonlite::fromJSON(linux_meta_path, simplifyVector = TRUE)
  linux_meta$commit <- paste(rep("b", 40L), collapse = "")
  .mfa_test_write_json(linux_meta, linux_meta_path)
  wrong_commit <- admission$moment_fugw_verify_admission(
    fixture$manifest, fixture$evidence, mode = "candidate",
    expected_commit = fixture$commit
  )
  expect_false(wrong_commit$passed)
  expect_true(
    "performance_commit_mismatch" %in% wrong_commit$technical_failure_reasons
  )
  linux_meta$commit <- fixture$commit
  .mfa_test_write_json(linux_meta, linux_meta_path)

  check_log <- file.path(
    fixture$release_trust[["Darwin"]], "rfugw.Rcheck", "00check.log"
  )
  original_check <- readLines(check_log, warn = FALSE)
  writeLines(c(original_check, "Status: 1 ERROR"), check_log)
  bad_check <- admission$moment_fugw_verify_admission(
    fixture$manifest, fixture$evidence, mode = "candidate",
    expected_commit = fixture$commit
  )
  expect_false(bad_check$passed)
  expect_true(
    "release_trust_r_cmd_check_failed" %in%
      bad_check$technical_failure_reasons
  )
  writeLines(original_check, check_log)

  image_path <- file.path(
    fixture$science[["Windows"]], "moment-fugw-image-evaluation.csv"
  )
  original_image <- readLines(image_path, warn = FALSE)
  image <- utils::read.csv(image_path, stringsAsFactors = FALSE)
  image$false_occluded_source_mass[
    image$method == "unbalanced_fugw" & image$overlap_percent == 80L
  ] <- 0.5
  utils::write.csv(image, image_path, row.names = FALSE)
  bad_science <- admission$moment_fugw_verify_admission(
    fixture$manifest, fixture$evidence, mode = "candidate",
    expected_commit = fixture$commit
  )
  expect_false(bad_science$passed)
  expect_true(
    "scientific_image_partial_overlap_80" %in%
      bad_science$technical_failure_reasons
  )
  writeLines(original_image, image_path)

  robustness_cases_path <- file.path(
    fixture$robustness[["Windows"]], "cases.csv"
  )
  original_robustness_cases <- readLines(robustness_cases_path, warn = FALSE)
  robustness_cases <- utils::read.csv(
    robustness_cases_path, stringsAsFactors = FALSE
  )
  robustness_cases$failures[[1L]] <- 1L
  utils::write.csv(
    robustness_cases, robustness_cases_path, row.names = FALSE
  )
  bad_robustness <- admission$moment_fugw_verify_admission(
    fixture$manifest, fixture$evidence, mode = "candidate",
    expected_commit = fixture$commit
  )
  expect_false(bad_robustness$passed)
  expect_true(any(c(
    "robustness_required_cases_failed_or_missing",
    "robustness_unexpected_test_outcome",
    "robustness_summary_counts_not_reproducible"
  ) %in% bad_robustness$technical_failure_reasons))
  writeLines(original_robustness_cases, robustness_cases_path)

  core_accuracy_path <- file.path(
    fixture$core[["Linux"]], "accuracy.csv"
  )
  original_core_accuracy <- readLines(core_accuracy_path, warn = FALSE)
  core_accuracy <- utils::read.csv(
    core_accuracy_path, stringsAsFactors = FALSE
  )
  core_accuracy$objective_relative_error[[1L]] <- 1
  utils::write.csv(core_accuracy, core_accuracy_path, row.names = FALSE)
  bad_core <- admission$moment_fugw_verify_admission(
    fixture$manifest, fixture$evidence, mode = "candidate",
    expected_commit = fixture$commit
  )
  expect_false(bad_core$passed)
  expect_true(
    "core_accuracy_csv_invalid" %in% bad_core$technical_failure_reasons
  )
  writeLines(original_core_accuracy, core_accuracy_path)

  core_multiscale_path <- file.path(
    fixture$core[["Darwin"]], "multiscale.csv"
  )
  original_core_multiscale <- readLines(core_multiscale_path, warn = FALSE)
  core_multiscale <- utils::read.csv(
    core_multiscale_path, stringsAsFactors = FALSE
  )
  core_multiscale$objective_relative_difference[[1L]] <- -1
  core_multiscale$plan_action_relative_difference[[2L]] <- 1
  utils::write.csv(
    core_multiscale, core_multiscale_path, row.names = FALSE
  )
  bad_multiscale <- admission$moment_fugw_verify_admission(
    fixture$manifest, fixture$evidence, mode = "candidate",
    expected_commit = fixture$commit
  )
  expect_false(bad_multiscale$passed)
  expect_true(
    "core_multiscale_csv_invalid" %in%
      bad_multiscale$technical_failure_reasons
  )
  writeLines(original_core_multiscale, core_multiscale_path)

  writeBin(charToRaw("tampered source artifact"), fixture$artifact)
  tampered <- admission$moment_fugw_verify_admission(
    fixture$manifest, fixture$evidence, mode = "candidate",
    expected_commit = fixture$commit
  )
  expect_false(tampered$passed)
  expect_true(
    "source_artifact_sha256_mismatch" %in%
      tampered$technical_failure_reasons
  )
  writeBin(charToRaw("exact canonical source artifact"), fixture$artifact)

  review$reviewer$implementation_contributor <- TRUE
  .mfa_test_write_json(review, review_path)
  conflicted <- admission$moment_fugw_verify_admission(
    fixture$manifest, fixture$evidence, review_path = review_path,
    mode = "promotion", expected_commit = fixture$commit
  )
  expect_false(conflicted$passed)
  expect_true(conflicted$technical_gates_passed)
  expect_false(conflicted$independent_review_passed)
  expect_true(
    "independent_reviewer_identity_or_conflict_invalid" %in%
      conflicted$unmet_promotion_gates
  )
})
