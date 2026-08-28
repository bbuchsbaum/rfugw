# Matrix-free Moment-FUGW benchmark helpers.
# This file is sourced by both the orchestrator and one-shot worker processes.

bench_moment_fugw_resource <- function(name) {
  installed <- system.file("bench", name, package = "rfugw")
  if (nzchar(installed) && file.exists(installed)) return(installed)
  candidates <- c(
    file.path("inst", "bench", name),
    file.path(dirname(sys.frame(1)$ofile %||% ""), name)
  )
  hit <- candidates[file.exists(candidates)]
  if (!length(hit)) {
    stop("Moment-FUGW benchmark resource not found: ", name, call. = FALSE)
  }
  normalizePath(hit[[1L]], mustWork = TRUE)
}

bench_moment_fugw_profile <- function(profile = c("smoke", "pr", "full")) {
  profile <- match.arg(profile)
  contract <- bench_read_moment_fugw_thresholds()
  out <- contract$profiles[[profile]]
  if (is.null(out)) stop("Unknown Moment-FUGW profile.", call. = FALSE)
  out$profile <- profile
  out$structure_ranks <- as.integer(contract$scaling$structure_ranks)
  out$heldout_batches <- as.integer(contract$scaling$heldout_batches)
  out$seed <- as.integer(contract$problem$seed)
  out$threads <- as.integer(contract$problem$threads)
  out$structure_rank <- as.integer(contract$problem$structure_rank)
  out$feature_dimension <- as.integer(contract$problem$feature_dimension)
  out$epsilon <- as.numeric(contract$problem$epsilon)
  out$reg_marginals <- as.numeric(contract$problem$reg_marginals)
  out$tol <- as.numeric(contract$problem$tol)
  out$tol_ot <- as.numeric(contract$problem$tol_ot)
  out
}

.bench_with_seed <- function(seed, code) {
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv)
  on.exit({
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  set.seed(as.integer(seed))
  force(code)
}

.bench_guard_cost_operator <- function(cost, label, audit) {
  original <- cost$block
  n_source <- cost$n_source
  n_target <- cost$n_target
  cost$block <- function(rows, columns) {
    elements <- as.double(length(rows)) * as.double(length(columns))
    audit$maximum_requested_block <- max(
      audit$maximum_requested_block, elements
    )
    audit$block_calls <- audit$block_calls + 1L
    if (elements >= as.double(n_source) * as.double(n_target)) {
      audit$full_block_calls <- audit$full_block_calls + 1L
      stop(
        "Prohibited full matrix request for ", label, ".",
        call. = FALSE
      )
    }
    original(rows, columns)
  }
  cost
}

bench_make_moment_fugw_problem <- function(
    n,
    structure_rank = 4L,
    feature_dimension = 8L,
    seed = 20260827L,
    guard_dense = TRUE) {
  n <- as.integer(n)
  structure_rank <- as.integer(structure_rank)
  feature_dimension <- as.integer(feature_dimension)
  if (!is.finite(n) || n < 2L || !is.finite(structure_rank) ||
      structure_rank < 1L || !is.finite(feature_dimension) ||
      feature_dimension < 1L) {
    stop("Invalid Moment-FUGW benchmark dimensions.", call. = FALSE)
  }
  generated <- .bench_with_seed(seed + 1009L * n + structure_rank, {
    source_coordinates <- matrix(rnorm(n * structure_rank), n, structure_rank)
    source_coordinates <- source_coordinates / sqrt(structure_rank)
    target_coordinates <- source_coordinates +
      matrix(rnorm(n * structure_rank, sd = 0.03), n, structure_rank)
    source_features <- matrix(rnorm(n * feature_dimension), n, feature_dimension)
    target_features <- source_features +
      matrix(rnorm(n * feature_dimension, sd = 0.05), n, feature_dimension)
    list(
      source_coordinates = source_coordinates,
      target_coordinates = target_coordinates,
      source_features = source_features,
      target_features = target_features
    )
  })
  audit <- new.env(parent = emptyenv())
  audit$block_calls <- 0L
  audit$full_block_calls <- 0L
  audit$maximum_requested_block <- 0
  Cx <- sqeuclidean_cost(generated$source_coordinates)
  Cy <- sqeuclidean_cost(generated$target_coordinates)
  M <- sqeuclidean_cost(
    generated$source_features, generated$target_features
  )
  if (isTRUE(guard_dense)) {
    Cx <- .bench_guard_cost_operator(Cx, "Cx", audit)
    Cy <- .bench_guard_cost_operator(Cy, "Cy", audit)
    M <- .bench_guard_cost_operator(M, "M", audit)
  }
  list(
    n_source = n,
    n_target = n,
    structure_rank = structure_rank,
    feature_dimension = feature_dimension,
    Cx = Cx,
    Cy = Cy,
    M = M,
    wx = rep(1 / n, n),
    wy = rep(1 / n, n),
    audit = audit,
    seed = as.integer(seed)
  )
}

bench_solve_moment_fugw <- function(
    problem,
    max_iter,
    max_iter_ot,
    tol,
    tol_ot,
    epsilon,
    reg_marginals,
    block_size) {
  fugw_factorized(
    Cx = problem$Cx,
    Cy = problem$Cy,
    wx = problem$wx,
    wy = problem$wy,
    M = problem$M,
    epsilon = epsilon,
    reg_marginals = reg_marginals,
    max_iter = as.integer(max_iter),
    max_iter_ot = as.integer(max_iter_ot),
    tol = tol,
    tol_ot = tol_ot,
    check_every = 1L,
    plan = "operator",
    block_size = as.integer(block_size),
    certify = TRUE
  )
}

bench_moment_fugw_dense_sentinel <- function(fit, problem) {
  memory <- fit$runtime_provenance$memory_contract %||% list()
  plans <- fit$plans %||% list()
  reasons <- character()
  if (!setequal(
    as.character(memory$prohibited_dense_allocations %||% character()),
    c("Cx", "Cy", "M", "P", "Q")
  )) reasons <- c(reasons, "missing_prohibited_allocation_contract")
  if (isTRUE(memory$dense_allocations_requested)) {
    reasons <- c(reasons, "dense_allocation_requested")
  }
  if (!.bench_moment_fugw_plan_is_implicit(plans$sample) ||
      !.bench_moment_fugw_plan_is_implicit(plans$feature)) {
    reasons <- c(reasons, "nonimplicit_returned_plan")
  }
  if ((problem$audit$full_block_calls %||% 0L) != 0L) {
    reasons <- c(reasons, "full_cost_block_requested")
  }
  full <- as.double(problem$n_source) * as.double(problem$n_target)
  max_tile <- as.double(memory$max_tile_elements %||% Inf)
  if (!is.finite(max_tile) || max_tile > full) {
    reasons <- c(reasons, "tile_exceeds_full_shape")
  }
  list(
    passed = !length(reasons),
    reasons = unique(reasons),
    full_block_calls = as.integer(problem$audit$full_block_calls %||% 0L),
    maximum_requested_block = as.double(
      problem$audit$maximum_requested_block %||% 0
    ),
    max_tile_elements = max_tile,
    full_matrix_elements = full
  )
}

.bench_rfugw_internal <- function(name) {
  getFromNamespace(name, "rfugw")
}

bench_recompute_moment_fugw_objective <- function(bundle) {
  fit <- bundle$fit
  problem <- bundle$problem
  source_terms <- .bench_rfugw_internal(".cost_factor_pair")(problem$Cx)
  target_terms <- .bench_rfugw_internal(".cost_factor_pair")(problem$Cy)
  .bench_rfugw_internal(".fugw_factorized_objective")(
    fit$implicit_plans$sample,
    fit$implicit_plans$feature,
    problem$M,
    fit$moments$feature,
    source_terms$left,
    target_terms$left,
    problem$wx,
    problem$wy,
    fit$reg_marginals,
    fit$regularization,
    fit$feature_weight,
    fit$structure_weight
  )
}

bench_moment_fugw_map_batch <- function(n, batch, seed) {
  index <- seq_len(as.double(n) * as.double(batch))
  matrix(
    sin(index * 0.017 + as.integer(seed) * 1e-6),
    nrow = as.integer(n), ncol = as.integer(batch)
  )
}

bench_loglog_exponent <- function(size, value) {
  keep <- is.finite(size) & size > 0 & is.finite(value) & value > 0
  size <- as.numeric(size[keep])
  value <- as.numeric(value[keep])
  if (length(unique(size)) < 2L) return(NA_real_)
  unname(stats::coef(stats::lm(log2(value) ~ log2(size)))[[2L]])
}

bench_detect_rss_backend <- function(sysname = Sys.info()[["sysname"]]) {
  if (!file.exists("/usr/bin/time")) return("unavailable")
  if (identical(sysname, "Darwin")) return("macos_time_l_bytes")
  if (identical(sysname, "Linux")) return("linux_time_v_kib")
  "unsupported"
}

bench_parse_peak_rss <- function(lines, sysname = Sys.info()[["sysname"]]) {
  lines <- as.character(lines)
  if (identical(sysname, "Darwin")) {
    hit <- grep("maximum resident set size", lines, value = TRUE)
    if (!length(hit)) return(NA_real_)
    value <- sub(
      "^[[:space:]]*([0-9]+)[[:space:]]+maximum resident set size.*$",
      "\\1", hit[[length(hit)]]
    )
    return(suppressWarnings(as.numeric(value)))
  }
  if (identical(sysname, "Linux")) {
    hit <- grep("Maximum resident set size [(]kbytes[)]", lines, value = TRUE)
    if (!length(hit)) return(NA_real_)
    value <- sub("^.*:[[:space:]]*([0-9]+)[[:space:]]*$", "\\1", hit[[1L]])
    return(suppressWarnings(as.numeric(value)) * 1024)
  }
  NA_real_
}

bench_run_fresh_process <- function(
    spec,
    artifact_dir,
    label,
    worker = bench_moment_fugw_resource("moment_fugw_worker.R")) {
  dir.create(artifact_dir, recursive = TRUE, showWarnings = FALSE)
  safe_label <- gsub("[^A-Za-z0-9_.-]", "_", label)
  spec_path <- file.path(artifact_dir, paste0(safe_label, ".spec.json"))
  result_path <- file.path(artifact_dir, paste0(safe_label, ".result.json"))
  stdout_path <- file.path(artifact_dir, paste0(safe_label, ".stdout.log"))
  time_path <- file.path(artifact_dir, paste0(safe_label, ".time.log"))
  jsonlite::write_json(spec, spec_path, auto_unbox = TRUE, pretty = TRUE,
                       null = "null")
  backend <- bench_detect_rss_backend()
  if (backend %in% c("unavailable", "unsupported")) {
    stop("Fresh-process RSS is unavailable on this platform.", call. = FALSE)
  }
  time_flag <- if (identical(Sys.info()[["sysname"]], "Darwin")) "-l" else "-v"
  rscript <- file.path(R.home("bin"), "Rscript")
  args <- c(
    time_flag,
    shQuote(rscript), shQuote(worker), shQuote(spec_path), shQuote(result_path)
  )
  env <- c(
    paste0("OMP_NUM_THREADS=", as.integer(spec$threads %||% 1L)),
    "OPENBLAS_NUM_THREADS=1", "MKL_NUM_THREADS=1",
    "VECLIB_MAXIMUM_THREADS=1"
  )
  status <- suppressWarnings(system2(
    "/usr/bin/time", args,
    stdout = stdout_path, stderr = time_path, env = env
  ))
  if (is.null(status)) status <- 0L
  time_lines <- if (file.exists(time_path)) {
    readLines(time_path, warn = FALSE)
  } else character()
  result <- if (file.exists(result_path)) {
    jsonlite::fromJSON(result_path, simplifyVector = FALSE)
  } else {
    list(
      valid = FALSE,
      error = paste0("worker_result_missing;process_status=", status)
    )
  }
  list(
    process_status = as.integer(status),
    rss_backend = backend,
    peak_rss_bytes = bench_parse_peak_rss(time_lines),
    result = result,
    artifacts = list(
      spec = spec_path,
      result = result_path,
      stdout = stdout_path,
      time = time_path
    )
  )
}

bench_validate_moment_fugw_design <- function() {
  contract <- bench_read_moment_fugw_thresholds()
  scaling <- contract$scaling
  checks <- c(
    four_full_sizes = length(scaling$domain_sizes) == 4L,
    reaches_5000 = max(scaling$domain_sizes) >= 5000,
    memory_exponent = identical(
      as.numeric(scaling$memory_doubling_exponent_max), 1.5
    ),
    runtime_envelope = identical(
      as.numeric(c(
        scaling$fixed_work_runtime_exponent_min,
        scaling$fixed_work_runtime_exponent_max
      )), c(1.8, 2.2)
    ),
    runtime_fit = identical(
      as.character(scaling$runtime_exponent_primary_fit),
      "largest_three_points"
    ) && isTRUE(scaling$runtime_exponent_all_points_required),
    rank_curve = identical(
      as.integer(scaling$structure_ranks), c(4L, 16L, 64L, 256L)
    ),
    heldout_curve = identical(
      as.integer(scaling$heldout_batches), c(1L, 10L, 100L, 1000L)
    ),
    r_visible_allocation =
      isTRUE(scaling$require_r_visible_allocation) &&
      setequal(
        as.character(scaling$r_visible_allocation_actions),
        c("solve", "objective", "apply", "adjoint")
      ),
    separate_convergence_curve = all(vapply(
      contract$profiles,
      function(profile) length(profile$converged_sizes) > 0L,
      logical(1)
    ))
  )
  list(passed = all(checks), checks = checks)
}

.bench_moment_base_spec <- function(config, package_mode, package_root) {
  list(
    package_mode = package_mode,
    package_root = normalizePath(package_root, mustWork = TRUE),
    protocol_path = bench_moment_fugw_resource("protocol.R"),
    benchmark_path = bench_moment_fugw_resource(
      "benchmark_moment_fugw.R"
    ),
    threads = as.integer(config$threads),
    seed = as.integer(config$seed),
    feature_dimension = as.integer(config$feature_dimension),
    epsilon = as.numeric(config$epsilon),
    reg_marginals = as.numeric(config$reg_marginals),
    tol = as.numeric(config$tol),
    tol_ot = as.numeric(config$tol_ot),
    block_size = as.integer(config$block_size)
  )
}

.bench_scalar <- function(x, default = NA) {
  if (is.null(x) || !length(x)) return(default)
  x[[1L]]
}

.bench_worker_row <- function(
    run, curve, phase, repetition, baseline_rss, profile) {
  result <- run$result %||% list()
  dense <- result$dense_sentinel %||% list()
  peak <- as.numeric(run$peak_rss_bytes %||% NA_real_)
  baseline_rss <- as.numeric(baseline_rss)
  increment <- if (is.finite(peak) && is.finite(baseline_rss)) {
    max(0, peak - baseline_rss)
  } else NA_real_
  raw_seconds <- as.numeric(result$operation_seconds %||% NA_real_)
  reportable <- isTRUE(result$timing_eligible) ||
    (identical(curve, "fixed_work_scaling") &&
      isTRUE(result$scaling_eligible)) ||
    (identical(curve, "operator_throughput") &&
      isTRUE(result$operator_eligible))
  data.frame(
    profile = profile,
    curve = curve,
    phase = phase,
    repetition = as.integer(repetition),
    action = as.character(result$action %||% NA_character_),
    method = as.character(result$method %||% "fugw_factorized"),
    evidence_class = as.character(result$evidence_class %||% NA_character_),
    n_source = as.integer(result$n_source %||% NA_integer_),
    n_target = as.integer(result$n_target %||% NA_integer_),
    structure_rank = as.integer(result$structure_rank %||% NA_integer_),
    feature_dimension = as.integer(
      result$feature_dimension %||% NA_integer_
    ),
    batch = as.integer(result$batch %||% NA_integer_),
    status = as.character(
      result$status %||% result$source_status %||% NA_character_
    ),
    process_status = as.integer(run$process_status %||% NA_integer_),
    infrastructure_error = isTRUE(result$infrastructure_error) ||
      !identical(as.integer(run$process_status %||% 1L), 0L),
    valid = isTRUE(result$valid),
    certified = isTRUE(result$certified) || isTRUE(result$source_certified),
    comparison_eligible = isTRUE(result$comparison_eligible),
    timing_eligible = isTRUE(result$timing_eligible),
    scaling_eligible = isTRUE(result$scaling_eligible),
    operator_eligible = isTRUE(result$operator_eligible),
    memory_contract_eligible = isTRUE(result$memory_contract_eligible) ||
      (identical(result$action, "solve") && isTRUE(result$scaling_eligible)),
    threshold_update_eligible = isTRUE(result$threshold_update_eligible),
    contract_status = as.character(result$contract_status %||% NA_character_),
    reject_reason = as.character(
      result$reject_reason %||% result$error %||% ""
    ),
    limitation_reason = as.character(result$limitation_reason %||% ""),
    iterations = as.integer(result$iterations %||% NA_integer_),
    inner_iterations = as.integer(result$inner_iterations %||% NA_integer_),
    residual = as.numeric(result$residual %||% NA_real_),
    setup_seconds_raw = as.numeric(result$setup_seconds %||% NA_real_),
    operation_seconds_raw = raw_seconds,
    end_to_end_seconds_raw = as.numeric(
      result$end_to_end_seconds %||% NA_real_
    ),
    reportable_seconds = if (reportable) raw_seconds else NA_real_,
    peak_rss_bytes = peak,
    baseline_rss_bytes = baseline_rss,
    incremental_peak_rss_bytes = increment,
    rss_backend = as.character(run$rss_backend %||% NA_character_),
    r_visible_allocation_bytes = as.numeric(
      result$r_visible_allocation_bytes %||% NA_real_
    ),
    dense_sentinel_passed = isTRUE(dense$passed),
    full_block_calls = as.integer(dense$full_block_calls %||% NA_integer_),
    max_tile_elements = as.numeric(dense$max_tile_elements %||% NA_real_),
    full_matrix_elements = as.numeric(
      dense$full_matrix_elements %||% NA_real_
    ),
    spec_artifact = as.character(run$artifacts$spec %||% NA_character_),
    result_artifact = as.character(run$artifacts$result %||% NA_character_),
    time_artifact = as.character(run$artifacts$time %||% NA_character_),
    stringsAsFactors = FALSE
  )
}

.bench_curve_medians <- function(rows, curve, action, value) {
  keep <- rows$curve == curve & rows$action == action &
    rows$phase == "measured" & rows$valid & is.finite(rows[[value]])
  selected <- rows[keep, , drop = FALSE]
  if (!nrow(selected)) {
    return(data.frame(n_source = integer(), value = numeric()))
  }
  aggregate(
    selected[[value]], list(n_source = selected$n_source), stats::median,
    na.rm = TRUE
  ) |>
    stats::setNames(c("n_source", "value"))
}

.bench_moment_fugw_release_hardware <- function(profile, package_mode, meta) {
  identical(profile, "full") &&
    identical(package_mode, "installed") &&
    !nzchar(as.character(meta$rfugw_fast_flags %||% "")[[1L]]) &&
    isTRUE(meta$git_provenance_complete) &&
    identical(meta$git_dirty, FALSE) &&
    length(meta$commit) == 1L &&
    !is.na(meta$commit[[1L]]) &&
    grepl("^[[:xdigit:]]{40}$", meta$commit[[1L]])
}

bench_profile_r_visible_allocation <- function(operation) {
  if (!is.function(operation)) {
    stop("`operation` must be a function.", call. = FALSE)
  }
  if (!isTRUE(capabilities("profmem"))) {
    stop("This R build does not support Rprofmem allocation profiling.",
         call. = FALSE)
  }
  path <- tempfile("moment-fugw-rprofmem-")
  profiling <- FALSE
  on.exit({
    if (profiling) utils::Rprofmem(NULL)
    unlink(path)
  }, add = TRUE)
  utils::Rprofmem(path, append = FALSE, threshold = 0)
  profiling <- TRUE
  value <- operation()
  utils::Rprofmem(NULL)
  profiling <- FALSE
  lines <- readLines(path, warn = FALSE)
  bytes <- suppressWarnings(as.numeric(sub("[[:space:]].*$", "", lines)))
  bytes <- bytes[is.finite(bytes)]
  total <- sum(bytes)
  if (!length(bytes) || !is.finite(total) || total <= 0) {
    stop("Rprofmem did not record a positive allocation total.", call. = FALSE)
  }
  list(value = value, bytes = total, events = length(bytes))
}

bench_run_moment_fugw_protocol <- function(
    profile = c("smoke", "pr", "full"),
    results_root = file.path("inst", "bench", "results", "moment_fugw"),
    package_mode = c("source", "installed"),
    package_root = ".") {
  profile <- match.arg(profile)
  package_mode <- match.arg(package_mode)
  bench_validate_moment_fugw_threshold_history()
  design <- bench_validate_moment_fugw_design()
  if (!isTRUE(design$passed)) {
    stop("Moment-FUGW benchmark design contract is invalid.", call. = FALSE)
  }
  config <- bench_moment_fugw_profile(profile)
  bench_pin_threads(config$threads)
  bench_archive_current(results_root)
  out_dir <- file.path(results_root, "current")
  raw_dir <- file.path(out_dir, "raw")
  fit_dir <- file.path(out_dir, "fits")
  dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(fit_dir, recursive = TRUE, showWarnings = FALSE)
  base_spec <- .bench_moment_base_spec(config, package_mode, package_root)
  baseline_spec <- utils::modifyList(base_spec, list(action = "baseline"))
  baseline <- bench_run_fresh_process(
    baseline_spec, raw_dir, "baseline"
  )
  if (!isTRUE(baseline$result$valid) ||
      !is.finite(baseline$peak_rss_bytes)) {
    stop("Fresh-process RSS baseline failed.", call. = FALSE)
  }
  baseline_rss <- baseline$peak_rss_bytes
  rows <- list()
  add_run <- function(run, curve, phase, repetition) {
    rows[[length(rows) + 1L]] <<- .bench_worker_row(
      run, curve, phase, repetition, baseline_rss, profile
    )
    invisible(run)
  }
  run_solve <- function(
      n, rank, max_iter, max_iter_ot, evidence_class, label, fit_path = "",
      worker_action = "solve") {
    spec <- utils::modifyList(base_spec, list(
      action = worker_action,
      n = as.integer(n),
      structure_rank = as.integer(rank),
      max_iter = as.integer(max_iter),
      max_iter_ot = as.integer(max_iter_ot),
      evidence_class = evidence_class,
      fit_path = fit_path
    ))
    bench_run_fresh_process(spec, raw_dir, label)
  }

  fixed_fits <- list()
  for (n in as.integer(config$domain_sizes)) {
    total <- as.integer(config$warmup + config$repetitions)
    for (index in seq_len(total)) {
      measured <- index > config$warmup
      repetition <- if (measured) index - config$warmup else index
      fit_path <- if (measured && repetition == 1L) {
        file.path(fit_dir, paste0("fixed-n", n, ".rds"))
      } else ""
      run <- run_solve(
        n, config$structure_rank, 1L, 1L, "fixed_work_scaling",
        paste0("fixed-n", n, "-", if (measured) "rep" else "warmup", repetition),
        fit_path
      )
      add_run(
        run, "fixed_work_scaling",
        if (measured) "measured" else "warmup", repetition
      )
      if (nzchar(fit_path) && file.exists(fit_path)) {
        fixed_fits[[as.character(n)]] <- fit_path
      }
    }
    allocation <- run_solve(
      n, config$structure_rank, 1L, 1L, "fixed_work_scaling",
      paste0("allocation-solve-n", n), "", worker_action = "allocation"
    )
    add_run(allocation, "r_visible_allocation", "measured", 1L)
    fit_path <- fixed_fits[[as.character(n)]] %||% ""
    if (nzchar(fit_path)) {
      for (action in c("objective", "apply", "adjoint")) {
        spec <- utils::modifyList(base_spec, list(
          action = action,
          fit_path = fit_path,
          batch = 1L,
          source_evidence_class = "fixed_work_scaling",
          measure_r_visible = TRUE
        ))
        run <- bench_run_fresh_process(
          spec, raw_dir, paste0("memory-", action, "-n", n)
        )
        add_run(run, "matrix_free_memory", "measured", 1L)
      }
    }
  }

  for (rank in as.integer(config$structure_ranks)) {
    run <- run_solve(
      config$rank_size, rank, 1L, 1L, "fixed_work_scaling",
      paste0("rank-", rank), ""
    )
    add_run(run, "rank_scaling", "measured", 1L)
  }

  certified_fit <- ""
  for (n in as.integer(config$converged_sizes)) {
    total <- as.integer(config$warmup + config$repetitions)
    for (index in seq_len(total)) {
      measured <- index > config$warmup
      repetition <- if (measured) index - config$warmup else index
      fit_path <- if (measured && repetition == 1L &&
          n == as.integer(config$operator_size)) {
        file.path(fit_dir, paste0("certified-n", n, ".rds"))
      } else ""
      run <- run_solve(
        n, config$structure_rank, config$max_iter, config$max_iter_ot,
        "certified_comparison",
        paste0("converged-n", n, "-", if (measured) "rep" else "warmup", repetition),
        fit_path
      )
      add_run(
        run, "end_to_end_convergence",
        if (measured) "measured" else "warmup", repetition
      )
      if (nzchar(fit_path) && file.exists(fit_path) &&
          isTRUE(run$result$certified)) certified_fit <- fit_path
    }
  }

  if (nzchar(certified_fit)) {
    for (batch in as.integer(config$heldout_batches)) {
      for (action in c("apply", "adjoint")) {
        spec <- utils::modifyList(base_spec, list(
          action = action,
          fit_path = certified_fit,
          batch = batch,
          source_evidence_class = "operator_throughput"
        ))
        run <- bench_run_fresh_process(
          spec, raw_dir, paste0(action, "-batch", batch)
        )
        add_run(run, "operator_throughput", "measured", 1L)
      }
    }
  }

  rows <- do.call(rbind, rows)
  fixed_time <- .bench_curve_medians(
    rows, "fixed_work_scaling", "solve", "reportable_seconds"
  )
  fixed_memory <- .bench_curve_medians(
    rows, "fixed_work_scaling", "solve", "incremental_peak_rss_bytes"
  )
  fixed_memory_absolute <- .bench_curve_medians(
    rows, "fixed_work_scaling", "solve", "peak_rss_bytes"
  )
  runtime_exponent_all <- bench_loglog_exponent(
    fixed_time$n_source, fixed_time$value
  )
  runtime_fit_count <- min(3L, nrow(fixed_time))
  runtime_exponent <- if (runtime_fit_count >= 2L) {
    bench_loglog_exponent(
      tail(fixed_time$n_source, runtime_fit_count),
      tail(fixed_time$value, runtime_fit_count)
    )
  } else NA_real_
  memory_exponent <- bench_loglog_exponent(
    fixed_memory$n_source, fixed_memory$value
  )
  absolute_memory_exponent <- bench_loglog_exponent(
    fixed_memory_absolute$n_source, fixed_memory_absolute$value
  )
  contract <- bench_read_moment_fugw_thresholds()
  full_sizes <- as.integer(contract$scaling$domain_sizes)
  measured_fixed <- rows[
    rows$curve == "fixed_work_scaling" & rows$phase == "measured", ,
    drop = FALSE
  ]
  completed_sizes <- sort(unique(
    measured_fixed$n_source[measured_fixed$valid]
  ))
  memory_rows <- rows[rows$curve %in% c(
    "fixed_work_scaling", "matrix_free_memory", "r_visible_allocation"
  ), , drop = FALSE]
  dense_contract_pass <- nrow(memory_rows) > 0L &&
    all(memory_rows$valid) && all(memory_rows$dense_sentinel_passed) &&
    all(memory_rows$full_block_calls == 0L)
  endpoint_rows <- rows[
    rows$curve == "end_to_end_convergence" & rows$phase == "measured" &
      rows$certified & rows$timing_eligible, , drop = FALSE
  ]
  largest_certified <- if (nrow(endpoint_rows)) {
    max(endpoint_rows$n_source)
  } else NA_integer_
  full_profile_complete <- identical(profile, "full") &&
    identical(completed_sizes, full_sizes)
  rank_complete <- setequal(
    rows$structure_rank[rows$curve == "rank_scaling" & rows$valid],
    as.integer(contract$scaling$structure_ranks)
  )
  operator_rows <- rows[
    rows$curve == "operator_throughput" & rows$valid, , drop = FALSE
  ]
  batch_complete <- setequal(
    operator_rows$batch, as.integer(contract$scaling$heldout_batches)
  ) && setequal(operator_rows$action, c("apply", "adjoint"))
  allocation_rows <- rows[
    rows$curve %in% c("r_visible_allocation", "matrix_free_memory") &
      rows$valid & is.finite(rows$r_visible_allocation_bytes) &
      rows$r_visible_allocation_bytes > 0,
    , drop = FALSE
  ]
  expected_allocations <- as.vector(outer(
    as.integer(config$domain_sizes),
    c("solve", "objective", "apply", "adjoint"),
    function(n, action) paste(n, action, sep = ":")
  ))
  observed_allocations <- paste(
    allocation_rows$n_source, allocation_rows$action, sep = ":"
  )
  r_visible_allocation_complete <- setequal(
    observed_allocations, expected_allocations
  )
  summary <- list(
    schema_version = 1L,
    profile = profile,
    package_mode = package_mode,
    admission_evaluable = FALSE,
    baseline_rss_bytes = baseline_rss,
    rss_backend = baseline$rss_backend,
    fixed_work_runtime_exponent = runtime_exponent,
    fixed_work_runtime_exponent_fit = "largest_three_points",
    fixed_work_runtime_exponent_all_points = runtime_exponent_all,
    expected_runtime_exponent = c(
      contract$scaling$fixed_work_runtime_exponent_min,
      contract$scaling$fixed_work_runtime_exponent_max
    ),
    incremental_memory_exponent = memory_exponent,
    absolute_peak_memory_exponent = absolute_memory_exponent,
    memory_exponent_max = contract$scaling$memory_doubling_exponent_max,
    completed_fixed_work_sizes = completed_sizes,
    configured_full_sizes = full_sizes,
    largest_completed_certified_endpoint = largest_certified,
    claim_ceiling_n = largest_certified,
    gates = list(
      design = design$passed,
      full_profile_complete = full_profile_complete,
      exact_commit_provenance = FALSE,
      clean_exact_source_tree = FALSE,
      release_hardware = FALSE,
      fixed_work_runtime_envelope = if (full_profile_complete) {
        is.finite(runtime_exponent) &&
          runtime_exponent >= contract$scaling$fixed_work_runtime_exponent_min &&
          runtime_exponent <= contract$scaling$fixed_work_runtime_exponent_max
      } else NA,
      near_linear_peak_memory = if (full_profile_complete) {
        is.finite(memory_exponent) &&
          memory_exponent < contract$scaling$memory_doubling_exponent_max
      } else NA,
      prohibited_dense_allocations_absent = dense_contract_pass,
      r_visible_allocation_complete = r_visible_allocation_complete,
      rank_curve_complete = rank_complete,
      operator_batch_curve_complete = batch_complete,
      larger_certified_endpoint = is.finite(largest_certified) &&
        largest_certified >= contract$problem$minimum_certified_endpoint
    ),
    release_baseline_admitted = FALSE,
    release_baseline_reason = paste0(
      "Candidate thresholds remain fail-closed until the full installed ",
      "conservative artifact and independent review are recorded in history."
    )
  )
  meta <- bench_capture_env(
    seed = config$seed,
    threads = config$threads,
    warmup = config$warmup,
    reps = config$repetitions,
    profile = paste("moment_fugw", profile, package_mode, sep = "_"),
    repository_root = package_root
  )
  history <- jsonlite::fromJSON(
    bench_moment_fugw_threshold_history_path(), simplifyVector = FALSE
  )
  meta$benchmark_schema_version <- 1L
  meta$package_mode <- package_mode
  meta$thresholds_md5 <- unname(tools::md5sum(
    bench_moment_fugw_thresholds_path()
  )[[1L]])
  meta$threshold_history_status <-
    history$entries[[length(history$entries)]]$admission_status
  meta$rss_backend <- baseline$rss_backend
  release_hardware_run <- .bench_moment_fugw_release_hardware(
    profile, package_mode, meta
  )
  summary$gates$exact_commit_provenance <-
    isTRUE(meta$git_provenance_complete)
  summary$gates$clean_exact_source_tree <-
    isTRUE(meta$git_provenance_complete) && identical(meta$git_dirty, FALSE)
  summary$gates$release_hardware <- release_hardware_run
  required_gates <- c(
    "design", "full_profile_complete", "exact_commit_provenance",
    "clean_exact_source_tree", "release_hardware",
    "fixed_work_runtime_envelope", "near_linear_peak_memory",
    "prohibited_dense_allocations_absent", "r_visible_allocation_complete",
    "rank_curve_complete", "operator_batch_curve_complete",
    "larger_certified_endpoint"
  )
  summary$admission_evaluable <- all(vapply(
    summary$gates[required_gates], isTRUE, logical(1)
  ))
  meta$r_visible_allocation_semantics <- paste0(
    "Required separate Rprofmem curve; never substituted for native peak ",
    "RSS. Fresh-process /usr/bin/time remains authoritative for peak memory."
  )
  utils::write.csv(rows, file.path(out_dir, "runs.csv"), row.names = FALSE)
  jsonlite::write_json(
    summary, file.path(out_dir, "summary.json"), auto_unbox = TRUE,
    pretty = TRUE, null = "null", na = "null"
  )
  meta$raw_artifacts <- sort(list.files(raw_dir, recursive = TRUE))
  jsonlite::write_json(
    meta, file.path(out_dir, "meta.json"), auto_unbox = TRUE,
    pretty = TRUE, null = "null", na = "null"
  )
  invisible(list(
    out_dir = out_dir,
    meta = meta,
    rows = rows,
    summary = summary
  ))
}
