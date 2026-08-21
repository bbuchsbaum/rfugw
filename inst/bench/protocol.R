# Shared correctness-controlled benchmark helpers.
# See inst/bench/PROTOCOL.md.

bench_thresholds_path <- function() {
  candidates <- c(
    "inst/bench/thresholds.json",
    file.path(system.file(package = "rfugw"), "bench", "thresholds.json")
  )
  hit <- candidates[file.exists(candidates)]
  if (!length(hit)) {
    stop("Could not find inst/bench/thresholds.json.", call. = FALSE)
  }
  hit[[1]]
}

bench_parse_suites <- function(x) {
  if (length(x) == 0L || (length(x) == 1L && !nzchar(x[[1]]))) {
    x <- Sys.getenv("RFUGW_PROTOCOL_SUITES", unset = "")
  }
  if (length(x) == 1L) {
    x <- strsplit(as.character(x), ",", fixed = TRUE)[[1]]
  }
  suites <- unique(trimws(as.character(x)))
  suites <- suites[nzchar(suites)]
  known <- c("linear", "fgw", "fugw", "semirelaxed", "partial", "ucoot", "sampled")
  if (!length(suites)) {
    return(known)
  }
  bad <- setdiff(suites, known)
  if (length(bad)) {
    stop("Unknown protocol suite(s): ", paste(bad, collapse = ", "), call. = FALSE)
  }
  suites
}

bench_read_thresholds <- function(method = NULL) {
  all <- jsonlite::fromJSON(bench_thresholds_path(), simplifyVector = TRUE)
  if (is.null(method)) {
    return(all)
  }
  if (!method %in% names(all)) {
    stop("Unknown benchmark method: ", method, call. = FALSE)
  }
  all[[method]]
}

bench_threshold_history_path <- function() {
  candidates <- c(
    "inst/bench/threshold-history.json",
    file.path(system.file(package = "rfugw"), "bench", "threshold-history.json")
  )
  hit <- candidates[file.exists(candidates)]
  if (!length(hit)) {
    stop("Could not find inst/bench/threshold-history.json.", call. = FALSE)
  }
  hit[[1]]
}

bench_validate_threshold_history <- function() {
  history <- jsonlite::fromJSON(bench_threshold_history_path(), simplifyVector = FALSE)
  entries <- history$entries %||% list()
  if (!length(entries)) {
    stop("Threshold history has no retained evidence entry.", call. = FALSE)
  }
  current <- unname(tools::md5sum(bench_thresholds_path())[[1]])
  recorded <- entries[[length(entries)]]$thresholds_md5 %||% ""
  if (!identical(current, recorded)) {
    stop(
      "thresholds.json changed without a matching reviewed evidence entry in threshold-history.json.",
      call. = FALSE
    )
  }
  if (!nzchar(entries[[length(entries)]]$evidence %||% "") ||
      !nzchar(entries[[length(entries)]]$review_requirement %||% "")) {
    stop("Latest threshold history entry lacks evidence or review requirements.", call. = FALSE)
  }
  invisible(TRUE)
}

bench_close <- function(actual, expected, atol, rtol) {
  actual <- as.numeric(actual)
  expected <- as.numeric(expected)
  all(is.finite(actual), is.finite(expected),
      abs(actual - expected) <= atol + rtol * abs(expected))
}

bench_find_description <- function() {
  if (file.exists("DESCRIPTION")) {
    return("DESCRIPTION")
  }
  installed <- system.file("DESCRIPTION", package = "rfugw")
  if (nzchar(installed) && file.exists(installed)) {
    return(installed)
  }
  for (up in c("..", file.path("..", ".."), file.path("..", "..", ".."))) {
    cand <- file.path(up, "DESCRIPTION")
    if (file.exists(cand)) {
      return(cand)
    }
  }
  stop("DESCRIPTION not found.", call. = FALSE)
}

bench_capture_env <- function(seed, threads, warmup, reps, profile = "conservative") {
  desc <- read.dcf(bench_find_description())
  commit <- tryCatch(
    system2("git", c("rev-parse", "HEAD"), stdout = TRUE, stderr = FALSE)[1],
    error = function(e) NA_character_
  )
  si <- as.list(Sys.info())
  runtime_provenance <- tryCatch(
    rfugw:::.solver_runtime_provenance(),
    error = function(e) list(error = conditionMessage(e))
  )
  r_config <- function(key) {
    tryCatch(
      system2(file.path(R.home("bin"), "R"), c("CMD", "config", key),
              stdout = TRUE, stderr = FALSE)[1],
      error = function(e) NA_character_
    )
  }
  list(
    package = unname(desc[1, "Package"]),
    version = unname(desc[1, "Version"]),
    commit = commit,
    r_version = paste(R.version$major, R.version$minor, sep = "."),
    r_platform = R.version$platform,
    compiler_cxx17 = r_config("CXX17"),
    compiler_cxx17std = r_config("CXX17STD"),
    blas = extSoftVersion()[["BLAS"]] %||% NA_character_,
    profile = profile,
    rfugw_fast_flags = Sys.getenv("RFUGW_FAST_FLAGS", ""),
    rfugw_extra_cxxflags = Sys.getenv("RFUGW_EXTRA_CXXFLAGS", ""),
    rfugw_openmp_flags = Sys.getenv("RFUGW_OPENMP_FLAGS", ""),
    omp_num_threads = Sys.getenv("OMP_NUM_THREADS", ""),
    openblas_num_threads = Sys.getenv("OPENBLAS_NUM_THREADS", ""),
    seed = as.integer(seed),
    threads = as.integer(threads),
    warmup = as.integer(warmup),
    reps = as.integer(reps),
    runtime_provenance = runtime_provenance,
    sysname = si$sysname,
    release = si$release,
    machine = si$machine,
    nodename = si$nodename,
    timestamp = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
  )
}

bench_pin_threads <- function(threads) {
  threads <- as.integer(threads)
  if (!is.finite(threads) || threads < 1L) threads <- 1L
  Sys.setenv(
    OMP_NUM_THREADS = as.character(threads),
    OPENBLAS_NUM_THREADS = "1",
    MKL_NUM_THREADS = "1",
    VECLIB_MAXIMUM_THREADS = "1"
  )
  threads
}

bench_make_problem <- function(kind = c("linear", "fgw", "ucoot"), n, seed) {
  kind <- match.arg(kind)
  set.seed(as.integer(seed) + as.integer(n))
  if (identical(kind, "ucoot")) {
    nf <- max(3L, as.integer(round(n * 0.6)))
    X <- matrix(rnorm(n * nf), n, nf)
    Y <- matrix(rnorm(n * nf), n, nf)
    return(list(kind = kind, n = n, X = X, Y = Y,
                wx_samp = rep(1 / n, n), wy_samp = rep(1 / n, n),
                wx_feat = rep(1 / nf, nf), wy_feat = rep(1 / nf, nf)))
  }
  if (identical(kind, "linear")) {
    X <- matrix(rnorm(n * 2L), n, 2L)
    Y <- matrix(rnorm(n * 2L), n, 2L)
    M <- as.matrix(stats::dist(rbind(X, Y)))[seq_len(n), n + seq_len(n)]
    M <- M / max(M)
    return(list(kind = kind, n = n, M = M, X = X, Y = Y,
                p = rep(1 / n, n), q = rep(1 / n, n)))
  }
  X <- matrix(rnorm(n * 2L), n, 2L)
  Y <- matrix(rnorm(n * 2L), n, 2L)
  F1 <- matrix(rnorm(n * 2L), n, 2L)
  F2 <- matrix(rnorm(n * 2L), n, 2L)
  C1 <- as.matrix(stats::dist(X)); C1 <- C1 / max(C1)
  C2 <- as.matrix(stats::dist(Y)); C2 <- C2 / max(C2)
  M <- as.matrix(stats::dist(rbind(F1, F2)))[seq_len(n), n + seq_len(n)]
  M <- M / max(M)
  list(kind = kind, n = n, M = M, C1 = C1, C2 = C2,
       p = rep(1 / n, n), q = rep(1 / n, n))
}

bench_check_quality <- function(
    method,
    result,
    problem,
    evidence_class = c("certified_comparison", "fixed_budget_performance")) {
  evidence_class <- match.arg(evidence_class)
  method_spec <- bench_read_thresholds(method)
  certified <- identical(evidence_class, "certified_comparison")
  spec <- if (certified) method_spec$quality else method_spec$performance_regression
  if (is.null(spec)) {
    return(list(
      valid = FALSE,
      certified = FALSE,
      comparison_eligible = FALSE,
      performance_regression_eligible = FALSE,
      evidence_class = evidence_class,
      reject_reason = paste0("no_", evidence_class, "_contract")
    ))
  }
  reasons <- character()
  status <- result$status
  if (is.null(status) || !length(status) || !nzchar(as.character(status)[1])) {
    reasons <- c(reasons, "missing_status")
  } else if (is.null(spec$status) || !status %in% spec$status) {
    reasons <- c(reasons, sprintf("status=%s", status))
  }
  if (certified) {
    if (!identical(status, "converged")) {
      reasons <- c(reasons, "not_converged")
    }
    required <- c(
      converged = result$converged,
      feasible = result$feasible,
      objective_consistent = result$objective_consistent,
      objective_components_consistent = result$objective_components_consistent
    )
    failed <- names(required)[vapply(required, function(x) !isTRUE(x), logical(1))]
    if (length(failed)) {
      reasons <- c(reasons, paste0("certificate_", failed))
    }
    if (!is.null(result$inner_converged) &&
        length(result$inner_converged) &&
        !is.na(result$inner_converged[[1]]) &&
        !isTRUE(result$inner_converged)) {
      reasons <- c(reasons, "certificate_inner_converged")
    }
  }
  residual <- result$residual %||% result$error %||% NA_real_
  residual <- as.numeric(residual)[1]
  if (!is.null(spec$residual_max) &&
      (!is.finite(residual) || residual > spec$residual_max)) {
    reasons <- c(reasons, sprintf("residual=%s", residual))
  }
  value <- tryCatch(rfugw_value(result), error = function(e) NA_real_)
  if (!is.finite(value) && !is.null(result$gw_dist_estimated)) {
    value <- as.numeric(result$gw_dist_estimated)[1]
  }
  if (!is.finite(value)) {
    reasons <- c(reasons, "nonfinite_objective")
  }
  if (method %in% c("ot_sinkhorn", "ot_emd")) {
    recon <- ot_linear_cost(problem$M, result)
    if (!bench_close(recon, value, spec$objective_atol, spec$objective_rtol)) {
      reasons <- c(reasons, "objective_mismatch")
    }
    ok_plan <- tryCatch({
      ot_validate_plan(result, problem$p, problem$q)
      TRUE
    }, error = function(e) FALSE)
    if (!ok_plan) reasons <- c(reasons, "invalid_plan")
  }
  if (identical(method, "ot_sinkhorn_unbalanced")) {
    recon <- ot_linear_cost(problem$M, result)
    if (!bench_close(recon, value, spec$objective_atol, spec$objective_rtol)) {
      reasons <- c(reasons, "objective_mismatch")
    }
    mass <- result$mass %||% sum(result$plan)
    if (!is.null(spec$mass_max) && is.finite(mass) && mass > spec$mass_max + 1e-8) {
      reasons <- c(reasons, sprintf("mass=%s", mass))
    }
  }
  if (identical(method, "fgw_entropic")) {
    recon <- ot_fgw_square(problem$M, problem$C1, problem$C2, result, alpha = 0.5)
    if (!bench_close(recon, value, spec$objective_atol, spec$objective_rtol)) {
      reasons <- c(reasons, "objective_mismatch")
    }
  }
  if (identical(method, "fugw_kl")) {
    if (!is.finite(value)) {
      reasons <- c(reasons, "nonfinite_objective")
    }
    if (is.null(result$pi_samp) || any(!is.finite(result$pi_samp))) {
      reasons <- c(reasons, "invalid_plan")
    }
  }
  if (identical(method, "sampled_gromov_wasserstein")) {
    plan <- result$plan %||% result
    if (!is.matrix(plan) || any(!is.finite(plan))) {
      reasons <- c(reasons, "invalid_plan")
    } else {
      if (max(abs(rowSums(plan) - problem$p)) > 1e-3 ||
          max(abs(colSums(plan) - problem$q)) > 1e-3) {
        reasons <- c(reasons, "invalid_plan")
      }
    }
    est <- result$gw_dist_estimated %||% NA_real_
    if (!is.finite(est)) {
      reasons <- c(reasons, "nonfinite_objective")
    }
  }
  if (method %in% c("semirelaxed_gromov_wasserstein",
                    "entropic_semirelaxed_gromov_wasserstein")) {
    recon <- ot_gw_square(problem$C1, problem$C2, result)
    if (!bench_close(recon, value, spec$objective_atol, spec$objective_rtol)) {
      reasons <- c(reasons, "objective_mismatch")
    }
    if (max(abs(rowSums(result$plan) - problem$p)) > 1e-8) {
      reasons <- c(reasons, "invalid_plan")
    }
  }
  if (method %in% c("partial_gromov_wasserstein",
                    "entropic_partial_gromov_wasserstein")) {
    recon <- ot_gw_square(problem$C1, problem$C2, result)
    if (!bench_close(recon, value, spec$objective_atol, spec$objective_rtol)) {
      reasons <- c(reasons, "objective_mismatch")
    }
    mass <- spec$mass %||% 0.7
    ok_plan <- tryCatch({
      ot_validate_plan(result, problem$p, problem$q, mass = mass, marginals = "partial")
      TRUE
    }, error = function(e) FALSE)
    if (!ok_plan) reasons <- c(reasons, "invalid_plan")
  }
  if (identical(method, "partial_fused_gromov_wasserstein")) {
    alpha <- spec$alpha %||% 0.5
    recon <- ot_fgw_square(problem$M, problem$C1, problem$C2, result, alpha = alpha)
    if (!bench_close(recon, value, spec$objective_atol, spec$objective_rtol)) {
      reasons <- c(reasons, "objective_mismatch")
    }
    mass <- spec$mass %||% 0.7
    ok_plan <- tryCatch({
      ot_validate_plan(result, problem$p, problem$q, mass = mass, marginals = "partial")
      TRUE
    }, error = function(e) FALSE)
    if (!ok_plan) reasons <- c(reasons, "invalid_plan")
  }
  if (method %in% c("unbalanced_co_optimal_transport",
                    "fused_unbalanced_across_spaces_divergence")) {
    if (!is.finite(value)) {
      reasons <- c(reasons, "nonfinite_objective")
    }
    if (is.null(result$pi_samp) || is.null(result$pi_feat) ||
        any(!is.finite(result$pi_samp)) || any(!is.finite(result$pi_feat))) {
      reasons <- c(reasons, "invalid_plan")
    }
  }
  if (method %in% c("semirelaxed_fused_gromov_wasserstein",
                    "entropic_semirelaxed_fused_gromov_wasserstein")) {
    alpha <- spec$alpha %||% 0.5
    recon <- ot_fgw_square(problem$M, problem$C1, problem$C2, result, alpha = alpha)
    if (!bench_close(recon, value, spec$objective_atol, spec$objective_rtol)) {
      reasons <- c(reasons, "objective_mismatch")
    }
    if (max(abs(rowSums(result$plan) - problem$p)) > 1e-8) {
      reasons <- c(reasons, "invalid_plan")
    }
  }
  valid <- !length(reasons)
  list(
    valid = valid,
    certified = certified && valid,
    comparison_eligible = certified && valid,
    performance_regression_eligible = !certified && valid,
    evidence_class = evidence_class,
    reject_reason = paste(unique(reasons), collapse = ";")
  )
}

`%||%` <- function(x, y) if (is.null(x)) y else x

bench_time_ms <- function(fn) {
  t0 <- proc.time()[["elapsed"]]
  val <- fn()
  elapsed <- (proc.time()[["elapsed"]] - t0) * 1000
  list(value = val, ms = elapsed)
}

bench_run_split <- function(prepare_fn, solve_fn, warmup = 1L, reps = 3L,
                            method, problem_for_quality = NULL,
                            evidence_class = "certified_comparison") {
  for (i in seq_len(warmup)) {
    dat <- prepare_fn()
    out <- solve_fn(dat)
    if (!is.null(method)) {
      q <- bench_check_quality(
        method, out, problem_for_quality %||% dat,
        evidence_class = evidence_class
      )
      if (!q$valid) {
        return(list(
          valid = FALSE,
          certified = FALSE,
          comparison_eligible = FALSE,
          performance_regression_eligible = FALSE,
          evidence_class = evidence_class,
          reject_reason = paste0("warmup:", q$reject_reason),
          prepare_ms = NA_real_,
          solve_ms = NA_real_,
          e2e_ms = NA_real_,
          mem_solve_bytes = NA_real_,
          mem_e2e_bytes = NA_real_,
          result = out
        ))
      }
    }
  }
  prep <- numeric(reps)
  solv <- numeric(reps)
  last <- NULL
  for (i in seq_len(reps)) {
    p <- bench_time_ms(prepare_fn)
    s <- bench_time_ms(function() solve_fn(p$value))
    prep[[i]] <- p$ms
    solv[[i]] <- s$ms
    last <- s$value
  }
  q <- if (!is.null(method)) {
    bench_check_quality(
      method, last, problem_for_quality %||% prepare_fn(),
      evidence_class = evidence_class
    )
  } else {
    list(
      valid = TRUE, certified = FALSE, comparison_eligible = FALSE,
      performance_regression_eligible = FALSE,
      evidence_class = "uncertified", reject_reason = ""
    )
  }
  mem_solve <- NA_real_
  mem_e2e <- NA_real_
  if (requireNamespace("bench", quietly = TRUE)) {
    b_solve <- bench::mark(solve_fn(prepare_fn()), iterations = 1L, memory = TRUE, check = FALSE)
    mem_solve <- as.numeric(b_solve$mem_alloc[[1]])
    b_e2e <- bench::mark({
      d <- prepare_fn()
      solve_fn(d)
    }, iterations = 1L, memory = TRUE, check = FALSE)
    mem_e2e <- as.numeric(b_e2e$mem_alloc[[1]])
  }
  list(
    valid = q$valid,
    certified = q$certified,
    comparison_eligible = q$comparison_eligible,
    performance_regression_eligible = q$performance_regression_eligible,
    evidence_class = q$evidence_class,
    reject_reason = q$reject_reason,
    prepare_ms = stats::median(prep),
    solve_ms = stats::median(solv),
    e2e_ms = stats::median(prep + solv),
    mem_solve_bytes = mem_solve,
    mem_e2e_bytes = mem_e2e,
    result = last
  )
}

bench_archive_current <- function(results_root = "inst/bench/results") {
  current <- file.path(results_root, "current")
  if (!dir.exists(current)) {
    return(invisible(FALSE))
  }
  scratch <- file.path(results_root, "scratch")
  dir.create(scratch, recursive = TRUE, showWarnings = FALSE)
  stamp <- format(Sys.time(), "%Y%m%dT%H%M%SZ", tz = "UTC")
  dest <- file.path(scratch, stamp)
  dir.create(dest, recursive = TRUE, showWarnings = FALSE)
  files <- list.files(current, full.names = TRUE)
  if (length(files)) {
    file.copy(files, dest, recursive = TRUE)
  }
  unlink(current, recursive = TRUE)
  invisible(TRUE)
}

bench_write_baseline <- function(meta, runs, quality, out_dir = "inst/bench/results/current") {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  jsonlite::write_json(meta, file.path(out_dir, "meta.json"), auto_unbox = TRUE, pretty = TRUE)
  utils::write.csv(runs, file.path(out_dir, "runs.csv"), row.names = FALSE)
  utils::write.csv(quality, file.path(out_dir, "quality.csv"), row.names = FALSE)
  invisible(out_dir)
}
