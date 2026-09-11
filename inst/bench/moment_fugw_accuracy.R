suppressPackageStartupMessages({
  rlib <- Sys.getenv("RFUGW_RLIB", unset = "")
  if (nzchar(rlib)) .libPaths(c(rlib, .libPaths()))
  if (identical(Sys.getenv("RFUGW_DEV_SOURCE", unset = "0"), "1")) {
    if (!requireNamespace("pkgload", quietly = TRUE)) {
      stop("RFUGW_DEV_SOURCE=1 requires pkgload.", call. = FALSE)
    }
    pkgload::load_all(".", quiet = TRUE)
  } else {
    library(rfugw)
  }
})

args <- commandArgs(trailingOnly = TRUE)
out_csv <- if (length(args) >= 1L) {
  args[[1L]]
} else {
  "inst/bench/moment-fugw-accuracy-evidence.csv"
}
profile <- if (length(args) >= 2L) match.arg(args[[2L]], c("quick", "full")) else "full"
base_seed <- if (length(args) >= 3L) as.integer(args[[3L]]) else 20260827L
if (!is.finite(base_seed)) stop("The base seed must be finite.", call. = FALSE)

Sys.setenv(
  OMP_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1"
)

oracle_candidates <- c(
  "inst/bench/moment_fugw_oracle.R",
  file.path(system.file(package = "rfugw"), "bench", "moment_fugw_oracle.R")
)
oracle_path <- oracle_candidates[file.exists(oracle_candidates)][[1L]]
source(oracle_path, local = TRUE)

relative_scalar_error <- function(actual, expected) {
  abs(actual - expected) / max(1e-12, abs(expected))
}

relative_matrix_error <- function(actual, expected) {
  sqrt(sum((actual - expected)^2)) /
    max(1e-12, sqrt(sum(expected^2)))
}

make_square_cost <- function(n, rank, seed) {
  set.seed(seed)
  factorized_cost(
    matrix(runif(n * rank, -.18, .18), n, rank),
    matrix(runif(n * rank, -.18, .18), n, rank),
    row = runif(n, .2, .45),
    column = runif(n, .2, .45),
    exact = TRUE,
    provenance = list(fixture_seed = seed, fixture_rank = rank)
  )
}

make_cross_cost <- function(ns, nt, rank, seed) {
  set.seed(seed)
  factorized_cost(
    matrix(runif(ns * rank, -.2, .2), ns, rank),
    matrix(runif(nt * rank, -.2, .2), nt, rank),
    row = runif(ns, .1, .3),
    column = runif(nt, .1, .3),
    exact = TRUE,
    provenance = list(fixture_seed = seed, fixture_rank = rank)
  )
}

design_candidates <- c(
  "inst/bench/moment-fugw-accuracy-design.csv",
  file.path(
    system.file(package = "rfugw"), "bench",
    "moment-fugw-accuracy-design.csv"
  )
)
design_candidates <- design_candidates[file.exists(design_candidates)]
if (!length(design_candidates)) {
  stop("Could not find the frozen Moment-FUGW accuracy design.", call. = FALSE)
}
case_table <- utils::read.csv(
  design_candidates[[1L]], stringsAsFactors = FALSE
)
design_columns <- c(
  "case_id", "seed_offset", "ns", "nt", "source_rank", "target_rank",
  "feature_rank", "block_size", "zero_weight", "near_zero_weight",
  "epsilon", "rho_source", "rho_target", "feature_weight",
  "structure_weight"
)
if (!identical(names(case_table), design_columns) || nrow(case_table) != 12L ||
    anyDuplicated(case_table$case_id) ||
    !identical(case_table$case_id, sprintf("mfugw-%02d", 1:12)) ||
    !is.logical(case_table$zero_weight) ||
    !is.logical(case_table$near_zero_weight) ||
    any(case_table$zero_weight & case_table$near_zero_weight) ||
    any(!is.finite(as.matrix(case_table[setdiff(
      design_columns, c("case_id", "zero_weight", "near_zero_weight")
    )]))) ||
    any(case_table[c(
      "ns", "nt", "source_rank", "target_rank", "feature_rank", "block_size"
    )] <= 0)) {
  stop("The frozen Moment-FUGW accuracy design is invalid.", call. = FALSE)
}
if (identical(profile, "quick")) case_table <- case_table[1:4, , drop = FALSE]

dense_moments <- function(plan, source_right, target_right) {
  source_mass <- rowSums(plan)
  target_mass <- colSums(plan)
  list(
    source_mass = source_mass,
    target_mass = target_mass,
    H = crossprod(source_right, plan %*% target_right),
    G_source = crossprod(source_right, source_right * source_mass),
    G_target = crossprod(target_right, target_right * target_mass),
    transported_mass = sum(plan)
  )
}

maximum_moment_error <- function(actual, expected) {
  max(
    relative_matrix_error(actual$source_mass, expected$source_mass),
    relative_matrix_error(actual$target_mass, expected$target_mass),
    relative_matrix_error(actual$H, expected$H),
    relative_matrix_error(actual$G_source, expected$G_source),
    relative_matrix_error(actual$G_target, expected$G_target),
    relative_scalar_error(actual$transported_mass, expected$transported_mass)
  )
}

failure_row <- function(spec, seed, message, dense_seconds = NA_real_,
                        implicit_seconds = NA_real_) {
  data.frame(
    spec, seed = seed, profile = profile,
    dense_status = "error", implicit_status = "error",
    dense_converged = FALSE, implicit_converged = FALSE,
    inner_certified = FALSE, outer_certified = FALSE,
    objective_relative_error = Inf,
    component_relative_error = Inf,
    direct_oracle_relative_error = Inf,
    action_relative_error = Inf,
    adjoint_relative_error = Inf,
    mass_relative_error = Inf,
    moment_relative_error = Inf,
    dense_seconds = dense_seconds,
    implicit_seconds = implicit_seconds,
    pass = FALSE,
    reject_reason = message,
    stringsAsFactors = FALSE
  )
}

run_case <- function(spec) {
  seed <- base_seed + spec$seed_offset
  Cx <- make_square_cost(spec$ns, spec$source_rank, seed + 1L)
  Cy <- make_square_cost(spec$nt, spec$target_rank, seed + 2L)
  M <- make_cross_cost(spec$ns, spec$nt, spec$feature_rank, seed + 3L)
  set.seed(seed + 4L)
  wx <- runif(spec$ns)
  wy <- runif(spec$nt)
  if (isTRUE(spec$zero_weight)) {
    wx[[1L]] <- 0
    wy[[spec$nt]] <- 0
  }
  if (isTRUE(spec$near_zero_weight)) {
    wx[[1L]] <- 1e-12
    wy[[spec$nt]] <- 1e-13
  }
  wx <- wx / sum(wx)
  wy <- wy / sum(wy)
  controls <- list(
    wx = wx,
    wy = wy,
    reg_marginals = c(spec$rho_source, spec$rho_target),
    epsilon = spec$epsilon,
    feature_weight = spec$feature_weight,
    structure_weight = spec$structure_weight,
    max_iter = 100L,
    tol = 5e-9,
    max_iter_ot = 10000L,
    tol_ot = 1e-11,
    rescale_plan = TRUE
  )

  dense_seconds <- NA_real_
  implicit_seconds <- NA_real_
  tryCatch({
    dense_time <- system.time({
      dense <- do.call(fugw_kl, c(list(
        Cx = cost_block(Cx), Cy = cost_block(Cy), M = cost_block(M),
        precision = "strict_double"
      ), controls))
    })
    dense_seconds <- unname(dense_time[["elapsed"]])
    implicit_time <- system.time({
      implicit <- do.call(fugw_factorized, c(list(
        Cx = Cx, Cy = Cy, M = M, block_size = spec$block_size
      ), controls))
    })
    implicit_seconds <- unname(implicit_time[["elapsed"]])
    sample <- transport_plan_materialize(implicit$plans$sample)
    feature <- transport_plan_materialize(implicit$plans$feature)

    objective_relative_error <- relative_scalar_error(
      implicit$fugw_cost, dense$fugw_cost
    )
    component_relative_error <- max(vapply(
      c("structure_unweighted", "feature_unweighted", "regularization", "total"),
      function(component) relative_scalar_error(
        implicit$objective_decomposition[[component]],
        dense$objective_decomposition[[component]]
      ),
      numeric(1)
    ))

    direct_oracle_relative_error <- NA_real_
    if (max(spec$ns, spec$nt) <= 8L) {
      direct <- moment_fugw_oracle(
        cost_block(Cx), cost_block(Cy), cost_block(M), wx, wy,
        sample, feature, controls$reg_marginals, controls$epsilon,
        controls$feature_weight, controls$structure_weight
      )
      direct_oracle_relative_error <- max(
        relative_scalar_error(implicit$fugw_cost, direct$fugw_cost),
        relative_scalar_error(
          implicit$objective_decomposition$structure_unweighted,
          direct$structure_unweighted
        ),
        relative_scalar_error(
          implicit$objective_decomposition$feature_unweighted,
          direct$feature_unweighted
        ),
        relative_scalar_error(
          implicit$objective_decomposition$regularization,
          direct$regularization
        )
      )
    }

    set.seed(seed + 5L)
    target_values <- matrix(rnorm(spec$nt * 3L), spec$nt, 3L)
    source_values <- matrix(rnorm(spec$ns * 3L), spec$ns, 3L)
    action_relative_error <- max(
      relative_matrix_error(
        transport_plan_apply(implicit$plans$sample, target_values),
        dense$pi_samp %*% target_values
      ),
      relative_matrix_error(
        transport_plan_apply(implicit$plans$feature, target_values),
        dense$pi_feat %*% target_values
      )
    )
    adjoint_relative_error <- max(
      relative_matrix_error(
        transport_plan_adjoint(implicit$plans$sample, source_values),
        t(dense$pi_samp) %*% source_values
      ),
      relative_matrix_error(
        transport_plan_adjoint(implicit$plans$feature, source_values),
        t(dense$pi_feat) %*% source_values
      )
    )
    mass_relative_error <- max(
      relative_scalar_error(sum(sample), sum(dense$pi_samp)),
      relative_scalar_error(sum(feature), sum(dense$pi_feat))
    )
    source_right <- rfugw:::.cost_factor_pair(Cx)$right
    target_right <- rfugw:::.cost_factor_pair(Cy)$right
    moment_relative_error <- max(
      maximum_moment_error(
        implicit$moments$sample,
        dense_moments(dense$pi_samp, source_right, target_right)
      ),
      maximum_moment_error(
        implicit$moments$feature,
        dense_moments(dense$pi_feat, source_right, target_right)
      )
    )

    dense_certified <- isTRUE(dense$converged) &&
      isTRUE(dense$inner_converged) &&
      isTRUE(dense$objective_consistent) &&
      isTRUE(dense$objective_components_consistent)
    implicit_certified <- isTRUE(implicit$converged) &&
      isTRUE(implicit$inner_uot_certified) &&
      isTRUE(implicit$certificate$outer_stationarity$certified) &&
      isTRUE(implicit$objective_consistent) &&
      isTRUE(implicit$objective_components_consistent)
    direct_ok <- is.na(direct_oracle_relative_error) ||
      direct_oracle_relative_error <= 1e-8
    passed <- dense_certified && implicit_certified &&
      objective_relative_error <= 1e-8 &&
      component_relative_error <= 1e-8 && direct_ok &&
      action_relative_error <= 1e-6 && adjoint_relative_error <= 1e-6 &&
      mass_relative_error <= 1e-7 && moment_relative_error <= 1e-7
    reasons <- c(
      if (!dense_certified) "dense_uncertified",
      if (!implicit_certified) "implicit_uncertified",
      if (objective_relative_error > 1e-8) "objective_error",
      if (component_relative_error > 1e-8) "component_error",
      if (!direct_ok) "direct_oracle_error",
      if (action_relative_error > 1e-6) "action_error",
      if (adjoint_relative_error > 1e-6) "adjoint_error",
      if (mass_relative_error > 1e-7) "mass_error",
      if (moment_relative_error > 1e-7) "moment_error"
    )

    data.frame(
      spec, seed = seed, profile = profile,
      dense_status = dense$status,
      implicit_status = implicit$status,
      dense_converged = isTRUE(dense$converged),
      implicit_converged = isTRUE(implicit$converged),
      inner_certified = isTRUE(implicit$inner_uot_certified),
      outer_certified = isTRUE(implicit$certificate$outer_stationarity$certified),
      objective_relative_error = objective_relative_error,
      component_relative_error = component_relative_error,
      direct_oracle_relative_error = direct_oracle_relative_error,
      action_relative_error = action_relative_error,
      adjoint_relative_error = adjoint_relative_error,
      mass_relative_error = mass_relative_error,
      moment_relative_error = moment_relative_error,
      dense_seconds = dense_seconds,
      implicit_seconds = implicit_seconds,
      pass = passed,
      reject_reason = paste(reasons, collapse = ";"),
      stringsAsFactors = FALSE
    )
  }, error = function(error) {
    failure_row(
      spec, seed, conditionMessage(error), dense_seconds, implicit_seconds
    )
  })
}

rows <- lapply(seq_len(nrow(case_table)), function(index) {
  result <- run_case(case_table[index, , drop = FALSE])
  cat(sprintf(
    "%s: %s (dense %.3fs, implicit %.3fs)\n",
    result$case_id, if (isTRUE(result$pass)) "PASS" else "FAIL",
    result$dense_seconds, result$implicit_seconds
  ))
  result
})
results <- do.call(rbind, rows)
dir.create(dirname(out_csv), recursive = TRUE, showWarnings = FALSE)
utils::write.csv(results, out_csv, row.names = FALSE)
cat("Wrote Moment-FUGW accuracy matrix:", out_csv, "\n")
print(results[c(
  "case_id", "ns", "nt", "dense_status", "implicit_status",
  "objective_relative_error", "action_relative_error",
  "moment_relative_error", "pass", "reject_reason"
)], row.names = FALSE)

if (any(!results$pass)) {
  stop("Moment-FUGW accuracy matrix contains failed or uncertified rows.",
       call. = FALSE)
}
