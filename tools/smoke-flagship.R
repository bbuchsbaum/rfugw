#!/usr/bin/env Rscript
# Certificate-bearing solver smoke used by sanitizer and serial-build CI.

library(rfugw)
if (!requireNamespace("jsonlite", quietly = TRUE)) {
  stop("Sanitizer smoke requires `jsonlite` for its evidence ledger.")
}

evidence_dir <- file.path(".gate", "numerical-trust")
dir.create(evidence_dir, recursive = TRUE, showWarnings = FALSE)
records <- list()

record_result <- function(name, result) {
  value <- rfugw_value(result)
  if (length(value) != 1L || !is.finite(value)) {
    stop(name, " did not return one finite certified value.", call. = FALSE)
  }
  if (!identical(result$status, "converged") || !isTRUE(result$converged)) {
    stop(name, " did not converge under its sanitizer smoke fixture.",
         call. = FALSE)
  }
  truth_fields <- c(
    "feasible", "objective_consistent", "inner_converged",
    "regularized_dual_consistent", "kkt_consistent",
    "stationarity_consistent", "nonnegative_within_tolerance"
  )
  present <- intersect(truth_fields, names(result))
  present <- present[vapply(result[present], function(x) {
    length(x) == 1L && !is.na(x)
  }, logical(1))]
  failed <- present[!vapply(result[present], isTRUE, logical(1))]
  if (length(failed)) {
    stop(
      name, " returned a failed certificate: ", paste(failed, collapse = ", "),
      call. = FALSE
    )
  }
  numeric_fields <- c(
    "residual", "row_residual", "col_residual", "feasibility_residual",
    "objective_residual", "regularized_duality_gap", "duality_gap",
    "primal_dual_gap", "kkt_residual", "fixed_point_residual",
    "max_inner_residual", "frank_wolfe_gap"
  )
  observed <- intersect(numeric_fields, names(result))
  records[[name]] <<- list(
    status = result$status,
    converged = result$converged,
    value = unname(value),
    certificates = as.list(result[present]),
    residuals = as.list(result[observed])
  )
  invisible(result)
}

M <- matrix(c(0, 0.2, 0.2, 0), 2L)
C <- M
p <- c(0.5, 0.5)

sink <- record_result("balanced_sinkhorn", ot_sinkhorn(
  M, p, p, epsilon = 0.2, method = "log", max_iter = 5000L, tol = 1e-9
))
stopifnot(isTRUE(sink$regularized_dual_consistent))

emd <- record_result("exact_emd", ot_emd(M, p, p))
stopifnot(abs(emd$duality_gap) <= emd$duality_gap_tolerance)

divergence <- record_result("sinkhorn_divergence", ot_sinkhorn_divergence(
  c(0, 1), c(0.2, 1.2), epsilon = 0.5, method = "log",
  max_iter = 5000L, tol = 1e-9
))
stopifnot(all(divergence$component_converged))

record_result("exact_partial_ot", ot_partial_emd(M, p, p, mass = 0.7))
record_result("entropic_partial_ot", ot_partial_sinkhorn(
  M, p, p, mass = 0.7, epsilon = 0.3, method = "log",
  max_iter = 10000L, tol = 1e-8
))
record_result("penalized_partial_ot", ot_partial_penalized(
  M, p, p, discard_penalty = 0.3
))

record_result("kl_unbalanced_ot", ot_sinkhorn_unbalanced(
  M, p, p, epsilon = 0.2, rho = 2, method = "log",
  max_iter = 5000L, tol = 1e-8
))
ti <- record_result("translation_invariant_kl_uot", ot_sinkhorn_unbalanced_ti(
  matrix(0.2, 1L), p = 2, q = 3, epsilon = 0.3, rho = c(1, 2),
  max_iter = 10000L, tol = 1e-10, plan = "dense"
))
stopifnot(abs(ti$primal_dual_gap) <= ti$primal_dual_gap_tolerance)

record_result("wasserstein_summary", ot_wasserstein_distance(
  c(0, 1), c(0.2, 1.2), p = 2
))
barycenter <- record_result("fixed_support_barycenter", ot_barycenter_weights(
  list(M, matrix(c(0.2, 0.8, 0.8, 0.2), 2L)),
  list(p, c(0.3, 0.7)), coefficients = c(0.4, 0.6),
  mode = "exact", tol = 1e-10
))
stopifnot(all(barycenter$component_status == "converged"))

record_result("entropic_fgw", fgw_entropic(
  M, C, C, p, p, epsilon = 0.2, max_iter = 100L, tol = 1e-7,
  sinkhorn_tol = 1e-8, precision = "strict_double"
))
record_result("exact_fgw", fgw_exact_cg(
  M, C, C, p, p, max_iter = 100L, tol_rel = 1e-8, tol_abs = 1e-8
))
record_result("fugw", fugw_kl(
  C, C, wx = p, wy = p, M = M, epsilon = 0.1, reg_marginals = 1,
  max_iter = 100L, tol = 1e-6, max_iter_ot = 5000L, tol_ot = 1e-8
))

manifest <- if (file.exists("release-artifact.json")) {
  jsonlite::read_json("release-artifact.json", simplifyVector = TRUE)
} else {
  NULL
}
ledger <- list(
  schema_version = 1L,
  package = "rfugw",
  scope = "release_sanitizer_smoke",
  source_commit = Sys.getenv("GITHUB_SHA", NA_character_),
  artifact_manifest = manifest,
  sanitizer = list(
    cxxflags = Sys.getenv("RFUGW_EXTRA_CXXFLAGS", NA_character_),
    libraries = Sys.getenv("RFUGW_EXTRA_LIBS", NA_character_),
    asan_options = Sys.getenv("ASAN_OPTIONS", NA_character_),
    ubsan_options = Sys.getenv("UBSAN_OPTIONS", NA_character_)
  ),
  results = records
)
ledger_path <- file.path(evidence_dir, "sanitizer-smoke.json")
jsonlite::write_json(
  ledger, ledger_path, auto_unbox = TRUE, pretty = TRUE, na = "null"
)
cat(sprintf(
  "Certificate-bearing sanitizer smoke passed for %d solver families; wrote %s.\n",
  length(records), ledger_path
))
