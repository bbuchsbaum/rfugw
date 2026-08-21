audit_lib <- Sys.getenv("RFUGW_AUDIT_LIB", unset = "")
if (nzchar(audit_lib)) {
  .libPaths(c(audit_lib, .libPaths()))
}

library(rfugw)

# Fractional counts are silently truncated.
truncated_count <- rfugw:::.validate_count(2.9, "max_iter")
cat("validate_count_2.9:", truncated_count, "\n")

# Partial validation accepts any total mass when `mass` is omitted.
p <- c(0.5, 0.5)
q <- c(0.5, 0.5)
zero_partial <- matrix(0, 2, 2)
partial_without_mass_accepted <- tryCatch({
  ot_validate_plan(zero_partial, p, q, marginals = "partial")
  TRUE
}, error = function(e) FALSE)
cat("partial_without_mass_accepted:", partial_without_mass_accepted, "\n")
cat("partial_without_mass_total:", sum(zero_partial), "\n")

# In relaxed mode supplied reference vectors are length-checked only.
relaxed_invalid_refs_accepted <- tryCatch({
  ot_validate_plan(
    diag(c(0.5, 0.5)),
    p = c(NA_real_, -1),
    q = c(Inf, 0),
    marginals = "relaxed"
  )
  TRUE
}, error = function(e) FALSE)
cat("relaxed_invalid_refs_accepted:", relaxed_invalid_refs_accepted, "\n")

# Symmetry uses a fixed absolute 1e-10 threshold with no relative component.
C_large <- matrix(c(0, 1e12, 1e12 + 1e-3, 0), 2, 2)
relative_asymmetry <- max(abs(C_large - t(C_large))) / max(abs(C_large))
absolute_symmetry_rejected <- inherits(
  try(rfugw:::.resolve_symmetric(TRUE, C_large, C_large), silent = TRUE),
  "try-error"
)
cat("large_cost_relative_asymmetry:", format(relative_asymmetry, digits = 16), "\n")
cat("large_cost_symmetric_true_rejected:", absolute_symmetry_rejected, "\n")

# A global additive cost shift should not change balanced entropic OT. The
# scaling kernel clips every shifted exponent and loses the cost contrast;
# the log-domain plan remains shift-invariant to numerical tolerance.
M <- matrix(c(0, 1, 1, 0), 2, 2)
shift <- 1e6
sink_args <- list(p = p, q = q, epsilon = 0.01, max_iter = 2000L, tol = 1e-12)
scaling_base <- do.call(ot_sinkhorn, c(list(M = M, method = "scaling"), sink_args))
scaling_shift <- do.call(ot_sinkhorn, c(list(M = M + shift, method = "scaling"), sink_args))
log_base <- do.call(ot_sinkhorn, c(list(M = M, method = "log"), sink_args))
log_shift <- do.call(ot_sinkhorn, c(list(M = M + shift, method = "log"), sink_args))
scaling_shift_error <- max(abs(scaling_base$plan - scaling_shift$plan))
log_shift_error <- max(abs(log_base$plan - log_shift$plan))
cat("scaling_additive_shift_plan_error:", format(scaling_shift_error, digits = 16), "\n")
cat("log_additive_shift_plan_error:", format(log_shift_error, digits = 16), "\n")
cat("scaling_shift_status:", scaling_shift$status, "\n")
cat("log_shift_status:", log_shift$status, "\n")

# The public unbalanced linear OT primitive has no log/scaling method switch.
cat("uot_has_method_argument:", "method" %in% names(formals(ot_sinkhorn_unbalanced)), "\n")

# UCOOT/across-space `sinkhorn_log` is currently an alias for scaling.
set.seed(2)
X <- matrix(rnorm(12), 4, 3)
Y <- matrix(rnorm(15), 5, 3)
ucoot_args <- list(
  X = X,
  Y = Y,
  reg_marginals = 2,
  epsilon = 0.1,
  max_iter = 10L,
  tol = 1e-7,
  max_iter_ot = 20L,
  tol_ot = 1e-8,
  log = TRUE
)
ucoot_scaling <- do.call(
  fused_unbalanced_across_spaces_divergence,
  c(ucoot_args, list(unbalanced_solver = "sinkhorn"))
)
ucoot_log_alias <- do.call(
  fused_unbalanced_across_spaces_divergence,
  c(ucoot_args, list(unbalanced_solver = "sinkhorn_log"))
)
cat("ucoot_sinkhorn_log_sample_identical:", identical(ucoot_scaling$pi_samp, ucoot_log_alias$pi_samp), "\n")
cat("ucoot_sinkhorn_log_feature_identical:", identical(ucoot_scaling$pi_feat, ucoot_log_alias$pi_feat), "\n")
cat("ucoot_sinkhorn_log_cost_identical:", identical(ucoot_scaling$ucoot_cost, ucoot_log_alias$ucoot_cost), "\n")

stopifnot(
  truncated_count == 2L,
  partial_without_mass_accepted,
  relaxed_invalid_refs_accepted,
  relative_asymmetry < 1e-12,
  absolute_symmetry_rejected,
  scaling_shift_error > 0.1,
  log_shift_error < 1e-7,
  !"method" %in% names(formals(ot_sinkhorn_unbalanced)),
  identical(ucoot_scaling$pi_samp, ucoot_log_alias$pi_samp),
  identical(ucoot_scaling$pi_feat, ucoot_log_alias$pi_feat),
  identical(ucoot_scaling$ucoot_cost, ucoot_log_alias$ucoot_cost)
)

sessionInfo()
