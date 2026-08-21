audit_lib <- Sys.getenv("RFUGW_AUDIT_LIB", unset = "")
if (nzchar(audit_lib)) {
  .libPaths(c(audit_lib, .libPaths()))
}

library(rfugw)

generalized_kl <- function(plan, p, q) {
  ref <- tcrossprod(p, q)
  if (any(plan > 0 & ref == 0)) {
    return(Inf)
  }
  z <- plan > 0
  sum(plan[z] * log(plan[z] / ref[z])) - sum(plan) + sum(ref)
}

plan_unbalanced <- matrix(c(0.30, 0.10, 0.00, 0.20), 2, 2)
p <- c(0.4, 0.6)
q <- c(0.7, 0.3)
reported <- ot_kl(plan_unbalanced, p, q)
expected <- generalized_kl(plan_unbalanced, p, q)
missing_mass_correction <- -sum(plan_unbalanced) + sum(tcrossprod(p, q))

cat("unbalanced_plan_mass:", sum(plan_unbalanced), "\n")
cat("reference_mass:", sum(tcrossprod(p, q)), "\n")
cat("reported_ot_kl:", format(reported, digits = 16), "\n")
cat("generalized_kl:", format(expected, digits = 16), "\n")
cat("missing_mass_correction:", format(missing_mass_correction, digits = 16), "\n")
cat("observed_difference:", format(expected - reported, digits = 16), "\n")

p_zero <- c(1, 0)
q_zero <- c(0.5, 0.5)
plan_outside_support <- matrix(c(0.4, 0.1, 0.4, 0.1), 2, 2)
reported_zero_support <- ot_kl(plan_outside_support, p_zero, q_zero)
expected_zero_support <- generalized_kl(plan_outside_support, p_zero, q_zero)

cat("zero_support_reported_finite:", is.finite(reported_zero_support), "\n")
cat("zero_support_reported:", format(reported_zero_support, digits = 16), "\n")
cat("zero_support_expected:", expected_zero_support, "\n")

invalid_negative <- suppressWarnings(ot_kl(plan_unbalanced, c(1.2, -0.2), q))
invalid_nonfinite <- suppressWarnings(ot_kl(plan_unbalanced, c(NA_real_, 1), q))
cat("negative_reference_result:", invalid_negative, "\n")
cat("nonfinite_reference_result:", invalid_nonfinite, "\n")

# Existing balanced tests cannot detect the missing mass terms because both
# the plan and product reference have unit mass.
balanced_plan <- diag(c(0.5, 0.5))
balanced_reported <- ot_kl(balanced_plan, c(0.5, 0.5), c(0.5, 0.5))
balanced_expected <- generalized_kl(balanced_plan, c(0.5, 0.5), c(0.5, 0.5))
cat("balanced_difference:", balanced_expected - balanced_reported, "\n")

stopifnot(
  abs((expected - reported) - missing_mass_correction) < 1e-12,
  missing_mass_correction > 0,
  is.finite(reported_zero_support),
  is.infinite(expected_zero_support),
  !isTRUE(all.equal(invalid_negative, expected)),
  is.na(invalid_nonfinite),
  abs(balanced_expected - balanced_reported) < 1e-12
)

sessionInfo()
