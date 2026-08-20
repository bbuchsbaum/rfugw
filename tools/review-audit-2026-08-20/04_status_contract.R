audit_lib <- Sys.getenv("RFUGW_AUDIT_LIB", unset = "")
if (nzchar(audit_lib)) {
  .libPaths(c(audit_lib, .libPaths()))
}

library(rfugw)

# Directly exercise the shared diagnostic helper. All values below are finite
# except the deliberately invalid objective and inner residual; the plan also
# contains negative entries and violates both marginals.
synthetic <- rfugw:::.attach_solver_diagnostics(
  out = list(
    plan = matrix(c(-1, 2, 2, -1), 2, 2),
    ot_dist = Inf,
    objective = Inf
  ),
  residual = 0,
  converged = TRUE,
  iterations = 1L,
  max_iter = 10L,
  p = c(0.5, 0.5),
  q = c(0.5, 0.5),
  inner_residual = Inf
)

cat("synthetic_status:", synthetic$status, "\n")
cat("synthetic_min_plan:", min(synthetic$plan), "\n")
cat("synthetic_row_residual:", synthetic$row_residual, "\n")
cat("synthetic_col_residual:", synthetic$col_residual, "\n")
cat("synthetic_inner_residual:", synthetic$inner_residual, "\n")
cat("synthetic_ot_dist:", synthetic$ot_dist, "\n")
cat("synthetic_objective:", synthetic$objective, "\n")

# A public FGW run can satisfy the outer plan-update tolerance while its final
# one-iteration Sinkhorn projection misses the requested marginal tolerance.
set.seed(8)
n <- 6L
C1 <- as.matrix(dist(matrix(rnorm(n * 2), n, 2))); C1 <- C1 / max(C1)
C2 <- as.matrix(dist(matrix(rnorm(n * 2), n, 2))); C2 <- C2 / max(C2)
M <- matrix(runif(n * n), n, n)
fgw_one_inner <- fgw_entropic(
  M, C1, C2,
  epsilon = 0.03,
  max_iter = 500L,
  tol = 1e-7,
  sinkhorn_max_iter = 1L,
  sinkhorn_tol = 1e-12,
  precision = "double",
  symmetric = TRUE
)

cat("fgw_status:", fgw_one_inner$status, "\n")
cat("fgw_outer_residual:", format(fgw_one_inner$residual, digits = 16), "\n")
cat("fgw_row_residual:", format(fgw_one_inner$row_residual, digits = 16), "\n")
cat("fgw_col_residual:", format(fgw_one_inner$col_residual, digits = 16), "\n")
cat("fgw_inner_residual:", fgw_one_inner$inner_residual, "\n")

# FUGW reports outer convergence but does not return an inner residual. A
# one-iteration inner budget and a well-solved budget lead to materially
# different plans and objectives while both are labelled converged.
set.seed(1)
nx <- 5L
ny <- 4L
Cx <- as.matrix(dist(matrix(rnorm(nx * 2), nx, 2))); Cx <- Cx / max(Cx)
Cy <- as.matrix(dist(matrix(rnorm(ny * 2), ny, 2))); Cy <- Cy / max(Cy)
M_fugw <- matrix(runif(nx * ny), nx, ny)
fugw_args <- list(
  Cx = Cx,
  Cy = Cy,
  epsilon = 0.1,
  reg_marginals = 1,
  alpha = 0.5,
  M = M_fugw,
  max_iter = 100L,
  tol = 1e-7,
  tol_ot = 1e-12,
  precision = "double"
)
fugw_one_inner <- do.call(fugw_kl, c(fugw_args, list(max_iter_ot = 1L)))
fugw_full_inner <- do.call(fugw_kl, c(fugw_args, list(max_iter_ot = 1000L)))

cat("fugw_one_inner_status:", fugw_one_inner$status, "\n")
cat("fugw_one_inner_residual:", fugw_one_inner$inner_residual, "\n")
cat("fugw_full_inner_status:", fugw_full_inner$status, "\n")
cat(
  "fugw_objective_gap:",
  format(abs(fugw_one_inner$fugw_cost - fugw_full_inner$fugw_cost), digits = 16),
  "\n"
)
cat(
  "fugw_plan_max_error:",
  format(max(abs(fugw_one_inner$pi_samp - fugw_full_inner$pi_samp)), digits = 16),
  "\n"
)

# The entropic partial adapter exposes neither status nor the final inner
# residual. Tight and one-step inner budgets yield materially different plans.
set.seed(5)
ns <- 5L
nt <- 4L
C1p <- as.matrix(dist(matrix(rnorm(ns * 2), ns, 2))); C1p <- C1p / max(C1p)
C2p <- as.matrix(dist(matrix(rnorm(nt * 2), nt, 2))); C2p <- C2p / max(C2p)
Mp <- matrix(runif(ns * nt), ns, nt)
partial_args <- list(
  M = Mp,
  C1 = C1p,
  C2 = C2p,
  p = rep(1 / ns, ns),
  q = rep(1 / nt, nt),
  reg = 0.1,
  m = 0.7,
  alpha = 0.5,
  numItermax = 50L,
  tol = 1e-7,
  symmetric = TRUE,
  log = TRUE,
  check_every = 1L,
  inner_tol = 1e-12
)
partial_one_inner <- do.call(
  entropic_partial_fused_gromov_wasserstein,
  c(partial_args, list(inner_max_iter = 1L))
)
partial_full_inner <- do.call(
  entropic_partial_fused_gromov_wasserstein,
  c(partial_args, list(inner_max_iter = 500L))
)

cat("partial_return_fields:", paste(names(partial_one_inner), collapse = ","), "\n")
cat(
  "partial_objective_gap:",
  format(abs(partial_one_inner$partial_fgw_dist - partial_full_inner$partial_fgw_dist), digits = 16),
  "\n"
)
cat(
  "partial_plan_max_error:",
  format(max(abs(partial_one_inner$plan - partial_full_inner$plan)), digits = 16),
  "\n"
)

stopifnot(
  identical(synthetic$status, "converged"),
  min(synthetic$plan) < 0,
  synthetic$row_residual > 0,
  !is.finite(synthetic$inner_residual),
  !is.finite(synthetic$ot_dist),
  !is.finite(synthetic$objective),
  identical(fgw_one_inner$status, "converged"),
  fgw_one_inner$row_residual > 1e-6,
  fgw_one_inner$row_residual > 1e4 * 1e-12,
  is.na(fgw_one_inner$inner_residual),
  identical(fugw_one_inner$status, "converged"),
  is.na(fugw_one_inner$inner_residual),
  abs(fugw_one_inner$fugw_cost - fugw_full_inner$fugw_cost) > 1e-4,
  !"status" %in% names(partial_one_inner),
  !"inner_residual" %in% names(partial_one_inner),
  max(abs(partial_one_inner$plan - partial_full_inner$plan)) > 1e-3
)

sessionInfo()
