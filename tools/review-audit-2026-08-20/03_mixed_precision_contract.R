audit_lib <- Sys.getenv("RFUGW_AUDIT_LIB", unset = "")
if (nzchar(audit_lib)) {
  .libPaths(c(audit_lib, .libPaths()))
}

library(rfugw)

make_problem <- function(n, seed) {
  set.seed(seed)
  C1 <- as.matrix(dist(matrix(rnorm(n * 3), n, 3)))
  C2 <- as.matrix(dist(matrix(rnorm(n * 3), n, 3)))
  list(
    C1 = C1 / max(C1),
    C2 = C2 / max(C2),
    M = matrix(runif(n * n), n, n)
  )
}

solve_problem <- function(d, precision, tol, sinkhorn_tol, max_iter = 500L) {
  fgw_entropic(
    d$M, d$C1, d$C2,
    epsilon = 0.1,
    max_iter = max_iter,
    tol = tol,
    sinkhorn_max_iter = 1000L,
    sinkhorn_tol = sinkhorn_tol,
    precision = precision,
    symmetric = TRUE
  )
}

small <- make_problem(8L, 1L)
tight_mixed <- solve_problem(small, "mixed", 1e-9, 1e-9)
floor_mixed <- solve_problem(small, "mixed", 1e-6, 1e-6)
strict_double <- solve_problem(small, "double", 1e-9, 1e-9)

default_precision <- fgw_entropic(
  small$M, small$C1, small$C2,
  epsilon = 0.1,
  max_iter = 500L,
  tol = 1e-9,
  sinkhorn_tol = 1e-9,
  symmetric = TRUE
)

cat("tight_mixed_status:", tight_mixed$status, "\n")
cat("tight_mixed_iterations:", tight_mixed$iterations, "of", tight_mixed$max_iter, "\n")
cat("tight_mixed_residual:", format(tight_mixed$residual, digits = 16), "\n")
cat("strict_double_status:", strict_double$status, "\n")
cat("strict_double_iterations:", strict_double$iterations, "\n")
cat("strict_double_residual:", format(strict_double$residual, digits = 16), "\n")
cat(
  "tight_vs_floor_plan_max_error:",
  format(max(abs(tight_mixed$plan - floor_mixed$plan)), digits = 16),
  "\n"
)
cat("default_equals_explicit_mixed:", identical(default_precision$plan, tight_mixed$plan), "\n")

provenance_fields <- c(
  "requested_tol", "effective_tol", "requested_inner_tol", "effective_inner_tol",
  "requested_precision", "effective_precision", "termination_reason"
)
cat(
  "missing_provenance_fields:",
  paste(setdiff(provenance_fields, names(tight_mixed)), collapse = ","),
  "\n"
)

# At 32x32 the documented public `precision = "double"` request is routed to
# the same mixed core as `precision = "mixed"` for scaling PGD.
large <- make_problem(32L, 32L)
large_double <- solve_problem(large, "double", 1e-9, 1e-9, max_iter = 100L)
large_mixed <- solve_problem(large, "mixed", 1e-9, 1e-9, max_iter = 100L)
cat("large_double_mixed_plan_identical:", identical(large_double$plan, large_mixed$plan), "\n")
cat("large_double_mixed_iterations:", large_double$iterations, large_mixed$iterations, "\n")
cat("large_double_backend_field:", large_double$backend, "\n")

stopifnot(
  tight_mixed$iterations < tight_mixed$max_iter,
  tight_mixed$residual > 1e-9,
  tight_mixed$residual < 1e-6,
  identical(tight_mixed$status, "max_iter"),
  identical(tight_mixed$plan, floor_mixed$plan),
  identical(default_precision$plan, tight_mixed$plan),
  identical(strict_double$status, "converged"),
  length(intersect(provenance_fields, names(tight_mixed))) == 0L,
  identical(large_double$plan, large_mixed$plan),
  identical(large_double$iterations, large_mixed$iterations)
)

sessionInfo()
