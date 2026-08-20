audit_lib <- Sys.getenv("RFUGW_AUDIT_LIB", unset = "")
if (nzchar(audit_lib)) {
  .libPaths(c(audit_lib, .libPaths()))
}

library(rfugw)

if (!requireNamespace("lpSolve", quietly = TRUE)) {
  stop("This probe requires lpSolve for an independent optimum.")
}

lp_reference <- function(M, p, q) {
  fit <- lpSolve::lp.transport(
    cost.mat = M,
    direction = "min",
    row.signs = rep("=", length(p)),
    row.rhs = p,
    col.signs = rep("=", length(q)),
    col.rhs = q,
    integers = integer(0)
  )
  if (fit$status != 0L) {
    stop("lpSolve failed with status ", fit$status)
  }
  list(plan = fit$solution, objective = sum(M * fit$solution))
}

set.seed(20260820)
false_certificates <- list()
max_objective_gap <- 0
n_cases <- 600L

for (case_id in seq_len(n_cases)) {
  nr <- sample(2:6, 1)
  nc <- sample(2:6, 1)
  total_units <- 20L
  p_raw <- as.vector(rmultinom(1, total_units, rep(1, nr)))
  q_raw <- as.vector(rmultinom(1, total_units, rep(1, nc)))
  p <- p_raw / total_units
  q <- q_raw / total_units
  M <- matrix(sample(0:9, nr * nc, replace = TRUE), nr, nc)

  out <- rfugw::ot_emd(M, p, q, max_iter = 20000L, tol = 1e-12)
  ref <- lp_reference(M, p, q)
  objective_gap <- abs(out$ot_dist - ref$objective)
  max_objective_gap <- max(max_objective_gap, objective_gap)
  row_residual <- max(abs(rowSums(out$plan) - p))
  col_residual <- max(abs(colSums(out$plan) - q))

  if (isTRUE(out$converged) &&
      (objective_gap > 1e-8 || row_residual > 1e-8 || col_residual > 1e-8)) {
    false_certificates[[length(false_certificates) + 1L]] <- list(
      case_id = case_id,
      M = M,
      p = p,
      q = q,
      out = out,
      reference = ref,
      objective_gap = objective_gap,
      row_residual = row_residual,
      col_residual = col_residual
    )
    break
  }
}

cat("random_or_degenerate_cases:", n_cases, "\n")
cat("false_certificates_found:", length(false_certificates), "\n")
cat("max_objective_gap:", format(max_objective_gap, digits = 16), "\n")
if (length(false_certificates)) {
  dput(false_certificates[[1]])
}

# A max-iteration exit is a useful control: this branch is currently reported
# as lp_failure, unlike the early internal-failure branches in the source.
M_control <- matrix(c(7, 1, 4, 3, 9, 0, 5, 2, 8, 6, 1, 4), 3, 4)
p_control <- c(0.2, 0.3, 0.5)
q_control <- c(0.1, 0.2, 0.3, 0.4)
limited <- ot_emd(M_control, p_control, q_control, max_iter = 1L)
cat("limited_status:", limited$status, "\n")
cat("limited_iterations:", limited$iterations, "\n")

# Partial exact FGW extracts the inner simplex plan without returning its
# convergence flag. With alpha=0 the outer problem exposes the LP direction
# directly, making a one-pivot and a fully solved direction comparable.
set.seed(174)
ns <- 4L
nt <- 5L
M <- matrix(runif(ns * nt), ns, nt)
C1 <- matrix(0, ns, ns)
C2 <- matrix(0, nt, nt)
p <- c(0.1, 0.2, 0.3, 0.4)
q <- c(0.08, 0.12, 0.2, 0.25, 0.35)
m <- 0.7
G0 <- (p %o% q) * m

partial_limited <- rfugw:::cpp_partial_fgw_exact_square(
  M, C1, C2, p, q, m, 0, FALSE, G0,
  1L, 0, 1L, 1L, 1e-12
)
partial_full <- rfugw:::cpp_partial_fgw_exact_square(
  M, C1, C2, p, q, m, 0, FALSE, G0,
  1L, 0, 1L, 20000L, 1e-12
)

cat("partial_return_fields:", paste(names(partial_limited), collapse = ","), "\n")
cat(
  "partial_limited_vs_full_objective_gap:",
  format(abs(partial_limited$objective - partial_full$objective), digits = 16),
  "\n"
)
cat(
  "partial_limited_vs_full_plan_max_error:",
  format(max(abs(partial_limited$plan - partial_full$plan)), digits = 16),
  "\n"
)

stopifnot(
  identical(limited$status, "lp_failure"),
  limited$iterations == 1L,
  !"lp_ok" %in% names(partial_limited),
  !"inner_iterations" %in% names(partial_limited),
  abs(partial_limited$objective - partial_full$objective) > 1e-8
)

sessionInfo()
