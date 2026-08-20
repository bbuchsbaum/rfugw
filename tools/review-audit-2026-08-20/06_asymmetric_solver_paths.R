audit_lib <- Sys.getenv("RFUGW_AUDIT_LIB", unset = "")
if (nzchar(audit_lib)) {
  .libPaths(c(audit_lib, .libPaths()))
}

library(rfugw)

if (!requireNamespace("lpSolve", quietly = TRUE)) {
  stop("This probe requires lpSolve.")
}

gw_brute <- function(C1, C2, G) {
  out <- 0
  for (i in seq_len(nrow(C1))) {
    for (j in seq_len(nrow(C1))) {
      for (k in seq_len(nrow(C2))) {
        for (l in seq_len(nrow(C2))) {
          out <- out + (C1[i, j] - C2[k, l])^2 * G[i, k] * G[j, l]
        }
      }
    }
  }
  out
}

gw_grad_brute <- function(C1, C2, G) {
  grad <- matrix(0, nrow(G), ncol(G))
  for (a in seq_len(nrow(G))) {
    for (b in seq_len(ncol(G))) {
      for (j in seq_len(nrow(G))) {
        for (l in seq_len(ncol(G))) {
          grad[a, b] <- grad[a, b] +
            (C1[a, j] - C2[b, l])^2 * G[j, l]
        }
      }
      for (i in seq_len(nrow(G))) {
        for (k in seq_len(ncol(G))) {
          grad[a, b] <- grad[a, b] +
            (C1[i, a] - C2[k, b])^2 * G[i, k]
        }
      }
    }
  }
  grad
}

transport_lp <- function(cost, p, q) {
  fit <- lpSolve::lp.transport(
    cost.mat = cost,
    direction = "min",
    row.signs = rep("=", length(p)),
    row.rhs = p,
    col.signs = rep("=", length(q)),
    col.rhs = q,
    integers = integer(0)
  )
  if (fit$status != 0L) stop("lpSolve status ", fit$status)
  fit$solution
}

line_minimize <- function(objective, G0, Gc) {
  delta <- Gc - G0
  f0 <- objective(G0)
  fh <- objective(G0 + 0.5 * delta)
  f1 <- objective(Gc)
  a <- 2 * (f1 - f0) - 4 * (fh - f0)
  b <- (f1 - f0) - a
  candidates <- c(0, 1)
  if (is.finite(a) && a > 0) {
    stationary <- -b / (2 * a)
    if (stationary > 0 && stationary < 1) {
      candidates <- c(candidates, stationary)
    }
  }
  values <- vapply(candidates, function(step) objective(G0 + step * delta), numeric(1))
  step <- candidates[[which.min(values)]]
  list(plan = G0 + step * delta, step = step, objective = min(values))
}

set.seed(91)
ns <- 4L
nt <- 3L
C1 <- matrix(runif(ns * ns, 0.05, 1.2), ns, ns); diag(C1) <- 0
C2 <- matrix(runif(nt * nt, 0.05, 1.3), nt, nt); diag(C2) <- 0
M <- matrix(runif(ns * nt), ns, nt)
p <- c(0.1, 0.2, 0.3, 0.4)
q <- c(0.2, 0.3, 0.5)
alpha <- 0.6
epsilon <- 0.2
G0 <- p %o% q

fgw_objective <- function(G) {
  (1 - alpha) * sum(M * G) + alpha * gw_brute(C1, C2, G)
}
correct_linearization <- (1 - alpha) * M + alpha * gw_grad_brute(C1, C2, G0)

# One asymmetric entropic FGW step must equal Sinkhorn applied to the exact
# gradient linearization.
entropic_native <- fgw_entropic(
  M, C1, C2,
  p = p,
  q = q,
  alpha = alpha,
  epsilon = epsilon,
  max_iter = 1L,
  tol = 1e-15,
  sinkhorn_max_iter = 5000L,
  sinkhorn_tol = 1e-12,
  init_plan = G0,
  precision = "double",
  sinkhorn_method = "scaling",
  symmetric = FALSE
)
entropic_expected <- ot_sinkhorn(
  correct_linearization,
  p,
  q,
  epsilon = epsilon,
  method = "scaling",
  max_iter = 5000L,
  tol = 1e-12
)

# One asymmetric exact FGW step is checked against an independent LP direction
# and an objective-only quadratic line search.
exact_direction <- transport_lp(correct_linearization, p, q)
exact_expected <- line_minimize(fgw_objective, G0, exact_direction)
exact_native <- rfugw:::cpp_fgw_exact_cg_square(
  M, C1, C2, p, q, alpha, FALSE,
  1L, 0, 0, 20000L, 1e-12, G0
)

# The semirelaxed feasible direction puts each row mass on a minimum-gradient
# column. Its native one-step result is checked against the same brute objective.
q0 <- rep(1 / nt, nt)
G0_sr <- p %o% q0
sr_grad <- (1 - alpha) * M + alpha * gw_grad_brute(C1, C2, G0_sr)
sr_direction <- matrix(0, ns, nt)
for (i in seq_len(ns)) {
  js <- which(sr_grad[i, ] == min(sr_grad[i, ]))
  sr_direction[i, js] <- p[[i]] / length(js)
}
sr_expected <- line_minimize(fgw_objective, G0_sr, sr_direction)
sr_native <- rfugw:::cpp_semirelaxed_fgw_exact_square(
  M, C1, C2, p, alpha, FALSE, G0_sr, 1L, 0, 0
)

# Repeat the exact one-step construction for partial FGW. The native path uses
# the shared defective gradient helper; the reference direction uses the brute
# derivative and an independent extended transport LP.
m <- 0.7
G0_partial <- (p %o% q) * m
partial_grad <- (1 - alpha) * M + alpha * gw_grad_brute(C1, C2, G0_partial)
dummy_a <- sum(q) - m
dummy_b <- sum(p) - m
partial_cost_ext <- matrix(0, ns + 1L, nt + 1L)
partial_cost_ext[seq_len(ns), seq_len(nt)] <- partial_grad
penalty <- max(partial_grad) + max(abs(partial_grad)) + 1
partial_cost_ext[ns + 1L, nt + 1L] <- penalty
partial_ext <- transport_lp(partial_cost_ext, c(p, dummy_a), c(q, dummy_b))
partial_direction <- partial_ext[seq_len(ns), seq_len(nt), drop = FALSE]
partial_expected <- line_minimize(fgw_objective, G0_partial, partial_direction)
partial_native <- rfugw:::cpp_partial_fgw_exact_square(
  M, C1, C2, p, q, m, alpha, FALSE, G0_partial,
  1L, 0, 1L, 20000L, 1e-12
)

metrics <- c(
  entropic_one_step_plan_error = max(abs(entropic_native$plan - entropic_expected$plan)),
  exact_one_step_plan_error = max(abs(exact_native$plan - exact_expected$plan)),
  exact_objective_error = abs(exact_native$fgw_dist - exact_expected$objective),
  semirelaxed_one_step_plan_error = max(abs(sr_native$plan - sr_expected$plan)),
  semirelaxed_objective_error = abs(sr_native$srfgw_dist - sr_expected$objective),
  partial_one_step_plan_error = max(abs(partial_native$plan - partial_expected$plan)),
  partial_native_minus_reference_objective = partial_native$objective - partial_expected$objective
)
print(metrics, digits = 16)

stopifnot(
  metrics[["entropic_one_step_plan_error"]] < 1e-9,
  metrics[["exact_one_step_plan_error"]] < 1e-9,
  metrics[["exact_objective_error"]] < 1e-9,
  metrics[["semirelaxed_one_step_plan_error"]] < 1e-9,
  metrics[["semirelaxed_objective_error"]] < 1e-9,
  metrics[["partial_one_step_plan_error"]] > 1e-4,
  metrics[["partial_native_minus_reference_objective"]] > 1e-8
)

sessionInfo()
