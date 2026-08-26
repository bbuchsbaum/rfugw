# Minimal downstream conformance fixture. This intentionally uses only exported
# rfugw functions and base R so a manifold-alignment client does not depend on
# native or package-private Sinkhorn state.
manifoldalign_ot_procrustes_step <- function(X, Y, epsilon, state = NULL) {
  stopifnot(is.matrix(X), is.matrix(Y), ncol(X) == ncol(Y))
  x2 <- rowSums(X^2)
  y2 <- rowSums(Y^2)
  cost <- outer(x2, y2, "+") - 2 * tcrossprod(X, Y)
  cost[cost < 0 & cost > -1e-12] <- 0
  fit <- rfugw::ot_sinkhorn(
    cost,
    epsilon = epsilon,
    method = "auto",
    max_iter = 10000L,
    tol = 1e-8,
    init_duals = state
  )
  coupling <- rfugw::rfugw_plan(fit)
  cross <- crossprod(X, coupling %*% Y)
  decomposition <- svd(cross)
  rotation <- decomposition$u %*% t(decomposition$v)
  list(
    rotation = rotation,
    result = fit,
    state = fit$dual_state
  )
}
