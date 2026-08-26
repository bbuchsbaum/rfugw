# Unregularized FGW via Conditional Gradient + LP Direction (square loss)

Conditional-gradient iterations with an exact linear-OT direction step.
The outer GW/FGW problem is non-convex; the returned plan is a
stationary point of this procedure, not a certified global minimizer.
Reported `fgw_dist` is the unregularized FGW objective.

## Usage

``` r
fgw_exact_cg(
  M,
  C1,
  C2,
  p = NULL,
  q = NULL,
  alpha = 0.5,
  symmetric = NULL,
  G0 = NULL,
  max_iter = 500L,
  tol_rel = 1e-09,
  tol_abs = 1e-09,
  lp_scale = 1e+06,
  lp_solver = c("cpp_transport", "lp_transport", "lp_matrix"),
  lp_max_iter = 20000L,
  lp_tol = 1e-12,
  feature_weight = NULL,
  structure_weight = NULL
)
```

## Arguments

- M:

  Cross-domain feature cost matrix (`ns x nt`).

- C1:

  Source structure cost matrix (`ns x ns`).

- C2:

  Target structure cost matrix (`nt x nt`).

- p:

  Source weights (default uniform). Renormalized to sum 1.

- q:

  Target weights (default uniform). Renormalized to sum 1.

- alpha:

  Trade-off between feature and structure terms, in `[0, 1]`. It is the
  structure share, so the feature share is `1 - alpha`.

- symmetric:

  `NULL` auto-detects symmetry of `C1` and `C2`. `TRUE` requires both
  costs to be symmetric within `1e-10`. `FALSE` uses the two-sided
  tensor.

- G0:

  Optional feasible initial coupling (`ns x nt`). Honored by every
  `lp_solver` backend.

- max_iter:

  Max outer FGW iterations.

- tol_rel:

  Relative stopping tolerance on objective change.

- tol_abs:

  Absolute stopping tolerance on objective change.

- lp_scale:

  Integer mass scaling used for LP marginals (lpSolve backends).

- lp_solver:

  LP backend for direction step: `"cpp_transport"` (default, C++
  transport-simplex backend), `"lp_transport"` (lpSolve), or
  `"lp_matrix"` (cached dense lpSolve constraints).

- lp_max_iter:

  Maximum iterations for the C++ transport-simplex LP backend.

- lp_tol:

  Optimality tolerance for the C++ transport-simplex LP backend.

- feature_weight, structure_weight:

  Explicit nonnegative aliases for the feature and structure shares.
  Supply both instead of `alpha`; they are normalized to sum to one.

## Value

A list with `plan`, `fgw_dist`, `iterations`, `error`, `rel_error`,
`loss_trace`, `status`, and `converged`.

## Examples

``` r
set.seed(1)
C1 <- as.matrix(dist(matrix(rnorm(8), 4, 2)))
C2 <- as.matrix(dist(matrix(rnorm(8), 4, 2)))
M <- matrix(runif(16), 4, 4)
out <- fgw_exact_cg(M, C1, C2, max_iter = 20L)
out$status
#> [1] "converged"
```
