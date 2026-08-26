# Entropic Fused Gromov-Wasserstein (square loss)

POT-style projected gradient iterations with Sinkhorn projections.
Reported `fgw_dist` is the unregularized FGW objective at the returned
plan; it does not add the entropic term.

## Usage

``` r
fgw_entropic(
  M,
  C1,
  C2,
  p = NULL,
  q = NULL,
  alpha = 0.5,
  epsilon = 0.1,
  max_iter = 1000L,
  tol = 1e-09,
  sinkhorn_max_iter = 1000L,
  sinkhorn_tol = 1e-09,
  init_plan = NULL,
  structure_rank = 0L,
  sinkhorn_method = c("auto", "scaling", "log"),
  precision = c("mixed", "double", "strict_double"),
  symmetric = NULL,
  solver = c("PGD", "PPA"),
  check_every = 1L,
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

- epsilon:

  Entropic regularization (must be positive).

- max_iter:

  Max outer FGW iterations.

- tol:

  Outer stopping tolerance on Frobenius norm of plan updates.

- sinkhorn_max_iter:

  Max Sinkhorn iterations per outer step.

- sinkhorn_tol:

  Sinkhorn stopping tolerance on marginal residual.

- init_plan:

  Optional initial coupling (`ns x nt`) used as a warm start.

- structure_rank:

  Optional low-rank approximation rank for structure matrices in the FGW
  tensor product (`0` disables approximation).

- sinkhorn_method:

  Sinkhorn variant: `"scaling"`, `"log"`, or `"auto"`. Auto evaluates a
  conservative initial-gradient/bound proxy and selects scaling only
  inside the precision-specific certified regime.

- precision:

  Numeric precision mode: `"mixed"` (default), `"double"`, or
  `"strict_double"` (`"mixed"` uses float inner iterations with double
  final objective evaluation). A mixed request tighter than `1e-6` is
  promoted and reported as `"strict_double"`. For larger problems,
  `"double"` may use a reported mixed-precision accelerated path;
  `"strict_double"` never dispatches to float arithmetic.

- symmetric:

  `NULL` auto-detects symmetry of `C1` and `C2`. `TRUE` requires both
  costs to be symmetric within `1e-10`. `FALSE` uses the two-sided
  tensor.

- solver:

  Either `"PGD"` or `"PPA"`.

- check_every:

  Evaluate outer stopping error every `check_every` iterations.

- feature_weight, structure_weight:

  Explicit nonnegative aliases for the feature and structure shares.
  Supply both instead of `alpha`; they are normalized to sum to one and
  the interpreted weights are returned.

## Value

A list with `plan`, `fgw_dist`, `iterations`, `error`, `residual`,
`status`, `converged`, and marginal residuals. See
`inst/solver-contract.md`.

## Examples

``` r
set.seed(1)
C1 <- as.matrix(dist(matrix(rnorm(8), 4, 2)))
C2 <- as.matrix(dist(matrix(rnorm(10), 5, 2)))
M <- matrix(runif(20), 4, 5)
out <- fgw_entropic(M, C1, C2, epsilon = 0.1, max_iter = 40L)
out$status
#> [1] "converged"
out$residual
#> [1] 3.0505e-10
```
