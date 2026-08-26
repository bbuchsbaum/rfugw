# Fused Unbalanced Gromov-Wasserstein (KL divergence, Sinkhorn inner solver)

POT-style BCD on two couplings with KL-unbalanced OT subproblems.

## Usage

``` r
fugw_kl(
  Cx,
  Cy,
  wx = NULL,
  wy = NULL,
  reg_marginals = c(10, 10),
  epsilon = 0.01,
  alpha = 0.5,
  M = NULL,
  init_pi = NULL,
  max_iter = 100L,
  tol = 1e-07,
  max_iter_ot = 500L,
  tol_ot = 1e-07,
  rescale_plan = TRUE,
  check_every = 1L,
  precision = c("double", "mixed", "strict_double"),
  feature_weight = NULL,
  structure_weight = NULL
)
```

## Arguments

- Cx:

  Source structure matrix.

- Cy:

  Target structure matrix.

- wx:

  Source weights (default uniform). Renormalized to sum 1.

- wy:

  Target weights (default uniform). Renormalized to sum 1.

- reg_marginals:

  Length-1 or length-2 marginal relaxation parameters.

- epsilon:

  Joint KL regularization parameter (entropic term).

- alpha:

  Legacy feature-cost coefficient, in `[0, 1]`; the structure
  coefficient is one. This convention differs from
  [`fgw_entropic()`](https://bbuchsbaum.github.io/rfugw/reference/fgw_entropic.md).

- M:

  Linear sample cost matrix (default `NULL`).

- init_pi:

  Optional initial sample coupling (`nx x ny`).

- max_iter:

  Max BCD iterations.

- tol:

  BCD stopping tolerance (`L1` delta on sample coupling).

- max_iter_ot:

  Max unbalanced Sinkhorn iterations per inner solve.

- tol_ot:

  Inner Sinkhorn tolerance.

- rescale_plan:

  Whether to rescale successive plans to equal mass.

- check_every:

  Evaluate BCD stopping criterion every `check_every` iterations.

- precision:

  Numeric precision mode for the C++ solver (`"double"`, `"mixed"`, or
  `"strict_double"`). Tight mixed requests are promoted to reported
  strict double rather than silently floored.

- feature_weight, structure_weight:

  Explicit nonnegative feature and structure coefficients. Supply both
  instead of `alpha`. Unlike FGW share aliases, FUGW coefficients are
  not normalized because their scale is meaningful relative to the KL
  penalties.

## Value

A list with `pi_samp`, `pi_feat`, `fugw_cost`, `linear_cost`,
`iterations`, `error`, `status`, and `converged`.

## Examples

``` r
set.seed(1)
Cx <- as.matrix(dist(matrix(rnorm(8), 4, 2)))
Cy <- as.matrix(dist(matrix(rnorm(10), 5, 2)))
out <- fugw_kl(Cx, Cy, epsilon = 0.05, max_iter = 20L)
out$status
#> [1] "inner_failure"
```
