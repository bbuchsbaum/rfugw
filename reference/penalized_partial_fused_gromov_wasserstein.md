# Penalized Variable-Mass Partial Fused Gromov-Wasserstein

Optimizes the square-loss FGW objective over nonnegative subcouplings
while choosing the transported mass. Unmatched source and target mass
each incur `discard_penalty` per unit. The outer problem is nonconvex: a
converged result certifies subcoupling feasibility, every exact linear
direction solve, objective and line-search consistency, and a
Frank-Wolfe stationarity gap, but not global optimality.

## Usage

``` r
penalized_partial_fused_gromov_wasserstein(
  M,
  C1,
  C2,
  discard_penalty,
  p = NULL,
  q = NULL,
  alpha = 0.5,
  G0 = NULL,
  numItermax = 500L,
  tol = 1e-08,
  inner_max_iter = 20000L,
  inner_tol = 1e-12,
  symmetric = NULL,
  trace = FALSE,
  verbose = FALSE
)
```

## Arguments

- M:

  Finite nonnegative source-by-target feature cost.

- C1, C2:

  Finite square within-domain structure matrices.

- discard_penalty:

  Finite nonnegative cost per unmatched unit on each marginal.
  Increasing it weakly favors more transported mass in each exact
  linearized subproblem.

- p, q:

  Finite nonnegative source and target measures. Defaults are uniform
  probability measures.

- alpha:

  Feature/structure tradeoff in `[0, 1]`.

- G0:

  Optional feasible nonnegative subcoupling. Its total mass is free; row
  and column sums must not exceed `p` and `q`.

- numItermax:

  Maximum Frank-Wolfe updates.

- tol:

  Relative Frank-Wolfe-gap tolerance.

- inner_max_iter:

  Maximum iterations for each certified penalized linear OT direction
  solve.

- inner_tol:

  Tolerance for each direction solve.

- symmetric:

  Whether to use the symmetric square-loss formula. By default it is
  inferred from `C1` and `C2`.

- trace:

  If `TRUE`, retain every accepted feasible iterate and its diagnostics
  in `debug_trace`.

- verbose:

  If `TRUE`, print iteration diagnostics.

## Value

An `rfugw_result` containing the plan, unrooted objective terms,
transported/discarded masses, feasibility and objective certificates,
Frank-Wolfe gap, line-search evidence, and nested exact-solver status.

## Details

The unrooted objective is

This differs from
[`partial_fused_gromov_wasserstein()`](https://bbuchsbaum.github.io/rfugw/reference/partial_fused_gromov_wasserstein.md),
which fixes `sum(G)`, and from
[`fugw_kl()`](https://bbuchsbaum.github.io/rfugw/reference/fugw_kl.md),
which uses generalized-KL marginal relaxation.

## Examples

``` r
C1 <- matrix(c(0, 1, 1, 0), 2, 2)
C2 <- matrix(c(0, 2, 2, 0), 2, 2)
M <- matrix(c(0.1, 2, 2, 0.2), 2, 2)
out <- penalized_partial_fused_gromov_wasserstein(
  M, C1, C2, discard_penalty = 1
)
out$transported_mass
#> [1] 1
out$frank_wolfe_gap
#> [1] 0
```
