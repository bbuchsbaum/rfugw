# Dense GW Plan Followed by Truncated SVD (Experimental)

Computes a dense entropic GW plan from sample distances, then compresses
that completed plan with truncated SVD. The name is deliberately
explicit: this is not an end-to-end low-rank or low-memory GW solver.

## Usage

``` r
dense_gromov_wasserstein_plan_svd(
  X_s,
  X_t,
  a = NULL,
  b = NULL,
  reg = 0,
  rank = NULL,
  rescale_cost = TRUE,
  stopThr = 1e-04,
  numItermax = 1000L,
  log = FALSE
)
```

## Arguments

- X_s:

  Source samples (`ns x d`).

- X_t:

  Target samples (`nt x d`).

- a:

  Optional source weights.

- b:

  Optional target weights.

- reg:

  Entropic regularization used for the inner GW solve.

- rank:

  Target low-rank factorization rank. Values below 1 error. Values above
  `min(ns, nt)` warn and clamp.

- rescale_cost:

  Whether to normalize structure costs to `[0, 1]`.

- stopThr:

  Outer GW stopping tolerance.

- numItermax:

  Outer GW max iterations.

- log:

  If `TRUE`, return diagnostics.

## Value

A list with `Q`, `R`, `g`, `representation`, `dense_plan_materialized`,
and byte counts for the dense structures and plan. If `log = TRUE`, it
also includes `plan`, `value_quad`, and `value`.

## Experimental

The solve materializes two dense square structure costs and a dense
rectangular plan before factorization. Peak memory is therefore not
bounded by the returned rank. Reconstruction error
`||T - T_r||_F / ||T||_F` decreases as rank increases up to
`min(ns, nt)`. See `inst/bench/sampled-budget-curves.md`.
