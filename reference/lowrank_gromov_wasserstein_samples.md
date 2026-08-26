# Deprecated pseudo-low-rank GW compatibility wrapper

`lowrank_gromov_wasserstein_samples()` was misnamed: it always
materializes and solves a dense GW problem before applying SVD. Use
[`dense_gromov_wasserstein_plan_svd()`](https://bbuchsbaum.github.io/rfugw/reference/dense_gromov_wasserstein_plan_svd.md)
for the same explicit operation. POT-shaped factorized-cost, Dykstra,
and initialization parameters were never implemented and now error when
changed from their legacy defaults.

## Usage

``` r
lowrank_gromov_wasserstein_samples(
  X_s,
  X_t,
  a = NULL,
  b = NULL,
  reg = 0,
  rank = NULL,
  alpha = 1e-10,
  gamma_init = "rescale",
  rescale_cost = TRUE,
  cost_factorized_Xs = NULL,
  cost_factorized_Xt = NULL,
  stopThr = 1e-04,
  numItermax = 1000L,
  stopThr_dykstra = 0.001,
  numItermax_dykstra = 10000L,
  seed_init = 49,
  warn = TRUE,
  warn_dykstra = FALSE,
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

- alpha:

  Deprecated unsupported POT-compatibility parameter.

- gamma_init:

  Deprecated unsupported POT-compatibility parameter.

- rescale_cost:

  Whether to normalize structure costs to `[0, 1]`.

- cost_factorized_Xs:

  Deprecated unsupported factorized source cost.

- cost_factorized_Xt:

  Deprecated unsupported factorized target cost.

- stopThr:

  Outer GW stopping tolerance.

- numItermax:

  Outer GW max iterations.

- stopThr_dykstra:

  Deprecated unsupported Dykstra tolerance.

- numItermax_dykstra:

  Deprecated unsupported Dykstra iteration budget.

- seed_init:

  Deprecated unsupported initialization seed.

- warn:

  Deprecated unsupported warning control.

- warn_dykstra:

  Deprecated unsupported Dykstra warning control.

- log:

  If `TRUE`, return diagnostics.

## Value

The result of
[`dense_gromov_wasserstein_plan_svd()`](https://bbuchsbaum.github.io/rfugw/reference/dense_gromov_wasserstein_plan_svd.md).
