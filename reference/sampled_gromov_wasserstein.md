# Sampled Gromov-Wasserstein (square loss)

POT-style stochastic GW estimator using sampled gradients and balanced
OT projection steps.

## Usage

``` r
sampled_gromov_wasserstein(
  C1,
  C2,
  p = NULL,
  q = NULL,
  loss_fun = "square_loss",
  nb_samples_grad = 100L,
  epsilon = 1,
  max_iter = 500L,
  log = FALSE,
  verbose = FALSE,
  random_state = NULL,
  sinkhorn_max_iter = 200L,
  sinkhorn_tol = 1e-09,
  lp_solver = c("lp_matrix", "lp_transport"),
  lp_scale = 1e+06
)
```

## Arguments

- C1:

  Source structure matrix.

- C2:

  Target structure matrix.

- p:

  Source weights (default uniform).

- q:

  Target weights (default uniform).

- loss_fun:

  Currently only `"square_loss"` is supported.

- nb_samples_grad:

  Number of sampled gradient points, or length-2 vector
  `(n_source_samples, n_target_samples)`. Values below 1 error. Source
  or target counts above `ns` / `nt` warn and clamp.

- epsilon:

  Entropic regularization for the OT projection step. If `<= 0`, exact
  LP projection is used.

- max_iter:

  Maximum stochastic iterations.

- log:

  If `TRUE`, return objective estimate and diagnostics.

- verbose:

  If `TRUE`, print iterative diagnostics.

- random_state:

  Optional seed.

- sinkhorn_max_iter:

  Sinkhorn iterations when `epsilon > 0`.

- sinkhorn_tol:

  Sinkhorn tolerance when `epsilon > 0`.

- lp_solver:

  LP backend used when `epsilon <= 0`.

- lp_scale:

  Integer scaling for LP marginals.

## Value

If `log = FALSE`, returns coupling matrix `T`. If `log = TRUE`, returns
a list with `plan`, `gw_dist_estimated`, and `iterations`.

## Experimental

Sampled GW is experimental. The certified 0.1 envelope is that a full
budget `(ns, nt)` is closer to dense
[`entropic_gromov_wasserstein()`](https://bbuchsbaum.github.io/rfugw/reference/entropic_gromov_wasserstein.md)
than a tiny budget such as `(2, 1)`, in square-loss GW and plan
Frobenius distance. Intermediate budgets are not certified as monotone.
See `inst/bench/sampled-budget-curves.md`.
