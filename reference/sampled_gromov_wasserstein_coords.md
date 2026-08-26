# Sampled Gromov-Wasserstein Directly from Coordinates

POT-style sampled GW solver that operates on point coordinates and
computes required structure-distance rows on the fly. This avoids
storing dense `n x n` distance matrices and is suitable after
sparse-graph diffusion embeddings.

## Usage

``` r
sampled_gromov_wasserstein_coords(
  X1,
  X2,
  p = NULL,
  q = NULL,
  metric = c("euclidean", "sqeuclidean"),
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
  lp_scale = 1e+06,
  use_cpp = TRUE,
  objective_samples = 2000L,
  sampling = c("stochastic", "deterministic")
)
```

## Arguments

- X1:

  Source coordinates (`ns x d`).

- X2:

  Target coordinates (`nt x d`).

- p:

  Source weights (default uniform).

- q:

  Target weights (default uniform).

- metric:

  Ground metric between coordinates (`"euclidean"` or `"sqeuclidean"`).

- loss_fun:

  Currently only `"square_loss"` is supported.

- nb_samples_grad:

  Number of sampled gradient points, or length-2 vector
  `(n_source_samples, n_target_samples)`. Values below 1 error. Source
  or target counts above `ns` / `nt` warn and clamp.

- epsilon:

  Entropic regularization for projection step.

- max_iter:

  Maximum stochastic iterations.

- log:

  If `TRUE`, return diagnostics.

- verbose:

  If `TRUE`, print iterative diagnostics.

- random_state:

  Optional seed.

- sinkhorn_max_iter:

  Sinkhorn iterations when `epsilon > 0`.

- sinkhorn_tol:

  Sinkhorn tolerance when `epsilon > 0`.

- lp_solver:

  LP backend when `epsilon <= 0`.

- lp_scale:

  Integer scaling for LP marginals.

- use_cpp:

  If `TRUE` and C++ kernel is available, use the optimized entropic
  coordinate-native kernel when `epsilon > 0`.

- objective_samples:

  Number of Monte Carlo samples for `gw_dist_estimated` when
  `log = TRUE`.

- sampling:

  Sampling policy for stochastic updates. `"stochastic"` uses weighted
  random sampling (POT-style). `"deterministic"` uses top-k probability
  selections for reproducible parity checks.

## Value

If `log = FALSE`, returns coupling matrix `T`. If `log = TRUE`, returns
a list with `plan`, `gw_dist_estimated`, and `iterations`.

## Experimental

Coordinate-native sampled GW is experimental. The certified 0.1 envelope
matches
[`sampled_gromov_wasserstein()`](https://bbuchsbaum.github.io/rfugw/reference/sampled_gromov_wasserstein.md):
a full budget `(ns, nt)` is closer to dense entropic GW than a tiny
budget. Intermediate budgets are not certified as monotone. Inputs scale
as `O(n d)` rather than `O(n^2)` structure costs. See
`inst/bench/sampled-budget-curves.md`.
