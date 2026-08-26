# Sampled GW from Sparse or Dense Similarity Graphs

Convenience wrapper that converts similarity graphs into diffusion
coordinates and solves sampled GW in coordinate space.

## Usage

``` r
sampled_gw_from_graphs(
  W1,
  W2,
  n_components = 16L,
  diffusion_time = 1,
  symmetrize = TRUE,
  self_loop = 1e-06,
  tol = 1e-12,
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
  return_embeddings = FALSE,
  sampling = c("stochastic", "deterministic")
)
```

## Arguments

- W1:

  Source similarity graph (`n1 x n1`) as dense matrix or sparse
  [`Matrix::sparseMatrix`](https://rdrr.io/pkg/Matrix/man/sparseMatrix.html).

- W2:

  Target similarity graph (`n2 x n2`) as dense matrix or sparse
  [`Matrix::sparseMatrix`](https://rdrr.io/pkg/Matrix/man/sparseMatrix.html).

- n_components:

  Number of diffusion coordinates.

- diffusion_time:

  Diffusion-time exponent.

- symmetrize:

  If `TRUE`, symmetrize each graph.

- self_loop:

  Optional diagonal regularization added before diffusion map.

- tol:

  Degree floor for graph validity checks.

- p:

  Source weights (default uniform).

- q:

  Target weights (default uniform).

- metric:

  Ground metric for coordinate GW.

- loss_fun:

  Currently only `"square_loss"` is supported.

- nb_samples_grad:

  Number of sampled gradient points, or length-2 vector
  `(n_source_samples, n_target_samples)`. Values below 1 error. Source
  or target counts above `ns` / `nt` warn and clamp.

- epsilon:

  Entropic regularization.

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

  Use C++ fast path for entropic coordinate solver when available.

- objective_samples:

  Monte Carlo samples for objective estimation.

- return_embeddings:

  If `TRUE`, include diffusion embeddings in output.

- sampling:

  Sampling policy passed to
  [`sampled_gromov_wasserstein_coords()`](https://bbuchsbaum.github.io/rfugw/reference/sampled_gromov_wasserstein_coords.md).

## Value

Coupling matrix or diagnostics list, mirroring
`sampled_gromov_wasserstein_coords`.

## Experimental

Graph-then-sampled GW is experimental. The certified 0.1 envelope is the
same quality-versus-budget claim as
[`sampled_gromov_wasserstein_coords()`](https://bbuchsbaum.github.io/rfugw/reference/sampled_gromov_wasserstein_coords.md),
plus that a sparse similarity graph plus `k` diffusion coordinates
stores less than two dense `n x n` structure costs. See
`inst/bench/sampled-budget-curves.md`.
