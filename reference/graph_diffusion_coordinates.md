# Diffusion Coordinates from a Similarity Graph

Build diffusion-map style coordinates from a dense or sparse nonnegative
similarity graph. This provides a memory-aware route from sparse kNN
graphs to geometry vectors that can be aligned with sampled GW without
materializing dense `n x n` structure matrices.

## Usage

``` r
graph_diffusion_coordinates(
  W,
  n_components = 16L,
  diffusion_time = 1,
  symmetrize = TRUE,
  self_loop = 0,
  tol = 1e-12,
  log = FALSE
)
```

## Arguments

- W:

  Similarity graph (`n x n`) as a dense numeric matrix or a sparse
  [`Matrix::sparseMatrix`](https://rdrr.io/pkg/Matrix/man/sparseMatrix.html).

- n_components:

  Number of diffusion coordinates to return.

- diffusion_time:

  Diffusion-time exponent applied to eigenvalues.

- symmetrize:

  If `TRUE`, symmetrize graph as `(W + t(W))/2`.

- self_loop:

  Optional diagonal value added to all nodes.

- tol:

  Degree floor used for validity checks.

- log:

  If `TRUE`, return diagnostics (`eigenvalues`, `degree`, `solver`).

## Value

If `log = FALSE`, returns a dense matrix of coordinates (`n x k`). If
`log = TRUE`, returns a list with `coords`, `eigenvalues`, `degree`, and
`solver`.
