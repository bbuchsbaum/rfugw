# Entropic Gromov-Wasserstein Barycenters

POT-compatible GW barycenter wrapper using the optimized FGW barycenter
core with zero features and `alpha = 1`.

## Usage

``` r
entropic_gromov_barycenters(
  N,
  Cs,
  ps = NULL,
  p = NULL,
  lambdas = NULL,
  loss_fun = "square_loss",
  epsilon = 0.1,
  symmetric = TRUE,
  max_iter = 1000L,
  tol = 1e-09,
  stop_criterion = c("barycenter", "loss"),
  warmstartT = FALSE,
  verbose = FALSE,
  log = FALSE,
  init_C = NULL,
  random_state = NULL,
  ...
)
```

## Arguments

- N:

  Number of barycenter nodes.

- Cs:

  List of structure matrices.

- ps:

  Optional list of source weights.

- p:

  Optional barycenter weights.

- lambdas:

  Optional barycenter sample weights.

- loss_fun:

  Currently only `"square_loss"`.

- epsilon:

  Entropic regularization.

- symmetric:

  Symmetry flag for inner GW solves.

- max_iter:

  Max outer barycenter iterations.

- tol:

  Outer stopping tolerance.

- stop_criterion:

  Currently only `"barycenter"` is supported.

- warmstartT:

  Warm-start inner GW solves.

- verbose:

  Print diagnostics.

- log:

  Return diagnostics/history.

- init_C:

  Optional initial barycenter structure.

- random_state:

  Optional seed for initialization.

- ...:

  Additional arguments. Unused extras are rejected when the solver uses
  `.reject_unused_dots()`; otherwise they are forwarded to the primary
  solver.

## Value

If `log = FALSE`, returns barycenter structure matrix `C`. If
`log = TRUE`, returns a list with `C`, `p`, `couplings`, `history`, and
`objective`.
