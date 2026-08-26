# Entropic Fused Gromov-Wasserstein Barycenters

POT-compatible alias for fixed-support entropic FGW barycenters.

## Usage

``` r
entropic_fused_gromov_barycenters(
  N,
  Ys,
  Cs,
  ps = NULL,
  p = NULL,
  lambdas = NULL,
  loss_fun = "square_loss",
  epsilon = 0.1,
  symmetric = TRUE,
  alpha = 0.5,
  max_iter = 1000L,
  tol = 1e-09,
  stop_criterion = c("barycenter", "loss"),
  warmstartT = FALSE,
  verbose = FALSE,
  log = FALSE,
  init_C = NULL,
  init_Y = NULL,
  fixed_structure = FALSE,
  fixed_features = FALSE,
  random_state = NULL,
  ...
)
```

## Arguments

- N:

  Number of barycenter nodes.

- Ys:

  List of feature matrices, one per sample (`ns x d`).

- Cs:

  List of structure matrices, one per sample (`ns x ns`).

- ps:

  Optional list of sample weights for each sample.

- p:

  Optional barycenter weights (`N`) (default uniform).

- lambdas:

  Optional sample weights across sets (default uniform).

- loss_fun:

  Inner loss, currently only `"square_loss"` is supported.

- epsilon:

  Entropic regularization used for inner FGW solves.

- symmetric:

  Assume symmetric structure matrices in inner FGW solves.

- alpha:

  FGW feature/structure trade-off.

- max_iter:

  Maximum outer barycenter iterations.

- tol:

  Stopping tolerance on barycenter updates.

- stop_criterion:

  Currently only `"barycenter"` is supported.

- warmstartT:

  If `TRUE`, warm-start inner FGW with previous couplings.

- verbose:

  If `TRUE`, print outer-iteration diagnostics.

- log:

  If `TRUE`, return full diagnostics; otherwise `Y` and `C`.

- init_C:

  Optional initial barycenter structure (`N x N`).

- init_Y:

  Optional initial barycenter features.

- fixed_structure:

  If `TRUE`, keep barycenter structure fixed at `init_C`.

- fixed_features:

  If `TRUE`, keep barycenter features fixed at `init_X`.

- random_state:

  Optional seed used only when `init_C` is `NULL`.

- ...:

  Additional arguments. Unused extras are rejected when the solver uses
  `.reject_unused_dots()`; otherwise they are forwarded to the primary
  solver.

## Value

If `log = FALSE`, returns a list with `Y` and `C`. If `log = TRUE`,
returns full diagnostics.
