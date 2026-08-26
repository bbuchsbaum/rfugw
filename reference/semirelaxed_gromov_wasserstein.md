# Semi-Relaxed Gromov-Wasserstein (non-entropic, square loss)

POT-compatible non-entropic semirelaxed GW solver via conditional
gradient.

## Usage

``` r
semirelaxed_gromov_wasserstein(
  C1,
  C2,
  p = NULL,
  loss_fun = "square_loss",
  symmetric = NULL,
  log = FALSE,
  G0 = NULL,
  max_iter = 10000L,
  tol_rel = 1e-09,
  tol_abs = 1e-09,
  random_state = 0,
  verbose = FALSE,
  ...
)
```

## Arguments

- C1:

  Source structure cost matrix (`ns x ns`).

- C2:

  Target structure cost matrix (`nt x nt`).

- p:

  Source weights (default uniform).

- loss_fun:

  Inner loss, currently only `"square_loss"` is supported.

- symmetric:

  If `NULL`, inferred from `C1` and `C2`; otherwise logical.

- log:

  If `TRUE`, include `loss_trace`.

- G0:

  Optional initial semirelaxed plan (`ns x nt`) satisfying
  `rowSums(G0) == p`.

- max_iter:

  Max mirror-descent iterations.

- tol_rel:

  Relative stopping tolerance.

- tol_abs:

  Absolute stopping tolerance.

- random_state:

  Ignored; retained for POT signature compatibility.

- verbose:

  If `TRUE`, print iterative error diagnostics.

- ...:

  Additional arguments. Unused extras are rejected when the solver uses
  `.reject_unused_dots()`; otherwise they are forwarded to the primary
  solver.

## Value

A list with `plan`, `q`, `srgw_dist`, `iterations`, `error`, and
`abs_error`.
