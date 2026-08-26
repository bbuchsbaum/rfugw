# Semi-Relaxed Fused Gromov-Wasserstein (non-entropic, square loss)

POT-compatible non-entropic semirelaxed FGW solver via conditional
gradient.

## Usage

``` r
semirelaxed_fused_gromov_wasserstein(
  M,
  C1,
  C2,
  p = NULL,
  loss_fun = "square_loss",
  symmetric = NULL,
  alpha = 0.5,
  G0 = NULL,
  log = FALSE,
  max_iter = 10000L,
  tol_rel = 1e-09,
  tol_abs = 1e-09,
  random_state = 0,
  verbose = FALSE,
  ...
)
```

## Arguments

- M:

  Cross-domain feature cost matrix (`ns x nt`).

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

- alpha:

  Feature/structure trade-off in `[0, 1]`.

- G0:

  Optional initial semirelaxed plan (`ns x nt`) with row marginals `p`.

- log:

  If `TRUE`, include `loss_trace`.

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

A list with `plan`, `q`, `srfgw_dist`, `lin_loss`, `quad_loss`,
`iterations`, `error`, and `abs_error`.
