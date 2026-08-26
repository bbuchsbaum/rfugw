# Partial Fused Gromov-Wasserstein (square loss)

POT-compatible exact partial FGW solver using conditional gradient.

## Usage

``` r
partial_fused_gromov_wasserstein(
  M,
  C1,
  C2,
  p = NULL,
  q = NULL,
  m = NULL,
  loss_fun = "square_loss",
  alpha = 0.5,
  nb_dummies = 1L,
  G0 = NULL,
  thres = 1,
  numItermax = 10000L,
  tol = 1e-08,
  symmetric = NULL,
  warn = TRUE,
  log = FALSE,
  verbose = FALSE,
  lp_solver = c("cpp_transport", "lp_matrix", "lp_transport"),
  lp_scale = 1e+06,
  ...
)
```

## Arguments

- M:

  Cross-domain feature cost matrix.

- C1:

  Source structure cost matrix (`ns x ns`).

- C2:

  Target structure cost matrix (`nt x nt`).

- p:

  Source weights (default uniform).

- q:

  Target weights (default uniform).

- m:

  Amount of mass to transport. Defaults to `min(sum(p), sum(q))`.

- loss_fun:

  Inner loss, currently only `"square_loss"` is supported.

- alpha:

  FGW tradeoff in `[0, 1]`.

- nb_dummies:

  Number of dummy nodes used in the partial OT linearized step.

- G0:

  Optional initial coupling (`ns x nt`) warm start.

- thres:

  Unused POT compatibility parameter.

- numItermax:

  Maximum CG iterations.

- tol:

  Relative stopping tolerance on the partial objective.

- symmetric:

  Assume `C1` and `C2` are symmetric if `TRUE`.

- warn:

  Ignored; retained for POT signature compatibility.

- log:

  If `TRUE`, return a list with diagnostics; otherwise the plan.

- verbose:

  If `TRUE`, print CG diagnostics.

- lp_solver:

  LP backend for the linearized partial OT step (`"cpp_transport"`
  default, `"lp_transport"`, `"lp_matrix"`).

- lp_scale:

  Integer scaling factor for LP marginal discretization.

- ...:

  Additional arguments. Unused extras are rejected when the solver uses
  `.reject_unused_dots()`; otherwise they are forwarded to the primary
  solver.

## Value

If `log = FALSE`, returns a coupling matrix. If `log = TRUE`, returns a
list with `plan`, `partial_fgw_dist`, `lin_loss`, `quad_loss`,
`iterations`, `error`, and `loss_trace`.
