# Entropic Partial Fused Gromov-Wasserstein (square loss)

POT-compatible entropic partial FGW solver.

## Usage

``` r
entropic_partial_fused_gromov_wasserstein(
  M,
  C1,
  C2,
  p = NULL,
  q = NULL,
  reg = 1,
  m = NULL,
  loss_fun = "square_loss",
  alpha = 0.5,
  G0 = NULL,
  numItermax = 1000L,
  tol = 1e-07,
  symmetric = NULL,
  log = FALSE,
  verbose = FALSE,
  check_every = 2L,
  inner_max_iter = 300L,
  inner_tol = 1e-12,
  method = c("auto", "scaling", "log")
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

- reg:

  Entropic regularization parameter (\>0).

- m:

  Amount of mass to transport. Defaults to `min(sum(p), sum(q))`.

- loss_fun:

  Inner loss, currently only `"square_loss"` is supported.

- alpha:

  FGW tradeoff in `[0, 1]`.

- G0:

  Optional initial coupling (`ns x nt`) warm start.

- numItermax:

  Maximum CG iterations.

- tol:

  Relative stopping tolerance on the partial objective.

- symmetric:

  Assume `C1` and `C2` are symmetric if `TRUE`.

- log:

  If `TRUE`, return a list with diagnostics; otherwise the plan.

- verbose:

  If `TRUE`, print CG diagnostics.

- check_every:

  Outer stopping check interval.

- inner_max_iter:

  Maximum iterations for inner entropic partial OT solve.

- inner_tol:

  Inner stopping tolerance for partial OT solve.

- method:

  Certified public partial-Sinkhorn backend: `"auto"`, bounded
  `"scaling"`, or genuine `"log"` Dykstra.

## Value

If `log = FALSE`, returns a coupling matrix. If `log = TRUE`, returns a
list with `plan`, `partial_fgw_dist`, `lin_loss`, `quad_loss`,
`iterations`, `error`, and `err_trace`.
