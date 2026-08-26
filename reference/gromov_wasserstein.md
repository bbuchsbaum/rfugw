# Unregularized Gromov-Wasserstein via Conditional Gradient (square loss)

POT-compatible unregularized GW wrapper using the FGW
conditional-gradient core with zero feature cost and `alpha = 1`. This
is not a global optimizer of the non-convex GW objective.

## Usage

``` r
gromov_wasserstein(
  C1,
  C2,
  p = NULL,
  q = NULL,
  loss_fun = "square_loss",
  symmetric = NULL,
  G0 = NULL,
  max_iter = 500L,
  tol_rel = 1e-09,
  tol_abs = 1e-09,
  lp_solver = c("cpp_transport", "lp_transport", "lp_matrix"),
  lp_max_iter = 20000L,
  lp_tol = 1e-12
)
```

## Arguments

- C1:

  Source structure cost matrix (`ns x ns`).

- C2:

  Target structure cost matrix (`nt x nt`).

- p:

  Source weights (default uniform).

- q:

  Target weights (default uniform).

- loss_fun:

  Inner loss, currently only `"square_loss"` is supported.

- symmetric:

  Assume `C1` and `C2` are symmetric if `TRUE`.

- G0:

  Optional initial coupling (`ns x nt`) warm start.

- max_iter:

  Max outer GW iterations.

- tol_rel:

  Relative stopping tolerance for exact CG.

- tol_abs:

  Absolute stopping tolerance for exact CG.

- lp_solver:

  LP backend (`"cpp_transport"`, `"lp_transport"`, `"lp_matrix"`).

- lp_max_iter:

  Maximum iterations for C++ transport simplex LP backend.

- lp_tol:

  Optimality tolerance for C++ transport simplex LP backend.

## Value

A list with `plan`, `gw_dist`, `iterations`, `error`, `loss_trace`.
