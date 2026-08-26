# Entropic Semirelaxed Fused Gromov-Wasserstein (square loss)

Mirror-descent solver in KL geometry with fixed source marginal and free
target marginal.

## Usage

``` r
entropic_semirelaxed_fused_gromov_wasserstein(
  M,
  C1,
  C2,
  p = NULL,
  loss_fun = "square_loss",
  symmetric = NULL,
  epsilon = 0.1,
  alpha = 0.5,
  G0 = NULL,
  max_iter = 10000L,
  tol = 1e-09,
  check_every = 10L,
  verbose = FALSE,
  precision = c("mixed", "double"),
  backend = c("cpp", "r")
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

- epsilon:

  Entropic regularization.

- alpha:

  Feature/structure trade-off in `[0, 1]`.

- G0:

  Optional initial semirelaxed plan (`ns x nt`) with row marginals `p`.

- max_iter:

  Max mirror-descent iterations.

- tol:

  Stopping tolerance on Frobenius norm of plan updates.

- check_every:

  Evaluate stopping criterion every `check_every` iterations.

- verbose:

  If `TRUE`, print iterative error diagnostics.

- precision:

  Numeric precision mode: `"mixed"` or `"double"`.

- backend:

  `"cpp"` or `"r"`.

## Value

A list with `plan`, `q`, `srfgw_dist`, `lin_loss`, `quad_loss`,
`iterations`, and `error`.
