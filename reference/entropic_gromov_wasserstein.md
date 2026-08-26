# Entropic Gromov-Wasserstein (square loss)

POT-compatible GW wrapper using the optimized FGW entropic core with a
zero feature cost and `alpha = 1`.

## Usage

``` r
entropic_gromov_wasserstein(
  C1,
  C2,
  p = NULL,
  q = NULL,
  loss_fun = "square_loss",
  epsilon = 0.1,
  symmetric = NULL,
  G0 = NULL,
  max_iter = 1000L,
  tol = 1e-09,
  solver = c("PGD", "PPA"),
  sinkhorn_max_iter = 1000L,
  sinkhorn_tol = 1e-09,
  sinkhorn_method = c("auto", "scaling", "log"),
  precision = c("mixed", "double", "strict_double"),
  check_every = 1L,
  structure_rank = 0L
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

- epsilon:

  Entropic regularization.

- symmetric:

  Assume `C1` and `C2` are symmetric if `TRUE`.

- G0:

  Optional initial coupling (`ns x nt`) warm start.

- max_iter:

  Max outer GW iterations.

- tol:

  Outer stopping tolerance on plan updates.

- solver:

  Either `"PGD"` or `"PPA"`.

- sinkhorn_max_iter:

  Max Sinkhorn iterations per outer step.

- sinkhorn_tol:

  Sinkhorn stopping tolerance.

- sinkhorn_method:

  Sinkhorn variant (`"scaling"` or `"log"`).

- precision:

  Numeric precision mode (`"mixed"`, `"double"`, or `"strict_double"`).

- check_every:

  Evaluate stopping criterion every `check_every` iterations.

- structure_rank:

  Optional low-rank rank for structure tensor product.

## Value

A list with `plan`, `gw_dist`, `iterations`, `error`.
