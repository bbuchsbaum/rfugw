# KL-unbalanced entropic optimal transport

KL-unbalanced entropic optimal transport

## Usage

``` r
ot_sinkhorn_unbalanced(
  M,
  p = NULL,
  q = NULL,
  epsilon = 0.05,
  rho = 10,
  max_iter = 500L,
  tol = 1e-07,
  init_plan = NULL,
  method = c("scaling", "log", "auto"),
  normalization = c("separate", "none", "joint")
)
```

## Arguments

- M:

  Cost matrix (`ns x nt`).

- p:

  Source finite nonnegative measure (default uniform).

- q:

  Target finite nonnegative measure (default uniform).

- epsilon:

  Positive entropic regularization.

- rho:

  Marginal generalized-KL penalties, length 1 or 2.

- max_iter:

  Maximum Sinkhorn iterations.

- tol:

  Marginal residual tolerance.

- init_plan:

  Optional warm start.

- method:

  `"scaling"`, genuine `"log"`, or `"auto"`. Auto uses the documented
  dynamic-range criterion from
  [`ot_sinkhorn()`](https://bbuchsbaum.github.io/rfugw/reference/ot_sinkhorn.md).

- normalization:

  Measure normalization policy. `"none"` preserves the supplied finite
  measures. `"joint"` applies one common factor so the mean
  source/target mass is one, preserving their mass ratio. `"separate"`
  normalizes each measure to mass one and is the backward-compatible
  default during the 0.1 transition.

## Value

An `rfugw_result`. `ot_dist` remains the unregularized transport term
for compatibility. `regularized_objective`, all three generalized-KL
terms, original/effective measures and masses, transported mass,
normalization policy, and a fixed-point certificate are also returned.

## Examples

``` r
M <- matrix(c(0, 1, 1, 0), 2, 2)
out <- ot_sinkhorn_unbalanced(M, epsilon = 0.1, rho = 2)
out$mass
#> [1] 0.9832372
```
