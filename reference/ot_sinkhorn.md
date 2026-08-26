# Balanced entropic optimal transport

Scaling or log-domain Sinkhorn for a linear cost. Reported `ot_dist` is
the unregularized `<M, plan>` cost. `regularized_objective` uses the
explicit convention `<M, plan> + epsilon * KL(plan || p %o% q)`; its
independently derived dual and gap are returned without changing the
legacy field.

## Usage

``` r
ot_sinkhorn(
  M,
  p = NULL,
  q = NULL,
  epsilon = 0.05,
  method = c("scaling", "log", "auto"),
  max_iter = 1000L,
  tol = 1e-09,
  init_plan = NULL,
  init_duals = NULL
)
```

## Arguments

- M:

  Cost matrix (`ns x nt`).

- p:

  Source weights (default uniform). Renormalized to sum 1.

- q:

  Target weights (default uniform). Renormalized to sum 1.

- epsilon:

  Positive entropic regularization.

- method:

  `"scaling"`, `"log"`, or `"auto"`. Auto uses scaling only when the
  maximum scaled exponent magnitude and span are at most 500; otherwise
  it selects the genuine log-domain implementation.

- max_iter:

  Maximum Sinkhorn iterations.

- tol:

  Marginal residual tolerance.

- init_plan:

  Optional strictly positive entropic plan warm start. It must have the
  requested shape and be zero outside any zero-weight support.

- init_duals:

  Optional reusable dual state, either a prior `rfugw_result` or a list
  with numeric `source` and `target` potentials. When both
  initialization forms are supplied, `init_duals` takes precedence after
  both inputs are validated.

## Value

An `rfugw_result` with `plan`, `ot_dist`, `status`, residuals, canonical
`source_potential` / `target_potential`, and reusable `dual_state`. The
gauge is fixed by a zero source-weighted mean. The regularized primal,
dual, product-reference KL, constant offset, and gap certificate are
also returned.

## Examples

``` r
M <- matrix(c(0, 1, 1, 0), 2, 2)
out <- ot_sinkhorn(M, epsilon = 0.1)
out$status
#> [1] "converged"
rfugw_value(out)
#> [1] 4.539787e-05
```
