# Exact balanced linear transport

Network-simplex / assignment backend. Reported `ot_dist` is `<M, plan>`.
An optimal result is accepted only when marginal feasibility, nonbasic
reduced costs, and the primal-dual gap pass the reported certificate
tolerances. Dual potentials and all certificate components are returned.

## Usage

``` r
ot_emd(M, p = NULL, q = NULL, max_iter = 20000L, tol = 1e-12)
```

## Arguments

- M:

  Cost matrix (`ns x nt`).

- p:

  Source weights (default uniform). Renormalized to sum 1.

- q:

  Target weights (default uniform). Renormalized to sum 1.

- max_iter:

  Maximum simplex iterations.

- tol:

  Optimality tolerance.

## Value

An `rfugw_result` with `plan`, `ot_dist`, `status`, exact termination
reason, source and target dual potentials, and optimality-certificate
diagnostics.

## Examples

``` r
M <- matrix(c(0, 2, 2, 0), 2, 2)
out <- ot_emd(M)
out$converged
#> [1] TRUE
ot_validate_plan(out, c(0.5, 0.5), c(0.5, 0.5))
```
