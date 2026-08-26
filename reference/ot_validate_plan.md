# Validate a transport plan

Checks shape, nonnegativity, total mass, and optional marginal
constraints. Accepts a raw matrix, sparse/edge/operator transport plan,
or an `rfugw_result`. Validation uses mass operations and does not
materialize a sparse or implicit plan.

## Usage

``` r
ot_validate_plan(
  plan,
  p = NULL,
  q = NULL,
  mass = NULL,
  marginals = c("balanced", "partial", "relaxed"),
  tol = 1e-08
)
```

## Arguments

- plan:

  Coupling matrix or result object.

- p:

  Optional source weights.

- q:

  Optional target weights.

- mass:

  Optional required total mass.

- marginals:

  `"balanced"` (row/col match `p`/`q`), `"partial"` (row/col do not
  exceed `p`/`q` and total mass matches), or `"relaxed"` (nonnegativity
  and finiteness only).

- tol:

  Residual tolerance.

## Value

The plan, invisibly, after validation.

## Examples

``` r
M <- matrix(c(0, 1, 1, 0), 2, 2)
out <- ot_emd(M)
ot_validate_plan(out, c(0.5, 0.5), c(0.5, 0.5))
ot_linear_cost(M, out)
#> [1] 0
ot_entropy(out)
#> [1] -0.6931472
```
