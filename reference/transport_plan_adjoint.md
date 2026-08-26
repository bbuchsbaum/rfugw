# Apply the Transport Adjoint

Computes `t(plan) %*% x` without materialization. The adjoint reverses
the linear action; it is not an inverse, a reverse conditional
distribution, or a barycentric reverse map.

## Usage

``` r
transport_plan_adjoint(plan, x)
```

## Arguments

- plan:

  Transport plan, operator, result, dense matrix, or sparse matrix.

- x:

  Numeric vector or matrix with one row per source.

## Value

Numeric vector or matrix with one row per target.
