# Apply a Transport Plan as a Linear Operator

Computes `plan %*% x` without materializing edge-list, sparse, or
implicit plans. This is the forward action from target-indexed values to
source rows; it is not row-normalized.

## Usage

``` r
transport_plan_apply(plan, x)
```

## Arguments

- plan:

  Transport plan, operator, result, dense matrix, or sparse matrix.

- x:

  Numeric vector or matrix with one row per target.

## Value

Numeric vector or matrix with one row per source.
