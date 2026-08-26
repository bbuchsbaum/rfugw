# Explicitly Materialize a Transport Plan

This is the only generic operation that intentionally allocates a dense
source-by-target matrix. Implicit operators must provide a materializer.

## Usage

``` r
transport_plan_materialize(plan)
```

## Arguments

- plan:

  Transport plan, operator, result, dense matrix, or sparse matrix.

## Value

Dense numeric matrix.
