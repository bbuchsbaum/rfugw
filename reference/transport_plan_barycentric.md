# Barycentric Projection Through a Transport Plan

Normalizes the forward or adjoint action by transported row/column mass.
Empty conditionals abstain as `NaN`, zero, or an error; they are never
replaced by a uniform match. Reverse barycentric projection is not an
inverse coupling.

## Usage

``` r
transport_plan_barycentric(
  plan,
  points,
  orientation = c("source_to_target", "target_to_source"),
  zero_mass = c("nan", "zero", "error")
)
```

## Arguments

- plan:

  Transport plan, operator, result, dense matrix, or sparse matrix.

- points:

  Destination point matrix.

- orientation:

  `"source_to_target"` or `"target_to_source"`.

- zero_mass:

  Empty-conditional abstention policy.

## Value

Projected point matrix.
