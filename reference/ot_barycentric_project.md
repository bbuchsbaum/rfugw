# Barycentric projection of a coupling

Maps source points through `plan` onto the target support, or the
reverse. Rows (or columns) with zero transported mass return `NaN` by
default.

## Usage

``` r
ot_barycentric_project(
  plan,
  points,
  orientation = c("source_to_target", "target_to_source"),
  zero_mass = c("nan", "zero")
)
```

## Arguments

- plan:

  Coupling or `rfugw_result`.

- points:

  Point matrix on the destination support (`nt x d` for
  `source_to_target`, `ns x d` for `target_to_source`).

- orientation:

  `"source_to_target"` or `"target_to_source"`.

- zero_mass:

  `"nan"` or `"zero"` for empty-mass rows.

## Value

Projected point matrix.

## Examples

``` r
plan <- diag(c(1, 0))
points <- matrix(c(0, 1, 2, 3), 2, 2, byrow = TRUE)
ot_barycentric_project(plan, points)
#>      [,1] [,2]
#> [1,]    0    1
#> [2,]  NaN  NaN
```
