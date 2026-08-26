# Entropic term of a plan

Computes `sum(plan * log(plan))` with the convention `0 log 0 = 0`.

## Usage

``` r
ot_entropy(plan)
```

## Arguments

- plan:

  Coupling or `rfugw_result`.

## Value

Numeric scalar.
