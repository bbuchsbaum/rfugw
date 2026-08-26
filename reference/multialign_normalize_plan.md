# Normalize a Coupling Matrix

Normalize a Coupling Matrix

## Usage

``` r
multialign_normalize_plan(
  P,
  mode = c("none", "row", "col", "sum"),
  eps = 1e-15
)
```

## Arguments

- P:

  Coupling/transport matrix.

- mode:

  Normalization mode.

- eps:

  Stabilizer for division.

## Value

Normalized coupling matrix.
