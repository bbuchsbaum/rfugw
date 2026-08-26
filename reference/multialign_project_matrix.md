# Project a Subject Matrix into Template Space

Computes `t(Pn) %*% R %*% Pn` where `Pn` is the normalized coupling.

## Usage

``` r
multialign_project_matrix(
  P,
  R,
  normalize = c("none", "row", "col", "sum"),
  eps = 1e-15
)
```

## Arguments

- P:

  Subject-to-template coupling (`n_subject x n_template`).

- R:

  Subject square matrix (`n_subject x n_subject`).

- normalize:

  Plan normalization mode passed to
  [`multialign_normalize_plan()`](https://bbuchsbaum.github.io/rfugw/reference/multialign_normalize_plan.md).

- eps:

  Stabilizer for normalization.

## Value

Projected template-space matrix.
