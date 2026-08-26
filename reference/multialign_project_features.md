# Project Template Features to Subject Nodes

Computes `Pn %*% E_template`, typically with row-normalized `Pn`.

## Usage

``` r
multialign_project_features(
  P,
  E_template,
  normalize = c("row", "none", "col", "sum"),
  eps = 1e-15
)
```

## Arguments

- P:

  Subject-to-template coupling (`n_subject x n_template`).

- E_template:

  Template feature/embedding matrix (`n_template x d`).

- normalize:

  Plan normalization mode.

- eps:

  Stabilizer for normalization.

## Value

Subject-space projected feature matrix (`n_subject x d`).
