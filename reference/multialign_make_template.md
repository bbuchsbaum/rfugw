# Build a Fixed Template Support for Multiset Alignment

Builds a fixed template support from pooled subject node features using
k-means, then defines template geometry as pairwise distances between
template feature centroids.

## Usage

``` r
multialign_make_template(
  subjects,
  k = NULL,
  feature_normalization = c("none", "zscore"),
  feature_cost_metric = c("euclidean", "sqeuclidean"),
  seed = 1L
)
```

## Arguments

- subjects:

  A list of subject sets, each with at least `C` and optional `F`, `w`,
  `id`.

- k:

  Number of template nodes. Defaults to rounded median subject size.

- feature_normalization:

  Feature normalization mode.

- feature_cost_metric:

  Feature distance metric.

- seed:

  Optional random seed for k-means reproducibility.

## Value

A template list with `C`, `F`, `w`, `id`.
