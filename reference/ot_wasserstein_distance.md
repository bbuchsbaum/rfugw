# Wasserstein p-distance with explicit root semantics

Calls the same certified transport path as
[`ot_wasserstein_cost()`](https://bbuchsbaum.github.io/rfugw/reference/ot_wasserstein_cost.md)
but makes the returned primary value `(minimum d^p cost)^(1/p)`. For
`p = 1` this is exactly the transport cost; no square root is applied.

## Usage

``` r
ot_wasserstein_distance(
  source = NULL,
  target = NULL,
  p = 2,
  source_weights = NULL,
  target_weights = NULL,
  cost = NULL,
  cost_power = NULL,
  metric = c("euclidean", "manhattan"),
  solver = c("exact", "sinkhorn"),
  epsilon = NULL,
  sinkhorn_method = NULL,
  max_iter = NULL,
  tol = NULL,
  init_plan = NULL,
  init_duals = NULL,
  allow_uncertified = FALSE
)
```

## Arguments

- source:

  Source support as a numeric vector or rows-by-features matrix.

- target:

  Target support with the same feature dimension.

- p:

  Positive Wasserstein ground-metric power.

- source_weights:

  Nonnegative source weights; normalized to probability.

- target_weights:

  Nonnegative target weights; normalized to probability.

- cost:

  Optional nonnegative source-by-target cost matrix. When supplied,
  `source` and `target` must be `NULL` and `cost_power` is required.

- cost_power:

  Positive power already represented by `cost`. The effective solver
  cost is `cost^(p / cost_power)`.

- metric:

  Raw-support metric, `"euclidean"` or `"manhattan"`.

- solver:

  `"exact"` for certified EMD or `"sinkhorn"` for the distinctly labeled
  cost of a certified entropic plan.

- epsilon:

  Entropic regularization for the Sinkhorn solver only.

- sinkhorn_method:

  `"auto"`, `"log"`, or `"scaling"` for Sinkhorn only.

- max_iter:

  Optional solver iteration limit; backend defaults apply.

- tol:

  Optional solver tolerance; backend defaults apply.

- init_plan:

  Optional Sinkhorn plan warm start.

- init_duals:

  Optional Sinkhorn dual warm state.

- allow_uncertified:

  If `FALSE` (default), fail instead of returning a value from a
  nonconverged solve. If `TRUE`, the result remains explicitly labeled
  uncertified.

## Value

An `rfugw_result` whose
[`rfugw_value()`](https://bbuchsbaum.github.io/rfugw/reference/rfugw_value.md)
is `wasserstein_p_distance`. The unrooted value remains available as
`wasserstein_p_cost`.

## Examples

``` r
out <- ot_wasserstein_distance(c(0, 2), c(1, 3), p = 1)
rfugw_value(out)
#> [1] 1
```
