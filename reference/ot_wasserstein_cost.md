# Wasserstein p-cost with explicit cost-power semantics

Solves classical balanced transport with ground cost `d^p` and returns a
certificate-rich result. Supply either raw supports or a nonnegative
cost matrix whose existing power is declared explicitly. A Sinkhorn
result is labeled as the cost of its certified entropic plan, not as
exact Wasserstein cost and not as a regularized objective.

## Usage

``` r
ot_wasserstein_cost(
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

An `rfugw_result`. `wasserstein_p_cost` is the transport cost,
`wasserstein_p_distance` is its `1 / p` root, and
[`rfugw_value()`](https://bbuchsbaum.github.io/rfugw/reference/rfugw_value.md)
returns the p-cost for this helper.

## Examples

``` r
out <- ot_wasserstein_cost(c(0, 2), c(1, 3), p = 2)
rfugw_value(out)
#> [1] 1
out$wasserstein_p_distance
#> [1] 1
```
