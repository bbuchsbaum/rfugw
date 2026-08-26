# Debiased Sinkhorn divergence with certified component solves

Computes the debiased quantity
`OTe(source, target) - 0.5 * OTe(source, source) - 0.5 * OTe(target, target)`,
where every `OTe` uses the same convention
`<cost, plan> + epsilon * KL(plan || source_weights %o% target_weights)`.
All three component solves must earn fresh feasibility, objective, and
regularized primal-dual certificates.

## Usage

``` r
ot_sinkhorn_divergence(
  source = NULL,
  target = NULL,
  p = 2,
  source_weights = NULL,
  target_weights = NULL,
  cost = NULL,
  source_cost = NULL,
  target_cost = NULL,
  cost_power = NULL,
  metric = c("euclidean", "manhattan"),
  epsilon = 0.05,
  method = c("auto", "log", "scaling"),
  max_iter = 1000L,
  tol = 1e-09,
  init_state = NULL,
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

- source_cost:

  Optional source self-cost matrix. Required with `cost`.

- target_cost:

  Optional target self-cost matrix. Required with `cost`.

- cost_power:

  Positive power already represented by `cost`. The effective solver
  cost is `cost^(p / cost_power)`.

- metric:

  Raw-support metric, `"euclidean"` or `"manhattan"`.

- epsilon:

  Entropic regularization for the Sinkhorn solver only.

- method:

  Balanced Sinkhorn backend: `"auto"`, `"log"`, or `"scaling"`.

- max_iter:

  Optional solver iteration limit; backend defaults apply.

- tol:

  Optional solver tolerance; backend defaults apply.

- init_state:

  Optional list with reusable `cross`, `source_self`, and/or
  `target_self` dual states or prior results.

- allow_uncertified:

  If `FALSE` (default), fail instead of returning a value from a
  nonconverged solve. If `TRUE`, the result remains explicitly labeled
  uncertified.

## Value

An `rfugw_result` whose primary value is `sinkhorn_divergence`.
`component_solves`, `component_status`, `component_values`,
`component_residuals`, and `component_runtime_seconds` retain the
complete three-solve provenance.

## Examples

``` r
x <- c(0, 1, 3)
y <- c(0.5, 2, 4)
out <- ot_sinkhorn_divergence(x, y, p = 2, epsilon = 1)
rfugw_value(out)
#> [1] 0.7535445
```
