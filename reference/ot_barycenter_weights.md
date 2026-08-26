# Fixed-support Wasserstein barycenter weights

Optimizes one probability vector on a user-supplied common target
support. Every `costs[[s]]` has source measure `measures[[s]]` on its
rows and the common barycenter support on its columns. This is a
linear-Wasserstein estimand, not a GW/FGW support-learning update.

## Usage

``` r
ot_barycenter_weights(
  costs,
  measures,
  coefficients = NULL,
  mode = c("exact", "regularized"),
  support_cost = NULL,
  epsilon = 0.05,
  init_weights = NULL,
  init_state = NULL,
  min_weight = 1e-12,
  max_iter = NULL,
  tol = NULL,
  step_size = 1,
  backtracking = 30L,
  armijo = 1e-04,
  sinkhorn_method = c("auto", "scaling", "log"),
  sinkhorn_max_iter = 5000L,
  sinkhorn_tol = 1e-09
)
```

## Arguments

- costs:

  Nonempty list of finite nonnegative source-by-common-support cost
  matrices. All matrices must have the same column count.

- measures:

  List of source probability weights, one per cost matrix. Finite
  nonnegative vectors with positive mass are normalized separately.

- coefficients:

  Optional nonnegative barycenter coefficients, normalized to sum one.
  Zero coefficients are allowed.

- mode:

  `"exact"` or `"regularized"`.

- support_cost:

  In regularized mode, the finite nonnegative symmetric cost matrix
  within the common support. It supplies the debiasing self term.

- epsilon:

  Positive regularization in cost units for regularized mode.

- init_weights:

  Optional deterministic starting weights. Zeros are accepted and
  projected onto the declared numerical interior.

- init_state:

  Optional certified `warm_state` or prior converged regularized
  barycenter result from the identical problem.

- min_weight:

  Numerical lower bound for regularized weights. Exact mode has no such
  bound and can return exact zeros.

- max_iter:

  Maximum joint-simplex iterations in regularized mode and component
  transport-simplex iterations in exact mode. `NULL` selects 500 and
  20,000, respectively.

- tol:

  KKT and simplex tolerance. `NULL` selects `1e-6` in regularized mode
  and `1e-10` in exact mode.

- step_size:

  Initial projected-gradient step size.

- backtracking:

  Maximum Armijo reductions per outer iteration.

- armijo:

  Armijo sufficient-decrease coefficient in `(0, 1)`.

- sinkhorn_method:

  `"auto"`, `"scaling"`, or `"log"`.

- sinkhorn_max_iter:

  Maximum iterations for each component Sinkhorn solve.

- sinkhorn_tol:

  Marginal tolerance for each component Sinkhorn solve.

## Value

An `rfugw_barycenter_result` with `weights`, component results,
objective convention/components, simplex and KKT certificates, complete
solve counts, warm state, status, and provenance.

## Details

Exact mode solves one joint linear program and then independently
certifies every fixed-weight component with
[`ot_emd()`](https://bbuchsbaum.github.io/rfugw/reference/ot_emd.md).
Regularized mode minimizes a semi-debiased Sinkhorn objective: the
weighted cross product-reference-KL objectives minus half the
common-support self objective. Source-self terms are constant in the
barycenter weights and intentionally omitted. A projected-gradient KKT
map, safeguarded Armijo descent, and every component certificate are
required for convergence.

## Examples

``` r
support <- matrix(1:3, ncol = 1)
C <- outer(1:3, 1:3, function(x, y) (x - y)^2)
C <- C / max(C)
out <- ot_barycenter_weights(
  list(C, C),
  list(c(0.2, 0.5, 0.3), c(0.4, 0.2, 0.4)),
  mode = "regularized",
  support_cost = C,
  epsilon = 0.1
)
out$weights
#> [1] 0.3018473 0.3480833 0.3500694
```
