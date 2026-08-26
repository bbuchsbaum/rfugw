# Linear optimal transport in rfugw

Linear optimal transport matches two weighted collections when the cost
of every source-to-target match is known. You supply a cost matrix and,
optionally, weights for the two collections. `rfugw` returns a coupling:
a nonnegative matrix whose entries say how much mass moves between each
pair.

This vignette builds one coupling, checks that it is trustworthy, and
uses it to map target coordinates back to the source. It then shows when
to choose an exact, entropic, or unbalanced solver.

``` r

library(rfugw)
```

## What are the inputs?

We will align two small point clouds in the same two-dimensional space.
Rows are points; `M[i, j]` is the Euclidean cost of matching source
point `i` to target point `j`. The weights `p` and `q` each sum to one.

``` r

set.seed(16)
X <- matrix(rnorm(16), ncol = 2)
Y <- matrix(rnorm(20), ncol = 2) + 0.4
M <- as.matrix(dist(rbind(X, Y)))[1:8, 9:18]
M <- M / max(M)
p <- rep(1 / nrow(X), nrow(X))
q <- rep(1 / nrow(Y), nrow(Y))
```

Linear OT requires a meaningful cross-domain cost. If your two supports
live in different metric spaces and such a cost is unavailable, use GW
or FGW instead; those methods compare within-domain structure.

## How do you compute a balanced plan?

[`ot_sinkhorn()`](https://bbuchsbaum.github.io/rfugw/reference/ot_sinkhorn.md)
is the usual starting point for a fast, dense coupling. Its `epsilon`
argument controls entropic regularization: larger values spread mass
more diffusely. The reported value is the unregularized transport cost
`sum(M * plan)`; it does not include the entropy penalty.

``` r

sink <- ot_sinkhorn(M, p, q, epsilon = 0.08, method = "auto")
data.frame(
  status = sink$status,
  cost = rfugw_value(sink),
  marginal_residual = sink$residual
)
#>      status      cost marginal_residual
#> 1 converged 0.2983975       2.41073e-10
```

The marginal residual measures how closely the row and column sums
reproduce `p` and `q`.
[`ot_validate_plan()`](https://bbuchsbaum.github.io/rfugw/reference/ot_validate_plan.md)
provides an explicit feasibility check; `print(sink)` and
`summary(sink)` offer compact diagnostics without printing the full
coupling.

When solving a sequence of nearby regularization values, reuse the
returned dual state instead of restarting from scratch. The state has a
stable gauge (the source-weighted potential mean is zero), works across
scaling and log backends, and contains no native pointer.

``` r

sink_finer <- ot_sinkhorn(
  M, p, q,
  epsilon = 0.05,
  method = "auto",
  init_duals = sink$dual_state
)
data.frame(
  status = sink_finer$status,
  initialization = sink_finer$initialization,
  iterations = sink_finer$iterations,
  residual = sink_finer$residual
)
#>      status initialization iterations     residual
#> 1 converged          duals         85 5.289409e-10
```

You may instead pass a previous entropic plan through `init_plan`. If
both are provided, both are validated and `init_duals` takes precedence.
A warm start never carries over convergence: the new solve must satisfy
its requested marginal certificate. Because `dual_state` is an ordinary
named list, it can be saved with
[`saveRDS()`](https://rdrr.io/r/base/readRDS.html) and reused in an
installed package from a fresh session.

## How do you use the coupling?

A barycentric projection replaces each source point by the weighted
average of the target points to which it sends mass.

``` r

projected <- ot_barycentric_project(sink, Y)
head(round(projected, 3), 4)
#>        [,1]   [,2]
#> [1,]  1.162  0.988
#> [2,]  0.371  0.550
#> [3,]  1.185  1.204
#> [4,] -0.911 -0.209
```

``` r

mean_nearest_distance <- function(A, B) {
  D <- as.matrix(dist(rbind(A, B)))[
    seq_len(nrow(A)), nrow(A) + seq_len(nrow(B)), drop = FALSE
  ]
  mean(apply(D, 1, min))
}
projection_distance <- c(
  before = mean_nearest_distance(X, Y),
  after = mean_nearest_distance(projected, Y)
)
```

``` r

data.frame(
  coordinates = names(projection_distance),
  mean_distance_to_nearest_target = unname(projection_distance)
)
#>   coordinates mean_distance_to_nearest_target
#> 1      before                       0.6350109
#> 2       after                       0.4004487
```

The result has one row per source point and the same coordinate columns
as `Y`. A source row with no transported mass cannot have a weighted
average and is returned as `NaN`; use `zero_mass = "zero"` if zeros are
preferable. In this example, the projected coordinates are closer to the
target support on average. That is an interpreted example result, not a
general guarantee for every cost or geometry.

## Which solver should you choose?

| Task | Function | Main trade-off |
|----|----|----|
| Balanced, exact linear OT | [`ot_emd()`](https://bbuchsbaum.github.io/rfugw/reference/ot_emd.md) | Sparse exact plan, best suited to small or moderate problems |
| Balanced, entropic linear OT | [`ot_sinkhorn()`](https://bbuchsbaum.github.io/rfugw/reference/ot_sinkhorn.md) | Faster dense plan; `epsilon` changes the solution |
| Exact transport of a specified fraction | [`ot_partial_emd()`](https://bbuchsbaum.github.io/rfugw/reference/ot_partial_emd.md) | Enforces a transported mass in `[0, 1]` |
| Entropic transport of a specified fraction | [`ot_partial_sinkhorn()`](https://bbuchsbaum.github.io/rfugw/reference/ot_partial_sinkhorn.md) | Dense fixed-mass subcoupling with certified scaling/log Dykstra |
| Exact transport with optimized mass | [`ot_partial_penalized()`](https://bbuchsbaum.github.io/rfugw/reference/ot_partial_penalized.md) | A linear discard penalty selects how much mass to transport |
| Marginals may change softly | [`ot_sinkhorn_unbalanced()`](https://bbuchsbaum.github.io/rfugw/reference/ot_sinkhorn_unbalanced.md) | KL penalties control the cost of creating or removing mass |
| Marginals change on sparse support | [`ot_sinkhorn_unbalanced_ti()`](https://bbuchsbaum.github.io/rfugw/reference/ot_sinkhorn_unbalanced_ti.md) | Certified TI potentials and operator application avoid a dense plan |

For example, the exact solver satisfies the same marginals without
entropic regularization.

``` r

exact <- ot_emd(M, p, q)
data.frame(
  converged = exact$converged,
  cost = rfugw_value(exact),
  nonzero_entries = sum(rfugw_plan(exact) > 1e-12)
)
#>   converged      cost nonzero_entries
#> 1      TRUE 0.2740385              16
```

The exact and entropic costs are useful summaries of their respective
plans, but comparing them alone is not a convergence test. Inspect the
solver status and validate the returned marginals.

## How do sparse plans remain operators?

Exact transport plans often have few nonzero edges. Converting those
edges to an `rfugw_transport_plan` preserves mass, linear application,
objectives, and barycentric projection without allocating a dense
source-by-target matrix.

``` r

edge_index <- which(exact$plan > 1e-12, arr.ind = TRUE)
exact_edges <- as_transport_plan(
  data.frame(
    source = edge_index[, 1],
    target = edge_index[, 2],
    weight = exact$plan[edge_index]
  ),
  n_source = nrow(exact$plan),
  n_target = ncol(exact$plan)
)
data.frame(
  representation = transport_plan_representation(exact_edges)$representation,
  stored_edges = nrow(exact_edges$data),
  transported_mass = transport_plan_mass(exact_edges),
  linear_cost = ot_linear_cost(M, exact_edges)
)
#>   representation stored_edges transported_mass linear_cost
#> 1      edge_list           16                1   0.2740385
```

[`transport_plan_apply()`](https://bbuchsbaum.github.io/rfugw/reference/transport_plan_apply.md)
computes the forward linear action;
[`transport_plan_adjoint()`](https://bbuchsbaum.github.io/rfugw/reference/transport_plan_adjoint.md)
computes the transpose action. The adjoint is not an inverse. Likewise,
target-to-source barycentric projection normalizes the adjoint by
transported column mass and should not be called an inverse map. Rows or
columns with no transported mass return `NaN` by default, an explicit
abstention rather than a fabricated uniform match. Materialization is
available only by calling
[`transport_plan_materialize()`](https://bbuchsbaum.github.io/rfugw/reference/transport_plan_materialize.md)
(or `rfugw_plan(x, materialize = TRUE)`) explicitly. Pruning reports
lost mass and removes convergence certificates from a result until the
changed plan is separately revalidated.

## What if only a fixed fraction should match?

When the overlap is known but not all mass should participate,
[`ot_partial_sinkhorn()`](https://bbuchsbaum.github.io/rfugw/reference/ot_partial_sinkhorn.md)
transports exactly the requested `mass` while treating the supplied
weights as row and column capacities. It uses the explicit regularized
objective

``` text
sum(M * plan) + epsilon * sum(plan * (log(plan) - 1)).
```

This counting-measure entropy convention differs from the
product-reference KL reported by balanced
[`ot_sinkhorn()`](https://bbuchsbaum.github.io/rfugw/reference/ot_sinkhorn.md),
so their regularized values should not be mixed.

``` r

partial_sink <- ot_partial_sinkhorn(
  M, p, q, mass = 0.7, epsilon = 0.08,
  method = "auto", tol = 1e-9
)
data.frame(
  status = partial_sink$status,
  backend = partial_sink$effective_sinkhorn_method,
  transported_mass = partial_sink$transported_mass,
  transport_term = partial_sink$transport_objective,
  entropy_term = partial_sink$weighted_entropy_minus_one,
  regularized_objective = rfugw_value(partial_sink),
  primal_dual_gap = partial_sink$duality_gap
)
#>      status backend transported_mass transport_term entropy_term
#> 1 converged scaling              0.7      0.1337304   -0.2468603
#>   regularized_objective primal_dual_gap
#> 1              -0.11313   -6.183588e-11
```

Auto dispatch uses bounded scaling when safe and genuine log Dykstra for
large cost-to-entropy range. A raw plan is not a valid warm start
because Dykstra also needs its correction state; pass a prior certified
result or its `warm_state` only for the identical problem.

## When is the result a Wasserstein distance?

A linear transport cost becomes a Wasserstein p-cost only after the
ground metric power is declared.
[`ot_wasserstein_cost()`](https://bbuchsbaum.github.io/rfugw/reference/ot_wasserstein_cost.md)
computes `sum(d^p * plan)`;
[`ot_wasserstein_distance()`](https://bbuchsbaum.github.io/rfugw/reference/ot_wasserstein_distance.md)
returns its `1 / p` root. This distinction prevents accidentally taking
a square root of W1 or treating a squared structural loss as a ground
distance.

For supports in a common coordinate system, supply the points directly.
Rows are observations and columns are features.

``` r

w2 <- ot_wasserstein_distance(
  source = X,
  target = Y,
  p = 2,
  source_weights = p,
  target_weights = q,
  metric = "euclidean",
  solver = "exact"
)
data.frame(
  p_cost = w2$wasserstein_p_cost,
  p_distance = rfugw_value(w2),
  power = w2$wasserstein_power,
  certification = w2$value_certification
)
#>     p_cost p_distance power             certification
#> 1 1.104775   1.051083     2 certified_exact_transport
```

If a cost matrix already exists, declare its current power. For example,
if `M2` stores squared Euclidean distances, use `cost_power = 2`. The
helper can then transform it consistently when the requested `p`
differs.

``` r

M_raw <- as.matrix(dist(rbind(X, Y)))[1:8, 9:18]
M2 <- M_raw^2
w2_from_cost <- ot_wasserstein_distance(
  cost = M2,
  cost_power = 2,
  p = 2,
  source_weights = p,
  target_weights = q
)
stopifnot(abs(rfugw_value(w2) - rfugw_value(w2_from_cost)) < 1e-10)
```

A user-supplied nonnegative matrix is accepted as a declared
dissimilarity, but `metric_certified` is `FALSE`: rfugw does not infer
symmetry or the triangle inequality from its values. The helpers
normalize each positive-mass weight vector to a probability measure and
retain original/effective mass provenance.

With `solver = "sinkhorn"`, the result is labeled
`estimate_kind = "entropic_plan"`. It is the unregularized p-cost of a
certified entropic plan, not the exact Wasserstein optimum and not the
entropy-regularized objective. Nonconverged solves fail by default.
`allow_uncertified = TRUE` exists only for explicit diagnostic
inspection and leaves `value_certified = FALSE`.

## What is the regularized objective and Sinkhorn divergence?

[`ot_sinkhorn()`](https://bbuchsbaum.github.io/rfugw/reference/ot_sinkhorn.md)
keeps `ot_dist = sum(M * plan)` for compatibility. It now also reports
the full regularized primal

``` text
sum(M * plan) + epsilon * KL(plan || p %o% q),
```

plus an independently evaluated dual, their gap, and the exact constant
offset from the alternative `sum(plan * (log(plan) - 1))` entropy
convention. This matters when comparing epsilon values: the linear term
alone is not the regularized optimum, and conventions that differ by
constants must not be mixed.

``` r

data.frame(
  linear_plan_cost = sink$ot_dist,
  regularized_primal = sink$regularized_objective,
  regularized_dual = sink$regularized_dual_objective,
  primal_dual_gap = sink$regularized_duality_gap
)
#>   linear_plan_cost regularized_primal regularized_dual primal_dual_gap
#> 1        0.2983975          0.3692867        0.3692867    -7.85309e-11
```

Entropic OT is biased: even a measure compared with itself generally has
a nonzero regularized value.
[`ot_sinkhorn_divergence()`](https://bbuchsbaum.github.io/rfugw/reference/ot_sinkhorn_divergence.md)
removes that self-bias by combining one cross solve and two self solves
under the same convention.

``` r

sinkhorn_div <- ot_sinkhorn_divergence(
  X, Y,
  p = 2,
  source_weights = p,
  target_weights = q,
  epsilon = 0.2,
  method = "auto",
  max_iter = 5000L,
  tol = 1e-7
)
data.frame(
  divergence = rfugw_value(sinkhorn_div),
  status = sinkhorn_div$status,
  cross = sinkhorn_div$component_values[["cross"]],
  source_self = sinkhorn_div$component_values[["source_self"]],
  target_self = sinkhorn_div$component_values[["target_self"]]
)
#>   divergence    status    cross source_self target_self
#> 1   1.005309 converged 1.406336   0.4087245   0.3933297
```

The divergence result retains the three complete component results,
statuses, residuals, iteration counts, and runtimes. It is certified
only when all three solves pass and the combined value is nonnegative
within the reported numerical tolerance. A nonconverged self solve is
never hidden behind a finite cross value.

## When should transported mass be optimized?

[`ot_partial_emd()`](https://bbuchsbaum.github.io/rfugw/reference/ot_partial_emd.md)
is appropriate when the transported mass is known in advance. When it is
part of the estimand,
[`ot_partial_penalized()`](https://bbuchsbaum.github.io/rfugw/reference/ot_partial_penalized.md)
instead balances transport cost against a linear penalty for mass left
unmatched on either side. The supplied weights remain finite measures;
they are not separately normalized.

``` r

penalized_low <- ot_partial_penalized(
  M, p, q, discard_penalty = 0.05
)
penalized_high <- ot_partial_penalized(
  M, p, q, discard_penalty = 0.5
)
data.frame(
  discard_penalty = c(0.05, 0.5),
  transported_mass = c(
    penalized_low$transported_mass,
    penalized_high$transported_mass
  ),
  transport_term = c(
    penalized_low$transport_term,
    penalized_high$transport_term
  ),
  discard_term = c(
    penalized_low$discard_penalty_term,
    penalized_high$discard_penalty_term
  )
)
#>   discard_penalty transported_mass transport_term discard_term
#> 1            0.05              0.1    0.006289643         0.09
#> 2            0.50              1.0    0.274038542         0.00
```

Each transported unit avoids one source and one target discard penalty,
so a larger penalty weakly favors more transport. The result reports
both discarded masses, the transport and penalty objective terms, and an
augmented exact-OT primal/dual certificate scaled back to the original
mass units. This is a linear variable-mass formulation, not the KL
marginal relaxation below.

## What if the marginals should change?

Contamination, missing support, or partial overlap can make fixed
marginals a poor model.
[`ot_sinkhorn_unbalanced()`](https://bbuchsbaum.github.io/rfugw/reference/ot_sinkhorn_unbalanced.md)
replaces equality constraints with KL penalties. Smaller `rho` permits
more marginal deviation; larger `rho` approaches the balanced problem.

``` r

unbalanced <- ot_sinkhorn_unbalanced(
  M, p, q, epsilon = 0.08, rho = 0.2, method = "auto",
  max_iter = 1000L, tol = 1e-7, normalization = "none"
)
data.frame(
  status = unbalanced$status,
  cost = rfugw_value(unbalanced),
  regularized_objective = unbalanced$regularized_objective,
  transported_mass = unbalanced$mass
)
#>      status      cost regularized_objective transported_mass
#> 1 converged 0.1116993             0.2387464        0.5026116
```

An unbalanced plan is not expected to reproduce `p` and `q` exactly. Its
total mass and relaxed marginal residuals are part of the result and
should be interpreted alongside the application-specific choice of
`rho`. Here `normalization = "none"` means the supplied finite measures,
including their absolute and unequal total masses, enter the
generalized-KL penalties unchanged. Use `"joint"` to apply one common
scale while preserving the mass ratio. The backward-compatible 0.1
default, `"separate"`, normalizes each side to one and should be
selected only when the inputs represent probability measures. The result
records original and effective masses and separates the transport term
(`ot_dist`) from the full `regularized_objective`.

When the allowed matches themselves are sparse, use the
translation-invariant path. Here each source has only two or three
admissible targets; the edge list can represent a genuine zero cost
without confusing it with a missing sparse- matrix entry.

``` r

ti_index <- do.call(rbind, lapply(seq_len(nrow(M)), function(i) {
  target <- unique(c(i, i + 1L, if (i == nrow(M)) ncol(M)))
  data.frame(source = i, target = target, cost = M[i, target])
}))
ti_uot <- ot_sinkhorn_unbalanced_ti(
  ti_index, p, q,
  epsilon = 0.08, rho = c(0.2, 0.3),
  n_source = nrow(M), n_target = ncol(M),
  plan = "operator", max_iter = 5000L, tol = 1e-8
)
data.frame(
  status = ti_uot$status,
  representation = ti_uot$plan_representation$representation,
  support_size = ti_uot$support_size,
  regularized_objective = rfugw_value(ti_uot),
  kkt_residual = ti_uot$kkt_residual,
  primal_dual_gap = ti_uot$primal_dual_gap
)
#>      status    representation support_size regularized_objective kkt_residual
#> 1 converged implicit_operator           17             0.3330139 7.674851e-08
#>   primal_dual_gap
#> 1    1.110223e-16
```

The operator applies the coupling or its adjoint from the stored cost
support and certified potentials. Request `plan = "sparse"` only when a
downstream consumer needs explicit weighted edges, or `plan = "dense"`
when the dense allocation is intentional. Unlike the older unbalanced
API’s legacy
[`rfugw_value()`](https://bbuchsbaum.github.io/rfugw/reference/rfugw_value.md)
convention, the TI result’s value is its full regularized objective. Its
support certificate reports both active-node coverage and a separate
exact max-flow answer to the hypothetical balanced-marginal problem;
unequal masses can make the latter false without making KL-UOT
infeasible.

## How does linear OT connect to FGW?

FGW combines a cross-domain feature cost `M` with structure matrices
`C1` and `C2`. `alpha = 0` uses only the linear feature term, while
`alpha = 1` uses only structure. Intermediate values blend the two.

``` r

C1 <- as.matrix(dist(X)); C1 <- C1 / max(C1)
C2 <- as.matrix(dist(Y)); C2 <- C2 / max(C2)
fgw <- fgw_entropic(
  M, C1, C2, alpha = 0.5, epsilon = 0.08, max_iter = 40L
)
data.frame(
  status = fgw$status,
  converged = fgw$converged,
  inner_status = fgw$inner_status,
  objective = fgw$fgw_dist
)
#>      status converged inner_status objective
#> 1 converged      TRUE    converged 0.1960542
```

[`ot_fgw_square()`](https://bbuchsbaum.github.io/rfugw/reference/ot_fgw_square.md)
recomputes the unregularized FGW objective from a plan. Its agreement
with the reported value checks the objective calculation; solver status
and marginal feasibility remain separate diagnostics.

## Where should you go next?

- [`vignette("solver-guide")`](https://bbuchsbaum.github.io/rfugw/articles/solver-guide.md)
  explains how to choose among GW, partial, semirelaxed, and FUGW
  formulations.
- [`vignette("barycenters")`](https://bbuchsbaum.github.io/rfugw/articles/barycenters.md)
  constructs a representative metric-measure space from several inputs.
- See
  [`?ot_sinkhorn`](https://bbuchsbaum.github.io/rfugw/reference/ot_sinkhorn.md),
  [`?ot_sinkhorn_divergence`](https://bbuchsbaum.github.io/rfugw/reference/ot_sinkhorn_divergence.md),
  [`?ot_emd`](https://bbuchsbaum.github.io/rfugw/reference/ot_emd.md),
  [`?ot_partial_sinkhorn`](https://bbuchsbaum.github.io/rfugw/reference/ot_partial_sinkhorn.md),
  [`?ot_partial_penalized`](https://bbuchsbaum.github.io/rfugw/reference/ot_partial_penalized.md),
  [`?ot_wasserstein_distance`](https://bbuchsbaum.github.io/rfugw/reference/ot_wasserstein_distance.md),
  and
  [`?ot_validate_plan`](https://bbuchsbaum.github.io/rfugw/reference/ot_validate_plan.md)
  for full argument and diagnostic contracts.
