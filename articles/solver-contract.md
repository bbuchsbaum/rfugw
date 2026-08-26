# A solver-client protocol for downstream libraries

The solver-client protocol gives an R package one stable route from an
explicit transport estimand to a certified coupling. The generic example
comes first here: DKGE and manifoldalign use the same boundary, but
neither package defines its semantics.

``` r

library(rfugw)
```

## What data and estimand are we declaring?

Suppose a client has seven source objects, five target objects, a
rectangular cross-cost matrix, and probability weights. Rows of the cost
are always source support; columns are target support. Our estimand is
balanced entropic linear transport, and the expected output is a 7-by-5
coupling with the supplied row and column marginals.

``` r

set.seed(31)
source_points <- matrix(rnorm(21), 7, 3)
target_points <- matrix(rnorm(15), 5, 3)
cost <- outer(rowSums(source_points^2), rowSums(target_points^2), "+") -
  2 * tcrossprod(source_points, target_points)
source_mass <- seq_len(7) / sum(seq_len(7))
target_mass <- rev(seq_len(5)) / sum(seq_len(5))

problem <- transport_problem_sinkhorn(
  cost,
  p = source_mass,
  q = target_mass,
  epsilon = 0.7,
  method = "auto",
  tol = 1e-9,
  mass_policy = "probability"
)
problem
#> <rfugw_transport_problem v1.0>
#>   estimand:   balanced_entropic_linear_ot
#>   solver:     ot_sinkhorn
#>   orientation:rows_source_columns_target
#>   shape:      7 x 5
#>   mass policy:probability
```

`mass_policy = "probability"` verifies rather than silently rescales.
Use `"normalize"` only when separate normalization is the intended
estimand; the problem then retains both original and effective masses.
KL-unbalanced problems instead require an explicit finite-measure,
joint, or separate policy.

## How does a client solve and inspect it?

``` r

fit <- transport_solve(problem)
plan <- rfugw_plan(fit)

data.frame(
  status = rfugw_status(fit),
  rows = nrow(plan),
  columns = ncol(plan),
  row_residual = rfugw_residuals(fit)$row_residual,
  column_residual = rfugw_residuals(fit)$col_residual
)
#>      status rows columns row_residual column_residual
#> 1 converged    7       5 9.880044e-10    5.551115e-17
```

The public accessors separate four concerns:

- [`rfugw_problem()`](https://bbuchsbaum.github.io/rfugw/reference/rfugw_problem.md)
  returns the declared estimand and requested controls;
- [`rfugw_masses()`](https://bbuchsbaum.github.io/rfugw/reference/rfugw_masses.md)
  separates input, effective, and transported mass;
- [`rfugw_provenance()`](https://bbuchsbaum.github.io/rfugw/reference/rfugw_provenance.md)
  records requested/effective backend controls, objective units,
  representation, status, and runtime;
- the established plan, value, status, and residual accessors remain the
  authoritative numerical boundary.

``` r

masses <- rfugw_masses(fit)
provenance <- rfugw_provenance(fit)
list(
  input_mass = c(
    source = masses$source$input_mass,
    target = masses$target$input_mass
  ),
  transported_mass = masses$transported_mass,
  estimand = provenance$estimand,
  backend = provenance$backend,
  objective_units = provenance$objective$units
)
#> $input_mass
#> source target 
#>      1      1 
#> 
#> $transported_mass
#> [1] 1
#> 
#> $estimand
#> [1] "balanced_entropic_linear_ot"
#> 
#> $backend
#> [1] "cpp_scaling"
#> 
#> $objective_units
#> [1] "cost_times_probability_mass"
```

## How is solver state reused safely?

[`rfugw_state()`](https://bbuchsbaum.github.io/rfugw/reference/rfugw_state.md)
returns opaque versioned state only from a certified result. Balanced
state contains canonical dual potentials, so the next solve can change
tolerance, epsilon, or switch between scaling and log backends. It is an
initialization—not a cached result—and the continuation must earn a new
certificate.

``` r

state <- rfugw_state(fit)
tighter_problem <- transport_problem_sinkhorn(
  cost,
  p = source_mass,
  q = target_mass,
  epsilon = 0.55,
  method = "log",
  tol = 1e-10,
  mass_policy = "probability"
)
continued <- transport_solve(tighter_problem, init_state = state)

data.frame(
  status = rfugw_status(continued),
  warm_started = continued$warm_started,
  state_accepted = continued$warm_start_accepted,
  iterations = continued$iterations,
  residual = continued$residual
)
#>      status warm_started state_accepted iterations   residual
#> 1 converged         TRUE           TRUE        110 2.9442e-11
```

Invalid dimensions, nonfinite potentials, changed support, incompatible
state types, and unsafe scaling starts fail before computation. A client
that wants an explicit fallback can request `state_policy = "cold"`; the
result then contains a structured rejection code and does not claim it
was warm-started. Partial-Sinkhorn state is stricter than balanced state
because its Dykstra corrections are bound to the identical problem and
effective backend.

State survives [`saveRDS()`](https://rdrr.io/r/base/readRDS.html) /
[`readRDS()`](https://rdrr.io/r/base/readRDS.html). The retained
fresh-process installed check is
`inst/bench/verify_solver_protocol_state_roundtrip.R`.

## Does the same boundary preserve a second formulation?

Exact transport uses a different constructor and objective certificate.
It does not accept iterative state.

``` r

exact_problem <- transport_problem_emd(
  cost,
  p = source_mass,
  q = target_mass,
  mass_policy = "probability"
)
exact <- transport_solve(exact_problem)
data.frame(
  estimand = rfugw_provenance(exact)$estimand,
  status = rfugw_status(exact),
  value = rfugw_value(exact),
  primal_dual_gap = exact$duality_gap
)
#>                   estimand    status    value primal_dual_gap
#> 1 exact_balanced_linear_ot converged 3.921256    4.440892e-16
```

The formulation-specific constructor is important: exact OT, partial OT,
KL-unbalanced OT, GW, FGW, and FUGW do not become equivalent because
their results share accessors.

## What does protocol v1 support?

``` r

transport_capabilities()[, c(
  "estimand", "protocol_constructor", "mass_policy", "warm_state",
  "plan_representation", "relational_costs"
)]
#>                                                    estimand
#> 1                               balanced_entropic_linear_ot
#> 2                                       sinkhorn_divergence
#> 3                                  exact_balanced_linear_ot
#> 4                        exact_fixed_mass_partial_linear_ot
#> 5                     entropic_fixed_mass_partial_linear_ot
#> 6                 penalized_variable_mass_partial_linear_ot
#> 7                          kl_unbalanced_entropic_linear_ot
#> 8             translation_invariant_kl_unbalanced_linear_ot
#> 9                             classical_wasserstein_summary
#> 10             fixed_support_wasserstein_barycenter_weights
#> 11                                       gromov_wasserstein
#> 12                                 fused_gromov_wasserstein
#> 13                    fixed_mass_partial_gromov_wasserstein
#> 14              fixed_mass_partial_fused_gromov_wasserstein
#> 15 penalized_variable_mass_partial_fused_gromov_wasserstein
#> 16                           semirelaxed_gromov_wasserstein
#> 17                     semirelaxed_fused_gromov_wasserstein
#> 18                      fused_unbalanced_gromov_wasserstein
#> 19                              across_spaces_unbalanced_ot
#> 20                          unbalanced_co_optimal_transport
#> 21              fixed_support_gromov_wasserstein_barycenter
#> 22                               multi_collection_alignment
#> 23                               sampled_gromov_wasserstein
#> 24              dense_plan_gromov_wasserstein_approximation
#>    protocol_constructor                               mass_policy warm_state
#> 1                  TRUE                     probability|normalize       TRUE
#> 2                 FALSE                               probability      FALSE
#> 3                  TRUE                     probability|normalize      FALSE
#> 4                  TRUE                     probability|normalize      FALSE
#> 5                  TRUE                     probability|normalize       TRUE
#> 6                 FALSE                            finite_measure      FALSE
#> 7                  TRUE finite_measure|joint|separate_probability      FALSE
#> 8                 FALSE                            finite_measure      FALSE
#> 9                 FALSE                               probability      FALSE
#> 10                FALSE                               probability       TRUE
#> 11                FALSE                               probability      FALSE
#> 12                FALSE                               probability      FALSE
#> 13                FALSE                               probability      FALSE
#> 14                FALSE                               probability      FALSE
#> 15                FALSE                            finite_measure      FALSE
#> 16                FALSE                  fixed_source_probability      FALSE
#> 17                FALSE                  fixed_source_probability      FALSE
#> 18                FALSE                            finite_measure      FALSE
#> 19                FALSE                            finite_measure      FALSE
#> 20                FALSE                            finite_measure      FALSE
#> 21                FALSE                               probability       TRUE
#> 22                FALSE                               probability       TRUE
#> 23                FALSE                               probability      FALSE
#> 24                FALSE                               probability      FALSE
#>            plan_representation relational_costs
#> 1                        dense            FALSE
#> 2                        dense            FALSE
#> 3                        dense            FALSE
#> 4                        dense            FALSE
#> 5                        dense            FALSE
#> 6                        dense            FALSE
#> 7                        dense            FALSE
#> 8  implicit|sparse_edges|dense            FALSE
#> 9                        dense            FALSE
#> 10                       dense            FALSE
#> 11                       dense             TRUE
#> 12                       dense             TRUE
#> 13                       dense             TRUE
#> 14                       dense             TRUE
#> 15                       dense             TRUE
#> 16                       dense             TRUE
#> 17                       dense             TRUE
#> 18                       dense             TRUE
#> 19                       dense             TRUE
#> 20                       dense             TRUE
#> 21                       dense             TRUE
#> 22                       dense             TRUE
#> 23                       dense             TRUE
#> 24 dense_input_low_rank_output             TRUE
```

A false `protocol_constructor` means the distinct named public solver
remains available but is not wrapped by protocol v1. Sparse and implicit
couplings are consumed through the transport-plan/operator API;
materialization stays an explicit action. Application-specific intensive
versus extensive normalization remains in the downstream library.

The installed generic fixture is
`inst/extdata/clients/generic_solver_protocol_fixture.R`. A DKGE-shaped
fixture and differential test reproduce its certified log-Sinkhorn plan,
objective, marginal error, failure behavior, and tighter-tolerance reuse
using only exported rfugw APIs on the client side. Manifoldalign
fixtures exercise the same principle for OT-Procrustes and sparse
translation-invariant UOT.
