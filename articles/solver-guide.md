# Choosing the right GW solver

Choosing a solver starts with the correspondence you want. Do both
distributions have fixed mass? Should every object participate? Do you
have features that can be compared across domains? The answers determine
the model; regularization and numerical controls come afterward.

This guide fits several formulations to the same two small metric
spaces. Each returns a soft correspondence, but the plans obey different
marginal constraints and their objective values are not directly
comparable.

Different numbers of source and target objects do not, by themselves,
require partial or unbalanced transport. Balanced transport handles
rectangular plans; choose a relaxed formulation only when the mass
assumptions call for one.

## What are we aligning?

The source and target contain different numbers of objects. Each has an
internal distance matrix, and each object has two features measured in a
common space.

``` r

library(rfugw)

source_xy <- rbind(c(0, 0), c(1, 0), c(1.2, 0.8), c(0.4, 1.3), c(-0.3, 0.6))
target_xy <- rbind(c(2.0, 1.0), c(2.5, 1.8), c(1.9, 2.4),
                   c(1.2, 2.0), c(1.0, 1.2), c(2.8, 2.6))
C1 <- as.matrix(dist(source_xy)); C1 <- C1 / max(C1)
C2 <- as.matrix(dist(target_xy)); C2 <- C2 / max(C2)
```

``` r

F1 <- cbind(source_xy[, 1] + source_xy[, 2], source_xy[, 2])
F2 <- cbind(target_xy[, 1] + target_xy[, 2] - 3, target_xy[, 2] - 1)
M <- as.matrix(dist(rbind(F1, F2)))[1:5, 6:11]
M <- M / max(M)
p <- rep(1 / nrow(C1), nrow(C1))
q <- rep(1 / nrow(C2), nrow(C2))
```

## Which family matches your assumptions?

Use the narrowest formulation that expresses the scientific problem:

| Your data and assumptions | Solver family |
|----|----|
| A direct cost between objects; no relational structure | [`ot_sinkhorn()`](https://bbuchsbaum.github.io/rfugw/reference/ot_sinkhorn.md) or [`ot_emd()`](https://bbuchsbaum.github.io/rfugw/reference/ot_emd.md) |
| Structure only; both marginals fixed | [`entropic_gromov_wasserstein()`](https://bbuchsbaum.github.io/rfugw/reference/entropic_gromov_wasserstein.md) |
| Structure plus shared features; both marginals fixed | [`fgw_entropic()`](https://bbuchsbaum.github.io/rfugw/reference/fgw_entropic.md) |
| Only a specified amount of mass should match | `entropic_partial_*()` or `partial_*()` |
| Unmatched mass has a linear cost and matched mass should be optimized | [`penalized_partial_fused_gromov_wasserstein()`](https://bbuchsbaum.github.io/rfugw/reference/penalized_partial_fused_gromov_wasserstein.md) |
| Source weights fixed; target weights learned | `entropic_semirelaxed_*()` or `semirelaxed_*()` |
| Both marginals may depart from reference weights | [`fugw_kl()`](https://bbuchsbaum.github.io/rfugw/reference/fugw_kl.md) |

Within the balanced GW and FGW families, the unregularized counterparts
[`gromov_wasserstein()`](https://bbuchsbaum.github.io/rfugw/reference/gromov_wasserstein.md)
and
[`fgw_exact_cg()`](https://bbuchsbaum.github.io/rfugw/reference/fgw_exact_cg.md)
use conditional gradient with exact linear-OT directions. The entropic
functions are usually the practical first attempt because their positive
regularization produces smoother plans and avoids an exact transport
solve at every outer iteration.

## When are both marginals fixed?

Balanced FGW requires every row and column sum to match the supplied
weights. Use it when both distributions represent complete probability
measures and all mass should be assigned.

``` r

balanced <- fgw_entropic(
  M, C1, C2, p = p, q = q,
  alpha = 0.5, epsilon = 0.08
)
balanced
#> <rfugw_result>
#>   formulation: fgw
#>   backend:     cpp_strict_double
#>   status:      converged
#>   value:       0.183826
#>   iterations:  30 / 1000
#>   residual:    7.624e-10
#>   row/col res: 4.161e-12 / 5.551e-17
```

Set `alpha = 1` or call `entropic_gromov_wasserstein(C1, C2)` when no
cross-domain features are available. Do not compare objective values
across different `alpha` or `epsilon` settings as though they were the
same quantity.

## When should only part of the mass match?

Partial GW is for overlap known through a total transported mass `m`.
Rows and columns may retain unmatched mass, but the plan must transport
exactly `m` and must not exceed either marginal.

``` r

partial <- entropic_partial_gromov_wasserstein(
  C1, C2, p = p, q = q,
  m = 0.7, reg = 0.10, log = TRUE
)
data.frame(
  status = rfugw_status(partial),
  objective = rfugw_value(partial),
  transported_mass = sum(rfugw_plan(partial))
)
#>      status  objective transported_mass
#> 1 converged 0.06493918              0.7
```

The exact
[`partial_gromov_wasserstein()`](https://bbuchsbaum.github.io/rfugw/reference/partial_gromov_wasserstein.md)
and
[`partial_fused_gromov_wasserstein()`](https://bbuchsbaum.github.io/rfugw/reference/partial_fused_gromov_wasserstein.md)
default to the built-in `lp_solver = "cpp_transport"`. The optional
`"lp_transport"` and `"lp_matrix"` backends require the suggested
package `lpSolve`.

Their entropic counterparts call the same certified
[`ot_partial_sinkhorn()`](https://bbuchsbaum.github.io/rfugw/reference/ot_partial_sinkhorn.md)
primitive available for direct-cost partial OT. Each outer projection
must pass fixed-mass capacity, objective, KKT, and duality certificates.
`method = "auto"` chooses bounded scaling or genuine log-domain Dykstra
from the actual linearized cost at every outer iteration; explicit
unsafe scaling fails.

Partial transport does not infer the scientifically correct overlap for
you: choose `m` from the meaning of the data, or treat it as a
sensitivity parameter and report how conclusions change.

When overlap is itself the estimand and unmatched source and target mass
each have a linear cost, use
[`penalized_partial_fused_gromov_wasserstein()`](https://bbuchsbaum.github.io/rfugw/reference/penalized_partial_fused_gromov_wasserstein.md).
This does not fix `m`: every iterate is a subcoupling, and the solver
chooses its mass while balancing feature, structure, and discard terms.

``` r

penalized_low <- penalized_partial_fused_gromov_wasserstein(
  M, C1, C2, discard_penalty = 0.05,
  p = p, q = q, alpha = 0.5
)
penalized_high <- penalized_partial_fused_gromov_wasserstein(
  M, C1, C2, discard_penalty = 0.4,
  p = p, q = q, alpha = 0.5
)
data.frame(
  discard_penalty = c(0.05, 0.4),
  transported_mass = c(
    penalized_low$transported_mass,
    penalized_high$transported_mass
  ),
  feature_term = c(penalized_low$feature_term, penalized_high$feature_term),
  structure_term = c(
    penalized_low$structure_term,
    penalized_high$structure_term
  ),
  discard_term = c(
    penalized_low$discard_penalty_term,
    penalized_high$discard_penalty_term
  )
)
#>   discard_penalty transported_mass feature_term structure_term discard_term
#> 1            0.05              0.5   0.02382383    0.005128658         0.05
#> 2            0.40              1.0   0.11471844    0.053865379         0.00
```

Each transported unit avoids two discard penalties, so the larger
penalty has the derived tendency toward more mass. The displayed
direction is a controlled fixture, not a global monotonicity theorem for
this nonconvex problem. A converged result certifies a stationary point
through a Frank-Wolfe gap; it does not certify the global FGW optimum.
[`fugw_kl()`](https://bbuchsbaum.github.io/rfugw/reference/fugw_kl.md)
below is a different estimand: it penalizes marginal changes with
generalized KL rather than a linear unmatched-mass term.

## When is the target distribution unknown?

Semirelaxed GW fixes the source marginal and learns the target marginal.
This fits problems where every source object must be represented but the
appropriate mass allocation over target objects is unknown.

``` r

semirelaxed <- entropic_semirelaxed_gromov_wasserstein(
  C1, C2, p = p, epsilon = 0.08, tol = 2e-8
)
data.frame(
  status = rfugw_status(semirelaxed),
  objective = rfugw_value(semirelaxed),
  target_mass = sum(rfugw_plan(semirelaxed))
)
#>      status  objective target_mass
#> 1 converged 0.04252806           1
```

This is not the same as partial transport: all source mass is
transported, but the solver decides how it is distributed across target
columns.

## When may both marginals move?

[`fugw_kl()`](https://bbuchsbaum.github.io/rfugw/reference/fugw_kl.md)
relaxes both marginal constraints using KL penalties. It is useful when
reference weights are meaningful but exact equality is too strong, for
example when collections contain unequal or unreliable mass.

``` r

unbalanced <- fugw_kl(
  C1, C2, wx = p, wy = q, M = M,
  alpha = 0.5, epsilon = 0.5,
  reg_marginals = c(1, 1),
  max_iter = 200L, tol = 1e-6,
  max_iter_ot = 1000L, tol_ot = 1e-6
)
unbalanced
#> <rfugw_result>
#>   formulation: fugw_kl
#>   backend:     cpp_double
#>   status:      converged
#>   value:       0.424622
#>   iterations:  5 / 200
#>   residual:    5.567e-07
```

Larger `reg_marginals` penalizes departure from `wx` and `wy` more
strongly; smaller values permit more creation, destruction, or
redistribution of mass. There is no universal setting: inspect the
transported mass and marginal deviations, and assess sensitivity on the
scale relevant to your application. The relatively smooth `epsilon`
above keeps this small teaching example’s nested solve well conditioned;
it is not a recommended value for every data set.

## How do regularized and exact solvers differ?

Entropic solvers require a positive `epsilon` (or `reg` for the partial
API). Regularization smooths the coupling and changes the optimized
objective. A smaller value is not automatically better: it may make the
inner transport problem harder and can expose numerical limitations.

Exact conditional-gradient variants optimize the unregularized objective
and often yield sparser directions. They still solve a non-convex outer
problem, so neither family certifies a global optimum. For asymmetric
structure costs, leave `symmetric = NULL` so the solver detects the
general path; setting `symmetric = TRUE` requires genuinely symmetric
inputs.

Start with documented defaults, normalize commensurate costs, and change
one control at a time. Use `sinkhorn_method = "log"` or stricter
precision only when diagnostics or the cost scale justify it; those are
numerical choices, not different scientific models.

## What should you check before using a plan?

For result objects, use the common interface:

``` r

rfugw_status(balanced)
#> [1] "converged"
rfugw_value(balanced)
#> [1] 0.1838261
rfugw_residuals(balanced)[c("feasible", "objective_consistent")]
#> $feasible
#> [1] TRUE
#> 
#> $objective_consistent
#> [1] TRUE
```

Check four things:

1.  The plan is finite, nonnegative, and non-degenerate.
2.  Its marginals satisfy the formulation you selected.
3.  Status and residuals support convergence rather than merely
    exhaustion of the iteration budget.
4.  The objective is interpreted only within its exact formulation and
    regularization convention.

[`ot_validate_plan()`](https://bbuchsbaum.github.io/rfugw/reference/ot_validate_plan.md)
checks balanced and partial marginal contracts. Unbalanced plans have no
fixed marginal target, so inspect their total mass, row and column sums,
and the diagnostics returned by
[`rfugw_residuals()`](https://bbuchsbaum.github.io/rfugw/reference/rfugw_residuals.md).

## Which APIs are experimental?

[`sampled_gromov_wasserstein()`](https://bbuchsbaum.github.io/rfugw/reference/sampled_gromov_wasserstein.md),
[`sampled_gromov_wasserstein_coords()`](https://bbuchsbaum.github.io/rfugw/reference/sampled_gromov_wasserstein_coords.md),
[`sampled_gw_from_graphs()`](https://bbuchsbaum.github.io/rfugw/reference/sampled_gw_from_graphs.md),
and
[`dense_gromov_wasserstein_plan_svd()`](https://bbuchsbaum.github.io/rfugw/reference/dense_gromov_wasserstein_plan_svd.md)
are experimental. Their names and signatures may change, and their
result status does not claim convergence certification. In particular,
the last function materializes two dense structure costs and a dense
plan before factorizing that plan. It is not an end-to-end low-rank or
low-memory GW solver. The older
[`lowrank_gromov_wasserstein_samples()`](https://bbuchsbaum.github.io/rfugw/reference/lowrank_gromov_wasserstein_samples.md)
name is deprecated; its formerly ignored POT-shaped factorized-cost,
Dykstra, seed, and warning arguments now error when supplied.

Use these only when their documented envelope matches your task. See
`inst/bench/sampled-budget-curves.md` for the package’s current evidence
and `inst/solver-contract.md` for detailed solver contracts.

## Where should you go next?

- [`vignette("rfugw")`](https://bbuchsbaum.github.io/rfugw/articles/rfugw.md)
  gives the shortest end-to-end FGW workflow.
- [`vignette("linear-ot")`](https://bbuchsbaum.github.io/rfugw/articles/linear-ot.md)
  covers direct-cost balanced and unbalanced OT.
- [`vignette("barycenters")`](https://bbuchsbaum.github.io/rfugw/articles/barycenters.md)
  covers representative metric spaces.
- [`vignette("multiset-alignment")`](https://bbuchsbaum.github.io/rfugw/articles/multiset-alignment.md)
  covers subject-to-template alignment.
