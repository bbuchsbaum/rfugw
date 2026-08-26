# Getting started with rfugw

Suppose you have two collections of objects—brain regions from two
subjects, for example—but their coordinates are not directly comparable.
You do know the distances *within* each collection, and you may have
features measured in both. `rfugw` turns those inputs into a transport
plan: a matrix of soft correspondences between the collections.

This vignette follows one small alignment from data to result. It
introduces the default entropic fused Gromov–Wasserstein (FGW) workflow
first, then shows how to omit features and how to inspect the result
safely.

## What data do you need?

Each row below is one object. The target is a reordered, rotated version
of the source, so ordinary coordinate distance between the two matrices
would be misleading.

``` r

library(rfugw)

source_xy <- rbind(A = c(0, 0), B = c(1, 0), C = c(1.4, 0.8),
                   D = c(0.7, 1.5), E = c(-0.3, 0.9), F = c(0.2, 0.4))
target_order <- c(4, 1, 6, 3, 5, 2)
```

``` r

theta <- pi / 3
rotation <- matrix(c(cos(theta), -sin(theta),
                     sin(theta),  cos(theta)), 2, 2)
target_xy <- source_xy[target_order, ] %*% rotation
target_xy <- sweep(target_xy, 2, c(4, -2), "+")
rownames(target_xy) <- rownames(source_xy)[target_order]
```

GW and FGW use three cost matrices:

- `C1` is the square structure cost within the source;
- `C2` is the square structure cost within the target;
- `M` is the rectangular feature cost from source rows to target rows.

Pairwise distances are a common choice for `C1` and `C2`. Here each
object also has a shared scalar feature; `M[i, j]` is the absolute
feature difference for a possible match.

``` r

C1 <- as.matrix(dist(source_xy))
C2 <- as.matrix(dist(target_xy))
C1 <- C1 / max(C1); C2 <- C2 / max(C2)

source_feature <- c(0.05, 0.25, 0.50, 0.70, 0.90, 0.38)
target_feature <- source_feature[target_order]
M <- abs(outer(source_feature, target_feature, "-"))
```

Scaling comparable costs to similar ranges makes `alpha` easier to
interpret. Cost matrices must be finite; structure matrices must be
square, while `M` must have `nrow(C1)` rows and `nrow(C2)` columns.

## How do you fit an alignment?

Use
[`fgw_entropic()`](https://bbuchsbaum.github.io/rfugw/reference/fgw_entropic.md)
when both structure and shared features should influence the match.
`alpha` is the structure share: `alpha = 1` uses structure alone, and
`alpha = 0` uses features alone. Positive `epsilon` smooths the plan;
larger values generally spread mass more broadly.

``` r

fit <- fgw_entropic(
  M, C1, C2, alpha = 0.5,
  epsilon = 0.08, sinkhorn_tol = 1e-8
)
fit
#> <rfugw_result>
#>   formulation: fgw
#>   backend:     cpp_strict_double
#>   status:      converged
#>   value:       0.0247043
#>   iterations:  16 / 1000
#>   residual:    9.018e-10
#>   row/col res: 1.613e-10 / 5.551e-17
```

The compact print method reports the formulation, backend, status,
objective, iteration count, and residual without printing the whole
plan. Optimization is non-convex, so this is a solver result, not a
certificate of the global optimum.

## What does the plan say?

Rows of the plan correspond to source objects and columns to target
objects. A balanced plan has the requested source and target weights as
its row and column sums. For a compact view, report the strongest target
for each source.

``` r

best_column <- apply(rfugw_plan(fit), 1, which.max)
matches <- data.frame(
  source = rownames(source_xy),
  target = rownames(target_xy)[best_column]
)
matches
#>   source target
#> 1      A      A
#> 2      B      B
#> 3      C      C
#> 4      D      D
#> 5      E      E
#> 6      F      F
```

In this constructed example, the recovered target labels agree with the
source labels. In real analyses there may be no one-to-one answer: keep
the full soft plan when its uncertainty or mass splitting matters.

## What if you have structure but no shared features?

Use
[`entropic_gromov_wasserstein()`](https://bbuchsbaum.github.io/rfugw/reference/entropic_gromov_wasserstein.md)
when only within-domain structure is meaningful. Its input is just the
two structure matrices.

``` r

gw_fit <- entropic_gromov_wasserstein(
  C1, C2, epsilon = 0.08, sinkhorn_tol = 1e-8
)
summary(gw_fit)
#> rfugw result summary
#>   fgw via cpp_strict_double: converged
#>   value=0.0251479  iterations=35/1000  plan=6 x 6
```

The result has the same central contract:
[`rfugw_plan()`](https://bbuchsbaum.github.io/rfugw/reference/rfugw_plan.md)
extracts the coupling,
[`rfugw_value()`](https://bbuchsbaum.github.io/rfugw/reference/rfugw_value.md)
extracts the documented objective, and
[`rfugw_status()`](https://bbuchsbaum.github.io/rfugw/reference/rfugw_status.md)
reports the stopping status.

Unregularized
[`gromov_wasserstein()`](https://bbuchsbaum.github.io/rfugw/reference/gromov_wasserstein.md)
and
[`fgw_exact_cg()`](https://bbuchsbaum.github.io/rfugw/reference/fgw_exact_cg.md)
are alternatives when you specifically want conditional-gradient
optimization with exact linear-OT directions. They can produce sharper
plans, but the outer GW/FGW problem remains non-convex.

## Which result fields should you trust?

Prefer the stable accessors over reaching into solver-specific lists:

- [`rfugw_plan()`](https://bbuchsbaum.github.io/rfugw/reference/rfugw_plan.md)
  returns the transport plan;
- [`rfugw_value()`](https://bbuchsbaum.github.io/rfugw/reference/rfugw_value.md)
  returns the documented objective;
- [`rfugw_status()`](https://bbuchsbaum.github.io/rfugw/reference/rfugw_status.md)
  returns the solver status;
- [`rfugw_residuals()`](https://bbuchsbaum.github.io/rfugw/reference/rfugw_residuals.md)
  returns feasibility and stopping diagnostics.

Always inspect status and residuals before interpreting a plan. A finite
plan can still come from a run that reached an iteration limit. Also
compare objective values only when the formulation, inputs, weights, and
regularization are the same; an entropically regularized value is not
interchangeable with an unregularized GW value.

Several solver families provide a documented scalar companion, such as
[`fgw_entropic2()`](https://bbuchsbaum.github.io/rfugw/reference/fgw_entropic2.md)
and
[`gromov_wasserstein2()`](https://bbuchsbaum.github.io/rfugw/reference/gromov_wasserstein2.md).
These return only the objective and are convenient when the coupling and
diagnostics are not needed.

## Where should you go next?

- [`vignette("solver-guide")`](https://bbuchsbaum.github.io/rfugw/articles/solver-guide.md)
  explains how marginal assumptions determine the solver family.
- [`vignette("linear-ot")`](https://bbuchsbaum.github.io/rfugw/articles/linear-ot.md)
  covers transport when you already have a direct cross-domain cost
  matrix.
- [`vignette("barycenters")`](https://bbuchsbaum.github.io/rfugw/articles/barycenters.md)
  constructs representative metric spaces.
- [`vignette("multiset-alignment")`](https://bbuchsbaum.github.io/rfugw/articles/multiset-alignment.md)
  aligns several collections to a template.
