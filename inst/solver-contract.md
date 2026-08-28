# rfugw solver contract and compatibility policy

This document is the supported contract for rfugw 0.1. Public solvers, aliases,
parameters, and diagnostics are interpreted against it. Implementation gaps
are either fixed in this release or linked to a tracked ticket.

Related tickets: `bd-01M05QY311B0WZEV6JD3F1FN90` (this contract),
`bd-01M05QY3VJTSR30MN5PGP7T1Y0` (diagnostics),
`bd-01M05QY6A7T8RX4D28A4D8XCZ7` (result class).

## Scope

rfugw is a focused optimal-transport library: linear OT, Gromov-Wasserstein
(GW), fused GW (FGW), unbalanced alignment, and a small set of scalable
variants. It does not pursue POT breadth.

Current global limits:

- Inner structural loss is `square_loss` only.
- There is no Python runtime dependency.
- Numerical claims are guarded by tests and, where advertised, POT fixtures.
- Symmetric and asymmetric partial GW/FGW use the same canonical square-loss
  tensor module. The general path is certified by O(n^4) loss/gradient,
  directional-derivative, permutation, exact-LP, and entropic-projection
  oracles. Passing `symmetric = TRUE` never overrides validation of asymmetric
  inputs.
- `inst/numerical-trust-charter.md` defines the test families and evidence
  required for a path to be classified as supported.
- `inst/numerical-path-matrix.csv` is the executable inventory of advertised
  backend, precision, adapter, start, threading, and approximation paths.

## Versioned solver-client protocol

The protocol-v1 boundary is for downstream libraries that need to construct a
transport estimand, solve it, retain certificates and provenance, and reuse
valid numerical state without reading solver-specific list fields.

1. Call a formulation-specific constructor. The current constructors are
   `transport_problem_sinkhorn()`, `transport_problem_emd()`,
   `transport_problem_partial_emd()`,
   `transport_problem_partial_sinkhorn()`, and
   `transport_problem_unbalanced()`. They do not accept a generic formulation
   string.
2. Call `transport_solve()`. The complete problem and any state are
   revalidated before a numerical backend is entered.
3. Consume `rfugw_plan()`, `rfugw_value()`, `rfugw_status()`, and
   `rfugw_residuals()`. `rfugw_problem()`, `rfugw_masses()`, and
   `rfugw_provenance()` expose stable cross-solver metadata.
4. For a certified iterative result, call `rfugw_state()`. The returned
   `rfugw_solver_state` is opaque, versioned, and serializable. A new solve
   must still earn its own convergence and feasibility certificate.

All problem objects declare `orientation = "rows_source_columns_target"`, the
scientific estimand, original and effective measures, explicit mass policy,
objective convention and units, and requested controls. Probability problems
require `mass_policy = "probability"` to verify unit sums or
`"normalize"` to request and record separate normalization. KL-unbalanced
problems require one of `"finite_measure"`, `"joint"`, or
`"separate_probability"`; the protocol has no ambiguous compatibility
default.

Balanced Sinkhorn state stores canonical dual potentials. It is compatible
across scaling, log, and auto dispatch and with positive-epsilon continuation,
subject to dimensions, active support, finiteness, and scaling-overflow
validation. Partial-Sinkhorn state stores Dykstra corrections and is bound to
the identical cost, measures, transported mass, epsilon, support, and
effective backend. Exact and current KL-unbalanced protocol problems do not
accept iterative state.

`transport_solve(..., state_policy = "error")` is fail-closed. With
`state_policy = "cold"`, an invalid state is deliberately rejected and the
result records `warm_started = FALSE`, `warm_start_accepted = FALSE`, and a
structured `warm_start_rejection` with `code` and `message`. State is never
treated as a cached answer or inherited convergence claim.

`transport_capabilities()` is the executable capability matrix. It keeps
balanced, exact, fixed-mass partial, variable-mass partial, KL-unbalanced,
translation-invariant UOT, GW, FGW, and FUGW as distinct rows. A false
`protocol_constructor` means the formulation remains available through its
named public solver and stable result accessors but has not been collapsed
into protocol v1. Sparse and implicit support remains explicit at the
transport-plan boundary; the protocol does not silently materialize it.

The installed generic client fixture in
`inst/extdata/clients/generic_solver_protocol_fixture.R` uses exported rfugw
APIs only. DKGE and manifoldalign are conformance case studies, not owners of
the contract: they retain their cost construction and intensive/extensive
scientific interpretation. The retained continuation and fresh-process RDS
evidence is in `inst/bench/solver-protocol-benchmark.md`.

## Compatibility policy

1. **Primary names are stable.** Native rfugw names (`fgw_entropic`,
   `fgw_exact_cg`, `fugw_kl`, and the linear-OT primitives added in Phase C)
   are the supported API.
2. **POT-style names are aliases.** They call the same backend and share
   diagnostics. They exist for familiarity, not as a second implementation.
3. **`*2` names return a scalar objective.** They call the corresponding
   solver and extract the documented objective field. They do not change
   formulation or stopping rules.
4. **Accepted parameters are operational, ignored, or unsupported.**
   - *Operational*: used by the solver.
   - *Explicitly ignored*: accepted for POT-shaped signatures, documented as
     ignored, and never silently change the result.
   - *Unsupported*: rejected before computation with the supported
     alternatives named in the error.
5. **Unregularized does not mean globally optimal.** Names such as
   `fgw_exact_cg`, `gromov_wasserstein`, and `fused_gromov_wasserstein` mean
   *unregularized conditional gradient with an exact linear-OT subproblem*.
   The outer GW/FGW problem is non-convex. A returned plan is a stationary
   point of that procedure, not a certified global minimizer.
6. **Weights are probability measures unless a solver explicitly accepts finite measures.**
   Most finite nonnegative weights with positive total mass are renormalized
   to sum 1. That renormalization is part of their contract
   (`bd-01M05QY426X8F52TTBEM34AXEY`). `ot_sinkhorn_unbalanced()` and
   `ot_sinkhorn_unbalanced_ti()` are exceptions. The former's
   `normalization = "none"` preserves finite nonnegative measures
   (including a zero measure), `"joint"` applies one common scale while
   retaining the source/target mass ratio, and `"separate"` retains the pre-0.1
   probability-normalized behavior. The backward-compatible 0.1 default is
   `"separate"`; callers using absolute mass must opt into `"none"` explicitly
   during this transition. The TI solver always preserves the supplied finite
   measures and requires positive total mass on both sides.
7. **Breaking changes require a NEWS entry** and, after 0.1, a minor version
   bump. Experimental APIs may change in a patch if labeled experimental.

## Export classification

### Primary solvers

| Function | Family | Maturity |
|---|---|---|
| `ot_sinkhorn` | Balanced entropic linear OT | Flagship |
| `ot_sinkhorn_divergence` | Debiased three-solve regularized OT divergence | Supported |
| `ot_emd` | Exact balanced linear OT | Flagship |
| `ot_partial_emd` | Exact nonnegative-cost partial linear OT | Supported |
| `ot_partial_sinkhorn` | Entropy-regularized fixed-mass partial linear OT | Supported |
| `ot_partial_penalized` | Exact variable-mass partial OT with linear discard penalty | Supported |
| `ot_sinkhorn_unbalanced` | KL-unbalanced entropic linear OT | Flagship |
| `ot_sinkhorn_unbalanced_ti` | Translation-invariant dense/sparse KL-unbalanced OT | Supported |
| `ot_wasserstein_cost`, `ot_wasserstein_distance` | Explicit classical Wasserstein summaries over exact or labeled entropic plans | Supported |
| `ot_barycenter_weights` | Fixed-support linear-Wasserstein barycenter weights; exact joint LP or regularized semi-debiased Sinkhorn objective | Supported |
| `fgw_entropic` | Entropic FGW | Flagship |
| `fgw_exact_cg` | Unregularized FGW (CG + LP) | Flagship |
| `fugw_kl` | Fused unbalanced GW (KL) | Flagship |
| `partial_gromov_wasserstein` | Partial GW (CG + partial LP) | Supported |
| `partial_fused_gromov_wasserstein` | Partial FGW | Supported |
| `penalized_partial_fused_gromov_wasserstein` | Variable-mass partial FGW with linear discard penalty | Supported |
| `entropic_partial_gromov_wasserstein` | Entropic partial GW | Supported |
| `entropic_partial_fused_gromov_wasserstein` | Entropic partial FGW | Supported |
| `semirelaxed_gromov_wasserstein` | Unregularized semirelaxed GW | Supported |
| `semirelaxed_fused_gromov_wasserstein` | Unregularized semirelaxed FGW | Supported |
| `entropic_semirelaxed_gromov_wasserstein` | Entropic semirelaxed GW | Supported |
| `entropic_semirelaxed_fused_gromov_wasserstein` | Entropic semirelaxed FGW | Supported |
| `fused_unbalanced_across_spaces_divergence` | Across-spaces unbalanced OT | Supported |
| `unbalanced_co_optimal_transport` | UCOOT (`reg_type = "independent"`) | Supported |
| `fgw_barycenters` | Fixed-support FGW barycenters | Supported |

### POT-name aliases

| Alias | Primary |
|---|---|
| `entropic_fused_gromov_wasserstein` | `fgw_entropic` |
| `entropic_gromov_wasserstein` | `fgw_entropic` with zero feature cost and `alpha = 1` |
| `fused_gromov_wasserstein` | `fgw_exact_cg` |
| `gromov_wasserstein` | `fgw_exact_cg` with zero feature cost and `alpha = 1` |
| `fused_unbalanced_gromov_wasserstein` | `fugw_kl` |
| `entropic_gromov_barycenters` | `fgw_barycenters` (structure only) |
| `entropic_fused_gromov_barycenters` | `fgw_barycenters` |
| `gromov_barycenters` | `entropic_gromov_barycenters` |

### Scalar `*2` aliases

`fgw_entropic2`, `fugw_kl2`, `entropic_fused_gromov_wasserstein2`,
`entropic_gromov_wasserstein2`, `fused_gromov_wasserstein2`,
`fused_unbalanced_gromov_wasserstein2`, `gromov_wasserstein2`,
`partial_gromov_wasserstein2`, `partial_fused_gromov_wasserstein2`,
`entropic_partial_gromov_wasserstein2`,
`entropic_partial_fused_gromov_wasserstein2`,
`entropic_semirelaxed_gromov_wasserstein2`,
`entropic_semirelaxed_fused_gromov_wasserstein2`,
`semirelaxed_gromov_wasserstein2`, `semirelaxed_fused_gromov_wasserstein2`,
`unbalanced_co_optimal_transport2`.

Each returns the documented objective field of its parent solver.

### Experimental / approximate

These are not flagship 0.1 solvers. The only certified claims are the
path-specific claims below. The admission review in
`inst/scalable-relational-ot-decision.md` promotes no scalable candidate.
Coordinate and graph sampled-GW inputs save structure-input storage but still
materialize a dense plan and at least five coupling-sized work matrices. The
new Moment-FUGW path has dense-oracle, matrix-free representation, four-size
scaling, multiscale-efficacy, partial-overlap image, and planted cortical
held-out-map candidate evidence. Those numerical, representation, geometry,
scientific, and release-evidence layers remain separate. Supported promotion
still requires a reviewed clean same-commit installed artifact and hosted
cross-platform replay.

| Function | Status | Certified envelope |
|---|---|---|
| `fugw_factorized` | Experimental | Solves the same joint-KL two-coupling objective as `fugw_kl()` for affine-bilinear structure and feature costs. Complete-support couplings are retained by TI-Sinkhorn potentials and applied in blocks. Small exact-factor fixtures certify objective components, both plan actions, masses, and `H`/Gram moments against dense `fugw_kl()`. The result separately reports inner UOT, outer stationarity, geometry representation, and complete-support certificates; global optimality is not claimed. Explicit `plan = "dense"` is the only normal dense-coupling allocation boundary. |
| `fugw_multiscale` | Experimental | Runs arbitrary-depth coarse-to-fine `fugw_domain()` hierarchies with per-level epsilon schedules. Parent maps prolong both potential pairs; a conditional-product lift preserves mass and recomputes fine-basis marginal, entropy, `H`, and Gram state. An independent tiny dense lift oracle guards the transfer. Only complete implicit support is implemented; `support = "adaptive"` fails closed pending an omitted-mass or reduced-cost certificate. |
| `sampled_gromov_wasserstein` | Experimental | A full budget `(ns, nt)` is closer to dense entropic GW than a tiny budget such as `(2, 1)`, in square-loss GW and plan Frobenius distance. Intermediate budgets are not certified as monotone. Budgets `< 1` error; source/target counts above `ns`/`nt` warn and clamp. |
| `sampled_gromov_wasserstein_coords` | Experimental | Same tiny-versus-full quality envelope. Inputs scale as `O(n d)` rather than two dense `n x n` structure costs. |
| `sampled_gw_from_graphs` | Experimental | Same envelope after diffusion coordinates. A sparse graph plus `k` embeddings stores less than two dense structure costs. |
| `dense_gromov_wasserstein_plan_svd` | Experimental | Post-hoc SVD of a dense GW plan. The solve first materializes two dense square costs and a dense rectangular plan, so its solve memory is `O(ns^2 + nt^2 + ns*nt)` and is not bounded by the output rank. Relative Frobenius reconstruction error decreases as rank grows up to `min(ns, nt)`. Rank `< 1` errors; rank above `min(ns, nt)` warns and clamps. |
| `lowrank_gromov_wasserstein_samples` | Deprecated | Compatibility wrapper for `dense_gromov_wasserstein_plan_svd()`. It warns on use. Former POT-shaped factorized-cost, Dykstra, seed, and warning parameters were never operational and are rejected whenever supplied. |

### Utilities

Plan validation and independent objectives (`bd-01M05QY759KA7W2VQHYSPEEJ1M`):
`ot_validate_plan`, `ot_linear_cost`, `ot_entropy`, `ot_kl`,
`ot_gw_square`, `ot_fgw_square`, `ot_barycentric_project`.

These evaluators are independent of solver internals. Reported
`ot_dist` / `fgw_dist` values must match `ot_linear_cost` /
`ot_fgw_square` on the returned plan.

`ot_kl(plan, p, q)` is generalized KL from `plan` to `p %o% q`: it
includes `-sum(plan) + sum(p %o% q)`. Zero-over-zero terms contribute zero,
while positive plan mass outside zero reference support returns `Inf`.
`ot_entropy()` remains the separate `sum(plan * log(plan))` functional.

Alignment helpers: `graph_diffusion_coordinates`, `multialign_fit`,
`multialign_make_template`, `multialign_normalize_plan`,
`multialign_project_features`, `multialign_project_matrix`.

### Transport plan and operator representations

`as_transport_plan()` accepts a dense matrix, a Matrix CSC/CSR-style sparse
matrix (canonicalized to CSC storage), or a 1-based edge list. Edge lists are
sorted by source then target; duplicates are summed by default or rejected on
request, and zero entries are removed. Missing rows and columns are valid and
carry zero transported mass. `transport_operator()` adds an implicit form with
mandatory forward and adjoint actions plus declared row/column masses.

The sparse-safe public surface is `transport_plan_shape()`,
`transport_plan_mass()`, `transport_plan_apply()`,
`transport_plan_adjoint()`, `transport_plan_barycentric()`, and
`transport_plan_representation()`. Forward application is `G %*% x`; the
adjoint is `t(G) %*% x`. Neither is a mathematical inverse. Reverse
barycentric projection is the adjoint action normalized by column mass, not a
reverse conditional unless that normalization is explicitly the intended
estimand. Zero-mass rows or columns abstain as `NaN`, zero, or an error; they
are never filled by a uniform match.

`transport_plan_materialize()` is the explicit dense allocation boundary.
Sparse validation, mass, application, barycentric projection, linear cost,
entropy, and product-reference KL operate on stored support. An implicit
operator can be validated from declared masses and applied without entries;
entrywise objectives require the caller to materialize it explicitly.

`transport_plan_prune()` reports removed mass. If it changes an
`rfugw_result`, status becomes `requires_revalidation`, convergence and
feasibility become false, objective certificates are cleared, and the stale
residual becomes infinite. Serialization preserves representation, callbacks,
metadata, and behavior. Existing dense results still return their base matrix
from `rfugw_plan()`; sparse/operator results return their representation unless
`materialize = TRUE` is explicitly requested.

### Translation-invariant sparse KL-UOT

`ot_sinkhorn_unbalanced_ti()` optimizes `<C,plan> + rho1 * KL(plan1|p) +
rho2 * KL(plan2|q) + epsilon * KL(plan|p %o% q)` without normalizing either
finite measure. Dense matrices mean complete support. A sparse Matrix means
its structural nonzeros; an edge list or manifoldalign-style CSR list can
represent zero or negative costs explicitly. Duplicate support edges are
rejected because summing costs has no valid transport meaning.

The default plan is an implicit forward/adjoint operator and does not store the
coupling. `plan = "sparse"` explicitly extracts canonical weighted edges;
`plan = "dense"` explicitly allocates the full plan. The translated `f` and
`g` potentials satisfy marginal KKT identities, while `fbar` and `gbar`
reconstruct the plan and may be shifted by equal and opposite constants.

KL-UOT is primal-feasible through the zero coupling even when exact marginals
cannot fit the sparse graph. The result therefore separates finite-potential
active-support coverage from a global exact balanced-marginal max-flow check.
Unequal total masses legitimately make the latter false without invalidating
UOT convergence. A convergence claim nevertheless requires finite potentials
and plan weights, an independently recomputed TI fixed point and KKT identity,
gauge invariance, and a matching generalized-KL primal/dual objective.

## Common diagnostic definitions

These names have one meaning across flagship solvers. Legacy field names
remain for compatibility.

| Term | Meaning |
|---|---|
| `value` / documented `*_dist` / `*_cost` | Objective of the advertised formulation, evaluated at the returned plan. Unless a field name includes entropy or KL, the value does **not** add the entropic regularizer. |
| `plan` | Coupling matrix or explicit transport-plan/operator representation. Balanced solvers satisfy the documented marginals up to residual. Dense legacy results remain base matrices. |
| `plan_representation`, `certificate_invalidated_by_pruning`, `pruning_lost_mass` | Representation/provenance fields for sparse, implicit, or pruned results. Positive pruning loss always invalidates prior convergence, feasibility, and objective certificates. |
| `iterations` | Outer iterations executed. |
| `inner_iterations` | Inner linear-OT / Sinkhorn iterations, when a nested solver exists. |
| `inner_residual`, `max_inner_residual` | Residual from the final required inner solve and the maximum across required inner solves. A final small value does not erase an earlier uncertified solve. |
| `inner_converged`, `inner_status` | Whether every formulation-required inner solve was certified, plus its final/specific status. Outer convergence cannot be `TRUE` when this field is `FALSE`. |
| `error` | Stopping residual used by that solver (see table below). Prefer `residual` for new code. |
| `residual` | Same numeric quantity as `error`, with a documented unit. |
| `rel_error` | Relative cost change, used by unregularized CG. |
| `row_residual`, `col_residual` | Balanced: `max(abs(rowSums(plan) - p))` and `max(abs(colSums(plan) - q))`. Partial: maximum positive row/column capacity violation. |
| `mass`, `mass_residual` | Returned mass and, for partial solvers, `abs(sum(plan) - requested_mass)`. Partial log results also expose `transported_mass_target` and whether it was defaulted. |
| `discarded_source_mass`, `discarded_target_mass`, `discard_penalty`, `discard_penalty_term` | Variable-mass partial OT reports the mass left on each side, the per-unit penalty, and their combined linear contribution. Increasing `discard_penalty` weakly favors transporting more mass. |
| `frank_wolfe_gap`, `frank_wolfe_gap_tolerance`, `stationarity_consistent` | Penalized partial FGW stationarity certificate. The nonnegative gap compares the returned iterate with a certified exact variable-mass linear direction; it certifies a stationary point, not a global optimum. |
| `line_search_residual`, `line_search_tolerance`, `line_search_consistent` | Maximum disagreement between the direct penalized-FGW objective and its exact one-dimensional quadratic along accepted directions. |
| `feasibility`, `feasibility_residual`, `feasibility_tolerance`, `feasible` | Formulation-specific certificate. Balanced solvers check both marginals; partial solvers check marginal inequalities and transported mass; semirelaxed solvers check the fixed source marginal; unbalanced solvers use the required scaling fixed-point certificate. |
| `normalization`, `original_*_mass`, `effective_*_mass`, `transported_mass` | Finite-measure provenance for KL-unbalanced linear OT. |
| `source_marginal_kl`, `target_marginal_kl`, `plan_product_kl`, `regularized_objective` | Independently recomputed generalized-KL objective terms for KL-unbalanced linear OT. `ot_dist` keeps its legacy meaning as the unregularized transport term. |
| `source_bar`, `target_bar`, `translation`, `source_potential`, `target_potential` | TI-UOT plan-reconstruction potentials, canonical translation, and translated marginal-KKT potentials. Equal/opposite shifts of the barred potentials leave the plan invariant. |
| `fixed_point_residual`, `kkt_residual`, `primal_dual_gap`, `gauge_residual` | Independent TI-UOT certificates. All must pass before convergence is reported. |
| `support_certificate` | Active-support coverage plus a global exact bipartite max-flow diagnostic. Balanced infeasibility is informational for unequal-mass UOT; uncovered active support prevents finite TI potentials. |
| `entropy_convention`, `entropy_reference`, `product_measure_kl`, `regularized_objective` | Balanced Sinkhorn uses `<M,T> + epsilon * KL(T || p %o% q)`. This does not change legacy `ot_dist = <M,T>`. |
| `regularized_dual_objective`, `regularized_duality_gap`, `regularized_duality_gap_tolerance`, `regularized_dual_consistent` | Dual certificate under the same product-reference convention. Canonical plan-reconstruction potentials are shifted by `-epsilon * log(weight)` on active support to obtain the regularized dual potentials. |
| `entropy_minus_one_objective`, `product_reference_offset` | Exposes the exact constant difference between the product-reference convention and `<M,T> + epsilon * sum(T * (log(T) - 1))`; constants are never hidden. |
| `entropy`, `entropy_minus_one`, `weighted_entropy_minus_one`, `partial_sinkhorn_objective` | Fixed-mass partial Sinkhorn uses counting-measure `<M,T> + epsilon * sum(T * (log(T) - 1))`. Its convention is not silently mixed with balanced product-reference KL. |
| `stationarity_residual`, `complementarity_residual`, `dual_feasibility_residual`, `kkt_tolerance` | Partial Sinkhorn capacity-dual, mass-dual, and primal-dual KKT certificate under its declared entropy convention. |
| `component_solves`, `component_status`, `component_values`, `component_residuals`, `component_runtime_seconds` | Full cross/source-self/target-self provenance for Sinkhorn divergence. A finite cross value cannot hide an uncertified self solve. |
| `barycenter_weights`, `coefficients`, `component_results`, `support_self_result` | Fixed-support Wasserstein barycenter output, input-measure coefficients, final certified cross solves, and the regularized common-support self correction. There is no single barycenter coupling. |
| `kkt_residual`, `kkt_tolerance`, `kkt_consistent`, `simplex_residual`, `objective_trace`, `objective_monotone` | Barycenter optimality and descent evidence. Exact mode uses joint-LP weight reduced costs and primal/dual agreement. Regularized mode uses the simplex projected-gradient map and safeguarded Armijo descent. |
| `wasserstein_p_cost`, `wasserstein_p_distance`, `wasserstein_power`, `value_root` | Explicit Wasserstein convention. The p-cost is `<d^p, plan>` and the distance is its `1/p` root. For `p = 1`, no square root is taken. |
| `estimate_kind`, `value_kind`, `value_certified`, `value_certification` | Distinguishes exact Wasserstein values, certified entropic-plan estimates, and explicitly requested uncertified plan inspection. |
| `input_cost_power`, `effective_cost_power`, `input_cost_scale`, `effective_cost_scale`, `metric_certified` | Provenance for raw metric supports or a user-supplied nonnegative cost with declared power. Arbitrary supplied costs do not acquire a metric certificate. |
| `objective_recomputed`, `objective_residual`, `objective_tolerance`, `objective_consistent` | Independent evaluation of the advertised objective, its discrepancy from the reported value, the comparison threshold, and the resulting certificate. Named objective components expose analogous `*_recomputed`, `*_residual`, and `*_consistent` fields. |
| `objective_components_consistent` | Whether every independently checked named objective component is finite and self-consistent. |
| `status` | Normally one of `converged`, `max_iter`, `inner_failure`, `infeasible`, `objective_mismatch`, `stationarity_failure`, `numerical_failure`, `lp_failure`. Exact linear OT preserves its more specific internal failure reason. |
| `converged` | `TRUE` only when `status == "converged"`: the plan is finite and nonnegative, formulation feasibility passes, reported and independently recomputed objectives agree, every required inner solve is certified, and the solver's stopping rule passes. Iteration count alone never establishes success. |
| `termination_reason` | Exact linear OT uses `optimal`, `max_iter`, `disconnected_basis`, `invalid_cycle`, `invalid_step`, `no_leaving_variable`, or `numerical_failure`. |
| `source_potential`, `target_potential` | Exact OT returns certificate duals whose weighted sum is `dual_objective`. Entropic balanced OT returns reusable potentials in the `weighted_source_mean_zero` gauge; together with `M` and `epsilon` they reconstruct the plan. |
| `dual_state`, `initialization`, `initialization_precedence` | Serializable balanced-Sinkhorn state and provenance. `init_duals` takes precedence over a simultaneously supplied valid `init_plan`; every new solve must earn a fresh feasibility certificate. |
| `primal_objective`, `dual_objective`, `duality_gap` | Exact-OT certificate values. The gap is primal minus dual. |
| `min_reduced_cost` | Minimum `M[i,j] - source_potential[i] - target_potential[j]` over nonbasic cells. |
| `feasibility_tolerance`, `reduced_cost_tolerance`, `duality_gap_tolerance` | The actual thresholds used to certify an exact transport result. |
| `requested_tol`, `effective_tol`, `requested_inner_tol`, `effective_inner_tol` | Caller requests and the thresholds actually used. They are always reported on precision-selecting flagship solvers. |
| `requested_precision`, `effective_precision`, `compute_precision` | Requested policy, selected backend policy, and arithmetic reported by the native implementation. |
| `backend_transition`, `automatic_backend_transition` | Names and flags any automatic precision/backend transition. `none` means no transition. |
| `warning_payload` | `NULL` for a certified result; otherwise a structured failure code and message suitable for callers that do not want emitted warnings. |

Flagship solvers return an `rfugw_result` list. `print()` / `summary()` show
diagnostics without dumping the plan. Accessors `rfugw_plan()`,
`rfugw_value()`, `rfugw_status()`, and `rfugw_residuals()` are the stable
downstream API. Legacy fields (`plan`, `fgw_dist`, `error`, ...) remain.
The object is a named list, so `saveRDS()` / `readRDS()` and copies preserve
fields and the S3 class.

`ot_barycenter_weights()` returns `rfugw_barycenter_result` because a
barycenter has one weight vector and several component couplings rather than a
single plan. `rfugw_value()`, `rfugw_status()`, and `rfugw_residuals()` remain
valid; component plans are retained in `component_results`.

For every public solver result, `converged == TRUE` therefore implies
`feasible == TRUE`, `objective_consistent == TRUE`,
`objective_components_consistent == TRUE`, and either no required nested
certificate or `inner_converged == TRUE`. A failed implication is represented
by its specific status and termination reason; it is never repaired by merely
reaching or avoiding an iteration limit.

### Residual units by family

| Family | `error` / `residual` |
|---|---|
| Exact balanced linear OT | Maximum of row residual, column residual, reduced-cost violation, and absolute primal-dual gap; `Inf` unless certified optimal |
| Entropic fixed-mass partial linear OT | Maximum checked Dykstra update, row/column capacity violation, and transported-mass residual; KKT and duality fields are separate mandatory certificates |
| Translation-invariant KL-UOT | Maximum native iterate change, independently recomputed TI fixed-point residual, and log-marginal KKT residual; primal-dual and gauge checks are separate mandatory certificates |
| Entropic GW/FGW | Frobenius norm of the outer plan update |
| Unregularized CG GW/FGW | Absolute objective change |
| FUGW / UCOOT | L1 change of the sample coupling |
| Partial CG | Relative objective change |
| Penalized variable-mass partial FGW | Certified Frank-Wolfe stationarity gap |
| Semirelaxed CG | Absolute or relative objective change (`abs_error`, `rel_error`) |

### Entropy-regularized fixed-mass partial linear OT

`ot_partial_sinkhorn()` keeps the exact subcoupling constraints from
`ot_partial_emd()` and minimizes

```text
<M,G> + epsilon * sum(G * (log(G) - 1)).
```

This is counting-measure entropy, not `KL(G || p %o% q)`. The
entropy-minus-one constant is retained because the public result reports the
complete regularized objective. Zero-weight rows and columns are removed from
the active product support and returned as exact zeros.

`method = "scaling"` is certified only when the shift-invariant cost range
divided by `epsilon` is at most 100. Unsafe explicit scaling fails;
`method = "auto"` selects genuine log-domain Dykstra. Convergence requires the
checked update, both capacity inequalities, transported mass, independent
objective, capacity-dual feasibility, complementarity, stationarity, and
primal-dual gap to pass. `warm_state` contains the plan representation and all
three Dykstra corrections; it is rejected unless certified and tied to the
identical cost, measures, mass, epsilon, active support, and effective backend.
Entropic partial GW and FGW reuse this solver for every required projection
and propagate its nested status.

### Penalized variable-mass partial linear OT

`ot_partial_penalized()` preserves the supplied finite source and target
measures while optimizing how much mass to transport. Its objective is

```text
sum(M * plan) + lambda * (sum(p) - sum(plan))
              + lambda * (sum(q) - sum(plan)).
```

Thus one transported unit avoids two discard penalties. Increasing `lambda`
weakly favors more transported mass; at an exact tie, multiple optimal masses
can exist. This linear capacity-constrained formulation is distinct from the
generalized-KL marginal relaxation controlled by `rho` in
`ot_sinkhorn_unbalanced()`.

The solver uses an exact symmetric dummy-node reduction: real-to-dummy and
dummy-to-real arcs cost `lambda`, while dummy-to-dummy costs zero. No big-M
constant is introduced. A common augmented-mass normalization is used only to
call the simplex exact-OT backend; the returned plan, primal and dual
objectives, gap, and discarded masses are scaled back to the original
finite-measure units.

### Penalized variable-mass partial FGW

`penalized_partial_fused_gromov_wasserstein()` optimizes over every
nonnegative plan satisfying `rowSums(plan) <= p` and
`colSums(plan) <= q`; no total mass is fixed. Its unrooted objective is

```text
(1 - alpha) * <M,G> + alpha * GW_square(C1,C2,G)
+ lambda * (sum(p) + sum(q) - 2 * sum(G)).
```

The mass penalty is included in the reported decomposition, full gradient,
exact linearized direction problem, and one-dimensional line search. Each
direction reuses `ot_partial_penalized()`. When the structural gradient has
negative entries, a common cost shift and half-shift of `lambda` produce an
exactly equivalent nonnegative-cost problem; both values are recorded for
every direction.

Convergence means subcoupling feasibility, direct objective and line-search
agreement, certified nested linear solves, and a Frank-Wolfe gap below its
reported tolerance. Because square-loss FGW is nonconvex, this is a stationary
point certificate rather than a global-optimality claim. Increasing `lambda`
has the derived more-transport direction for every exact linearized subproblem,
but separate nonconvex runs need not be globally monotone across local
stationary points. Fixed-mass partial FGW and generalized-KL FUGW remain
distinct estimands.

## Parameter policy for flagship solvers

### Operational on `ot_sinkhorn`

`M`, `p`, `q`, `epsilon`, `method`, `max_iter`, `tol`, `init_plan`, and
`init_duals`. Entropic plan starts must be finite, nonnegative, strictly
positive on the active product support, and zero outside zero-weight support.
Dual starts may be a prior `rfugw_result` or a list with finite `source` and
`target` vectors of the requested lengths. When both are supplied, both are
validated and duals take precedence. Scaling and log backends return the same
dual convention up to numerical tolerance: the source potential has weighted
mean zero, the target receives the opposite gauge shift, and inactive support
coordinates use zero by convention. Plan reconstruction applies on the active
product support; the returned plan is exactly zero elsewhere. Returned state never carries converged
status into a new call.

The balanced regularized objective is
`<M,T> + epsilon * KL(T || p %o% q)`. Its dual is evaluated under that same
reference-measure convention and must pass the reported gap tolerance before a
converged result remains certified. `ot_sinkhorn_divergence()` combines one
cross and two self objectives with coefficients `1, -1/2, -1/2`. It exposes
all three results and fails closed when any component is uncertified or the
combined value is negative beyond its numerical tolerance.

### Operational on `fgw_entropic` / `entropic_*gromov_wasserstein`

`M`, `C1`, `C2`, `p`, `q`, `alpha`, `epsilon`, `max_iter`, `tol`,
`sinkhorn_max_iter`, `sinkhorn_tol`, `init_plan` / `G0`, `structure_rank`,
`sinkhorn_method` (`scaling`, `log`, `auto`), `precision` (`mixed`, `double`,
`strict_double`),
`symmetric`, `solver` (`PGD`, `PPA`), `check_every`.

The float tolerance boundary is `1e-6`. A `mixed` request with either outer
or inner tolerance below that boundary is promoted to `strict_double`; the
requested tolerance is preserved. `strict_double` never enters a float path.
For scaling-domain PGD problems at least 32 by 32, `double` may select the
reported `mixed_accelerated` backend only when both tolerances are at least
`1e-6`. Log Sinkhorn, PPA, tight requests, and smaller problems remain double.
The result records the transition and the C++-reported arithmetic path.
`sinkhorn_method = "auto"` computes a conservative scaled-cost criterion. In
double precision, scaling is eligible only when both the maximum exponent
magnitude and scaled span are at most 500 (50 for mixed/float arithmetic).
Auto records the requested/effective method, threshold, metric, transition, and
reason. It selects genuine log-domain Sinkhorn outside that regime. Explicit
scaling outside the same regime errors and never silently returns a clipped
kernel as the requested problem.

`loss_fun` is accepted only as `"square_loss"`; any other value is
unsupported. In particular, `"kl_loss"` is deliberately rejected because the
directed loss is infinite on ordinary zero-diagonal distance costs unless an
estimand-changing logarithm floor is introduced. The prototype, POT convention,
near-zero differential, risks, and reconsideration gate are recorded in
`inst/kl-structural-loss-decision.md`.

### Operational on `fgw_exact_cg` / `gromov_wasserstein` / `fused_gromov_wasserstein`

`M`, `C1`, `C2`, `p`, `q`, `alpha`, `symmetric`, `G0`, `max_iter`,
`tol_rel`, `tol_abs`, `lp_solver`, `lp_max_iter`, `lp_tol`.
`lp_scale` is operational only for the lpSolve backends.

`G0` is a feasible warm start on **every** backend, including
`cpp_transport` (`bd-01M05QY3DS3E3YHQQ0DHWMADDX`). Invalid shape,
negativity, non-finite values, or marginal violations fail before
computation.

### Operational on `fugw_kl`

`Cx`, `Cy`, `wx`, `wy`, `reg_marginals`, `epsilon`, `alpha`, `M`,
`init_pi`, `max_iter`, `tol`, `max_iter_ot`, `tol_ot`, `rescale_plan`,
`check_every`, `precision` (`mixed`, `double`, `strict_double`). Tight mixed
requests are promoted to strict double with no tolerance floor.

### UCOOT / across-spaces (`bd-01M05QY37G4BNPY5WNH9K387A4`)

Supported:

- `divergence = "kl"`
- `unbalanced_solver = "sinkhorn"`; the former `"sinkhorn_log"` scaling alias
  is deprecated and errors rather than claiming log-domain behavior
- `reg_type = "joint"` or `"independent"`
- `epsilon > 0` (default `1e-2`)

Unsupported and rejected before computation:

- `divergence = "l2"`
- `unbalanced_solver = "mm"` or `"lbfgsb"`

Explicitly ignored:

- `init_duals` (POT-shaped; unused)
- `...` extra unused arguments (rejected if present after this contract)

### Explicitly ignored POT-compat parameters

| Parameter | Where | Policy |
|---|---|---|
| `thres`, `warn` | Partial solvers | Ignored |
| `random_state` | Unregularized semirelaxed CG | Ignored |
| `alpha`, `gamma_init`, factorized-cost, Dykstra, seed, and warning arguments | `lowrank_gromov_wasserstein_samples` | Unsupported and rejected when supplied; migrate to `dense_gromov_wasserstein_plan_svd()` |

## Symmetry (`bd-01M05QY3MGVG3YQS0DTRG1PR7K`)

- Default is auto-detect: `symmetric = NULL` means
  `max(abs(C - t(C))) <= 1e-10 + 1e-12 * max(abs(C), 1)` on both structure
  costs. The absolute-plus-relative rule is stable across matrix scales.
- `symmetric = TRUE` is a correctness claim. If the costs fail the
  tolerance, the call errors. It is not a silent fast-path override.
- `symmetric = FALSE` always uses the two-sided tensor.
- The symmetric fast path and the general path must agree on symmetric
  inputs.

All public iteration, rank, sampling-budget, dummy-node, and check-interval
counts must be exact finite integers in range. Fractional values are rejected;
they are never truncated before validation.

## Warm starts

| Solver | Parameter | Contract |
|---|---|---|
| Entropic GW/FGW | `init_plan` / `G0` | Nonempty plans are used; they need not be exactly feasible because Sinkhorn projects. Non-finite or negative entries fail. |
| Unregularized GW/FGW | `G0` | Must be feasible (shape, finite, nonnegative, marginals). Consumed by C++ and lpSolve backends. |
| FUGW | `init_pi` | Used as the sample/feature start. |
| UCOOT / across-spaces | `init_pi` | List with `pi_samp` and `pi_feat`. Consumed as the BCD start. Inner Sinkhorn scalings are warm-started across BCD steps, with a recorded fallback if the residual guard rejects them. |
| Partial | `G0` | Must satisfy partial mass constraints. |
| Semirelaxed | `G0` | Must satisfy source (row) marginals. |

## Threading, RNG, and memory

OpenMP is optional. Single-subject flagship solvers are serial: changing
`OMP_NUM_THREADS` must not change their plans or objectives. The only
parallel kernels are the C++ batched paths used by `multialign_fit()`
(`n_threads`) and `cpp_feature_cost_batch()`.

- `n_threads` is subject-parallel. It does not thread a single Sinkhorn
  or FGW solve.
- `n_threads = 1` is the serial fallback and must match `n_threads > 1`
  within `1e-10` on the same data and stopping rules.
- When `n_threads > 1`, BLAS is pinned to one thread unless
  `RFUGW_PIN_BLAS_THREADS=0`. Running OpenMP subjects and a threaded
  BLAS together is oversubscription and is not a certified speed path.
- Speed claims for threading apply only to those batch kernels. See
  `inst/bench/benchmark_thread_scaling.R` and
  `inst/bench/threading-memory.md`.
- Hosted ASan/UBSan builds disable OpenMP. They certify serial memory
  safety, not data races. Threaded correctness includes 1/2/4 equivalence;
  the nightly pure-C++ OpenMP harness runs under ThreadSanitizer for shared
  caches, disjoint writes, exception capture, and nested suppression. See
  `inst/thread-safety.md`.

RNG and deterministic sampling:

- Flagship solvers are deterministic given the inputs. They do not draw
  random numbers.
- `sampled_gromov_wasserstein()` is stochastic. `random_state` (or
  `set.seed`) makes a run reproducible.
- `sampled_gromov_wasserstein_coords(sampling = "deterministic")` uses
  top-k probability indices instead of weighted sampling. C++ and R
  paths must match. This is for parity checks, not a quality claim.

Memory:

- Returned plans are dense.
- `structure_knn` changes the structure metric by filling far entries
  with `max(C)`. It does **not** reduce storage; the matrix stays
  `n x n`.
- `use_cpp_feature_fused = TRUE` avoids materializing `M_list` in R.
  That is the certified densification reduction. Evidence:
  `inst/bench/threading-memory.md`.

## Known discrepancies and linked tickets

| Issue | Ticket |
|---|---|
| UCOOT advertised unsupported solvers / unusable default epsilon | `bd-01M05QY37G4BNPY5WNH9K387A4` |
| Exact CG `G0` ignored on `cpp_transport` | `bd-01M05QY3DS3E3YHQQ0DHWMADDX` |
| Flagship solvers default `symmetric = TRUE` without validation | `bd-01M05QY3MGVG3YQS0DTRG1PR7K` |
| Missing `status` / `converged` / residual fields | `bd-01M05QY3VJTSR30MN5PGP7T1Y0` |
| Incomplete validation of alpha, ranks, degenerate mass | `bd-01M05QY426X8F52TTBEM34AXEY` |
| “Exact” wording that can be read as global optimality | `bd-01M05QY49XVBDCEX94N7W6T9M7` |
| Uniform result class | `bd-01M05QY6A7T8RX4D28A4D8XCZ7` |
| Linear-OT certification vs analytic/POT | `bd-01M05QY7BST72M59M3SJEV8BV9` |

## Evidence policy

Solver claims require regression evidence: invariant tests, reference
differentials where a reference exists, and answer-quality-controlled
benchmarks for performance claims. Speed evidence must follow
`inst/bench/PROTOCOL.md`. See `CONTRIBUTING.md`.
