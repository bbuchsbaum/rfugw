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
6. **Weights are probability measures unless a solver is unbalanced.**
   Finite nonnegative weights with positive total mass are renormalized to
   sum 1. That renormalization is part of the contract
   (`bd-01M05QY426X8F52TTBEM34AXEY`).
7. **Breaking changes require a NEWS entry** and, after 0.1, a minor version
   bump. Experimental APIs may change in a patch if labeled experimental.

## Export classification

### Primary solvers

| Function | Family | Maturity |
|---|---|---|
| `ot_sinkhorn` | Balanced entropic linear OT | Flagship |
| `ot_emd` | Exact balanced linear OT | Flagship |
| `ot_partial_emd` | Exact nonnegative-cost partial linear OT | Supported |
| `ot_sinkhorn_unbalanced` | KL-unbalanced entropic linear OT | Flagship |
| `fgw_entropic` | Entropic FGW | Flagship |
| `fgw_exact_cg` | Unregularized FGW (CG + LP) | Flagship |
| `fugw_kl` | Fused unbalanced GW (KL) | Flagship |
| `partial_gromov_wasserstein` | Partial GW (CG + partial LP) | Supported |
| `partial_fused_gromov_wasserstein` | Partial FGW | Supported |
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
quality-versus-budget envelope in `inst/bench/sampled-budget-curves.md`.

| Function | Status | Certified envelope |
|---|---|---|
| `sampled_gromov_wasserstein` | Experimental | A full budget `(ns, nt)` is closer to dense entropic GW than a tiny budget such as `(2, 1)`, in square-loss GW and plan Frobenius distance. Intermediate budgets are not certified as monotone. Budgets `< 1` error; source/target counts above `ns`/`nt` warn and clamp. |
| `sampled_gromov_wasserstein_coords` | Experimental | Same tiny-versus-full quality envelope. Inputs scale as `O(n d)` rather than two dense `n x n` structure costs. |
| `sampled_gw_from_graphs` | Experimental | Same envelope after diffusion coordinates. A sparse graph plus `k` embeddings stores less than two dense structure costs. |
| `lowrank_gromov_wasserstein_samples` | Experimental / misnamed | Post-hoc SVD of a dense GW plan. Relative Frobenius reconstruction error decreases as rank grows up to `min(ns, nt)`. Rank `< 1` errors; rank above `min(ns, nt)` warns and clamps. Lifecycle: `bd-01M05QY9GKDDTB3CXTKAJE0E8G`. |

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

## Common diagnostic definitions

These names have one meaning across flagship solvers. Legacy field names
remain for compatibility.

| Term | Meaning |
|---|---|
| `value` / documented `*_dist` / `*_cost` | Objective of the advertised formulation, evaluated at the returned plan. Unless a field name includes entropy or KL, the value does **not** add the entropic regularizer. |
| `plan` | Coupling matrix. Balanced solvers satisfy the documented marginals up to residual. |
| `iterations` | Outer iterations executed. |
| `inner_iterations` | Inner linear-OT / Sinkhorn iterations, when a nested solver exists. |
| `inner_residual`, `max_inner_residual` | Residual from the final required inner solve and the maximum across required inner solves. A final small value does not erase an earlier uncertified solve. |
| `inner_converged`, `inner_status` | Whether every formulation-required inner solve was certified, plus its final/specific status. Outer convergence cannot be `TRUE` when this field is `FALSE`. |
| `error` | Stopping residual used by that solver (see table below). Prefer `residual` for new code. |
| `residual` | Same numeric quantity as `error`, with a documented unit. |
| `rel_error` | Relative cost change, used by unregularized CG. |
| `row_residual`, `col_residual` | Balanced: `max(abs(rowSums(plan) - p))` and `max(abs(colSums(plan) - q))`. Partial: maximum positive row/column capacity violation. |
| `mass`, `mass_residual` | Returned mass and, for partial solvers, `abs(sum(plan) - requested_mass)`. Partial log results also expose `transported_mass_target` and whether it was defaulted. |
| `feasibility`, `feasibility_residual`, `feasibility_tolerance`, `feasible` | Formulation-specific certificate. Balanced solvers check both marginals; partial solvers check marginal inequalities and transported mass; semirelaxed solvers check the fixed source marginal; unbalanced solvers use the required scaling fixed-point certificate. |
| `objective_recomputed`, `objective_residual`, `objective_tolerance`, `objective_consistent` | Independent evaluation of the advertised objective, its discrepancy from the reported value, the comparison threshold, and the resulting certificate. Named objective components expose analogous `*_recomputed`, `*_residual`, and `*_consistent` fields. |
| `objective_components_consistent` | Whether every independently checked named objective component is finite and self-consistent. |
| `status` | Normally one of `converged`, `max_iter`, `inner_failure`, `infeasible`, `objective_mismatch`, `numerical_failure`, `lp_failure`. Exact linear OT preserves its more specific internal failure reason. |
| `converged` | `TRUE` only when `status == "converged"`: the plan is finite and nonnegative, formulation feasibility passes, reported and independently recomputed objectives agree, every required inner solve is certified, and the solver's stopping rule passes. Iteration count alone never establishes success. |
| `termination_reason` | Exact linear OT uses `optimal`, `max_iter`, `disconnected_basis`, `invalid_cycle`, `invalid_step`, `no_leaving_variable`, or `numerical_failure`. |
| `source_potential`, `target_potential` | Dual potentials for exact balanced linear OT. Their weighted sum is `dual_objective`. |
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
| Entropic GW/FGW | Frobenius norm of the outer plan update |
| Unregularized CG GW/FGW | Absolute objective change |
| FUGW / UCOOT | L1 change of the sample coupling |
| Partial CG | Relative objective change |
| Semirelaxed CG | Absolute or relative objective change (`abs_error`, `rel_error`) |

## Parameter policy for flagship solvers

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
unsupported.

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
| `gamma_init`, factorized-cost and Dykstra arguments | `lowrank_gromov_wasserstein_samples` | Ignored; see experimental notice |

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
| Pseudo-low-rank GW name | `bd-01M05QY9GKDDTB3CXTKAJE0E8G` |
| KL structural loss not supported | `bd-01M05QY99ZQHFQ710DVAR7APX8` |

## Evidence policy

Solver claims require regression evidence: invariant tests, reference
differentials where a reference exists, and answer-quality-controlled
benchmarks for performance claims. Speed evidence must follow
`inst/bench/PROTOCOL.md`. See `CONTRIBUTING.md`.
