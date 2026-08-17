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
| `error` | Stopping residual used by that solver (see table below). Prefer `residual` for new code. |
| `residual` | Same numeric quantity as `error`, with a documented unit. |
| `rel_error` | Relative cost change, used by unregularized CG. |
| `row_residual`, `col_residual` | `max(abs(rowSums(plan) - p))` and `max(abs(colSums(plan) - q))` after any documented renormalization. Unbalanced solvers report mass and KL residuals instead. |
| `status` | One of `converged`, `max_iter`, `numerical_failure`, `lp_failure`. |
| `converged` | `TRUE` only when `status == "converged"`. Iteration limit and NaN/Inf breakdown are never success. |

Flagship solvers return an `rfugw_result` list. `print()` / `summary()` show
diagnostics without dumping the plan. Accessors `rfugw_plan()`,
`rfugw_value()`, `rfugw_status()`, and `rfugw_residuals()` are the stable
downstream API. Legacy fields (`plan`, `fgw_dist`, `error`, ...) remain.
The object is a named list, so `saveRDS()` / `readRDS()` and copies preserve
fields and the S3 class.

### Residual units by family

| Family | `error` / `residual` |
|---|---|
| Entropic GW/FGW | Frobenius norm of the outer plan update |
| Unregularized CG GW/FGW | Absolute objective change |
| FUGW / UCOOT | L1 change of the sample coupling |
| Partial CG | Relative objective change |
| Semirelaxed CG | Absolute or relative objective change (`abs_error`, `rel_error`) |

## Parameter policy for flagship solvers

### Operational on `fgw_entropic` / `entropic_*gromov_wasserstein`

`M`, `C1`, `C2`, `p`, `q`, `alpha`, `epsilon`, `max_iter`, `tol`,
`sinkhorn_max_iter`, `sinkhorn_tol`, `init_plan` / `G0`, `structure_rank`,
`sinkhorn_method` (`scaling`, `log`), `precision` (`mixed`, `double`),
`symmetric`, `solver` (`PGD`, `PPA`), `check_every`.

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
`check_every`, `precision`.

### UCOOT / across-spaces (`bd-01M05QY37G4BNPY5WNH9K387A4`)

Supported:

- `divergence = "kl"`
- `unbalanced_solver = "sinkhorn"` ( `"sinkhorn_log"` is an accepted alias
  of the same scaling-domain implementation until a true log-domain path
  exists)
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
  `max(abs(C - t(C))) <= 1e-10` on both structure costs.
- `symmetric = TRUE` is a correctness claim. If the costs fail the
  tolerance, the call errors. It is not a silent fast-path override.
- `symmetric = FALSE` always uses the two-sided tensor.
- The symmetric fast path and the general path must agree on symmetric
  inputs.

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
  safety, not data races. Threaded correctness is the 1-vs-N
  equivalence test. There is no TSan job.

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
