# rfugw review verification report

Date: 2026-08-20  
Reviewed source commit: `41b41ed713e68fe36a69ad088f34385e78f5699e`  
Audit checkout: `5aef129cdb765a51f9b56b25b66fd2fd1a5f6c16`

## Executive disposition

The review's central judgment is supported. The asymmetric partial gradient,
mixed-precision tolerance contract, convergence/status contract, partial inner
solver propagation, and public KL evaluator defects were reproduced in fresh R
processes. The claimed boundary around the symmetric/general flagship code is
also supported: independent one-step probes matched standard entropic FGW,
exact FGW, and semirelaxed FGW on asymmetric inputs while partial FGW failed.

One qualification is important. The exact transport simplex does conflate all
early exits in its Boolean status, but no false optimality certificate was found
in 600 seeded random and degenerate small problems. Reachability of those
internal failure branches through valid public inputs remains inconclusive. The
partial exact solver's loss of the inner status is independently reproduced and
material.

One additional release blocker was found: the exact reviewed tarball fails
`R CMD check` because a benchmark-protocol test uses a source-tree-relative
path that is invalid in the check layout.

## Provenance and baseline

- The current checkout was clean before audit artifacts were added.
- No package source, test, documentation, or workflow file differs between the
  reviewed commit and the audit checkout; all intervening changes are `.mote/`
  tracker metadata.
- Clean source installation succeeded with R 4.5.1, Rcpp 1.1.2,
  RcppArmadillo 15.4.2.1, Homebrew clang 20.1.8, and Apple arm64.
- The source-tree test suite passed with 616 or more passing expectations and
  two expected skips.
- The repository numerical accuracy gate passed all 28 configured checks.
- Hosted GitHub Actions were not rerun or treated as verified by this local
  audit.

### Isolated R CMD check

The exact reviewed commit was exported with `git archive`, built in `/tmp`, and
checked with `LC_ALL=C LANG=C R CMD check --as-cran --no-manual`.

Final status: **1 ERROR, 0 WARNING, 2 NOTE**.

- Error: `tests/testthat/test-bench-protocol.R` sources
  `inst/bench/protocol.R` through `testthat::test_path("..", "..", ...)`.
  That path exists in a source checkout but not in the installed/check test
  layout. The check stopped after 616 passes, one test failure, one test warning,
  and two skips.
- Note: new submission and development version `0.0.1.9000`.
- Note: `CONTRIBUTING.md` and `_pkgdown.yml` are non-standard top-level files.

An earlier check attempt failed at DESCRIPTION metadata because of invalid host
locale variables. The controlled-locale rerun eliminated that host artifact.

## Verified findings

| ID | Disposition | Evidence | Priority |
|---|---|---|---|
| P0-01 asymmetric partial gradient | **Reproduced** | Objective matched brute force to `2.22e-16`, but native gradient max error was `1.17`; native directional derivative error was `1.33` while brute derivative error was `3.90e-10`. | P0 |
| P0-02 mixed tolerance floor | **Reproduced** | Requested `1e-9` and requested `1e-6` mixed runs produced bit-identical plans. | P0 |
| P0-03 false mixed termination reason | **Reproduced** | Mixed run stopped at iteration 18 of 500 with residual `7.63e-7` and returned `status = "max_iter"`; strict double converged at `6.86e-10`. | P0 |
| Hidden double-to-mixed dispatch | **Reproduced** | At 32x32, public `precision = "double"` and `"mixed"` produced bit-identical plans and iteration counts; result backend was only `"cpp"`. | P0 contract defect |
| P0-04 simplex early-exit conflation | **Confirmed in implementation; runtime reachability inconclusive** | Every early `break` shares `converged = (it < max_iter)`. No invalid certificate was found in 600 differential cases; maximum objective gap from `lpSolve` was `2.66e-15`. | P0 hardening boundary |
| P0-05 partial exact discards LP status | **Reproduced** | Returned fields contain no `lp_ok` or inner iterations. One-pivot and full inner solves differed by `0.259` in objective and `0.2265` in plan. | P0 |
| P0-06 formulation-aware convergence | **Reproduced** | Public FGW reported convergence with row infeasibility `1.35e-5` against requested inner tolerance `1e-12`; inner residual was `NA`. | P0 |
| FUGW inner certification | **Reproduced** | One-step and well-solved inner budgets both reported converged; objective gap `0.00322`, plan max error `0.0115`, both inner residuals `NA`. | P0 |
| Entropic partial inner certification | **Reproduced** | Adapter exposes neither status nor inner residual; inner budgets changed objective by `0.0289` and plan by `0.0863`. | P0 |
| Central result helper validation | **Reproduced** | A negative, marginal-infeasible plan with infinite `ot_dist`, generic `objective`, and inner residual was labelled converged. | P0 |
| P0-07 generalized KL mass terms | **Reproduced** | Unbalanced example omitted exactly `0.4 = -mass(plan) + mass(reference)` and returned a negative value where generalized KL was positive. | P0 |
| P0-08 KL support rule | **Reproduced** | Positive mass outside zero reference support returned finite `137.5`; mathematical result is `Inf`. Invalid reference weights are not rejected. | P0 |
| P0-09 asymmetric coverage gap | **Confirmed** | Partial fixtures are symmetric. Existing asymmetric coverage is concentrated in auto-detection and semirelaxed comparisons, without independent gradient laws. | P0 gate |
| Standard asymmetric FGW | **Verified on independent one-step probes** | Entropic and exact one-step plans matched independent brute-gradient constructions within `8.33e-17` and `5.55e-17`. | Keep regression |
| Asymmetric semirelaxed FGW | **Verified on independent one-step probe** | Native plan matched the brute-gradient row-wise direction exactly; objective error `5.55e-17`. | Keep regression |
| Asymmetric partial solver effect | **Reproduced** | Native one-step plan differed by `0.20` and objective was `0.00105` worse than the correct independent step. | P0 |
| Scaling clipping | **Reproduced** | Adding a constant `1e6` changed the scaling plan by `0.25`; log-domain plan changed by `2.54e-9`. | P1 |
| True log-domain UOT gap | **Confirmed** | `ot_sinkhorn_unbalanced()` has no method selector; UCOOT `sinkhorn` and `sinkhorn_log` outputs were bit-identical. | P1 |
| Fractional counts | **Reproduced** | `.validate_count(2.9)` returns integer `2`. | P2 |
| Partial mass validation | **Confirmed behavior; contract needs clarification** | A zero-mass partial plan is accepted when `mass` is omitted. No intended target mass can be inferred without that argument. | P3 |
| Reference-vector validation | **Reproduced** | Relaxed validation accepts supplied vectors containing `NA`, `Inf`, and negative values because only length is checked. | P2 |
| Absolute-only symmetry | **Reproduced** | `symmetric = TRUE` rejected costs with relative asymmetry `9.77e-16` because absolute difference exceeded `1e-10`. | P2 |
| Environment-dependent behavior | **Confirmed in source** | Warm starts, adaptive tolerances, residual floors/caps, matvec dispatch, sampled precision, and autotuning read `RFUGW_*` variables; solver results do not record them. | P1 |
| Monolithic native core | **Confirmed** | `src/fgw_core.cpp` is 8,115 lines and contains all major solver families, diagnostics, precision variants, batching, and environment controls. | P1 architecture |
| Structure-rank setup | **Confirmed** | Both precision paths call full `arma::svd_econ()` before truncation. | P2/documentation |
| Low-rank samples naming | **Confirmed and already documented** | Function solves dense entropic GW, materializes `T`, then calls R `svd(T)`. | P2/documentation |
| Graph/dense behavior | **Confirmed and already documented** | Returned plans and structure matrices remain dense; kNN fills far entries rather than changing storage. | P2/documentation |
| OpenMP Rcpp boundary | **Confirmed in source** | Worker loops call `Rcpp::NumericMatrix::begin()`, `nrow()`, `ncol()`, and vector accessors inside parallel regions. | P1 hardening |
| Benchmark threshold interpretation | **Confirmed** | FUGW and UOT quality rules admit `max_iter`; sampled objectives allow absolute and relative tolerance `1`. These are regression gates, not convergence certificates. | P1 semantics |
| `alpha` semantics | **Confirmed** | FGW uses `(1-alpha) * M + alpha * GW`; FUGW constructs both linear costs as `(alpha/2) * M`. | P2 API |
| Release check path | **New finding, reproduced** | Exact reviewed tarball fails at `test-bench-protocol.R:1` under R CMD check. | P1 release blocker |

## Why the existing gates passed

The existing tests and accuracy gate are useful and should remain. They do not
exercise the failure laws above:

- Structural fixtures are overwhelmingly symmetric.
- Objective equality cannot detect the asymmetric target-marginal orientation
  error because `q' A q == q' A' q`.
- Mixed/double proximity thresholds do not test requested/effective tolerance
  truthfulness.
- Outer stopping tests do not require inner residual or final feasibility.
- Exact-OT tests compare plans/objectives but do not demand dual/reduced-cost
  certificates or branch-specific termination reasons.
- Balanced KL tests cancel the missing generalized mass terms.

## Dependency-aware remediation plan

### Phase 0: freeze and establish gates

1. Temporarily mark asymmetric partial GW/FGW experimental or reject
   `symmetric = FALSE` in exact and entropic partial entry points.
2. Fix the R CMD check test path using an installed-resource lookup with a
   source-checkout fallback.
3. Add the seven audit reproducers here as failing regression tests before
   changing solver code.
4. Define release evidence channels explicitly: local check, hosted CI,
   benchmark artifact, and publication record.

Acceptance: exact reviewed failures reproduce; corrected branch starts with a
clean local package check and the new tests fail for the intended reasons.

### Phase 1: canonical square-loss algebra

1. Introduce a small native module with explicit functions for forward tensor,
   reverse tensor, scalar objective, and full gradient.
2. Correct the forward target term to `(C2 % C2) * qG`; retain the transpose for
   the reverse term.
3. Route partial exact and entropic solvers through the canonical implementation.
4. Remove independent duplicated formulas only after equivalence tests cover
   float/double, symmetric/general, balanced/partial, and semirelaxed paths.

Acceptance:

- Seeded finite-difference and directional-derivative laws pass on random
  asymmetric inputs.
- Exact and entropic partial paths agree with independent small-problem
  references.
- Symmetric and general paths agree on symmetric inputs.
- Existing POT fixtures and accuracy gates remain within their declared bounds.

### Phase 2: exact transport termination and certificates

1. Replace the Boolean with an explicit termination enum: optimal, max-iter,
   disconnected basis, invalid cycle, invalid step, no leaving variable, and
   numerical failure.
2. Return final row/column residuals, potentials, minimum nonbasic reduced cost,
   primal objective, dual objective, and duality gap.
3. Treat only an optimal, feasible, gap-certified result as successful.
4. Make partial transport return the complete inner report; propagate failure
   through partial FGW and exact FGW.
5. Add deterministic branch tests or an instrumented native test seam for
   branches that valid public inputs rarely reach.

Acceptance: `converged` implies primal feasibility, reduced-cost optimality, and
duality-gap bounds; low inner iteration budgets can never be reported as solved.

### Phase 3: truthful precision and tolerance semantics

1. Decide the public contract for tight mixed requests: reject, promote to
   strict double, or expose the relaxed floor. Do not silently change it.
2. Add `strict_double` if `double` is intended to retain float acceleration.
3. Return requested/effective outer tolerance, inner tolerance, precision,
   backend, approximation, and termination reason from native kernels.
4. Evaluate status against the effective rule that actually stopped the solver.
5. Apply the same schema to batched/multialign paths.

Acceptance: an early tolerance stop is never called max-iteration; results make
the 32x32 dispatch and all floors observable; contract tests cover thresholds.

### Phase 4: formulation-aware result certification

1. Define a common report structure but attach formulation-specific certificate
   requirements.
2. Always validate finite and nonnegative plans and finite public objectives.
3. Balanced solvers require row and column feasibility; partial solvers require
   upper bounds and transported mass; nested solvers require all final inner
   solves to pass; unbalanced solvers require their fixed-point certificate.
4. Preserve distinct reasons for max iterations, stagnation, inner failure,
   infeasibility, and numerical failure.
5. Recompute the documented objective independently at the returned plan.

Acceptance: a package-wide property test establishes
`converged => all formulation-relevant certificates pass`.

### Phase 5: repair the public KL evaluator

1. Make `ot_kl()` compute generalized KL by default, including mass terms.
2. Return `Inf` when positive plan mass lies outside reference support.
3. Validate reference vectors as finite and nonnegative.
4. If compatibility requires the old expression, expose it under an explicit
   relative-entropy-term name rather than overloading KL.

Acceptance: balanced, unbalanced, zero-support, zero-mass, and invalid-reference
tests pass against an independent implementation.

### Phase 6: adversarial test and release matrix

Add mandatory test families for every supported solver:

1. Gradient laws and asymmetric differential fixtures.
2. Convergence implication and objective-recomputation laws.
3. Exact transport primal-dual certificates and degenerate bases.
4. Zero supports, zero weights, singletons, and partial mass boundaries.
5. Additive shifts, scale extremes, negative costs, and small epsilon.
6. Mixed requested/effective contract and environment provenance.
7. Random seeded small-problem differentials against brute force, POT, and
   `lpSolve`.
8. Threaded native tests whose worker regions touch only owned C++ data.

Acceptance: local R CMD check, Linux/macOS/Windows checks, sanitizer jobs,
linear-OT certification, flagship accuracy, and answer-quality-controlled
benchmarks all pass on the same release commit. Hosted artifacts must be linked;
local success is not hosted proof.

### Phase 7: numerical and API hardening

1. Add `method = "auto"` for balanced Sinkhorn using dynamic range and requested
   tolerance; implement true log-domain unbalanced Sinkhorn before advertising
   `sinkhorn_log` elsewhere.
2. Reject fractional counts and validate all supplied reference vectors.
3. Use combined absolute/relative symmetry checking.
4. Clarify partial validation: require `mass` when mass certification is wanted,
   or document that it checks only upper bounds when omitted.
5. Add semantic aliases such as `structure_weight` and `feature_weight`; record
   the interpreted objective decomposition.
6. Record all environment-derived solver settings or replace them with explicit
   controls.

### Phase 8: architecture and threading

After behavioral contracts are locked:

1. Split linear OT, GW algebra, solver families, approximations, batching, and
   diagnostics into separately testable native modules.
2. Convert Rcpp inputs to owned or pointer/shape plain-C++ views before entering
   OpenMP regions; do not call Rcpp methods from workers.
3. Keep POT adapters thin over canonical native result-producing functions.
4. Preserve the existing benchmark protocol while labeling full-SVD and dense
   setup costs explicitly.

### Phase 9: breadth only after hardening

Do not gate 0.1 on new breadth. After Phases 0-8, consider public partial linear
OT, regularized objective decompositions, dual outputs, Sinkhorn divergence,
and fixed-support Wasserstein barycenters. GPU, autodiff, arbitrary losses, and
full POT parity remain reasonable non-goals.

## Release decision

A trust-oriented 0.1 should remain blocked on Phases 0-6. A narrower alpha can
remain useful if asymmetric partial methods are disabled or prominently marked,
mixed precision is described as approximate, and `converged` is not presented
as a complete numerical certificate until the report contract is repaired.
