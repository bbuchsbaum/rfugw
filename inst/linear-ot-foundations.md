# Post-certification linear-OT foundation decisions

This assessment is deliberately narrower than POT. It covers five candidate
additions after the focused 0.1 certification gates; it does not imply broad
API or algorithm parity.

| Candidate | Decision | Evidence and boundary |
|---|---|---|
| Exact partial linear OT | Accepted as `ot_partial_emd()` | Reuses the certified exact simplex through a one-dummy reduction. The public contract requires nonnegative costs, normalized weights, mass in `[0,1]`, partial marginal/mass feasibility, and the augmented primal/dual certificate. Analytic boundaries, an independent `lpSolve` oracle, permutation and scaling laws, invalid-input adversaries, regression tests, and a local baseline cover it. |
| Genuine log-domain KL-unbalanced Sinkhorn | Accepted as the existing `ot_sinkhorn_unbalanced(method="log"/"auto")` | The log backend is not a scaling alias. Moderate regimes match the independent scaling implementation, permutation equivariance is tested, unsafe dynamic range forces log and rejects explicit scaling, certificate telemetry is retained, and the log/auto benchmark records its cost. |
| Entropic partial linear OT | Accepted as `ot_partial_sinkhorn()` | The fixed-mass counting-measure entropy objective is explicit. Moderate scaling and genuine log Dykstra agree; unsafe scaling fails and auto selects log. Independent dual optimization, boundary/zero-support/constant/tiny-entropy fixtures, KKT and primal-dual certificates, backend-bound warm states, and partial GW/FGW parity cover the public path. |
| Sinkhorn divergence | Accepted as `ot_sinkhorn_divergence()` | `ot_sinkhorn()` retains the unregularized linear term `<M,G>` and separately certifies the regularized primal `<M,G> + epsilon KL(G || p %o% q)` against its matching dual. The divergence combines one cross and two self values under that convention, retains all component results, supports reusable dual state, and fails closed on component or nonnegativity failure. Independent 2x2 optimization, constant-offset, identity, symmetry, permutation, joint cost-epsilon scaling, zero-support, tiny-epsilon, and failure-propagation tests cover it. |
| Fixed-support Wasserstein barycenter weights | Accepted as `ot_barycenter_weights()` | Exact mode solves the joint fixed-support LP and independently re-certifies every component with `ot_emd()`. Regularized mode optimizes a semi-debiased product-reference-KL Sinkhorn objective using certified dual gradients, safeguarded descent, a simplex projected-gradient KKT map, reusable component state, and an explicit numerical weight floor. Independent joint-LP/convex-optimizer, finite-difference, identity, permutation, zero-weight, duplicate-support, coefficient, tiny-entropy, warm-state, and benchmark evidence guard it. Existing GW/FGW barycenters remain distinct support/structure estimands. |

## Accepted exact-partial contract

For nonnegative `M`, probability vectors `p` and `q`, and transported mass
`m`, `ot_partial_emd()` solves

`min_G <M,G>` subject to `G >= 0`, `G 1 <= p`, `G^T 1 <= q`, and
`sum(G) = m`.

It augments each side by one dummy coordinate with mass `1-m`; a penalty larger
than every real cost makes augmented dummy-to-dummy mass zero at optimum, so
the real block transports exactly `m`. The returned `source_potential` and
`target_potential` include the dummy coordinate because they certify the
equivalent augmented LP. `mass_certified`, `feasible`,
`objective_consistent`, reduced-cost feasibility, and duality gap must all pass
before `status="converged"`.

## Accepted entropy-regularized partial contract

For finite `M`, probability vectors `p` and `q`, requested mass `m`, and
positive `epsilon`, `ot_partial_sinkhorn()` solves

`<M,G> + epsilon * sum(G * (log(G) - 1))`

on the same fixed-mass subcoupling feasible set. Zero entries contribute zero.
This is the counting-measure entropy-minus-one convention, not the
product-reference KL convention used by balanced `ot_sinkhorn()` diagnostics.
The `-epsilon * m` term is constant for fixed mass but is reported rather than
silently dropped.

The scaling Dykstra backend is limited to shift-invariant dynamic range 100.
An explicit unsafe request errors; `method = "auto"` selects the genuine
log-domain implementation. A result is converged only when its update,
capacity and mass residuals pass and independently reconstructed capacity
duals, mass dual, stationarity, complementarity, dual feasibility, objective,
and primal-dual gap agree. Warm state includes all three Dykstra corrections
and is accepted only for an identical problem and effective backend.
Entropic partial GW and FGW now call this public solver for every linearized
projection and propagate its certification.

The representative local benchmark is in
`inst/bench/linear-foundations-baseline.csv`. It records certificates before
timing, setup/solve/end-to-end time, allocation size, status, objective, mass
residual, precision/thread provenance, commit, and environment. It is local
evidence only, not hosted or publication evidence.

The dedicated scaling/log/warm/dispatch evidence is in
`inst/bench/partial-sinkhorn-baseline.csv` with interpretation in
`inst/bench/partial-sinkhorn-benchmark.md`.
