# Post-certification linear-OT foundation decisions

This assessment is deliberately narrower than POT. It covers four candidate
additions after the focused 0.1 certification gates; it does not imply broad
API or algorithm parity.

| Candidate | Decision | Evidence and boundary |
|---|---|---|
| Exact partial linear OT | Accepted as `ot_partial_emd()` | Reuses the certified exact simplex through a one-dummy reduction. The public contract requires nonnegative costs, normalized weights, mass in `[0,1]`, partial marginal/mass feasibility, and the augmented primal/dual certificate. Analytic boundaries, an independent `lpSolve` oracle, permutation and scaling laws, invalid-input adversaries, regression tests, and a local baseline cover it. |
| Genuine log-domain KL-unbalanced Sinkhorn | Accepted as the existing `ot_sinkhorn_unbalanced(method="log"/"auto")` | The log backend is not a scaling alias. Moderate regimes match the independent scaling implementation, permutation equivariance is tested, unsafe dynamic range forces log and rejects explicit scaling, certificate telemetry is retained, and the log/auto benchmark records its cost. |
| Entropic partial linear OT | Deferred | The existing internal Dykstra projection is scaling-domain only. There is no genuine log-domain Dykstra implementation or complete convergence certificate. Public exposure would contradict the new robust log/auto contract; it remains an internal primitive for the explicitly bounded partial-GW path. |
| Sinkhorn divergence | Deferred | `ot_sinkhorn()` intentionally reports the unregularized linear term `<M,G>`, not a certified regularized primal/dual value. Debiasing that field would define the wrong functional. Admission requires a stable regularized-OT convention, self-cost oracle, positivity/identity laws, and dual/primal certificates. |
| Fixed-support Wasserstein barycenter weights | Deferred | Existing certified primitives return couplings but not the dual-gradient contract needed for a trustworthy fixed-support weight optimizer. Admission requires simplex-constrained optimality/KKT diagnostics, permutation and identical-input laws, adversarial zero weights, and quality-versus-iteration evidence. Existing FGW barycenters solve a different structural problem and are not reused under a misleading name. |

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

The representative local benchmark is in
`inst/bench/linear-foundations-baseline.csv`. It records certificates before
timing, setup/solve/end-to-end time, allocation size, status, objective, mass
residual, precision/thread provenance, commit, and environment. It is local
evidence only, not hosted or publication evidence.
