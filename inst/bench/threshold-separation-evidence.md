# Threshold evidence: certification versus fixed-budget performance

The previous threshold schema allowed `max_iter` for unbalanced Sinkhorn and
FUGW quality rows, while sampled GW had no required status. That made a loose
fixed-budget run look interchangeable with a converged solver result. The
schema-v2 boundary removes that ambiguity:

- certified comparisons require `status == "converged"` plus formulation
  feasibility, independent objective consistency, named component consistency,
  and any required inner certificate;
- fixed-budget FUGW accepts only its declared `max_iter` sentinel and is labeled
  `fixed_budget_performance`, `certified = FALSE`, and
  `comparison_eligible = FALSE`;
- sampled GW remains `experimental`; its three budget rows report exact-plan
  objective error versus budget and never claim convergence;
- log-domain unbalanced Sinkhorn has a converged certification configuration;
  any future `max_iter` row must use its separate performance-only contract.

On the local deterministic n=6/n=8 fixtures, the formerly loose FUGW settings
produced `inner_failure` after outer convergence pressure, whereas a 10-outer-
iteration budget produced the intended `max_iter` sentinel (residual about
`0.005` on n=6, below its performance-only ceiling of `0.01`). Log-domain
unbalanced Sinkhorn at tolerance `1e-6` produced `converged`; sampled GW reports
`experimental`. These observations justify separating the evidence classes;
they do not establish portable runtime thresholds.

Historical CSVs produced before schema v2 retain their meaning but cannot be
promoted into the certified series. A threshold edit changes the MD5 recorded
in `threshold-history.json`; the protocol refuses to run until a new retained
evidence entry and review requirement accompany that edit.

## Nightly boundary correction (2026-08-21)

The deterministic nightly fixture at seed 20260816 exposed three boundary
settings that did not satisfy their declared evidence class. The FUGW n=8 run
converged before its ten-iteration fixed budget and therefore could not be
classified as a `max_iter` performance sentinel; a nine-iteration budget keeps
both n=8 and n=12 rows fixed-budget, with outer residuals approximately
`6.9e-6` and `7.0e-3`. Entropic semirelaxed GW at n=12 required 675 iterations
to reach its unchanged `1e-6` tolerance, so its outer budget is now 800.

UCOOT's unbalanced Sinkhorn loop previously accepted heuristic stopping tests
before checking the fixed-point residual used by its public certificate.
After requiring candidate stops to satisfy that reported residual, the
unchanged regularization and `1e-7` inner tolerance certify both UCOOT methods
at n=8 and n=12 with a 5000-iteration safety budget. Observed worst inner
residuals were below `1e-7`; all four outer solves converged below `1e-6`.
These are deterministic correctness budgets, not portable runtime claims.
