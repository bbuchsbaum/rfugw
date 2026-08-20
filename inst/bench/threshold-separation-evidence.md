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
