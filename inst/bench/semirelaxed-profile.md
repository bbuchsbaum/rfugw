# Semi-relaxed hot-path profile

Evidence for `bd-01M05QY815DMGYFZZFHXYQKYCY`. Times are median
milliseconds from `Rscript inst/bench/profile_semirelaxed.R` on the
conservative `-O2` profile, one thread, seed `20260816`, 40 CG
iterations. Problems come from `bench_make_problem("fgw", n, seed)`.

## Dominant work

On both exact and fast paths the per-iteration cost is:

1. GEMM for `C1 %*% G %*% t(2 C2)` (and the transpose pair when
   `C1`/`C2` are asymmetric)
2. Row-min Frank–Wolfe direction
3. Linesearch dots against `Acur` and the `C2^2 q` marginal

The R reference core also rebuilds `C1^2`, `C2^2`, and the full tensor
from scratch every iteration. That is the main exact-path waste.

## Timings

| n | path | median ms |
|---|---|---|
| 96 | R core, symmetric FGW | 4.5 |
| 96 | C++ exact, symmetric FGW | 1.0 |
| 96 | C++ fast double | 1.0 |
| 96 | C++ fast mixed | 0.5 |
| 96 | R core, asymmetric GW | 11.5 |
| 96 | C++ exact, asymmetric GW | 1.5 |
| 96 | C++ fast, empty `M` | 1.0 |
| 96 | C++ fast, explicit zero `M` | 1.5 |
| 64 | R core, asymmetric GW | 6.0 |
| 64 | C++ exact, asymmetric GW | 1.0 |

n = 32 is at timer resolution (0–1.5 ms) and is not used as evidence.

## What improved

- Asymmetric unregularized SR now uses `cpp_semirelaxed_fgw_exact_square`
  instead of the R CG core. At n = 96 that is about 8×.
- Exact CG caches `C1^2` / `C2^2` and updates `Acur` incrementally.
  Symmetric exact now matches the fast double path at these sizes.
- GW wrappers pass an empty feature cost when `alpha = 1` instead of
  allocating `ns × nt` zeros.

## Floor

Skipping the zero feature-cost matrix is real but GEMM-dominated.
At n ≤ 64 the empty-`M` vs zero-`M` difference is inside the 0.5 ms
clock. At n = 96 the fast GW path is 1.0 ms vs 1.5 ms. Mixed precision
is the larger symmetric win once `ns * nt ≥ 4000`.

Protocol rows for `semirelaxed_gromov_wasserstein`,
`semirelaxed_fused_gromov_wasserstein`, and
`entropic_semirelaxed_gromov_wasserstein` live in
`inst/bench/thresholds.json` and `inst/bench/run_protocol.R`.
Those n = 12 / 16 smokes are quality guards, not speed claims.
