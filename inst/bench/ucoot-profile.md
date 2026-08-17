# UCOOT / across-spaces hot-path profile

Evidence for `bd-01M05QY8EM21TNSPDCMJ55YS42`. Times are from
`Rscript inst/bench/profile_ucoot.R` on the conservative `-O2`
profile, one thread, seed `20260816`, 2 reps, 30 BCD iterations.
Problems come from `bench_make_problem("ucoot", n, seed)`:
sample size `n`, feature size `round(0.6 n)`.

Sample and feature phases are timed separately inside
`cpp_ucoot_kl` (`feat_ms`, `samp_ms`).

## C++ BCD vs R core

| n | path | median ms |
|---|---|---|
| 16 | C++ independent | 0.5 |
| 16 | R independent | 5.5 |
| 24 | C++ independent | 2.5 |
| 24 | R independent | 11.5 |
| 32 | C++ independent | 3.5 |
| 32 | R independent | 13.5 |
| 16 | C++ joint | 0.5 |
| 24 | C++ joint | 2.5 |
| 32 | C++ joint | 2.0 |

The C++ path is faster at every profiled size (about 4–11×).
There is no R-core crossover in this range. The R core remains
as a reference fallback.

## Sample vs feature phases

| n | feature ms | sample ms |
|---|---|---|
| 16 | 0.27 | 0.35 |
| 24 | 0.81 | 1.12 |
| 32 | 1.02 | 2.38 |

The sample phase is the larger of the two once `n_sample > n_feature`,
which is the representative rectangular layout. Both phases reuse
`SinkhornUnbalancedWorkspace` buffers (`K`, scalings, matvecs) and
in-place UOT cost workspaces (`uot_cost`, `scratch`, marginals).

## Warm starts

Inner Sinkhorn scalings are warm-started across BCD steps. The first
step cannot reuse scalings. Later steps either keep the previous
`(u, v)` or record `inner_warm_fallback_*` when the residual guard
rejects them. Tests require that after iteration 1 every feature
step is either warm-started or a recorded fallback. `init_pi` is
consumed as the BCD start.

## POT differentials

`inst/extdata/fixtures/ucoot_pot_modes.json` (POT 0.9.6.post1,
seed 20260816) covers both `reg_type = "independent"` (UCOOT) and
`"joint"` (across-spaces) on the same 7×5 / 6×4 problem. Public
solvers stay within 2e-2 of POT plans and `ucoot_cost`.

Protocol rows for both methods live in `inst/bench/thresholds.json`
and `inst/bench/run_protocol.R` (n = 8 / 12). A conservative `-O2`
protocol run with seed `20260816` accepted all four UCOOT rows.
Those smokes are quality guards, not speed claims.
