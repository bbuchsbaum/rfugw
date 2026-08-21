# rfugw benchmark protocol

This is the only accepted way to generate speed evidence for rfugw 0.1
(`bd-01M05QY7SXKJR408SD2MPEJC1T`). A faster wrong answer is not a
benchmark. Later Phase D tickets (semirelaxed, partial, UCOOT, sampled
quality, threading, flagship CI) must use this protocol.

## Entry point

From the package root, after installing the conservative profile:

```
Rscript inst/bench/run_protocol.R
```

Optional arguments: `reps` `seed` `out_dir` `threads` `suites`.

`suites` is a comma-separated list from
`linear,fgw,fugw,semirelaxed,partial,ucoot,sampled`. Empty means all.
`RFUGW_PROTOCOL_SCALE` is `full` (default), `pr`, or `nightly` and
selects the size grid.

Writes:

- `inst/bench/results/current/meta.json` — environment and versions
- `inst/bench/results/current/runs.csv` — one row per (suite, method, n)
- `inst/bench/results/current/quality.csv` — quality checks for those rows

Every row declares one of two evidence classes:

- `certified_comparison`: eligible for answer-quality and timing comparison
  only when `status == "converged"` and every formulation certificate passes;
- `fixed_budget_performance`: a runtime sentinel at a declared iteration or
  approximation budget. It is always `certified = FALSE` and
  `comparison_eligible = FALSE`, even when its regression checks pass.

`max_iter`, `inner_failure`, and `experimental` are never aliases for
`converged`. FUGW's loose fixed-budget row and sampled GW's budget curve belong
only to the second class.

Scratch under `inst/bench/results/scratch/` is not a baseline. The entry
point archives any previous `current/` tree into `scratch/` before
writing. `inst/bench/results/` is gitignored.

Certified baselines use the conservative compile profile in
`inst/build-profiles.md`. `RFUGW_FAST_FLAGS` runs are labeled
`profile=fast` and cannot replace a conservative baseline.

## Shared comparison contract

Every compared method on a row-group must share:

| Field | Source |
|---|---|
| Problem data | `bench_make_problem()` with the same `kind`, `n`, and `seed` |
| Seed | CLI / `meta$seed` |
| Threads | `OMP_NUM_THREADS` and BLAS pinned to `meta$threads` |
| Stopping | `inst/bench/thresholds.json` `stop` block for that method |
| Quality | same file, `quality` block |

Do not mix seeds, sizes, `epsilon`, `max_iter`, or thread counts inside
one comparison.

## Timing and memory

Times are milliseconds.

| Column | Meaning |
|---|---|
| `prepare_ms` / `setup_ms` | Cost/structure construction only |
| `solve_ms` | Solver call only (median over timed reps after warmup) |
| `e2e_ms` | `prepare_ms + solve_ms` |
| `mem_solve_bytes` | Peak allocation during the solver call when `bench` is available; otherwise `NA` |
| `mem_e2e_bytes` | Peak allocation of prepare+solve when available |

Warmup is required and is not a timed rep. Default: 1 warmup, 3 timed
reps. Report the median of timed reps.

## Quality gate

A certified-comparison run is eligible only if every check in its `quality`
block passes:

- `status` is exactly `converged`; `converged`, `feasible`,
  `objective_consistent`, and `objective_components_consistent` are true
- every required nested solver has `inner_converged == TRUE`
- residuals / mass within thresholds
- independent objective (`ot_linear_cost`, `ot_fgw_square`,
  `ot_gw_square`, …) matches the reported value within
  `atol + rtol * |expected|`
- semirelaxed methods also require `rowSums(plan) == p`

Invalid runs are written with `valid = FALSE` and a `reject_reason`.
They must not be used as a baseline, a speedup, or a CI threshold
update. `run_protocol.R` exits nonzero if any intended baseline row is
invalid.

Fixed-budget rows use `performance_regression`, never `quality`, and cannot set
`certified` or `comparison_eligible`. Their time caps detect regressions at the
same budget; they do not support numerical or convergence claims.

## Approximation curves

Sampled GW emits three rows per problem with `(nb_p, nb_q)` budgets `(2,1)`,
`(4,2)`, and `(8,4)`. `quality.csv` records the exact-CG reference objective,
absolute error, and relative error for each budget; `runs.csv` records the
corresponding setup, solve, end-to-end, and allocation measures. A speed-only
sampled row is inadmissible.

## Recorded environment

`meta.json` always includes:

- package version, git commit, R version and platform
- compiler, BLAS, and `RFUGW_*` flag environment
- `Sys.info()` hardware fields
- `OMP_NUM_THREADS`, BLAS thread env
- seed, warmup count, timed repetitions
- timestamp (UTC)

Rows separately record requested/effective/compute precision and
requested/effective thread counts. Environment metadata stays in `meta.json`;
answer quality stays in `quality.csv`; timing and allocation stay in
`runs.csv`.

## Existing scripts

`benchmark_suite.R`, `benchmark_sparse_sampled.R`, and the nightly
guard remain. New comparisons must go through `protocol.R` helpers so
quality and metadata stay uniform. Do not treat a CSV from those older
scripts as a 0.1 baseline unless it was regenerated through this
protocol.

## Hosted CI gates

Quality-controlled hosted gates use this protocol, not the older
suite CSVs.

| Gate | Workflow | Suites | Scale |
|---|---|---|---|
| Fast PR | `.github/workflows/flagship-gate.yml` | FGW, FUGW, semirelaxed | `pr` |
| Nightly | `.github/workflows/flagship-gate-nightly.yml` | FGW, FUGW, semirelaxed, partial, UCOOT, sampled, plus a 1-vs-2 thread smoke | `nightly` |

The sparse-sampled POT speed gates stay separate. They are not a
substitute for these quality gates.

Entry point:

```
inst/bench/run_flagship_gate.sh inst/bench/results/current \
  fgw,fugw,semirelaxed pr 1 20260816 1
```

Artifacts in the output directory: `meta.json`, `runs.csv`,
`quality.csv`, `gate_report.md`. Nightly also writes `threads.csv`.
Hosted jobs upload that tree.

### Regression criteria

- Quality is the hard gate. A row is invalid unless
  `inst/bench/thresholds.json` checks pass. Invalid rows cannot update
  a cap or count as a pass.
- Time caps in `inst/bench/ci_time_caps.json` are median `solve_ms`
  after warmup. A row fails only if that median exceeds the cap. Caps
  are slack for GitHub-hosted runner noise, not tight speed claims.
- Changing `thresholds.json` or a time cap requires a protocol artifact from
  the same commit in the PR, a retained entry in `threshold-history.json`, and
  an explicit reviewer note. The protocol verifies the recorded threshold MD5,
  so silent threshold edits are rejected. Schema-v1 CSVs remain historical
  performance data and are never reclassified as schema-v2 certification.

### Exit codes

`run_flagship_gate.sh` and `gate_protocol.R` distinguish failures:

| Code | Meaning |
|---|---|
| 0 | Quality and caps passed |
| 1 | Solver / quality / time-cap regression |
| 2 | Infrastructure: missing artifacts, unreadable JSON, protocol did not write |
