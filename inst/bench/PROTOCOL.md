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

Optional arguments: `reps` `seed` `out_dir` `threads`.

Writes:

- `inst/bench/results/current/meta.json` — environment and versions
- `inst/bench/results/current/runs.csv` — one row per (suite, method, n)
- `inst/bench/results/current/quality.csv` — quality checks for those rows

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
| `prepare_ms` | Cost/structure construction only |
| `solve_ms` | Solver call only (median over timed reps after warmup) |
| `e2e_ms` | `prepare_ms + solve_ms` |
| `mem_solve_bytes` | Peak allocation during the solver call when `bench` is available; otherwise `NA` |
| `mem_e2e_bytes` | Peak allocation of prepare+solve when available |

Warmup is required and is not a timed rep. Default: 1 warmup, 3 timed
reps. Report the median of timed reps.

## Quality gate

A run is valid only if every check in its `quality` block passes:

- `status` is in the allowed set (usually `converged`)
- residuals / mass within thresholds
- independent objective (`ot_linear_cost`, `ot_fgw_square`,
  `ot_gw_square`, …) matches the reported value within
  `atol + rtol * |expected|`
- semirelaxed methods also require `rowSums(plan) == p`

Invalid runs are written with `valid = FALSE` and a `reject_reason`.
They must not be used as a baseline, a speedup, or a CI threshold
update. `run_protocol.R` exits nonzero if any intended baseline row is
invalid.

## Recorded environment

`meta.json` always includes:

- package version, git commit, R version
- compiler / `RFUGW_*` flag environment
- `Sys.info()` hardware fields
- `OMP_NUM_THREADS`, BLAS thread env
- seed, warmup count, timed repetitions
- timestamp (UTC)

## Existing scripts

`benchmark_suite.R`, `benchmark_sparse_sampled.R`, and the nightly
guard remain. New comparisons must go through `protocol.R` helpers so
quality and metadata stay uniform. Do not treat a CSV from those older
scripts as a 0.1 baseline unless it was regenerated through this
protocol.
