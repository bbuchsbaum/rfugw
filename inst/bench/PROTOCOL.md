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

`mem_solve_bytes` and `mem_e2e_bytes` are R-visible allocations. They are not
native peak memory and must never be labeled RSS. The Moment-FUGW protocol
below uses a fresh process and the operating system's resource accounting for
authoritative peak RSS.

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

## Moment-FUGW admission protocol

Factorized and multiscale FUGW use structured certificates rather than the
legacy literal `status == "converged"` check. Run their dedicated protocol from
the package root:

```
Rscript inst/bench/run_moment_fugw_benchmark.R smoke
Rscript inst/bench/run_moment_fugw_benchmark.R full \
  inst/bench/results/moment_fugw installed
Rscript inst/bench/run_moment_fugw_core_certification.R local \
  inst/bench/results/moment_fugw_core source
Rscript inst/bench/run_moment_fugw_robustness.R local \
  inst/bench/results/moment_fugw_robustness source
```

The arguments are `profile`, `results_root`, and `package_mode`. `package_mode`
is `source` for development and `installed` for a release candidate. Only a
full, installed, conservative-profile run can become release evidence.
The runner also requires an exact 40-hex Git commit and a successfully queried,
clean repository status. A failed or unavailable Git command records unknown
provenance and fails the release gate; it is never interpreted as an empty,
clean status. Expected build outputs are gitignored, so a canonical release job
can preserve this check while installing and measuring its verified tarball.

The versioned design and gates live in
`inst/bench/moment-fugw-thresholds.json`. Its history is
`inst/bench/moment-fugw-threshold-history.json`. The initial contract is
deliberately `candidate`: a full artifact and independent review must be added
to history before any row can have `threshold_update_eligible = TRUE`.

The versioned core contract is
`inst/bench/moment-fugw-core-contract.json`. Its release profile regenerates
the complete 12-case dense/direct-oracle matrix (through 128 by 127, including
exact-zero and positive near-zero weights), the frozen ten-seed multiscale
efficacy study, and 27 named oracle, final-KKT, geometry, hierarchy, and
multiscale cases. It retains raw CSVs, binds every participating test/helper/
executable hash, binds the full case specification in
`moment-fugw-accuracy-design.csv`, and pins BLAS and OpenMP to one thread. A
summary boolean is
not evidence by itself: cross-platform admission independently recomputes the
thresholds and counts from those raw artifacts.

The kernel admission profile is `inst/bench/profile_moment_fugw.R`. It times
the whole fixed-work solve and separately profiles TI sweeps, plan statistics,
moment reduction, plan difference, objective recomputation, both final KKT
audits, apply, and adjoint. Retained before/after medians and fresh-process RSS
checks live in `inst/bench/moment-fugw-kernel-evidence.csv`; the CSV is
candidate source-tree evidence, not a release baseline. Its gate requires at
least 15 percent end-to-end improvement at both n = 1000 and n = 2000, no more
than 5 percent regression at retained smaller endpoints, and no absolute
peak-RSS regression.

The protocol writes:

- `current/meta.json`: commit, compiler, BLAS, hardware, seeds, pinned threads,
  profile, threshold checksum, RSS backend, and the raw-artifact manifest;
- `current/runs.csv`: every warmup, measurement, rejected result, raw time, and
  eligibility decision;
- `current/summary.json`: fitted slopes, largest certified endpoint, claim
  ceiling, and admission gates;
- `current/raw/`: the JSON input, JSON result, stdout, and raw resource report
  from every fresh worker process.

The historical retained candidate installed-artifact replay is under
`inst/bench/evidence/moment-fugw-full-installed-candidate-20260827/`. Its small
archive contains all 228 pre-allocation-audit raw worker receipts plus the
unmodified 56-row summary, metadata, and run table; cached fit objects are
excluded because they are execution scratch. It predates the mandatory
R-visible allocation curve and therefore cannot satisfy the current contract.
The current full profile has 60 rows and 244 raw receipts: one package-load
baseline plus 60 measured or warmup workers, each retaining its spec, result,
stdout, and native resource report. Historical replay and current contract are
kept distinct. Neither can become a release baseline without exact clean
same-commit provenance.

The amended installed candidate is retained separately under
`moment-fugw-full-installed-rprofmem-candidate-20260827/`. Its exact checked
tarball completed all 60 rows, all 244 raw receipts, the native-RSS and
R-visible-allocation gates, and a 512-node certified endpoint. Its primary
runtime exponent is 1.943 and incremental RSS exponent is 1.324. It remains
non-release evidence because the tarball was built from a dirty working tree;
the attestation records that boundary rather than treating the reported parent
commit as the source identity.

Previous `current/` artifacts move to `scratch/`. They remain diagnostic and do
not silently become baselines.

### Structured quality gate

A `fugw_factorized` certified row requires all of the following:

- status `converged_stationary` or `converged_approximate_stationary`;
- final inner UOT, consecutive plan/moment stationarity, and both final
  post-rescale block KKT audits certified;
- equal coupling mass, objective/component consistency, feasibility, and
  numerical finiteness;
- certified geometry and an exact feature cost for the benchmark fixture;
- complete full implicit support with zero omitted-kernel bound;
- implicit, non-materialized sample and feature plans;
- the memory contract naming `Cx`, `Cy`, `M`, `P`, and `Q` as prohibited dense
  allocations.

`fugw_multiscale` additionally requires every level's inner and outer
certificate plus the computed hierarchy-transfer and geometry audits. A
`converged_final_level_only`, geometry-uncertified, hierarchy-uncertified,
nonfinite, or otherwise nonstationary result remains in `runs.csv` with a
reason, while its reportable timing is `NA`.

Historical inner failures may be retained in a successful result; the gate
uses the final recovered two-sided certificate. Global optimality is never
claimed.

### Evidence classes and curves

Moment-FUGW keeps four claims separate:

- `fixed_work_scaling`: exactly one outer update and one iteration per UOT
  block, used only to characterize full-support scaling. These rows are
  `scaling_eligible` but never certified, comparison eligible, timing-comparison
  eligible, or threshold-update eligible;
- `end_to_end_convergence`: independent certified solves, used for actual
  converged runtime and the largest supported endpoint;
- `rank_scaling`: bilinear structure ranks 4, 16, 64, and 256 at fixed domain
  size;
- `operator_throughput`: apply and adjoint batches 1, 10, 100, and 1000 from a
  certified serialized implicit operator.

The full fixed-rank size curve is 625, 1250, 2500, and 5000. Its primary
asymptotic fit uses the largest three points and must fall in 1.8--2.2,
reflecting the intended full-support quadratic arithmetic. The all-four-point
fit is also mandatory and retained as an overhead-sensitive diagnostic; it is
never silently discarded. This is not a claim of subquadratic runtime.
End-to-end convergence is fitted and reported separately.

### Multiscale transfer efficacy

`inst/bench/benchmark_moment_fugw_multiscale.R` freezes a ten-seed planted
rigid-permutation hierarchy. Every transferred and cold fine solve shares the
same final epsilon, objective weights, marginal penalties, geometry, stopping
tolerances, tile size, and full-support certificate. The retained
`inst/bench/moment-fugw-multiscale-evidence.csv` records both work counts,
objectives, plan-action differences, residuals, and a held-out map score.

Admission requires at least 80 percent of cases to reduce fine-level inner
work, a median reduction of at least 20 percent, transferred objectives no
worse than 1e-4 relative to certified cold solves, and held-out score
noninferiority within one percent. The frozen candidate currently has fewer
fine-level inner iterations in all ten cases and a 34.6 percent median
reduction. This is efficacy evidence for the planted fixture, not a general
runtime or release-hardware claim.

Every level also declares either `certified_endpoint` or
`budgeted_warm_start`. The default requires all levels to be certified. A
successful finest solve following explicitly budgeted intermediate work gets
the distinct status `converged_final_level_with_budgeted_warm_starts`; it is
not relabeled as all-level stationary. Potential prolongation is gauge
recentered, while epsilon changes preserve dual potentials in cost units.

### Planted image and cortical validation

The scientific harness is governed by the frozen
`inst/bench/moment-fugw-validation-protocol.json`. Tuning and evaluation seeds
are disjoint, every requested row is retained, and solver failures are data
rather than silently omitted cases. Run it from the package root with:

```
Rscript inst/bench/run_moment_fugw_scientific_validation.R tuning \
  inst/bench/results/moment_fugw_scientific/tuning source
Rscript inst/bench/run_moment_fugw_scientific_validation.R evaluation \
  inst/bench/results/moment_fugw_scientific/evaluation installed
```

The nightly numerical-trust workflow replays the tuning split against the
installed package and uploads its CSVs and JSON summary. The frozen evaluation
split is reserved for threshold-independent candidate or release evidence.
The shell launcher with the same basename is only a convenience wrapper around
the direct R entry point.

The scientific summary separates task success from release admissibility.
Image, surface, and geometry-rank gates may pass while
`admission_evaluable` remains false. Release admissibility additionally
requires the installed evaluation split, an exact Git commit, and a
successfully queried clean status; failed Git commands record `git_dirty` as
unknown rather than false.

The image suite plants deformation, nonuniform weights, feature noise, and
100, 80, 60, and 40 percent overlap. Its partial-overlap gate requires
unbalanced FUGW to have strictly lower median false occluded-source mass than a
balanced FGW comparator at every partial-overlap level, while preserving
endpoint accuracy and recovered overlap mass. The retained evaluation has all
40 requested rows and passes all three paired gates.

The cortical suite uses 36 deterministic landmarks derived from the left
fsaverage6 pial/inflated mesh. Its train and held-out maps are planted
geodesic radial-basis channels; this is a controlled alignment benchmark, not
evidence about observed cortical physiology. The exact-factor evaluation has
five certified rows with held-out correlation 0.99942--0.99949. A separate
16-landmark dense fixture gates every objective component and sample-plan
action at tolerance 1e-6; its retained maximum error is 2.27e-11. The rank
curve reports geometry representation error independently of numerical solver
and coupling discrepancy.

Runtime validation reads retained hexadecimal-double CSV derivatives rather
than loading a legacy serialized surface object. The 16- and 36-landmark
semantic fingerprints, resource MD5s, raw-surface MD5s, source repository
commit, and deterministic derivation are recorded in
`moment-fugw-fsaverage6-fixture-provenance.json`. The maintainer-only
`generate_moment_fugw_cortical_fixture.R` regenerates those resources from raw
FreeSurfer triangle surfaces and fails if CSV round-trip changes the semantic
fingerprint. This replay hardening changed no estimand, threshold, seed, or
evaluation result.

### Native peak memory

Each action runs in a new `Rscript` process. On macOS the backend parses byte
units from `/usr/bin/time -l`; on Linux it converts the KiB units from
`/usr/bin/time -v` to bytes. The raw report is always retained. A matching
fresh-process package-load baseline is subtracted for the primary memory-growth
curve; absolute peak RSS and its slope are also reported.

The full profile passes the memory gate only if the incremental peak-RSS
doubling exponent is below 1.5 and every solve, objective recomputation, apply,
and adjoint action passes the dense-allocation sentinel. A runtime memory
contract alone is not enough: the cost callbacks abort on a full block request,
and the worker verifies that neither returned coupling was materialized.

Peak RSS includes native allocations. A separate fresh solve and each
objective/apply/adjoint worker are also profiled with base `Rprofmem`; the
reported `r_visible_allocation_bytes` is cumulative R-visible allocation, not
peak memory and never a substitute for native RSS. All four actions at all four
full-profile sizes are mandatory. Missing R allocation or OS RSS is an
infrastructure failure, never a zero or an inferred measurement.

### Cross-platform candidate and promotion receipts

`tools/numerical-trust/verify-moment-fugw-admission.R` recomputes the admission
decision from the canonical source manifest, both full installed native-RSS
artifacts, all three installed frozen scientific artifacts, and all three
installed core and metamorphic/adversarial robustness artifacts. The core
replay is governed by `moment-fugw-core-contract.json`; it recomputes all 12
dense/oracle rows, exact-zero and positive near-zero coverage, all ten
multiscale work/accuracy rows, and all 27 named AC1--AC5 cases from raw CSVs.
The robustness replay
is governed by `moment-fugw-robustness-contract.json`; it binds the source-test
hashes, recomputes case outcomes from `cases.csv`, requires all 15 declared
AC8 properties, and verifies one-thread BLAS/OpenMP settings on Darwin, Linux,
and Windows. The verifier also checks exact commit and clean-tree provenance,
artifact hashes, benchmark design rows, retained raw-worker JSON, R-visible
allocation coverage, runtime and memory slopes recomputed from `runs.csv`,
certified endpoint, frozen protocol checksum, partial-overlap comparator metrics
recomputed from CSV, and dense surface parity. The three hosted
numerical-trust ledgers, canonical tarball copies, and successful `R CMD check`
logs are commit- and digest-bound. Its SHA-256 evidence fingerprint is
deterministic over all retained layers.

Candidate mode can establish only
`technical_candidate_passed_review_pending`. Promotion mode additionally
requires a separately authored review bound to that evidence fingerprint, with
independent numerical, performance, and scientific decisions. The review
schema is in `moment-fugw-independent-review-template.json`. Even a passing
promotion receipt records `supported = false`; changing public status remains
a separate reviewed source change. Global optimality, adaptive support, GPU
scale, and 100k-node performance remain unclaimed.

### Claim ceiling

`summary.json` records `largest_completed_certified_endpoint` and copies it to
`claim_ceiling_n`. Documentation and release notes must stop there. Reaching
5000 under a one-iteration scaling contract does not certify a converged 5000
by 5000 alignment. A larger endpoint counts for release only on documented
installed conservative hardware; source-mode smoke evidence cannot promote the
candidate thresholds. A dirty worktree also fails the release-hardware gate,
even when the package version and reported Git SHA look correct.

## Approximation curves

Sampled GW emits three rows per problem with `(nb_p, nb_q)` budgets `(2,1)`,
`(4,2)`, and `(8,4)`. `quality.csv` records the exact-CG reference objective,
absolute error, and relative error for each budget; `runs.csv` records the
corresponding setup, solve, end-to-end, and allocation measures. A speed-only
sampled row is inadmissible.

## Recorded environment

`meta.json` always includes:

- package version, git commit, dirty-worktree flag and status-entry count,
  R version, and platform
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

The numerical-trust release matrix additionally installs the one canonical
source artifact, runs the full Moment-FUGW native-RSS profile on macOS and
Linux, and replays the frozen scientific evaluation on macOS, Linux, and
Windows. Those hosted artifacts still require independent review before the
candidate threshold history can be promoted. A dependent aggregation job
downloads the five Moment-FUGW artifacts and first emits the fail-closed
technical candidate receipt.

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
