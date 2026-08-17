# Contributing to rfugw

rfugw is a focused optimal-transport library. Changes should make a solver
more correct, more honest, faster under a quality gate, or easier to use
without expanding scope toward a POT clone.

## Compatibility

Follow `inst/solver-contract.md`.

- Do not silently ignore a newly accepted parameter.
- Do not advertise a choice the implementation rejects.
- Unregularized GW/FGW is a conditional-gradient procedure, not a global
  optimizer. Do not write “exact solution” without that qualification.
- Primary names are stable. POT-style names are aliases.
- After 0.1, breaking changes need a NEWS entry and a minor version bump.
  Experimental APIs must stay labeled experimental.

## Tests

Every behavior change needs a regression test.

- Flagship solvers: validation, warm starts, symmetry, convergence status,
  and at least one invariant (nonnegativity, marginals or mass, finite
  objective).
- Reference comparisons: record fixture provenance and the POT version when
  a fixture is generated.
- Do not skip a test to make CI green. If a path is optional, test both
  present and absent dependency behavior.

Run the package tests from the repo root:

```r
devtools::test()
```

## Benchmarks

Speed claims are invalid unless the run meets the documented answer-quality
threshold. The canonical protocol is `inst/bench/PROTOCOL.md`.

```
Rscript inst/bench/run_protocol.R
```

- Comparisons must share data, seed, threads, stopping rules, and thresholds
  from `inst/bench/thresholds.json`.
- Solver-only and end-to-end time (and peak memory) are separate columns.
- Invalid-quality rows cannot become baselines or speedups.
- Record commit, compiler flags, thread counts, seed, and hardware
  (`results/current/meta.json`).
- Threshold changes need reviewed evidence, not a silent edit.
- Scratch under `inst/bench/results/scratch/` is not a baseline.

## Generated documentation

Rd files are generated from roxygen comments. After changing exported
signatures or `@return` fields:

```r
devtools::document()
```

Do not hand-edit `man/*.Rd` or `R/RcppExports.R` / `src/RcppExports.cpp`.
After C++ signature changes run `Rcpp::compileAttributes()`.
`NAMESPACE` is maintained by hand; generate Rd with
`roxygen2::roxygenise(roclets = c("rd", "collate"))`.

Vignettes must rebuild from a clean library with declared dependencies
only. Guard optional examples on Suggests packages.

CI helpers (not in the tarball; `tools/` is Rbuildignored):

- `tools/check-generated.R` — Rcpp attributes and Rd coverage
- `tools/check-docs.R` — pkgdown reference coverage of exports
- `tools/smoke-flagship.R` — serial / sanitizer smoke suite

## Evidence for 0.1

Release-blocking work is tracked in mote. Close a ticket only when its
acceptance criteria have tests or recorded check artifacts. Local checks,
hosted CI, and publication evidence stay separate.

Local release-quality gate (builds the clean tarball, docs/Rcpp
consistency, R CMD check, and a skip inventory):

```
Rscript tools/release-gate.R
```

Writes `.gate/latest.json` and `.gate/latest.txt`. Hosted evidence is
GitHub Actions (`R-CMD-check`, `sanitizer`, `pkgdown`). A local pass is
not publication evidence.
