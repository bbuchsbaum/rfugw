# rfugw vignette gauntlet

Benchmark: [Introduction to dplyr](https://dplyr.tidyverse.org/articles/dplyr.html)

## Quality bar

- Begin with the reader's problem, data, and intended result.
- Introduce one coherent workflow before variants or tuning.
- Keep code compact and output purposeful.
- Support computational claims with executable, quiet checks.
- Match prose to the exported API and current runtime behavior.
- Render every affected vignette and inspect the resulting HTML.

## Progress

- [x] Inventory the repository and protect pre-existing work.
- [x] Fetch and inspect the live benchmark.
- [x] Map the reader journey across the five-vignette suite.
- [x] Complete the first rewrite and independent critique of each vignette.
- [x] Build a source tarball and rebuild the complete suite under `R CMD check`.
- [x] Address convergence, diagnostic-contract, copyability, and compatibility findings.
- [x] Pass the final independent comparison for all five vignettes.

## Current focus

Complete. Every visible chunk sequence runs without hidden setup, every
showcased solver state is either certified or explicitly qualified, and the
independent exit critic selected each rfugw vignette over the benchmark.

## Evidence so far

- Clean temporary package installation: passed.
- Source tarball creation, including all five vignettes: passed.
- `R CMD check --no-manual` on the tarball under an explicit UTF-8 locale:
  `Status: OK`.
- Initial shell locale was unusable and caused a metadata-stage check failure;
  this cleared under the installed `en_US.UTF-8` locale.
- Product-controlled browser: unavailable. The only connected browser is a
  user Chrome profile, which this task is not authorized to use.
- Final source-tarball `R CMD check --no-manual`: `Status: OK`.
- Final independent benchmark result: rfugw won all five binary comparisons.
- Final browser-process audit: no automated top-level browser processes found.

## Discovery

`rfugw` aligns weighted objects using linear OT, GW/FGW, unbalanced variants,
barycenters, and multiset template alignment, with the main solvers backed by
C++ kernels.

The central teaching surface is `ot_sinkhorn()`, `ot_emd()`,
`ot_sinkhorn_unbalanced()`, `gromov_wasserstein()`,
`entropic_gromov_wasserstein()`, `fgw_exact_cg()`, `fgw_entropic()`,
`fugw_kl()`, `fgw_barycenters()`, and `multialign_fit()`.

There is no packaged teaching dataset. The current guides use small seeded
point clouds; `inst/extdata/fixtures/` contains frozen validation fixtures, not
reader-facing data. Relevant optional packages are `Matrix`, `RSpectra`,
`lpSolve`, `knitr`, and `rmarkdown`.

Reader journey:

1. `rfugw`: understand structure matrices, feature costs, and couplings.
2. `linear-ot`: solve the simpler shared-coordinate case.
3. `solver-guide`: choose balanced, partial, semirelaxed, or unbalanced GW/FGW.
4. `barycenters`: learn a representative metric-measure object.
5. `multiset-alignment`: align many inputs through a shared template.
