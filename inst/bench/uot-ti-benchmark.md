# Translation-invariant sparse KL-UOT benchmark

Run `Rscript inst/bench/benchmark_uot_ti.R` from the repository root. The
machine-readable result is `uot-ti-baseline.csv`.

The benchmark solves the same 300 by 300 ring-support problem in two ways. The
dense reference assigns cost 25 to forbidden edges; at epsilon 0.3 those edges
have negligible recovered mass. The sparse solve represents only the 2,400
allowed edges. Both use unequal finite measures, the same KL penalties, the
implicit operator result, and full independent certification.

Reported phases are setup, native solve, independent certification, operator
application, and explicit materialization. Peak resident memory is sampled
from a forked child process with `ps`; an unavailable sample is recorded as
`NA`, never as zero. Because RSS includes the R runtime and inherited package
state, the table also reports deterministic peak problem-storage estimates from
the input, CSR plus CSC working arrays, potentials, signals, outputs, and (only
for that phase) the dense materialization. These estimates exclude allocator
and interpreter overhead and are labeled as estimates. Materialization is an
intentional dense allocation boundary, not part of sparse solve or apply.

Answer quality is the maximum absolute plan difference and absolute regularized
objective difference against the dense large-cost reference, alongside each
solver's KKT residual and primal-dual gap. This benchmark supports only this
specified regime; it does not make a universal sparse speedup claim.
