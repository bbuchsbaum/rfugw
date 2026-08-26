# Entropic partial Sinkhorn benchmark

Local evidence was collected on 2026-08-24 with
`benchmark_partial_sinkhorn.R`; the machine-readable result is
`partial-sinkhorn-baseline.csv`. This is performance and dispatch evidence, not
a substitute for the mathematical tests.

All four cases were certified for feasibility, independent objective
recomputation, KKT conditions, and primal-dual gap. On the 12 by 10 moderate
fixture, scaling and log Dykstra both took 130 checked cycles and returned the
same residual (`8.56e-10`). Scaling used about 4.36 MB of recorded allocations
and 22 ms; log used 6.65 MB and 116 ms. Resuming the certified log state took
one cycle, 0.17 MB, and 1 ms. These single local timings are descriptive, not a
cross-platform speed guarantee.

The adversarial auto case had shift-invariant dynamic range 1249.7, above the
declared scaling threshold 100, and selected the genuine log backend. It
certified in 130 cycles with residual `7.99e-10`. The benchmark records the
requested/effective backend, warm-start flag, iterations, feasibility and
objective residuals, duality gap, runtime, allocation count/bytes, seed,
platform, R version, and source commit.
