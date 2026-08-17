# Partial GW / FGW hot-path profile

Evidence for `bd-01M05QY883AD27TYRRPMV7WKGG`. Times are median
milliseconds from `Rscript inst/bench/profile_partial.R` on the
conservative `-O2` profile, one thread, seed `20260816`, 3 reps.
Exact rows use 30 CG iterations. Entropic rows use 40 outer
iterations. Problems come from `bench_make_problem("fgw", n, seed)`
with mass `m = 0.7`.

Exact (native LP vs lpSolve) and entropic (C++ outer vs R outer)
are profiled separately because they have different bottlenecks.

## Exact path

Each CG iteration solves a dummy-extended balanced transport LP
of size `(n + nb_dummies)^2`. That LP dominates the runtime.

| n | path | median ms |
|---|---|---|
| 16 | `cpp_transport` | 6 |
| 16 | `lp_matrix` | 7 |
| 24 | `cpp_transport` | 21 |
| 24 | `lp_matrix` | 8 |
| 32 | `cpp_transport` | 98 |
| 32 | `lp_matrix` | 32 |

n = 16 is inside the 1 ms clock. From n = 24 the native simplex
is slower than lpSolve (about 3× at n = 32). The default stays
`cpp_transport` so exact partial has no extra dependency.
Callers who have `lpSolve` and care about exact-path speed can
pass `lp_solver = "lp_matrix"`.

What improved on this path:

- The default CG loop, tensor, linesearch, and dummy-extended LP
  now run in `cpp_partial_fgw_exact_square` instead of R + lpSolve
- GW wrappers pass an empty feature cost when `alpha = 1`
- Mass `m` and row/column upper bounds stay within tolerance
  against the lpSolve CG core

## Entropic path

The inner Dykstra-like projection was already C++. The remaining
work is the outer gradient and the R call overhead.

| n | path | median ms |
|---|---|---|
| 64 | C++ outer, empty `M` | 0 |
| 64 | R outer | 1 |
| 64 | C++ outer, explicit zero `M` | 1 |
| 96 | C++ outer, empty `M` | 1 |
| 96 | R outer | 2 |
| 96 | C++ outer, explicit zero `M` | 2 |

n ≤ 48 is at timer resolution and is not used as evidence.

What improved on this path:

- The outer loop now runs in `cpp_partial_fgw_entropic_square`
- Empty `M` at `alpha = 1` avoids an `ns × nt` GEMM term
- The inner `q2` / `q3` updates now use the pre-multiply kernels,
  matching POT and the R reference. The post-multiply form overflowed
  at high mass (for example n = 6, `m = 0.8`)

## POT evidence

`inst/extdata/fixtures/partial_pot_multisize.json` (POT 0.9.6.post1,
seed 20260816) covers n ∈ {6, 6, 8, 10} and m ∈ {0.5, 0.8, 0.6, 0.7}.
Entropic partial GW matches POT to about 1e-8. POT's entropic
partial FGW gradient adds the scalar `(1 - alpha) * sum(G * M)`
instead of `(1 - alpha) * M`, so FGW is checked against
`ot_fgw_square`, not POT's FGW log value.

Protocol rows for `partial_gromov_wasserstein`,
`partial_fused_gromov_wasserstein`, and
`entropic_partial_gromov_wasserstein` live in
`inst/bench/thresholds.json` and `inst/bench/run_protocol.R`.
Those n = 8 / 10 smokes are quality guards, not speed claims.
A conservative `-O2` protocol run with seed `20260816` accepted
all six partial rows.
