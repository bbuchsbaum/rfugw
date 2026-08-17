# Sampled and approximate GW: certified envelope

Evidence for `bd-01M05QY8NJRGY4NDT2YQ7ZHQ4H`. These methods stay
experimental. The only 0.1 claims are the tiny-versus-full budget
quality envelope, rank reconstruction, and input-memory scaling below.
They are not certified substitutes for dense
`entropic_gromov_wasserstein()`.

Regenerate numbers with:

```sh
RFUGW_RLIB=/path/to/.tmp-lib Rscript inst/bench/sampled_budget_curves.R \
  inst/bench/results 20260816 5
```

Numbers below are from a conservative `-O2` local install, one thread,
seed `20260816`. CSV copies live in `inst/bench/results/` (gitignored).

## Sampled GW versus dense

Problem: `n = 16`, 2-D Gaussian clouds, costs scaled to `[0, 1]`.
Reference is dense `entropic_gromov_wasserstein()` at `epsilon = 0.1`
(square-loss GW `2529.045`). Each budget uses 5 seeds and 40 outer
iterations. Gap is `|ot_gw_square(sampled) - ot_gw_square(dense)|`.
Frobenius is `||T_sampled - T_dense||_F`.

| nb_p | nb_q | median gap | median Frobenius | median GW | median ms |
|---|---|---|---|---|---|
| 2 | 1 | 1.83 | 0.192 | 2530.877 | 2 |
| 4 | 1 | 2.20 | 0.217 | 2531.242 | 1 |
| 8 | 2 | 1.94 | 0.198 | 2530.986 | 2 |
| 16 | 16 | 0.92 | 0.112 | 2529.960 | 2 |

Certified claim: a full budget `(ns, nt)` is closer to the dense
reference than a tiny budget such as `(2, 1)`, in both square-loss GW
and plan Frobenius distance. At `n = 16` the full budget cut the median
gap from 1.83 to 0.92 and the median Frobenius from 0.192 to 0.112.
The full-budget run is deterministic across seeds.

Intermediate budgets are **not** certified as monotone. `(4, 1)` and
`(8, 2)` sit inside the seed scatter of `(2, 1)`. Tiny usable budgets
still return valid couplings.

Times at `n = 16` are at timer resolution and are not a speed claim.
Speed versus POT is a separate sparse-sampled performance gate.

Unusable budgets:

- `nb_samples_grad < 1`, non-finite values, or a vector that is not
  length 1 or 2: error
- length-2 counts above `ns` or `nt`: warning and clamp
- a scalar larger than `ns` still uses the POT-style remap
  `(ns, nb %/% ns)`; if that target count exceeds `nt`, it warns and
  clamps

## Rank approximation

`lowrank_gromov_wasserstein_samples()` factorizes a dense entropic GW
plan by truncated SVD. It is not a factorized solver.

| rank | relative Frobenius error | SVD ms |
|---|---|---|
| 1 | 0.562 | 7 |
| 2 | 0.210 | 0 |
| 4 | 0.035 | 0 |
| 8 | 0.0014 | 0 |
| 16 | ~0 | 0 |

Certified claim: reconstruction error decreases as rank grows up to
`min(ns, nt)`. Rank `< 1` errors. Rank above `min(ns, nt)` warns and
clamps. SVD time at `n = 16` is not a speed claim.

## Memory scaling

The returned plan is always dense (`n x n`). The certified input-memory
claim is only about structure storage. Bytes are `object.size()` of the
inputs, not peak allocator traffic.

| n | dense C1+C2 | coords X1+X2 | sparse 6-NN graphs | 8-D embeddings | plan |
|---|---|---|---|---|---|
| 32 | 25552 | 1968 | 6656 | 2264 | 8408 |
| 64 | 82896 | 3504 | 10176 | 4312 | 32984 |
| 96 | 173008 | 5040 | 13624 | 6360 | 73944 |

- Dense sampled GW stores two `n x n` costs plus the plan
- Coordinate-native sampled GW stores two `n x d` point clouds plus the
  plan
- Graph sampled GW stores sparse similarities, `k` diffusion
  coordinates, and the plan

At `n = 96`, coordinates are about 3% of two dense costs, and a sparse
graph plus 8-D embeddings is about 12%. The plan still grows as `n^2`
and already exceeds the coordinate inputs at `n = 32`.

## What is not claimed

- Sampled GW is not a certified match to POT or to dense GW at a fixed
  intermediate budget
- Quality does not improve at every budget increment
- Speed versus POT is out of scope here
- `lowrank_gromov_wasserstein_samples()` is not a low-rank solver
