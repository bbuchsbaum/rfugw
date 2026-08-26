# KL structural loss decision

## Decision: defer

rfugw 0.1 supports only `loss_fun = "square_loss"`. It deliberately rejects
`"kl_loss"` before computation. This is an estimand decision, not an absent
string dispatch: adding an implicit numerical floor would silently define a
different directed loss for the zero-diagonal distance matrices used by the
package.

The decision was evaluated for `bd-01M05QY99ZQHFQ710DVAR7APX8`. The retained
prototype is `inst/numerical-trust/kl-structural-loss-prototype.R`; reproduce its differential
and timing evidence with:

```sh
Rscript inst/bench/evaluate_kl_structural_loss.R \
  inst/bench/kl-structural-loss-baseline.csv
```

## Estimand and domain

The directed structural loss is

```text
L(a, b) = a log(a / b) - a + b.
```

It is finite for `a = 0, b >= 0` under the continuous convention, but is
infinite when `a > 0, b = 0`. It is asymmetric. Its natural use case is the
comparison of strictly positive similarity, kernel, rate, or intensity
matrices whose ratio has a scientific interpretation. Ordinary metric cost
matrices instead have a zero diagonal. Under a positive coupling, their GW
four-index objective generally contains positive source distances compared
with zero target distances and is therefore infinite.

POT 0.9.7 exposes this loss through the factorization
`f1(a) = a log(a) - a`, `f2(b) = b`, `h1(a) = a`, `h2(b) = log(b)` and adds
`1e-18` inside logarithms in its implementation. The rfugw prototype reproduces
that factorization, but calls the constant `log_floor` because it changes the
objective rather than merely stabilizing an algebraically equivalent
calculation.

## Differential and numerical evidence

The prototype compares an independent O(ns^2 nt^2) enumeration with the
factorized matrix expression. On strictly positive and near-zero fixtures they
agree to floating-point tolerance for every checked floor. On zero-diagonal
metric costs the unfloored objective is infinite; positive floors make it
finite, but its value moves materially across `1e-18`, `1e-15`, and `1e-12`.
The retained CSV records the objective differential, the floor, timings, R
version, and platform.

This evidence supports the rejection boundary, not a performance claim. The
small timings only confirm that the factorization is executable; they do not
establish an rfugw backend baseline.

## Risks and compatibility

- A hidden floor changes objective values and gradients, especially near zero.
- The directed loss does not inherit the symmetry of square structural loss.
- Cost rescaling is not innocuous because it interacts with the logarithmic
  ratio and any floor.
- rfugw's native objective, gradient, line search, independent oracle, and
  convergence certificates currently implement square loss as one coherent
  algebra. Enabling only a POT-shaped argument would fracture that contract.
- Existing POT familiarity is not sufficient compatibility: identical loss
  names would still disagree if floor, domain, or objective constants differ.

## Reconsideration gate

A future implementation requires a separate positive-structure contract with
either strict positivity or an explicit user-visible floor in the estimand. It
must add an independent four-index oracle, analytic gradient checks, native/R
and POT differentials under the same convention, near-zero adversarial tests,
scale-law documentation, objective-component certificates, and retained
quality-gated benchmarks. Until that complete path exists, all public GW, FGW,
partial, semirelaxed, sampled, and barycenter entry points reject `"kl_loss"`.
