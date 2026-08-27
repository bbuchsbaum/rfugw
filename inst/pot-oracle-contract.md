# Live POT differential certification

`rfugw` uses the Python Optimal Transport package (POT) as a live differential
oracle, not as its sole source of numerical truth.  The gate combines POT
outputs with independently recomputed objectives, marginal and mass
feasibility, solver certificates, and formulation-specific invariants.  A POT
result cannot override those checks merely because two plans look similar.

## Replaying the gate

The baseline environment is pinned in
`tools/pot-oracle/requirements-baseline.txt`.  From a source checkout with
`rfugw` already installed, run:

```sh
python tools/pot-oracle/generate.py \
  --profile=pr --channel=baseline \
  --output=.gate/pot-oracle/pot.json
Rscript tools/pot-oracle/compare.R \
  --input=.gate/pot-oracle/pot.json \
  --output=.gate/pot-oracle/rfugw.json \
  --scope=pr
```

The Python receipt records the POT, NumPy, SciPy, and Python versions, the
platform, timestamp, profile, channel, seed, function-specific POT entry point,
parameters, captured diagnostics, per-case input and output SHA-256 digests,
and an aggregate case digest.  It is generated evidence, not a checked-in
golden fixture.  The R receipt adds the installed `rfugw` version and R session,
every tolerance result, strict and monitored result counts, and all active
exceptions.

## Profiles and channels

| Lane | Profile | POT channel | Authority |
|---|---|---|---|
| Pull request and main push | `pr`, one seed | exact baseline pins | required |
| Nightly baseline | `nightly`, three seeds including asymmetric structures | exact baseline pins | required |
| Nightly compatibility | `nightly`, three seeds | bounded latest POT/NumPy/SciPy | informational until promoted |
| Tagged release | `release`, three seeds | exact baseline pins | required against the canonical source artifact |

The baseline channel refuses to run with any POT version other than
`0.9.7.post1`.  The latest channel is allowed to move only within the bounded
requirements file and produces its own versioned receipt.  Latest-channel
failures do not silently change the release oracle; promotion requires review
of the diff and an explicit baseline-pin update.

## Certified formulation matrix

Each seed covers the following families through public, function-specific POT
and installed `rfugw` APIs:

| Family | Strict evidence |
|---|---|
| Balanced linear OT | EMD plan and cost; Sinkhorn scaling/log plans and product-KL objectives; KL-UOT plan, mass, and objective |
| Partial linear OT | exact and entropic plans/objectives; penalized plan, transported mass, and discard-penalty objective |
| Balanced GW/FGW | exact-CG and entropic plans and independently named square-loss objectives, with explicit initial plans |
| Partial GW/FGW | exact and entropic partial-GW; exact partial-FGW; transported-mass semantics |
| Semirelaxed GW/FGW | exact and entropic plans and objectives |
| Translation-invariant KL-UOT | ordinary POT UOT plan/objective plus independent generalized-KL recomputation |
| Sinkhorn divergence | all three component plans, component values, and debiased value under the product-reference KL convention |
| FUGW | sample/feature couplings and the KL-Sinkhorn FUGW objective |
| UCOOT | sample/feature couplings and the independent-regularization UCOOT objective |
| Sampled GW | dense POT reference, `rfugw` feasibility, and full-budget quality versus a six-seed tiny-budget distribution |
| Entropic GW/FGW barycenters | one-step couplings, off-diagonal structure update, and the closed-form feature update derived from POT couplings |

Strict calls fail on warnings, stdout, or stderr unless the exact path is one
of the named exceptions below.  Tolerances are attached to each receipt case;
they are not inferred from a single favorable run.  Plan equality is required
only for deterministic formulations with aligned initialization and stopping
controls.  Stochastic sampled-GW plans are monitors, not parity targets.

## Named POT exceptions

These exceptions are narrow and remain in every comparison receipt:

1. **Entropic partial FGW gradient.** POT 0.9.7.post1 adds a scalar feature
   term to the entropic partial-FGW gradient where the feature-cost matrix is
   required.  POT plan parity is monitored.  The `rfugw` square-loss objective,
   its independent recomputation, feasibility, and gradient-law tests are
   authoritative.
2. **Translation-invariant KL-UOT with unequal marginal penalties.** POT's
   specialized translation-invariant path does not reach the ordinary POT UOT
   objective for the certified unequal-penalty case.  Ordinary POT UOT plus the
   independent generalized-KL objective is authoritative; the specialized plan
   remains a monitor and may not undercut the certified objective.
3. **Sampled GW.** POT's sampled API does not expose inner Sinkhorn convergence
   controls, can emit repeated convergence warnings, cannot reliably exercise
   the full `(n, n)` budget once a row becomes sparse, and is not uniformly
   monotone across seeds.  POT stochastic plans and its `(8, 2)` quality result
   are monitors.  The strict `rfugw` claim is limited to feasibility and a true
   `(n, n)` budget beating the median of six `(2, 1)` runs in objective gap and
   Frobenius distance to the dense reference.
4. **Entropic barycenter return semantics.** POT retains positive entropic
   self-costs on the learned structure diagonal.  In the fused routine it also
   updates `X` but returns the stale initial `Y` variable.  `rfugw`'s metric
   invariant `diag(C) = 0`, off-diagonal parity, coupling parity, and the
   closed-form feature update computed from POT couplings are authoritative.

The current certification makes no POT parity claim for low-rank wrappers that
perform a dense SVD, for alternative structural losses, or for stochastic plan
identity.  Such paths require a new estimand and acceptance rule before they
can become release gates.
