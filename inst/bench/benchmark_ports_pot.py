#!/usr/bin/env python3
"""POT reference benchmarks for newly ported API families."""

from __future__ import annotations

import csv
import statistics
import sys
import time
from pathlib import Path

import numpy as np
import ot


def make_pair_cost(rng: np.random.RandomState, n: int, d: int = 3):
    x1 = rng.normal(size=(n, d))
    x2 = rng.normal(size=(n, d))
    c1 = ot.dist(x1, x1)
    c2 = ot.dist(x2, x2)
    c1 = c1 / c1.max()
    c2 = c2 / c2.max()
    return c1, c2


def make_feature_cost(rng: np.random.RandomState, ns: int, nt: int, d: int = 2):
    f1 = rng.normal(size=(ns, d))
    f2 = rng.normal(size=(nt, d))
    m = ot.dist(f1, f2)
    m = m / m.max()
    return m


def make_ucoot_problem(rng: np.random.RandomState, ns: int, nf: int, nt: int, mf: int):
    x = rng.normal(size=(ns, nf))
    y = rng.normal(size=(nt, mf))
    return x, y


def time_many(fn, runs: int):
    timings = []
    for _ in range(runs):
        t0 = time.perf_counter()
        fn()
        t1 = time.perf_counter()
        timings.append((t1 - t0) * 1000.0)
    return min(timings), statistics.median(timings)


def bench_ports(iters: int, seed: int):
    rng = np.random.RandomState(seed)
    rows: list[dict[str, float | int | str]] = []

    for n in [25, 35]:
        ns = n
        nt = n + 4
        c1_all, c2_all = make_pair_cost(rng, max(ns, nt))
        c1 = c1_all[:ns, :ns]
        c2 = c2_all[:nt, :nt]
        m = make_feature_cost(rng, ns, nt)
        p = ot.unif(ns)
        q = ot.unif(nt)

        def _epgw():
            ot.gromov.entropic_partial_gromov_wasserstein(
                c1,
                c2,
                p=p,
                q=q,
                reg=0.2,
                m=0.7,
                loss_fun="square_loss",
                numItermax=80,
                tol=1e-7,
                log=False,
            )

        mn, med = time_many(_epgw, iters)
        rows.append(
            {
                "suite": "partial",
                "n": n,
                "method": "pot_entropic_partial_gw",
                "min_ms": mn,
                "median_ms": med,
                "iter_per_sec": 1000.0 / med if med > 0 else float("nan"),
                "runs": iters,
            }
        )

        def _epfgw():
            ot.gromov.entropic_partial_fused_gromov_wasserstein(
                m,
                c1,
                c2,
                p=p,
                q=q,
                reg=0.2,
                m=0.7,
                loss_fun="square_loss",
                alpha=0.6,
                numItermax=80,
                tol=1e-7,
                log=False,
            )

        mn, med = time_many(_epfgw, iters)
        rows.append(
            {
                "suite": "partial",
                "n": n,
                "method": "pot_entropic_partial_fgw",
                "min_ms": mn,
                "median_ms": med,
                "iter_per_sec": 1000.0 / med if med > 0 else float("nan"),
                "runs": iters,
            }
        )

        def _pgw():
            ot.gromov.partial_gromov_wasserstein(
                c1,
                c2,
                p=p,
                q=q,
                m=0.7,
                nb_dummies=5,
                loss_fun="square_loss",
                numItermax=120,
                tol=1e-8,
                log=False,
            )

        try:
            mn, med = time_many(_pgw, iters)
            rows.append(
                {
                    "suite": "partial",
                    "n": n,
                    "method": "pot_partial_gw_exact",
                    "min_ms": mn,
                    "median_ms": med,
                    "iter_per_sec": 1000.0 / med if med > 0 else float("nan"),
                    "runs": iters,
                }
            )
        except Exception:
            rows.append(
                {
                    "suite": "partial",
                    "n": n,
                    "method": "pot_partial_gw_exact_error",
                    "min_ms": float("nan"),
                    "median_ms": float("nan"),
                    "iter_per_sec": float("nan"),
                    "runs": iters,
                }
            )

        def _pfgw():
            ot.gromov.partial_fused_gromov_wasserstein(
                m,
                c1,
                c2,
                p=p,
                q=q,
                m=0.7,
                nb_dummies=5,
                loss_fun="square_loss",
                alpha=0.5,
                numItermax=120,
                tol=1e-8,
                log=False,
            )

        try:
            mn, med = time_many(_pfgw, iters)
            rows.append(
                {
                    "suite": "partial",
                    "n": n,
                    "method": "pot_partial_fgw_exact",
                    "min_ms": mn,
                    "median_ms": med,
                    "iter_per_sec": 1000.0 / med if med > 0 else float("nan"),
                    "runs": iters,
                }
            )
        except Exception:
            rows.append(
                {
                    "suite": "partial",
                    "n": n,
                    "method": "pot_partial_fgw_exact_error",
                    "min_ms": float("nan"),
                    "median_ms": float("nan"),
                    "iter_per_sec": float("nan"),
                    "runs": iters,
                }
            )

    for n in [80, 140, 220]:
        ns = n
        nt = round(n * 0.8)
        c1_all, c2_all = make_pair_cost(rng, max(ns, nt))
        c1 = c1_all[:ns, :ns]
        c2 = c2_all[:nt, :nt]
        m = make_feature_cost(rng, ns, nt)
        p = ot.unif(ns)

        def _srgw():
            ot.gromov.semirelaxed_gromov_wasserstein(
                c1,
                c2,
                p=p,
                loss_fun="square_loss",
                max_iter=200,
                tol_rel=1e-9,
                tol_abs=1e-9,
                log=False,
            )

        mn, med = time_many(_srgw, iters)
        rows.append(
            {
                "suite": "semirelaxed",
                "n": n,
                "method": "pot_semirelaxed_gw",
                "min_ms": mn,
                "median_ms": med,
                "iter_per_sec": 1000.0 / med if med > 0 else float("nan"),
                "runs": iters,
            }
        )

        def _srfgw():
            ot.gromov.semirelaxed_fused_gromov_wasserstein(
                m,
                c1,
                c2,
                p=p,
                loss_fun="square_loss",
                alpha=0.55,
                max_iter=200,
                tol_rel=1e-9,
                tol_abs=1e-9,
                log=False,
            )

        mn, med = time_many(_srfgw, iters)
        rows.append(
            {
                "suite": "semirelaxed",
                "n": n,
                "method": "pot_semirelaxed_fgw",
                "min_ms": mn,
                "median_ms": med,
                "iter_per_sec": 1000.0 / med if med > 0 else float("nan"),
                "runs": iters,
            }
        )

    for n in [90, 160, 260]:
        c1, c2 = make_pair_cost(rng, n)
        p = ot.unif(n)
        q = ot.unif(n)

        def _sampled():
            ot.gromov.sampled_gromov_wasserstein(
                c1,
                c2,
                p,
                q,
                loss_fun=lambda a, b: (a - b) ** 2,
                nb_samples_grad=(16, 2),
                epsilon=0.1,
                max_iter=120,
                random_state=seed,
                log=False,
            )

        mn, med = time_many(_sampled, iters)
        rows.append(
            {
                "suite": "sampled",
                "n": n,
                "method": "pot_sampled_gw",
                "min_ms": mn,
                "median_ms": med,
                "iter_per_sec": 1000.0 / med if med > 0 else float("nan"),
                "runs": iters,
            }
        )

    for n in [40, 70, 110]:
        ns = n
        nf = max(10, round(n * 0.6))
        nt = n - 5
        mf = max(9, round(n * 0.5))
        x, y = make_ucoot_problem(rng, ns, nf, nt, mf)

        def _ucoot():
            ot.gromov.unbalanced_co_optimal_transport(
                x,
                y,
                reg_marginals=(10, 8),
                epsilon=(0.05, 0.03),
                divergence="kl",
                unbalanced_solver="sinkhorn",
                alpha=(0.0, 0.0),
                max_iter=40,
                tol=1e-7,
                max_iter_ot=200,
                tol_ot=1e-7,
                log=False,
            )

        mn, med = time_many(_ucoot, iters)
        rows.append(
            {
                "suite": "ucoot",
                "n": n,
                "method": "pot_ucoot_kl",
                "min_ms": mn,
                "median_ms": med,
                "iter_per_sec": 1000.0 / med if med > 0 else float("nan"),
                "runs": iters,
            }
        )

    return rows


def main():
    iters = int(sys.argv[1]) if len(sys.argv) >= 2 else 3
    out_csv = Path(sys.argv[2]) if len(sys.argv) >= 3 else Path("rfugw/inst/bench/results/benchmark_ports_pot_latest.csv")
    seed = int(sys.argv[3]) if len(sys.argv) >= 4 else 123

    rows = bench_ports(iters=iters, seed=seed)
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    with out_csv.open("w", newline="") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=["suite", "n", "method", "min_ms", "median_ms", "iter_per_sec", "runs"],
        )
        writer.writeheader()
        writer.writerows(rows)

    print("== POT benchmark ports ==")
    for r in rows:
        print(
            f"{r['suite']:11s} n={int(r['n']):3d} {r['method']:28s} "
            f"median_ms={r['median_ms']:.3f}"
        )
    print(f"\nWrote POT benchmark CSV: {out_csv}")


if __name__ == "__main__":
    main()
