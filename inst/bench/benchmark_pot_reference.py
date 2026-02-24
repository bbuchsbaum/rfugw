#!/usr/bin/env python3
"""Reference POT benchmarks on the same synthetic workloads as rfugw bench suite."""

from __future__ import annotations

import csv
import statistics
import sys
import time
from pathlib import Path

import numpy as np
import ot


def make_fgw_problem(rng: np.random.RandomState, n: int, d_struct: int = 3, d_feat: int = 5):
    x1 = rng.normal(size=(n, d_struct))
    x2 = rng.normal(size=(n, d_struct))
    y1 = rng.normal(size=(n, d_feat))
    y2 = rng.normal(size=(n, d_feat))

    c1 = ot.dist(x1, x1)
    c2 = ot.dist(x2, x2)
    c1 = c1 / c1.max()
    c2 = c2 / c2.max()

    m = ot.dist(y1, y2)
    m = m / m.max()
    p = ot.unif(n)
    q = ot.unif(n)
    return c1, c2, m, p, q


def make_fugw_problem(rng: np.random.RandomState, n: int):
    x1 = rng.normal(size=(n, 3))
    x2 = rng.normal(size=(n, 3))
    cx = ot.dist(x1, x1)
    cy = ot.dist(x2, x2)
    cx = cx / cx.max()
    cy = cy / cy.max()
    m = np.ones((n, n), dtype=np.float64)
    idx = np.arange(n)
    m[idx, idx[::-1]] = 0.0
    wx = ot.unif(n)
    wy = ot.unif(n)
    return cx, cy, m, wx, wy


def time_many(fn, runs: int):
    timings = []
    for _ in range(runs):
        t0 = time.perf_counter()
        fn()
        t1 = time.perf_counter()
        timings.append((t1 - t0) * 1000.0)
    return min(timings), statistics.median(timings)


def bench_pot(iters: int, seed: int):
    rng = np.random.RandomState(seed)
    rows = []

    for n in [40, 80, 120, 160]:
        c1, c2, m, p, q = make_fgw_problem(rng, n)

        def _run_fgw():
            ot.gromov.entropic_fused_gromov_wasserstein(
                m,
                c1,
                c2,
                p,
                q,
                loss_fun="square_loss",
                epsilon=0.05,
                alpha=0.5,
                max_iter=200,
                tol=1e-9,
                solver="PGD",
                numItermax=500,
                stopThr=1e-9,
                log=False,
            )

        mn, med = time_many(_run_fgw, iters)
        rows.append(
            {
                "suite": "fgw",
                "n": n,
                "method": "pot_fgw_pgd",
                "min_ms": mn,
                "median_ms": med,
                "iter_per_sec": 1000.0 / med if med > 0 else float("nan"),
                "runs": iters,
            }
        )

        if n <= 80:
            def _run_fgw_exact():
                ot.gromov.fused_gromov_wasserstein(
                    m,
                    c1,
                    c2,
                    p,
                    q,
                    loss_fun="square_loss",
                    alpha=0.5,
                    armijo=False,
                    max_iter=120,
                    tol_rel=1e-9,
                    tol_abs=1e-9,
                    log=False,
                )

            mn, med = time_many(_run_fgw_exact, iters)
            rows.append(
                {
                    "suite": "fgw",
                    "n": n,
                    "method": "pot_fgw_exact",
                    "min_ms": mn,
                    "median_ms": med,
                    "iter_per_sec": 1000.0 / med if med > 0 else float("nan"),
                    "runs": iters,
                }
            )

    for n in [30, 50, 70, 90]:
        cx, cy, m, wx, wy = make_fugw_problem(rng, n)

        def _run_fugw():
            ot.gromov.fused_unbalanced_gromov_wasserstein(
                cx,
                cy,
                wx=wx,
                wy=wy,
                reg_marginals=(100, 50),
                epsilon=1e-2,
                divergence="kl",
                unbalanced_solver="sinkhorn",
                alpha=0.5,
                M=m,
                max_iter=80,
                tol=1e-8,
                max_iter_ot=300,
                tol_ot=1e-8,
                log=False,
            )

        mn, med = time_many(_run_fugw, iters)
        rows.append(
            {
                "suite": "fugw",
                "n": n,
                "method": "pot_fugw_kl",
                "min_ms": mn,
                "median_ms": med,
                "iter_per_sec": 1000.0 / med if med > 0 else float("nan"),
                "runs": iters,
            }
        )

    return rows


def main():
    iters = int(sys.argv[1]) if len(sys.argv) >= 2 else 3
    out_csv = Path(sys.argv[2]) if len(sys.argv) >= 3 else Path("rfugw/inst/bench/results/pot_benchmark_latest.csv")
    seed = int(sys.argv[3]) if len(sys.argv) >= 4 else 42

    rows = bench_pot(iters=iters, seed=seed)
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    with out_csv.open("w", newline="") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=["suite", "n", "method", "min_ms", "median_ms", "iter_per_sec", "runs"],
        )
        writer.writeheader()
        writer.writerows(rows)

    print("== POT reference benchmark ==")
    for r in rows:
        print(
            f"{r['suite']:4s} n={int(r['n']):3d} {r['method']:12s} "
            f"median_ms={r['median_ms']:.3f} iter_per_sec={r['iter_per_sec']:.3f}"
        )
    print(f"\nWrote POT benchmark CSV: {out_csv}")


if __name__ == "__main__":
    main()
