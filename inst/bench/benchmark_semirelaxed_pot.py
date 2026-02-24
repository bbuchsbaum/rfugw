#!/usr/bin/env python3
"""Reference POT benchmarks for semirelaxed GW/FGW workloads."""

from __future__ import annotations

import csv
import statistics
import sys
import time
from pathlib import Path

import numpy as np
import ot


def make_problem(rng: np.random.RandomState, n: int, d_struct: int = 3, d_feat: int = 5):
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
    return c1, c2, m, p


def time_many(fn, runs: int):
    vals = []
    for _ in range(runs):
        t0 = time.perf_counter()
        fn()
        t1 = time.perf_counter()
        vals.append((t1 - t0) * 1000.0)
    return min(vals), statistics.median(vals)


def bench(iters: int, seed: int):
    rng = np.random.RandomState(seed)
    rows = []

    for n in [80, 120, 180, 240]:
        c1, c2, m, p = make_problem(rng, n)

        warm_srfgw = ot.gromov.entropic_semirelaxed_fused_gromov_wasserstein(
            m,
            c1,
            c2,
            p=p,
            epsilon=0.05,
            alpha=0.5,
            symmetric=True,
            max_iter=200,
            tol=1e-9,
            log=False,
        )
        warm_srfgw_obj = ot.gromov.semirelaxed_fused_gromov_wasserstein2(
            m,
            c1,
            c2,
            p=p,
            alpha=0.5,
            T=warm_srfgw,
        )

        def _run_srfgw():
            ot.gromov.entropic_semirelaxed_fused_gromov_wasserstein(
                m,
                c1,
                c2,
                p=p,
                epsilon=0.05,
                alpha=0.5,
                symmetric=True,
                max_iter=200,
                tol=1e-9,
                log=False,
            )

        mn, med = time_many(_run_srfgw, iters)
        rows.append(
            {
                "suite": "srfgw",
                "n": n,
                "method": "pot_srfgw",
                "min_ms": mn,
                "median_ms": med,
                "iter_per_sec": 1000.0 / med if med > 0 else float("nan"),
                "warm_obj": warm_srfgw_obj,
                "runs": iters,
            }
        )

        warm_srgw = ot.gromov.entropic_semirelaxed_gromov_wasserstein(
            c1,
            c2,
            p=p,
            epsilon=0.05,
            symmetric=True,
            max_iter=200,
            tol=1e-9,
            log=False,
        )
        warm_srgw_obj = ot.gromov.semirelaxed_gromov_wasserstein2(
            c1,
            c2,
            p=p,
            T=warm_srgw,
        )

        def _run_srgw():
            ot.gromov.entropic_semirelaxed_gromov_wasserstein(
                c1,
                c2,
                p=p,
                epsilon=0.05,
                symmetric=True,
                max_iter=200,
                tol=1e-9,
                log=False,
            )

        mn, med = time_many(_run_srgw, iters)
        rows.append(
            {
                "suite": "srgw",
                "n": n,
                "method": "pot_srgw",
                "min_ms": mn,
                "median_ms": med,
                "iter_per_sec": 1000.0 / med if med > 0 else float("nan"),
                "warm_obj": warm_srgw_obj,
                "runs": iters,
            }
        )

    return rows


def main():
    iters = int(sys.argv[1]) if len(sys.argv) >= 2 else 3
    out_csv = Path(sys.argv[2]) if len(sys.argv) >= 3 else Path("rfugw/inst/bench/results/pot_benchmark_semirelaxed_latest.csv")
    seed = int(sys.argv[3]) if len(sys.argv) >= 4 else 42

    rows = bench(iters=iters, seed=seed)
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    with out_csv.open("w", newline="") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=["suite", "n", "method", "min_ms", "median_ms", "iter_per_sec", "warm_obj", "runs"],
        )
        writer.writeheader()
        writer.writerows(rows)

    print("== POT semirelaxed benchmark ==")
    for r in rows:
        print(
            f"{r['suite']:5s} n={int(r['n']):3d} {r['method']:10s} "
            f"median_ms={r['median_ms']:.3f} iter_per_sec={r['iter_per_sec']:.3f} "
            f"warm_obj={r['warm_obj']:.8f}"
        )
    print(f"\nWrote POT benchmark CSV: {out_csv}")


if __name__ == "__main__":
    main()
