#!/usr/bin/env python3
"""Benchmark POT sampled GW in sparse benchmark comparison modes."""

from __future__ import annotations

import json
import statistics
import sys
import time

import numpy as np
import ot


def _ru_maxrss_bytes() -> int:
    try:
        import platform
        import resource

        rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        if platform.system().lower() == "darwin":
            return int(rss)
        return int(rss) * 1024
    except Exception:
        return 0


def _time_many(fn, runs: int) -> tuple[float, float]:
    ts = []
    for _ in range(runs):
        t0 = time.perf_counter()
        fn()
        t1 = time.perf_counter()
        ts.append((t1 - t0) * 1000.0)
    return min(ts), statistics.median(ts)


def _run_sampled(c1: np.ndarray, c2: np.ndarray, seed: int, with_log: bool):
    n1 = c1.shape[0]
    n2 = c2.shape[0]
    p = ot.unif(n1)
    q = ot.unif(n2)
    return ot.gromov.sampled_gromov_wasserstein(
        c1,
        c2,
        p,
        q,
        loss_fun=lambda a, b: (a - b) ** 2,
        nb_samples_grad=(16, 2),
        epsilon=0.1,
        max_iter=120,
        random_state=seed,
        log=with_log,
    )


def main() -> int:
    if len(sys.argv) < 4:
        raise SystemExit(
            "Usage: benchmark_sparse_sampled_pot.py <mode: cost|coords> <a_csv> <b_csv> [iters] [seed]"
        )

    mode = str(sys.argv[1]).strip().lower()
    a = np.loadtxt(sys.argv[2], delimiter=",")
    b = np.loadtxt(sys.argv[3], delimiter=",")
    iters = int(sys.argv[4]) if len(sys.argv) >= 5 else 3
    seed = int(sys.argv[5]) if len(sys.argv) >= 6 else 123

    if mode not in {"cost", "coords"}:
        raise SystemExit(f"Unsupported mode: {mode}")

    if mode == "cost":
        c1 = a
        c2 = b

        def _run():
            _run_sampled(c1, c2, seed=seed, with_log=False)

        out_log = _run_sampled(c1, c2, seed=seed, with_log=True)
        method = "pot_sampled_dense_solver_from_cost"
    else:
        e1 = a
        e2 = b

        def _run():
            c1 = ot.dist(e1, e1)
            c2 = ot.dist(e2, e2)
            _run_sampled(c1, c2, seed=seed, with_log=False)

        c1 = ot.dist(e1, e1)
        c2 = ot.dist(e2, e2)
        out_log = _run_sampled(c1, c2, seed=seed, with_log=True)
        method = "pot_sampled_dense_end_to_end_from_coords"

    rss_before = _ru_maxrss_bytes()
    min_ms, median_ms = _time_many(_run, iters)
    rss_after = _ru_maxrss_bytes()
    gw_obj = float(out_log[1]["gw_dist_estimated"])

    out = {
        "method": method,
        "min_ms": float(min_ms),
        "median_ms": float(median_ms),
        "mem_bytes": int(max(0, rss_after - rss_before)),
        "gw_obj": gw_obj,
        "runs": int(iters),
    }
    print(json.dumps(out))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
