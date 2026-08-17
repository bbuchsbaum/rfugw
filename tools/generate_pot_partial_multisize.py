#!/usr/bin/env python3
"""Compact POT objectives for partial GW/FGW at several sizes and masses."""
import json
import pathlib
import numpy as np
import ot

ROOT = pathlib.Path(__file__).resolve().parents[1]
out_path = ROOT / "inst" / "extdata" / "fixtures" / "partial_pot_multisize.json"

rng = np.random.RandomState(20260816)
cases = []
for n, m, reg in [(6, 0.5, 0.25), (6, 0.8, 0.25), (8, 0.6, 0.2), (10, 0.7, 0.2)]:
    Xs = rng.normal(size=(n, 3))
    Xt = rng.normal(size=(n, 3))
    Fs = rng.normal(size=(n, 2))
    Ft = rng.normal(size=(n, 2))
    C1 = ot.dist(Xs, Xs)
    C2 = ot.dist(Xt, Xt)
    C1 = C1 / C1.max()
    C2 = C2 / C2.max()
    M = ot.dist(Fs, Ft)
    M = M / M.max()
    p = ot.unif(n)
    q = ot.unif(n)
    T_gw, log_gw = ot.gromov.entropic_partial_gromov_wasserstein(
        C1,
        C2,
        p=p,
        q=q,
        reg=reg,
        m=m,
        loss_fun="square_loss",
        numItermax=60,
        tol=1e-7,
        log=True,
    )
    T_fgw, log_fgw = ot.gromov.entropic_partial_fused_gromov_wasserstein(
        M,
        C1,
        C2,
        p=p,
        q=q,
        reg=reg,
        m=m,
        loss_fun="square_loss",
        alpha=0.5,
        numItermax=60,
        tol=1e-7,
        log=True,
    )
    cases.append(
        {
            "n": n,
            "m": m,
            "reg": reg,
            "alpha": 0.5,
            "numItermax": 60,
            "tol": 1e-7,
            "C1": C1.tolist(),
            "C2": C2.tolist(),
            "M": M.tolist(),
            "p": p.tolist(),
            "q": q.tolist(),
            "epgw_dist": float(
                log_gw.get("partial_gw_dist", log_gw.get("gw_dist", np.nan))
            ),
            "epfgw_dist": float(
                log_fgw.get("partial_fgw_dist", log_fgw.get("fgw_dist", np.nan))
            ),
        }
    )

out_path.write_text(
    json.dumps(
        {
            "pot_version": ot.__version__,
            "seed": 20260816,
            "cases": cases,
        },
        indent=2,
    )
)
print("wrote", out_path, "n_cases", len(cases))
