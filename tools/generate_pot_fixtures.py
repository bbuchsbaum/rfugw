#!/usr/bin/env python3
import json
import numpy as np
import pathlib

ROOT = pathlib.Path(__file__).resolve().parents[2]

import ot  # noqa: E402

rng = np.random.RandomState(7)
out_dir = ROOT / "rfugw" / "inst" / "extdata" / "fixtures"
out_dir.mkdir(parents=True, exist_ok=True)


def tolist(x):
    if isinstance(x, np.ndarray):
        return x.tolist()
    return x

# FGW entropic fixture
ns, nt = 9, 11
Xs = rng.normal(size=(ns, 3))
Xt = rng.normal(size=(nt, 3))
Ys = rng.normal(size=(ns, 2))
Yt = rng.normal(size=(nt, 2))

C1 = ot.dist(Xs, Xs)
C2 = ot.dist(Xt, Xt)
M = ot.dist(Ys, Yt)

C1 = C1 / C1.max()
C2 = C2 / C2.max()
M = M / M.max()

p = ot.unif(ns)
q = ot.unif(nt)

fgw_params = {
    "alpha": 0.4,
    "epsilon": 0.05,
    "max_iter": 300,
    "tol": 1e-10,
    "solver": "PGD",
    "sinkhorn_numItermax": 5000,
    "sinkhorn_stopThr": 1e-12,
}

T_fgw, log_fgw = ot.gromov.entropic_fused_gromov_wasserstein(
    M,
    C1,
    C2,
    p,
    q,
    loss_fun="square_loss",
    epsilon=fgw_params["epsilon"],
    alpha=fgw_params["alpha"],
    max_iter=fgw_params["max_iter"],
    tol=fgw_params["tol"],
    solver=fgw_params["solver"],
    numItermax=fgw_params["sinkhorn_numItermax"],
    stopThr=fgw_params["sinkhorn_stopThr"],
    log=True,
)

fgw_fixture = {
    "inputs": {
        "M": tolist(M),
        "C1": tolist(C1),
        "C2": tolist(C2),
        "p": tolist(p),
        "q": tolist(q),
    },
    "params": fgw_params,
    "outputs": {
        "plan": tolist(T_fgw),
        "fgw_dist": float(log_fgw["fgw_dist"]),
    },
}

(out_dir / "fgw_entropic_square_fixture.json").write_text(json.dumps(fgw_fixture))

# FGW exact-CG fixture (square loss)
n_exact = 8
Xe = rng.normal(size=(n_exact, 3))
Ye = rng.normal(size=(n_exact, 2))
perm = np.arange(n_exact)[::-1]

C1e = ot.dist(Xe, Xe)
C2e = C1e[perm][:, perm]
Me = ot.dist(Ye, Ye[perm])
C1e = C1e / C1e.max()
C2e = C2e / C2e.max()
Me = Me / Me.max()
pe = ot.unif(n_exact)
qe = ot.unif(n_exact)

fgw_exact_params = {
    "alpha": 0.5,
    "max_iter": 120,
    "tol_rel": 1e-9,
    "tol_abs": 1e-9,
    "armijo": False,
}

T_exact, log_exact = ot.gromov.fused_gromov_wasserstein(
    Me,
    C1e,
    C2e,
    pe,
    qe,
    loss_fun="square_loss",
    alpha=fgw_exact_params["alpha"],
    armijo=fgw_exact_params["armijo"],
    max_iter=fgw_exact_params["max_iter"],
    tol_rel=fgw_exact_params["tol_rel"],
    tol_abs=fgw_exact_params["tol_abs"],
    log=True,
)

fgw_exact_fixture = {
    "inputs": {
        "M": tolist(Me),
        "C1": tolist(C1e),
        "C2": tolist(C2e),
        "p": tolist(pe),
        "q": tolist(qe),
    },
    "params": fgw_exact_params,
    "outputs": {
        "plan": tolist(T_exact),
        "fgw_dist": float(log_exact["fgw_dist"]),
    },
}

(out_dir / "fgw_exact_square_fixture.json").write_text(json.dumps(fgw_exact_fixture))

# FUGW fixture (KL + Sinkhorn)
nx, ny = 7, 8
Xs = rng.normal(size=(nx, 3))
Xt = rng.normal(size=(ny, 3))
Cx = ot.dist(Xs, Xs)
Cy = ot.dist(Xt, Xt)
Cx = Cx / Cx.max()
Cy = Cy / Cy.max()

wx = ot.unif(nx)
wy = ot.unif(ny)
M_lin = np.ones((nx, ny))
for i in range(min(nx, ny)):
    M_lin[i, ny - 1 - i] = 0

fugw_params = {
    "reg_marginals": [100.0, 50.0],
    "epsilon": 1e-2,
    "alpha": 0.5,
    "max_iter": 60,
    "tol": 1e-9,
    "max_iter_ot": 500,
    "tol_ot": 1e-9,
    "unbalanced_solver": "sinkhorn",
    "divergence": "kl",
}

pi_s, pi_f, log_fugw = ot.gromov.fused_unbalanced_gromov_wasserstein(
    Cx,
    Cy,
    wx=wx,
    wy=wy,
    reg_marginals=tuple(fugw_params["reg_marginals"]),
    epsilon=fugw_params["epsilon"],
    divergence=fugw_params["divergence"],
    unbalanced_solver=fugw_params["unbalanced_solver"],
    alpha=fugw_params["alpha"],
    M=M_lin,
    init_pi=None,
    max_iter=fugw_params["max_iter"],
    tol=fugw_params["tol"],
    max_iter_ot=fugw_params["max_iter_ot"],
    tol_ot=fugw_params["tol_ot"],
    log=True,
)

fugw2 = ot.gromov.fused_unbalanced_gromov_wasserstein2(
    Cx,
    Cy,
    wx=wx,
    wy=wy,
    reg_marginals=tuple(fugw_params["reg_marginals"]),
    epsilon=fugw_params["epsilon"],
    divergence=fugw_params["divergence"],
    unbalanced_solver=fugw_params["unbalanced_solver"],
    alpha=fugw_params["alpha"],
    M=M_lin,
    init_pi=None,
    max_iter=fugw_params["max_iter"],
    tol=fugw_params["tol"],
    max_iter_ot=fugw_params["max_iter_ot"],
    tol_ot=fugw_params["tol_ot"],
    log=False,
)

fugw_fixture = {
    "inputs": {
        "Cx": tolist(Cx),
        "Cy": tolist(Cy),
        "wx": tolist(wx),
        "wy": tolist(wy),
        "M": tolist(M_lin),
    },
    "params": fugw_params,
    "outputs": {
        "pi_samp": tolist(pi_s),
        "pi_feat": tolist(pi_f),
        "fugw_cost": float(fugw2),
        "fugw_cost_log": float(log_fugw["fugw_cost"]),
    },
}

(out_dir / "fugw_kl_sinkhorn_fixture.json").write_text(json.dumps(fugw_fixture))

print("wrote fixtures to", out_dir)
