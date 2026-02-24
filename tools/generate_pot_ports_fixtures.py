#!/usr/bin/env python3
import json
import pathlib
import numpy as np
import ot

ROOT = pathlib.Path(__file__).resolve().parents[1]
out_path = ROOT / "inst" / "extdata" / "fixtures" / "pot_api_ports_fixture.json"
out_path.parent.mkdir(parents=True, exist_ok=True)

rng = np.random.RandomState(123)

def tolist(x):
    if isinstance(x, np.ndarray):
        return x.tolist()
    return x

fixture = {}

# Non-entropic semirelaxed GW/FGW fixtures
ns, nt = 10, 8
Xs = rng.normal(size=(ns, 3))
Xt = rng.normal(size=(nt, 3))
Fs = rng.normal(size=(ns, 2))
Ft = rng.normal(size=(nt, 2))

C1 = ot.dist(Xs, Xs)
C2 = ot.dist(Xt, Xt)
C1 = C1 / C1.max()
C2 = C2 / C2.max()
M = ot.dist(Fs, Ft)
M = M / M.max()
p = ot.unif(ns)

T_srgw, log_srgw = ot.gromov.semirelaxed_gromov_wasserstein(
    C1,
    C2,
    p=p,
    loss_fun="square_loss",
    max_iter=200,
    tol_rel=1e-9,
    tol_abs=1e-9,
    log=True,
)

T_srfgw, log_srfgw = ot.gromov.semirelaxed_fused_gromov_wasserstein(
    M,
    C1,
    C2,
    p=p,
    loss_fun="square_loss",
    alpha=0.55,
    max_iter=200,
    tol_rel=1e-9,
    tol_abs=1e-9,
    log=True,
)

fixture["semirelaxed"] = {
    "inputs": {
        "C1": tolist(C1),
        "C2": tolist(C2),
        "M": tolist(M),
        "p": tolist(p),
    },
    "params": {
        "alpha": 0.55,
        "max_iter": 200,
        "tol_rel": 1e-9,
        "tol_abs": 1e-9,
    },
    "outputs": {
        "T_srgw": tolist(T_srgw),
        "srgw_dist": float(log_srgw["srgw_dist"]),
        "T_srfgw": tolist(T_srfgw),
        "srfgw_dist": float(log_srfgw["srfgw_dist"]),
    },
}

# Entropic partial GW/FGW fixtures
ns, nt = 9, 11
Xs = rng.normal(size=(ns, 3))
Xt = rng.normal(size=(nt, 3))
Fs = rng.normal(size=(ns, 2))
Ft = rng.normal(size=(nt, 2))

C1 = ot.dist(Xs, Xs)
C2 = ot.dist(Xt, Xt)
C1 = C1 / C1.max()
C2 = C2 / C2.max()
M = ot.dist(Fs, Ft)
M = M / M.max()
p = ot.unif(ns)
q = ot.unif(nt)
mass = 0.7

T_epgw, log_epgw = ot.gromov.entropic_partial_gromov_wasserstein(
    C1,
    C2,
    p=p,
    q=q,
    reg=0.2,
    m=mass,
    loss_fun="square_loss",
    numItermax=80,
    tol=1e-7,
    log=True,
)

T_epfgw, log_epfgw = ot.gromov.entropic_partial_fused_gromov_wasserstein(
    M,
    C1,
    C2,
    p=p,
    q=q,
    reg=0.2,
    m=mass,
    loss_fun="square_loss",
    alpha=0.6,
    numItermax=80,
    tol=1e-7,
    log=True,
)

fixture["entropic_partial"] = {
    "inputs": {
        "C1": tolist(C1),
        "C2": tolist(C2),
        "M": tolist(M),
        "p": tolist(p),
        "q": tolist(q),
    },
    "params": {
        "reg": 0.2,
        "m": mass,
        "alpha": 0.6,
        "numItermax": 80,
        "tol": 1e-7,
    },
    "outputs": {
        "T_epgw": tolist(T_epgw),
        "epgw_dist": float(log_epgw["partial_gw_dist"]),
        "T_epfgw": tolist(T_epfgw),
        "epfgw_dist": float(log_epfgw["partial_fgw_dist"]),
    },
}

# UCOOT fixture
nx_samp, nx_feat = 7, 5
ny_samp, ny_feat = 6, 4
X = rng.normal(size=(nx_samp, nx_feat))
Y = rng.normal(size=(ny_samp, ny_feat))

pi_samp, pi_feat, log_u = ot.gromov.unbalanced_co_optimal_transport(
    X,
    Y,
    reg_marginals=(10, 8),
    epsilon=(0.05, 0.03),
    divergence="kl",
    unbalanced_solver="sinkhorn",
    alpha=(0.0, 0.0),
    max_iter=40,
    tol=1e-7,
    max_iter_ot=200,
    tol_ot=1e-7,
    log=True,
)

fixture["ucoot"] = {
    "inputs": {
        "X": tolist(X),
        "Y": tolist(Y),
    },
    "params": {
        "reg_marginals": [10.0, 8.0],
        "epsilon": [0.05, 0.03],
        "max_iter": 40,
        "tol": 1e-7,
        "max_iter_ot": 200,
        "tol_ot": 1e-7,
    },
    "outputs": {
        "pi_samp": tolist(pi_samp),
        "pi_feat": tolist(pi_feat),
        "ucoot_cost": float(log_u["ucoot_cost"]),
    },
}

# Sampled GW fixture (deterministic with fixed random_state)
n = 16
X1 = rng.normal(size=(n, 3))
X2 = rng.normal(size=(n, 3))
C1 = ot.dist(X1, X1)
C2 = ot.dist(X2, X2)
C1 = C1 / C1.max()
C2 = C2 / C2.max()
p = ot.unif(n)
q = ot.unif(n)

T_samp, log_samp = ot.gromov.sampled_gromov_wasserstein(
    C1,
    C2,
    p,
    q,
    loss_fun=lambda a, b: (a - b) ** 2,
    nb_samples_grad=(8, 2),
    epsilon=0.1,
    max_iter=80,
    random_state=77,
    log=True,
)

fixture["sampled"] = {
    "inputs": {
        "C1": tolist(C1),
        "C2": tolist(C2),
        "p": tolist(p),
        "q": tolist(q),
    },
    "params": {
        "nb_samples_grad": [8, 2],
        "epsilon": 0.1,
        "max_iter": 80,
        "random_state": 77,
    },
    "outputs": {
        "plan": tolist(T_samp),
        "gw_dist_estimated": float(log_samp["gw_dist_estimated"]),
    },
}

out_path.write_text(json.dumps(fixture))
print(f"wrote {out_path}")
