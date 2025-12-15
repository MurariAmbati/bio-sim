from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp


@dataclass(frozen=True)
class ModelParams:
    antigen: float = 0.8
    cd28: float = 0.5
    lck: float = 1.0
    shp1: float = 0.6
    strength: float = 1.0
    decay: float = 1.0


def _hill(x: float, k: float = 0.5, n: float = 2.0) -> float:
    x = max(0.0, float(x))
    return (x**n) / (k**n + x**n)


def _ode(t: float, y: np.ndarray, p: ModelParams) -> np.ndarray:
    """compact tcr signaling ode.

    state (all in [0, 1] scale, interpreted as activities):
    0: proximal (lck/zap70/lat lump)
    1: plcγ
    2: ip3
    3: ca (cytosolic)
    4: ras
    5: erk
    6: pkcθ
    7: ikk
    8: nfkb
    """

    proximal, plc, ip3, ca, ras, erk, pkc, ikk, nfkb = y

    # inputs
    stim = _hill(p.antigen, k=0.35, n=3.0)
    costim = _hill(p.cd28, k=0.5, n=2.0)

    # negative regulation (shp1) acts more strongly on proximal activation
    inh = 1.0 / (1.0 + 1.5 * max(0.0, p.shp1))

    # global scales
    s = max(0.0, p.strength)
    d = max(0.01, p.decay)

    # proximal activation: antigen + lck, damped by shp1
    k_prox_on = 2.0 * s * max(0.0, p.lck) * stim * (0.6 + 0.4 * costim) * inh
    k_prox_off = 1.2 * d
    d_proximal = k_prox_on * (1.0 - proximal) - k_prox_off * proximal

    # plcγ downstream of proximal
    k_plc_on = 2.2 * s * _hill(proximal, k=0.35, n=2.0)
    k_plc_off = 1.4 * d
    d_plc = k_plc_on * (1.0 - plc) - k_plc_off * plc

    # ip3 production from plcγ
    k_ip3_on = 2.0 * s * _hill(plc, k=0.4, n=2.0)
    k_ip3_off = 1.8 * d
    d_ip3 = k_ip3_on * (1.0 - ip3) - k_ip3_off * ip3

    # calcium release via ip3 + decay/pump
    k_ca_on = 3.0 * s * _hill(ip3, k=0.45, n=3.0)
    k_ca_off = 2.2 * d
    d_ca = k_ca_on * (1.0 - ca) - k_ca_off * ca

    # ras activation via proximal + costim
    k_ras_on = 2.0 * s * _hill(proximal, k=0.35, n=2.0) * (0.7 + 0.3 * costim)
    k_ras_off = 1.6 * d
    d_ras = k_ras_on * (1.0 - ras) - k_ras_off * ras

    # erk downstream of ras
    k_erk_on = 2.4 * s * _hill(ras, k=0.4, n=2.0)
    k_erk_off = 1.7 * d
    d_erk = k_erk_on * (1.0 - erk) - k_erk_off * erk

    # pkcθ activation depends on proximal and calcium (represents dag+ca synergy)
    k_pkc_on = 2.0 * s * _hill(proximal, k=0.35, n=2.0) * (0.4 + 0.6 * _hill(ca, k=0.4, n=2.0))
    k_pkc_off = 1.5 * d
    d_pkc = k_pkc_on * (1.0 - pkc) - k_pkc_off * pkc

    # ikk downstream of pkcθ with some costim contribution
    k_ikk_on = 2.2 * s * _hill(pkc, k=0.4, n=2.0) * (0.7 + 0.3 * costim)
    k_ikk_off = 1.4 * d
    d_ikk = k_ikk_on * (1.0 - ikk) - k_ikk_off * ikk

    # nfkb activation downstream of ikk (nuclear activity proxy)
    k_nf_on = 2.0 * s * _hill(ikk, k=0.4, n=2.0)
    k_nf_off = 1.0 * d
    d_nfkb = k_nf_on * (1.0 - nfkb) - k_nf_off * nfkb

    return np.array([d_proximal, d_plc, d_ip3, d_ca, d_ras, d_erk, d_pkc, d_ikk, d_nfkb], dtype=float)


def simulate(t: np.ndarray, p: ModelParams) -> pd.DataFrame:
    t = np.asarray(t, dtype=float)
    if t.ndim != 1 or t.size < 2:
        raise ValueError("t must be a 1d array with at least 2 points")

    y0 = np.zeros(9, dtype=float)

    sol = solve_ivp(
        fun=lambda tt, yy: _ode(tt, yy, p),
        t_span=(float(t[0]), float(t[-1])),
        y0=y0,
        t_eval=t,
        method="LSODA",
        rtol=1e-6,
        atol=1e-9,
    )
    if not sol.success:
        raise RuntimeError(f"solver failed: {sol.message}")

    cols = [
        "proximal",
        "plcg",
        "ip3",
        "ca",
        "ras",
        "erk",
        "pkctheta",
        "ikk",
        "nfkb",
    ]

    df = pd.DataFrame(sol.y.T, columns=cols)
    df.insert(0, "t_min", sol.t)
    return df
