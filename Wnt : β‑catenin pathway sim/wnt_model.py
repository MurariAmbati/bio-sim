from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.integrate import solve_ivp


@dataclass(frozen=True)
class WntParams:
    k_prod: float = 0.25
    k_stab: float = 0.6
    k_base: float = 0.05
    k_deg: float = 0.9

    k_syn: float = 0.8
    k_turn: float = 0.25

    k_inhib: float = 1.2
    K: float = 0.7
    n: float = 3.0


def wnt_ode(t: float, y: np.ndarray, p: WntParams, wnt: float) -> np.ndarray:
    """compact wnt / β‑catenin ode.

    state:
      y[0] = beta  (effective cytosolic/total β‑catenin)
      y[1] = dc    (active destruction complex)

    wnt in [0, 1] increases beta stability and decreases dc synthesis.
    beta nonlinearly inhibits dc activity (a minimal feedback that can yield multistability).
    """

    beta = max(float(y[0]), 0.0)
    dc = max(float(y[1]), 0.0)
    wnt = float(np.clip(wnt, 0.0, 1.0))

    hill = (beta**p.n) / (p.K**p.n + beta**p.n + 1e-12)

    d_beta = p.k_prod + p.k_stab * wnt - p.k_deg * dc * beta - p.k_base * beta
    d_dc = p.k_syn * (1.0 - wnt) - p.k_turn * dc - p.k_inhib * hill * dc

    return np.array([d_beta, d_dc], dtype=float)


def simulate(
    p: WntParams,
    wnt: float,
    y0: tuple[float, float],
    t_end: float,
    n_points: int = 1200,
    method: str = "LSODA",
) -> tuple[np.ndarray, np.ndarray]:
    t_eval = np.linspace(0.0, float(t_end), int(n_points))

    sol = solve_ivp(
        fun=lambda t, y: wnt_ode(t, y, p, wnt),
        t_span=(0.0, float(t_end)),
        y0=np.array(y0, dtype=float),
        t_eval=t_eval,
        method=method,
        rtol=1e-7,
        atol=1e-9,
    )
    if not sol.success:
        raise RuntimeError(sol.message)

    return sol.t, sol.y


def endpoints_from_multiple_inits(
    p: WntParams,
    wnt: float,
    t_end: float,
    beta_inits: np.ndarray,
    dc_inits: np.ndarray,
) -> np.ndarray:
    endpoints = []
    for b0 in beta_inits:
        for d0 in dc_inits:
            _, y = simulate(p=p, wnt=wnt, y0=(float(b0), float(d0)), t_end=t_end, n_points=600)
            endpoints.append([float(y[0, -1]), float(y[1, -1])])
    return np.array(endpoints, dtype=float)


def cluster_endpoints(endpoints: np.ndarray, tol: float = 2e-2) -> np.ndarray:
    """greedy clustering by euclidean distance; returns unique centers."""
    if endpoints.size == 0:
        return endpoints

    centers: list[np.ndarray] = []
    for pt in endpoints:
        if not centers:
            centers.append(pt)
            continue
        dists = np.linalg.norm(np.stack(centers) - pt, axis=1)
        if float(np.min(dists)) > tol:
            centers.append(pt)

    return np.stack(centers)
