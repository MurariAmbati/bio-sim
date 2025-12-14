from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd

from .core import CoreMapping
from .integrators import EulerIntegrator, IntegratorConfig, OdeSystem, RK4Integrator
from .math_utils import clamp01, hill
from .params import OdeParams


@dataclass(frozen=True)
class NetworkState:
    """State container for the extended network."""

    ptch: float
    smo: float
    gli_a: float
    gli_r: float
    hhip: float

    def as_vector(self) -> np.ndarray:
        return np.array([self.ptch, self.smo, self.gli_a, self.gli_r, self.hhip], dtype=float)

    @staticmethod
    def from_vector(y: np.ndarray) -> "NetworkState":
        return NetworkState(ptch=float(y[0]), smo=float(y[1]), gli_a=float(y[2]), gli_r=float(y[3]), hhip=float(y[4]))


class ExtendedNetworkSystem(OdeSystem):
    """Vector field for the extended HH → PTCH ⟂ SMO → GLI(A/R) network with HHIP sink."""

    def __init__(
        self,
        hh: float,
        params: OdeParams,
        ptch_ko: bool = False,
        smo_ko: bool = False,
        gli_ko: bool = False,
    ) -> None:
        self.hh = clamp01(hh)
        self.p = params
        self.ptch_ko = bool(ptch_ko)
        self.smo_ko = bool(smo_ko)
        self.gli_ko = bool(gli_ko)
        self.mapping = CoreMapping(params.core)

    def _hh_eff(self, hhip: float) -> float:
        return float(np.clip(self.hh / (1.0 + max(0.0, self.p.hhip_strength) * max(0.0, hhip)), 0.0, 1.0))

    def f(self, t: float, y: np.ndarray) -> np.ndarray:
        # unpack
        P, S, A, R, H = [float(v) for v in y]

        hh_eff = self._hh_eff(H)

        fb = float(hill(A, n=self.p.fb_n, k=self.p.fb_k))
        p_prod = self.p.k_p_prod * (1.0 + self.p.fb_strength * fb)
        dP = p_prod - self.p.k_p_deg * P

        h_act = float(hill(A, n=self.p.hhip_n, k=self.p.hhip_k))
        dH = self.p.k_h_prod * h_act - self.p.k_h_deg * H

        P_eff = 0.0 if self.ptch_ko else P
        smo_target, _gli_ss = self.mapping.transduction(hh_eff, P_eff)
        smo_target = float(np.asarray(smo_target))

        dS = self.p.k_s_act * (smo_target - S) - self.p.k_s_deg * S

        a_target = float(hill(S, n=self.p.core.hill_n, k=self.p.core.hill_k))
        r_target = float(1.0 - a_target)
        dA = self.p.k_a_act * (a_target - A) - self.p.k_a_deg * A
        dR = self.p.k_r_prod * (r_target - R) - self.p.k_r_deg * R

        if self.ptch_ko:
            dP = 0.0
        if self.smo_ko:
            dS = 0.0
        if self.gli_ko:
            dA = 0.0
            dR = 0.0

        return np.array([dP, dS, dA, dR, dH], dtype=float)


class NetworkSimulator:
    """High-level simulator built on reusable integrators."""

    def __init__(self, params: OdeParams) -> None:
        self.params = params

    def simulate_extended(
        self,
        hh: float,
        t_end: float,
        dt: float,
        method: str = "rk4",
        ptch_ko: bool = False,
        smo_ko: bool = False,
        gli_ko: bool = False,
        y0: NetworkState | None = None,
    ) -> pd.DataFrame:
        y0 = y0 or NetworkState(ptch=0.35, smo=0.05, gli_a=0.05, gli_r=0.25, hhip=0.10)
        system = ExtendedNetworkSystem(hh=hh, params=self.params, ptch_ko=ptch_ko, smo_ko=smo_ko, gli_ko=gli_ko)
        cfg = IntegratorConfig(dt=float(dt), t_end=float(t_end))

        method = method.strip().lower()
        if method == "euler":
            integrator = EulerIntegrator(cfg)
        else:
            integrator = RK4Integrator(cfg)

        t, y = integrator.solve(system, y0.as_vector())

        P = np.clip(y[:, 0], 0.0, 3.0)
        S = np.clip(y[:, 1], 0.0, 1.0)
        A = np.clip(y[:, 2], 0.0, 1.0)
        R = np.clip(y[:, 3], 0.0, 1.0)
        H = np.clip(y[:, 4], 0.0, 3.0)

        gli_net = np.clip(A - self.params.r_weight * R, 0.0, 1.0)
        return pd.DataFrame({"t": t, "ptch": P, "smo": S, "gliA": A, "gliR": R, "gli_net": gli_net, "hhip": H})
