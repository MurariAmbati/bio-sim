from __future__ import annotations

import numpy as np
import pandas as pd

from .core import CoreMapping
from .math_utils import clamp01, hill
from .params import OdeParams


def hh_effective(hh: float, hhip: float, strength: float) -> float:
    hh = clamp01(hh)
    hhip = max(0.0, float(hhip))
    strength = max(0.0, float(strength))
    return float(np.clip(hh / (1.0 + strength * hhip), 0.0, 1.0))


class FeedbackDynamicsModel:
    """ODE-level dynamics for hedgehog pathway modules."""

    def __init__(self, params: OdeParams) -> None:
        self.p = params
        self.mapping = CoreMapping(params.core)

    def simulate_simple(self, hh: float, t_end: float, dt: float, ptch_ko: bool, smo_ko: bool, gli_ko: bool) -> pd.DataFrame:
        hh = clamp01(hh)
        dt = max(1e-4, float(dt))
        t_end = max(dt, float(t_end))

        steps = int(np.ceil(t_end / dt))
        t = np.linspace(0.0, steps * dt, steps + 1)

        P = np.zeros_like(t)
        S = np.zeros_like(t)
        G = np.zeros_like(t)
        P[0] = 0.35
        S[0] = 0.05
        G[0] = 0.05

        for i in range(steps):
            Pi, Si, Gi = P[i], S[i], G[i]

            fb = float(hill(Gi, n=self.p.fb_n, k=self.p.fb_k))
            p_prod = self.p.k_p_prod * (1.0 + self.p.fb_strength * fb)
            dP = p_prod - self.p.k_p_deg * Pi

            P_eff = 0.0 if ptch_ko else Pi
            smo_target, gli_target = self.mapping.transduction(hh, P_eff)
            smo_target = float(np.asarray(smo_target))
            gli_target = float(np.asarray(gli_target))

            dS = self.p.k_s_act * (smo_target - Si) - self.p.k_s_deg * Si
            dG = self.p.k_g_act * (gli_target - Gi) - self.p.k_g_deg * Gi

            Pn = Pi + dt * dP
            Sn = Si + dt * dS
            Gn = Gi + dt * dG

            if ptch_ko:
                Pn = 0.0
            if smo_ko:
                Sn = 0.0
            if gli_ko:
                Gn = 0.0

            P[i + 1] = float(np.clip(Pn, 0.0, 3.0))
            S[i + 1] = float(np.clip(Sn, 0.0, 1.0))
            G[i + 1] = float(np.clip(Gn, 0.0, 1.0))

        return pd.DataFrame({"t": t, "ptch": P, "smo": S, "gli": G})

    def simulate_extended(self, hh: float, t_end: float, dt: float, ptch_ko: bool, smo_ko: bool, gli_ko: bool) -> pd.DataFrame:
        hh = clamp01(hh)
        dt = max(1e-4, float(dt))
        t_end = max(dt, float(t_end))

        steps = int(np.ceil(t_end / dt))
        t = np.linspace(0.0, steps * dt, steps + 1)

        P = np.zeros_like(t)
        S = np.zeros_like(t)
        A = np.zeros_like(t)
        R = np.zeros_like(t)
        H = np.zeros_like(t)

        P[0] = 0.35
        S[0] = 0.05
        A[0] = 0.05
        R[0] = 0.25
        H[0] = 0.10

        for i in range(steps):
            Pi, Si, Ai, Ri, Hi = P[i], S[i], A[i], R[i], H[i]

            hh_eff = hh_effective(hh, Hi, self.p.hhip_strength)

            fb = float(hill(Ai, n=self.p.fb_n, k=self.p.fb_k))
            p_prod = self.p.k_p_prod * (1.0 + self.p.fb_strength * fb)
            dP = p_prod - self.p.k_p_deg * Pi

            h_act = float(hill(Ai, n=self.p.hhip_n, k=self.p.hhip_k))
            dH = self.p.k_h_prod * h_act - self.p.k_h_deg * Hi

            P_eff = 0.0 if ptch_ko else Pi
            smo_target, _gli_ss = self.mapping.transduction(hh_eff, P_eff)
            smo_target = float(np.asarray(smo_target))

            dS = self.p.k_s_act * (smo_target - Si) - self.p.k_s_deg * Si

            a_target = float(hill(Si, n=self.p.core.hill_n, k=self.p.core.hill_k))
            r_target = float(1.0 - a_target)
            dA = self.p.k_a_act * (a_target - Ai) - self.p.k_a_deg * Ai
            dR = self.p.k_r_prod * (r_target - Ri) - self.p.k_r_deg * Ri

            Pn = Pi + dt * dP
            Sn = Si + dt * dS
            An = Ai + dt * dA
            Rn = Ri + dt * dR
            Hn = Hi + dt * dH

            if ptch_ko:
                Pn = 0.0
            if smo_ko:
                Sn = 0.0
            if gli_ko:
                An = 0.0
                Rn = 0.0

            P[i + 1] = float(np.clip(Pn, 0.0, 3.0))
            S[i + 1] = float(np.clip(Sn, 0.0, 1.0))
            A[i + 1] = float(np.clip(An, 0.0, 1.0))
            R[i + 1] = float(np.clip(Rn, 0.0, 1.0))
            H[i + 1] = float(np.clip(Hn, 0.0, 3.0))

        gli_net = np.clip(A - self.p.r_weight * R, 0.0, 1.0)
        return pd.DataFrame({"t": t, "ptch": P, "smo": S, "gliA": A, "gliR": R, "gli_net": gli_net, "hhip": H})
