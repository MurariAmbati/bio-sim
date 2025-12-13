from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Optional, Tuple

import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp

from .input import Stimulus
from .model import Parameters, SPECIES, pack_state, rhs, unpack_state


@dataclass(frozen=True)
class ODEResult:
    t: np.ndarray
    y: np.ndarray

    def to_frame(self) -> pd.DataFrame:
        df = pd.DataFrame(self.y.T, columns=list(SPECIES))
        df.insert(0, "t", self.t)
        return df


def simulate_ode(
    p: Parameters,
    x0: Dict[str, float],
    stim: Stimulus,
    t_span: Tuple[float, float] = (0.0, 120.0),
    n_points: int = 1201,
    rtol: float = 1e-7,
    atol: float = 1e-9,
) -> ODEResult:
    y0 = pack_state(x0)
    t_eval = np.linspace(float(t_span[0]), float(t_span[1]), int(n_points))

    sol = solve_ivp(
        fun=lambda t, y: rhs(t, y, p=p, stim=stim),
        t_span=(float(t_span[0]), float(t_span[1])),
        y0=y0,
        t_eval=t_eval,
        method="LSODA",
        rtol=rtol,
        atol=atol,
    )

    if not sol.success:
        raise RuntimeError(f"ODE solve failed: {sol.message}")

    y = sol.y.copy()
    # Numerical solvers can introduce tiny negative values near zero; clip only those.
    tiny = 1e-12
    y[(y < 0.0) & (y > -tiny)] = 0.0

    return ODEResult(t=sol.t, y=y)


# -------------------------
# Stochastic (SSA) mode
# -------------------------

@dataclass(frozen=True)
class SSAResult:
    t: np.ndarray
    y: np.ndarray

    def to_frame(self) -> pd.DataFrame:
        df = pd.DataFrame(self.y, columns=list(SPECIES))
        df.insert(0, "t", self.t)
        return df


def _clip_nonneg_int(x: np.ndarray) -> np.ndarray:
    x = np.asarray(x)
    x = np.where(x < 0, 0, x)
    return x.astype(int)


def simulate_ssa(
    p: Parameters,
    x0: Dict[str, float],
    stim: Stimulus,
    t_span: Tuple[float, float] = (0.0, 120.0),
    max_events: int = 200_000,
    volume_fL: float = 1.0,
    seed: Optional[int] = None,
) -> SSAResult:
    """Single-cell stochastic simulation (Gillespie SSA).

    This SSA uses a simplified mass-action reaction network that approximates
    the ODE cascade structure. It is intended for variability exploration.

    The ODE concentrations are converted to molecule counts via volume.
    """

    rng = np.random.default_rng(seed)

    # concentration (µM) -> molecules: N = conc_uM * 1e-6 (mol/L) * NA * V(L)
    NA = 6.02214076e23
    V_L = float(volume_fL) * 1e-15
    scale = NA * V_L * 1e-6

    species = list(SPECIES)
    idx = {s: i for i, s in enumerate(species)}

    y = np.zeros(len(species), dtype=int)
    for s in species:
        y[idx[s]] = int(max(0.0, x0.get(s, 0.0) * scale))

    t0, t1 = float(t_span[0]), float(t_span[1])
    t = t0

    ts = [t]
    ys = [y.copy()]

    def conc(name: str) -> float:
        return y[idx[name]] / scale

    # For SSA we define propensities with mass-action approximations.
    # Reactions (stoichiometry is in molecule counts):
    # 1) R + L -> R*  (L(t) acts as time-varying factor)
    # 2) R* -> R
    # 3) R* -> degraded (R*)
    # 4) RasGDP -> RasGTP (catalyzed by R*, inhibited by ERKPP)
    # 5) RasGTP -> RasGDP
    # 6) Raf -> Raf* (catalyzed by RasGTP)
    # 7) Raf* -> Raf
    # 8-9) MEK -> MEKP -> MEKPP (catalyzed by Raf*)
    # 10-11) MEKPP -> MEKP -> MEK (basal)
    # 12-13) ERK -> ERKP -> ERKPP (catalyzed by MEKPP)
    # 14-15) ERKPP -> ERKP -> ERK (basal + DUSP)
    # 16) DUSP synthesis (Hill in ERKPP)
    # 17) DUSP degradation

    def hill(x: float, K: float, n: float) -> float:
        x = max(x, 0.0)
        K = max(K, 1e-12)
        n = max(n, 1e-12)
        return (x**n) / (K**n + x**n)

    for _ in range(int(max_events)):
        if t >= t1:
            break

        L = stim.value(t)

        R = y[idx["R"]]
        R_star = y[idx["R_star"]]
        RasGDP = y[idx["RasGDP"]]
        RasGTP = y[idx["RasGTP"]]
        Raf = y[idx["Raf"]]
        Raf_star = y[idx["Raf_star"]]
        MEK = y[idx["MEK"]]
        MEK_P = y[idx["MEK_P"]]
        MEK_PP = y[idx["MEK_PP"]]
        ERK = y[idx["ERK"]]
        ERK_P = y[idx["ERK_P"]]
        ERK_PP = y[idx["ERK_PP"]]
        DUSP = y[idx["DUSP"]]

        # feedback gate uses concentrations for smoothness
        fb = 1.0 / (1.0 + p.fb_Ras_alpha * hill(conc("ERK_PP"), p.fb_Ras_K, p.fb_Ras_n))

        a = np.zeros(17, dtype=float)

        a[0] = p.k_R_on * L * R
        a[1] = p.k_R_off * R_star
        a[2] = p.k_R_deg * R_star

        # catalytic props: scale by catalysts count and substrate count
        a[3] = fb * p.V_Ras_act * R_star * (RasGDP / max(p.K_Ras_act * scale, 1.0))
        a[4] = p.k_Ras_gdp * RasGTP

        a[5] = p.k_Raf_act * RasGTP * (Raf / max(p.K_Raf_act * scale, 1.0))
        a[6] = p.k_Raf_deact * Raf_star

        a[7] = p.V_MEK_p1 * Raf_star * (MEK / max(p.K_MEK_p1 * scale, 1.0))
        a[8] = p.V_MEK_p2 * Raf_star * (MEK_P / max(p.K_MEK_p2 * scale, 1.0))
        a[9] = p.V_MEK_dp2 * (MEK_PP / max(p.K_MEK_dp2 * scale, 1.0))
        a[10] = p.V_MEK_dp1 * (MEK_P / max(p.K_MEK_dp1 * scale, 1.0))

        a[11] = p.V_ERK_p1 * MEK_PP * (ERK / max(p.K_ERK_p1 * scale, 1.0))
        a[12] = p.V_ERK_p2 * MEK_PP * (ERK_P / max(p.K_ERK_p2 * scale, 1.0))

        Vdp1 = (p.V_ERK_dp1_basal + p.V_ERK_dp1_dusp * (DUSP / scale))
        Vdp2 = (p.V_ERK_dp2_basal + p.V_ERK_dp2_dusp * (DUSP / scale))
        a[13] = Vdp2 * (ERK_PP / max(p.K_ERK_dp2 * scale, 1.0))
        a[14] = Vdp1 * (ERK_P / max(p.K_ERK_dp1 * scale, 1.0))

        a[15] = (p.dusp_syn0 + p.dusp_syn_max * hill(conc("ERK_PP"), p.dusp_K, p.dusp_n)) * scale
        a[16] = p.dusp_deg * DUSP

        a0 = float(np.sum(a))
        if not np.isfinite(a0) or a0 <= 0.0:
            # nothing happens; jump to end
            t = t1
            ts.append(t)
            ys.append(y.copy())
            break

        r1, r2 = rng.random(2)
        dt = -np.log(max(r1, 1e-15)) / a0
        t_next = t + dt

        # record at event time
        t = t_next

        # choose reaction
        threshold = r2 * a0
        cum = 0.0
        mu = 0
        for i in range(len(a)):
            cum += a[i]
            if cum >= threshold:
                mu = i
                break

        # apply stoichiometry
        if mu == 0 and R > 0:
            y[idx["R"]] -= 1
            y[idx["R_star"]] += 1
        elif mu == 1 and R_star > 0:
            y[idx["R_star"]] -= 1
            y[idx["R"]] += 1
        elif mu == 2 and R_star > 0:
            y[idx["R_star"]] -= 1
        elif mu == 3 and RasGDP > 0:
            y[idx["RasGDP"]] -= 1
            y[idx["RasGTP"]] += 1
        elif mu == 4 and RasGTP > 0:
            y[idx["RasGTP"]] -= 1
            y[idx["RasGDP"]] += 1
        elif mu == 5 and Raf > 0:
            y[idx["Raf"]] -= 1
            y[idx["Raf_star"]] += 1
        elif mu == 6 and Raf_star > 0:
            y[idx["Raf_star"]] -= 1
            y[idx["Raf"]] += 1
        elif mu == 7 and MEK > 0:
            y[idx["MEK"]] -= 1
            y[idx["MEK_P"]] += 1
        elif mu == 8 and MEK_P > 0:
            y[idx["MEK_P"]] -= 1
            y[idx["MEK_PP"]] += 1
        elif mu == 9 and MEK_PP > 0:
            y[idx["MEK_PP"]] -= 1
            y[idx["MEK_P"]] += 1
        elif mu == 10 and MEK_P > 0:
            y[idx["MEK_P"]] -= 1
            y[idx["MEK"]] += 1
        elif mu == 11 and ERK > 0:
            y[idx["ERK"]] -= 1
            y[idx["ERK_P"]] += 1
        elif mu == 12 and ERK_P > 0:
            y[idx["ERK_P"]] -= 1
            y[idx["ERK_PP"]] += 1
        elif mu == 13 and ERK_PP > 0:
            y[idx["ERK_PP"]] -= 1
            y[idx["ERK_P"]] += 1
        elif mu == 14 and ERK_P > 0:
            y[idx["ERK_P"]] -= 1
            y[idx["ERK"]] += 1
        elif mu == 15:
            y[idx["DUSP"]] += 1
        elif mu == 16 and DUSP > 0:
            y[idx["DUSP"]] -= 1

        y = _clip_nonneg_int(y)
        ts.append(t)
        ys.append(y.copy())

    t_arr = np.asarray(ts, dtype=float)
    y_arr_counts = np.asarray(ys, dtype=int)

    # counts -> concentrations for reporting (match ODE columns)
    y_arr = y_arr_counts / scale

    return SSAResult(t=t_arr, y=y_arr)
