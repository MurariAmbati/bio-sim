from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Tuple

import numpy as np

from .input import Stimulus


SPECIES: Tuple[str, ...] = (
    # receptor + Ras
    "R",
    "R_star",
    "RasGDP",
    "RasGTP",
    # Raf
    "Raf",
    "Raf_star",
    # MEK phosphorylation states
    "MEK",
    "MEK_P",
    "MEK_PP",
    # ERK phosphorylation states
    "ERK",
    "ERK_P",
    "ERK_PP",
    # inducible phosphatase (negative feedback)
    "DUSP",
)


@dataclass
class Parameters:
    # receptor activation
    k_R_on: float = 1.0
    k_R_off: float = 0.2
    k_R_deg: float = 0.01

    # Ras activation/inactivation
    V_Ras_act: float = 1.5
    K_Ras_act: float = 0.5
    k_Ras_gdp: float = 1.0  # GTP hydrolysis (RasGTP -> RasGDP)

    # ERK->Ras negative feedback (reduces Ras activation)
    fb_Ras_alpha: float = 2.0
    fb_Ras_K: float = 0.4
    fb_Ras_n: float = 2.0

    # Raf activation/inactivation
    k_Raf_act: float = 2.0
    K_Raf_act: float = 0.3
    k_Raf_deact: float = 0.4

    # MEK phosphorylation by Raf*
    V_MEK_p1: float = 3.0
    K_MEK_p1: float = 0.5
    V_MEK_p2: float = 3.0
    K_MEK_p2: float = 0.5

    # MEK dephosphorylation (basal)
    V_MEK_dp1: float = 1.5
    K_MEK_dp1: float = 0.5
    V_MEK_dp2: float = 1.5
    K_MEK_dp2: float = 0.5

    # ERK phosphorylation by MEK_PP
    V_ERK_p1: float = 5.0
    K_ERK_p1: float = 0.4
    V_ERK_p2: float = 5.0
    K_ERK_p2: float = 0.4

    # ERK dephosphorylation: basal + DUSP-dependent
    V_ERK_dp1_basal: float = 1.0
    V_ERK_dp2_basal: float = 1.0
    V_ERK_dp1_dusp: float = 3.0
    V_ERK_dp2_dusp: float = 3.0
    K_ERK_dp1: float = 0.3
    K_ERK_dp2: float = 0.3

    # DUSP induction and decay
    dusp_syn0: float = 0.02
    dusp_syn_max: float = 0.4
    dusp_K: float = 0.25
    dusp_n: float = 2.0
    dusp_deg: float = 0.05

    # totals / pools
    Ras_tot: float = 1.0
    Raf_tot: float = 1.0
    MEK_tot: float = 1.0
    ERK_tot: float = 1.0


def default_parameters() -> Parameters:
    return Parameters()


def default_initial_conditions(p: Parameters) -> Dict[str, float]:
    return {
        "R": 1.0,
        "R_star": 0.0,
        "RasGDP": p.Ras_tot,
        "RasGTP": 0.0,
        "Raf": p.Raf_tot,
        "Raf_star": 0.0,
        "MEK": p.MEK_tot,
        "MEK_P": 0.0,
        "MEK_PP": 0.0,
        "ERK": p.ERK_tot,
        "ERK_P": 0.0,
        "ERK_PP": 0.0,
        "DUSP": 0.05,
    }


def _hill(x: float, K: float, n: float) -> float:
    x = max(x, 0.0)
    K = max(K, 1e-12)
    n = max(n, 1e-12)
    return (x**n) / (K**n + x**n)


def _mm(V: float, S: float, K: float) -> float:
    S = max(S, 0.0)
    K = max(K, 1e-12)
    return V * S / (K + S)


def pack_state(x0: Dict[str, float]) -> np.ndarray:
    return np.array([float(x0[s]) for s in SPECIES], dtype=float)


def unpack_state(y: np.ndarray) -> Dict[str, float]:
    return {s: float(y[i]) for i, s in enumerate(SPECIES)}


def rhs(t: float, y: np.ndarray, p: Parameters, stim: Stimulus) -> np.ndarray:
    """ODE right-hand side.

    State variables are concentrations (arbitrary consistent units).
    """

    s = unpack_state(y)

    R = s["R"]
    R_star = s["R_star"]
    RasGDP = s["RasGDP"]
    RasGTP = s["RasGTP"]
    Raf = s["Raf"]
    Raf_star = s["Raf_star"]
    MEK = s["MEK"]
    MEK_P = s["MEK_P"]
    MEK_PP = s["MEK_PP"]
    ERK = s["ERK"]
    ERK_P = s["ERK_P"]
    ERK_PP = s["ERK_PP"]
    DUSP = s["DUSP"]

    L = stim.value(t)

    # receptor
    v_R_on = p.k_R_on * L * R
    v_R_off = p.k_R_off * R_star
    v_R_deg = p.k_R_deg * R_star

    # ERK -> Ras negative feedback gate (reduces Ras activation)
    fb = 1.0 / (1.0 + p.fb_Ras_alpha * _hill(ERK_PP, p.fb_Ras_K, p.fb_Ras_n))

    # Ras GDP->GTP activated by receptor
    v_Ras_act = fb * R_star * _mm(p.V_Ras_act, RasGDP, p.K_Ras_act)
    v_Ras_gdp = p.k_Ras_gdp * RasGTP

    # Raf activation by RasGTP
    v_Raf_act = RasGTP * _mm(p.k_Raf_act, Raf, p.K_Raf_act)
    v_Raf_deact = p.k_Raf_deact * Raf_star

    # MEK dual phosphorylation by Raf*
    v_MEK_p1 = Raf_star * _mm(p.V_MEK_p1, MEK, p.K_MEK_p1)
    v_MEK_p2 = Raf_star * _mm(p.V_MEK_p2, MEK_P, p.K_MEK_p2)

    # MEK dephosphorylation
    v_MEK_dp1 = _mm(p.V_MEK_dp1, MEK_P, p.K_MEK_dp1)
    v_MEK_dp2 = _mm(p.V_MEK_dp2, MEK_PP, p.K_MEK_dp2)

    # ERK dual phosphorylation by MEK_PP
    v_ERK_p1 = MEK_PP * _mm(p.V_ERK_p1, ERK, p.K_ERK_p1)
    v_ERK_p2 = MEK_PP * _mm(p.V_ERK_p2, ERK_P, p.K_ERK_p2)

    # ERK dephosphorylation (basal + DUSP-dependent)
    V_dp1 = p.V_ERK_dp1_basal + p.V_ERK_dp1_dusp * DUSP
    V_dp2 = p.V_ERK_dp2_basal + p.V_ERK_dp2_dusp * DUSP
    v_ERK_dp1 = _mm(V_dp1, ERK_P, p.K_ERK_dp1)
    v_ERK_dp2 = _mm(V_dp2, ERK_PP, p.K_ERK_dp2)

    # DUSP dynamics induced by ppERK
    v_dusp_syn = p.dusp_syn0 + p.dusp_syn_max * _hill(ERK_PP, p.dusp_K, p.dusp_n)
    v_dusp_deg = p.dusp_deg * DUSP

    dy = {k: 0.0 for k in SPECIES}

    dy["R"] += -v_R_on + v_R_off
    dy["R_star"] += v_R_on - v_R_off - v_R_deg

    dy["RasGDP"] += -v_Ras_act + v_Ras_gdp
    dy["RasGTP"] += v_Ras_act - v_Ras_gdp

    dy["Raf"] += -v_Raf_act + v_Raf_deact
    dy["Raf_star"] += v_Raf_act - v_Raf_deact

    dy["MEK"] += -v_MEK_p1 + v_MEK_dp1
    dy["MEK_P"] += v_MEK_p1 - v_MEK_p2 - v_MEK_dp1 + v_MEK_dp2
    dy["MEK_PP"] += v_MEK_p2 - v_MEK_dp2

    dy["ERK"] += -v_ERK_p1 + v_ERK_dp1
    dy["ERK_P"] += v_ERK_p1 - v_ERK_p2 - v_ERK_dp1 + v_ERK_dp2
    dy["ERK_PP"] += v_ERK_p2 - v_ERK_dp2

    dy["DUSP"] += v_dusp_syn - v_dusp_deg

    return np.array([dy[s] for s in SPECIES], dtype=float)


def validate_conservation(y: np.ndarray, p: Parameters) -> Dict[str, float]:
    s = unpack_state(y)
    return {
        "Ras_tot": s["RasGDP"] + s["RasGTP"],
        "Raf_tot": s["Raf"] + s["Raf_star"],
        "MEK_tot": s["MEK"] + s["MEK_P"] + s["MEK_PP"],
        "ERK_tot": s["ERK"] + s["ERK_P"] + s["ERK_PP"],
        "R_tot": s["R"] + s["R_star"],
    }


def species_index() -> Dict[str, int]:
    return {s: i for i, s in enumerate(SPECIES)}


def species_list() -> List[str]:
    return list(SPECIES)
