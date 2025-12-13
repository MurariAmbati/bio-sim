from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, Dict, Tuple

import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp


@dataclass(frozen=True)
class ODEParams:
    # ligand/receptor binding
    kon: float = 1.0
    koff: float = 0.2

    # jak activation/inactivation
    k_jak_on: float = 1.0
    k_jak_off: float = 0.5

    # stat phosphorylation / dephosphorylation
    k_p: float = 2.0
    k_dp: float = 0.6

    # dimerization / dissociation
    k_dim: float = 1.0
    k_undim: float = 0.3

    # nuclear import/export
    k_imp: float = 0.8
    k_exp: float = 0.5

    # socs transcription/decay (negative feedback)
    k_socs_txn: float = 0.8
    k_socs_deg: float = 0.4

    # socs inhibition strength on jak activation
    k_socs_inhib: float = 2.0


STATE_NAMES = (
    "lr",  # ligand-receptor complex
    "jak",  # active jak (0..rt)
    "pstat",  # phosphorylated stat
    "dimer",  # stat dimer
    "nstat",  # nuclear stat (proxy for transcriptional output)
    "socs",  # negative feedback protein
)


def default_initial_conditions() -> Dict[str, float]:
    return {
        "lr": 0.0,
        "jak": 0.0,
        "pstat": 0.0,
        "dimer": 0.0,
        "nstat": 0.0,
        "socs": 0.0,
    }


def cytokine_input(
    t: float,
    *,
    level: float,
    pulse: bool,
    pulse_end: float,
) -> float:
    if not pulse:
        return float(level)
    return float(level) if t <= pulse_end else 0.0


def rhs(
    t: float,
    y: np.ndarray,
    *,
    params: ODEParams,
    cytokine: Callable[[float], float],
    totals: Dict[str, float],
    socs_knockout: bool,
) -> np.ndarray:
    lr, jak, pstat, dimer, nstat, socs = y

    # totals / conserved pools (simple proxies)
    rt = float(totals.get("rt", 1.0))  # receptor total
    stat_t = float(totals.get("stat_t", 1.0))

    ligand = float(cytokine(t))

    # socs inhibition reduces jak activation rate
    socs_effect = 1.0
    if not socs_knockout:
        socs_effect = 1.0 / (1.0 + params.k_socs_inhib * socs)

    # lr complex
    d_lr = params.kon * ligand * (rt - lr) - params.koff * lr

    # jak activation driven by lr, inhibited by socs
    d_jak = params.k_jak_on * lr * socs_effect * (rt - jak) - params.k_jak_off * jak

    # stat phosphorylation depends on jak and available stat
    stat_free = max(stat_t - pstat - 2.0 * dimer - nstat, 0.0)
    d_pstat = params.k_p * jak * stat_free - params.k_dp * pstat - 2.0 * params.k_dim * (pstat**2) + 2.0 * params.k_undim * dimer

    # dimer dynamics (pstat dimerizes)
    d_dimer = params.k_dim * (pstat**2) - params.k_undim * dimer - params.k_imp * dimer

    # nuclear stat from imported dimers, exports/decays (proxy)
    d_nstat = params.k_imp * dimer - params.k_exp * nstat

    # socs transcription induced by nuclear stat
    d_socs = (0.0 if socs_knockout else params.k_socs_txn * nstat) - params.k_socs_deg * socs

    return np.array([d_lr, d_jak, d_pstat, d_dimer, d_nstat, d_socs], dtype=float)


def simulate_ode(
    *,
    params: ODEParams = ODEParams(),
    t_end: float = 50.0,
    n_steps: int = 800,
    cytokine_level: float = 1.0,
    pulse: bool = False,
    pulse_end: float = 10.0,
    socs_knockout: bool = False,
    totals: Dict[str, float] | None = None,
    y0: Dict[str, float] | None = None,
) -> Tuple[pd.DataFrame, Dict[str, float]]:
    totals = dict(totals or {"rt": 1.0, "stat_t": 1.0})
    y0_dict = default_initial_conditions()
    if y0:
        y0_dict.update(y0)

    y0_vec = np.array([y0_dict[name] for name in STATE_NAMES], dtype=float)

    t_eval = np.linspace(0.0, float(t_end), int(n_steps))

    cfun = lambda t: cytokine_input(
        t,
        level=float(cytokine_level),
        pulse=bool(pulse),
        pulse_end=float(pulse_end),
    )

    sol = solve_ivp(
        lambda t, y: rhs(
            t,
            y,
            params=params,
            cytokine=cfun,
            totals=totals,
            socs_knockout=bool(socs_knockout),
        ),
        t_span=(0.0, float(t_end)),
        y0=y0_vec,
        t_eval=t_eval,
        method="RK45",
        rtol=1e-6,
        atol=1e-9,
    )

    if not sol.success:
        raise RuntimeError(f"ode integration failed: {sol.message}")

    df = pd.DataFrame(sol.y.T, columns=list(STATE_NAMES))
    df.insert(0, "t", sol.t)

    meta = {
        "t_end": float(t_end),
        "n_steps": int(n_steps),
        "cytokine_level": float(cytokine_level),
        "pulse": bool(pulse),
        "pulse_end": float(pulse_end),
        "socs_knockout": bool(socs_knockout),
        **{f"total_{k}": float(v) for k, v in totals.items()},
    }

    return df, meta
