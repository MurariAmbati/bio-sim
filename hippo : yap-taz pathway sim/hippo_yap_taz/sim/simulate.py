from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Dict, Optional

import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp

from hippo_yap_taz.sim.model import (
    InitialConditions,
    ModelParams,
    ModelSettings,
    STATE_NAMES,
    pack_state,
    rhs,
    unpack_state,
)


@dataclass(frozen=True)
class SimulationSettings:
    stiffness: float
    density: float
    t_end: float
    n_points: int
    method: str


@dataclass(frozen=True)
class SimulationResult:
    settings: SimulationSettings
    params: ModelParams
    initial_conditions: InitialConditions
    t: np.ndarray
    y: np.ndarray
    df: pd.DataFrame
    steady: Dict[str, float]


def _validate_settings(stiffness: float, density: float, t_end: float, n_points: int, method: str) -> None:
    if not (0.0 <= stiffness <= 1.0):
        raise ValueError("stiffness must be in [0, 1]")
    if not (0.0 <= density <= 1.0):
        raise ValueError("density must be in [0, 1]")
    if t_end <= 0:
        raise ValueError("t_end must be > 0")
    if n_points < 5:
        raise ValueError("n_points must be >= 5")
    if method not in {"BDF", "LSODA", "RK45"}:
        raise ValueError("method must be BDF, LSODA, or RK45")


def run_simulation(
    *,
    stiffness: float,
    density: float,
    t_end: float = 300.0,
    n_points: int = 1200,
    method: str = "BDF",
    params: Optional[ModelParams] = None,
    initial_conditions: Optional[InitialConditions] = None,
) -> SimulationResult:
    _validate_settings(stiffness, density, t_end, n_points, method)

    params = params or ModelParams()
    initial_conditions = initial_conditions or InitialConditions()

    settings = ModelSettings(stiffness=float(stiffness), density=float(density))
    sim_settings = SimulationSettings(
        stiffness=float(stiffness),
        density=float(density),
        t_end=float(t_end),
        n_points=int(n_points),
        method=str(method),
    )

    y0 = pack_state(initial_conditions)
    t_eval = np.linspace(0.0, float(t_end), int(n_points))

    sol = solve_ivp(
        fun=lambda t, y: rhs(t, y, settings, params),
        t_span=(0.0, float(t_end)),
        y0=y0,
        t_eval=t_eval,
        method=str(method),
        rtol=1e-6,
        atol=1e-9,
        max_step=float(t_end) / 200.0,
    )

    if not sol.success:
        raise RuntimeError(f"ode solver failed: {sol.message}")

    y = np.asarray(sol.y, dtype=float).T
    df = pd.DataFrame(y, columns=STATE_NAMES)
    df.insert(0, "t", np.asarray(sol.t, dtype=float))

    steady_vec = y[-1, :]
    steady = unpack_state(steady_vec)
    steady.update(
        {
            "yap_total": float(steady["yap_c"] + steady["yap_p"] + steady["yap_n"]),
            "yap_nuclear_fraction": float(
                steady["yap_n"]
                / max(1e-12, steady["yap_c"] + steady["yap_p"] + steady["yap_n"])
            ),
        }
    )

    return SimulationResult(
        settings=sim_settings,
        params=params,
        initial_conditions=initial_conditions,
        t=np.asarray(sol.t, dtype=float),
        y=y,
        df=df,
        steady=steady,
    )


def run_to_steady_state(
    *,
    stiffness: float,
    density: float,
    t_end: float = 800.0,
    n_points: int = 2000,
    method: str = "BDF",
    params: Optional[ModelParams] = None,
    initial_conditions: Optional[InitialConditions] = None,
) -> Dict[str, float]:
    res = run_simulation(
        stiffness=stiffness,
        density=density,
        t_end=t_end,
        n_points=n_points,
        method=method,
        params=params,
        initial_conditions=initial_conditions,
    )
    return dict(res.steady)


def params_as_dict(p: ModelParams) -> Dict[str, float]:
    return dict(asdict(p))
