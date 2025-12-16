from __future__ import annotations

from dataclasses import asdict
from typing import Dict, Optional, Tuple

import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp

from .model import ModelState, SPECIES, outcome_scores, pack_initial_conditions, rhs
from .params import InitialConditions, ModelParams


def simulate(
    params: ModelParams,
    ic: InitialConditions,
    t_end_s: float = 3600.0,
    n_points: int = 600,
    method: str = "LSODA",
    rtol: float = 1e-6,
    atol: float = 1e-9,
) -> ModelState:
    t_eval = np.linspace(0.0, float(t_end_s), int(n_points))
    y0 = pack_initial_conditions(asdict(ic))

    sol = solve_ivp(
        fun=lambda t, y: rhs(t, y, params),
        t_span=(0.0, float(t_end_s)),
        y0=y0,
        t_eval=t_eval,
        method=method,
        rtol=rtol,
        atol=atol,
    )
    if not sol.success:
        raise RuntimeError(f"ODE solve failed: {sol.message}")

    return ModelState(t=sol.t, y=sol.y)


def to_dataframe(state: ModelState) -> pd.DataFrame:
    df = pd.DataFrame({"t_s": state.t})
    for i, name in enumerate(SPECIES):
        df[name] = state.y[i, :]
    return df


def outcomes_from_state(state: ModelState) -> Dict[str, float]:
    y_final = state.y[:, -1]
    return outcome_scores(y_final)
