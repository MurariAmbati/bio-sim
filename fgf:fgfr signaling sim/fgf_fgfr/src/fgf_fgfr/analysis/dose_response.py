from __future__ import annotations

from dataclasses import replace
from typing import Dict, Iterable, List, Tuple

import numpy as np
import pandas as pd

from ..params import InitialConditions, ModelParams
from ..simulate import outcomes_from_state, simulate


def dose_response(
    doses_nM: Iterable[float],
    base_params: ModelParams,
    base_ic: InitialConditions,
    t_end_s: float = 3600.0,
) -> pd.DataFrame:
    rows: List[Dict[str, float]] = []
    for dose in doses_nM:
        ic = replace(base_ic, L=float(dose))
        state = simulate(base_params, ic, t_end_s=t_end_s)
        out = outcomes_from_state(state)
        row = {"dose_nM": float(dose), **out}
        rows.append(row)
    return pd.DataFrame(rows)
