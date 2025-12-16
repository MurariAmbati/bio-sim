from __future__ import annotations

from dataclasses import replace
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd

from ..params import InitialConditions, ModelParams
from ..simulate import outcomes_from_state, simulate


PARAM_FIELDS = [
    "kon",
    "koff",
    "k_dim",
    "k_p",
    "k_dp",
    "k_adapt_on",
    "k_adapt_off",
    "k_erk_act",
    "k_akt_act",
    "k_plc_act",
    "k_int",
]


def local_sensitivity(
    base_params: ModelParams,
    base_ic: InitialConditions,
    outcome: str = "proliferation",
    rel_step: float = 0.1,
    t_end_s: float = 3600.0,
) -> pd.DataFrame:
    # Finite difference on log-parameters: d log(outcome) / d log(param)
    base_state = simulate(base_params, base_ic, t_end_s=t_end_s)
    base_out = outcomes_from_state(base_state)
    y0 = float(base_out[outcome])
    if y0 <= 0:
        y0 = 1e-12

    rows: List[Dict[str, float]] = []
    for field in PARAM_FIELDS:
        p0 = float(getattr(base_params, field))
        if p0 <= 0:
            continue
        p1 = p0 * (1.0 + rel_step)
        p2 = p0 * (1.0 - rel_step)

        p_plus = replace(base_params, **{field: p1})
        p_minus = replace(base_params, **{field: p2})

        out_plus = float(outcomes_from_state(simulate(p_plus, base_ic, t_end_s=t_end_s))[outcome])
        out_minus = float(outcomes_from_state(simulate(p_minus, base_ic, t_end_s=t_end_s))[outcome])

        out_plus = max(out_plus, 1e-12)
        out_minus = max(out_minus, 1e-12)

        # symmetric log-derivative
        sens = (np.log(out_plus) - np.log(out_minus)) / (np.log(p1) - np.log(p2))
        rows.append({"param": field, "elasticity": float(sens)})

    df = pd.DataFrame(rows).sort_values("elasticity", key=lambda s: s.abs(), ascending=False)
    return df
