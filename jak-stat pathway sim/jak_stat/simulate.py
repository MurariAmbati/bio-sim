from __future__ import annotations

from dataclasses import asdict

from .model_ode import ODEParams, simulate_ode
from .model_logical import LogicalParams, simulate_logical


def demo():
    df_ode, meta_ode = simulate_ode(params=ODEParams(), t_end=40.0, cytokine_level=1.0)
    df_log, meta_log = simulate_logical(params=LogicalParams(), steps=25, cytokine_on=True)

    return {
        "ode": {"meta": meta_ode, "tail": df_ode.tail(3).to_dict(orient="records"), "params": asdict(ODEParams())},
        "logical": {"meta": meta_log, "tail": df_log.tail(3).to_dict(orient="records")},
    }
