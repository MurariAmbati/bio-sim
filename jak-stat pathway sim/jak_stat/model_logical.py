from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Tuple

import pandas as pd


NODES = (
    "cytokine",
    "lr",
    "jak",
    "pstat",
    "dimer",
    "nstat",
    "socs",
)


@dataclass(frozen=True)
class LogicalParams:
    # if True, socs node is forced off and does not inhibit
    socs_knockout: bool = False


def _bool(x: int | bool) -> int:
    return 1 if bool(x) else 0


def step(state: Dict[str, int], *, params: LogicalParams) -> Dict[str, int]:
    cytokine = _bool(state.get("cytokine", 0))
    lr = _bool(state.get("lr", 0))
    jak = _bool(state.get("jak", 0))
    pstat = _bool(state.get("pstat", 0))
    dimer = _bool(state.get("dimer", 0))
    nstat = _bool(state.get("nstat", 0))
    socs = _bool(state.get("socs", 0))

    # wiring: cytokine -> lr -> jak -> pstat -> dimer -> nstat -> socs
    # feedback: socs inhibits jak activation

    next_lr = cytokine

    if params.socs_knockout:
        next_socs = 0
        next_jak = lr
    else:
        next_socs = nstat
        next_jak = 1 if (lr and not socs) else 0

    next_pstat = jak
    next_dimer = pstat
    next_nstat = dimer

    return {
        "cytokine": cytokine,
        "lr": next_lr,
        "jak": next_jak,
        "pstat": next_pstat,
        "dimer": next_dimer,
        "nstat": next_nstat,
        "socs": next_socs,
    }


def simulate_logical(
    *,
    params: LogicalParams = LogicalParams(),
    steps: int = 30,
    cytokine_on: bool = True,
    pulse: bool = False,
    pulse_end_step: int = 10,
    initial: Dict[str, int] | None = None,
) -> Tuple[pd.DataFrame, Dict[str, int | bool]]:
    state = {name: 0 for name in NODES}
    if initial:
        for k, v in initial.items():
            if k in state:
                state[k] = _bool(v)

    rows: List[Dict[str, int]] = []
    for k in range(int(steps) + 1):
        if not pulse:
            state["cytokine"] = _bool(cytokine_on)
        else:
            state["cytokine"] = 1 if (k <= int(pulse_end_step) and bool(cytokine_on)) else 0

        rows.append({"step": k, **{n: int(state[n]) for n in NODES}})
        state = step(state, params=params)

    df = pd.DataFrame(rows)
    meta = {
        "steps": int(steps),
        "cytokine_on": bool(cytokine_on),
        "pulse": bool(pulse),
        "pulse_end_step": int(pulse_end_step),
        "socs_knockout": bool(params.socs_knockout),
    }
    return df, meta
