from __future__ import annotations

from dataclasses import asdict, replace
from typing import Dict, Tuple

from ..params import ContextName, InitialConditions, ModelParams, preset


def default_context() -> ContextName:
    return "angiogenesis"


def default_model() -> Tuple[ModelParams, InitialConditions]:
    pr = preset(default_context(), ligand_nM=10.0)
    return pr.params, pr.ic


def apply_knockouts(params: ModelParams, ko_erk: bool, ko_akt: bool, ko_plc: bool, ko_stat: bool) -> ModelParams:
    # Knockouts implemented as strong inhibitors (keeps model stable).
    mek_inhib = 1.0 if ko_erk else params.mek_inhib
    pi3k_inhib = 1.0 if ko_akt else params.pi3k_inhib
    plc_inhib = 1.0 if ko_plc else params.plc_inhib
    stat_inhib = 1.0 if ko_stat else params.stat_inhib
    return replace(
        params,
        mek_inhib=mek_inhib,
        pi3k_inhib=pi3k_inhib,
        plc_inhib=plc_inhib,
        stat_inhib=stat_inhib,
    )
