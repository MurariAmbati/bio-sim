from __future__ import annotations

from dataclasses import dataclass
from typing import Literal


ContextName = Literal["development", "angiogenesis", "tissue_repair"]


@dataclass(frozen=True)
class ModelParams:
    # Ligand-receptor binding
    kon: float = 1.0e-3  # 1/(nM*s)
    koff: float = 2.0e-3  # 1/s

    # Dimerization / activation
    k_dim: float = 5.0e-4  # 1/(nM*s)
    k_undim: float = 1.0e-3  # 1/s
    k_p: float = 8.0e-3  # 1/s  (phosphorylation)
    k_dp: float = 2.5e-3  # 1/s  (dephosphorylation)

    # Adapter recruitment (FRS2/GRB2/SOS lumped)
    k_adapt_on: float = 8.0e-4
    k_adapt_off: float = 2.0e-3

    # MAPK cascade (lumped Ras->Raf->MEK->ERK)
    k_erk_act: float = 9.0e-3
    k_erk_deact: float = 3.0e-3

    # PI3K/AKT branch (lumped)
    k_akt_act: float = 6.0e-3
    k_akt_deact: float = 2.5e-3

    # PLCγ/Ca2+ branch (lumped)
    k_plc_act: float = 4.0e-3
    k_plc_deact: float = 2.0e-3

    # STAT branch (optional lump)
    k_stat_act: float = 2.0e-3
    k_stat_deact: float = 1.5e-3

    # Receptor trafficking (lumped internalization/return)
    k_int: float = 8.0e-4
    k_rec: float = 3.0e-4
    k_deg: float = 1.0e-4

    # Inhibitors (0..1, 1 = full block)
    fgfr_inhib: float = 0.0
    mek_inhib: float = 0.0
    pi3k_inhib: float = 0.0
    plc_inhib: float = 0.0
    stat_inhib: float = 0.0


@dataclass(frozen=True)
class InitialConditions:
    # Species are in nM-equivalent (relative units ok)
    L: float = 10.0
    R: float = 50.0
    LR: float = 0.0
    D: float = 0.0
    Rp: float = 0.0
    A: float = 0.0

    ERKp: float = 0.0
    AKTp: float = 0.0
    PLCa: float = 0.0
    STATp: float = 0.0

    Rint: float = 0.0


@dataclass(frozen=True)
class ContextPreset:
    name: ContextName
    description: str
    params: ModelParams
    ic: InitialConditions


def preset(context: ContextName, ligand_nM: float = 10.0) -> ContextPreset:
    if context == "development":
        return ContextPreset(
            name=context,
            description="Higher mitogenic bias; moderate survival; patterned ERK dynamics.",
            params=ModelParams(k_erk_act=1.1e-2, k_akt_act=5.0e-3, k_int=7.0e-4),
            ic=InitialConditions(L=ligand_nM, R=60.0),
        )
    if context == "angiogenesis":
        return ContextPreset(
            name=context,
            description="Migration/angiogenic bias; stronger PLC/Ca2+ contribution.",
            params=ModelParams(k_plc_act=5.5e-3, k_erk_act=8.0e-3, k_int=9.0e-4),
            ic=InitialConditions(L=ligand_nM, R=45.0),
        )
    if context == "tissue_repair":
        return ContextPreset(
            name=context,
            description="Survival + proliferation; stronger PI3K/AKT tone.",
            params=ModelParams(k_akt_act=8.0e-3, k_erk_act=9.0e-3, k_int=8.0e-4),
            ic=InitialConditions(L=ligand_nM, R=55.0),
        )
    raise ValueError(f"Unknown context: {context}")
