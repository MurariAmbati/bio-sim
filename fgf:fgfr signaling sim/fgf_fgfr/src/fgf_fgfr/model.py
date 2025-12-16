from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Tuple

import numpy as np

from .params import ModelParams


SPECIES: List[str] = [
    "L",
    "R",
    "LR",
    "D",
    "Rp",
    "A",
    "ERKp",
    "AKTp",
    "PLCa",
    "STATp",
    "Rint",
]


@dataclass(frozen=True)
class ModelState:
    t: np.ndarray
    y: np.ndarray  # shape (n_species, n_time)

    def as_dict(self) -> Dict[str, np.ndarray]:
        return {name: self.y[i, :] for i, name in enumerate(SPECIES)}


def rhs(t: float, y: np.ndarray, p: ModelParams) -> np.ndarray:
    # y order: SPECIES
    L, R, LR, D, Rp, A, ERKp, AKTp, PLCa, STATp, Rint = y

    # Inhibitory scalars
    fgfr_scale = 1.0 - float(np.clip(p.fgfr_inhib, 0.0, 1.0))
    mek_scale = 1.0 - float(np.clip(p.mek_inhib, 0.0, 1.0))
    pi3k_scale = 1.0 - float(np.clip(p.pi3k_inhib, 0.0, 1.0))
    plc_scale = 1.0 - float(np.clip(p.plc_inhib, 0.0, 1.0))
    stat_scale = 1.0 - float(np.clip(p.stat_inhib, 0.0, 1.0))

    # Ligand binding
    bind = p.kon * L * R
    unbind = p.koff * LR

    # Dimerization from LR
    dim = p.k_dim * LR * LR
    undim = p.k_undim * D

    # Activation/phosphorylation of dimer
    phos = fgfr_scale * p.k_p * D
    dephos = p.k_dp * Rp

    # Adapter recruitment depends on active receptor Rp
    adapt_on = p.k_adapt_on * Rp
    adapt_off = p.k_adapt_off * A

    # Downstream activation: simple first-order driven by A
    erk_act = mek_scale * p.k_erk_act * A
    erk_deact = p.k_erk_deact * ERKp

    akt_act = pi3k_scale * p.k_akt_act * A
    akt_deact = p.k_akt_deact * AKTp

    plc_act = plc_scale * p.k_plc_act * A
    plc_deact = p.k_plc_deact * PLCa

    stat_act = stat_scale * p.k_stat_act * A
    stat_deact = p.k_stat_deact * STATp

    # Trafficking: internalize active + inactive receptor (lumped from R + LR + D + Rp)
    # This is intentionally simple; it mainly provides a desensitization timescale.
    rec_pool = R + LR + 2.0 * D + Rp
    internalize = p.k_int * rec_pool
    recycle = p.k_rec * Rint
    degrade = p.k_deg * Rint

    dL = -bind + unbind
    dR = -bind + unbind - internalize + recycle
    dLR = bind - unbind - 2.0 * dim + 2.0 * undim
    dD = dim - undim - phos
    dRp = phos - dephos
    dA = adapt_on - adapt_off

    dERKp = erk_act - erk_deact
    dAKTp = akt_act - akt_deact
    dPLCa = plc_act - plc_deact
    dSTATp = stat_act - stat_deact

    dRint = internalize - recycle - degrade

    return np.array([dL, dR, dLR, dD, dRp, dA, dERKp, dAKTp, dPLCa, dSTATp, dRint], dtype=float)


def outcome_scores(y: np.ndarray) -> Dict[str, float]:
    # Heuristic "phenotype" scores from terminal signaling state.
    # Uses final values and a bit of cross-talk.
    L, R, LR, D, Rp, A, ERKp, AKTp, PLCa, STATp, Rint = y

    proliferation = 0.65 * ERKp + 0.25 * STATp + 0.10 * AKTp
    survival = 0.75 * AKTp + 0.15 * ERKp + 0.10 * STATp
    migration = 0.55 * PLCa + 0.35 * ERKp + 0.10 * AKTp
    angiogenesis = 0.45 * ERKp + 0.35 * PLCa + 0.20 * AKTp

    return {
        "proliferation": float(max(0.0, proliferation)),
        "survival": float(max(0.0, survival)),
        "migration": float(max(0.0, migration)),
        "angiogenesis": float(max(0.0, angiogenesis)),
    }


def species_index(name: str) -> int:
    try:
        return SPECIES.index(name)
    except ValueError as exc:
        raise KeyError(name) from exc


def pack_initial_conditions(ic: Dict[str, float]) -> np.ndarray:
    y0 = np.zeros(len(SPECIES), dtype=float)
    for k, v in ic.items():
        y0[species_index(k)] = float(v)
    return y0


def unpack(y: np.ndarray) -> Dict[str, float]:
    return {SPECIES[i]: float(y[i]) for i in range(len(SPECIES))}
