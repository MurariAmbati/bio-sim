from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Dict, List

import numpy as np


STATE_NAMES: List[str] = [
    "mst",  # mst1/2 activity (upstream hippo kinases)
    "lats",  # lats1/2 activity
    "yap_c",  # unphosphorylated cytoplasmic yap/taz
    "yap_p",  # phosphorylated cytoplasmic yap/taz
    "yap_n",  # nuclear yap/taz (active)
    "gene",  # representative yap/tead target gene expression
    "actin",  # actin tension / rhoa-like mechanosensing state
]


@dataclass(frozen=True)
class ModelSettings:
    stiffness: float = 0.7  # 0..1
    density: float = 0.2  # 0..1 (contact inhibition)


@dataclass(frozen=True)
class ModelParams:
    # mst module
    k_m_act: float = 0.35
    k_m_ci: float = 0.50
    k_m_deact: float = 0.55

    # lats module
    k_l_act: float = 1.10
    k_l_deact: float = 0.80

    # yap shuttling / phosphorylation
    k_y_syn: float = 0.25
    k_y_deg_c: float = 0.02
    k_y_deg_p: float = 0.03
    k_y_deg_n: float = 0.02

    k_phos: float = 2.20
    k_dephos: float = 0.65

    k_imp: float = 0.85
    k_exp: float = 0.65

    # transcriptional response (yap/tead)
    k_txn: float = 1.00
    k_gene_deg: float = 0.15
    hill_k: float = 0.35
    hill_n: float = 2.2

    # actin mechanosensing state
    k_a_act: float = 1.20
    k_a_inhib_mst: float = 0.85
    k_a_decay: float = 0.35

    # coupling: actin inhibits mst/lats activation and promotes nuclear import
    actin_inhib_mst: float = 0.85
    actin_inhib_lats: float = 0.70
    actin_promote_import: float = 0.80

    eps: float = 1e-9


@dataclass(frozen=True)
class InitialConditions:
    mst: float = 0.15
    lats: float = 0.10
    yap_c: float = 0.60
    yap_p: float = 0.25
    yap_n: float = 0.15
    gene: float = 0.10
    actin: float = 0.50


def _clip01(x: float) -> float:
    return float(min(1.0, max(0.0, x)))


def _hill(x: float, k: float, n: float, eps: float) -> float:
    x = max(0.0, float(x))
    k = max(eps, float(k))
    n = max(eps, float(n))
    num = x**n
    den = k**n + num
    return float(num / max(eps, den))


def pack_state(ic: InitialConditions) -> np.ndarray:
    return np.array(
        [ic.mst, ic.lats, ic.yap_c, ic.yap_p, ic.yap_n, ic.gene, ic.actin],
        dtype=float,
    )


def unpack_state(y: np.ndarray) -> Dict[str, float]:
    return {name: float(y[i]) for i, name in enumerate(STATE_NAMES)}


def rhs(t: float, y: np.ndarray, settings: ModelSettings, p: ModelParams) -> np.ndarray:
    mst, lats, yap_c, yap_p, yap_n, gene, actin = (float(v) for v in y)

    stiffness = _clip01(settings.stiffness)
    density = _clip01(settings.density)

    # actin dynamics: stiffness builds tension; mst activity relaxes tension (contact inhibition)
    d_actin = (
        p.k_a_act * stiffness * (1.0 - actin)
        - p.k_a_inhib_mst * mst * actin
        - p.k_a_decay * actin
    )

    # mst activation: promoted by density (contact inhibition) and low actin tension
    mst_drive = p.k_m_act * (1.0 - p.actin_inhib_mst * actin) + p.k_m_ci * density
    d_mst = mst_drive * (1.0 - mst) - p.k_m_deact * mst

    # lats activation: driven by mst; inhibited by actin tension
    lats_drive = p.k_l_act * mst * (1.0 - p.actin_inhib_lats * actin)
    d_lats = lats_drive * (1.0 - lats) - p.k_l_deact * lats

    # yap pool bookkeeping
    total_yap = max(p.eps, yap_c + yap_p + yap_n)

    # lats phosphorylates cytoplasmic yap; dephos returns to unphos pool
    phos = p.k_phos * lats * yap_c
    dephos = p.k_dephos * yap_p

    # nuclear import favored by actin tension; export is baseline
    # import is also reduced when most of cytoplasmic yap is phosphorylated
    cyt_total = max(p.eps, yap_c + yap_p)
    frac_phos_cyt = yap_p / cyt_total
    import_gate = (1.0 - frac_phos_cyt) * (1.0 + p.actin_promote_import * actin)
    imp = p.k_imp * import_gate * yap_c
    exp = p.k_exp * yap_n

    d_yap_c = p.k_y_syn - phos + dephos - imp + exp - p.k_y_deg_c * yap_c
    d_yap_p = phos - dephos - p.k_y_deg_p * yap_p
    d_yap_n = imp - exp - p.k_y_deg_n * yap_n

    # transcriptional output driven by nuclear yap via hill nonlinearity
    d_gene = p.k_txn * _hill(yap_n / total_yap, p.hill_k, p.hill_n, p.eps) - p.k_gene_deg * gene

    return np.array([d_mst, d_lats, d_yap_c, d_yap_p, d_yap_n, d_gene, d_actin], dtype=float)


def to_dict(settings: ModelSettings, params: ModelParams, ic: InitialConditions) -> Dict[str, Dict[str, float]]:
    return {
        "settings": asdict(settings),
        "params": asdict(params),
        "initial_conditions": asdict(ic),
    }
