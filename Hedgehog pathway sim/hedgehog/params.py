from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class CoreParams:
    # receptor + transduction
    k_hh: float = 6.0
    alpha_ptch: float = 8.0

    # gli activation nonlinearity
    hill_n: float = 3.0
    hill_k: float = 0.45


@dataclass(frozen=True)
class OdeParams:
    # timescales (1/time)
    k_p_prod: float = 0.35
    k_p_deg: float = 0.18
    k_s_act: float = 1.8
    k_s_deg: float = 1.2
    k_g_act: float = 2.0
    k_g_deg: float = 0.9

    # feedback: gli -> ptch
    fb_strength: float = 1.2
    fb_k: float = 0.35
    fb_n: float = 2.0

    # extended layer: hhip ligand sink induced by gliA
    hhip_strength: float = 0.9
    hhip_k: float = 0.35
    hhip_n: float = 2.0
    k_h_prod: float = 0.25
    k_h_deg: float = 0.12

    # gli processing (activator vs repressor)
    k_a_act: float = 2.2
    k_a_deg: float = 0.8
    k_r_prod: float = 1.2
    k_r_deg: float = 0.7
    r_weight: float = 0.9

    # transduction mapping
    core: CoreParams = CoreParams()


@dataclass(frozen=True)
class SpatialParams:
    # ligand diffusion / decay in 1d tissue (dimensionless)
    D: float = 0.03
    k_decay: float = 0.02
    source_hh: float = 1.0
    source_width: float = 0.08  # fraction of domain

    # optional: hhip sink (produced by gli) shapes hh gradient
    sink_strength: float = 0.8
    k_hhip_prod: float = 0.25
    k_hhip_deg: float = 0.10
    hhip_k: float = 0.35
    hhip_n: float = 2.0

    # receptor/transduction + gli nonlinearity
    core: CoreParams = CoreParams()

    # gli relaxation
    tau_gli: float = 6.0
