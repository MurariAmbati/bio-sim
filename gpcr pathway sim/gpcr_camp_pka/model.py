from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, Dict, List, Tuple

import numpy as np


STATE_NAMES: List[str] = [
    # receptor states
    "R",
    "RL",
    "Rstar",
    "Rp",
    "RArr",
    "Rint",
    # arrestin pool
    "Arr",
    # GRK
    "GRK0",
    "GRKstar",
    # G protein
    "Gs_trimer",
    "Gs_GTP",
    "Gs_aGDP",
    "Gbg",
    # AC pools
    "AC",
    "AC_act",
    "AC_p",
    "AC_p_act",
    # PDE pools
    "PDE_mem",
    "PDE_mem_p",
    "PDE_cyt",
    "PDE_cyt_p",
    # cAMP compartments
    "cAMP_mem",
    "cAMP_cyt",
    # PKA activation ladder
    "PKA0",
    "PKA1",
    "PKA2",
    "PKAstar",
    # downstream
    "CREB",
    "pCREB",
]


def ligand_pulse(t: float, ligand_cfg: Dict[str, float]) -> float:
    """Simple pulse: baseline outside [t_on, t_off), pulse value inside."""
    L0 = float(ligand_cfg.get("L_baseline", 0.0))
    Lp = float(ligand_cfg.get("L_pulse", 0.0))
    t_on = float(ligand_cfg.get("t_on", 0.0))
    t_off = float(ligand_cfg.get("t_off", 0.0))
    return Lp if (t_on <= t < t_off) else L0


def _mm_rate(kcat: float, enzyme: float, substrate: float, Km: float) -> float:
    # Safe Michaelis-Menten-like rate (no division by 0)
    return kcat * enzyme * substrate / (Km + max(substrate, 0.0) + 1e-12)


def ode_rhs(
    t: float,
    y: np.ndarray,
    p: Dict,
    ligand_fn: Callable[[float], float],
) -> np.ndarray:
    """Right-hand side for the GPCR→cAMP/PKA ODE model."""

    idx = {name: i for i, name in enumerate(STATE_NAMES)}

    # unpack state
    (
        R,
        RL,
        Rstar,
        Rp,
        RArr,
        Rint,
        Arr,
        GRK0,
        GRKstar,
        Gs_trimer,
        Gs_GTP,
        Gs_aGDP,
        Gbg,
        AC,
        AC_act,
        AC_p,
        AC_p_act,
        PDE_mem,
        PDE_mem_p,
        PDE_cyt,
        PDE_cyt_p,
        cAMP_mem,
        cAMP_cyt,
        PKA0,
        PKA1,
        PKA2,
        PKAstar,
        CREB,
        pCREB,
    ) = y

    # parameters
    rec = p["receptor"]
    gs = p["Gs"]
    ac = p["AC"]
    camp = p["cAMP"]
    pde = p["PDE"]
    pka = p["PKA"]
    creb_p = p["CREB"]

    L = float(ligand_fn(t))

    # --- receptor binding & activation ---
    v_bind = rec["k_on"] * R * L
    v_unbind = rec["k_off"] * RL
    v_act = rec["k_act"] * RL
    v_deact = rec["k_deact"] * Rstar

    # --- GRK activation by Gbg (proxy for membrane recruitment) ---
    v_grk_act = rec["k_grk_act"] * GRK0 * Gbg
    v_grk_deact = rec["k_grk_deact"] * GRKstar

    # --- receptor phosphorylation to Rp (homologous + heterologous) ---
    v_r_phos_grk = rec["k_grk_phos"] * GRKstar * Rstar
    v_r_phos_pka = rec["k_pka_r_phos"] * PKAstar * Rstar
    v_r_dephos = rec["k_r_dephos"] * Rp

    # --- arrestin binding/internalization/recycling ---
    v_arr_on = rec["k_arr_on"] * Rp * Arr
    v_arr_off = rec["k_arr_off"] * RArr
    v_int = rec["k_int"] * RArr
    v_rec = rec["k_rec"] * Rint

    # --- Gs cycle driven by active receptor ---
    # R* catalyzes activation of Gs trimer: Gs_trimer -> Gs_GTP + Gbg
    v_gef = gs["k_gef"] * Rstar * Gs_trimer
    v_gtpase = gs["k_gtpase"] * Gs_GTP
    v_reassoc = gs["k_reassoc"] * Gs_aGDP * Gbg

    # --- AC activation by Gs_GTP ---
    v_ac_bind = ac["k_bind"] * AC * Gs_GTP
    v_ac_unbind = ac["k_unbind"] * AC_act

    v_acp_bind = ac["k_bind_p"] * AC_p * Gs_GTP
    v_acp_unbind = ac["k_unbind_p"] * AC_p_act

    # --- PKA phosphorylation of AC pools ---
    v_ac_phos = ac["k_pka_ac"] * PKAstar * AC
    v_ac_dephos = ac["k_pp_ac"] * AC_p

    v_acact_phos = ac["k_pka_ac"] * PKAstar * AC_act
    v_acact_dephos = ac["k_pp_ac"] * AC_p_act

    # --- cAMP production near membrane (by AC active pools) ---
    v_camp_prod = ac["k_cat"] * AC_act + ac["k_cat_p"] * AC_p_act

    # --- cAMP degradation (membrane and cytosol PDEs) ---
    Km = pde["Km"]
    v_pde_mem = _mm_rate(pde["k_cat_mem"], PDE_mem, cAMP_mem, Km) + _mm_rate(
        pde["k_cat_mem_p"], PDE_mem_p, cAMP_mem, Km
    )
    v_pde_cyt = _mm_rate(pde["k_cat_cyt"], PDE_cyt, cAMP_cyt, Km) + _mm_rate(
        pde["k_cat_cyt_p"], PDE_cyt_p, cAMP_cyt, Km
    )

    # --- cAMP diffusion between membrane microdomain and cytosol ---
    v_diff = camp["k_diff"] * (cAMP_mem - cAMP_cyt)

    # --- PKA activation ladder via cAMP binding ---
    v_pka01 = pka["k_on1"] * PKA0 * cAMP_cyt
    v_pka10 = pka["k_off1"] * PKA1

    v_pka12 = pka["k_on2"] * PKA1 * cAMP_cyt
    v_pka21 = pka["k_off2"] * PKA2

    v_release = pka["k_release"] * PKA2
    v_inact = pka["k_inact"] * PKAstar

    # --- PKA phosphorylation of PDE pools (feedback) ---
    v_pde_mem_phos = pde["k_pka_pde"] * PKAstar * PDE_mem
    v_pde_mem_dephos = pde["k_pp_pde"] * PDE_mem_p

    v_pde_cyt_phos = pde["k_pka_pde"] * PKAstar * PDE_cyt
    v_pde_cyt_dephos = pde["k_pp_pde"] * PDE_cyt_p

    # --- downstream: CREB phosphorylation ---
    v_creb_phos = creb_p["k_pka"] * PKAstar * CREB
    v_creb_dephos = creb_p["k_pp"] * pCREB

    dydt = np.zeros_like(y)

    # Receptor
    dydt[idx["R"]] = -v_bind + v_unbind + v_deact + v_r_dephos + v_rec
    dydt[idx["RL"]] = v_bind - v_unbind - v_act
    dydt[idx["Rstar"]] = v_act - v_deact - v_r_phos_grk - v_r_phos_pka
    dydt[idx["Rp"]] = v_r_phos_grk + v_r_phos_pka - v_r_dephos - v_arr_on + v_arr_off
    dydt[idx["RArr"]] = v_arr_on - v_arr_off - v_int
    dydt[idx["Rint"]] = v_int - v_rec

    # Arrestin pool
    dydt[idx["Arr"]] = -v_arr_on + v_arr_off + v_int

    # GRK
    dydt[idx["GRK0"]] = -v_grk_act + v_grk_deact
    dydt[idx["GRKstar"]] = v_grk_act - v_grk_deact

    # G protein
    dydt[idx["Gs_trimer"]] = -v_gef + v_reassoc
    dydt[idx["Gs_GTP"]] = v_gef - v_gtpase - v_ac_bind + v_ac_unbind - v_acp_bind + v_acp_unbind
    dydt[idx["Gs_aGDP"]] = v_gtpase - v_reassoc
    # GRK activation is modeled as recruitment/activation by Gβγ, not stoichiometric consumption of Gβγ.
    dydt[idx["Gbg"]] = v_gef - v_reassoc

    # AC pools
    dydt[idx["AC"]] = -v_ac_bind - v_ac_phos + v_ac_unbind + v_ac_dephos
    dydt[idx["AC_act"]] = v_ac_bind - v_ac_unbind - v_acact_phos + v_acact_dephos
    dydt[idx["AC_p"]] = -v_acp_bind + v_ac_phos + v_acp_unbind - v_ac_dephos
    dydt[idx["AC_p_act"]] = v_acp_bind - v_acp_unbind + v_acact_phos - v_acact_dephos

    # PDE pools
    dydt[idx["PDE_mem"]] = -v_pde_mem_phos + v_pde_mem_dephos
    dydt[idx["PDE_mem_p"]] = v_pde_mem_phos - v_pde_mem_dephos

    dydt[idx["PDE_cyt"]] = -v_pde_cyt_phos + v_pde_cyt_dephos
    dydt[idx["PDE_cyt_p"]] = v_pde_cyt_phos - v_pde_cyt_dephos

    # cAMP
    dydt[idx["cAMP_mem"]] = v_camp_prod - v_pde_mem - v_diff

    # cytosolic cAMP includes binding to PKA ladder
    dydt[idx["cAMP_cyt"]] = (
        +v_diff
        - v_pde_cyt
        # PKA binding/unbinding stoichiometry: consumes/releases free cAMP
        - v_pka01
        + v_pka10
        - v_pka12
        + v_pka21
    )

    # PKA ladder
    dydt[idx["PKA0"]] = -v_pka01 + v_pka10 + v_inact
    dydt[idx["PKA1"]] = v_pka01 - v_pka10 - v_pka12 + v_pka21
    dydt[idx["PKA2"]] = v_pka12 - v_pka21 - v_release
    dydt[idx["PKAstar"]] = v_release - v_inact

    # Downstream
    dydt[idx["CREB"]] = -v_creb_phos + v_creb_dephos
    dydt[idx["pCREB"]] = v_creb_phos - v_creb_dephos

    return dydt


# ---------------- Stochastic (SSA) ----------------

@dataclass(frozen=True)
class Reaction:
    name: str
    # stoichiometry updates: list of (species_index, delta)
    updates: Tuple[Tuple[int, int], ...]
    propensity: Callable[[float, np.ndarray, Dict], float]


def build_ssa_reactions(p: Dict, ligand_fn: Callable[[float], float]) -> Tuple[List[str], List[Reaction]]:
    """A reduced Gillespie SSA system approximating the same pathway.

    Notes:
    - Uses mass-action propensities only.
    - Replaces MM PDE rates with linear degradation: k_deg * PDE * cAMP.
    - Produces cAMP as a birth from active AC complexes.

    This is intended as a *single-molecule flavored* option; it won't perfectly match the ODE.
    """

    names = STATE_NAMES
    i = {n: k for k, n in enumerate(names)}

    rec = p["receptor"]
    gs = p["Gs"]
    ac = p["AC"]
    camp = p["cAMP"]
    pde = p["PDE"]
    pka = p["PKA"]

    # choose linearized degradation rates for SSA
    k_deg_mem = float(pde["k_cat_mem"]) / float(pde["Km"])
    k_deg_cyt = float(pde["k_cat_cyt"]) / float(pde["Km"])

    rxns: List[Reaction] = []

    # R + L -> RL (treat L as time-varying constant; not consumed)
    rxns.append(
        Reaction(
            "bind",
            updates=((i["R"], -1), (i["RL"], +1)),
            propensity=lambda t, x, pp: rec["k_on"] * x[i["R"]] * ligand_fn(t),
        )
    )
    rxns.append(
        Reaction(
            "unbind",
            updates=((i["R"], +1), (i["RL"], -1)),
            propensity=lambda t, x, pp: rec["k_off"] * x[i["RL"]],
        )
    )
    rxns.append(
        Reaction(
            "act",
            updates=((i["RL"], -1), (i["Rstar"], +1)),
            propensity=lambda t, x, pp: rec["k_act"] * x[i["RL"]],
        )
    )
    rxns.append(
        Reaction(
            "deact",
            updates=((i["Rstar"], -1), (i["R"], +1)),
            propensity=lambda t, x, pp: rec["k_deact"] * x[i["Rstar"]],
        )
    )

    # Gs activation by R*: Gs_trimer -> Gs_GTP + Gbg
    rxns.append(
        Reaction(
            "Gs_GEF",
            updates=((i["Gs_trimer"], -1), (i["Gs_GTP"], +1), (i["Gbg"], +1)),
            propensity=lambda t, x, pp: gs["k_gef"] * x[i["Rstar"]] * x[i["Gs_trimer"]],
        )
    )
    rxns.append(
        Reaction(
            "Gs_GTPase",
            updates=((i["Gs_GTP"], -1), (i["Gs_aGDP"], +1)),
            propensity=lambda t, x, pp: gs["k_gtpase"] * x[i["Gs_GTP"]],
        )
    )
    rxns.append(
        Reaction(
            "reassoc",
            updates=((i["Gs_aGDP"], -1), (i["Gbg"], -1), (i["Gs_trimer"], +1)),
            propensity=lambda t, x, pp: gs["k_reassoc"] * x[i["Gs_aGDP"]] * x[i["Gbg"]],
        )
    )

    # AC activation by Gs_GTP
    rxns.append(
        Reaction(
            "AC_bind",
            updates=((i["AC"], -1), (i["AC_act"], +1)),
            propensity=lambda t, x, pp: ac["k_bind"] * x[i["AC"]] * x[i["Gs_GTP"]],
        )
    )
    rxns.append(
        Reaction(
            "AC_unbind",
            updates=((i["AC"], +1), (i["AC_act"], -1)),
            propensity=lambda t, x, pp: ac["k_unbind"] * x[i["AC_act"]],
        )
    )

    # cAMP birth from active AC (mem compartment)
    rxns.append(
        Reaction(
            "cAMP_prod",
            updates=((i["cAMP_mem"], +1),),
            propensity=lambda t, x, pp: ac["k_cat"] * x[i["AC_act"]],
        )
    )

    # linear PDE degradation
    rxns.append(
        Reaction(
            "cAMP_deg_mem",
            updates=((i["cAMP_mem"], -1),),
            propensity=lambda t, x, pp: k_deg_mem * x[i["PDE_mem"]] * x[i["cAMP_mem"]],
        )
    )
    rxns.append(
        Reaction(
            "cAMP_deg_cyt",
            updates=((i["cAMP_cyt"], -1),),
            propensity=lambda t, x, pp: k_deg_cyt * x[i["PDE_cyt"]] * x[i["cAMP_cyt"]],
        )
    )

    # diffusion as two reactions (mem->cyt and cyt->mem)
    rxns.append(
        Reaction(
            "diff_mem_to_cyt",
            updates=((i["cAMP_mem"], -1), (i["cAMP_cyt"], +1)),
            propensity=lambda t, x, pp: camp["k_diff"] * x[i["cAMP_mem"]],
        )
    )
    rxns.append(
        Reaction(
            "diff_cyt_to_mem",
            updates=((i["cAMP_cyt"], -1), (i["cAMP_mem"], +1)),
            propensity=lambda t, x, pp: camp["k_diff"] * x[i["cAMP_cyt"]],
        )
    )

    # minimal PKA activation via cAMP binding ladder
    rxns.append(
        Reaction(
            "PKA0_bind",
            updates=((i["PKA0"], -1), (i["PKA1"], +1), (i["cAMP_cyt"], -1)),
            propensity=lambda t, x, pp: pka["k_on1"] * x[i["PKA0"]] * x[i["cAMP_cyt"]],
        )
    )
    rxns.append(
        Reaction(
            "PKA1_bind",
            updates=((i["PKA1"], -1), (i["PKA2"], +1), (i["cAMP_cyt"], -1)),
            propensity=lambda t, x, pp: pka["k_on2"] * x[i["PKA1"]] * x[i["cAMP_cyt"]],
        )
    )
    rxns.append(
        Reaction(
            "PKA_release",
            updates=((i["PKA2"], -1), (i["PKAstar"], +1)),
            propensity=lambda t, x, pp: pka["k_release"] * x[i["PKA2"]],
        )
    )
    rxns.append(
        Reaction(
            "PKA_inact",
            updates=((i["PKAstar"], -1), (i["PKA0"], +1)),
            propensity=lambda t, x, pp: pka["k_inact"] * x[i["PKAstar"]],
        )
    )

    return names, rxns
