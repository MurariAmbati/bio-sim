from __future__ import annotations

import copy
from pathlib import Path
from typing import Dict, Tuple

import matplotlib.pyplot as plt
import numpy as np
import streamlit as st
import yaml
from scipy.integrate import solve_ivp

from model import STATE_NAMES, ligand_pulse, ode_rhs


def _get_param(p: Dict, path: Tuple[str, ...]) -> float:
    cur = p
    for key in path:
        cur = cur[key]
    return float(cur)


def _set_param(p: Dict, path: Tuple[str, ...], value: float) -> None:
    cur = p
    for key in path[:-1]:
        cur = cur[key]
    cur[path[-1]] = float(value)


@st.cache_data
def load_params(path: str) -> Dict:
    with Path(path).open("r", encoding="utf-8") as f:
        return yaml.safe_load(f)


def build_y0(p: Dict) -> np.ndarray:
    init = p["init"]
    return np.array([float(init[name]) for name in STATE_NAMES], dtype=float)


def run_ode(p: Dict) -> Tuple[np.ndarray, np.ndarray]:
    lig_cfg = p["ligand"]
    lig_fn = lambda t: ligand_pulse(t, lig_cfg)

    y0 = build_y0(p)
    t0 = float(p["sim"]["t0"])
    t_end = float(p["sim"]["t_end"])
    max_step = float(p["sim"].get("max_step", 0.5))

    rhs = lambda t, y: ode_rhs(t, y, p, lig_fn)

    sol = solve_ivp(
        rhs,
        t_span=(t0, t_end),
        y0=y0,
        method="LSODA",
        max_step=max_step,
        rtol=1e-7,
        atol=1e-9,
    )
    if not sol.success:
        raise RuntimeError(sol.message)

    return sol.t, sol.y


def _plot_ode(t: np.ndarray, Y: np.ndarray, p: Dict) -> plt.Figure:
    i = {n: k for k, n in enumerate(STATE_NAMES)}
    lig_cfg = p["ligand"]
    lig_fn = lambda tt: ligand_pulse(tt, lig_cfg)
    L = np.array([lig_fn(tt) for tt in t])

    fig, axes = plt.subplots(4, 1, figsize=(10, 10), sharex=True)

    axes[0].plot(t, L)
    axes[0].set_ylabel("ligand")

    axes[1].plot(t, Y[i["Rstar"], :], label="r*")
    axes[1].plot(t, Y[i["RArr"], :], label="r:arr")
    axes[1].plot(t, Y[i["Rint"], :], label="r_int")
    axes[1].set_ylabel("receptor")
    axes[1].legend(loc="best")

    axes[2].plot(t, Y[i["cAMP_mem"], :], label="camp_mem")
    axes[2].plot(t, Y[i["cAMP_cyt"], :], label="camp_cyt")
    axes[2].set_ylabel("camp")
    axes[2].legend(loc="best")

    axes[3].plot(t, Y[i["PKAstar"], :], label="pka*")
    axes[3].plot(t, Y[i["pCREB"], :], label="pcreb")
    axes[3].set_ylabel("pka / creb")
    axes[3].set_xlabel("time (s)")
    axes[3].legend(loc="best")

    fig.tight_layout()
    return fig


def _plot_ssa(t: np.ndarray, X: np.ndarray) -> plt.Figure:
    i = {n: k for k, n in enumerate(STATE_NAMES)}
    fig, ax = plt.subplots(1, 1, figsize=(10, 4))
    ax.step(t, X[:, i["cAMP_mem"]], where="post", label="camp_mem")
    ax.step(t, X[:, i["cAMP_cyt"]], where="post", label="camp_cyt")
    ax.step(t, X[:, i["PKAstar"]], where="post", label="pka*")
    ax.set_xlabel("time (s)")
    ax.set_ylabel("counts")
    ax.legend(loc="best")
    fig.tight_layout()
    return fig


def main() -> None:
    st.set_page_config(page_title="gpcr → camp/pka", layout="wide")
    st.title("gpcr → camp / pka")

    p_base = load_params("params.yaml")

    with st.sidebar:
        mode = st.radio("simulation", ["ode", "ssa"], index=0)
        seed = st.number_input("ssa seed", min_value=0, max_value=10_000_000, value=1, step=1)

        st.divider()
        st.subheader("params")
        st.caption("15 sliders")

        # exactly 15 sliders
        with st.expander("ligand", expanded=True):
            L_pulse = st.slider("l_pulse", 0.0, 1.0, _get_param(p_base, ("ligand", "L_pulse")), 0.01)
            t_on = st.slider(
                "t_on (s)",
                0.0,
                float(p_base["sim"]["t_end"]),
                _get_param(p_base, ("ligand", "t_on")),
                1.0,
            )
            t_off = st.slider(
                "t_off (s)",
                0.0,
                float(p_base["sim"]["t_end"]),
                _get_param(p_base, ("ligand", "t_off")),
                1.0,
            )

        with st.expander("receptor", expanded=False):
            k_act = st.slider("k_act", 0.0, 3.0, _get_param(p_base, ("receptor", "k_act")), 0.01)
            k_deact = st.slider("k_deact", 0.0, 2.0, _get_param(p_base, ("receptor", "k_deact")), 0.01)
            k_grk_phos = st.slider(
                "k_grk_phos", 0.0, 5.0, _get_param(p_base, ("receptor", "k_grk_phos")), 0.05
            )

        with st.expander("gs / ac", expanded=False):
            k_gef = st.slider("gs k_gef", 0.0, 6.0, _get_param(p_base, ("Gs", "k_gef")), 0.05)
            k_cat = st.slider("ac k_cat", 0.0, 10.0, _get_param(p_base, ("AC", "k_cat")), 0.05)
            k_cat_p = st.slider("ac k_cat_p", 0.0, 12.0, _get_param(p_base, ("AC", "k_cat_p")), 0.05)

        with st.expander("camp / pde", expanded=False):
            k_diff = st.slider("camp k_diff", 0.0, 2.0, _get_param(p_base, ("cAMP", "k_diff")), 0.01)
            Km = st.slider("pde km", 0.1, 10.0, _get_param(p_base, ("PDE", "Km")), 0.05)
            k_cat_cyt = st.slider(
                "pde k_cat_cyt", 0.0, 12.0, _get_param(p_base, ("PDE", "k_cat_cyt")), 0.05
            )
            k_pka_pde = st.slider(
                "pde k_pka_pde", 0.0, 2.0, _get_param(p_base, ("PDE", "k_pka_pde")), 0.01
            )

        with st.expander("pka / creb", expanded=False):
            k_on1 = st.slider("pka k_on1", 0.0, 10.0, _get_param(p_base, ("PKA", "k_on1")), 0.05)
            k_release = st.slider(
                "pka k_release", 0.0, 5.0, _get_param(p_base, ("PKA", "k_release")), 0.05
            )
            creb_k_pka = st.slider(
                "creb k_pka", 0.0, 3.0, _get_param(p_base, ("CREB", "k_pka")), 0.01
            )

        with st.expander("ssa", expanded=False):
            ssa_scale = st.slider(
                "ssa.scale",
                50.0,
                2000.0,
                float(p_base.get("ssa", {}).get("scale", 400.0)),
                10.0,
            )

        run = st.button("run", type="primary")

    # apply overrides on a copy
    p = copy.deepcopy(p_base)

    if t_off <= t_on:
        t_off = min(float(p["sim"]["t_end"]), t_on + 1.0)

    _set_param(p, ("ligand", "L_pulse"), L_pulse)
    _set_param(p, ("ligand", "t_on"), t_on)
    _set_param(p, ("ligand", "t_off"), t_off)

    _set_param(p, ("receptor", "k_act"), k_act)
    _set_param(p, ("receptor", "k_deact"), k_deact)
    _set_param(p, ("receptor", "k_grk_phos"), k_grk_phos)

    _set_param(p, ("Gs", "k_gef"), k_gef)
    _set_param(p, ("AC", "k_cat"), k_cat)
    _set_param(p, ("AC", "k_cat_p"), k_cat_p)

    _set_param(p, ("cAMP", "k_diff"), k_diff)
    _set_param(p, ("PDE", "Km"), Km)
    _set_param(p, ("PDE", "k_cat_cyt"), k_cat_cyt)
    _set_param(p, ("PDE", "k_pka_pde"), k_pka_pde)

    _set_param(p, ("PKA", "k_on1"), k_on1)
    _set_param(p, ("PKA", "k_release"), k_release)
    _set_param(p, ("CREB", "k_pka"), creb_k_pka)

    p.setdefault("ssa", {})
    p["ssa"]["scale"] = float(ssa_scale)

    st.caption("base: params.yaml (with slider overrides)")

    if not run:
        st.info("adjust sliders and click run")
        return

    if mode == "ode":
        with st.spinner("integrating ode..."):
            t, Y = run_ode(p)

        i = {n: k for k, n in enumerate(STATE_NAMES)}
        st.metric("peak camp_cyt", f"{float(Y[i['cAMP_cyt'], :].max()):.3f}")
        st.metric("peak pka*", f"{float(Y[i['PKAstar'], :].max()):.3f}")
        st.pyplot(_plot_ode(t, Y, p), clear_figure=True)

    else:
        from run_ssa import build_x0, gillespie  # noqa: wps433

        with st.spinner("running ssa..."):
            t0 = float(p["sim"]["t0"])
            t_end = float(p["sim"]["t_end"])
            x0 = build_x0(p)
            traj = gillespie(t0=t0, t_end=t_end, x0=x0, p=p, seed=int(seed))

        st.metric("ssa events", f"{traj.x.shape[0] - 1}")
        st.pyplot(_plot_ssa(traj.t, traj.x), clear_figure=True)


if __name__ == "__main__":
    main()
