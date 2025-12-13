from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd
import streamlit as st

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from mapk_erk import Stimulus, default_initial_conditions, default_parameters, simulate_ode, simulate_ssa


st.set_page_config(page_title="MAPK/ERK simulator", layout="wide")

st.title("MAPK / ERK pathway simulator")

with st.sidebar:
    st.header("Pages")
    page = st.radio("Go to", ["Simulator", "Interpretation & Use"], label_visibility="collapsed")

    st.header("Simulation")
    mode = st.selectbox("Mode", ["ODE (deterministic)", "SSA (stochastic, single-cell)"])

    st.subheader("Time")
    t_end = st.slider("End time", min_value=10.0, max_value=300.0, value=120.0, step=5.0)

    st.subheader("Stimulus")
    stim_kind = st.selectbox("Stimulus type", ["step", "pulse", "exp"])
    stim_amp = st.slider("Amplitude", min_value=0.0, max_value=5.0, value=1.0, step=0.1)
    stim_ton = st.slider("t_on", min_value=0.0, max_value=200.0, value=0.0, step=1.0)
    stim_toff = st.slider("t_off (pulse)", min_value=0.0, max_value=200.0, value=30.0, step=1.0)
    stim_tau = st.slider("tau (exp)", min_value=0.5, max_value=60.0, value=10.0, step=0.5)

    st.subheader("Key feedback")
    fb_alpha = st.slider("ERK→Ras feedback strength (alpha)", 0.0, 10.0, 2.0, 0.1)
    dusp_syn_max = st.slider("DUSP inducible synthesis", 0.0, 1.5, 0.4, 0.02)

    st.subheader("Stochastic options")
    volume_fL = st.slider("Cell volume (fL)", min_value=0.1, max_value=10.0, value=1.0, step=0.1)
    seed = st.number_input("Random seed", min_value=0, max_value=10_000_000, value=1, step=1)


def _render_interpretation_page() -> None:
    st.subheader("How to interpret the outputs")

    c1, c2 = st.columns(2)
    with c1:
        st.info(
            "**ppERK (ERK-PP)** is the main readout: it’s commonly used as a proxy for ERK pathway activity.\n\n"
            "- **Peak** ppERK: how strongly the pathway activates.\n"
            "- **AUC** ppERK: how much total signaling occurred over time.\n"
            "- **Transient vs sustained**: whether signaling adapts back down or stays high depends on feedback strength and stimulus shape."
        )

        st.info(
            "**Stimulus shapes**:\n\n"
            "- **step**: sustained ligand after `t_on`.\n"
            "- **pulse**: ligand between `t_on` and `t_off`.\n"
            "- **exp**: ligand turns on then decays with time constant `tau`."
        )

    with c2:
        st.success(
            "**Negative feedback in this model**:\n\n"
            "- **ERK→Ras feedback (alpha)**: higher values suppress Ras activation as ppERK rises, encouraging adaptation.\n"
            "- **DUSP induction**: ppERK induces DUSP, which increases ERK dephosphorylation and adds delayed shutoff."
        )

        st.warning(
            "**Important modeling note**:\n\n"
            "The **ODE mode** is a deterministic, population-average view. The **SSA mode** is a simplified single-cell approximation (mass-action SSA) and is best used to explore qualitative variability and noise-driven differences, not to match every ODE rate law exactly."
        )

    st.subheader("What this simulator can help with")
    c3, c4 = st.columns(2)
    with c3:
        st.info(
            "- Exploring how **dose** and **timing** of ligand change ppERK dynamics\n"
            "- Testing how feedback strength shifts **peak**, **duration**, and **adaptation**\n"
            "- Comparing pulse vs sustained inputs for **cell-fate-like** signaling patterns"
        )
    with c4:
        st.info(
            "- Sensitivity intuition: which layers (Ras/Raf/MEK/ERK) bottleneck responses\n"
            "- Single-cell intuition: how stochasticity can create **cell-to-cell variability**\n"
            "- Designing in silico perturbations (e.g., increase DUSP induction to mimic stronger phosphatase response)"
        )


if page == "Interpretation & Use":
    _render_interpretation_page()
    st.stop()

p = default_parameters()
# update a couple of parameters from UI
p.fb_Ras_alpha = float(fb_alpha)
p.dusp_syn_max = float(dusp_syn_max)

x0 = default_initial_conditions(p)
stim = Stimulus(kind=stim_kind, amplitude=float(stim_amp), t_on=float(stim_ton), t_off=float(stim_toff), tau=float(stim_tau))

run = st.button("Run simulation", type="primary")

if run:
    with st.spinner("Simulating..."):
        if mode.startswith("ODE"):
            res = simulate_ode(p=p, x0=x0, stim=stim, t_span=(0.0, float(t_end)), n_points=int(max(201, t_end * 10)))
            df = res.to_frame()
        else:
            res = simulate_ssa(
                p=p,
                x0=x0,
                stim=stim,
                t_span=(0.0, float(t_end)),
                volume_fL=float(volume_fL),
                seed=int(seed),
            )
            df = res.to_frame()

    st.success("Done")

    c1, c2 = st.columns([2, 1])
    with c1:
        st.subheader("ERK phosphorylation")
        plot_df = df[["t", "ERK", "ERK_P", "ERK_PP"]]
        st.line_chart(plot_df, x="t")

    with c2:
        st.subheader("Summary")
        peak = float(df["ERK_PP"].max())
        auc = float(((df["ERK_PP"].values[:-1] + df["ERK_PP"].values[1:]) / 2 * (df["t"].values[1:] - df["t"].values[:-1])).sum())
        st.metric("Peak ppERK", f"{peak:.4g}")
        st.metric("AUC ppERK", f"{auc:.4g}")

    st.subheader("All species")
    st.dataframe(df, use_container_width=True)

    csv = df.to_csv(index=False).encode("utf-8")
    st.download_button("Download CSV", data=csv, file_name="mapk_erk_simulation.csv", mime="text/csv")

else:
    st.caption("Adjust stimulus/feedback in the sidebar, then click Run simulation.")
