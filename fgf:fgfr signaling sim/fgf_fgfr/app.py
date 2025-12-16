from __future__ import annotations

import sys
from pathlib import Path
from dataclasses import replace

import numpy as np
import streamlit as st

APP_DIR = Path(__file__).resolve().parent
SRC_DIR = APP_DIR / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from fgf_fgfr.analysis.dose_response import dose_response
from fgf_fgfr.analysis.sensitivity import local_sensitivity
from fgf_fgfr.params import ContextName, InitialConditions, ModelParams, preset
from fgf_fgfr.simulate import outcomes_from_state, simulate, to_dataframe
from fgf_fgfr.viz.plots import (
    dose_response_plot,
    outcome_bar,
    sensitivity_bar,
    timeseries_plot,
)


st.set_page_config(page_title="FGF/FGFR Signaling", layout="wide")

st.title("FGF/FGFR signaling: interactive simulator")

with st.sidebar:
    st.header("Context")
    context: ContextName = st.selectbox(
        "Biological context",
        options=["development", "angiogenesis", "tissue_repair"],
        index=1,
    )
    ligand = st.slider("FGF dose (nM)", 0.0, 200.0, 10.0, 1.0)
    receptor = st.slider("FGFR abundance (nM-eq)", 0.0, 200.0, 50.0, 1.0)

    st.divider()
    st.header("Perturbations")
    fgfr_inhib = st.slider("FGFR inhibitor", 0.0, 1.0, 0.0, 0.05)
    mek_inhib = st.slider("MEK inhibitor", 0.0, 1.0, 0.0, 0.05)
    pi3k_inhib = st.slider("PI3K inhibitor", 0.0, 1.0, 0.0, 0.05)
    plc_inhib = st.slider("PLCγ inhibitor", 0.0, 1.0, 0.0, 0.05)
    stat_inhib = st.slider("STAT inhibitor", 0.0, 1.0, 0.0, 0.05)

    st.divider()
    st.header("Simulation")
    t_end_min = st.slider("Duration (min)", 1, 240, 60, 1)
    n_points = st.slider("Resolution (points)", 200, 2000, 700, 50)

pr = preset(context, ligand_nM=float(ligand))
params = replace(
    pr.params,
    fgfr_inhib=float(fgfr_inhib),
    mek_inhib=float(mek_inhib),
    pi3k_inhib=float(pi3k_inhib),
    plc_inhib=float(plc_inhib),
    stat_inhib=float(stat_inhib),
)
ic = replace(pr.ic, L=float(ligand), R=float(receptor))

colA, colB = st.columns([1.2, 1.0], gap="large")

with colA:
    st.subheader("Pathway activity (time series)")
    try:
        state = simulate(params, ic, t_end_s=float(t_end_min) * 60.0, n_points=int(n_points))
    except Exception as e:
        st.error(str(e))
        st.stop()

    df = to_dataframe(state)

    fig = timeseries_plot(df, ["Rp", "A", "ERKp", "AKTp", "PLCa", "STATp"], title="Core signaling")
    st.plotly_chart(fig, use_container_width=True)

    with st.expander("Raw simulation data"):
        st.dataframe(df, use_container_width=True)
        st.download_button(
            "Download CSV",
            data=df.to_csv(index=False).encode("utf-8"),
            file_name="fgf_fgfr_timeseries.csv",
            mime="text/csv",
        )

with colB:
    st.subheader("Predicted outcomes")
    outcomes = outcomes_from_state(state)
    st.plotly_chart(outcome_bar(outcomes), use_container_width=True)

    st.caption(
        "Outcome scores are heuristic mappings from terminal pathway signals (ERK/AKT/PLC/STAT)."
    )

st.divider()

colC, colD = st.columns([1.0, 1.0], gap="large")

with colC:
    st.subheader("Dose-response")
    y_name = st.selectbox("Outcome to scan", options=["proliferation", "survival", "migration", "angiogenesis"], index=0)
    doses = np.unique(np.round(np.geomspace(0.1, 200.0, 18), 3))
    dr = dose_response(doses, params, ic, t_end_s=float(t_end_min) * 60.0)
    st.plotly_chart(dose_response_plot(dr, y=y_name), use_container_width=True)

with colD:
    st.subheader("Sensitivity")
    outcome = st.selectbox("Outcome for sensitivity", options=["proliferation", "survival", "migration", "angiogenesis"], index=0, key="sens")
    sens = local_sensitivity(params, ic, outcome=outcome, rel_step=0.1, t_end_s=float(t_end_min) * 60.0)
    st.plotly_chart(sensitivity_bar(sens, top_n=10), use_container_width=True)

st.divider()

st.subheader("Model notes")
st.write(
    "This app is a compact mechanistic toy-model capturing canonical FGF/FGFR branches (MAPK, PI3K/AKT, PLCγ/Ca2+, STAT) "
    "plus receptor desensitization via internalization. It’s intended for interactive exploration and hypothesis sketching."
)

st.code(
    """Key species (ODE state):
- L: free ligand (FGF)
- R: surface receptor
- LR: ligand-bound receptor
- D: dimerized complex
- Rp: active/phosphorylated receptor
- A: recruited adaptor complex (FRS2/GRB2/SOS lump)
- ERKp, AKTp, PLCa, STATp: active downstream signals
- Rint: internalized receptor pool
""",
    language="text",
)
