from __future__ import annotations

import numpy as np
import pandas as pd
import streamlit as st

from tcr_signaling.model import ModelParams, simulate
from tcr_signaling.pathway import build_pathway_graph
from tcr_signaling.plots import plot_network, plot_timeseries


st.set_page_config(page_title="tcr signaling", layout="wide")

st.title("tcr signaling")
st.caption("simplified tcr→nf-κb / mapk / ca²⁺ dynamics with interactive parameters")

with st.sidebar:
    st.header("inputs")
    antigen = st.slider("antigen (tcr stimulus)", 0.0, 1.0, 0.8, 0.01)
    cd28 = st.slider("cd28 costimulation", 0.0, 1.0, 0.5, 0.01)

    st.header("regulators")
    lck = st.slider("lck activity", 0.0, 2.0, 1.0, 0.01)
    shp1 = st.slider("shp1 inhibition", 0.0, 2.0, 0.6, 0.01)

    st.header("time")
    t_end = st.slider("t end (min)", 1.0, 60.0, 30.0, 1.0)
    n_points = st.slider("points", 100, 2000, 600, 50)

    st.header("kinetics")
    strength = st.slider("global signal strength", 0.1, 5.0, 1.0, 0.05)
    decay = st.slider("global decay", 0.1, 5.0, 1.0, 0.05)

params = ModelParams(
    antigen=antigen,
    cd28=cd28,
    lck=lck,
    shp1=shp1,
    strength=strength,
    decay=decay,
)

col_a, col_b = st.columns([1, 1])

with col_a:
    st.subheader("pathway")
    g = build_pathway_graph()
    st.plotly_chart(plot_network(g), width="stretch")

with col_b:
    st.subheader("dynamics")
    t = np.linspace(0.0, float(t_end), int(n_points))
    df = simulate(t, params)
    st.plotly_chart(plot_timeseries(df), width="stretch")
    st.dataframe(df.tail(10), width="stretch")

st.divider()
with st.expander("notes", expanded=False):
    st.write(
        "this model is intentionally compact: it captures canonical information flow (tcr→zap70/lat→plcγ/calcineurin, ras→erk, pkcθ/ikk→nf-κb) "
        "with a few tunable inputs and negative regulation. it is not a curated quantitative model and should not be used for prediction."
    )
