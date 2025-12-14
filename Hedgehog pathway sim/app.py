from __future__ import annotations

import numpy as np
import pandas as pd
import streamlit as st

from hedgehog.core import CoreMapping
from hedgehog.network import NetworkSimulator
from hedgehog.ode import FeedbackDynamicsModel
from hedgehog.params import CoreParams, OdeParams, SpatialParams
from hedgehog.sweeps import gli_surface
from hedgehog.tissue import Tissue1DModel
from hedgehog.viz import heatmap_figure, phase_figure, profiles_figure, timeseries_figure


st.set_page_config(page_title="hedgehog", page_icon="🧬", layout="centered")


@st.cache_data(show_spinner=False)
def steady_state(hh: float, ptch_level: float, core: CoreParams) -> dict[str, float]:
    return CoreMapping(core).steady_state(hh, ptch_level)


@st.cache_data(show_spinner=False)
def simulate_feedback_ode(
    hh: float,
    p: OdeParams,
    t_end: float,
    dt: float,
    ptch_ko: bool,
    smo_ko: bool,
    gli_ko: bool,
) -> pd.DataFrame:
    return FeedbackDynamicsModel(p).simulate_simple(hh, t_end=t_end, dt=dt, ptch_ko=ptch_ko, smo_ko=smo_ko, gli_ko=gli_ko)


@st.cache_data(show_spinner=False)
def simulate_extended_network_integrated(
    hh: float,
    p: OdeParams,
    t_end: float,
    dt: float,
    method: str,
    ptch_ko: bool,
    smo_ko: bool,
    gli_ko: bool,
) -> pd.DataFrame:
    return NetworkSimulator(p).simulate_extended(
        hh=hh,
        t_end=t_end,
        dt=dt,
        method=method,
        ptch_ko=ptch_ko,
        smo_ko=smo_ko,
        gli_ko=gli_ko,
    )


@st.cache_data(show_spinner=False)
def simulate_spatial(
    p: SpatialParams,
    nx: int,
    t_end: float,
    dt: float,
    ptch_level: float,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    return Tissue1DModel(p).simulate(nx=nx, t_end=t_end, dt=dt, ptch_level=ptch_level)


@st.cache_data(show_spinner=False)
def compute_gli_surface(
    core: CoreParams,
    ptch_min: float,
    ptch_max: float,
    hh_points: int,
    ptch_points: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    mapping = CoreMapping(core)
    return gli_surface(mapping, ptch_min=ptch_min, ptch_max=ptch_max, hh_points=hh_points, ptch_points=ptch_points)


st.title("hedgehog pathway")
st.caption("developmental patterning · stem cell maintenance · cancer")

st.write(
    "advanced qualitative sandbox with three layers: steady-state mapping, "
    "feedback dynamics (gli → ptch), and a 1d tissue simulation (diffusion + readout). "
    "all values are dimensionless and not fitted to a dataset."
)

with st.sidebar:
    mode = st.radio("mode", ["steady state", "feedback dynamics", "1d tissue"], index=1)

    st.divider()
    st.header("core mapping")
    hh = st.slider("hh ligand", 0.0, 1.0, 0.35, 0.01)
    ptch_level = st.slider("baseline ptch level", 0.0, 2.5, 0.8, 0.05)
    k_hh = st.slider("hh potency (k_hh)", 0.1, 20.0, 6.0, 0.1)
    alpha_ptch = st.slider("ptch repression strength (alpha)", 0.1, 20.0, 8.0, 0.1)
    hill_n = st.slider("gli hill n", 1.0, 6.0, 3.0, 0.5)
    hill_k = st.slider("gli half-activation k", 0.05, 0.95, 0.45, 0.01)

core = CoreParams(k_hh=k_hh, alpha_ptch=alpha_ptch, hill_n=hill_n, hill_k=hill_k)
mapping = CoreMapping(core)

if mode == "steady state":
    ss = steady_state(hh, ptch_level, core)

    c1, c2, c3 = st.columns(3)
    c1.metric("ptch repression", f"{ss['ptch_repression']:.3f}")
    c2.metric("smo", f"{ss['smo']:.3f}")
    c3.metric("gli* (ss)", f"{ss['gli_ss']:.3f}")

    dose = mapping.dose_response(ptch_level, points=201)
    fig = timeseries_figure(
        dose,
        x="hh",
        y_cols=["smo", "gli_ss"],
        height=320,
        x_title="hh",
        y_title="activity (0..1)",
    )
    st.subheader("dose response")
    st.plotly_chart(fig, use_container_width=True)

    with st.expander("gli surface (hh, ptch)", expanded=False):
        ptch_min, ptch_max = st.slider("ptch range", 0.0, 5.0, (0.0, 2.5), 0.1)
        hh_points = st.select_slider("hh points", options=[41, 61, 81, 101], value=81)
        ptch_points = st.select_slider("ptch points", options=[41, 61, 81], value=61)

        hh_grid, ptch_grid, gli = compute_gli_surface(core, ptch_min, ptch_max, hh_points, ptch_points)
        fig_s = heatmap_figure(gli, hh_grid, ptch_grid, colorscale="Viridis", height=340, x_title="hh", y_title="ptch")
        st.plotly_chart(fig_s, use_container_width=True)

elif mode == "feedback dynamics":
    with st.sidebar:
        st.divider()
        st.header("ode dynamics")
        t_end = st.slider("t_end", 5.0, 300.0, 80.0, 5.0)
        dt = st.select_slider("dt", options=[0.01, 0.02, 0.05, 0.1, 0.2], value=0.05)
        st.caption("smaller dt is slower but smoother")

        st.divider()
        st.header("integrator")
        method = st.selectbox("method", ["rk4", "euler"], index=0)

        st.divider()
        st.header("feedback")
        fb_strength = st.slider("gli→ptch strength", 0.0, 4.0, 1.2, 0.05)
        fb_k = st.slider("feedback k", 0.05, 0.95, 0.35, 0.01)
        fb_n = st.slider("feedback n", 1.0, 6.0, 2.0, 0.5)

        st.divider()
        st.header("more detail")
        extended = st.checkbox("gli a/r + hhip sink", value=True)
        hhip_strength = st.slider("hhip sink strength", 0.0, 3.0, 0.9, 0.05)
        r_weight = st.slider("gliR weight", 0.0, 1.5, 0.9, 0.05)

        st.divider()
        st.header("perturbations")
        ptch_ko = st.checkbox("ptch knockout", value=False)
        smo_ko = st.checkbox("smo knockout", value=False)
        gli_ko = st.checkbox("gli knockout", value=False)

    ode_p = OdeParams(
        fb_strength=fb_strength,
        fb_k=fb_k,
        fb_n=fb_n,
        hhip_strength=hhip_strength,
        r_weight=r_weight,
        core=core,
    )

    if extended:
        df = simulate_extended_network_integrated(
            hh,
            ode_p,
            t_end=t_end,
            dt=dt,
            method=method,
            ptch_ko=ptch_ko,
            smo_ko=smo_ko,
            gli_ko=gli_ko,
        )

        c1, c2, c3 = st.columns(3)
        c1.metric("ptch (final)", f"{df['ptch'].iloc[-1]:.3f}")
        c2.metric("gli_net (final)", f"{df['gli_net'].iloc[-1]:.3f}")
        c3.metric("hhip (final)", f"{df['hhip'].iloc[-1]:.3f}")

        fig = timeseries_figure(
            df,
            x="t",
            y_cols=["ptch", "smo", "gliA", "gliR", "gli_net", "hhip"],
            height=390,
            x_title="time",
            y_title="level / activity",
        )
        st.subheader("network dynamics")
        st.plotly_chart(fig, use_container_width=True)

        fig2 = phase_figure(df, x="ptch", y="gli_net", height=320, x_title="ptch", y_title="gli_net")
        st.subheader("phase view")
        st.plotly_chart(fig2, use_container_width=True)

    else:
        df = simulate_feedback_ode(hh, ode_p, t_end=t_end, dt=dt, ptch_ko=ptch_ko, smo_ko=smo_ko, gli_ko=gli_ko)

        c1, c2, c3 = st.columns(3)
        c1.metric("ptch (final)", f"{df['ptch'].iloc[-1]:.3f}")
        c2.metric("smo (final)", f"{df['smo'].iloc[-1]:.3f}")
        c3.metric("gli (final)", f"{df['gli'].iloc[-1]:.3f}")

        fig = timeseries_figure(df, x="t", y_cols=["ptch", "smo", "gli"], height=360, x_title="time", y_title="level / activity")
        st.subheader("feedback dynamics")
        st.plotly_chart(fig, use_container_width=True)

        fig2 = phase_figure(df, x="ptch", y="gli", height=320, x_title="ptch", y_title="gli")
        st.subheader("phase view")
        st.plotly_chart(fig2, use_container_width=True)

else:
    with st.sidebar:
        st.divider()
        st.header("tissue simulation")
        nx = st.slider("nx", 51, 301, 151, 10)
        t_end = st.slider("t_end", 10.0, 500.0, 160.0, 10.0)
        dt = st.select_slider("dt", options=[0.01, 0.02, 0.05, 0.1, 0.2], value=0.05)

        st.divider()
        st.header("ligand transport")
        D = st.slider("diffusion D", 0.0, 0.2, 0.03, 0.005)
        k_decay = st.slider("decay k", 0.0, 0.2, 0.02, 0.005)
        source_hh = st.slider("source hh", 0.0, 1.0, 1.0, 0.01)
        source_width = st.slider("source width", 0.01, 0.25, 0.08, 0.01)

        st.divider()
        st.header("hhip sink")
        sink_strength = st.slider("sink strength", 0.0, 3.0, 0.8, 0.05)
        k_hhip_prod = st.slider("hhip prod", 0.0, 1.0, 0.25, 0.01)
        k_hhip_deg = st.slider("hhip deg", 0.0, 1.0, 0.10, 0.01)

        st.divider()
        st.header("gli relaxation")
        tau_gli = st.slider("tau_gli", 0.5, 50.0, 6.0, 0.5)

    sp = SpatialParams(
        D=D,
        k_decay=k_decay,
        source_hh=source_hh,
        source_width=source_width,
        sink_strength=sink_strength,
        k_hhip_prod=k_hhip_prod,
        k_hhip_deg=k_hhip_deg,
        core=core,
        tau_gli=tau_gli,
    )
    space_time, final = simulate_spatial(sp, nx=nx, t_end=t_end, dt=dt, ptch_level=ptch_level)

    st.subheader("final profiles")
    fig = profiles_figure(final, x="x", y_cols=["hh", "gli", "hhip"], height=340, x_title="position (x)", y_title="level / activity")
    st.plotly_chart(fig, use_container_width=True)

    st.subheader("space-time heatmaps")
    hh_piv = space_time.pivot(index="t", columns="x", values="hh")
    gli_piv = space_time.pivot(index="t", columns="x", values="gli")
    hhip_piv = space_time.pivot(index="t", columns="x", values="hhip")

    fig_hh = heatmap_figure(hh_piv.values, hh_piv.columns.values, hh_piv.index.values, colorscale="Viridis", height=300, x_title="x", y_title="t")
    st.caption("hh")
    st.plotly_chart(fig_hh, use_container_width=True)

    fig_g = heatmap_figure(gli_piv.values, gli_piv.columns.values, gli_piv.index.values, colorscale="Plasma", height=300, x_title="x", y_title="t")
    st.caption("gli")
    st.plotly_chart(fig_g, use_container_width=True)

    fig_h = heatmap_figure(hhip_piv.values, hhip_piv.columns.values, hhip_piv.index.values, colorscale="Cividis", height=300, x_title="x", y_title="t")
    st.caption("hhip")
    st.plotly_chart(fig_h, use_container_width=True)


with st.expander("model summary", expanded=False):
    st.write(
        "- steady state: hh reduces effective ptch repression; ptch represses smo; smo activates gli via hill\n"
        "- feedback ode: gli upregulates ptch (negative feedback) with tunable strength\n"
        "- extended ode: gli split into activator/repressor; gliA induces hhip which reduces effective hh\n"
        "- 1d tissue: hh diffuses/decays from a localized source; gli relaxes to local steady state; hhip adds a sink"
    )
