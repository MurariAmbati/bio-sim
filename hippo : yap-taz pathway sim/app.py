import streamlit as st

from hippo_yap_taz.sim.simulate import run_simulation
from hippo_yap_taz.viz.plots import (
    make_network_figure,
    make_phase_plot,
    make_sensitivity_figure,
    make_steady_state_heatmap,
    make_timeseries_figure,
)


st.set_page_config(page_title="hippo / yap–taz pathway sim", layout="wide")

st.title("hippo / yap–taz pathway")
st.caption("organ size control · proliferation vs apoptosis · mechano-sensing")

with st.sidebar:
    st.header("controls")
    scenario = st.selectbox(
        "scenario",
        [
            "baseline",
            "soft substrate (low stiffness)",
            "stiff substrate (high stiffness)",
            "high density (contact inhibition)",
        ],
        index=0,
    )

    stiffness = st.slider("substrate stiffness (0=soft, 1=stiff)", 0.0, 1.0, 0.7, 0.01)
    density = st.slider("cell density / contact inhibition (0..1)", 0.0, 1.0, 0.2, 0.01)

    st.divider()
    t_end = st.number_input("simulation time", min_value=10.0, max_value=2_000.0, value=300.0, step=10.0)
    n_points = st.number_input("output points", min_value=200, max_value=10_000, value=1200, step=100)
    method = st.selectbox("ode method", ["BDF", "LSODA", "RK45"], index=0)

    st.divider()
    grid_n = st.slider("heatmap grid size", 10, 60, 25, 1)

# scenario presets (light touch; sliders still apply)
if scenario == "soft substrate (low stiffness)":
    stiffness = min(stiffness, 0.25)
elif scenario == "stiff substrate (high stiffness)":
    stiffness = max(stiffness, 0.85)
elif scenario == "high density (contact inhibition)":
    density = max(density, 0.75)

result = run_simulation(
    stiffness=float(stiffness),
    density=float(density),
    t_end=float(t_end),
    n_points=int(n_points),
    method=str(method),
)

col_a, col_b = st.columns([2, 1])
with col_a:
    st.subheader("time series")
    st.plotly_chart(make_timeseries_figure(result), use_container_width=True)

with col_b:
    st.subheader("network state")
    st.plotly_chart(make_network_figure(result), use_container_width=True)

col_c, col_d = st.columns(2)
with col_c:
    st.subheader("phase portrait")
    st.plotly_chart(make_phase_plot(result), use_container_width=True)

with col_d:
    st.subheader("local sensitivity (steady nuclear yap)")
    st.plotly_chart(make_sensitivity_figure(result), use_container_width=True)

st.subheader("steady-state map")
st.plotly_chart(
    make_steady_state_heatmap(
        base_settings=result.settings,
        grid_n=int(grid_n),
    ),
    use_container_width=True,
)
