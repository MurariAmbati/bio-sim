from __future__ import annotations

import numpy as np
import pandas as pd
import plotly.graph_objects as go
import streamlit as st

from wnt_model import WntParams, cluster_endpoints, endpoints_from_multiple_inits, simulate, wnt_ode


st.set_page_config(page_title="wnt_beta_catenin", layout="wide")

st.title("wnt / β‑catenin pathway · ode")
st.caption("development · stem cell regulation · cancer · slug: wnt_beta_catenin")


with st.sidebar:
    st.subheader("inputs")
    wnt = st.slider("wnt input", min_value=0.0, max_value=1.0, value=0.35, step=0.01)

    st.subheader("initial state")
    beta0 = st.number_input("beta(0)", min_value=0.0, value=0.2, step=0.05)
    dc0 = st.number_input("dc(0)", min_value=0.0, value=1.0, step=0.05)

    st.subheader("time")
    t_end = st.number_input("t_end", min_value=5.0, value=80.0, step=5.0)
    n_points = st.slider("points", min_value=200, max_value=5000, value=1200, step=100)

    st.subheader("parameters")
    col1, col2 = st.columns(2)
    with col1:
        k_prod = st.number_input("k_prod", min_value=0.0, value=0.25, step=0.05)
        k_stab = st.number_input("k_stab", min_value=0.0, value=0.6, step=0.05)
        k_base = st.number_input("k_base", min_value=0.0, value=0.05, step=0.01)
        k_deg = st.number_input("k_deg", min_value=0.0, value=0.9, step=0.05)
    with col2:
        k_syn = st.number_input("k_syn", min_value=0.0, value=0.8, step=0.05)
        k_turn = st.number_input("k_turn", min_value=0.0, value=0.25, step=0.05)
        k_inhib = st.number_input("k_inhib", min_value=0.0, value=1.2, step=0.05)
        K = st.number_input("K", min_value=1e-6, value=0.7, step=0.05)

    n = st.slider("hill n", min_value=1.0, max_value=6.0, value=3.0, step=0.5)

    st.subheader("advanced")
    show_vector_field = st.checkbox("phase-plane vector field", value=True)
    do_scan = st.checkbox("steady-state scan", value=True)

    scan_steps = st.slider("scan steps", min_value=15, max_value=120, value=55, step=5)
    init_grid = st.slider("scan init grid", min_value=2, max_value=8, value=4, step=1)
    cluster_tol = st.slider("cluster tol", min_value=0.005, max_value=0.2, value=0.02, step=0.005)


p = WntParams(
    k_prod=k_prod,
    k_stab=k_stab,
    k_base=k_base,
    k_deg=k_deg,
    k_syn=k_syn,
    k_turn=k_turn,
    k_inhib=k_inhib,
    K=K,
    n=float(n),
)


@st.cache_data(show_spinner=False)
def _simulate_cached(p_dict: dict, wnt: float, beta0: float, dc0: float, t_end: float, n_points: int):
    p_obj = WntParams(**p_dict)
    t, y = simulate(p=p_obj, wnt=wnt, y0=(beta0, dc0), t_end=t_end, n_points=n_points)
    return t, y


@st.cache_data(show_spinner=False)
def _scan_cached(p_dict: dict, t_end: float, steps: int, init_grid: int, cluster_tol: float):
    p_obj = WntParams(**p_dict)

    wnts = np.linspace(0.0, 1.0, int(steps))

    beta_inits = np.linspace(0.0, 3.0, int(init_grid))
    dc_inits = np.linspace(0.0, 3.0, int(init_grid))

    rows = []
    for u in wnts:
        endpoints = endpoints_from_multiple_inits(p_obj, float(u), t_end=float(t_end), beta_inits=beta_inits, dc_inits=dc_inits)
        centers = cluster_endpoints(endpoints, tol=float(cluster_tol))
        for c in centers:
            rows.append({"wnt": float(u), "beta_ss": float(c[0]), "dc_ss": float(c[1]), "n_ss": int(centers.shape[0])})

    return pd.DataFrame(rows)


t, y = _simulate_cached(p.__dict__, float(wnt), float(beta0), float(dc0), float(t_end), int(n_points))

beta = y[0]
dc = y[1]


df = pd.DataFrame({"t": t, "beta": beta, "dc": dc})

col_left, col_right = st.columns([1.2, 1])

with col_left:
    st.subheader("time series")
    fig = go.Figure()
    fig.add_trace(go.Scatter(x=df["t"], y=df["beta"], mode="lines", name="beta"))
    fig.add_trace(go.Scatter(x=df["t"], y=df["dc"], mode="lines", name="dc", yaxis="y2"))

    fig.update_layout(
        height=420,
        margin=dict(l=10, r=10, t=25, b=10),
        xaxis_title="time",
        yaxis=dict(title="beta"),
        yaxis2=dict(title="dc", overlaying="y", side="right"),
        legend=dict(orientation="h"),
    )
    st.plotly_chart(fig, use_container_width=True)

    st.download_button(
        "download simulation csv",
        data=df.to_csv(index=False).encode("utf-8"),
        file_name="wnt_beta_catenin_sim.csv",
        mime="text/csv",
    )

with col_right:
    st.subheader("phase plane")

    fig2 = go.Figure()

    if show_vector_field:
        b_min = 0.0
        b_max = max(3.0, float(np.max(beta)) * 1.05)
        d_min = 0.0
        d_max = max(3.0, float(np.max(dc)) * 1.05)

        bb = np.linspace(b_min, b_max, 18)
        dd = np.linspace(d_min, d_max, 18)
        BB, DD = np.meshgrid(bb, dd)

        U = np.zeros_like(BB)
        V = np.zeros_like(DD)
        for i in range(BB.shape[0]):
            for j in range(BB.shape[1]):
                dv = wnt_ode(0.0, np.array([BB[i, j], DD[i, j]]), p, float(wnt))
                U[i, j] = dv[0]
                V[i, j] = dv[1]

        mag = np.sqrt(U**2 + V**2) + 1e-12
        U2 = U / mag
        V2 = V / mag

        fig2.add_trace(
            go.Cone(
                x=BB.flatten(),
                y=DD.flatten(),
                z=np.zeros(BB.size),
                u=U2.flatten(),
                v=V2.flatten(),
                w=np.zeros(BB.size),
                sizemode="absolute",
                sizeref=0.35,
                anchor="tail",
                showscale=False,
                name="vector field",
            )
        )

    fig2.add_trace(go.Scatter(x=beta, y=dc, mode="lines", name="trajectory"))
    fig2.add_trace(go.Scatter(x=[beta[0]], y=[dc[0]], mode="markers", name="start"))
    fig2.add_trace(go.Scatter(x=[beta[-1]], y=[dc[-1]], mode="markers", name="end"))

    fig2.update_layout(
        height=420,
        margin=dict(l=10, r=10, t=25, b=10),
        scene=dict(
            xaxis_title="beta",
            yaxis_title="dc",
            zaxis=dict(visible=False),
        ),
        showlegend=True,
    )

    st.plotly_chart(fig2, use_container_width=True)


st.subheader("steady-state scan (wnt → endpoints)")
if do_scan:
    with st.spinner("scanning…"):
        sdf = _scan_cached(p.__dict__, float(t_end), int(scan_steps), int(init_grid), float(cluster_tol))

    if sdf.empty:
        st.info("no scan data")
    else:
        fig3 = go.Figure()
        fig3.add_trace(go.Scatter(x=sdf["wnt"], y=sdf["beta_ss"], mode="markers", name="beta_ss"))
        fig3.update_layout(
            height=360,
            margin=dict(l=10, r=10, t=25, b=10),
            xaxis_title="wnt",
            yaxis_title="beta steady state",
        )
        st.plotly_chart(fig3, use_container_width=True)

        fig4 = go.Figure()
        fig4.add_trace(go.Scatter(x=sdf["wnt"], y=sdf["n_ss"], mode="lines+markers", name="# distinct endpoints"))
        fig4.update_layout(
            height=220,
            margin=dict(l=10, r=10, t=25, b=10),
            xaxis_title="wnt",
            yaxis_title="count",
        )
        st.plotly_chart(fig4, use_container_width=True)

        st.download_button(
            "download scan csv",
            data=sdf.to_csv(index=False).encode("utf-8"),
            file_name="wnt_beta_catenin_scan.csv",
            mime="text/csv",
        )
else:
    st.caption("scan disabled")
