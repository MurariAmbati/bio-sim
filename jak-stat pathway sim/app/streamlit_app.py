from __future__ import annotations

import streamlit as st

from jak_stat.model_ode import ODEParams, simulate_ode
from jak_stat.model_logical import LogicalParams, NODES, simulate_logical
from jak_stat.viz import logical_heatmap, timeseries_plot


st.set_page_config(page_title="jak_stat", layout="wide")

st.title("jak_stat")
st.caption("jak–stat pathway • cytokine and growth factor signaling")

with st.sidebar:
    st.header("model")
    model = st.selectbox("type", ["ode", "logical"], index=0)

    st.header("scenario")
    scenario = st.selectbox("case", ["baseline", "pulse", "socs knockout"], index=0)

    st.header("input")
    if model == "ode":
        cytokine_level = st.slider("cytokine level", 0.0, 5.0, 1.0, 0.1)
        t_end = st.slider("t end", 10.0, 200.0, 60.0, 5.0)
        n_steps = st.slider("n steps", 200, 3000, 900, 50)
        pulse_end = st.slider("pulse end", 1.0, float(t_end), 10.0, 1.0)

        st.header("parameters")
        kon = st.slider("kon", 0.0, 5.0, 1.0, 0.05)
        koff = st.slider("koff", 0.0, 2.0, 0.2, 0.02)
        k_p = st.slider("k_p", 0.0, 5.0, 2.0, 0.05)
        k_dp = st.slider("k_dp", 0.0, 2.0, 0.6, 0.02)
        k_socs_txn = st.slider("k_socs_txn", 0.0, 2.0, 0.8, 0.02)
        k_socs_deg = st.slider("k_socs_deg", 0.0, 2.0, 0.4, 0.02)
        k_socs_inhib = st.slider("k_socs_inhib", 0.0, 10.0, 2.0, 0.1)

    else:
        steps = st.slider("steps", 5, 100, 30, 1)
        cytokine_on = st.checkbox("cytokine on", value=True)
        pulse_end_step = st.slider("pulse end step", 0, int(steps), min(10, int(steps)), 1)


socs_knockout = scenario == "socs knockout"
pulse = scenario == "pulse"

if model == "ode":
    params = ODEParams(
        kon=kon,
        koff=koff,
        k_p=k_p,
        k_dp=k_dp,
        k_socs_txn=k_socs_txn,
        k_socs_deg=k_socs_deg,
        k_socs_inhib=k_socs_inhib,
    )

    df, meta = simulate_ode(
        params=params,
        t_end=t_end,
        n_steps=n_steps,
        cytokine_level=cytokine_level,
        pulse=pulse,
        pulse_end=pulse_end,
        socs_knockout=socs_knockout,
    )

    left, right = st.columns([2, 1])

    with left:
        st.plotly_chart(
            timeseries_plot(
                df,
                x="t",
                y=["lr", "jak", "pstat", "dimer", "nstat", "socs"],
                title="ode trajectories",
            ),
            use_container_width=True,
        )

    with right:
        st.subheader("notes")
        st.write("- lr: ligand–receptor complex")
        st.write("- jak: active jak")
        st.write("- nstat: nuclear stat proxy")
        st.write("- socs: negative feedback")
        st.subheader("meta")
        st.json(meta)
        st.subheader("data")
        st.dataframe(df.head(12), use_container_width=True)

else:
    params = LogicalParams(socs_knockout=socs_knockout)

    df, meta = simulate_logical(
        params=params,
        steps=steps,
        cytokine_on=cytokine_on,
        pulse=pulse,
        pulse_end_step=pulse_end_step,
    )

    left, right = st.columns([2, 1])

    with left:
        st.plotly_chart(
            logical_heatmap(df, nodes=NODES, title="logical state (0/1)"),
            use_container_width=True,
        )

    with right:
        st.subheader("meta")
        st.json(meta)
        st.subheader("data")
        st.dataframe(df.head(16), use_container_width=True)
