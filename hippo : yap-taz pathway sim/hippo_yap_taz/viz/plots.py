from __future__ import annotations

from dataclasses import asdict
from typing import Dict, Iterable, List, Tuple

import networkx as nx
import numpy as np
import plotly.graph_objects as go

from hippo_yap_taz.sim.model import ModelParams
from hippo_yap_taz.sim.simulate import SimulationResult, run_to_steady_state


def _series(df, name: str) -> np.ndarray:
    return np.asarray(df[name].to_numpy(dtype=float))


def make_timeseries_figure(res: SimulationResult) -> go.Figure:
    df = res.df

    fig = go.Figure()
    fig.add_trace(go.Scatter(x=df["t"], y=_series(df, "mst"), name="mst", mode="lines"))
    fig.add_trace(go.Scatter(x=df["t"], y=_series(df, "lats"), name="lats", mode="lines"))
    fig.add_trace(go.Scatter(x=df["t"], y=_series(df, "actin"), name="actin", mode="lines"))
    fig.add_trace(go.Scatter(x=df["t"], y=_series(df, "yap_c"), name="yap_c", mode="lines"))
    fig.add_trace(go.Scatter(x=df["t"], y=_series(df, "yap_p"), name="yap_p", mode="lines"))
    fig.add_trace(go.Scatter(x=df["t"], y=_series(df, "yap_n"), name="yap_n", mode="lines"))
    fig.add_trace(go.Scatter(x=df["t"], y=_series(df, "gene"), name="target_gene", mode="lines"))

    fig.update_layout(
        height=420,
        xaxis_title="time",
        yaxis_title="state",
        legend_title="species",
        margin=dict(l=40, r=20, t=20, b=40),
    )

    subtitle = (
        f"stiffness={res.settings.stiffness:.2f} · density={res.settings.density:.2f} · "
        f"nuclear_yap_fraction={res.steady['yap_nuclear_fraction']:.3f}"
    )
    fig.update_layout(title=dict(text=subtitle, x=0.01, xanchor="left"))

    return fig


def make_phase_plot(res: SimulationResult) -> go.Figure:
    df = res.df
    fig = go.Figure()

    fig.add_trace(
        go.Scatter(
            x=_series(df, "lats"),
            y=_series(df, "yap_n"),
            mode="lines",
            name="trajectory",
        )
    )

    fig.add_trace(
        go.Scatter(
            x=[float(df["lats"].iloc[-1])],
            y=[float(df["yap_n"].iloc[-1])],
            mode="markers",
            marker=dict(size=10),
            name="steady",
        )
    )

    fig.update_layout(
        height=380,
        xaxis_title="lats activity",
        yaxis_title="nuclear yap",
        margin=dict(l=40, r=20, t=20, b=40),
    )
    return fig


def _build_graph() -> nx.DiGraph:
    g = nx.DiGraph()
    g.add_nodes_from(["stiffness", "density", "actin", "mst", "lats", "yap_p", "yap_n", "gene"])

    # edges represent activation (positive) or inhibition (negative)
    g.add_edge("stiffness", "actin", sign=+1)
    g.add_edge("mst", "actin", sign=-1)

    g.add_edge("density", "mst", sign=+1)
    g.add_edge("actin", "mst", sign=-1)

    g.add_edge("mst", "lats", sign=+1)
    g.add_edge("actin", "lats", sign=-1)

    g.add_edge("lats", "yap_p", sign=+1)
    g.add_edge("yap_p", "yap_n", sign=-1)
    g.add_edge("actin", "yap_n", sign=+1)

    g.add_edge("yap_n", "gene", sign=+1)

    return g


def _node_value(res: SimulationResult, node: str) -> float:
    if node in {"stiffness", "density"}:
        return float(getattr(res.settings, node))
    if node == "yap_p":
        return float(res.steady["yap_p"])
    return float(res.steady.get(node, 0.0))


def make_network_figure(res: SimulationResult) -> go.Figure:
    g = _build_graph()

    pos = nx.spring_layout(g, seed=2)

    node_x = []
    node_y = []
    node_text = []
    node_color = []

    for n in g.nodes:
        x, y = pos[n]
        node_x.append(x)
        node_y.append(y)
        val = _node_value(res, n)
        node_color.append(val)
        node_text.append(f"{n}: {val:.3g}")

    edge_x = []
    edge_y = []
    edge_color = []
    for u, v, data in g.edges(data=True):
        x0, y0 = pos[u]
        x1, y1 = pos[v]
        edge_x += [x0, x1, None]
        edge_y += [y0, y1, None]
        edge_color.append("#444" if data.get("sign", 1) > 0 else "#999")

    # edges as a single trace (uniform style). plotly doesn't support per-segment color easily,
    # so we encode sign via dash.
    edge_trace = go.Scatter(
        x=edge_x,
        y=edge_y,
        line=dict(width=2, color="#666", dash="solid"),
        hoverinfo="none",
        mode="lines",
        name="edges",
    )

    # mark inhibitory edges with a second dashed layer
    inh_edge_x = []
    inh_edge_y = []
    for u, v, data in g.edges(data=True):
        if data.get("sign", 1) < 0:
            x0, y0 = pos[u]
            x1, y1 = pos[v]
            inh_edge_x += [x0, x1, None]
            inh_edge_y += [y0, y1, None]

    inh_edge_trace = go.Scatter(
        x=inh_edge_x,
        y=inh_edge_y,
        line=dict(width=3, color="#999", dash="dash"),
        hoverinfo="none",
        mode="lines",
        name="inhibition",
    )

    node_trace = go.Scatter(
        x=node_x,
        y=node_y,
        mode="markers+text",
        text=list(g.nodes),
        textposition="bottom center",
        hovertext=node_text,
        hoverinfo="text",
        marker=dict(
            showscale=True,
            colorscale="Viridis",
            size=18,
            color=node_color,
            colorbar=dict(title="state"),
            line=dict(width=1, color="#222"),
        ),
        name="nodes",
    )

    fig = go.Figure(data=[edge_trace, inh_edge_trace, node_trace])
    fig.update_layout(
        height=420,
        margin=dict(l=10, r=10, t=10, b=10),
        xaxis=dict(visible=False),
        yaxis=dict(visible=False),
        showlegend=False,
    )
    return fig


def make_steady_state_heatmap(*, base_settings, grid_n: int = 25) -> go.Figure:
    # base_settings is SimulationSettings
    grid_n = int(max(5, grid_n))

    stiffness_vals = np.linspace(0.0, 1.0, grid_n)
    density_vals = np.linspace(0.0, 1.0, grid_n)

    z = np.zeros((grid_n, grid_n), dtype=float)

    # keep this moderately fast; rely on short integration + caching higher up if desired
    for i, dens in enumerate(density_vals):
        for j, stiff in enumerate(stiffness_vals):
            steady = run_to_steady_state(
                stiffness=float(stiff),
                density=float(dens),
                t_end=max(400.0, float(base_settings.t_end)),
                n_points=1400,
                method=str(base_settings.method),
            )
            z[i, j] = float(steady["yap_nuclear_fraction"])

    fig = go.Figure(
        data=
        [
            go.Heatmap(
                z=z,
                x=stiffness_vals,
                y=density_vals,
                colorscale="Viridis",
                colorbar=dict(title="nuclear yap fraction"),
            )
        ]
    )

    fig.update_layout(
        height=480,
        xaxis_title="substrate stiffness",
        yaxis_title="cell density / contact inhibition",
        margin=dict(l=50, r=20, t=20, b=50),
    )

    return fig


def make_sensitivity_figure(res: SimulationResult) -> go.Figure:
    # local (one-at-a-time) sensitivity around current params
    p0 = res.params
    p0_dict = asdict(p0)

    # only a curated subset of params
    keys = [
        "k_m_ci",
        "k_m_act",
        "k_l_act",
        "k_phos",
        "k_dephos",
        "k_imp",
        "k_exp",
        "actin_inhib_mst",
        "actin_inhib_lats",
        "actin_promote_import",
        "k_txn",
    ]

    base = float(res.steady["yap_nuclear_fraction"])
    rel = []

    for k in keys:
        v = float(p0_dict[k])
        if v == 0:
            rel.append(0.0)
            continue

        p_plus = ModelParams(**{**p0_dict, k: v * 1.10})
        p_minus = ModelParams(**{**p0_dict, k: v * 0.90})

        s_plus = run_to_steady_state(
            stiffness=float(res.settings.stiffness),
            density=float(res.settings.density),
            t_end=max(500.0, float(res.settings.t_end)),
            n_points=1600,
            method=str(res.settings.method),
            params=p_plus,
        )
        s_minus = run_to_steady_state(
            stiffness=float(res.settings.stiffness),
            density=float(res.settings.density),
            t_end=max(500.0, float(res.settings.t_end)),
            n_points=1600,
            method=str(res.settings.method),
            params=p_minus,
        )

        # symmetric relative sensitivity: (f+ - f-)/(2*0.1*v) * v / f = (f+ - f-)/(0.2*f)
        rel_sens = (float(s_plus["yap_nuclear_fraction"]) - float(s_minus["yap_nuclear_fraction"])) / max(1e-12, 0.2 * base)
        rel.append(rel_sens)

    fig = go.Figure(
        data=[
            go.Bar(
                x=keys,
                y=rel,
            )
        ]
    )
    fig.update_layout(
        height=380,
        xaxis_title="parameter",
        yaxis_title="relative sensitivity of nuclear yap fraction",
        margin=dict(l=40, r=20, t=20, b=80),
    )
    return fig
