from __future__ import annotations

import math

import networkx as nx
import pandas as pd
import plotly.graph_objects as go


def plot_network(g: nx.DiGraph) -> go.Figure:
    # deterministic layout
    pos = nx.spring_layout(g, seed=7, k=0.8)

    edge_x = []
    edge_y = []
    for u, v in g.edges():
        x0, y0 = pos[u]
        x1, y1 = pos[v]
        edge_x += [x0, x1, None]
        edge_y += [y0, y1, None]

    edge_trace = go.Scatter(
        x=edge_x,
        y=edge_y,
        line=dict(width=1),
        hoverinfo="none",
        mode="lines",
    )

    node_x = []
    node_y = []
    labels = []
    for n in g.nodes():
        x, y = pos[n]
        node_x.append(x)
        node_y.append(y)
        labels.append(n)

    node_trace = go.Scatter(
        x=node_x,
        y=node_y,
        mode="markers+text",
        text=labels,
        textposition="top center",
        hoverinfo="text",
        marker=dict(size=12),
    )

    fig = go.Figure(data=[edge_trace, node_trace])
    fig.update_layout(
        margin=dict(l=10, r=10, t=10, b=10),
        xaxis=dict(visible=False),
        yaxis=dict(visible=False),
        showlegend=False,
    )
    return fig


def plot_timeseries(df: pd.DataFrame) -> go.Figure:
    fig = go.Figure()

    for col in ["nfkb", "erk", "ca"]:
        if col in df.columns:
            fig.add_trace(go.Scatter(x=df["t_min"], y=df[col], mode="lines", name=col))

    fig.update_layout(
        xaxis_title="time (min)",
        yaxis_title="activity (0–1)",
        margin=dict(l=10, r=10, t=10, b=10),
        legend=dict(orientation="h"),
    )
    fig.update_yaxes(range=[-0.05, 1.05])
    return fig
