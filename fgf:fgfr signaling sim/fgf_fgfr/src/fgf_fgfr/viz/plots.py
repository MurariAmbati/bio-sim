from __future__ import annotations

from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
import plotly.graph_objects as go


def timeseries_plot(df: pd.DataFrame, y_cols: Sequence[str], title: str) -> go.Figure:
    fig = go.Figure()
    for c in y_cols:
        fig.add_trace(go.Scatter(x=df["t_s"], y=df[c], mode="lines", name=c))
    fig.update_layout(
        title=title,
        xaxis_title="Time (s)",
        yaxis_title="Level (a.u. / nM-eq)",
        legend_title="Species",
        height=420,
        margin=dict(l=20, r=20, t=50, b=20),
    )
    return fig


def outcome_bar(outcomes: Dict[str, float], title: str = "Outcome scores") -> go.Figure:
    keys = list(outcomes.keys())
    vals = [outcomes[k] for k in keys]
    fig = go.Figure(go.Bar(x=keys, y=vals))
    fig.update_layout(
        title=title,
        xaxis_title="Outcome",
        yaxis_title="Score (a.u.)",
        height=360,
        margin=dict(l=20, r=20, t=50, b=20),
    )
    return fig


def dose_response_plot(df: pd.DataFrame, y: str) -> go.Figure:
    fig = go.Figure(go.Scatter(x=df["dose_nM"], y=df[y], mode="lines+markers"))
    fig.update_layout(
        title=f"Dose-response: {y}",
        xaxis_title="FGF dose (nM)",
        yaxis_title="Outcome (a.u.)",
        height=380,
        margin=dict(l=20, r=20, t=50, b=20),
    )
    return fig


def sensitivity_bar(df: pd.DataFrame, top_n: int = 10) -> go.Figure:
    d = df.head(int(top_n)).iloc[::-1]
    fig = go.Figure(go.Bar(x=d["elasticity"], y=d["param"], orientation="h"))
    fig.update_layout(
        title=f"Local sensitivity (top {top_n})",
        xaxis_title="Elasticity (d log outcome / d log param)",
        yaxis_title="Parameter",
        height=420,
        margin=dict(l=20, r=20, t=50, b=20),
    )
    return fig
