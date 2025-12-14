from __future__ import annotations

import pandas as pd
import plotly.graph_objects as go


def timeseries_figure(
    df: pd.DataFrame,
    x: str,
    y_cols: list[str],
    height: int,
    x_title: str,
    y_title: str,
) -> go.Figure:
    fig = go.Figure()
    for col in y_cols:
        fig.add_trace(go.Scatter(x=df[x], y=df[col], mode="lines", name=col))
    fig.update_layout(
        height=height,
        margin=dict(l=10, r=10, t=10, b=10),
        xaxis_title=x_title,
        yaxis_title=y_title,
        legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="left", x=0.0),
    )
    return fig


def phase_figure(df: pd.DataFrame, x: str, y: str, height: int, x_title: str, y_title: str) -> go.Figure:
    fig = go.Figure()
    fig.add_trace(go.Scatter(x=df[x], y=df[y], mode="lines", name="trajectory"))
    fig.update_layout(
        height=height,
        margin=dict(l=10, r=10, t=10, b=10),
        xaxis_title=x_title,
        yaxis_title=y_title,
        showlegend=False,
    )
    return fig


def profiles_figure(df: pd.DataFrame, x: str, y_cols: list[str], height: int, x_title: str, y_title: str) -> go.Figure:
    fig = go.Figure()
    for col in y_cols:
        fig.add_trace(go.Scatter(x=df[x], y=df[col], mode="lines", name=col))
    fig.update_layout(
        height=height,
        margin=dict(l=10, r=10, t=10, b=10),
        xaxis_title=x_title,
        yaxis_title=y_title,
        legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="left", x=0.0),
    )
    return fig


def heatmap_figure(z, x, y, colorscale: str, height: int, x_title: str, y_title: str) -> go.Figure:
    fig = go.Figure(data=go.Heatmap(z=z, x=x, y=y, colorscale=colorscale))
    fig.update_layout(height=height, margin=dict(l=10, r=10, t=10, b=10), xaxis_title=x_title, yaxis_title=y_title)
    return fig
