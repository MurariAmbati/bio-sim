from __future__ import annotations

from typing import Iterable, Sequence

import pandas as pd
import plotly.express as px


def timeseries_plot(df: pd.DataFrame, *, x: str, y: Sequence[str], title: str):
    long_df = df.melt(id_vars=[x], value_vars=list(y), var_name="species", value_name="value")
    fig = px.line(long_df, x=x, y="value", color="species", title=title)
    fig.update_layout(legend_title_text="")
    return fig


def logical_heatmap(df: pd.DataFrame, *, nodes: Iterable[str], title: str):
    nodes = list(nodes)
    plot_df = df[["step", *nodes]].set_index("step")
    fig = px.imshow(
        plot_df.T,
        aspect="auto",
        color_continuous_scale="Blues",
        origin="lower",
        title=title,
        labels={"x": "step", "y": "node", "color": "state"},
    )
    return fig
