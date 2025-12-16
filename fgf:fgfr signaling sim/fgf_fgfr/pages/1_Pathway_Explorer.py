from __future__ import annotations

import sys
from pathlib import Path
import math
from typing import Dict, Tuple

import networkx as nx
import numpy as np
import plotly.graph_objects as go
import streamlit as st
from streamlit_plotly_events import plotly_events

APP_DIR = Path(__file__).resolve().parents[1]
SRC_DIR = APP_DIR / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from fgf_fgfr.viz.pathway_graph import pathway_graph


st.set_page_config(page_title="Pathway Explorer", layout="wide")

st.title("FGF/FGFR pathway explorer")

st.write(
    "Tap/click nodes to inspect what they represent. This is a canonical pathway map used by the simulator (lumped model)."
)


g = pathway_graph()

# Layout (spring) in 2D
pos = nx.spring_layout(g, seed=2, k=0.9)

node_ids = list(g.nodes())
xs = [pos[n][0] for n in node_ids]
ys = [pos[n][1] for n in node_ids]
labels = [g.nodes[n].get("label", n) for n in node_ids]

# edges as segments
edge_x = []
edge_y = []
for u, v in g.edges():
    edge_x += [pos[u][0], pos[v][0], None]
    edge_y += [pos[u][1], pos[v][1], None]

edge_trace = go.Scatter(
    x=edge_x,
    y=edge_y,
    mode="lines",
    line=dict(width=1),
    hoverinfo="none",
)

node_trace = go.Scatter(
    x=xs,
    y=ys,
    mode="markers+text",
    text=labels,
    textposition="top center",
    marker=dict(size=18),
    customdata=node_ids,
    hovertemplate="%{text}<extra></extra>",
)

fig = go.Figure(data=[edge_trace, node_trace])
fig.update_layout(
    height=650,
    margin=dict(l=20, r=20, t=20, b=20),
    xaxis=dict(visible=False),
    yaxis=dict(visible=False),
)

selected = plotly_events(fig, click_event=True, hover_event=False, select_event=False, override_height=650)

if selected:
    point = selected[0]
    curve = point.get("curveNumber", 1)
    idx = point.get("pointIndex")
    if idx is not None:
        node_id = node_ids[int(idx)]
        st.subheader(f"Selected: {g.nodes[node_id].get('label', node_id)}")
        st.write({"id": node_id, **g.nodes[node_id]})

        preds = list(g.predecessors(node_id))
        succs = list(g.successors(node_id))
        c1, c2 = st.columns(2)
        with c1:
            st.caption("Upstream")
            st.write([g.nodes[p].get("label", p) for p in preds] or "(none)")
        with c2:
            st.caption("Downstream")
            st.write([g.nodes[s].get("label", s) for s in succs] or "(none)")
else:
    st.info("Tap a node to inspect it.")
