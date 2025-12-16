from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Tuple

import networkx as nx


@dataclass(frozen=True)
class Node:
    id: str
    label: str
    kind: str  # ligand, receptor, adapter, kinase, outcome


def pathway_graph() -> nx.DiGraph:
    # High-level canonical FGF/FGFR signaling with major branches.
    g = nx.DiGraph()

    nodes = [
        Node("FGF", "FGF", "ligand"),
        Node("FGFR", "FGFR", "receptor"),
        Node("FGFR_p", "FGFR-P", "receptor"),
        Node("FRS2", "FRS2/GRB2/SOS", "adapter"),
        Node("RAS", "RAS", "kinase"),
        Node("RAF", "RAF", "kinase"),
        Node("MEK", "MEK", "kinase"),
        Node("ERK", "ERK", "kinase"),
        Node("PI3K", "PI3K", "kinase"),
        Node("AKT", "AKT", "kinase"),
        Node("PLCg", "PLCγ", "kinase"),
        Node("Ca", "Ca2+", "kinase"),
        Node("STAT", "STAT", "kinase"),
        Node("PROLIF", "Proliferation", "outcome"),
        Node("SURV", "Survival", "outcome"),
        Node("MIG", "Migration", "outcome"),
        Node("ANGIO", "Angiogenesis", "outcome"),
    ]

    for n in nodes:
        g.add_node(n.id, label=n.label, kind=n.kind)

    edges = [
        ("FGF", "FGFR"),
        ("FGFR", "FGFR_p"),
        ("FGFR_p", "FRS2"),
        ("FRS2", "RAS"),
        ("RAS", "RAF"),
        ("RAF", "MEK"),
        ("MEK", "ERK"),
        ("FRS2", "PI3K"),
        ("PI3K", "AKT"),
        ("FGFR_p", "PLCg"),
        ("PLCg", "Ca"),
        ("FGFR_p", "STAT"),
        ("ERK", "PROLIF"),
        ("AKT", "SURV"),
        ("Ca", "MIG"),
        ("ERK", "ANGIO"),
        ("Ca", "ANGIO"),
        ("AKT", "ANGIO"),
        ("ERK", "MIG"),
        ("AKT", "PROLIF"),
        ("STAT", "PROLIF"),
        ("STAT", "SURV"),
    ]

    g.add_edges_from(edges)
    return g


def node_palette(kind: str) -> str:
    # Return semantic palette names (actual colors applied by Plotly defaults).
    # We keep these as categories so the UI can map them consistently.
    return {
        "ligand": "ligand",
        "receptor": "receptor",
        "adapter": "adapter",
        "kinase": "kinase",
        "outcome": "outcome",
    }.get(kind, "other")
