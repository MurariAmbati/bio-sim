from __future__ import annotations

import networkx as nx


def build_pathway_graph() -> nx.DiGraph:
    g = nx.DiGraph()

    nodes = [
        "tcr",
        "lck",
        "zap70/lat",
        "plcγ",
        "ip3",
        "ca²⁺",
        "ras",
        "raf",
        "mek",
        "erk",
        "pkcθ",
        "ikk",
        "nf-κb",
        "shp1",
        "cd28",
    ]
    g.add_nodes_from(nodes)

    # core flow
    g.add_edge("tcr", "zap70/lat")
    g.add_edge("lck", "zap70/lat")
    g.add_edge("cd28", "zap70/lat")

    g.add_edge("zap70/lat", "plcγ")
    g.add_edge("plcγ", "ip3")
    g.add_edge("ip3", "ca²⁺")

    g.add_edge("zap70/lat", "ras")
    g.add_edge("ras", "raf")
    g.add_edge("raf", "mek")
    g.add_edge("mek", "erk")

    # nfkb arm (via pkcθ→ikk)
    g.add_edge("zap70/lat", "pkcθ")
    g.add_edge("ca²⁺", "pkcθ")
    g.add_edge("pkcθ", "ikk")
    g.add_edge("ikk", "nf-κb")

    # negative regulator
    g.add_edge("shp1", "zap70/lat")

    return g
