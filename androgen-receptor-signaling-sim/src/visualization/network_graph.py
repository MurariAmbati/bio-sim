"""
network graph visualization for ar signaling pathway
"""

import numpy as np
import matplotlib.pyplot as plt
import networkx as nx
from typing import Optional, Dict, List, Tuple
import plotly.graph_objects as go
from plotly.subplots import make_subplots

from ..models.ar_pathway import ArPathwayModel


def plot_reaction_network(model: ArPathwayModel,
                          layout: str = 'hierarchical',
                          show_rates: bool = True,
                          node_size: int = 3000,
                          figsize: tuple = (16, 12),
                          save: Optional[str] = None,
                          show: bool = True,
                          interactive: bool = False) -> plt.Figure:
    """
    visualize reaction network as a directed graph
    
    args:
        model: ar pathway model
        layout: graph layout ('hierarchical', 'spring', 'circular', 'kamada_kawai')
        show_rates: whether to show rate constants on edges
        node_size: size of nodes
        figsize: figure size
        save: path to save figure
        show: whether to display
        interactive: create interactive plotly version
        
    returns:
        matplotlib figure or plotly figure
    """
    if interactive:
        return _plot_interactive_network(model, show_rates, save)
    
    # create directed graph
    G = nx.DiGraph()
    
    # define node categories
    species_categories = {
        'ligands': ['testosterone', 'dht'],
        'cytoplasmic_ar': ['ar_free_cytoplasm', 'ar_testosterone_cytoplasm', 
                          'ar_dht_cytoplasm'],
        'nuclear_ar': ['ar_free_nucleus', 'ar_testosterone_nucleus', 
                      'ar_dht_nucleus'],
        'dimers': ['ar_testosterone_dimer_nucleus', 'ar_dht_dimer_nucleus'],
        'dna_bound': ['ar_testosterone_dna', 'ar_dht_dna'],
        'coregulators': ['coactivator_free', 'corepressor_free',
                        'ar_dna_coactivator', 'ar_dna_corepressor'],
        'transcripts': ['mrna_psa', 'mrna_klk2', 'mrna_tmprss2', 
                       'mrna_nkx31', 'mrna_fkbp5'],
        'proteins': ['protein_psa'],
        'dna': ['are_sites_free']
    }
    
    # add nodes with categories
    node_colors = {
        'ligands': '#F18F01',
        'cytoplasmic_ar': '#2E86AB',
        'nuclear_ar': '#A23B72',
        'dimers': '#C73E1D',
        'dna_bound': '#6A4C93',
        'coregulators': '#52B788',
        'transcripts': '#E63946',
        'proteins': '#06D6A0',
        'dna': '#FFB703'
    }
    
    for category, species_list in species_categories.items():
        for species in species_list:
            G.add_node(species, category=category, color=node_colors[category])
    
    # add edges from reactions
    reactions = model.get_reaction_list()
    
    for rxn in reactions:
        for reactant in rxn['reactants']:
            for product in rxn['products']:
                if reactant != product:  # avoid self-loops for catalytic reactions
                    if G.has_edge(reactant, product):
                        G[reactant][product]['weight'] += 1
                    else:
                        G.add_edge(reactant, product, 
                                 weight=1, 
                                 rate=rxn['rate'],
                                 reaction=rxn['name'])
    
    # create layout
    if layout == 'hierarchical':
        pos = _hierarchical_layout(G, species_categories)
    elif layout == 'spring':
        pos = nx.spring_layout(G, k=2, iterations=50, seed=42)
    elif layout == 'circular':
        pos = nx.circular_layout(G)
    elif layout == 'kamada_kawai':
        pos = nx.kamada_kawai_layout(G)
    else:
        pos = nx.spring_layout(G, seed=42)
    
    # plot
    fig, ax = plt.subplots(figsize=figsize)
    
    # draw nodes by category
    for category, species_list in species_categories.items():
        nodes_in_cat = [n for n in species_list if n in G.nodes()]
        if nodes_in_cat:
            nx.draw_networkx_nodes(G, pos, 
                                 nodelist=nodes_in_cat,
                                 node_color=node_colors[category],
                                 node_size=node_size,
                                 alpha=0.9,
                                 ax=ax,
                                 label=category.replace('_', ' '))
    
    # draw edges
    nx.draw_networkx_edges(G, pos, 
                          edge_color='gray',
                          arrows=True,
                          arrowsize=20,
                          arrowstyle='->',
                          width=2,
                          alpha=0.6,
                          ax=ax,
                          connectionstyle='arc3,rad=0.1')
    
    # draw labels
    labels = {node: node.replace('_', '\n') for node in G.nodes()}
    nx.draw_networkx_labels(G, pos, labels, 
                           font_size=7,
                           font_weight='bold',
                           ax=ax)
    
    # add rate constants if requested
    if show_rates:
        edge_labels = {(u, v): d['rate'] 
                      for u, v, d in G.edges(data=True) if 'rate' in d}
        # only show subset to avoid clutter
        edge_labels = dict(list(edge_labels.items())[:20])
        nx.draw_networkx_edge_labels(G, pos, edge_labels, 
                                    font_size=6,
                                    ax=ax)
    
    ax.set_title('AR Signaling Reaction Network', 
                fontsize=18, fontweight='bold', pad=20)
    ax.legend(loc='upper left', bbox_to_anchor=(1.0, 1.0), 
             frameon=True, shadow=True)
    ax.axis('off')
    
    plt.tight_layout()
    
    if save:
        plt.savefig(save, dpi=300, bbox_inches='tight')
    
    if show:
        plt.show()
    else:
        plt.close()
    
    return fig


def _hierarchical_layout(G: nx.DiGraph, 
                        species_categories: Dict[str, List[str]]) -> Dict:
    """
    create hierarchical layout based on biological flow
    
    args:
        G: networkx graph
        species_categories: species grouped by category
        
    returns:
        position dictionary
    """
    pos = {}
    
    # define y-levels (flow from top to bottom)
    levels = {
        'ligands': 5,
        'cytoplasmic_ar': 4,
        'nuclear_ar': 3,
        'dimers': 2.5,
        'dna': 2,
        'dna_bound': 1.5,
        'coregulators': 1.5,
        'transcripts': 0.5,
        'proteins': 0
    }
    
    for category, species_list in species_categories.items():
        y = levels.get(category, 2)
        n_species = len([s for s in species_list if s in G.nodes()])
        
        for i, species in enumerate(species_list):
            if species in G.nodes():
                x = (i - n_species/2) * 2
                pos[species] = (x, y)
    
    return pos


def _plot_interactive_network(model: ArPathwayModel,
                              show_rates: bool = True,
                              save: Optional[str] = None) -> go.Figure:
    """
    create interactive plotly network visualization
    
    args:
        model: ar pathway model
        show_rates: show rate constants
        save: path to save html
        
    returns:
        plotly figure
    """
    # create networkx graph
    G = nx.DiGraph()
    
    # add reactions as edges
    reactions = model.get_reaction_list()
    
    for rxn in reactions:
        for reactant in rxn['reactants']:
            for product in rxn['products']:
                if reactant != product:
                    G.add_edge(reactant, product, 
                             rate=rxn['rate'],
                             reaction=rxn['name'])
    
    # layout
    pos = nx.spring_layout(G, k=3, iterations=50, seed=42)
    
    # create edge trace
    edge_x = []
    edge_y = []
    
    for edge in G.edges():
        x0, y0 = pos[edge[0]]
        x1, y1 = pos[edge[1]]
        edge_x.extend([x0, x1, None])
        edge_y.extend([y0, y1, None])
    
    edge_trace = go.Scatter(
        x=edge_x, y=edge_y,
        line=dict(width=2, color='#888'),
        hoverinfo='none',
        mode='lines',
        showlegend=False
    )
    
    # create node trace
    node_x = []
    node_y = []
    node_text = []
    node_color = []
    
    for node in G.nodes():
        x, y = pos[node]
        node_x.append(x)
        node_y.append(y)
        
        # node info
        neighbors = list(G.neighbors(node))
        node_info = f'{node}<br>→ {len(neighbors)} reactions'
        node_text.append(node_info)
        
        # color by type
        if 'mrna' in node:
            node_color.append('#E63946')
        elif 'protein' in node:
            node_color.append('#06D6A0')
        elif 'ar' in node and 'dna' in node:
            node_color.append('#6A4C93')
        elif 'ar' in node:
            node_color.append('#2E86AB')
        else:
            node_color.append('#F18F01')
    
    node_trace = go.Scatter(
        x=node_x, y=node_y,
        mode='markers+text',
        hoverinfo='text',
        text=[node.replace('_', '<br>') for node in G.nodes()],
        textposition='top center',
        textfont=dict(size=8),
        hovertext=node_text,
        marker=dict(
            showscale=False,
            color=node_color,
            size=30,
            line=dict(width=2, color='white')
        ),
        showlegend=False
    )
    
    # create figure
    fig = go.Figure(data=[edge_trace, node_trace],
                   layout=go.Layout(
                       title='<b>AR Signaling Reaction Network</b><br>(Interactive)',
                       titlefont=dict(size=20),
                       showlegend=False,
                       hovermode='closest',
                       margin=dict(b=20, l=5, r=5, t=60),
                       xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
                       yaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
                       plot_bgcolor='white',
                       height=800
                   ))
    
    if save:
        fig.write_html(save)
    
    return fig


def plot_pathway_flow(model: ArPathwayModel,
                     state_values: Dict[str, float],
                     figsize: tuple = (14, 10),
                     save: Optional[str] = None,
                     show: bool = True) -> plt.Figure:
    """
    visualize pathway with flux magnitudes
    
    args:
        model: ar pathway model
        state_values: current state values for flux calculation
        figsize: figure size
        save: path to save
        show: whether to display
        
    returns:
        matplotlib figure
    """
    fig, ax = plt.subplots(figsize=figsize)
    
    # create simplified flow diagram
    # (this would be a custom schematic showing key steps with widths proportional to flux)
    
    # for now, create a simple sankey-style diagram
    from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
    
    # define positions for key components
    positions = {
        'testosterone': (2, 8),
        'dht': (2, 6),
        'ar_cytoplasm': (5, 7),
        'ar_nucleus': (8, 7),
        'ar_dna': (11, 7),
        'mrna': (14, 7),
        'protein': (17, 7)
    }
    
    colors = {
        'testosterone': '#F18F01',
        'dht': '#C73E1D',
        'ar_cytoplasm': '#2E86AB',
        'ar_nucleus': '#A23B72',
        'ar_dna': '#6A4C93',
        'mrna': '#E63946',
        'protein': '#06D6A0'
    }
    
    # draw boxes
    for name, (x, y) in positions.items():
        value = state_values.get(name, 0)
        box = FancyBboxPatch((x-0.5, y-0.3), 1, 0.6,
                            boxstyle="round,pad=0.1",
                            facecolor=colors[name],
                            edgecolor='black',
                            linewidth=2,
                            alpha=0.7)
        ax.add_patch(box)
        ax.text(x, y, f'{name}\n{value:.0f}',
               ha='center', va='center',
               fontsize=9, fontweight='bold')
    
    # draw arrows (simplified)
    arrow_specs = [
        ('testosterone', 'ar_cytoplasm'),
        ('dht', 'ar_cytoplasm'),
        ('ar_cytoplasm', 'ar_nucleus'),
        ('ar_nucleus', 'ar_dna'),
        ('ar_dna', 'mrna'),
        ('mrna', 'protein')
    ]
    
    for start, end in arrow_specs:
        x1, y1 = positions[start]
        x2, y2 = positions[end]
        
        arrow = FancyArrowPatch((x1+0.5, y1), (x2-0.5, y2),
                              arrowstyle='->,head_width=0.4,head_length=0.8',
                              color='black',
                              linewidth=3,
                              alpha=0.6)
        ax.add_patch(arrow)
    
    ax.set_xlim(0, 20)
    ax.set_ylim(0, 10)
    ax.set_aspect('equal')
    ax.axis('off')
    ax.set_title('AR Signaling Pathway Flow', fontsize=18, fontweight='bold')
    
    plt.tight_layout()
    
    if save:
        plt.savefig(save, dpi=300, bbox_inches='tight')
    
    if show:
        plt.show()
    else:
        plt.close()
    
    return fig


def plot_network_statistics(model: ArPathwayModel,
                           figsize: tuple = (14, 6),
                           save: Optional[str] = None,
                           show: bool = True) -> plt.Figure:
    """
    analyze and visualize network topology statistics
    
    args:
        model: ar pathway model
        figsize: figure size
        save: path to save
        show: whether to display
        
    returns:
        matplotlib figure
    """
    # create graph
    G = nx.DiGraph()
    reactions = model.get_reaction_list()
    
    for rxn in reactions:
        for reactant in rxn['reactants']:
            for product in rxn['products']:
                if reactant != product:
                    G.add_edge(reactant, product)
    
    fig, axes = plt.subplots(1, 3, figsize=figsize)
    
    # degree distribution
    degrees = [G.degree(n) for n in G.nodes()]
    in_degrees = [G.in_degree(n) for n in G.nodes()]
    out_degrees = [G.out_degree(n) for n in G.nodes()]
    
    axes[0].hist([in_degrees, out_degrees], label=['in-degree', 'out-degree'],
                bins=range(max(degrees)+2), alpha=0.7)
    axes[0].set_xlabel('degree', fontweight='bold')
    axes[0].set_ylabel('count', fontweight='bold')
    axes[0].set_title('Degree Distribution', fontweight='bold')
    axes[0].legend()
    axes[0].grid(True, alpha=0.3)
    
    # node importance (betweenness centrality)
    centrality = nx.betweenness_centrality(G)
    nodes = list(centrality.keys())
    values = list(centrality.values())
    
    # top 10
    top_indices = np.argsort(values)[-10:]
    top_nodes = [nodes[i] for i in top_indices]
    top_values = [values[i] for i in top_indices]
    
    axes[1].barh(range(len(top_nodes)), top_values, color='#2E86AB')
    axes[1].set_yticks(range(len(top_nodes)))
    axes[1].set_yticklabels([n.replace('_', ' ') for n in top_nodes], fontsize=8)
    axes[1].set_xlabel('betweenness centrality', fontweight='bold')
    axes[1].set_title('Top 10 Hub Species', fontweight='bold')
    axes[1].grid(True, alpha=0.3, axis='x')
    
    # network metrics
    metrics = {
        'nodes': G.number_of_nodes(),
        'edges': G.number_of_edges(),
        'density': nx.density(G),
        'avg degree': np.mean(degrees),
    }
    
    # try to calculate clustering (if possible for directed graph)
    try:
        metrics['clustering'] = nx.average_clustering(G.to_undirected())
    except:
        metrics['clustering'] = 0.0
    
    axes[2].axis('off')
    metric_text = '\n'.join([f'{k}: {v:.3f}' for k, v in metrics.items()])
    axes[2].text(0.5, 0.5, metric_text,
                ha='center', va='center',
                fontsize=14, fontweight='bold',
                bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    axes[2].set_title('Network Metrics', fontweight='bold')
    
    plt.tight_layout()
    
    if save:
        plt.savefig(save, dpi=300, bbox_inches='tight')
    
    if show:
        plt.show()
    else:
        plt.close()
    
    return fig
