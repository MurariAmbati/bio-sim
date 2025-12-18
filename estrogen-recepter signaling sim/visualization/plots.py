"""
visualization utilities for er signaling
"""

import numpy as np
import pandas as pd
import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots
from typing import Dict, List, Optional


class PathwayVisualizer:
    """
    create interactive visualizations of er signaling pathways
    """
    
    @staticmethod
    def plot_receptor_dynamics(df: pd.DataFrame) -> go.Figure:
        """plot er-alpha and er-beta dynamics"""
        fig = make_subplots(
            rows=2, cols=2,
            subplot_titles=(
                'receptor concentration',
                'transcriptional activity',
                'receptor state distribution',
                'phosphorylation status'
            )
        )
        
        # concentration
        fig.add_trace(
            go.Scatter(x=df['time_hours'], y=df['er_alpha_concentration'],
                      name='er-α', line=dict(color='#e74c3c', width=2)),
            row=1, col=1
        )
        fig.add_trace(
            go.Scatter(x=df['time_hours'], y=df['er_beta_concentration'],
                      name='er-β', line=dict(color='#3498db', width=2)),
            row=1, col=1
        )
        
        # transcriptional activity
        fig.add_trace(
            go.Scatter(x=df['time_hours'], y=df['er_alpha_transcriptional_activity'],
                      name='er-α activity', line=dict(color='#e74c3c', width=2),
                      showlegend=False),
            row=1, col=2
        )
        
        # state distribution (final timepoint)
        states = ['unbound', 'ligand_bound', 'dimerized', 'dna_bound']
        # approximate state distribution
        final_activity = df['er_alpha_transcriptional_activity'].iloc[-1]
        state_values = [
            max(0, 1 - final_activity),
            0.3 * final_activity,
            0.3 * final_activity,
            0.4 * final_activity
        ]
        
        fig.add_trace(
            go.Bar(x=states, y=state_values, marker_color='#e74c3c',
                  showlegend=False),
            row=2, col=1
        )
        
        # phosphorylation
        if 'er_alpha_phosphorylation' in df.columns:
            fig.add_trace(
                go.Scatter(x=df['time_hours'], y=df['er_alpha_phosphorylation'],
                          name='phosphorylation', line=dict(color='#9b59b6', width=2),
                          showlegend=False),
                row=2, col=2
            )
        
        fig.update_xaxes(title_text="time (hours)", row=2, col=1)
        fig.update_xaxes(title_text="time (hours)", row=2, col=2)
        fig.update_xaxes(title_text="time (hours)", row=1, col=1)
        fig.update_xaxes(title_text="time (hours)", row=1, col=2)
        
        fig.update_yaxes(title_text="concentration (nm)", row=1, col=1)
        fig.update_yaxes(title_text="activity", row=1, col=2)
        fig.update_yaxes(title_text="fraction", row=2, col=1)
        fig.update_yaxes(title_text="phosphorylation", row=2, col=2)
        
        fig.update_layout(height=800, title_text="estrogen receptor dynamics")
        
        return fig
    
    @staticmethod
    def plot_target_genes(df: pd.DataFrame) -> go.Figure:
        """plot target gene expression"""
        genes = ['pr', 'cyclin_d1', 'bcl2', 'tff1', 'vegf', 'c_myc']
        
        fig = make_subplots(
            rows=3, cols=2,
            subplot_titles=[gene.upper() for gene in genes]
        )
        
        positions = [(1, 1), (1, 2), (2, 1), (2, 2), (3, 1), (3, 2)]
        colors = ['#e74c3c', '#3498db', '#2ecc71', '#f39c12', '#9b59b6', '#1abc9c']
        
        for (gene, pos, color) in zip(genes, positions, colors):
            mrna_col = f'gene_{gene}_mrna'
            protein_col = f'gene_{gene}_protein'
            
            if mrna_col in df.columns:
                fig.add_trace(
                    go.Scatter(x=df['time_hours'], y=df[mrna_col],
                              name='mrna', line=dict(color=color, dash='dot'),
                              legendgroup=gene, showlegend=(pos == (1, 1))),
                    row=pos[0], col=pos[1]
                )
            
            if protein_col in df.columns:
                fig.add_trace(
                    go.Scatter(x=df['time_hours'], y=df[protein_col],
                              name='protein', line=dict(color=color, width=2),
                              legendgroup=gene, showlegend=(pos == (1, 1))),
                    row=pos[0], col=pos[1]
                )
        
        fig.update_xaxes(title_text="time (hours)")
        fig.update_yaxes(title_text="level (au)")
        fig.update_layout(height=1000, title_text="er target gene expression")
        
        return fig
    
    @staticmethod
    def plot_signaling_cascades(df: pd.DataFrame) -> go.Figure:
        """plot non-genomic signaling pathways"""
        fig = make_subplots(
            rows=2, cols=2,
            subplot_titles=(
                'mapk/erk pathway',
                'pi3k/akt pathway',
                'calcium signaling',
                'nitric oxide production'
            )
        )
        
        # mapk/erk
        fig.add_trace(
            go.Scatter(x=df['time_hours'], y=df['nongenomic_mapk_active'],
                      name='mapk', line=dict(color='#e74c3c', width=2)),
            row=1, col=1
        )
        fig.add_trace(
            go.Scatter(x=df['time_hours'], y=df['nongenomic_erk_active'],
                      name='erk', line=dict(color='#c0392b', width=2, dash='dot')),
            row=1, col=1
        )
        
        # pi3k/akt
        fig.add_trace(
            go.Scatter(x=df['time_hours'], y=df['nongenomic_pi3k_active'],
                      name='pi3k', line=dict(color='#3498db', width=2),
                      showlegend=False),
            row=1, col=2
        )
        fig.add_trace(
            go.Scatter(x=df['time_hours'], y=df['nongenomic_akt_active'],
                      name='akt', line=dict(color='#2980b9', width=2, dash='dot'),
                      showlegend=False),
            row=1, col=2
        )
        
        # calcium
        fig.add_trace(
            go.Scatter(x=df['time_hours'], y=df['nongenomic_calcium_cytosol'],
                      name='cytosolic ca²⁺', line=dict(color='#2ecc71', width=2),
                      showlegend=False),
            row=2, col=1
        )
        
        # nitric oxide
        fig.add_trace(
            go.Scatter(x=df['time_hours'], y=df['nongenomic_no_concentration'],
                      name='no', line=dict(color='#9b59b6', width=2),
                      showlegend=False),
            row=2, col=2
        )
        
        fig.update_xaxes(title_text="time (hours)")
        fig.update_yaxes(title_text="activity", row=1, col=1)
        fig.update_yaxes(title_text="activity", row=1, col=2)
        fig.update_yaxes(title_text="[ca²⁺] (μm)", row=2, col=1)
        fig.update_yaxes(title_text="[no] (μm)", row=2, col=2)
        
        fig.update_layout(height=800, title_text="non-genomic signaling pathways")
        
        return fig
    
    @staticmethod
    def plot_cellular_outcomes(df: pd.DataFrame) -> go.Figure:
        """plot cellular phenotype outcomes"""
        fig = make_subplots(
            rows=2, cols=2,
            subplot_titles=(
                'proliferation rate',
                'cell cycle progression',
                'survival signaling',
                'cumulative divisions'
            )
        )
        
        # proliferation
        fig.add_trace(
            go.Scatter(x=df['time_hours'], y=df['cell_proliferation_rate'],
                      line=dict(color='#e74c3c', width=2),
                      fill='tozeroy', showlegend=False),
            row=1, col=1
        )
        
        # cell cycle
        fig.add_trace(
            go.Scatter(x=df['time_hours'], y=df['cell_cell_cycle_progress'],
                      line=dict(color='#3498db', width=2), showlegend=False),
            row=1, col=2
        )
        
        # survival
        fig.add_trace(
            go.Scatter(x=df['time_hours'], y=df['cell_survival_signal'],
                      name='survival', line=dict(color='#2ecc71', width=2)),
            row=2, col=1
        )
        fig.add_trace(
            go.Scatter(x=df['time_hours'], y=df['cell_apoptosis_index'],
                      name='apoptosis', line=dict(color='#e67e22', width=2, dash='dot')),
            row=2, col=1
        )
        
        # divisions
        fig.add_trace(
            go.Scatter(x=df['time_hours'], y=df['cell_division_count'],
                      line=dict(color='#9b59b6', width=2),
                      mode='lines+markers', showlegend=False),
            row=2, col=2
        )
        
        fig.update_xaxes(title_text="time (hours)")
        fig.update_yaxes(title_text="rate (au)", row=1, col=1)
        fig.update_yaxes(title_text="progress", row=1, col=2)
        fig.update_yaxes(title_text="signal", row=2, col=1)
        fig.update_yaxes(title_text="divisions", row=2, col=2)
        
        fig.update_layout(height=800, title_text="cellular outcomes")
        
        return fig
    
    @staticmethod
    def plot_dose_response(dose_response_df: pd.DataFrame, 
                          response_column: str,
                          title: str = "dose-response curve") -> go.Figure:
        """plot dose-response relationship"""
        fig = go.Figure()
        
        fig.add_trace(go.Scatter(
            x=dose_response_df['dose'],
            y=dose_response_df[response_column],
            mode='lines+markers',
            line=dict(color='#e74c3c', width=3),
            marker=dict(size=10, color='#c0392b')
        ))
        
        fig.update_xaxes(type="log", title="dose (mg)")
        fig.update_yaxes(title=response_column)
        fig.update_layout(title=title, height=500)
        
        return fig
    
    @staticmethod
    def plot_heatmap(data_dict: Dict[str, pd.DataFrame],
                    column: str,
                    title: str = "heatmap") -> go.Figure:
        """create heatmap comparing conditions"""
        # extract time series for each condition
        matrix = []
        labels = []
        
        for label, df in data_dict.items():
            if column in df.columns:
                matrix.append(df[column].values)
                labels.append(label)
        
        if not matrix:
            return go.Figure()
        
        # get time points from first dataframe
        time_hours = list(data_dict.values())[0]['time_hours'].values
        
        fig = go.Figure(data=go.Heatmap(
            z=matrix,
            x=time_hours,
            y=labels,
            colorscale='RdYlBu_r',
            colorbar=dict(title=column)
        ))
        
        fig.update_xaxes(title="time (hours)")
        fig.update_layout(title=title, height=400)
        
        return fig


class NetworkVisualizer:
    """
    visualize er signaling network topology
    """
    
    @staticmethod
    def create_pathway_network() -> go.Figure:
        """create interactive network diagram of er pathways"""
        
        # define nodes
        nodes = {
            # ligands
            'estradiol': {'x': 0, 'y': 5, 'color': '#e74c3c', 'size': 20},
            
            # receptors
            'er_alpha': {'x': 2, 'y': 5, 'color': '#e74c3c', 'size': 30},
            'er_beta': {'x': 2, 'y': 3, 'color': '#3498db', 'size': 25},
            
            # nuclear
            'nucleus': {'x': 4, 'y': 4, 'color': '#95a5a6', 'size': 15},
            'ere': {'x': 6, 'y': 4, 'color': '#7f8c8d', 'size': 15},
            
            # target genes
            'cyclin_d1': {'x': 8, 'y': 5, 'color': '#2ecc71', 'size': 20},
            'bcl2': {'x': 8, 'y': 3, 'color': '#2ecc71', 'size': 20},
            
            # non-genomic
            'membrane_er': {'x': 2, 'y': 7, 'color': '#f39c12', 'size': 20},
            'mapk': {'x': 4, 'y': 8, 'color': '#9b59b6', 'size': 20},
            'pi3k': {'x': 4, 'y': 6, 'color': '#1abc9c', 'size': 20},
            'akt': {'x': 6, 'y': 6, 'color': '#16a085', 'size': 20},
            
            # outcomes
            'proliferation': {'x': 10, 'y': 6, 'color': '#e67e22', 'size': 25},
            'survival': {'x': 10, 'y': 4, 'color': '#27ae60', 'size': 25},
        }
        
        # define edges
        edges = [
            ('estradiol', 'er_alpha'),
            ('estradiol', 'er_beta'),
            ('er_alpha', 'nucleus'),
            ('er_beta', 'nucleus'),
            ('nucleus', 'ere'),
            ('ere', 'cyclin_d1'),
            ('ere', 'bcl2'),
            ('estradiol', 'membrane_er'),
            ('membrane_er', 'mapk'),
            ('membrane_er', 'pi3k'),
            ('pi3k', 'akt'),
            ('mapk', 'proliferation'),
            ('cyclin_d1', 'proliferation'),
            ('akt', 'survival'),
            ('bcl2', 'survival'),
        ]
        
        # create figure
        fig = go.Figure()
        
        # add edges
        for source, target in edges:
            fig.add_trace(go.Scatter(
                x=[nodes[source]['x'], nodes[target]['x']],
                y=[nodes[source]['y'], nodes[target]['y']],
                mode='lines',
                line=dict(color='#bdc3c7', width=2),
                hoverinfo='none',
                showlegend=False
            ))
        
        # add nodes
        for name, props in nodes.items():
            fig.add_trace(go.Scatter(
                x=[props['x']],
                y=[props['y']],
                mode='markers+text',
                marker=dict(size=props['size'], color=props['color'],
                          line=dict(color='white', width=2)),
                text=name.replace('_', ' '),
                textposition='top center',
                name=name,
                showlegend=False
            ))
        
        fig.update_layout(
            title="estrogen receptor signaling network",
            xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
            yaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
            height=600,
            plot_bgcolor='white'
        )
        
        return fig
