"""
Streamlit App for PPAR Signaling Simulation
Interactive visualization of PPAR-α, PPAR-γ, and PPAR-δ pathways
"""

import streamlit as st
import numpy as np
import pandas as pd
import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots
import seaborn as sns
import matplotlib.pyplot as plt

from ppar_model import PPARSignalingModel, PPARParameters, compare_isoform_specificity


# Page configuration
st.set_page_config(
    page_title="PPAR Signaling Simulator",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded"
)

# Custom CSS
st.markdown("""
    <style>
    .main-header {
        font-size: 3rem;
        font-weight: bold;
        text-align: center;
        color: #1f77b4;
        margin-bottom: 0.5rem;
    }
    .sub-header {
        font-size: 1.2rem;
        text-align: center;
        color: #666;
        margin-bottom: 2rem;
    }
    .metric-box {
        background-color: #f0f2f6;
        padding: 1rem;
        border-radius: 0.5rem;
        margin: 0.5rem 0;
    }
    .info-box {
        background-color: #e8f4f8;
        padding: 1rem;
        border-radius: 0.5rem;
        border-left: 4px solid #1f77b4;
        margin: 1rem 0;
    }
    </style>
""", unsafe_allow_html=True)


def create_header():
    """Create app header"""
    st.markdown('<div class="main-header">PPAR Signaling Simulator</div>', unsafe_allow_html=True)
    st.markdown(
        '<div class="sub-header">peroxisome proliferator-activated receptors: α, γ, δ isoforms</div>',
        unsafe_allow_html=True
    )
    
    st.write("""
**about ppar signaling:** ppars are nuclear receptor transcription factors that regulate 
lipid metabolism, glucose homeostasis, and inflammatory responses. the three isoforms have 
distinct but overlapping functions:

- **ppar-α:** fatty acid oxidation, ketogenesis (liver, muscle, heart)
- **ppar-γ:** adipogenesis, insulin sensitivity (adipose tissue)
- **ppar-δ:** fatty acid oxidation, anti-inflammatory (ubiquitous)
    """)


def plot_ppar_dynamics(df: pd.DataFrame, title: str = "PPAR Signaling Dynamics"):
    """Plot main PPAR signaling dynamics"""
    fig = make_subplots(
        rows=3, cols=2,
        subplot_titles=(
            'PPAR Isoform Activation',
            'Gene Expression (mRNA)',
            'Protein Levels',
            'Metabolic Outputs',
            'Nuclear Receptor Dynamics',
            'Health Metrics'
        ),
        vertical_spacing=0.12,
        horizontal_spacing=0.1
    )
    
    # 1. PPAR-RXR-DNA complexes (active transcription)
    for ppar, color in [('PPAR-α-RXR-DNA', '#ff7f0e'), 
                        ('PPAR-γ-RXR-DNA', '#2ca02c'), 
                        ('PPAR-δ-RXR-DNA', '#d62728')]:
        fig.add_trace(
            go.Scatter(x=df['Time'], y=df[ppar], name=ppar.split('-')[0],
                      line=dict(color=color, width=2.5),
                      legendgroup='ppar'),
            row=1, col=1
        )
    
    # 2. mRNA levels
    mrna_cols = ['mRNA_FAO_α', 'mRNA_IS_γ', 'mRNA_AI_δ']
    mrna_names = ['fatty acid oxidation', 'insulin sensitivity', 'anti-inflammatory']
    colors_mrna = ['#ff7f0e', '#2ca02c', '#d62728']
    
    for mrna, name, color in zip(mrna_cols, mrna_names, colors_mrna):
        fig.add_trace(
            go.Scatter(x=df['Time'], y=df[mrna], name=name,
                      line=dict(color=color, width=2, dash='dot'),
                      legendgroup='mrna'),
            row=1, col=2
        )
    
    # 3. Protein levels
    proteins = ['Protein_FAO', 'Protein_IS', 'Protein_AI']
    protein_names = ['fatty acid oxidation', 'glucose uptake', 'anti-inflammatory']
    
    for prot, name, color in zip(proteins, protein_names, colors_mrna):
        fig.add_trace(
            go.Scatter(x=df['Time'], y=df[prot], name=name,
                      line=dict(color=color, width=2.5),
                      legendgroup='protein'),
            row=2, col=1
        )
    
    # 4. Metabolic outputs
    metabolic = [
        ('Lipid_Accumulation', 'lipid accumulation', '#e377c2'),
        ('Insulin_Sensitivity', 'insulin sensitivity', '#7f7f7f'),
        ('Inflammatory_State', 'inflammation', '#bcbd22')
    ]
    
    for var, name, color in metabolic:
        fig.add_trace(
            go.Scatter(x=df['Time'], y=df[var], name=name,
                      line=dict(color=color, width=2.5),
                      legendgroup='metabolic'),
            row=2, col=2
        )
    
    # 5. Nuclear receptor dynamics
    rxr_data = [
        ('RXR_free', 'free RXR', '#17becf'),
        ('PPAR-α-RXR', 'PPAR-α-RXR', '#ff7f0e'),
        ('PPAR-γ-RXR', 'PPAR-γ-RXR', '#2ca02c'),
    ]
    
    for var, name, color in rxr_data:
        fig.add_trace(
            go.Scatter(x=df['Time'], y=df[var], name=name,
                      line=dict(color=color, width=2),
                      legendgroup='rxr'),
            row=3, col=1
        )
    
    # 6. Overall health metric
    fig.add_trace(
        go.Scatter(x=df['Time'], y=df['Metabolic_Health'], 
                  name='metabolic health score',
                  line=dict(color='#1f77b4', width=3),
                  fill='tozeroy', fillcolor='rgba(31, 119, 180, 0.2)',
                  legendgroup='health'),
        row=3, col=2
    )
    
    # Update axes
    fig.update_xaxes(title_text="time (h)", row=3, col=1)
    fig.update_xaxes(title_text="time (h)", row=3, col=2)
    fig.update_yaxes(title_text="concentration", row=1, col=1)
    fig.update_yaxes(title_text="mrna level", row=1, col=2)
    fig.update_yaxes(title_text="protein level", row=2, col=1)
    fig.update_yaxes(title_text="state", row=2, col=2)
    fig.update_yaxes(title_text="concentration", row=3, col=1)
    fig.update_yaxes(title_text="health score", row=3, col=2)
    
    fig.update_layout(
        height=1200,
        title_text=title,
        title_font_size=20,
        showlegend=True,
        template='plotly_white',
        legend=dict(orientation="v", yanchor="top", y=1, xanchor="left", x=1.02)
    )
    
    return fig


def plot_dose_response(results: dict, metric: str):
    """Plot dose-response curves"""
    doses = []
    values = []
    
    for dose_key, df in results.items():
        dose = float(dose_key.split('_')[1])
        steady_state = df.iloc[-len(df)//10:]
        value = steady_state[metric].mean()
        doses.append(dose)
        values.append(value)
    
    fig = go.Figure()
    fig.add_trace(go.Scatter(
        x=doses, y=values,
        mode='lines+markers',
        line=dict(color='#1f77b4', width=3),
        marker=dict(size=10, color='#ff7f0e')
    ))
    
    fig.update_layout(
        title=f'dose-response curve: {metric}',
        xaxis_title='ligand concentration',
        yaxis_title=metric.replace('_', ' ').lower(),
        template='plotly_white',
        height=400
    )
    
    return fig


def plot_comparison_heatmap(comparison_df: pd.DataFrame):
    """Plot heatmap comparing isoform activities"""
    # Prepare data
    metrics = ['PPAR-α Activity', 'PPAR-γ Activity', 'PPAR-δ Activity',
               'Fatty Acid Oxidation', 'Insulin Sensitivity', 'Anti-inflammatory']
    
    data = comparison_df[metrics].values
    isoforms = comparison_df['Isoform'].values
    
    fig = go.Figure(data=go.Heatmap(
        z=data.T,
        x=isoforms,
        y=metrics,
        colorscale='YlOrRd',
        text=np.round(data.T, 2),
        texttemplate='%{text}',
        textfont={"size": 12},
        colorbar=dict(title="activity level")
    ))
    
    fig.update_layout(
        title='isoform-specific activity profiles',
        xaxis_title='ppar isoform',
        yaxis_title='',
        height=500,
        template='plotly_white'
    )
    
    return fig


def plot_pathway_network():
    """Create network visualization of PPAR signaling"""
    fig = go.Figure()
    
    # Node positions (manually arranged for clarity)
    nodes = {
        'Ligands': (0, 3),
        'PPAR-α': (2, 4), 'PPAR-γ': (2, 3), 'PPAR-δ': (2, 2),
        'RXR': (4, 3),
        'PPAR-RXR': (6, 3),
        'DNA': (8, 3),
        'FAO Genes': (10, 4), 'IS Genes': (10, 3), 'AI Genes': (10, 2),
        'Lipid Metabolism': (12, 4),
        'Insulin Sensitivity': (12, 3),
        'Inflammation': (12, 2)
    }
    
    # Add edges
    edges = [
        ('Ligands', 'PPAR-α'), ('Ligands', 'PPAR-γ'), ('Ligands', 'PPAR-δ'),
        ('PPAR-α', 'PPAR-RXR'), ('PPAR-γ', 'PPAR-RXR'), ('PPAR-δ', 'PPAR-RXR'),
        ('RXR', 'PPAR-RXR'),
        ('PPAR-RXR', 'DNA'),
        ('DNA', 'FAO Genes'), ('DNA', 'IS Genes'), ('DNA', 'AI Genes'),
        ('FAO Genes', 'Lipid Metabolism'),
        ('IS Genes', 'Insulin Sensitivity'),
        ('AI Genes', 'Inflammation')
    ]
    
    # Draw edges
    for edge in edges:
        x0, y0 = nodes[edge[0]]
        x1, y1 = nodes[edge[1]]
        fig.add_trace(go.Scatter(
            x=[x0, x1], y=[y0, y1],
            mode='lines',
            line=dict(color='gray', width=2),
            hoverinfo='none',
            showlegend=False
        ))
    
    # Draw nodes
    node_colors = ['#ff7f0e', '#ff7f0e', '#2ca02c', '#d62728', '#17becf', 
                   '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22',
                   '#ff9896', '#c5b0d5', '#c49c94']
    
    for (name, pos), color in zip(nodes.items(), node_colors):
        fig.add_trace(go.Scatter(
            x=[pos[0]], y=[pos[1]],
            mode='markers+text',
            marker=dict(size=30, color=color),
            text=name,
            textposition='middle center',
            textfont=dict(size=9, color='white', family='Arial Black'),
            hoverinfo='text',
            hovertext=name,
            showlegend=False
        ))
    
    fig.update_layout(
        title='ppar signaling pathway network',
        showlegend=False,
        xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
        yaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
        height=600,
        template='plotly_white',
        hovermode='closest'
    )
    
    return fig


def main():
    """Main app function"""
    create_header()
    
    # Sidebar controls
    st.sidebar.header("simulation parameters")
    
    sim_mode = st.sidebar.selectbox(
        "simulation mode",
        ["basic dynamics", "drug response", "isoform comparison", "pathway network"]
    )
    
    if sim_mode == "basic dynamics":
        st.sidebar.subheader("ligand concentrations")
        
        ligand_alpha = st.sidebar.slider(
            "ppar-α ligand (fibrates)",
            0.0, 5.0, 0.0, 0.1,
            help="e.g., fenofibrate, gemfibrozil"
        )
        
        ligand_gamma = st.sidebar.slider(
            "ppar-γ ligand (tzds)",
            0.0, 5.0, 0.0, 0.1,
            help="e.g., pioglitazone, rosiglitazone"
        )
        
        ligand_delta = st.sidebar.slider(
            "ppar-δ ligand",
            0.0, 5.0, 0.0, 0.1,
            help="e.g., gw501516, gw0742"
        )
        
        st.sidebar.subheader("simulation settings")
        
        t_max = st.sidebar.slider("simulation time (h)", 10, 200, 100, 10)
        
        # Advanced parameters
        with st.sidebar.expander("advanced parameters"):
            ppar_alpha_base = st.slider("ppar-α expression", 0.1, 2.0, 1.0, 0.1)
            ppar_gamma_base = st.slider("ppar-γ expression", 0.1, 2.0, 1.0, 0.1)
            ppar_delta_base = st.slider("ppar-δ expression", 0.1, 2.0, 1.0, 0.1)
            free_fatty_acids = st.slider("free fatty acids", 0.0, 3.0, 1.0, 0.1)
            inflammatory_cytokines = st.slider("inflammatory cytokines", 0.0, 2.0, 0.5, 0.1)
        
        if st.sidebar.button("run simulation", type="primary"):
            with st.spinner("simulating ppar signaling..."):
                # Create model with custom parameters
                params = PPARParameters(
                    ppar_alpha_base=ppar_alpha_base,
                    ppar_gamma_base=ppar_gamma_base,
                    ppar_delta_base=ppar_delta_base,
                    free_fatty_acids=free_fatty_acids,
                    inflammatory_cytokines=inflammatory_cytokines
                )
                
                model = PPARSignalingModel(params)
                df = model.simulate(
                    (0, t_max), 
                    1000,
                    ligand_alpha,
                    ligand_gamma,
                    ligand_delta
                )
                
                # Store in session state
                st.session_state.simulation_df = df
                st.session_state.model = model
        
        # Display results
        if 'simulation_df' in st.session_state:
            df = st.session_state.simulation_df
            model = st.session_state.model
            
            # Main dynamics plot
            st.plotly_chart(plot_ppar_dynamics(df), use_container_width=True)
            
            # Metrics
            st.subheader("steady-state metrics")
            metrics = model.get_steady_state_metrics(df)
            
            col1, col2, col3, col4 = st.columns(4)
            
            with col1:
                st.metric("ppar-α activity", f"{metrics['PPAR-α Activity']:.3f}")
                st.metric("fatty acid oxidation", f"{metrics['Fatty Acid Oxidation']:.3f}")
            
            with col2:
                st.metric("ppar-γ activity", f"{metrics['PPAR-γ Activity']:.3f}")
                st.metric("insulin sensitivity", f"{metrics['Insulin Sensitivity']:.3f}")
            
            with col3:
                st.metric("ppar-δ activity", f"{metrics['PPAR-δ Activity']:.3f}")
                st.metric("anti-inflammatory", f"{metrics['Anti-inflammatory']:.3f}")
            
            with col4:
                st.metric("lipid accumulation", f"{metrics['Lipid Accumulation']:.3f}")
                st.metric("metabolic health", f"{metrics['Metabolic Health']:.3f}")
            
            # Download data
            st.subheader("export data")
            csv = df.to_csv(index=False)
            st.download_button(
                label="download simulation data (csv)",
                data=csv,
                file_name="ppar_simulation.csv",
                mime="text/csv"
            )
    
    elif sim_mode == "drug response":
        st.subheader("dose-response analysis")
        
        drug_type = st.sidebar.selectbox(
            "drug type",
            ["alpha", "gamma", "delta", "pan"],
            format_func=lambda x: {
                "alpha": "ppar-α agonist (fibrate)",
                "gamma": "ppar-γ agonist (tzd)",
                "delta": "ppar-δ agonist",
                "pan": "pan-ppar agonist"
            }[x]
        )
        
        dose_min = st.sidebar.number_input("min dose", 0.0, 10.0, 0.0, 0.5)
        dose_max = st.sidebar.number_input("max dose", 0.0, 10.0, 5.0, 0.5)
        n_doses = st.sidebar.slider("number of doses", 3, 15, 8)
        
        if st.sidebar.button("run dose-response", type="primary"):
            with st.spinner("running dose-response analysis..."):
                model = PPARSignalingModel()
                dose_range = np.linspace(dose_min, dose_max, n_doses)
                results = model.simulate_drug_response(drug_type, dose_range)
                
                st.session_state.dose_results = results
                st.session_state.drug_type = drug_type
        
        if 'dose_results' in st.session_state:
            results = st.session_state.dose_results
            
            # Select metric to plot
            metric = st.selectbox(
                "select metric",
                ["Insulin_Sensitivity", "Lipid_Accumulation", "Inflammatory_State",
                 "Metabolic_Health", "Protein_FAO", "Protein_IS"]
            )
            
            col1, col2 = st.columns(2)
            
            with col1:
                st.plotly_chart(plot_dose_response(results, metric), use_container_width=True)
            
            with col2:
                # Time course for selected doses
                selected_dose = st.select_slider(
                    "select dose for time course",
                    options=list(results.keys()),
                    format_func=lambda x: f"dose: {x.split('_')[1]}"
                )
                
                df_selected = results[selected_dose]
                
                fig = go.Figure()
                fig.add_trace(go.Scatter(
                    x=df_selected['Time'],
                    y=df_selected[metric],
                    mode='lines',
                    line=dict(color='#1f77b4', width=2.5)
                ))
                fig.update_layout(
                    title=f'{metric} time course',
                    xaxis_title='time (h)',
                    yaxis_title=metric.replace('_', ' ').lower(),
                    template='plotly_white',
                    height=400
                )
                st.plotly_chart(fig, use_container_width=True)
    
    elif sim_mode == "isoform comparison":
        st.subheader("isoform-specific activity analysis")
        
        if st.sidebar.button("compare isoforms", type="primary"):
            with st.spinner("analyzing isoform specificity..."):
                comparison_df = compare_isoform_specificity()
                st.session_state.comparison_df = comparison_df
        
        if 'comparison_df' in st.session_state:
            comparison_df = st.session_state.comparison_df
            
            st.plotly_chart(plot_comparison_heatmap(comparison_df), use_container_width=True)
            
            st.subheader("detailed comparison")
            st.dataframe(comparison_df.set_index('Isoform'), use_container_width=True)
            
            # Bar chart comparison
            fig = go.Figure()
            metrics = ['Fatty Acid Oxidation', 'Insulin Sensitivity', 'Anti-inflammatory']
            colors = ['#ff7f0e', '#2ca02c', '#d62728']
            
            for metric, color in zip(metrics, colors):
                fig.add_trace(go.Bar(
                    name=metric,
                    x=comparison_df['Isoform'],
                    y=comparison_df[metric],
                    marker_color=color
                ))
            
            fig.update_layout(
                title='primary functions by isoform',
                xaxis_title='ppar isoform',
                yaxis_title='activity level',
                barmode='group',
                template='plotly_white',
                height=500
            )
            
            st.plotly_chart(fig, use_container_width=True)
    
    else:  # pathway network
        st.subheader("ppar signaling pathway")
        st.plotly_chart(plot_pathway_network(), use_container_width=True)
        
        st.markdown("""
        ### pathway components
        
        **ligand binding:**
        - fatty acids, eicosanoids (endogenous)
        - fibrates (ppar-α agonists)
        - thiazolidinediones (ppar-γ agonists)
        - synthetic agonists (ppar-δ)
        
        **heterodimer formation:**
        - ppar isoforms form obligate heterodimers with rxr (retinoid x receptor)
        - rxr availability can limit ppar signaling
        
        **target genes:**
        - **ppar-α:** cpt1a, acox1, hmgcs2 (lipid oxidation, ketogenesis)
        - **ppar-γ:** glut4, fabp4, cebpa (glucose uptake, adipogenesis)
        - **ppar-δ:** pdk4, angptl4, cd36 (fatty acid metabolism, inflammation)
        
        **metabolic effects:**
        - improved lipid profiles
        - enhanced insulin sensitivity
        - reduced inflammation
        - increased energy expenditure
        """)
    
    # Footer
    st.sidebar.markdown("---")
    st.sidebar.markdown("""
    **model features:**
    - 28 state variables
    - 3 ppar isoforms
    - ligand binding kinetics
    - rxr heterodimerization
    - dna binding & transcription
    - metabolic outputs
    - feedback regulation
    """)


if __name__ == "__main__":
    main()
