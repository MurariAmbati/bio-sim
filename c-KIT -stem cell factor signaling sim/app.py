import streamlit as st
import numpy as np
import pandas as pd
import plotly.graph_objects as go
import plotly.express as px
from scipy.integrate import odeint
from dataclasses import dataclass
from typing import Tuple, Dict, List
import matplotlib.pyplot as plt

# Page configuration
st.set_page_config(
    page_title="c-KIT/SCF Signaling Simulation",
    page_icon="⚗",
    layout="wide",
    initial_sidebar_state="expanded"
)

# Custom CSS
st.markdown("""
<style>
    .main-header {
        font-size: 2.5rem;
        color: #2c3e50;
        font-weight: bold;
        margin-bottom: 0.5rem;
    }
    .sub-header {
        font-size: 1.2rem;
        color: #7f8c8d;
        margin-bottom: 2rem;
    }
    .metric-card {
        background-color: #f8f9fa;
        padding: 1rem;
        border-radius: 0.5rem;
        border-left: 4px solid #3498db;
    }
</style>
""", unsafe_allow_html=True)

@dataclass
class SignalingParameters:
    """Parameters for c-KIT/SCF signaling pathway"""
    # Receptor binding and activation
    k_bind: float = 0.1  # SCF-cKIT binding rate
    k_unbind: float = 0.05  # Unbinding rate
    k_dim: float = 0.2  # Dimerization rate
    k_undim: float = 0.01  # Dimer dissociation rate
    k_phos: float = 0.15  # Phosphorylation rate
    k_dephos: float = 0.08  # Dephosphorylation rate
    
    # PI3K/AKT pathway
    k_pi3k_act: float = 0.12  # PI3K activation
    k_pi3k_deact: float = 0.06  # PI3K deactivation
    k_akt_act: float = 0.18  # AKT activation
    k_akt_deact: float = 0.09  # AKT deactivation
    
    # MAPK/ERK pathway
    k_ras_act: float = 0.14  # RAS activation
    k_ras_deact: float = 0.07  # RAS deactivation
    k_mek_act: float = 0.16  # MEK activation
    k_mek_deact: float = 0.08  # MEK deactivation
    k_erk_act: float = 0.20  # ERK activation
    k_erk_deact: float = 0.10  # ERK deactivation
    
    # JAK/STAT pathway
    k_jak_act: float = 0.13  # JAK activation
    k_jak_deact: float = 0.065  # JAK deactivation
    k_stat_act: float = 0.17  # STAT activation
    k_stat_deact: float = 0.085  # STAT deactivation
    
    # Downstream effects
    k_prolif: float = 0.05  # Cell proliferation rate
    k_diff: float = 0.03  # Cell differentiation rate
    k_survival: float = 0.04  # Cell survival rate
    k_apop: float = 0.02  # Apoptosis rate
    
    # Degradation and internalization
    k_deg: float = 0.04  # General degradation rate
    k_intern: float = 0.06  # Receptor internalization


class CKITSignalingModel:
    """
    Comprehensive model of c-KIT/SCF signaling pathway
    
    Components:
    - SCF (Stem Cell Factor)
    - c-KIT receptor (monomers and dimers)
    - PI3K/AKT pathway
    - MAPK/ERK pathway
    - JAK/STAT pathway
    - Cell fate outcomes
    """
    
    def __init__(self, params: SignalingParameters):
        self.params = params
        
    def model_equations(self, state: np.ndarray, t: float, scf_input: float) -> List[float]:
        """
        Differential equations for c-KIT/SCF signaling
        
        State variables:
        0: SCF (free)
        1: cKIT (free receptor)
        2: SCF-cKIT complex
        3: cKIT dimer (activated)
        4: cKIT-p (phosphorylated)
        5: PI3K active
        6: AKT active
        7: RAS active
        8: MEK active
        9: ERK active
        10: JAK active
        11: STAT active
        12: Proliferation signal
        13: Differentiation signal
        14: Survival signal
        15: Apoptosis signal
        """
        
        SCF, cKIT, SCF_cKIT, cKIT_dim, cKIT_p, PI3K_a, AKT_a, RAS_a, MEK_a, ERK_a, JAK_a, STAT_a, Prolif, Diff, Surv, Apop = state
        
        p = self.params
        
        # Receptor binding and activation
        d_SCF = -p.k_bind * SCF * cKIT + p.k_unbind * SCF_cKIT + scf_input - p.k_deg * SCF
        d_cKIT = -p.k_bind * SCF * cKIT + p.k_unbind * SCF_cKIT - p.k_deg * cKIT + 0.1
        d_SCF_cKIT = p.k_bind * SCF * cKIT - p.k_unbind * SCF_cKIT - p.k_dim * SCF_cKIT
        d_cKIT_dim = p.k_dim * SCF_cKIT - p.k_undim * cKIT_dim - p.k_phos * cKIT_dim
        d_cKIT_p = p.k_phos * cKIT_dim - p.k_dephos * cKIT_p - p.k_intern * cKIT_p
        
        # PI3K/AKT pathway
        d_PI3K_a = p.k_pi3k_act * cKIT_p * (1 - PI3K_a) - p.k_pi3k_deact * PI3K_a
        d_AKT_a = p.k_akt_act * PI3K_a * (1 - AKT_a) - p.k_akt_deact * AKT_a
        
        # MAPK/ERK pathway
        d_RAS_a = p.k_ras_act * cKIT_p * (1 - RAS_a) - p.k_ras_deact * RAS_a
        d_MEK_a = p.k_mek_act * RAS_a * (1 - MEK_a) - p.k_mek_deact * MEK_a
        d_ERK_a = p.k_erk_act * MEK_a * (1 - ERK_a) - p.k_erk_deact * ERK_a
        
        # JAK/STAT pathway
        d_JAK_a = p.k_jak_act * cKIT_p * (1 - JAK_a) - p.k_jak_deact * JAK_a
        d_STAT_a = p.k_stat_act * JAK_a * (1 - STAT_a) - p.k_stat_deact * STAT_a
        
        # Cell fate outcomes
        d_Prolif = p.k_prolif * (ERK_a + AKT_a) - p.k_deg * Prolif
        d_Diff = p.k_diff * (STAT_a + 0.5 * ERK_a) - p.k_deg * Diff
        d_Surv = p.k_survival * (AKT_a + STAT_a) - p.k_apop * Surv
        d_Apop = p.k_apop * (1 - AKT_a) * (1 - Surv) - p.k_deg * Apop
        
        return [d_SCF, d_cKIT, d_SCF_cKIT, d_cKIT_dim, d_cKIT_p, d_PI3K_a, d_AKT_a, 
                d_RAS_a, d_MEK_a, d_ERK_a, d_JAK_a, d_STAT_a, d_Prolif, d_Diff, d_Surv, d_Apop]
    
    def simulate(self, initial_conditions: np.ndarray, time_points: np.ndarray, 
                 scf_profile: np.ndarray) -> pd.DataFrame:
        """Run simulation with given initial conditions and SCF profile"""
        
        solutions = []
        for i, t in enumerate(time_points[:-1]):
            sol = odeint(self.model_equations, initial_conditions, 
                        [time_points[i], time_points[i+1]], 
                        args=(scf_profile[i],))
            solutions.append(sol[-1])
            initial_conditions = sol[-1]
        
        solutions = np.array(solutions)
        
        df = pd.DataFrame(solutions, columns=[
            'SCF', 'cKIT', 'SCF-cKIT', 'cKIT_dimer', 'cKIT_p', 
            'PI3K_active', 'AKT_active', 'RAS_active', 'MEK_active', 'ERK_active',
            'JAK_active', 'STAT_active', 'Proliferation', 'Differentiation', 
            'Survival', 'Apoptosis'
        ])
        df['Time'] = time_points[:-1]
        
        return df


def create_scf_profile(profile_type: str, time_points: np.ndarray, 
                       amplitude: float = 1.0, duration: float = 50.0) -> np.ndarray:
    """Generate different SCF stimulation profiles"""
    
    if profile_type == "Constant":
        return np.ones_like(time_points) * amplitude
    elif profile_type == "Pulse":
        profile = np.zeros_like(time_points)
        mask = (time_points > 10) & (time_points < 10 + duration)
        profile[mask] = amplitude
        return profile
    elif profile_type == "Oscillatory":
        return amplitude * (0.5 + 0.5 * np.sin(2 * np.pi * time_points / duration))
    elif profile_type == "Ramp":
        return amplitude * np.minimum(time_points / duration, 1.0)
    elif profile_type == "Decay":
        return amplitude * np.exp(-time_points / duration)
    else:
        return np.zeros_like(time_points)


def plot_receptor_dynamics(df: pd.DataFrame) -> go.Figure:
    """Visualize receptor binding and activation dynamics"""
    
    fig = go.Figure()
    
    fig.add_trace(go.Scatter(x=df['Time'], y=df['cKIT'], 
                            name='Free c-KIT', line=dict(color='#3498db', width=2)))
    fig.add_trace(go.Scatter(x=df['Time'], y=df['SCF-cKIT'], 
                            name='SCF-cKIT Complex', line=dict(color='#e74c3c', width=2)))
    fig.add_trace(go.Scatter(x=df['Time'], y=df['cKIT_dimer'], 
                            name='c-KIT Dimer', line=dict(color='#f39c12', width=2)))
    fig.add_trace(go.Scatter(x=df['Time'], y=df['cKIT_p'], 
                            name='Phospho-cKIT', line=dict(color='#2ecc71', width=2)))
    
    fig.update_layout(
        title="Receptor Dynamics",
        xaxis_title="Time (min)",
        yaxis_title="Concentration (AU)",
        template="plotly_white",
        hovermode='x unified',
        height=400
    )
    
    return fig


def plot_signaling_pathways(df: pd.DataFrame) -> go.Figure:
    """Visualize activation of major signaling pathways"""
    
    fig = go.Figure()
    
    # PI3K/AKT pathway
    fig.add_trace(go.Scatter(x=df['Time'], y=df['PI3K_active'], 
                            name='PI3K', line=dict(color='#9b59b6', width=2, dash='solid')))
    fig.add_trace(go.Scatter(x=df['Time'], y=df['AKT_active'], 
                            name='AKT', line=dict(color='#8e44ad', width=2, dash='dot')))
    
    # MAPK/ERK pathway
    fig.add_trace(go.Scatter(x=df['Time'], y=df['RAS_active'], 
                            name='RAS', line=dict(color='#e67e22', width=2, dash='solid')))
    fig.add_trace(go.Scatter(x=df['Time'], y=df['MEK_active'], 
                            name='MEK', line=dict(color='#d35400', width=2, dash='dot')))
    fig.add_trace(go.Scatter(x=df['Time'], y=df['ERK_active'], 
                            name='ERK', line=dict(color='#c0392b', width=2, dash='dash')))
    
    # JAK/STAT pathway
    fig.add_trace(go.Scatter(x=df['Time'], y=df['JAK_active'], 
                            name='JAK', line=dict(color='#16a085', width=2, dash='solid')))
    fig.add_trace(go.Scatter(x=df['Time'], y=df['STAT_active'], 
                            name='STAT', line=dict(color='#1abc9c', width=2, dash='dot')))
    
    fig.update_layout(
        title="Signaling Pathway Activation",
        xaxis_title="Time (min)",
        yaxis_title="Activation Level",
        template="plotly_white",
        hovermode='x unified',
        height=400
    )
    
    return fig


def plot_cell_fate(df: pd.DataFrame) -> go.Figure:
    """Visualize cell fate outcomes"""
    
    fig = go.Figure()
    
    fig.add_trace(go.Scatter(x=df['Time'], y=df['Proliferation'], 
                            name='Proliferation', fill='tonexty',
                            line=dict(color='#27ae60', width=2)))
    fig.add_trace(go.Scatter(x=df['Time'], y=df['Differentiation'], 
                            name='Differentiation', fill='tonexty',
                            line=dict(color='#3498db', width=2)))
    fig.add_trace(go.Scatter(x=df['Time'], y=df['Survival'], 
                            name='Survival', fill='tonexty',
                            line=dict(color='#f39c12', width=2)))
    fig.add_trace(go.Scatter(x=df['Time'], y=df['Apoptosis'], 
                            name='Apoptosis', fill='tonexty',
                            line=dict(color='#e74c3c', width=2)))
    
    fig.update_layout(
        title="Cell Fate Outcomes",
        xaxis_title="Time (min)",
        yaxis_title="Signal Strength",
        template="plotly_white",
        hovermode='x unified',
        height=400
    )
    
    return fig


def plot_pathway_heatmap(df: pd.DataFrame) -> go.Figure:
    """Create heatmap of pathway activation over time"""
    
    pathways = ['PI3K_active', 'AKT_active', 'RAS_active', 'MEK_active', 
                'ERK_active', 'JAK_active', 'STAT_active']
    pathway_names = ['PI3K', 'AKT', 'RAS', 'MEK', 'ERK', 'JAK', 'STAT']
    
    heatmap_data = df[pathways].T.values
    
    fig = go.Figure(data=go.Heatmap(
        z=heatmap_data,
        x=df['Time'],
        y=pathway_names,
        colorscale='RdYlGn',
        colorbar=dict(title="Activation")
    ))
    
    fig.update_layout(
        title="Pathway Activation Heatmap",
        xaxis_title="Time (min)",
        yaxis_title="Pathway Component",
        height=350
    )
    
    return fig


def plot_phase_portrait(df: pd.DataFrame, x_var: str, y_var: str) -> go.Figure:
    """Create phase portrait of two variables"""
    
    fig = go.Figure()
    
    # Create color gradient based on time
    colors = df['Time']
    
    fig.add_trace(go.Scatter(
        x=df[x_var], 
        y=df[y_var],
        mode='markers+lines',
        marker=dict(
            size=5,
            color=colors,
            colorscale='Viridis',
            showscale=True,
            colorbar=dict(title="Time (min)")
        ),
        line=dict(color='rgba(100,100,100,0.3)', width=1),
        name='Trajectory'
    ))
    
    # Mark start and end points
    fig.add_trace(go.Scatter(
        x=[df[x_var].iloc[0]], 
        y=[df[y_var].iloc[0]],
        mode='markers',
        marker=dict(size=12, color='green', symbol='star'),
        name='Start'
    ))
    
    fig.add_trace(go.Scatter(
        x=[df[x_var].iloc[-1]], 
        y=[df[y_var].iloc[-1]],
        mode='markers',
        marker=dict(size=12, color='red', symbol='x'),
        name='End'
    ))
    
    fig.update_layout(
        title=f"Phase Portrait: {y_var} vs {x_var}",
        xaxis_title=x_var.replace('_', ' ').title(),
        yaxis_title=y_var.replace('_', ' ').title(),
        template="plotly_white",
        height=400
    )
    
    return fig


def main():
    """Main Streamlit application"""
    
    # Header
    st.markdown('<div class="main-header">c-KIT / SCF Signaling Simulation</div>', 
                unsafe_allow_html=True)
    st.markdown('<div class="sub-header">Hematopoietic and Mast Cell Regulation</div>', 
                unsafe_allow_html=True)
    
    # Sidebar controls
    st.sidebar.header("Simulation Parameters")
    
    # SCF stimulation profile
    st.sidebar.subheader("SCF Stimulation")
    scf_profile_type = st.sidebar.selectbox(
        "Profile Type",
        ["Constant", "Pulse", "Oscillatory", "Ramp", "Decay"]
    )
    
    scf_amplitude = st.sidebar.slider("SCF Amplitude", 0.0, 5.0, 1.0, 0.1)
    scf_duration = st.sidebar.slider("Duration/Period (min)", 10.0, 200.0, 50.0, 10.0)
    
    # Simulation time
    st.sidebar.subheader("Time Course")
    sim_duration = st.sidebar.slider("Simulation Duration (min)", 50, 500, 200, 50)
    time_steps = st.sidebar.slider("Time Steps", 100, 2000, 500, 100)
    
    # Advanced parameters
    with st.sidebar.expander("Advanced Parameters"):
        st.subheader("Receptor Dynamics")
        k_bind = st.slider("Binding Rate", 0.01, 0.5, 0.1, 0.01)
        k_phos = st.slider("Phosphorylation Rate", 0.01, 0.5, 0.15, 0.01)
        
        st.subheader("Pathway Activation")
        k_akt = st.slider("AKT Activation", 0.01, 0.5, 0.18, 0.01)
        k_erk = st.slider("ERK Activation", 0.01, 0.5, 0.20, 0.01)
        k_stat = st.slider("STAT Activation", 0.01, 0.5, 0.17, 0.01)
    
    # Create parameters object with user inputs
    params = SignalingParameters(
        k_bind=k_bind,
        k_phos=k_phos,
        k_akt_act=k_akt,
        k_erk_act=k_erk,
        k_stat_act=k_stat
    )
    
    # Run simulation button
    if st.sidebar.button("Run Simulation", type="primary"):
        
        with st.spinner("Running simulation..."):
            # Initialize model
            model = CKITSignalingModel(params)
            
            # Time points
            time_points = np.linspace(0, sim_duration, time_steps)
            
            # SCF profile
            scf_profile = create_scf_profile(scf_profile_type, time_points, 
                                            scf_amplitude, scf_duration)
            
            # Initial conditions
            initial_conditions = np.array([
                0.5,  # SCF
                1.0,  # cKIT
                0.0,  # SCF-cKIT
                0.0,  # cKIT_dimer
                0.0,  # cKIT_p
                0.0,  # PI3K_active
                0.0,  # AKT_active
                0.0,  # RAS_active
                0.0,  # MEK_active
                0.0,  # ERK_active
                0.0,  # JAK_active
                0.0,  # STAT_active
                0.0,  # Proliferation
                0.0,  # Differentiation
                0.0,  # Survival
                0.0   # Apoptosis
            ])
            
            # Run simulation
            results_df = model.simulate(initial_conditions, time_points, scf_profile)
            
            # Store in session state
            st.session_state['results'] = results_df
            st.session_state['scf_profile'] = scf_profile
            st.session_state['time_points'] = time_points
    
    # Display results if available
    if 'results' in st.session_state:
        df = st.session_state['results']
        scf_profile = st.session_state['scf_profile']
        time_points = st.session_state['time_points']
        
        # Summary metrics
        st.header("Summary Metrics")
        col1, col2, col3, col4 = st.columns(4)
        
        with col1:
            st.metric("Peak ERK Activation", f"{df['ERK_active'].max():.3f}",
                     delta=f"+{(df['ERK_active'].max() - df['ERK_active'].iloc[0]):.3f}")
        
        with col2:
            st.metric("Peak AKT Activation", f"{df['AKT_active'].max():.3f}",
                     delta=f"+{(df['AKT_active'].max() - df['AKT_active'].iloc[0]):.3f}")
        
        with col3:
            st.metric("Final Proliferation", f"{df['Proliferation'].iloc[-1]:.3f}",
                     delta=f"+{df['Proliferation'].iloc[-1]:.3f}")
        
        with col4:
            st.metric("Final Differentiation", f"{df['Differentiation'].iloc[-1]:.3f}",
                     delta=f"+{df['Differentiation'].iloc[-1]:.3f}")
        
        # Tabs for different visualizations
        tab1, tab2, tab3, tab4, tab5 = st.tabs([
            "Receptor Dynamics", 
            "Signaling Pathways", 
            "Cell Fate",
            "Heatmap",
            "Phase Space"
        ])
        
        with tab1:
            st.plotly_chart(plot_receptor_dynamics(df), use_container_width=True)
            
            # SCF input profile
            fig_scf = go.Figure()
            fig_scf.add_trace(go.Scatter(x=time_points[:-1], y=scf_profile[:-1],
                                        fill='tozeroy', line=dict(color='#e74c3c', width=2)))
            fig_scf.update_layout(title="SCF Input Profile", xaxis_title="Time (min)",
                                 yaxis_title="SCF Concentration", template="plotly_white",
                                 height=300)
            st.plotly_chart(fig_scf, use_container_width=True)
        
        with tab2:
            st.plotly_chart(plot_signaling_pathways(df), use_container_width=True)
            
            # Pathway comparison
            col1, col2 = st.columns(2)
            with col1:
                st.subheader("Pathway Peak Times")
                peak_times = {
                    'PI3K/AKT': df['Time'][df['AKT_active'].idxmax()],
                    'MAPK/ERK': df['Time'][df['ERK_active'].idxmax()],
                    'JAK/STAT': df['Time'][df['STAT_active'].idxmax()]
                }
                for pathway, time in peak_times.items():
                    st.write(f"**{pathway}**: {time:.1f} min")
            
            with col2:
                st.subheader("Pathway Peak Amplitudes")
                peak_amps = {
                    'PI3K/AKT': df['AKT_active'].max(),
                    'MAPK/ERK': df['ERK_active'].max(),
                    'JAK/STAT': df['STAT_active'].max()
                }
                for pathway, amp in peak_amps.items():
                    st.write(f"**{pathway}**: {amp:.3f}")
        
        with tab3:
            st.plotly_chart(plot_cell_fate(df), use_container_width=True)
            
            # Final outcomes
            st.subheader("Cell Fate Analysis")
            col1, col2 = st.columns(2)
            
            with col1:
                fate_df = pd.DataFrame({
                    'Outcome': ['Proliferation', 'Differentiation', 'Survival', 'Apoptosis'],
                    'Level': [
                        df['Proliferation'].iloc[-1],
                        df['Differentiation'].iloc[-1],
                        df['Survival'].iloc[-1],
                        df['Apoptosis'].iloc[-1]
                    ]
                })
                fig_pie = px.pie(fate_df, values='Level', names='Outcome',
                                title='Final Cell Fate Distribution',
                                color_discrete_sequence=px.colors.qualitative.Set2)
                st.plotly_chart(fig_pie, use_container_width=True)
            
            with col2:
                fig_bar = px.bar(fate_df, x='Outcome', y='Level',
                                title='Cell Fate Signal Strengths',
                                color='Outcome',
                                color_discrete_sequence=px.colors.qualitative.Set2)
                st.plotly_chart(fig_bar, use_container_width=True)
        
        with tab4:
            st.plotly_chart(plot_pathway_heatmap(df), use_container_width=True)
            
            # Correlation matrix
            st.subheader("Pathway Correlations")
            pathways = ['PI3K_active', 'AKT_active', 'RAS_active', 'MEK_active', 
                       'ERK_active', 'JAK_active', 'STAT_active']
            corr_matrix = df[pathways].corr()
            
            fig_corr = go.Figure(data=go.Heatmap(
                z=corr_matrix.values,
                x=['PI3K', 'AKT', 'RAS', 'MEK', 'ERK', 'JAK', 'STAT'],
                y=['PI3K', 'AKT', 'RAS', 'MEK', 'ERK', 'JAK', 'STAT'],
                colorscale='RdBu',
                zmid=0,
                text=corr_matrix.values.round(2),
                texttemplate='%{text}',
                textfont={"size": 10}
            ))
            fig_corr.update_layout(title="Pathway Correlation Matrix", height=400)
            st.plotly_chart(fig_corr, use_container_width=True)
        
        with tab5:
            col1, col2 = st.columns(2)
            
            with col1:
                x_var = st.selectbox("X Variable", df.columns[1:-1], index=5)
            with col2:
                y_var = st.selectbox("Y Variable", df.columns[1:-1], index=9)
            
            st.plotly_chart(plot_phase_portrait(df, x_var, y_var), use_container_width=True)
            
            # Multiple phase portraits
            st.subheader("Key Phase Portraits")
            col1, col2 = st.columns(2)
            
            with col1:
                st.plotly_chart(plot_phase_portrait(df, 'AKT_active', 'Proliferation'), 
                              use_container_width=True)
            
            with col2:
                st.plotly_chart(plot_phase_portrait(df, 'ERK_active', 'Differentiation'), 
                              use_container_width=True)
        
        # Data export
        st.header("Export Results")
        col1, col2 = st.columns(2)
        
        with col1:
            csv = df.to_csv(index=False)
            st.download_button(
                label="Download CSV",
                data=csv,
                file_name="ckit_scf_simulation_results.csv",
                mime="text/csv"
            )
        
        with col2:
            st.download_button(
                label="Download JSON",
                data=df.to_json(),
                file_name="ckit_scf_simulation_results.json",
                mime="application/json"
            )
    
    else:
        # Welcome message
        st.info("Configure parameters in the sidebar and click 'Run Simulation' to begin.")
        
        st.header("About c-KIT/SCF Signaling")
        
        col1, col2 = st.columns(2)
        
        with col1:
            st.subheader("Biological Context")
            st.write("""
            The c-KIT receptor tyrosine kinase and its ligand Stem Cell Factor (SCF) 
            play crucial roles in:
            - **Hematopoiesis**: Development of blood cells
            - **Mast cell function**: Immune responses and allergic reactions
            - **Gametogenesis**: Reproductive cell development
            - **Melanogenesis**: Pigmentation
            """)
        
        with col2:
            st.subheader("Key Pathways")
            st.write("""
            This simulation models three major signaling cascades:
            - **PI3K/AKT**: Cell survival and growth
            - **MAPK/ERK**: Proliferation and differentiation
            - **JAK/STAT**: Gene transcription and cell fate
            """)
        
        st.subheader("Model Features")
        st.write("""
        - Receptor binding, dimerization, and phosphorylation
        - Signal transduction through multiple pathways
        - Cell fate outcome integration
        - Dynamic SCF stimulation profiles
        - Comprehensive visualization suite
        """)


if __name__ == "__main__":
    main()
