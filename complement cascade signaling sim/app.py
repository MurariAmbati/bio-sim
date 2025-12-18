"""
Streamlit visualization app for complement cascade simulation
Interactive exploration of opsonization, lysis, and inflammation
"""

import streamlit as st
import numpy as np
import pandas as pd
import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots
from complement_cascade import ComplementCascade, ComplementPathway


# Page configuration
st.set_page_config(
    page_title="complement cascade signaling",
    page_icon="⚡",
    layout="wide",
    initial_sidebar_state="expanded"
)

# Custom CSS for minimalistic design
st.markdown("""
    <style>
    .main {background-color: #ffffff;}
    h1 {color: #2c3e50; font-weight: 300; letter-spacing: -1px;}
    h2 {color: #34495e; font-weight: 300;}
    h3 {color: #7f8c8d; font-weight: 300;}
    .stMetric {background-color: #f8f9fa; padding: 10px; border-radius: 5px;}
    .stMetric label {color: #000000 !important;}
    .stMetric [data-testid="stMetricValue"] {color: #000000 !important;}
    </style>
""", unsafe_allow_html=True)


def create_pathway_diagram():
    """Create interactive pathway diagram"""
    fig = go.Figure()
    
    # Classical pathway (left)
    fig.add_trace(go.Scatter(
        x=[1, 1, 1, 2], y=[5, 4, 3, 2],
        mode='markers+text+lines',
        name='classical',
        text=['C1qrs', 'C4', 'C2', 'C4b2a'],
        textposition='middle right',
        marker=dict(size=20, color='#3498db'),
        line=dict(color='#3498db', width=2)
    ))
    
    # Alternative pathway (middle)
    fig.add_trace(go.Scatter(
        x=[3, 3, 3, 3], y=[5, 4, 3, 2],
        mode='markers+text+lines',
        name='alternative',
        text=['C3(H2O)', 'Factor B', 'Factor D', 'C3bBb'],
        textposition='middle center',
        marker=dict(size=20, color='#e74c3c'),
        line=dict(color='#e74c3c', width=2)
    ))
    
    # Lectin pathway (right)
    fig.add_trace(go.Scatter(
        x=[5, 5, 5, 4], y=[5, 4, 3, 2],
        mode='markers+text+lines',
        name='lectin',
        text=['MBL', 'MASP-2', 'C4/C2', 'C4b2a'],
        textposition='middle left',
        marker=dict(size=20, color='#2ecc71'),
        line=dict(color='#2ecc71', width=2)
    ))
    
    # C3 convertase convergence
    fig.add_trace(go.Scatter(
        x=[3], y=[1],
        mode='markers+text',
        name='C3 convertase',
        text=['C3 → C3b + C3a'],
        textposition='bottom center',
        marker=dict(size=30, color='#9b59b6', symbol='star')
    ))
    
    # Terminal pathway
    fig.add_trace(go.Scatter(
        x=[3, 3, 3], y=[0, -1, -2],
        mode='markers+text+lines',
        name='terminal',
        text=['C5 convertase', 'C5b-6-7-8', 'MAC (C9n)'],
        textposition='bottom center',
        marker=dict(size=20, color='#e67e22'),
        line=dict(color='#e67e22', width=3)
    ))
    
    fig.update_layout(
        title='complement cascade pathways',
        showlegend=True,
        height=500,
        xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
        yaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
        plot_bgcolor='white',
        font=dict(size=10, color='black')
    )
    
    return fig


def plot_time_series(df: pd.DataFrame):
    """Plot concentration time series"""
    fig = make_subplots(
        rows=2, cols=2,
        subplot_titles=('complement proteins', 'anaphylatoxins', 'effector outcomes', 'inflammation'),
        vertical_spacing=0.15,
        horizontal_spacing=0.1
    )
    
    # Complement proteins
    fig.add_trace(go.Scatter(x=df['time'], y=df['C3'], name='C3', line=dict(color='#3498db')), row=1, col=1)
    fig.add_trace(go.Scatter(x=df['time'], y=df['C3b'], name='C3b', line=dict(color='#e74c3c')), row=1, col=1)
    fig.add_trace(go.Scatter(x=df['time'], y=df['C5'], name='C5', line=dict(color='#2ecc71')), row=1, col=1)
    fig.add_trace(go.Scatter(x=df['time'], y=df['MAC'], name='MAC', line=dict(color='#9b59b6')), row=1, col=1)
    
    # Anaphylatoxins
    fig.add_trace(go.Scatter(x=df['time'], y=df['C3a'], name='C3a', line=dict(color='#f39c12')), row=1, col=2)
    fig.add_trace(go.Scatter(x=df['time'], y=df['C5a'], name='C5a', line=dict(color='#e74c3c', width=3)), row=1, col=2)
    
    # Effector outcomes
    fig.add_trace(go.Scatter(x=df['time'], y=df['opsonized_count'], name='opsonized', 
                            line=dict(color='#3498db'), fill='tozeroy'), row=2, col=1)
    fig.add_trace(go.Scatter(x=df['time'], y=df['lysed_count'], name='lysed', 
                            line=dict(color='#e74c3c'), fill='tozeroy'), row=2, col=1)
    
    # Inflammation
    fig.add_trace(go.Scatter(x=df['time'], y=df['inflammation'], name='inflammation', 
                            line=dict(color='#e67e22', width=3), fill='tozeroy'), row=2, col=2)
    
    fig.update_xaxes(title_text="time", row=2, col=1)
    fig.update_xaxes(title_text="time", row=2, col=2)
    fig.update_yaxes(title_text="concentration (μg/mL)", row=1, col=1)
    fig.update_yaxes(title_text="concentration (μg/mL)", row=1, col=2)
    fig.update_yaxes(title_text="pathogen count", row=2, col=1)
    fig.update_yaxes(title_text="response level", row=2, col=2)
    
    fig.update_layout(
        height=700,
        showlegend=True,
        template='plotly_white',
        font=dict(size=10)
    )
    
    return fig


def plot_pathogen_status(sim: ComplementCascade):
    """Visualize pathogen opsonization and lysis status"""
    pathogens_data = []
    for p in sim.state.pathogens:
        pathogens_data.append({
            'pathogen': p.pathogen_id,
            'C3b deposited': p.c3b_deposited,
            'MAC complexes': p.mac_complexes,
            'opsonized': 'yes' if p.opsonized else 'no',
            'lysed': 'yes' if p.lysed else 'no'
        })
    
    df_pathogens = pd.DataFrame(pathogens_data)
    
    fig = make_subplots(
        rows=1, cols=2,
        subplot_titles=('C3b deposition', 'MAC insertion'),
        specs=[[{"type": "bar"}, {"type": "bar"}]]
    )
    
    fig.add_trace(
        go.Bar(x=df_pathogens['pathogen'], y=df_pathogens['C3b deposited'], 
               name='C3b', marker_color='#3498db'),
        row=1, col=1
    )
    
    fig.add_trace(
        go.Bar(x=df_pathogens['pathogen'], y=df_pathogens['MAC complexes'], 
               name='MAC', marker_color='#e74c3c'),
        row=1, col=2
    )
    
    fig.update_layout(height=400, showlegend=False, template='plotly_white')
    fig.update_yaxes(title_text="count", row=1, col=1)
    fig.update_yaxes(title_text="count", row=1, col=2)
    
    return fig, df_pathogens


def plot_component_heatmap(df: pd.DataFrame):
    """Create heatmap of component concentrations over time"""
    components = ['C3', 'C3b', 'C5', 'MAC', 'C3a', 'C5a']
    
    # Normalize each component to 0-1 range for better visualization
    heatmap_data = []
    for comp in components:
        values = df[comp].values
        if values.max() > 0:
            normalized = (values - values.min()) / (values.max() - values.min())
        else:
            normalized = values
        heatmap_data.append(normalized)
    
    fig = go.Figure(data=go.Heatmap(
        z=heatmap_data,
        x=df['time'],
        y=components,
        colorscale='RdYlBu_r',
        colorbar=dict(title='normalized<br>concentration')
    ))
    
    fig.update_layout(
        title='component activity over time',
        xaxis_title='time',
        yaxis_title='component',
        height=400,
        template='plotly_white'
    )
    
    return fig


def plot_3d_phase_space(df: pd.DataFrame):
    """3D phase space plot"""
    fig = go.Figure(data=[go.Scatter3d(
        x=df['C3b'],
        y=df['MAC'],
        z=df['inflammation'],
        mode='lines+markers',
        marker=dict(
            size=3,
            color=df['time'],
            colorscale='Viridis',
            showscale=True,
            colorbar=dict(title='time')
        ),
        line=dict(color='darkgray', width=2)
    )])
    
    fig.update_layout(
        title='phase space: opsonization → lysis → inflammation',
        scene=dict(
            xaxis_title='C3b (opsonization)',
            yaxis_title='MAC (lysis)',
            zaxis_title='inflammation'
        ),
        height=600,
        template='plotly_white'
    )
    
    return fig


def main():
    # Header
    st.title("complement cascade signaling")
    st.markdown("*opsonization • lysis • inflammation*")
    
    # Sidebar controls
    st.sidebar.header("simulation parameters")
    
    duration = st.sidebar.slider("duration", min_value=10.0, max_value=200.0, value=50.0, step=10.0)
    n_pathogens = st.sidebar.slider("pathogen count", min_value=1, max_value=20, value=5)
    pathogen_size = st.sidebar.slider("pathogen surface area", min_value=1.0, max_value=50.0, value=10.0)
    
    st.sidebar.subheader("pathway activation")
    activate_classical = st.sidebar.checkbox("classical pathway", value=True)
    classical_intensity = st.sidebar.slider("classical intensity", 0.0, 3.0, 1.5, 0.1) if activate_classical else 0.0
    
    activate_alternative = st.sidebar.checkbox("alternative pathway", value=True)
    
    activate_lectin = st.sidebar.checkbox("lectin pathway", value=False)
    lectin_intensity = st.sidebar.slider("lectin intensity", 0.0, 3.0, 1.0, 0.1) if activate_lectin else 0.0
    
    # Run simulation button
    if st.sidebar.button("run simulation", type="primary"):
        with st.spinner("simulating complement cascade..."):
            # Initialize simulation
            sim = ComplementCascade()
            
            # Add pathogens
            for i in range(n_pathogens):
                sim.add_pathogen(f"pathogen_{i}", surface_area=pathogen_size)
            
            # Activate pathways
            if activate_classical:
                sim.activate_classical_pathway(intensity=classical_intensity)
            if activate_lectin:
                sim.activate_lectin_pathway(intensity=lectin_intensity)
            if activate_alternative:
                sim.activate_alternative_pathway()
            
            # Run simulation
            sim.run(duration=duration)
            
            # Store in session state
            st.session_state['sim'] = sim
            st.session_state['df'] = sim.get_history_dataframe()
    
    # Display pathway diagram
    st.subheader("cascade architecture")
    st.plotly_chart(create_pathway_diagram(), use_container_width=True)
    
    # Display results if simulation has been run
    if 'sim' in st.session_state and 'df' in st.session_state:
        sim = st.session_state['sim']
        df = st.session_state['df']
        
        st.divider()
        
        # Metrics
        st.subheader("outcomes")
        col1, col2, col3, col4, col5 = st.columns(5)
        
        with col1:
            st.metric("C3b peak", f"{df['C3b'].max():.1f}")
        with col2:
            st.metric("MAC formed", f"{df['MAC'].iloc[-1]:.1f}")
        with col3:
            st.metric("opsonized", f"{int(df['opsonized_count'].iloc[-1])}/{n_pathogens}")
        with col4:
            st.metric("lysed", f"{int(df['lysed_count'].iloc[-1])}/{n_pathogens}")
        with col5:
            st.metric("inflammation", f"{df['inflammation'].iloc[-1]:.1f}")
        
        st.divider()
        
        # Time series plots
        st.subheader("temporal dynamics")
        st.plotly_chart(plot_time_series(df), use_container_width=True)
        
        st.divider()
        
        # Pathogen status
        st.subheader("pathogen targeting")
        fig_pathogen, df_pathogens = plot_pathogen_status(sim)
        
        col1, col2 = st.columns([2, 1])
        with col1:
            st.plotly_chart(fig_pathogen, use_container_width=True)
        with col2:
            st.dataframe(df_pathogens, use_container_width=True, hide_index=True)
        
        st.divider()
        
        # Component heatmap
        st.subheader("component activity")
        st.plotly_chart(plot_component_heatmap(df), use_container_width=True)
        
        st.divider()
        
        # 3D phase space
        st.subheader("phase space trajectory")
        st.plotly_chart(plot_3d_phase_space(df), use_container_width=True)
        
        st.divider()
        
        # Data export
        st.subheader("export data")
        csv = df.to_csv(index=False)
        st.download_button(
            label="download simulation data (csv)",
            data=csv,
            file_name="complement_cascade_simulation.csv",
            mime="text/csv"
        )
        
        # Raw data viewer
        with st.expander("view raw data"):
            st.dataframe(df, use_container_width=True)
    
    else:
        st.info("configure parameters and click 'run simulation' to begin")
    
    # Footer
    st.divider()
    st.markdown("""
        <div style='text-align: center; color: #7f8c8d; font-size: 0.9em;'>
        complement cascade signaling simulation • classical • alternative • lectin pathways
        </div>
    """, unsafe_allow_html=True)


if __name__ == "__main__":
    main()
