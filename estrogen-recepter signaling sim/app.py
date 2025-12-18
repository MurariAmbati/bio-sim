"""
streamlit web application for estrogen receptor signaling simulation
"""

import streamlit as st
import pandas as pd
import numpy as np
import plotly.graph_objects as go
from models import LigandLibrary, CellType
from simulation import (
    ERSimulation, 
    DoseResponseSimulation, 
    TimeCourseSimulation,
    TissueSpecificSimulation
)
from visualization import PathwayVisualizer, NetworkVisualizer

# page config
st.set_page_config(
    page_title="estrogen receptor signaling",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded"
)

# custom css
st.markdown("""
<style>
    .main-header {
        font-size: 2.5rem;
        font-weight: 600;
        color: #2c3e50;
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
        border-left: 4px solid #e74c3c;
    }
</style>
""", unsafe_allow_html=True)


def main():
    """main application"""
    
    # header
    st.markdown('<p class="main-header">estrogen receptor signaling simulation</p>', 
                unsafe_allow_html=True)
    st.markdown('<p class="sub-header">comprehensive mechanistic model of er-α and er-β signaling pathways</p>', 
                unsafe_allow_html=True)
    
    # sidebar navigation
    st.sidebar.title("navigation")
    page = st.sidebar.radio(
        "select analysis",
        ["pathway overview", "time course simulation", "dose-response", 
         "tissue specificity", "ligand comparison", "network visualization"]
    )
    
    if page == "pathway overview":
        show_pathway_overview()
    elif page == "time course simulation":
        show_timecourse()
    elif page == "dose-response":
        show_dose_response()
    elif page == "tissue specificity":
        show_tissue_specificity()
    elif page == "ligand comparison":
        show_ligand_comparison()
    elif page == "network visualization":
        show_network()


def show_pathway_overview():
    """display pathway overview and biological context"""
    
    st.header("estrogen receptor signaling pathways")
    
    col1, col2 = st.columns(2)
    
    with col1:
        st.subheader("genomic pathway")
        st.markdown("""
        **classical nuclear signaling (hours timescale)**
        
        1. **ligand binding**: estradiol binds to er-α or er-β
        2. **receptor dimerization**: formation of homodimers or heterodimers
        3. **nuclear translocation**: receptor complex enters nucleus
        4. **dna binding**: dimers bind to estrogen response elements (eres)
        5. **transcription**: recruitment of coactivators and rna polymerase ii
        6. **target gene expression**: synthesis of proteins regulating:
           - cell proliferation (cyclin d1, c-myc)
           - survival (bcl-2)
           - differentiation (progesterone receptor)
           - metabolism (tff1, vegf)
        
        **key features:**
        - high specificity via ere sequences
        - integration with chromatin remodeling
        - tissue-specific coregulator expression
        - feedback regulation
        """)
    
    with col2:
        st.subheader("non-genomic pathway")
        st.markdown("""
        **rapid membrane signaling (minutes timescale)**
        
        1. **membrane er activation**: gpr30/gper and membrane-localized ers
        2. **kinase cascade activation**:
           - mapk/erk pathway → cell proliferation
           - pi3k/akt pathway → cell survival
           - src kinase → multiple downstream effects
        3. **calcium signaling**: rapid ca²⁺ release from er stores
        4. **nitric oxide production**: enos activation via akt
        5. **crosstalk with genomic pathway**: kinases phosphorylate nuclear ers
        
        **key features:**
        - rapid onset (seconds to minutes)
        - amplification via kinase cascades
        - integration with growth factor signaling
        - modulation of genomic pathway
        """)
    
    st.markdown("---")
    
    # receptor subtypes
    st.subheader("receptor subtypes: er-α vs er-β")
    
    col1, col2, col3 = st.columns(3)
    
    with col1:
        st.metric("er-α", "595 amino acids", "chromosome 6")
        st.markdown("""
        **primary locations:**
        - uterus
        - breast
        - ovary
        - bone
        
        **functions:**
        - proliferation
        - reproduction
        - bone maintenance
        """)
    
    with col2:
        st.metric("er-β", "530 amino acids", "chromosome 14")
        st.markdown("""
        **primary locations:**
        - ovary
        - prostate
        - cardiovascular
        - cns
        
        **functions:**
        - anti-proliferative
        - differentiation
        - neuroprotection
        """)
    
    with col3:
        st.metric("homology", "97%", "dna binding domain")
        st.markdown("""
        **key differences:**
        - ligand binding affinity
        - tissue distribution
        - coregulator preference
        - target gene selectivity
        
        **ratio effects:**
        - α:β balance determines response
        - β can antagonize α activity
        """)
    
    st.markdown("---")
    
    # ligand classes
    st.subheader("ligand pharmacology")
    
    ligand_data = {
        'ligand': ['estradiol', 'tamoxifen', 'raloxifene', 'fulvestrant', 'genistein'],
        'class': ['agonist', 'serm', 'serm', 'serd', 'phytoestrogen'],
        'er-α affinity (nm)': [0.1, 10, 5, 0.5, 100],
        'er-β affinity (nm)': [0.5, 20, 10, 1.0, 20],
        'breast effect': ['agonist', 'antagonist', 'antagonist', 'antagonist', 'partial agonist'],
        'bone effect': ['agonist', 'agonist', 'agonist', 'antagonist', 'agonist'],
    }
    
    df_ligands = pd.DataFrame(ligand_data)
    st.dataframe(df_ligands, use_container_width=True)
    
    st.info("""
    **serms (selective estrogen receptor modulators)**: tissue-specific agonist/antagonist activity
    
    **serds (selective estrogen receptor degraders)**: promote receptor degradation, pure antagonists
    
    **phytoestrogens**: plant-derived compounds with weak estrogenic activity, er-β selective
    """)


def show_timecourse():
    """time course simulation interface"""
    
    st.header("time course simulation")
    
    # parameters
    col1, col2 = st.columns([1, 2])
    
    with col1:
        st.subheader("simulation parameters")
        
        # cell type
        cell_type = st.selectbox(
            "cell type",
            [CellType.BREAST_CANCER_MCF7, CellType.BREAST_CANCER_T47D,
             CellType.OSTEOBLAST, CellType.ENDOMETRIAL, CellType.ENDOTHELIAL],
            format_func=lambda x: x.replace('_', ' ').title()
        )
        
        # ligand
        ligands = LigandLibrary.list_ligands()
        ligand = st.selectbox("ligand", ligands, format_func=str.title)
        
        # dose
        dose = st.slider("dose (mg)", 0.01, 100.0, 10.0, step=0.01,
                        format="%.2f", help="oral dose in milligrams")
        
        # duration
        duration = st.slider("simulation duration (hours)", 1, 96, 48)
        
        # run button
        run_sim = st.button("run simulation", type="primary", use_container_width=True)
        
        if run_sim:
            st.session_state['run_timecourse'] = True
            st.session_state['tc_params'] = {
                'cell_type': cell_type,
                'ligand': ligand,
                'dose': dose,
                'duration': duration
            }
    
    with col2:
        if st.session_state.get('run_timecourse'):
            params = st.session_state['tc_params']
            
            with st.spinner("running simulation..."):
                # run simulation
                sim = ERSimulation(params['cell_type'], params['duration'])
                sim.add_treatment(params['ligand'], params['dose'])
                df = sim.run()
                
                # store results
                st.session_state['tc_results'] = df
                st.session_state['tc_sim'] = sim
            
            st.success("simulation completed!")
    
    # display results
    if 'tc_results' in st.session_state:
        df = st.session_state['tc_results']
        
        st.markdown("---")
        
        # summary metrics
        st.subheader("summary statistics")
        
        stats = st.session_state['tc_sim'].get_summary_statistics(df)
        
        col1, col2, col3, col4 = st.columns(4)
        col1.metric("er-α peak (nm)", f"{stats['er_alpha_max']:.2f}")
        col2.metric("transcription auc", f"{stats['transcription_auc']:.1f}")
        col3.metric("cell divisions", int(stats['final_divisions']))
        col4.metric("survival signal", f"{stats['survival_signal_mean']:.3f}")
        
        st.markdown("---")
        
        # visualizations
        tabs = st.tabs(["receptors", "target genes", "signaling", "cellular outcomes"])
        
        with tabs[0]:
            fig = PathwayVisualizer.plot_receptor_dynamics(df)
            st.plotly_chart(fig, use_container_width=True)
        
        with tabs[1]:
            fig = PathwayVisualizer.plot_target_genes(df)
            st.plotly_chart(fig, use_container_width=True)
        
        with tabs[2]:
            fig = PathwayVisualizer.plot_signaling_cascades(df)
            st.plotly_chart(fig, use_container_width=True)
        
        with tabs[3]:
            fig = PathwayVisualizer.plot_cellular_outcomes(df)
            st.plotly_chart(fig, use_container_width=True)
        
        # download data
        st.markdown("---")
        csv = df.to_csv(index=False)
        st.download_button(
            label="download simulation data (csv)",
            data=csv,
            file_name=f"er_simulation_{params['ligand']}_{params['dose']}mg.csv",
            mime="text/csv"
        )


def show_dose_response():
    """dose-response analysis"""
    
    st.header("dose-response analysis")
    
    col1, col2 = st.columns([1, 2])
    
    with col1:
        st.subheader("parameters")
        
        cell_type = st.selectbox(
            "cell type",
            [CellType.BREAST_CANCER_MCF7, CellType.BREAST_CANCER_T47D],
            format_func=lambda x: x.replace('_', ' ').title(),
            key='dr_cell'
        )
        
        ligand = st.selectbox(
            "ligand",
            LigandLibrary.list_ligands(),
            format_func=str.title,
            key='dr_ligand'
        )
        
        st.markdown("**dose range**")
        dose_min = st.number_input("minimum (mg)", 0.001, 10.0, 0.01, format="%.3f")
        dose_max = st.number_input("maximum (mg)", 1.0, 1000.0, 100.0)
        n_doses = st.slider("number of doses", 5, 20, 10)
        
        duration = st.slider("duration (hours)", 12, 72, 24, key='dr_duration')
        
        if st.button("run dose-response", type="primary", use_container_width=True):
            doses = np.logspace(np.log10(dose_min), np.log10(dose_max), n_doses)
            
            with st.spinner("running simulations..."):
                dr_sim = DoseResponseSimulation(cell_type)
                dr_df = dr_sim.run_dose_response(ligand, doses, duration)
                
                st.session_state['dr_results'] = dr_df
                st.session_state['dr_sim'] = dr_sim
            
            st.success("analysis completed!")
    
    with col2:
        if 'dr_results' in st.session_state:
            dr_df = st.session_state['dr_results']
            
            # response metric selection
            response_options = [
                'transcription_auc', 'er_alpha_mean', 'cyclin_d1_peak',
                'proliferation_rate_mean', 'mapk_peak', 'akt_peak'
            ]
            
            response = st.selectbox(
                "response metric",
                response_options,
                format_func=lambda x: x.replace('_', ' ').title()
            )
            
            # plot dose-response curve
            fig = PathwayVisualizer.plot_dose_response(dr_df, response)
            st.plotly_chart(fig, use_container_width=True)
            
            # calculate ec50
            ec50 = st.session_state['dr_sim'].calculate_ec50(dr_df, response)
            st.metric("ec50", f"{ec50:.3f} mg")
            
            # show data table
            st.dataframe(dr_df, use_container_width=True)


def show_tissue_specificity():
    """tissue-specific effects analysis"""
    
    st.header("tissue specificity analysis")
    
    st.markdown("""
    analyze how serms (selective estrogen receptor modulators) exhibit different effects
    across tissue types due to variations in:
    - er-α/er-β expression ratios
    - coregulator expression profiles
    - chromatin accessibility
    - cellular context
    """)
    
    col1, col2 = st.columns([1, 2])
    
    with col1:
        ligand = st.selectbox(
            "select ligand",
            ['tamoxifen', 'raloxifene', 'estradiol'],
            format_func=str.title,
            key='tissue_ligand'
        )
        
        dose = st.slider("dose (mg)", 0.1, 50.0, 10.0, key='tissue_dose')
        
        tissues = st.multiselect(
            "tissue types",
            [CellType.BREAST_CANCER_MCF7, CellType.OSTEOBLAST, 
             CellType.ENDOMETRIAL, CellType.ENDOTHELIAL],
            default=[CellType.BREAST_CANCER_MCF7, CellType.OSTEOBLAST],
            format_func=lambda x: x.replace('_', ' ').title()
        )
        
        if st.button("run analysis", type="primary", use_container_width=True):
            if len(tissues) < 2:
                st.error("select at least 2 tissue types")
            else:
                with st.spinner("simulating across tissues..."):
                    tissue_sim = TissueSpecificSimulation()
                    results = tissue_sim.run_tissue_comparison(ligand, dose, tissues)
                    
                    st.session_state['tissue_results'] = results
                    st.session_state['tissue_sim'] = tissue_sim
                
                st.success("analysis completed!")
    
    with col2:
        if 'tissue_results' in st.session_state:
            results = st.session_state['tissue_results']
            
            # selectivity metrics
            st.subheader("tissue selectivity indices")
            
            metrics = ['er_alpha_transcriptional_activity', 'gene_cyclin_d1_protein',
                      'nongenomic_mapk_active', 'cell_proliferation_rate']
            
            selected_metric = st.selectbox(
                "metric for selectivity",
                metrics,
                format_func=lambda x: x.replace('_', ' ').replace('gene ', '').title()
            )
            
            selectivity = st.session_state['tissue_sim'].calculate_selectivity_index(
                results, selected_metric, list(results.keys())[0]
            )
            
            # display selectivity
            for tissue, index in selectivity.items():
                st.metric(
                    tissue.replace('_', ' ').title(),
                    f"{index:.2f}x",
                    f"{'agonist' if index > 0.5 else 'antagonist'} activity"
                )
            
            st.markdown("---")
            
            # heatmap comparison
            st.subheader("temporal dynamics across tissues")
            
            fig = PathwayVisualizer.plot_heatmap(
                results, 
                selected_metric,
                f"{selected_metric} across tissues"
            )
            st.plotly_chart(fig, use_container_width=True)


def show_ligand_comparison():
    """compare multiple ligands"""
    
    st.header("ligand comparison")
    
    col1, col2 = st.columns([1, 2])
    
    with col1:
        cell_type = st.selectbox(
            "cell type",
            [CellType.BREAST_CANCER_MCF7, CellType.OSTEOBLAST],
            format_func=lambda x: x.replace('_', ' ').title(),
            key='comp_cell'
        )
        
        available_ligands = LigandLibrary.list_ligands()
        selected_ligands = st.multiselect(
            "select ligands to compare",
            available_ligands,
            default=['estradiol', 'tamoxifen', 'raloxifene'],
            format_func=str.title
        )
        
        dose = st.slider("dose (mg)", 0.1, 50.0, 10.0, key='comp_dose')
        duration = st.slider("duration (hours)", 12, 72, 24, key='comp_duration')
        
        if st.button("compare ligands", type="primary", use_container_width=True):
            if len(selected_ligands) < 2:
                st.error("select at least 2 ligands")
            else:
                with st.spinner("running comparisons..."):
                    tc_sim = TimeCourseSimulation(cell_type)
                    results = tc_sim.compare_ligands(selected_ligands, dose, duration)
                    
                    st.session_state['comp_results'] = results
                
                st.success("comparison completed!")
    
    with col2:
        if 'comp_results' in st.session_state:
            results = st.session_state['comp_results']
            
            # metric selection
            metrics = [
                'er_alpha_transcriptional_activity',
                'gene_cyclin_d1_protein',
                'gene_bcl2_protein',
                'nongenomic_mapk_active',
                'cell_proliferation_rate'
            ]
            
            metric = st.selectbox(
                "select metric",
                metrics,
                format_func=lambda x: x.replace('_', ' ').replace('gene ', '').title()
            )
            
            # create comparison plot
            fig = go.Figure()
            
            colors = ['#e74c3c', '#3498db', '#2ecc71', '#f39c12', '#9b59b6']
            
            for i, (ligand, df) in enumerate(results.items()):
                if metric in df.columns:
                    fig.add_trace(go.Scatter(
                        x=df['time_hours'],
                        y=df[metric],
                        name=ligand.title(),
                        line=dict(color=colors[i % len(colors)], width=2)
                    ))
            
            fig.update_layout(
                title=f"{metric.replace('_', ' ').replace('gene ', '').title()} comparison",
                xaxis_title="time (hours)",
                yaxis_title=metric.replace('_', ' ').replace('gene ', ''),
                height=500
            )
            
            st.plotly_chart(fig, use_container_width=True)
            
            # summary table
            st.subheader("peak responses")
            
            summary_data = []
            for ligand, df in results.items():
                if metric in df.columns:
                    summary_data.append({
                        'ligand': ligand.title(),
                        'peak': df[metric].max(),
                        'auc': np.trapz(df[metric], df['time_hours']),
                        'time_to_peak': df.loc[df[metric].idxmax(), 'time_hours']
                    })
            
            summary_df = pd.DataFrame(summary_data)
            st.dataframe(summary_df, use_container_width=True)


def show_network():
    """display pathway network"""
    
    st.header("signaling network visualization")
    
    st.markdown("""
    interactive network diagram showing major components of estrogen receptor signaling:
    - ligand-receptor interactions
    - genomic and non-genomic pathways
    - signaling cascades
    - cellular outcomes
    """)
    
    fig = NetworkVisualizer.create_pathway_network()
    st.plotly_chart(fig, use_container_width=True)
    
    # pathway descriptions
    col1, col2 = st.columns(2)
    
    with col1:
        st.subheader("genomic pathway components")
        with st.expander("estrogen response elements (eres)"):
            st.markdown("""
            **consensus sequence**: GGTCAnnnTGACC
            
            palindromic dna motifs recognized by er dimers. variations in sequence
            and spacing affect binding affinity and transcriptional output.
            """)
        
        with st.expander("coregulators"):
            st.markdown("""
            **coactivators**:
            - src-1, src-2, src-3 (p160 family)
            - cbp/p300 (histone acetyltransferases)
            - mediator complex
            
            **corepressors**:
            - ncor1, ncor2 (smrt)
            - recruited by antagonist-bound er
            - recruit histone deacetylases (hdacs)
            """)
    
    with col2:
        st.subheader("non-genomic pathway components")
        with st.expander("membrane receptors"):
            st.markdown("""
            **gpr30/gper**: g-protein coupled estrogen receptor
            - distinct from classical ers
            - activates camp and calcium signaling
            - mediates rapid estrogen effects
            
            **membrane-localized er-α/β**:
            - palmitoylation targets ers to membrane
            - interact with growth factor receptors
            - activate kinase cascades
            """)
        
        with st.expander("kinase cascades"):
            st.markdown("""
            **mapk/erk pathway**:
            ras → raf → mek → erk → transcription factors
            
            **pi3k/akt pathway**:
            pi3k → pip3 → pdk1 → akt → mtor, foxo
            
            **crosstalk**:
            - kinases phosphorylate nuclear ers
            - modulates ligand-independent activation
            - integrates er with growth signals
            """)


if __name__ == "__main__":
    main()
