import streamlit as st
import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import plotly.express as px
from thyroid_receptor import (
    ThyroidReceptorModel, ThyroidReceptorParameters,
    simulate_hypothyroidism, simulate_hyperthyroidism,
    simulate_tr_mutation, simulate_development
)


st.set_page_config(page_title="thyroid receptor signaling", layout="wide")


def create_time_series_plot(t, solution, indices, labels, title):
    """create multi-trace time series plot"""
    fig = go.Figure()
    
    colors = px.colors.qualitative.Set2
    
    for idx, label in zip(indices, labels):
        fig.add_trace(go.Scatter(
            x=t, y=solution[:, idx],
            mode='lines',
            name=label,
            line=dict(width=2)
        ))
    
    fig.update_layout(
        title=dict(text=title, font=dict(size=14)),
        xaxis_title="time (hours)",
        yaxis_title="concentration (au)",
        template="plotly_white",
        hovermode='x unified',
        showlegend=True,
        height=400,
        margin=dict(l=50, r=50, t=50, b=50)
    )
    
    return fig


def create_phase_plot(solution, idx1, idx2, label1, label2, title):
    """create phase space plot"""
    fig = go.Figure()
    
    fig.add_trace(go.Scatter(
        x=solution[:, idx1],
        y=solution[:, idx2],
        mode='lines',
        line=dict(width=2, color='#2E86AB'),
        showlegend=False
    ))
    
    # mark start and end
    fig.add_trace(go.Scatter(
        x=[solution[0, idx1]],
        y=[solution[0, idx2]],
        mode='markers',
        marker=dict(size=10, color='green', symbol='circle'),
        name='start'
    ))
    
    fig.add_trace(go.Scatter(
        x=[solution[-1, idx1]],
        y=[solution[-1, idx2]],
        mode='markers',
        marker=dict(size=10, color='red', symbol='square'),
        name='end'
    ))
    
    fig.update_layout(
        title=dict(text=title, font=dict(size=14)),
        xaxis_title=label1,
        yaxis_title=label2,
        template="plotly_white",
        height=400,
        margin=dict(l=50, r=50, t=50, b=50)
    )
    
    return fig


def create_heatmap(data, x_labels, y_labels, title):
    """create heatmap visualization"""
    fig = go.Figure(data=go.Heatmap(
        z=data,
        x=x_labels,
        y=y_labels,
        colorscale='Viridis',
        hoverongaps=False
    ))
    
    fig.update_layout(
        title=dict(text=title, font=dict(size=14)),
        template="plotly_white",
        height=400,
        margin=dict(l=50, r=50, t=50, b=50)
    )
    
    return fig


def create_comparison_plot(data_dict, variable_idx, ylabel, title):
    """create comparison plot across conditions"""
    fig = go.Figure()
    
    for label, (t, solution) in data_dict.items():
        fig.add_trace(go.Scatter(
            x=t, y=solution[:, variable_idx],
            mode='lines',
            name=label,
            line=dict(width=2)
        ))
    
    fig.update_layout(
        title=dict(text=title, font=dict(size=14)),
        xaxis_title="time (hours)",
        yaxis_title=ylabel,
        template="plotly_white",
        hovermode='x unified',
        height=400,
        margin=dict(l=50, r=50, t=50, b=50)
    )
    
    return fig


def main():
    st.title("thyroid hormone receptor signaling")
    
    # sidebar controls
    st.sidebar.header("simulation parameters")
    
    sim_mode = st.sidebar.selectbox(
        "simulation mode",
        ["normal physiology", "hypothyroidism", "hyperthyroidism", 
         "receptor mutation", "development", "parameter exploration"]
    )
    
    # parameter controls
    st.sidebar.subheader("model parameters")
    
    k_syn_t4 = st.sidebar.slider("t4 synthesis rate", 0.1, 2.0, 0.8, 0.1)
    k_syn_t3 = st.sidebar.slider("t3 synthesis rate", 0.05, 1.0, 0.2, 0.05)
    k_d2 = st.sidebar.slider("type 2 deiodinase activity", 0.1, 1.5, 0.6, 0.1)
    k_bind = st.sidebar.slider("receptor binding rate", 0.2, 2.0, 1.2, 0.1)
    k_trans = st.sidebar.slider("transcription activation", 0.5, 4.0, 2.0, 0.1)
    
    sim_time = st.sidebar.slider("simulation time (hours)", 50, 500, 200, 50)
    
    # create model with custom parameters
    params = ThyroidReceptorParameters()
    params.k_syn_t4 = k_syn_t4
    params.k_syn_t3 = k_syn_t3
    params.k_d2 = k_d2
    params.k_bind_t3 = k_bind
    params.k_trans_act = k_trans
    
    model = ThyroidReceptorModel(params)
    
    # run simulation based on mode
    if sim_mode == "normal physiology":
        st.header("normal thyroid hormone signaling")
        
        t, solution = model.simulate((0, sim_time))
        
        col1, col2 = st.columns(2)
        
        with col1:
            fig1 = create_time_series_plot(
                t, solution, [1, 2, 3], ['t4', 't3', 'rt3'],
                "thyroid hormone levels"
            )
            st.plotly_chart(fig1, use_container_width=True)
            
            fig3 = create_time_series_plot(
                t, solution, [4, 5, 6, 7],
                ['tr free', 'tr-t3', 'tr-t3-coactivator', 'tr-corepressor'],
                "receptor states"
            )
            st.plotly_chart(fig3, use_container_width=True)
        
        with col2:
            fig2 = create_time_series_plot(
                t, solution, [0], ['tsh'],
                "tsh levels"
            )
            st.plotly_chart(fig2, use_container_width=True)
            
            fig4 = create_time_series_plot(
                t, solution, [9, 11],
                ['metabolic enzymes', 'metabolic rate'],
                "metabolic effects"
            )
            st.plotly_chart(fig4, use_container_width=True)
        
        # phase plots
        st.subheader("phase space dynamics")
        col3, col4 = st.columns(2)
        
        with col3:
            fig5 = create_phase_plot(
                solution, 2, 5, 't3', 'tr-t3 complex',
                "t3 vs receptor binding"
            )
            st.plotly_chart(fig5, use_container_width=True)
        
        with col4:
            fig6 = create_phase_plot(
                solution, 6, 8, 'activated receptor', 'target mrna',
                "receptor activation vs transcription"
            )
            st.plotly_chart(fig6, use_container_width=True)
        
        # metrics
        st.subheader("physiological metrics")
        metrics = model.calculate_metrics(solution)
        
        col5, col6, col7, col8 = st.columns(4)
        col5.metric("mean t3", f"{metrics['mean_t3']:.2f}")
        col6.metric("t4/t3 ratio", f"{metrics['t4_t3_ratio']:.1f}")
        col7.metric("receptor occupancy", f"{metrics['receptor_occupancy']:.2%}")
        col8.metric("metabolic rate", f"{metrics['mean_metabolic_rate']:.2f}")
    
    elif sim_mode == "hypothyroidism":
        st.header("hypothyroid state simulation")
        
        severity = st.slider("disease severity", 0.0, 1.0, 0.5, 0.1)
        
        # normal vs hypothyroid
        t_normal, sol_normal = model.simulate((0, sim_time))
        t_hypo, sol_hypo = simulate_hypothyroidism(
            ThyroidReceptorModel(params), severity
        )
        
        col1, col2 = st.columns(2)
        
        with col1:
            fig1 = create_comparison_plot(
                {'normal': (t_normal, sol_normal), 
                 'hypothyroid': (t_hypo, sol_hypo)},
                2, 't3 concentration', "t3 levels comparison"
            )
            st.plotly_chart(fig1, use_container_width=True)
            
            fig3 = create_comparison_plot(
                {'normal': (t_normal, sol_normal),
                 'hypothyroid': (t_hypo, sol_hypo)},
                11, 'metabolic rate', "metabolic rate comparison"
            )
            st.plotly_chart(fig3, use_container_width=True)
        
        with col2:
            fig2 = create_comparison_plot(
                {'normal': (t_normal, sol_normal),
                 'hypothyroid': (t_hypo, sol_hypo)},
                0, 'tsh concentration', "tsh levels comparison"
            )
            st.plotly_chart(fig2, use_container_width=True)
            
            fig4 = create_comparison_plot(
                {'normal': (t_normal, sol_normal),
                 'hypothyroid': (t_hypo, sol_hypo)},
                10, 'developmental factor', "developmental impact"
            )
            st.plotly_chart(fig4, use_container_width=True)
        
        # metrics comparison
        st.subheader("condition comparison")
        metrics_normal = model.calculate_metrics(sol_normal)
        metrics_hypo = model.calculate_metrics(sol_hypo)
        
        col5, col6, col7 = st.columns(3)
        col5.metric("t3 change", f"{metrics_hypo['mean_t3']:.2f}",
                   f"{(metrics_hypo['mean_t3']-metrics_normal['mean_t3']):.2f}")
        col6.metric("metabolic rate change", f"{metrics_hypo['mean_metabolic_rate']:.2f}",
                   f"{(metrics_hypo['mean_metabolic_rate']-metrics_normal['mean_metabolic_rate']):.2f}")
        col7.metric("tsh change", f"{metrics_hypo['mean_tsh']:.2f}",
                   f"{(metrics_hypo['mean_tsh']-metrics_normal['mean_tsh']):.2f}")
    
    elif sim_mode == "hyperthyroidism":
        st.header("hyperthyroid state simulation")
        
        severity = st.slider("disease severity", 0.0, 1.0, 0.5, 0.1)
        
        t_normal, sol_normal = model.simulate((0, sim_time))
        t_hyper, sol_hyper = simulate_hyperthyroidism(
            ThyroidReceptorModel(params), severity
        )
        
        col1, col2 = st.columns(2)
        
        with col1:
            fig1 = create_comparison_plot(
                {'normal': (t_normal, sol_normal),
                 'hyperthyroid': (t_hyper, sol_hyper)},
                2, 't3 concentration', "t3 levels comparison"
            )
            st.plotly_chart(fig1, use_container_width=True)
            
            fig3 = create_comparison_plot(
                {'normal': (t_normal, sol_normal),
                 'hyperthyroid': (t_hyper, sol_hyper)},
                11, 'metabolic rate', "metabolic rate comparison"
            )
            st.plotly_chart(fig3, use_container_width=True)
        
        with col2:
            fig2 = create_comparison_plot(
                {'normal': (t_normal, sol_normal),
                 'hyperthyroid': (t_hyper, sol_hyper)},
                0, 'tsh concentration', "tsh suppression"
            )
            st.plotly_chart(fig2, use_container_width=True)
            
            fig4 = create_comparison_plot(
                {'normal': (t_normal, sol_normal),
                 'hyperthyroid': (t_hyper, sol_hyper)},
                9, 'metabolic enzymes', "metabolic enzyme induction"
            )
            st.plotly_chart(fig4, use_container_width=True)
    
    elif sim_mode == "receptor mutation":
        st.header("thyroid receptor mutation simulation")
        
        binding_defect = st.slider("binding defect severity", 0.0, 1.0, 0.7, 0.1)
        
        t_normal, sol_normal = model.simulate((0, sim_time))
        t_mut, sol_mut = simulate_tr_mutation(
            ThyroidReceptorModel(params), binding_defect
        )
        
        col1, col2 = st.columns(2)
        
        with col1:
            fig1 = create_comparison_plot(
                {'normal': (t_normal, sol_normal),
                 'mutant': (t_mut, sol_mut)},
                5, 'tr-t3 complex', "receptor binding comparison"
            )
            st.plotly_chart(fig1, use_container_width=True)
            
            fig3 = create_comparison_plot(
                {'normal': (t_normal, sol_normal),
                 'mutant': (t_mut, sol_mut)},
                8, 'target mrna', "transcriptional response"
            )
            st.plotly_chart(fig3, use_container_width=True)
        
        with col2:
            fig2 = create_comparison_plot(
                {'normal': (t_normal, sol_normal),
                 'mutant': (t_mut, sol_mut)},
                2, 't3 concentration', "compensatory t3 elevation"
            )
            st.plotly_chart(fig2, use_container_width=True)
            
            fig4 = create_comparison_plot(
                {'normal': (t_normal, sol_normal),
                 'mutant': (t_mut, sol_mut)},
                11, 'metabolic rate', "metabolic dysfunction"
            )
            st.plotly_chart(fig4, use_container_width=True)
    
    elif sim_mode == "development":
        st.header("developmental thyroid signaling")
        
        t_dev, sol_dev = simulate_development(ThyroidReceptorModel(params))
        
        col1, col2 = st.columns(2)
        
        with col1:
            fig1 = create_time_series_plot(
                t_dev, sol_dev, [2, 1], ['t3', 't4'],
                "hormone levels during development"
            )
            st.plotly_chart(fig1, use_container_width=True)
            
            fig3 = create_time_series_plot(
                t_dev, sol_dev, [6, 8], ['activated receptor', 'target mrna'],
                "transcriptional program"
            )
            st.plotly_chart(fig3, use_container_width=True)
        
        with col2:
            fig2 = create_time_series_plot(
                t_dev, sol_dev, [10], ['developmental factor'],
                "developmental progression"
            )
            st.plotly_chart(fig2, use_container_width=True)
            
            fig4 = create_time_series_plot(
                t_dev, sol_dev, [11], ['metabolic rate'],
                "metabolic maturation"
            )
            st.plotly_chart(fig4, use_container_width=True)
        
        st.subheader("developmental stages")
        st.write("early (0-50h): low hormone, receptor priming")
        st.write("mid (50-100h): hormone surge, rapid development")
        st.write("late (100-150h): mature hormone levels, homeostasis")
    
    elif sim_mode == "parameter exploration":
        st.header("parameter sensitivity analysis")
        
        # parameter sweep
        param_name = st.selectbox(
            "parameter to explore",
            ["k_syn_t4", "k_d2", "k_bind_t3", "k_trans_act"]
        )
        
        param_range = st.slider(
            "parameter range (fold change)",
            0.1, 3.0, (0.5, 2.0), 0.1
        )
        
        n_values = 5
        param_values = np.linspace(param_range[0], param_range[1], n_values)
        
        results = {}
        for pval in param_values:
            temp_params = ThyroidReceptorParameters()
            setattr(temp_params, param_name, pval * getattr(params, param_name))
            temp_model = ThyroidReceptorModel(temp_params)
            t, sol = temp_model.simulate((0, 200))
            results[f"{pval:.2f}x"] = (t, sol)
        
        col1, col2 = st.columns(2)
        
        with col1:
            fig1 = create_comparison_plot(
                results, 2, 't3 concentration',
                f"t3 response to {param_name}"
            )
            st.plotly_chart(fig1, use_container_width=True)
            
            fig3 = create_comparison_plot(
                results, 11, 'metabolic rate',
                f"metabolic rate vs {param_name}"
            )
            st.plotly_chart(fig3, use_container_width=True)
        
        with col2:
            fig2 = create_comparison_plot(
                results, 5, 'tr-t3 complex',
                f"receptor binding vs {param_name}"
            )
            st.plotly_chart(fig2, use_container_width=True)
            
            fig4 = create_comparison_plot(
                results, 8, 'target mrna',
                f"transcription vs {param_name}"
            )
            st.plotly_chart(fig4, use_container_width=True)
        
        # sensitivity heatmap
        st.subheader("parameter sensitivity matrix")
        
        variables_of_interest = [2, 5, 8, 11]  # t3, tr-t3, mrna, metabolism
        var_labels = ['t3', 'tr-t3', 'mrna', 'metabolism']
        
        sensitivity_matrix = np.zeros((len(var_labels), n_values))
        
        for i, pval in enumerate(param_values):
            temp_params = ThyroidReceptorParameters()
            setattr(temp_params, param_name, pval * getattr(params, param_name))
            temp_model = ThyroidReceptorModel(temp_params)
            t, sol = temp_model.simulate((0, 200))
            
            for j, var_idx in enumerate(variables_of_interest):
                sensitivity_matrix[j, i] = np.mean(sol[-100:, var_idx])
        
        # normalize
        sensitivity_matrix = sensitivity_matrix / sensitivity_matrix[:, n_values//2:n_values//2+1]
        
        fig5 = create_heatmap(
            sensitivity_matrix,
            [f"{pv:.2f}x" for pv in param_values],
            var_labels,
            f"normalized response to {param_name}"
        )
        st.plotly_chart(fig5, use_container_width=True)
    
    # footer
    st.sidebar.markdown("---")
    st.sidebar.caption("thyroid hormone receptor signaling model")
    st.sidebar.caption("comprehensive mechanistic simulation")


if __name__ == "__main__":
    main()
