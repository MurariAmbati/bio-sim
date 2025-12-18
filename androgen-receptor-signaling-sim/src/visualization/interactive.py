"""
interactive dashboard for ar signaling exploration
"""

import dash
from dash import dcc, html, Input, Output, State
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import numpy as np

from ..models.ar_pathway import ArPathwayModel, PathwayState
from ..models.parameters import ParameterSet, get_drug_parameters
from ..simulation.simulator import Simulator


def create_dashboard(port: int = 8050):
    """
    create and launch interactive dashboard
    
    args:
        port: port to run dashboard on
    """
    app = dash.Dash(__name__)
    
    app.layout = html.Div([
        html.H1('AR Signaling Interactive Dashboard',
               style={'textAlign': 'center', 'color': '#2E86AB'}),
        
        html.Div([
            html.Div([
                html.H3('Simulation Parameters'),
                
                html.Label('Testosterone Concentration (nM):'),
                dcc.Slider(id='testosterone-slider',
                          min=-2, max=2, step=0.1, value=0,
                          marks={i: f'10^{i}' for i in range(-2, 3)}),
                
                html.Label('DHT Concentration (nM):'),
                dcc.Slider(id='dht-slider',
                          min=-2, max=2, step=0.1, value=0,
                          marks={i: f'10^{i}' for i in range(-2, 3)}),
                
                html.Label('AR Expression Level:'),
                dcc.Slider(id='ar-expression-slider',
                          min=0.1, max=10, step=0.1, value=1,
                          marks={i: f'{i}x' for i in range(0, 11, 2)}),
                
                html.Label('Simulation Time (seconds):'),
                dcc.Input(id='time-input', type='number', 
                         value=10000, min=100, max=100000),
                
                html.Label('Drug Treatment:'),
                dcc.Dropdown(
                    id='drug-dropdown',
                    options=[
                        {'label': 'None (Control)', 'value': 'none'},
                        {'label': 'Enzalutamide', 'value': 'enzalutamide'},
                        {'label': 'Apalutamide', 'value': 'apalutamide'},
                        {'label': 'Bicalutamide', 'value': 'bicalutamide'},
                        {'label': 'Abiraterone', 'value': 'abiraterone'},
                        {'label': 'Finasteride', 'value': 'finasteride'}
                    ],
                    value='none'
                ),
                
                html.Label('Drug Concentration (μM):'),
                dcc.Slider(id='drug-conc-slider',
                          min=-1, max=2, step=0.1, value=1,
                          marks={i: f'10^{i}' for i in range(-1, 3)}),
                
                html.Button('Run Simulation', id='run-button',
                           style={'marginTop': '20px', 'fontSize': '16px'}),
                
            ], style={'width': '30%', 'display': 'inline-block', 
                     'verticalAlign': 'top', 'padding': '20px'}),
            
            html.Div([
                html.H3('Simulation Results'),
                dcc.Loading(
                    id="loading",
                    type="default",
                    children=html.Div(id='output-graphs')
                )
            ], style={'width': '68%', 'display': 'inline-block', 
                     'verticalAlign': 'top', 'padding': '20px'})
        ])
    ])
    
    @app.callback(
        Output('output-graphs', 'children'),
        Input('run-button', 'n_clicks'),
        State('testosterone-slider', 'value'),
        State('dht-slider', 'value'),
        State('ar-expression-slider', 'value'),
        State('time-input', 'value'),
        State('drug-dropdown', 'value'),
        State('drug-conc-slider', 'value')
    )
    def update_simulation(n_clicks, t_log, dht_log, ar_expr, 
                         sim_time, drug_name, drug_conc_log):
        if n_clicks is None:
            return html.Div('Click "Run Simulation" to start')
        
        # convert log values
        testosterone = 10**t_log * 1e-9  # to M
        dht = 10**dht_log * 1e-9
        drug_conc = 10**drug_conc_log * 1e-6  # to M
        
        # create model
        params = ParameterSet()
        params.testosterone_external = testosterone
        params.dht_external = dht
        params.ar_total *= ar_expr
        
        # add drug if selected
        drug = None
        if drug_name != 'none':
            drug = get_drug_parameters(drug_name, drug_conc)
        
        model = ArPathwayModel(parameters=params, drug=drug)
        
        # run simulation
        simulator = Simulator(model, time_end=sim_time, dt=sim_time/1000)
        result = simulator.run()
        
        # create plots
        fig = make_subplots(
            rows=3, cols=2,
            subplot_titles=('AR Compartment Distribution',
                          'DNA-Bound AR',
                          'Target Gene Expression (mRNA)',
                          'PSA Protein',
                          'Ligand-Bound AR',
                          'Transcriptional Activity'),
            specs=[[{"secondary_y": False}, {"secondary_y": False}],
                   [{"secondary_y": False}, {"secondary_y": False}],
                   [{"secondary_y": False}, {"secondary_y": False}]]
        )
        
        time = result.time
        
        # ar compartment distribution
        ar_cyto = (result.get_species('ar_free_cytoplasm') +
                  result.get_species('ar_testosterone_cytoplasm') +
                  result.get_species('ar_dht_cytoplasm'))
        ar_nuc = (result.get_species('ar_free_nucleus') +
                 result.get_species('ar_testosterone_nucleus') +
                 result.get_species('ar_dht_nucleus'))
        
        fig.add_trace(go.Scatter(x=time, y=ar_cyto, name='Cytoplasm',
                                line=dict(color='#2E86AB', width=2)),
                     row=1, col=1)
        fig.add_trace(go.Scatter(x=time, y=ar_nuc, name='Nucleus',
                                line=dict(color='#A23B72', width=2)),
                     row=1, col=1)
        
        # dna-bound ar
        ar_t_dna = result.get_species('ar_testosterone_dna')
        ar_dht_dna = result.get_species('ar_dht_dna')
        
        fig.add_trace(go.Scatter(x=time, y=ar_t_dna, name='Testosterone',
                                line=dict(color='#F18F01', width=2)),
                     row=1, col=2)
        fig.add_trace(go.Scatter(x=time, y=ar_dht_dna, name='DHT',
                                line=dict(color='#C73E1D', width=2)),
                     row=1, col=2)
        
        # mrna
        mrna_psa = result.get_species('mrna_psa')
        mrna_klk2 = result.get_species('mrna_klk2')
        
        fig.add_trace(go.Scatter(x=time, y=mrna_psa, name='PSA',
                                line=dict(color='#E63946', width=2)),
                     row=2, col=1)
        fig.add_trace(go.Scatter(x=time, y=mrna_klk2, name='KLK2',
                                line=dict(color='#F77F00', width=2)),
                     row=2, col=1)
        
        # protein
        protein_psa = result.get_species('protein_psa')
        
        fig.add_trace(go.Scatter(x=time, y=protein_psa, name='PSA Protein',
                                line=dict(color='#06D6A0', width=3)),
                     row=2, col=2)
        
        # ligand-bound ar
        ar_t_bound = (result.get_species('ar_testosterone_cytoplasm') +
                     result.get_species('ar_testosterone_nucleus'))
        ar_dht_bound = (result.get_species('ar_dht_cytoplasm') +
                       result.get_species('ar_dht_nucleus'))
        
        fig.add_trace(go.Scatter(x=time, y=ar_t_bound, name='AR:T',
                                line=dict(color='#F18F01', width=2)),
                     row=3, col=1)
        fig.add_trace(go.Scatter(x=time, y=ar_dht_bound, name='AR:DHT',
                                line=dict(color='#C73E1D', width=2)),
                     row=3, col=1)
        
        # transcriptional activity (dna-bound with coactivator)
        coactivator_complex = result.get_species('ar_dna_coactivator')
        
        fig.add_trace(go.Scatter(x=time, y=coactivator_complex,
                                name='Active Transcription',
                                line=dict(color='#6A4C93', width=3)),
                     row=3, col=2)
        
        # update layout
        fig.update_xaxes(title_text="Time (s)")
        fig.update_yaxes(title_text="Molecules")
        
        fig.update_layout(
            height=900,
            showlegend=True,
            title_text=f"AR Signaling Dynamics - {drug_name.title() if drug_name != 'none' else 'Control'}",
            title_font_size=18
        )
        
        # summary statistics
        steady_state = result.get_steady_state()
        if steady_state:
            summary = html.Div([
                html.H4('Steady State Values:'),
                html.P(f'Nuclear AR: {steady_state.ar_free_nucleus:.0f} molecules'),
                html.P(f'DNA-Bound AR: {(steady_state.ar_testosterone_dna + steady_state.ar_dht_dna):.0f} molecules'),
                html.P(f'PSA mRNA: {steady_state.mrna_psa:.0f} molecules'),
                html.P(f'PSA Protein: {steady_state.protein_psa:.0f} molecules'),
            ], style={'marginTop': '20px', 'padding': '15px',
                     'backgroundColor': '#f0f0f0', 'borderRadius': '5px'})
        else:
            summary = html.P('Steady state not reached', 
                           style={'marginTop': '20px', 'color': 'red'})
        
        return html.Div([
            dcc.Graph(figure=fig),
            summary
        ])
    
    return app


def launch_dashboard(port: int = 8050, debug: bool = True):
    """
    launch the interactive dashboard
    
    args:
        port: port to run on
        debug: enable debug mode
    """
    app = create_dashboard(port)
    print(f'\nLaunching AR Signaling Dashboard...')
    print(f'Open browser to: http://localhost:{port}\n')
    app.run_server(debug=debug, port=port)


if __name__ == '__main__':
    launch_dashboard()
