"""
bifurcation analysis for ar signaling
"""

import numpy as np
import matplotlib.pyplot as plt
from typing import Optional, List, Tuple
from scipy.optimize import fsolve

from ..models.ar_pathway import ArPathwayModel, PathwayState
from ..simulation.simulator import Simulator


def bifurcation_analysis(model: ArPathwayModel,
                        parameter_name: str,
                        param_range: Tuple[float, float],
                        n_points: int = 50,
                        output_species: str = 'protein_psa') -> dict:
    """
    perform bifurcation analysis varying a single parameter
    
    args:
        model: ar pathway model
        parameter_name: parameter to vary
        param_range: (min, max) range for parameter
        n_points: number of parameter values to test
        output_species: species to track as output
        
    returns:
        dictionary with bifurcation data
    """
    param_values = np.linspace(param_range[0], param_range[1], n_points)
    
    steady_states = []
    stable = []
    
    for param_val in param_values:
        # update parameter
        model.params.update(**{parameter_name: param_val})
        
        # run simulation to find steady state
        simulator = Simulator(model, time_end=20000, dt=20)
        result = simulator.run()
        
        # check for steady state
        data = result.get_species(output_species)
        ss_value = np.mean(data[-int(len(data)*0.1):])
        cv = np.std(data[-int(len(data)*0.1):]) / (ss_value + 1e-10)
        
        steady_states.append(ss_value)
        stable.append(cv < 0.01)  # stable if low variation
    
    return {
        'parameter': parameter_name,
        'parameter_values': param_values,
        'steady_states': np.array(steady_states),
        'stability': np.array(stable),
        'output_species': output_species
    }


def plot_bifurcation_diagram(bifurcation_data: dict,
                            figsize: tuple = (10, 6),
                            save: Optional[str] = None,
                            show: bool = True) -> plt.Figure:
    """
    plot bifurcation diagram
    
    args:
        bifurcation_data: output from bifurcation_analysis
        figsize: figure size
        save: path to save
        show: whether to display
        
    returns:
        matplotlib figure
    """
    fig, ax = plt.subplots(figsize=figsize)
    
    param_vals = bifurcation_data['parameter_values']
    steady_states = bifurcation_data['steady_states']
    stability = bifurcation_data['stability']
    
    # plot stable branches
    stable_mask = stability
    ax.plot(param_vals[stable_mask], steady_states[stable_mask],
           'o', color='#2E86AB', markersize=5, label='stable', alpha=0.7)
    
    # plot unstable branches
    unstable_mask = ~stability
    if np.any(unstable_mask):
        ax.plot(param_vals[unstable_mask], steady_states[unstable_mask],
               'o', color='red', markersize=5, label='unstable', alpha=0.7)
    
    ax.set_xlabel(bifurcation_data['parameter'].replace('_', ' '),
                 fontsize=12, fontweight='bold')
    ax.set_ylabel(f"{bifurcation_data['output_species']} (steady state)",
                 fontsize=12, fontweight='bold')
    ax.set_title('bifurcation diagram', fontsize=14, fontweight='bold')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    
    if save:
        plt.savefig(save, dpi=300, bbox_inches='tight')
    
    if show:
        plt.show()
    else:
        plt.close()
    
    return fig


def phase_plane_analysis(model: ArPathwayModel,
                        species1: str,
                        species2: str,
                        n_trajectories: int = 10,
                        figsize: tuple = (10, 8),
                        save: Optional[str] = None,
                        show: bool = True) -> plt.Figure:
    """
    create phase plane plot for two species
    
    args:
        model: ar pathway model
        species1: first species
        species2: second species
        n_trajectories: number of trajectories with different initial conditions
        figsize: figure size
        save: path to save
        show: whether to display
        
    returns:
        matplotlib figure
    """
    fig, ax = plt.subplots(figsize=figsize)
    
    colors = plt.cm.viridis(np.linspace(0, 1, n_trajectories))
    
    for i in range(n_trajectories):
        # vary initial conditions
        initial_state = model.get_initial_state()
        
        # perturb initial state
        state_array = initial_state.to_array()
        state_array = state_array * (1 + 0.5 * (np.random.rand(len(state_array)) - 0.5))
        initial_state = PathwayState.from_array(state_array)
        
        # simulate
        simulator = Simulator(model, time_end=10000, dt=10)
        result = simulator.run(initial_state=initial_state)
        
        # plot trajectory
        x = result.get_species(species1)
        y = result.get_species(species2)
        
        ax.plot(x, y, color=colors[i], linewidth=1.5, alpha=0.7)
        ax.plot(x[0], y[0], 'o', color=colors[i], markersize=8, 
               markeredgecolor='black', markeredgewidth=1)  # start
        ax.plot(x[-1], y[-1], 's', color=colors[i], markersize=8,
               markeredgecolor='black', markeredgewidth=1)  # end
    
    ax.set_xlabel(species1.replace('_', ' '), fontsize=12, fontweight='bold')
    ax.set_ylabel(species2.replace('_', ' '), fontsize=12, fontweight='bold')
    ax.set_title('phase plane analysis', fontsize=14, fontweight='bold')
    ax.grid(True, alpha=0.3)
    
    # add legend
    from matplotlib.lines import Line2D
    legend_elements = [
        Line2D([0], [0], marker='o', color='w', markerfacecolor='gray',
              markersize=8, label='initial state'),
        Line2D([0], [0], marker='s', color='w', markerfacecolor='gray',
              markersize=8, label='final state')
    ]
    ax.legend(handles=legend_elements)
    
    plt.tight_layout()
    
    if save:
        plt.savefig(save, dpi=300, bbox_inches='tight')
    
    if show:
        plt.show()
    else:
        plt.close()
    
    return fig
