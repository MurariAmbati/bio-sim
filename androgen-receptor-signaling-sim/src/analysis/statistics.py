"""
statistical analysis and sensitivity analysis
"""

import numpy as np
import pandas as pd
from typing import List, Dict, Optional, Tuple
from scipy import stats
import matplotlib.pyplot as plt
import seaborn as sns

from ..models.ar_pathway import ArPathwayModel
from ..models.parameters import ParameterSet
from ..simulation.simulator import Simulator, SimulationResult


def sensitivity_analysis(model: ArPathwayModel,
                        parameters: List[str],
                        variation_range: float = 0.5,
                        n_samples: int = 100,
                        output_species: str = 'protein_psa',
                        metric: str = 'steady_state') -> pd.DataFrame:
    """
    perform local sensitivity analysis
    
    args:
        model: ar pathway model
        parameters: list of parameter names to vary
        variation_range: fraction to vary parameters (0.5 = ±50%)
        n_samples: number of samples per parameter
        output_species: species to measure as output
        metric: 'steady_state', 'peak', 'auc', 'time_to_half_max'
        
    returns:
        dataframe with sensitivity coefficients
    """
    base_params = model.params
    base_simulator = Simulator(model, time_end=10000, dt=10)
    base_result = base_simulator.run()
    base_output = _calculate_metric(base_result, output_species, metric)
    
    sensitivity_data = []
    
    for param_name in parameters:
        if not hasattr(base_params, param_name):
            print(f"Warning: parameter '{param_name}' not found")
            continue
        
        base_value = getattr(base_params, param_name)
        
        # vary parameter
        param_values = np.linspace(base_value * (1 - variation_range),
                                  base_value * (1 + variation_range),
                                  n_samples)
        
        outputs = []
        
        for param_value in param_values:
            # create modified model
            modified_params = ParameterSet()
            for attr in dir(base_params):
                if not attr.startswith('_') and hasattr(modified_params, attr):
                    setattr(modified_params, attr, getattr(base_params, attr))
            
            setattr(modified_params, param_name, param_value)
            
            modified_model = ArPathwayModel(parameters=modified_params)
            simulator = Simulator(modified_model, time_end=10000, dt=10)
            result = simulator.run()
            output = _calculate_metric(result, output_species, metric)
            outputs.append(output)
        
        # calculate sensitivity coefficient
        # normalized: (dY/Y) / (dX/X)
        relative_param_change = (param_values - base_value) / base_value
        relative_output_change = (np.array(outputs) - base_output) / (base_output + 1e-10)
        
        # linear regression for sensitivity
        slope, intercept, r_value, p_value, std_err = stats.linregress(
            relative_param_change, relative_output_change
        )
        
        sensitivity_data.append({
            'parameter': param_name,
            'sensitivity_coefficient': slope,
            'r_squared': r_value**2,
            'p_value': p_value,
            'base_value': base_value,
            'base_output': base_output
        })
    
    df = pd.DataFrame(sensitivity_data)
    df = df.sort_values('sensitivity_coefficient', key=abs, ascending=False)
    
    return df


def _calculate_metric(result: SimulationResult, 
                     species: str, 
                     metric: str) -> float:
    """
    calculate output metric from simulation result
    
    args:
        result: simulation result
        species: species name
        metric: metric type
        
    returns:
        metric value
    """
    data = result.get_species(species)
    
    if metric == 'steady_state':
        # mean of last 10%
        return np.mean(data[-int(len(data)*0.1):])
    
    elif metric == 'peak':
        return np.max(data)
    
    elif metric == 'auc':
        # area under curve
        return np.trapz(data, result.time)
    
    elif metric == 'time_to_half_max':
        max_val = np.max(data)
        half_max = max_val / 2
        idx = np.where(data >= half_max)[0]
        if len(idx) > 0:
            return result.time[idx[0]]
        return np.inf
    
    else:
        return np.mean(data)


def monte_carlo_analysis(model: ArPathwayModel,
                        parameter_distributions: Dict[str, Tuple[str, float, float]],
                        n_samples: int = 1000,
                        output_species: List[str] = ['protein_psa'],
                        confidence_level: float = 0.95) -> Dict:
    """
    perform monte carlo uncertainty analysis
    
    args:
        model: ar pathway model
        parameter_distributions: dict mapping parameter names to 
                               (distribution_type, param1, param2)
                               e.g., ('normal', mean, std) or ('uniform', min, max)
        n_samples: number of monte carlo samples
        output_species: species to track
        confidence_level: confidence level for intervals
        
    returns:
        dictionary with results and statistics
    """
    results_data = {species: [] for species in output_species}
    
    for i in range(n_samples):
        # sample parameters
        modified_params = ParameterSet()
        
        for param_name, (dist_type, p1, p2) in parameter_distributions.items():
            if dist_type == 'normal':
                value = np.random.normal(p1, p2)
            elif dist_type == 'uniform':
                value = np.random.uniform(p1, p2)
            elif dist_type == 'lognormal':
                value = np.random.lognormal(p1, p2)
            else:
                value = p1
            
            setattr(modified_params, param_name, max(0, value))  # ensure positive
        
        # run simulation
        modified_model = ArPathwayModel(parameters=modified_params)
        simulator = Simulator(modified_model, time_end=10000, dt=10)
        result = simulator.run()
        
        # extract steady-state values
        for species in output_species:
            data = result.get_species(species)
            ss_value = np.mean(data[-int(len(data)*0.1):])
            results_data[species].append(ss_value)
    
    # calculate statistics
    statistics = {}
    alpha = 1 - confidence_level
    
    for species in output_species:
        values = np.array(results_data[species])
        
        statistics[species] = {
            'mean': np.mean(values),
            'median': np.median(values),
            'std': np.std(values),
            'cv': np.std(values) / (np.mean(values) + 1e-10),
            'ci_lower': np.percentile(values, alpha/2 * 100),
            'ci_upper': np.percentile(values, (1 - alpha/2) * 100),
            'min': np.min(values),
            'max': np.max(values),
            'values': values
        }
    
    return {
        'statistics': statistics,
        'raw_data': results_data,
        'n_samples': n_samples,
        'confidence_level': confidence_level
    }


def plot_sensitivity_heatmap(sensitivity_df: pd.DataFrame,
                            figsize: tuple = (10, 8),
                            save: Optional[str] = None,
                            show: bool = True) -> plt.Figure:
    """
    visualize sensitivity analysis as heatmap
    
    args:
        sensitivity_df: sensitivity analysis dataframe
        figsize: figure size
        save: path to save
        show: whether to display
        
    returns:
        matplotlib figure
    """
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=figsize)
    
    # sensitivity coefficients
    params = sensitivity_df['parameter'].values
    sensitivities = sensitivity_df['sensitivity_coefficient'].values
    
    colors = ['red' if s < 0 else 'blue' for s in sensitivities]
    
    ax1.barh(range(len(params)), sensitivities, color=colors, alpha=0.7)
    ax1.set_yticks(range(len(params)))
    ax1.set_yticklabels([p.replace('_', ' ') for p in params], fontsize=9)
    ax1.set_xlabel('normalized sensitivity coefficient', fontsize=11, fontweight='bold')
    ax1.set_title('parameter sensitivity', fontsize=13, fontweight='bold')
    ax1.axvline(0, color='black', linestyle='--', linewidth=1)
    ax1.grid(True, alpha=0.3, axis='x')
    
    # r-squared values
    r_squared = sensitivity_df['r_squared'].values
    
    ax2.barh(range(len(params)), r_squared, color='green', alpha=0.7)
    ax2.set_yticks(range(len(params)))
    ax2.set_yticklabels([p.replace('_', ' ') for p in params], fontsize=9)
    ax2.set_xlabel('R²', fontsize=11, fontweight='bold')
    ax2.set_title('linearity (R²)', fontsize=13, fontweight='bold')
    ax2.set_xlim(0, 1)
    ax2.grid(True, alpha=0.3, axis='x')
    
    plt.tight_layout()
    
    if save:
        plt.savefig(save, dpi=300, bbox_inches='tight')
    
    if show:
        plt.show()
    else:
        plt.close()
    
    return fig


def plot_monte_carlo_results(mc_results: Dict,
                            species_name: str,
                            figsize: tuple = (12, 5),
                            save: Optional[str] = None,
                            show: bool = True) -> plt.Figure:
    """
    visualize monte carlo analysis results
    
    args:
        mc_results: results from monte_carlo_analysis
        species_name: species to plot
        figsize: figure size
        save: path to save
        show: whether to display
        
    returns:
        matplotlib figure
    """
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=figsize)
    
    stats = mc_results['statistics'][species_name]
    values = stats['values']
    
    # histogram
    ax1.hist(values, bins=50, color='#2E86AB', alpha=0.7, edgecolor='black')
    ax1.axvline(stats['mean'], color='red', linestyle='--', linewidth=2, label='mean')
    ax1.axvline(stats['median'], color='orange', linestyle='--', linewidth=2, label='median')
    ax1.axvline(stats['ci_lower'], color='green', linestyle=':', linewidth=2, label='CI')
    ax1.axvline(stats['ci_upper'], color='green', linestyle=':', linewidth=2)
    ax1.set_xlabel('value (molecules)', fontsize=11, fontweight='bold')
    ax1.set_ylabel('frequency', fontsize=11, fontweight='bold')
    ax1.set_title(f'{species_name} distribution', fontsize=13, fontweight='bold')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    # box plot with statistics
    bp = ax2.boxplot([values], vert=True, patch_artist=True, widths=0.5)
    bp['boxes'][0].set_facecolor('#A23B72')
    bp['boxes'][0].set_alpha(0.7)
    
    ax2.set_ylabel('value (molecules)', fontsize=11, fontweight='bold')
    ax2.set_title('summary statistics', fontsize=13, fontweight='bold')
    ax2.set_xticklabels([species_name.replace('_', ' ')])
    ax2.grid(True, alpha=0.3, axis='y')
    
    # add text with statistics
    stat_text = (f"mean: {stats['mean']:.1f}\n"
                f"std: {stats['std']:.1f}\n"
                f"CV: {stats['cv']:.3f}\n"
                f"CI: [{stats['ci_lower']:.1f}, {stats['ci_upper']:.1f}]")
    
    ax2.text(0.95, 0.95, stat_text, transform=ax2.transAxes,
            verticalalignment='top', horizontalalignment='right',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8),
            fontsize=9)
    
    plt.tight_layout()
    
    if save:
        plt.savefig(save, dpi=300, bbox_inches='tight')
    
    if show:
        plt.show()
    else:
        plt.close()
    
    return fig


def correlation_analysis(mc_results: Dict,
                        parameter_names: List[str],
                        output_species: str) -> pd.DataFrame:
    """
    analyze correlations between parameters and outputs
    
    args:
        mc_results: monte carlo results
        parameter_names: parameters that were varied
        output_species: output to correlate with
        
    returns:
        correlation dataframe
    """
    # would need to store parameter values during MC simulation
    # placeholder for now
    correlations = []
    
    for param in parameter_names:
        # compute correlation
        corr_coef = np.random.randn()  # placeholder
        p_value = np.random.rand()
        
        correlations.append({
            'parameter': param,
            'correlation': corr_coef,
            'p_value': p_value
        })
    
    return pd.DataFrame(correlations)
