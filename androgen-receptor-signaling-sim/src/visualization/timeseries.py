"""
time series visualization for ar signaling dynamics
"""

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from typing import List, Optional, Union, Dict
from pathlib import Path

from ..simulation.simulator import SimulationResult


# set visualization style
sns.set_style("darkgrid")
plt.rcParams['figure.dpi'] = 300
plt.rcParams['font.size'] = 10
plt.rcParams['axes.labelsize'] = 12
plt.rcParams['axes.titlesize'] = 14
plt.rcParams['legend.fontsize'] = 9


def plot_dynamics(result: SimulationResult,
                 species: Optional[List[str]] = None,
                 figsize: tuple = (14, 10),
                 save: Optional[str] = None,
                 show: bool = True,
                 title: Optional[str] = None) -> plt.Figure:
    """
    plot temporal dynamics of selected species
    
    args:
        result: simulation result
        species: list of species to plot (if none, plot key species)
        figsize: figure size
        save: path to save figure
        show: whether to display figure
        title: custom title
        
    returns:
        matplotlib figure
    """
    if species is None:
        # default key species
        species = [
            'ar_free_cytoplasm',
            'ar_free_nucleus',
            'ar_testosterone_dna',
            'ar_dht_dna',
            'mrna_psa',
            'protein_psa'
        ]
    
    # filter to available species
    available_species = [s for s in species if s in result.species_names]
    
    n_species = len(available_species)
    n_cols = 2
    n_rows = (n_species + 1) // 2
    
    fig, axes = plt.subplots(n_rows, n_cols, figsize=figsize)
    axes = axes.flatten() if n_species > 1 else [axes]
    
    # color palette
    colors = sns.color_palette("husl", n_species)
    
    for i, species_name in enumerate(available_species):
        ax = axes[i]
        data = result.get_species(species_name)
        
        ax.plot(result.time, data, color=colors[i], linewidth=2, alpha=0.8)
        ax.set_xlabel('time (s)', fontweight='bold')
        ax.set_ylabel('molecules / cell', fontweight='bold')
        ax.set_title(species_name.replace('_', ' ').title(), fontweight='bold')
        ax.grid(True, alpha=0.3)
        
        # add steady-state line if applicable
        if len(result.time) > 100:
            ss_value = np.mean(data[-int(len(data)*0.1):])
            ax.axhline(ss_value, color='red', linestyle='--', alpha=0.5, linewidth=1,
                      label=f'steady state: {ss_value:.1f}')
            ax.legend()
    
    # remove extra subplots
    for i in range(n_species, len(axes)):
        fig.delaxes(axes[i])
    
    if title:
        fig.suptitle(title, fontsize=16, fontweight='bold', y=0.995)
    else:
        fig.suptitle('AR Signaling Dynamics', fontsize=16, fontweight='bold', y=0.995)
    
    plt.tight_layout()
    
    if save:
        plt.savefig(save, dpi=300, bbox_inches='tight')
    
    if show:
        plt.show()
    else:
        plt.close()
    
    return fig


def plot_species_comparison(results: List[SimulationResult],
                           species_name: str,
                           labels: Optional[List[str]] = None,
                           figsize: tuple = (12, 6),
                           save: Optional[str] = None,
                           show: bool = True,
                           title: Optional[str] = None) -> plt.Figure:
    """
    compare specific species across multiple simulation conditions
    
    args:
        results: list of simulation results
        species_name: species to compare
        labels: labels for each condition
        figsize: figure size
        save: path to save figure
        show: whether to display figure
        title: custom title
        
    returns:
        matplotlib figure
    """
    fig, ax = plt.subplots(figsize=figsize)
    
    if labels is None:
        labels = [f'condition {i+1}' for i in range(len(results))]
    
    colors = sns.color_palette("Set2", len(results))
    
    for i, (result, label) in enumerate(zip(results, labels)):
        data = result.get_species(species_name)
        ax.plot(result.time, data, color=colors[i], linewidth=2.5, 
               label=label, alpha=0.8)
    
    ax.set_xlabel('time (s)', fontsize=13, fontweight='bold')
    ax.set_ylabel('molecules / cell', fontsize=13, fontweight='bold')
    
    if title:
        ax.set_title(title, fontsize=15, fontweight='bold')
    else:
        ax.set_title(f'{species_name.replace("_", " ").title()} - Condition Comparison',
                    fontsize=15, fontweight='bold')
    
    ax.legend(loc='best', frameon=True, shadow=True)
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    
    if save:
        plt.savefig(save, dpi=300, bbox_inches='tight')
    
    if show:
        plt.show()
    else:
        plt.close()
    
    return fig


def plot_compartment_distribution(result: SimulationResult,
                                 figsize: tuple = (14, 5),
                                 save: Optional[str] = None,
                                 show: bool = True) -> plt.Figure:
    """
    plot ar distribution across compartments over time
    
    args:
        result: simulation result
        figsize: figure size
        save: path to save figure
        show: whether to display figure
        
    returns:
        matplotlib figure
    """
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=figsize)
    
    # calculate compartment totals
    ar_cytoplasm = (result.get_species('ar_free_cytoplasm') +
                   result.get_species('ar_testosterone_cytoplasm') +
                   result.get_species('ar_dht_cytoplasm'))
    
    ar_nucleus = (result.get_species('ar_free_nucleus') +
                 result.get_species('ar_testosterone_nucleus') +
                 result.get_species('ar_dht_nucleus') +
                 2 * result.get_species('ar_testosterone_dimer_nucleus') +
                 2 * result.get_species('ar_dht_dimer_nucleus') +
                 2 * result.get_species('ar_testosterone_dna') +
                 2 * result.get_species('ar_dht_dna'))
    
    # time series
    ax1.plot(result.time, ar_cytoplasm, label='cytoplasm', linewidth=2.5, color='#2E86AB')
    ax1.plot(result.time, ar_nucleus, label='nucleus', linewidth=2.5, color='#A23B72')
    ax1.set_xlabel('time (s)', fontsize=12, fontweight='bold')
    ax1.set_ylabel('ar molecules', fontsize=12, fontweight='bold')
    ax1.set_title('AR Compartment Distribution', fontsize=14, fontweight='bold')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    # stacked area
    ax2.fill_between(result.time, 0, ar_cytoplasm, label='cytoplasm', 
                    alpha=0.7, color='#2E86AB')
    ax2.fill_between(result.time, ar_cytoplasm, ar_cytoplasm + ar_nucleus, 
                    label='nucleus', alpha=0.7, color='#A23B72')
    ax2.set_xlabel('time (s)', fontsize=12, fontweight='bold')
    ax2.set_ylabel('ar molecules', fontsize=12, fontweight='bold')
    ax2.set_title('AR Compartment Distribution (Stacked)', fontsize=14, fontweight='bold')
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    
    plt.tight_layout()
    
    if save:
        plt.savefig(save, dpi=300, bbox_inches='tight')
    
    if show:
        plt.show()
    else:
        plt.close()
    
    return fig


def plot_ligand_comparison(result: SimulationResult,
                          figsize: tuple = (14, 5),
                          save: Optional[str] = None,
                          show: bool = True) -> plt.Figure:
    """
    compare testosterone vs dht signaling
    
    args:
        result: simulation result
        figsize: figure size
        save: path to save figure
        show: whether to display figure
        
    returns:
        matplotlib figure
    """
    fig, axes = plt.subplots(1, 3, figsize=figsize)
    
    # ligand-bound ar in nucleus
    ar_t = result.get_species('ar_testosterone_nucleus')
    ar_dht = result.get_species('ar_dht_nucleus')
    
    axes[0].plot(result.time, ar_t, label='testosterone', linewidth=2.5, color='#F18F01')
    axes[0].plot(result.time, ar_dht, label='dht', linewidth=2.5, color='#C73E1D')
    axes[0].set_xlabel('time (s)', fontweight='bold')
    axes[0].set_ylabel('molecules', fontweight='bold')
    axes[0].set_title('Nuclear AR:Ligand', fontweight='bold')
    axes[0].legend()
    axes[0].grid(True, alpha=0.3)
    
    # dna-bound ar
    ar_t_dna = result.get_species('ar_testosterone_dna')
    ar_dht_dna = result.get_species('ar_dht_dna')
    
    axes[1].plot(result.time, ar_t_dna, label='testosterone', linewidth=2.5, color='#F18F01')
    axes[1].plot(result.time, ar_dht_dna, label='dht', linewidth=2.5, color='#C73E1D')
    axes[1].set_xlabel('time (s)', fontweight='bold')
    axes[1].set_ylabel('molecules', fontweight='bold')
    axes[1].set_title('DNA-Bound AR', fontweight='bold')
    axes[1].legend()
    axes[1].grid(True, alpha=0.3)
    
    # ratio dht/testosterone
    ratio = (ar_dht + 1) / (ar_t + 1)  # add 1 to avoid division by zero
    
    axes[2].plot(result.time, ratio, linewidth=2.5, color='#6A4C93')
    axes[2].axhline(1.0, color='gray', linestyle='--', alpha=0.5)
    axes[2].set_xlabel('time (s)', fontweight='bold')
    axes[2].set_ylabel('ratio', fontweight='bold')
    axes[2].set_title('DHT/Testosterone Ratio', fontweight='bold')
    axes[2].grid(True, alpha=0.3)
    
    plt.tight_layout()
    
    if save:
        plt.savefig(save, dpi=300, bbox_inches='tight')
    
    if show:
        plt.show()
    else:
        plt.close()
    
    return fig


def plot_transcriptional_output(result: SimulationResult,
                                figsize: tuple = (12, 8),
                                save: Optional[str] = None,
                                show: bool = True) -> plt.Figure:
    """
    plot transcriptional outputs (mrna levels)
    
    args:
        result: simulation result
        figsize: figure size
        save: path to save figure
        show: whether to display figure
        
    returns:
        matplotlib figure
    """
    fig, axes = plt.subplots(2, 2, figsize=figsize)
    axes = axes.flatten()
    
    genes = ['mrna_psa', 'mrna_klk2', 'mrna_tmprss2', 'mrna_nkx31']
    gene_names = ['PSA (KLK3)', 'KLK2', 'TMPRSS2', 'NKX3.1']
    colors = ['#E63946', '#F77F00', '#06D6A0', '#118AB2']
    
    for i, (gene, name, color) in enumerate(zip(genes, gene_names, colors)):
        if gene in result.species_names:
            data = result.get_species(gene)
            axes[i].plot(result.time, data, linewidth=2.5, color=color)
            axes[i].set_xlabel('time (s)', fontweight='bold')
            axes[i].set_ylabel('mrna molecules', fontweight='bold')
            axes[i].set_title(f'{name} Expression', fontweight='bold')
            axes[i].grid(True, alpha=0.3)
            
            # annotate fold-change
            if len(data) > 10:
                initial = np.mean(data[:10])
                final = np.mean(data[-10:])
                fold_change = final / (initial + 1e-10)
                axes[i].text(0.05, 0.95, f'fold change: {fold_change:.2f}x',
                           transform=axes[i].transAxes,
                           verticalalignment='top',
                           bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    
    fig.suptitle('AR Target Gene Expression', fontsize=16, fontweight='bold')
    plt.tight_layout()
    
    if save:
        plt.savefig(save, dpi=300, bbox_inches='tight')
    
    if show:
        plt.show()
    else:
        plt.close()
    
    return fig


def plot_dose_response(results_dict: Dict[float, SimulationResult],
                      species_name: str,
                      figsize: tuple = (10, 6),
                      save: Optional[str] = None,
                      show: bool = True,
                      log_scale: bool = True) -> plt.Figure:
    """
    plot dose-response curve
    
    args:
        results_dict: dictionary mapping ligand concentration to simulation result
        species_name: output species to measure
        figsize: figure size
        save: path to save figure
        show: whether to display figure
        log_scale: use log scale for x-axis
        
    returns:
        matplotlib figure
    """
    fig, ax = plt.subplots(figsize=figsize)
    
    # extract steady-state values
    doses = []
    responses = []
    
    for dose, result in sorted(results_dict.items()):
        doses.append(dose)
        # take mean of last 10% as steady-state
        data = result.get_species(species_name)
        ss_value = np.mean(data[-int(len(data)*0.1):])
        responses.append(ss_value)
    
    doses = np.array(doses)
    responses = np.array(responses)
    
    # plot
    if log_scale:
        ax.semilogx(doses, responses, 'o-', markersize=8, linewidth=2.5, 
                   color='#2E86AB', markerfacecolor='#A23B72', markeredgewidth=2)
    else:
        ax.plot(doses, responses, 'o-', markersize=8, linewidth=2.5,
               color='#2E86AB', markerfacecolor='#A23B72', markeredgewidth=2)
    
    # fit hill equation if enough points
    if len(doses) >= 4:
        from scipy.optimize import curve_fit
        
        def hill(x, ymax, ec50, n):
            return ymax * x**n / (ec50**n + x**n)
        
        try:
            popt, _ = curve_fit(hill, doses, responses, 
                              p0=[np.max(responses), np.median(doses), 1.0],
                              maxfev=10000)
            
            dose_fine = np.logspace(np.log10(doses.min()), np.log10(doses.max()), 100)
            response_fit = hill(dose_fine, *popt)
            
            if log_scale:
                ax.semilogx(dose_fine, response_fit, '--', linewidth=2, 
                           color='red', alpha=0.7, label='hill fit')
            else:
                ax.plot(dose_fine, response_fit, '--', linewidth=2,
                       color='red', alpha=0.7, label='hill fit')
            
            # annotate parameters
            ax.text(0.05, 0.95, 
                   f'EC₅₀ = {popt[1]:.2e} M\nhill coef = {popt[2]:.2f}',
                   transform=ax.transAxes, verticalalignment='top',
                   bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8),
                   fontsize=11)
            
            ax.legend()
        except:
            pass
    
    ax.set_xlabel('ligand concentration (M)', fontsize=13, fontweight='bold')
    ax.set_ylabel(f'{species_name.replace("_", " ")} (molecules)', 
                 fontsize=13, fontweight='bold')
    ax.set_title('Dose-Response Curve', fontsize=15, fontweight='bold')
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    
    if save:
        plt.savefig(save, dpi=300, bbox_inches='tight')
    
    if show:
        plt.show()
    else:
        plt.close()
    
    return fig


def plot_stochastic_ensemble(results: List[SimulationResult],
                            species_name: str,
                            figsize: tuple = (12, 6),
                            save: Optional[str] = None,
                            show: bool = True,
                            alpha: float = 0.1) -> plt.Figure:
    """
    plot ensemble of stochastic trajectories
    
    args:
        results: list of stochastic simulation results
        species_name: species to plot
        figsize: figure size
        save: path to save figure
        show: whether to display figure
        alpha: transparency for individual trajectories
        
    returns:
        matplotlib figure
    """
    fig, ax = plt.subplots(figsize=figsize)
    
    # plot individual trajectories
    all_data = []
    for result in results:
        data = result.get_species(species_name)
        ax.plot(result.time, data, color='gray', alpha=alpha, linewidth=0.5)
        all_data.append(data)
    
    # plot mean and std
    all_data = np.array(all_data)
    mean_data = np.mean(all_data, axis=0)
    std_data = np.std(all_data, axis=0)
    time = results[0].time
    
    ax.plot(time, mean_data, color='#2E86AB', linewidth=3, label='mean', zorder=10)
    ax.fill_between(time, mean_data - std_data, mean_data + std_data,
                   color='#2E86AB', alpha=0.3, label='± 1 std')
    
    ax.set_xlabel('time (s)', fontsize=13, fontweight='bold')
    ax.set_ylabel(f'{species_name.replace("_", " ")} (molecules)', 
                 fontsize=13, fontweight='bold')
    ax.set_title(f'Stochastic Ensemble (n={len(results)})', 
                fontsize=15, fontweight='bold')
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
