import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from typing import Dict, List, Tuple
import pandas as pd


def analyze_steady_state(solution: np.ndarray, 
                         equilibration_fraction: float = 0.5) -> Dict[str, float]:
    """analyze steady state values from simulation
    
    args:
        solution: simulation output array
        equilibration_fraction: fraction of simulation to consider equilibrated
    
    returns:
        dictionary of steady state values
    """
    
    n_points = solution.shape[0]
    start_idx = int(n_points * equilibration_fraction)
    
    steady_state = {}
    labels = [
        'tsh', 't4', 't3', 'rt3', 'tr_free', 'tr_t3', 
        'tr_t3_coreg', 'tr_corep', 'target_mrna', 
        'metabolic_enzyme', 'developmental_factor', 
        'metabolic_rate', 'd2', 'd3'
    ]
    
    for i, label in enumerate(labels):
        steady_state[label] = np.mean(solution[start_idx:, i])
        steady_state[f'{label}_std'] = np.std(solution[start_idx:, i])
    
    return steady_state


def calculate_response_time(t: np.ndarray, 
                           solution: np.ndarray, 
                           variable_idx: int,
                           threshold: float = 0.63) -> float:
    """calculate time to reach threshold of final value
    
    args:
        t: time array
        solution: simulation output
        variable_idx: index of variable to analyze
        threshold: fraction of final value (default: 1-1/e)
    
    returns:
        response time in hours
    """
    
    initial = solution[0, variable_idx]
    final = solution[-1, variable_idx]
    target = initial + threshold * (final - initial)
    
    # find crossing
    if final > initial:
        crossing_idx = np.where(solution[:, variable_idx] >= target)[0]
    else:
        crossing_idx = np.where(solution[:, variable_idx] <= target)[0]
    
    if len(crossing_idx) > 0:
        return t[crossing_idx[0]]
    else:
        return np.nan


def calculate_oscillation_metrics(t: np.ndarray,
                                  solution: np.ndarray,
                                  variable_idx: int) -> Dict[str, float]:
    """calculate oscillation properties
    
    returns:
        amplitude, period, and damping metrics
    """
    
    signal = solution[:, variable_idx]
    mean_val = np.mean(signal)
    
    # find peaks
    from scipy.signal import find_peaks
    
    peaks, _ = find_peaks(signal - mean_val)
    troughs, _ = find_peaks(-(signal - mean_val))
    
    metrics = {}
    
    if len(peaks) > 1:
        # period
        periods = np.diff(t[peaks])
        metrics['mean_period'] = np.mean(periods)
        metrics['period_std'] = np.std(periods)
        
        # amplitude
        peak_amplitudes = signal[peaks] - mean_val
        metrics['mean_amplitude'] = np.mean(np.abs(peak_amplitudes))
        
        # damping (ratio of successive peaks)
        if len(peak_amplitudes) > 1:
            damping_ratios = peak_amplitudes[1:] / peak_amplitudes[:-1]
            metrics['damping_ratio'] = np.mean(damping_ratios)
    else:
        metrics['mean_period'] = np.nan
        metrics['mean_amplitude'] = np.nan
        metrics['damping_ratio'] = np.nan
    
    return metrics


def calculate_feedback_strength(solution: np.ndarray,
                                tsh_idx: int = 0,
                                hormone_idx: int = 2) -> float:
    """calculate feedback regulation strength
    
    args:
        solution: simulation output
        tsh_idx: index of tsh variable
        hormone_idx: index of hormone variable
    
    returns:
        correlation coefficient (negative indicates negative feedback)
    """
    
    tsh = solution[:, tsh_idx]
    hormone = solution[:, hormone_idx]
    
    correlation = np.corrcoef(tsh, hormone)[0, 1]
    
    return correlation


def perform_sensitivity_analysis(base_params: dict,
                                 param_name: str,
                                 param_range: Tuple[float, float],
                                 n_points: int = 20) -> Dict[str, np.ndarray]:
    """perform sensitivity analysis for a parameter
    
    args:
        base_params: baseline parameter dictionary
        param_name: parameter to vary
        param_range: (min, max) range for parameter
        n_points: number of points to sample
    
    returns:
        dictionary with parameter values and output metrics
    """
    
    from thyroid_receptor import ThyroidReceptorModel, ThyroidReceptorParameters
    
    param_values = np.linspace(param_range[0], param_range[1], n_points)
    
    results = {
        'param_values': param_values,
        'mean_t3': np.zeros(n_points),
        'mean_t4': np.zeros(n_points),
        'metabolic_rate': np.zeros(n_points),
        'receptor_occupancy': np.zeros(n_points),
        'tsh': np.zeros(n_points),
    }
    
    for i, pval in enumerate(param_values):
        params = ThyroidReceptorParameters()
        setattr(params, param_name, pval)
        
        model = ThyroidReceptorModel(params)
        t, solution = model.simulate((0, 200))
        
        # extract metrics
        steady_state = analyze_steady_state(solution)
        results['mean_t3'][i] = steady_state['t3']
        results['mean_t4'][i] = steady_state['t4']
        results['metabolic_rate'][i] = steady_state['metabolic_rate']
        results['tsh'][i] = steady_state['tsh']
        
        # receptor occupancy
        tr_bound = steady_state['tr_t3']
        tr_total = steady_state['tr_free'] + steady_state['tr_t3']
        results['receptor_occupancy'][i] = tr_bound / tr_total if tr_total > 0 else 0
    
    return results


def create_dataframe_from_solution(t: np.ndarray, 
                                   solution: np.ndarray) -> pd.DataFrame:
    """convert simulation output to pandas dataframe
    
    args:
        t: time array
        solution: simulation output array
    
    returns:
        dataframe with time and all variables
    """
    
    labels = [
        'tsh', 't4', 't3', 'rt3', 'tr_free', 'tr_t3',
        'tr_t3_coreg', 'tr_corep', 'target_mrna',
        'metabolic_enzyme', 'developmental_factor',
        'metabolic_rate', 'd2', 'd3'
    ]
    
    data = {'time': t}
    for i, label in enumerate(labels):
        data[label] = solution[:, i]
    
    return pd.DataFrame(data)


def calculate_dose_response(hormone_doses: np.ndarray,
                           param_name: str = 'k_syn_t3') -> Dict[str, np.ndarray]:
    """calculate dose-response relationship
    
    args:
        hormone_doses: array of relative hormone doses
        param_name: parameter to modify for dosing
    
    returns:
        dictionary with doses and response metrics
    """
    
    from thyroid_receptor import ThyroidReceptorModel, ThyroidReceptorParameters
    
    n_doses = len(hormone_doses)
    
    results = {
        'doses': hormone_doses,
        't3_levels': np.zeros(n_doses),
        'receptor_activation': np.zeros(n_doses),
        'transcription': np.zeros(n_doses),
        'metabolic_response': np.zeros(n_doses),
    }
    
    for i, dose in enumerate(hormone_doses):
        params = ThyroidReceptorParameters()
        setattr(params, param_name, dose * getattr(params, param_name))
        
        model = ThyroidReceptorModel(params)
        t, solution = model.simulate((0, 200))
        
        steady_state = analyze_steady_state(solution)
        results['t3_levels'][i] = steady_state['t3']
        results['receptor_activation'][i] = steady_state['tr_t3_coreg']
        results['transcription'][i] = steady_state['target_mrna']
        results['metabolic_response'][i] = steady_state['metabolic_rate']
    
    return results


def calculate_ec50(doses: np.ndarray, response: np.ndarray) -> float:
    """calculate ec50 from dose-response data
    
    args:
        doses: dose values
        response: response values
    
    returns:
        ec50 value
    """
    
    from scipy.optimize import curve_fit
    
    # hill equation
    def hill_equation(x, ec50, hill, top, bottom):
        return bottom + (top - bottom) / (1 + (ec50 / x) ** hill)
    
    try:
        # fit
        popt, _ = curve_fit(
            hill_equation, doses, response,
            p0=[np.median(doses), 1.0, np.max(response), np.min(response)],
            maxfev=10000
        )
        return popt[0]
    except:
        return np.nan


def compare_conditions(solutions: Dict[str, np.ndarray],
                      time_points: np.ndarray) -> pd.DataFrame:
    """compare multiple simulation conditions
    
    args:
        solutions: dictionary of condition names to solution arrays
        time_points: time array
    
    returns:
        dataframe with comparative metrics
    """
    
    comparison_data = []
    
    for condition, solution in solutions.items():
        steady_state = analyze_steady_state(solution)
        
        comparison_data.append({
            'condition': condition,
            't3': steady_state['t3'],
            't4': steady_state['t4'],
            't4_t3_ratio': steady_state['t4'] / steady_state['t3'] if steady_state['t3'] > 0 else np.nan,
            'tsh': steady_state['tsh'],
            'receptor_bound': steady_state['tr_t3'],
            'receptor_free': steady_state['tr_free'],
            'transcription': steady_state['target_mrna'],
            'metabolic_rate': steady_state['metabolic_rate'],
            'developmental_factor': steady_state['developmental_factor'],
        })
    
    return pd.DataFrame(comparison_data)


def calculate_energy_landscape(solution: np.ndarray,
                               var1_idx: int,
                               var2_idx: int,
                               n_bins: int = 50) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """calculate energy landscape for phase space
    
    args:
        solution: simulation output
        var1_idx: first variable index
        var2_idx: second variable index
        n_bins: number of bins for histogram
    
    returns:
        x edges, y edges, energy landscape
    """
    
    var1 = solution[:, var1_idx]
    var2 = solution[:, var2_idx]
    
    # 2d histogram
    hist, x_edges, y_edges = np.histogram2d(
        var1, var2, bins=n_bins, density=True
    )
    
    # convert to energy: E = -log(P)
    energy = np.zeros_like(hist)
    nonzero = hist > 0
    energy[nonzero] = -np.log(hist[nonzero])
    energy[~nonzero] = np.max(energy[nonzero]) + 1
    
    return x_edges, y_edges, energy.T


def export_data(t: np.ndarray, 
               solution: np.ndarray,
               filename: str = "simulation_data.csv"):
    """export simulation data to csv
    
    args:
        t: time array
        solution: simulation output
        filename: output filename
    """
    
    df = create_dataframe_from_solution(t, solution)
    df.to_csv(filename, index=False)
    print(f"data exported to {filename}")


def calculate_system_stability(solution: np.ndarray,
                               window_size: int = 100) -> Dict[str, float]:
    """calculate system stability metrics
    
    args:
        solution: simulation output
        window_size: window for calculating variance
    
    returns:
        stability metrics
    """
    
    n_vars = solution.shape[1]
    stability = {}
    
    # calculate coefficient of variation for each variable
    for i in range(n_vars):
        var_data = solution[-window_size:, i]
        mean_val = np.mean(var_data)
        std_val = np.std(var_data)
        
        if mean_val > 0:
            cv = std_val / mean_val
        else:
            cv = np.nan
        
        stability[f'cv_var{i}'] = cv
    
    # overall stability score (lower is more stable)
    valid_cvs = [v for v in stability.values() if not np.isnan(v)]
    stability['overall_stability'] = np.mean(valid_cvs) if valid_cvs else np.nan
    
    return stability
