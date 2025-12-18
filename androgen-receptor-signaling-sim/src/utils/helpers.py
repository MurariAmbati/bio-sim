"""
helper functions
"""

import numpy as np
from typing import List, Tuple


def molecules_to_concentration(n_molecules: float,
                               volume_liters: float) -> float:
    """
    convert molecule count to molar concentration
    
    args:
        n_molecules: number of molecules
        volume_liters: volume in liters
        
    returns:
        concentration in M
    """
    avogadro = 6.022e23
    return n_molecules / (avogadro * volume_liters)


def concentration_to_molecules(concentration_M: float,
                               volume_liters: float) -> float:
    """
    convert molar concentration to molecule count
    
    args:
        concentration_M: concentration in M
        volume_liters: volume in liters
        
    returns:
        number of molecules
    """
    avogadro = 6.022e23
    return concentration_M * avogadro * volume_liters


def calculate_kd(k_on: float, k_off: float) -> float:
    """
    calculate dissociation constant from rate constants
    
    args:
        k_on: association rate constant
        k_off: dissociation rate constant
        
    returns:
        kd in same units as k_off/k_on
    """
    return k_off / k_on


def calculate_half_life(k_degradation: float) -> float:
    """
    calculate half-life from first-order degradation rate
    
    args:
        k_degradation: degradation rate constant (s^-1)
        
    returns:
        half-life in seconds
    """
    return np.log(2) / k_degradation


def calculate_steady_state_ratio(k_forward: float, 
                                 k_reverse: float) -> float:
    """
    calculate steady-state ratio for reversible reaction
    
    args:
        k_forward: forward rate constant
        k_reverse: reverse rate constant
        
    returns:
        equilibrium constant
    """
    return k_forward / k_reverse


def calculate_ec50_from_ic50(ic50: float, 
                             ligand_conc: float,
                             ligand_kd: float) -> float:
    """
    convert ic50 to ec50 using cheng-prusoff equation
    
    args:
        ic50: ic50 value
        ligand_conc: concentration of competing ligand
        ligand_kd: kd of competing ligand
        
    returns:
        ec50 value
    """
    return ic50 / (1 + ligand_conc / ligand_kd)


def time_to_units(time_seconds: float, 
                 unit: str = 'minutes') -> float:
    """
    convert time in seconds to other units
    
    args:
        time_seconds: time in seconds
        unit: 'minutes', 'hours', 'days'
        
    returns:
        time in requested units
    """
    conversions = {
        'seconds': 1,
        'minutes': 60,
        'hours': 3600,
        'days': 86400
    }
    
    if unit not in conversions:
        raise ValueError(f"unknown unit: {unit}")
    
    return time_seconds / conversions[unit]


def smooth_trajectory(data: np.ndarray, 
                     window_size: int = 5) -> np.ndarray:
    """
    smooth noisy trajectory data using moving average
    
    args:
        data: input data array
        window_size: size of smoothing window
        
    returns:
        smoothed data
    """
    if window_size < 2:
        return data
    
    kernel = np.ones(window_size) / window_size
    smoothed = np.convolve(data, kernel, mode='same')
    
    # fix edges
    for i in range(window_size // 2):
        smoothed[i] = np.mean(data[:2*i+1])
        smoothed[-(i+1)] = np.mean(data[-(2*i+1):])
    
    return smoothed


def find_peaks(data: np.ndarray, 
              threshold: float = 0.5) -> List[int]:
    """
    find peaks in time series data
    
    args:
        data: input data
        threshold: relative threshold (0-1)
        
    returns:
        list of peak indices
    """
    threshold_value = threshold * (np.max(data) - np.min(data)) + np.min(data)
    
    peaks = []
    for i in range(1, len(data) - 1):
        if (data[i] > data[i-1] and 
            data[i] > data[i+1] and 
            data[i] > threshold_value):
            peaks.append(i)
    
    return peaks


def calculate_auc(time: np.ndarray, 
                 values: np.ndarray) -> float:
    """
    calculate area under curve using trapezoidal rule
    
    args:
        time: time points
        values: values at time points
        
    returns:
        area under curve
    """
    return np.trapz(values, time)


def normalize_to_control(data: np.ndarray, 
                        control: np.ndarray) -> np.ndarray:
    """
    normalize data to control condition
    
    args:
        data: experimental data
        control: control data
        
    returns:
        normalized data (fold change)
    """
    return data / (control + 1e-10)  # avoid division by zero
