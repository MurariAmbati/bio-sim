"""
data handling utilities
"""

import pandas as pd
import numpy as np
import json
from pathlib import Path
from typing import Optional, Dict

from ..simulation.simulator import SimulationResult
from ..models.parameters import ParameterSet


def save_simulation_result(result: SimulationResult,
                          filepath: str,
                          format: str = 'csv'):
    """
    save simulation result to file
    
    args:
        result: simulation result
        filepath: output file path
        format: 'csv', 'json', or 'hdf5'
    """
    filepath = Path(filepath)
    filepath.parent.mkdir(parents=True, exist_ok=True)
    
    if format == 'csv':
        df = result.to_dataframe()
        df.to_csv(filepath, index=False)
    
    elif format == 'json':
        data = {
            'time': result.time.tolist(),
            'species_names': result.species_names,
            'states': result.states.tolist()
        }
        with open(filepath, 'w') as f:
            json.dump(data, f, indent=2)
    
    elif format == 'hdf5':
        df = result.to_dataframe()
        df.to_hdf(filepath, key='simulation', mode='w')
    
    else:
        raise ValueError(f"unknown format: {format}")


def load_simulation_result(filepath: str,
                          format: Optional[str] = None) -> SimulationResult:
    """
    load simulation result from file
    
    args:
        filepath: input file path
        format: file format (auto-detected if none)
        
    returns:
        simulation result
    """
    filepath = Path(filepath)
    
    if format is None:
        format = filepath.suffix[1:]  # remove '.'
    
    if format == 'csv':
        df = pd.read_csv(filepath)
        time = df['time'].values
        species_names = [col for col in df.columns if col != 'time']
        states = df[species_names].values
        
    elif format == 'json':
        with open(filepath, 'r') as f:
            data = json.load(f)
        time = np.array(data['time'])
        species_names = data['species_names']
        states = np.array(data['states'])
    
    elif format in ['hdf5', 'h5']:
        df = pd.read_hdf(filepath, key='simulation')
        time = df['time'].values
        species_names = [col for col in df.columns if col != 'time']
        states = df[species_names].values
    
    else:
        raise ValueError(f"unknown format: {format}")
    
    # create parameter set (placeholder)
    params = ParameterSet()
    
    return SimulationResult(time, states, species_names, params)


def export_for_plotting(result: SimulationResult,
                       species: list,
                       filepath: str):
    """
    export selected species for external plotting
    
    args:
        result: simulation result
        species: list of species to export
        filepath: output csv file
    """
    df = pd.DataFrame({'time': result.time})
    
    for sp in species:
        if sp in result.species_names:
            df[sp] = result.get_species(sp)
    
    df.to_csv(filepath, index=False)
    print(f"exported {len(species)} species to {filepath}")


def batch_save_results(results_dict: Dict[str, SimulationResult],
                      output_dir: str,
                      format: str = 'csv'):
    """
    save multiple simulation results
    
    args:
        results_dict: dictionary mapping names to results
        output_dir: output directory
        format: file format
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    for name, result in results_dict.items():
        safe_name = name.replace(' ', '_').replace('/', '_')
        filepath = output_dir / f"{safe_name}.{format}"
        save_simulation_result(result, str(filepath), format=format)
    
    print(f"saved {len(results_dict)} results to {output_dir}")
