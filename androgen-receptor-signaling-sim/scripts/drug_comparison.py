#!/usr/bin/env python
"""
drug treatment comparison simulation
"""

import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))

from src.models.ar_pathway import ArPathwayModel
from src.models.parameters import ParameterSet, get_drug_parameters
from src.simulation.simulator import Simulator
from src.visualization.timeseries import plot_species_comparison, plot_transcriptional_output

import matplotlib.pyplot as plt
import numpy as np


def main():
    """compare different drug treatments"""
    
    print("=" * 60)
    print("AR SIGNALING - DRUG TREATMENT COMPARISON")
    print("=" * 60)
    
    # drugs to test
    drugs = ['control', 'enzalutamide', 'apalutamide', 'bicalutamide', 'abiraterone']
    drug_concentration = 10e-6  # 10 μM
    
    results = {}
    
    for drug_name in drugs:
        print(f"\nSimulating: {drug_name.upper()}")
        
        # create model
        params = ParameterSet()
        params.testosterone_external = 1e-9
        params.dht_external = 1e-9
        
        if drug_name == 'control':
            drug = None
        else:
            drug = get_drug_parameters(drug_name, drug_concentration)
            print(f"  Drug concentration: {drug_concentration*1e6:.1f} μM")
        
        model = ArPathwayModel(parameters=params, drug=drug)
        
        # simulate
        simulator = Simulator(model, time_end=10000, dt=10)
        result = simulator.run()
        results[drug_name] = result
        
        # report steady state
        ss = result.get_steady_state()
        if ss:
            print(f"  PSA mRNA: {ss.mrna_psa:.0f} molecules")
            print(f"  PSA protein: {ss.protein_psa:.0f} molecules")
            
            # calculate inhibition
            if drug_name != 'control':
                control_ss = results['control'].get_steady_state()
                inhibition = (1 - ss.protein_psa / control_ss.protein_psa) * 100
                print(f"  Inhibition: {inhibition:.1f}%")
    
    # create output directory
    os.makedirs('plots', exist_ok=True)
    
    # generate comparison plots
    print("\n" + "=" * 60)
    print("Generating comparison plots...")
    
    # compare PSA protein
    plot_species_comparison(list(results.values()),
                           'protein_psa',
                           labels=list(results.keys()),
                           title='PSA Protein - Drug Comparison',
                           save='plots/drug_comparison_psa.png',
                           show=False)
    
    # compare DNA-bound AR
    plot_species_comparison(list(results.values()),
                           'ar_dht_dna',
                           labels=list(results.keys()),
                           title='DNA-Bound AR - Drug Comparison',
                           save='plots/drug_comparison_dna.png',
                           show=False)
    
    # create dose-response curves
    print("Generating dose-response curves...")
    
    for drug_name in ['enzalutamide', 'apalutamide']:
        print(f"  {drug_name}...")
        
        concentrations = np.logspace(-10, -5, 10)  # 0.1 nM to 10 μM
        dose_results = {}
        
        for conc in concentrations:
            params = ParameterSet()
            drug = get_drug_parameters(drug_name, conc)
            model = ArPathwayModel(parameters=params, drug=drug)
            simulator = Simulator(model, time_end=10000, dt=10)
            result = simulator.run()
            dose_results[conc] = result
        
        from src.visualization.timeseries import plot_dose_response
        plot_dose_response(dose_results,
                          'protein_psa',
                          title=f'{drug_name.title()} Dose-Response',
                          save=f'plots/dose_response_{drug_name}.png',
                          show=False)
    
    print("\n" + "=" * 60)
    print("ANALYSIS COMPLETE!")
    print("Plots saved to: plots/")
    print("=" * 60)


if __name__ == '__main__':
    main()
