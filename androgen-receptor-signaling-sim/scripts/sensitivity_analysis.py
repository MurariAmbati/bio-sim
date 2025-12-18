#!/usr/bin/env python
"""
parameter sensitivity analysis
"""

import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))

from src.models.ar_pathway import ArPathwayModel
from src.models.parameters import ParameterSet
from src.analysis.statistics import sensitivity_analysis, plot_sensitivity_heatmap

import pandas as pd


def main():
    """perform parameter sensitivity analysis"""
    
    print("=" * 60)
    print("AR SIGNALING - SENSITIVITY ANALYSIS")
    print("=" * 60)
    
    # create base model
    params = ParameterSet()
    model = ArPathwayModel(parameters=params)
    
    # parameters to analyze
    parameters_to_test = [
        'k_nuclear_import',
        'k_nuclear_export',
        'k_dimerization',
        'k_dna_on',
        'k_dna_off',
        'k_transcription_activated',
        'k_mrna_degradation',
        'k_ar_degradation',
        'k_5alpha_reductase',
        'k_coactivator_on',
        'ar_total',
        'num_are_sites'
    ]
    
    print(f"\nAnalyzing sensitivity of {len(parameters_to_test)} parameters...")
    print("Output species: protein_psa")
    print("Metric: steady state")
    
    # run sensitivity analysis
    print("\nRunning simulations...")
    sensitivity_df = sensitivity_analysis(
        model,
        parameters=parameters_to_test,
        variation_range=0.5,
        n_samples=50,
        output_species='protein_psa',
        metric='steady_state'
    )
    
    # display results
    print("\n" + "=" * 60)
    print("SENSITIVITY RESULTS (sorted by absolute sensitivity):")
    print("=" * 60)
    print()
    print(sensitivity_df.to_string(index=False))
    
    # save results
    os.makedirs('results', exist_ok=True)
    sensitivity_df.to_csv('results/sensitivity_analysis.csv', index=False)
    print("\nResults saved to: results/sensitivity_analysis.csv")
    
    # create visualization
    os.makedirs('plots', exist_ok=True)
    print("\nGenerating sensitivity heatmap...")
    plot_sensitivity_heatmap(sensitivity_df,
                            save='plots/sensitivity_heatmap.png',
                            show=False)
    
    print("\n" + "=" * 60)
    print("ANALYSIS COMPLETE!")
    print("=" * 60)
    
    # interpret results
    print("\nKEY FINDINGS:")
    top_3 = sensitivity_df.head(3)
    for idx, row in top_3.iterrows():
        param = row['parameter']
        coef = row['sensitivity_coefficient']
        direction = "increases" if coef > 0 else "decreases"
        print(f"  - {param}: {abs(coef):.3f} (PSA {direction} with parameter increase)")


if __name__ == '__main__':
    main()
