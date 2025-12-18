#!/usr/bin/env python
"""
basic ar signaling simulation script
"""

import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))

from src.models.ar_pathway import ArPathwayModel
from src.models.parameters import ParameterSet
from src.simulation.simulator import Simulator
from src.visualization.timeseries import plot_dynamics, plot_compartment_distribution
from src.visualization.pathway_diagram import plot_pathway_schematic
from src.visualization.network_graph import plot_reaction_network

import matplotlib.pyplot as plt


def main():
    """run basic ar signaling simulation"""
    
    print("=" * 60)
    print("AR SIGNALING SIMULATION")
    print("=" * 60)
    
    # create model
    print("\n1. Initializing AR pathway model...")
    params = ParameterSet()
    params.testosterone_external = 1e-9  # 1 nM
    params.dht_external = 1e-9  # 1 nM
    
    model = ArPathwayModel(parameters=params)
    print(f"   - Total AR: {params.ar_total} molecules")
    print(f"   - Testosterone: {params.testosterone_external*1e9:.1f} nM")
    print(f"   - DHT: {params.dht_external*1e9:.1f} nM")
    
    # create simulator
    print("\n2. Setting up simulator...")
    simulator = Simulator(model, time_end=10000, dt=10)
    print(f"   - Simulation time: {simulator.time_end} seconds")
    print(f"   - Time step: {simulator.dt} seconds")
    
    # run simulation
    print("\n3. Running simulation...")
    result = simulator.run()
    print(f"   - Simulation complete!")
    print(f"   - Generated {len(result.time)} time points")
    
    # check steady state
    ss = result.get_steady_state()
    if ss:
        print("\n4. Steady state reached:")
        print(f"   - Nuclear AR: {ss.ar_free_nucleus:.0f} molecules")
        print(f"   - DNA-bound AR: {ss.ar_testosterone_dna + ss.ar_dht_dna:.0f} molecules")
        print(f"   - PSA mRNA: {ss.mrna_psa:.0f} molecules")
        print(f"   - PSA protein: {ss.protein_psa:.0f} molecules")
    else:
        print("\n4. Steady state not reached")
    
    # create output directory
    os.makedirs('plots', exist_ok=True)
    
    # generate visualizations
    print("\n5. Generating visualizations...")
    
    # time series
    print("   - Creating time series plots...")
    plot_dynamics(result, 
                 save='plots/basic_dynamics.png',
                 show=False,
                 title='AR Signaling Dynamics - Basic Simulation')
    
    # compartment distribution
    print("   - Creating compartment distribution...")
    plot_compartment_distribution(result,
                                 save='plots/compartment_distribution.png',
                                 show=False)
    
    # pathway schematic
    print("   - Creating pathway schematic...")
    plot_pathway_schematic(model, state=ss, show_values=True,
                          save='plots/pathway_schematic.png',
                          show=False)
    
    # reaction network
    print("   - Creating reaction network...")
    plot_reaction_network(model,
                         layout='hierarchical',
                         save='plots/reaction_network.png',
                         show=False)
    
    print("\n" + "=" * 60)
    print("SIMULATION COMPLETE!")
    print("Plots saved to: plots/")
    print("=" * 60)


if __name__ == '__main__':
    main()
