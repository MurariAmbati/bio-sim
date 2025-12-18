"""
example usage of thyroid hormone receptor signaling model
"""

from thyroid_receptor import (
    ThyroidReceptorModel,
    ThyroidReceptorParameters,
    simulate_hypothyroidism,
    simulate_hyperthyroidism,
    simulate_tr_mutation,
    simulate_development
)
from analysis import (
    analyze_steady_state,
    calculate_response_time,
    compare_conditions,
    create_dataframe_from_solution
)
import matplotlib.pyplot as plt
import numpy as np


def example_basic_simulation():
    """basic simulation of normal thyroid signaling"""
    print("running basic simulation...")
    
    # create model with default parameters
    model = ThyroidReceptorModel()
    
    # simulate for 200 hours
    t, solution = model.simulate((0, 200))
    
    # calculate metrics
    metrics = model.calculate_metrics(solution)
    print(f"\nphysiological metrics:")
    print(f"  mean t3: {metrics['mean_t3']:.2f} au")
    print(f"  mean t4: {metrics['mean_t4']:.2f} au")
    print(f"  t4/t3 ratio: {metrics['t4_t3_ratio']:.1f}")
    print(f"  receptor occupancy: {metrics['receptor_occupancy']:.1%}")
    print(f"  metabolic rate: {metrics['mean_metabolic_rate']:.2f} au")
    
    # plot hormone levels
    plt.figure(figsize=(12, 4))
    
    plt.subplot(1, 3, 1)
    plt.plot(t, solution[:, 1], label='t4', linewidth=2)
    plt.plot(t, solution[:, 2], label='t3', linewidth=2)
    plt.xlabel('time (hours)')
    plt.ylabel('concentration (au)')
    plt.title('thyroid hormones')
    plt.legend()
    plt.grid(alpha=0.3)
    
    plt.subplot(1, 3, 2)
    plt.plot(t, solution[:, 5], label='tr-t3', linewidth=2)
    plt.plot(t, solution[:, 6], label='tr-t3-coactivator', linewidth=2)
    plt.xlabel('time (hours)')
    plt.ylabel('concentration (au)')
    plt.title('receptor activation')
    plt.legend()
    plt.grid(alpha=0.3)
    
    plt.subplot(1, 3, 3)
    plt.plot(t, solution[:, 11], label='metabolic rate', linewidth=2)
    plt.xlabel('time (hours)')
    plt.ylabel('rate (au)')
    plt.title('metabolic effects')
    plt.legend()
    plt.grid(alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('basic_simulation.png', dpi=150)
    print("\nplot saved as basic_simulation.png")


def example_disease_comparison():
    """compare normal and disease states"""
    print("\ncomparing disease states...")
    
    # create base model
    params = ThyroidReceptorParameters()
    model = ThyroidReceptorModel(params)
    
    # run simulations
    t_normal, sol_normal = model.simulate((0, 200))
    t_hypo, sol_hypo = simulate_hypothyroidism(
        ThyroidReceptorModel(ThyroidReceptorParameters()), severity=0.6
    )
    t_hyper, sol_hyper = simulate_hyperthyroidism(
        ThyroidReceptorModel(ThyroidReceptorParameters()), severity=0.6
    )
    
    # compare conditions
    conditions = {
        'normal': sol_normal,
        'hypothyroid': sol_hypo,
        'hyperthyroid': sol_hyper
    }
    
    comparison = compare_conditions(conditions, t_normal)
    print("\ncomparison of conditions:")
    print(comparison.to_string(index=False))
    
    # plot comparison
    plt.figure(figsize=(12, 8))
    
    plt.subplot(2, 2, 1)
    plt.plot(t_normal, sol_normal[:, 2], label='normal', linewidth=2)
    plt.plot(t_hypo, sol_hypo[:, 2], label='hypothyroid', linewidth=2)
    plt.plot(t_hyper, sol_hyper[:, 2], label='hyperthyroid', linewidth=2)
    plt.xlabel('time (hours)')
    plt.ylabel('t3 concentration (au)')
    plt.title('t3 levels')
    plt.legend()
    plt.grid(alpha=0.3)
    
    plt.subplot(2, 2, 2)
    plt.plot(t_normal, sol_normal[:, 0], label='normal', linewidth=2)
    plt.plot(t_hypo, sol_hypo[:, 0], label='hypothyroid', linewidth=2)
    plt.plot(t_hyper, sol_hyper[:, 0], label='hyperthyroid', linewidth=2)
    plt.xlabel('time (hours)')
    plt.ylabel('tsh concentration (au)')
    plt.title('tsh levels')
    plt.legend()
    plt.grid(alpha=0.3)
    
    plt.subplot(2, 2, 3)
    plt.plot(t_normal, sol_normal[:, 11], label='normal', linewidth=2)
    plt.plot(t_hypo, sol_hypo[:, 11], label='hypothyroid', linewidth=2)
    plt.plot(t_hyper, sol_hyper[:, 11], label='hyperthyroid', linewidth=2)
    plt.xlabel('time (hours)')
    plt.ylabel('metabolic rate (au)')
    plt.title('metabolic rate')
    plt.legend()
    plt.grid(alpha=0.3)
    
    plt.subplot(2, 2, 4)
    plt.plot(t_normal, sol_normal[:, 10], label='normal', linewidth=2)
    plt.plot(t_hypo, sol_hypo[:, 10], label='hypothyroid', linewidth=2)
    plt.plot(t_hyper, sol_hyper[:, 10], label='hyperthyroid', linewidth=2)
    plt.xlabel('time (hours)')
    plt.ylabel('developmental factor (au)')
    plt.title('developmental effects')
    plt.legend()
    plt.grid(alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('disease_comparison.png', dpi=150)
    print("plot saved as disease_comparison.png")


def example_parameter_sensitivity():
    """explore parameter sensitivity"""
    print("\nparameter sensitivity analysis...")
    
    base_params = ThyroidReceptorParameters()
    
    # vary t4 synthesis rate
    k_syn_values = np.linspace(0.4, 1.6, 10)
    
    metrics_list = []
    for k_syn in k_syn_values:
        params = ThyroidReceptorParameters()
        params.k_syn_t4 = k_syn
        model = ThyroidReceptorModel(params)
        
        t, solution = model.simulate((0, 200))
        metrics = model.calculate_metrics(solution)
        metrics_list.append(metrics)
    
    # plot sensitivity
    plt.figure(figsize=(12, 4))
    
    plt.subplot(1, 3, 1)
    plt.plot(k_syn_values, [m['mean_t3'] for m in metrics_list], 
             'o-', linewidth=2, markersize=6)
    plt.xlabel('k_syn_t4')
    plt.ylabel('mean t3 (au)')
    plt.title('t3 sensitivity')
    plt.grid(alpha=0.3)
    
    plt.subplot(1, 3, 2)
    plt.plot(k_syn_values, [m['receptor_occupancy'] for m in metrics_list],
             'o-', linewidth=2, markersize=6)
    plt.xlabel('k_syn_t4')
    plt.ylabel('receptor occupancy')
    plt.title('receptor binding sensitivity')
    plt.grid(alpha=0.3)
    
    plt.subplot(1, 3, 3)
    plt.plot(k_syn_values, [m['mean_metabolic_rate'] for m in metrics_list],
             'o-', linewidth=2, markersize=6)
    plt.xlabel('k_syn_t4')
    plt.ylabel('metabolic rate (au)')
    plt.title('metabolic sensitivity')
    plt.grid(alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('parameter_sensitivity.png', dpi=150)
    print("plot saved as parameter_sensitivity.png")


def example_developmental_trajectory():
    """simulate developmental progression"""
    print("\ndevelopmental simulation...")
    
    model = ThyroidReceptorModel()
    t_dev, sol_dev = simulate_development(model)
    
    # analyze developmental stages
    steady_early = analyze_steady_state(sol_dev[:300, :])
    steady_mid = analyze_steady_state(sol_dev[300:600, :])
    steady_late = analyze_steady_state(sol_dev[600:, :])
    
    print("\ndevelopmental stages:")
    print(f"early (0-50h):")
    print(f"  t3: {steady_early['t3']:.2f} au")
    print(f"  developmental factor: {steady_early['developmental_factor']:.2f} au")
    print(f"mid (50-100h):")
    print(f"  t3: {steady_mid['t3']:.2f} au")
    print(f"  developmental factor: {steady_mid['developmental_factor']:.2f} au")
    print(f"late (100-150h):")
    print(f"  t3: {steady_late['t3']:.2f} au")
    print(f"  developmental factor: {steady_late['developmental_factor']:.2f} au")
    
    # plot
    plt.figure(figsize=(12, 4))
    
    plt.subplot(1, 3, 1)
    plt.plot(t_dev, sol_dev[:, 2], linewidth=2)
    plt.axvline(50, color='gray', linestyle='--', alpha=0.5)
    plt.axvline(100, color='gray', linestyle='--', alpha=0.5)
    plt.xlabel('time (hours)')
    plt.ylabel('t3 concentration (au)')
    plt.title('t3 during development')
    plt.grid(alpha=0.3)
    
    plt.subplot(1, 3, 2)
    plt.plot(t_dev, sol_dev[:, 10], linewidth=2, color='#E63946')
    plt.axvline(50, color='gray', linestyle='--', alpha=0.5)
    plt.axvline(100, color='gray', linestyle='--', alpha=0.5)
    plt.xlabel('time (hours)')
    plt.ylabel('developmental factor (au)')
    plt.title('developmental progression')
    plt.grid(alpha=0.3)
    
    plt.subplot(1, 3, 3)
    plt.plot(t_dev, sol_dev[:, 11], linewidth=2, color='#F77F00')
    plt.axvline(50, color='gray', linestyle='--', alpha=0.5)
    plt.axvline(100, color='gray', linestyle='--', alpha=0.5)
    plt.xlabel('time (hours)')
    plt.ylabel('metabolic rate (au)')
    plt.title('metabolic maturation')
    plt.grid(alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('developmental_trajectory.png', dpi=150)
    print("plot saved as developmental_trajectory.png")


if __name__ == "__main__":
    print("thyroid hormone receptor signaling examples\n" + "="*50)
    
    # run examples
    example_basic_simulation()
    example_disease_comparison()
    example_parameter_sensitivity()
    example_developmental_trajectory()
    
    print("\n" + "="*50)
    print("all examples completed")
