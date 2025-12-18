"""
Example usage scripts for PPAR signaling model
"""

import numpy as np
import matplotlib.pyplot as plt
from ppar_model import PPARSignalingModel, PPARParameters, compare_isoform_specificity


def example_basic_simulation():
    """Basic simulation with PPAR-γ agonist"""
    print("Example 1: Basic PPAR-γ activation (TZD treatment)")
    print("-" * 60)
    
    model = PPARSignalingModel()
    
    # Simulate TZD treatment (PPAR-γ agonist)
    df = model.simulate(
        t_span=(0, 100),
        n_points=1000,
        ligand_alpha=0.0,
        ligand_gamma=1.5,  # TZD dose
        ligand_delta=0.0
    )
    
    # Print steady-state metrics
    metrics = model.get_steady_state_metrics(df)
    print("\nSteady-state metrics:")
    for key, value in metrics.items():
        print(f"  {key}: {value:.4f}")
    
    print("\n")


def example_fibrate_treatment():
    """Simulate fibrate treatment for dyslipidemia"""
    print("Example 2: Fibrate treatment (PPAR-α activation)")
    print("-" * 60)
    
    # Create model with elevated free fatty acids
    params = PPARParameters(
        free_fatty_acids=2.5,  # Elevated FFA
        inflammatory_cytokines=0.8
    )
    model = PPARSignalingModel(params)
    
    # Simulate before and after fibrate treatment
    df_baseline = model.simulate((0, 100), 500, 0.0, 0.0, 0.0)
    df_treatment = model.simulate((0, 100), 500, 2.0, 0.0, 0.0)
    
    print("\nLipid accumulation:")
    print(f"  Baseline: {df_baseline['Lipid_Accumulation'].iloc[-1]:.4f}")
    print(f"  With fibrate: {df_treatment['Lipid_Accumulation'].iloc[-1]:.4f}")
    
    print("\nFatty acid oxidation:")
    print(f"  Baseline: {df_baseline['Protein_FAO'].iloc[-1]:.4f}")
    print(f"  With fibrate: {df_treatment['Protein_FAO'].iloc[-1]:.4f}")
    
    print("\n")


def example_dose_response():
    """Dose-response analysis"""
    print("Example 3: Dose-response curve")
    print("-" * 60)
    
    model = PPARSignalingModel()
    doses = np.linspace(0, 5, 10)
    
    # Test PPAR-γ agonist
    results = model.simulate_drug_response('gamma', doses, (0, 100), 300)
    
    print("\nInsulin sensitivity vs dose:")
    for dose_key, df in results.items():
        dose = float(dose_key.split('_')[1])
        is_value = df['Insulin_Sensitivity'].iloc[-50:].mean()
        print(f"  Dose {dose:.2f}: IS = {is_value:.4f}")
    
    print("\n")


def example_isoform_comparison():
    """Compare isoform specificity"""
    print("Example 4: Isoform-specific effects")
    print("-" * 60)
    
    comparison_df = compare_isoform_specificity()
    
    print("\nIsoform activity profiles:")
    print(comparison_df.to_string(index=False))
    
    print("\n")


def example_pan_agonist():
    """Simulate pan-PPAR agonist"""
    print("Example 5: Pan-PPAR agonist")
    print("-" * 60)
    
    model = PPARSignalingModel()
    
    # Pan agonist activates all isoforms
    df = model.simulate(
        t_span=(0, 100),
        n_points=1000,
        ligand_alpha=1.0,
        ligand_gamma=1.0,
        ligand_delta=1.0
    )
    
    metrics = model.get_steady_state_metrics(df)
    
    print("\nCombined effects:")
    print(f"  Fatty acid oxidation: {metrics['Fatty Acid Oxidation']:.4f}")
    print(f"  Insulin sensitivity: {metrics['Insulin Sensitivity']:.4f}")
    print(f"  Anti-inflammatory: {metrics['Anti-inflammatory']:.4f}")
    print(f"  Metabolic health score: {metrics['Metabolic Health']:.4f}")
    
    print("\n")


def example_disease_modeling():
    """Model metabolic syndrome"""
    print("Example 6: Metabolic syndrome simulation")
    print("-" * 60)
    
    # Create disease state
    disease_params = PPARParameters(
        free_fatty_acids=3.0,  # Elevated
        glucose_concentration=10.0,  # Hyperglycemia
        inflammatory_cytokines=2.0,  # Chronic inflammation
        ppar_alpha_base=0.6,  # Reduced in disease
        ppar_gamma_base=0.7
    )
    
    model_disease = PPARSignalingModel(disease_params)
    
    # Untreated
    df_untreated = model_disease.simulate((0, 100), 500, 0.0, 0.0, 0.0)
    
    # Treated with PPAR-γ agonist
    df_treated = model_disease.simulate((0, 100), 500, 0.0, 2.0, 0.0)
    
    print("\nMetabolic health score:")
    print(f"  Untreated: {df_untreated['Metabolic_Health'].iloc[-1]:.4f}")
    print(f"  PPAR-γ agonist: {df_treated['Metabolic_Health'].iloc[-1]:.4f}")
    
    print("\nInsulin sensitivity:")
    print(f"  Untreated: {df_untreated['Insulin_Sensitivity'].iloc[-1]:.4f}")
    print(f"  PPAR-γ agonist: {df_treated['Insulin_Sensitivity'].iloc[-1]:.4f}")
    
    print("\nInflammatory state:")
    print(f"  Untreated: {df_untreated['Inflammatory_State'].iloc[-1]:.4f}")
    print(f"  PPAR-γ agonist: {df_treated['Inflammatory_State'].iloc[-1]:.4f}")
    
    print("\n")


def example_time_course_plot():
    """Generate time course plot"""
    print("Example 7: Time course visualization")
    print("-" * 60)
    
    model = PPARSignalingModel()
    
    # Simulate PPAR-α activation
    df = model.simulate((0, 100), 1000, 2.0, 0.0, 0.0)
    
    # Create plot
    fig, axes = plt.subplots(2, 2, figsize=(12, 8))
    
    # PPAR-RXR-DNA complexes
    axes[0, 0].plot(df['Time'], df['PPAR-α-RXR-DNA'], 'r-', label='PPAR-α')
    axes[0, 0].plot(df['Time'], df['PPAR-γ-RXR-DNA'], 'g-', label='PPAR-γ')
    axes[0, 0].plot(df['Time'], df['PPAR-δ-RXR-DNA'], 'b-', label='PPAR-δ')
    axes[0, 0].set_xlabel('Time (h)')
    axes[0, 0].set_ylabel('DNA-bound PPAR-RXR')
    axes[0, 0].legend()
    axes[0, 0].set_title('Nuclear Receptor Activation')
    
    # Protein expression
    axes[0, 1].plot(df['Time'], df['Protein_FAO'], 'orange', label='FAO')
    axes[0, 1].plot(df['Time'], df['Protein_IS'], 'green', label='Insulin sensitivity')
    axes[0, 1].plot(df['Time'], df['Protein_AI'], 'purple', label='Anti-inflammatory')
    axes[0, 1].set_xlabel('Time (h)')
    axes[0, 1].set_ylabel('Protein Level')
    axes[0, 1].legend()
    axes[0, 1].set_title('Target Protein Expression')
    
    # Metabolic outputs
    axes[1, 0].plot(df['Time'], df['Lipid_Accumulation'], 'r-', label='Lipid accumulation')
    axes[1, 0].plot(df['Time'], df['Insulin_Sensitivity'], 'g-', label='Insulin sensitivity')
    axes[1, 0].set_xlabel('Time (h)')
    axes[1, 0].set_ylabel('State')
    axes[1, 0].legend()
    axes[1, 0].set_title('Metabolic Outputs')
    
    # Health score
    axes[1, 1].plot(df['Time'], df['Metabolic_Health'], 'b-', linewidth=2)
    axes[1, 1].fill_between(df['Time'], df['Metabolic_Health'], alpha=0.3)
    axes[1, 1].set_xlabel('Time (h)')
    axes[1, 1].set_ylabel('Health Score')
    axes[1, 1].set_title('Overall Metabolic Health')
    
    plt.tight_layout()
    plt.savefig('ppar_time_course.png', dpi=300, bbox_inches='tight')
    print("Plot saved as 'ppar_time_course.png'")
    print("\n")


def main():
    """Run all examples"""
    print("\n" + "="*60)
    print("PPAR SIGNALING MODEL - EXAMPLE USAGE")
    print("="*60 + "\n")
    
    example_basic_simulation()
    example_fibrate_treatment()
    example_dose_response()
    example_isoform_comparison()
    example_pan_agonist()
    example_disease_modeling()
    example_time_course_plot()
    
    print("="*60)
    print("All examples completed!")
    print("="*60)


if __name__ == "__main__":
    main()
