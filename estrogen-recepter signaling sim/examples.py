"""
example usage of estrogen receptor simulation
"""

import numpy as np
import pandas as pd
from models import CellType, LigandLibrary
from simulation import (
    ERSimulation,
    DoseResponseSimulation,
    TimeCourseSimulation,
    TissueSpecificSimulation
)
from visualization import PathwayVisualizer


def example_1_basic_timecourse():
    """basic time course simulation with estradiol"""
    print("\n=== example 1: basic time course ===\n")
    
    # create simulation for mcf-7 breast cancer cells
    sim = ERSimulation(
        cell_type=CellType.BREAST_CANCER_MCF7,
        duration_hours=48.0,
        dt_minutes=1.0
    )
    
    # add estradiol treatment
    sim.add_treatment('estradiol', dose_mg=1.0)
    
    # run simulation
    print("running 48-hour simulation...")
    df = sim.run(verbose=True)
    
    # get summary statistics
    stats = sim.get_summary_statistics(df)
    
    print(f"\nresults:")
    print(f"  peak er-α concentration: {stats['er_alpha_max']:.2f} nm")
    print(f"  transcriptional auc: {stats['transcription_auc']:.1f}")
    print(f"  cyclin d1 peak: {stats['cyclin_d1_peak']:.2f}")
    print(f"  mapk peak activation: {stats['mapk_peak']:.3f}")
    print(f"  final cell divisions: {int(stats['final_divisions'])}")
    
    return df


def example_2_serm_comparison():
    """compare tamoxifen and raloxifene"""
    print("\n=== example 2: serm comparison ===\n")
    
    tc_sim = TimeCourseSimulation(CellType.BREAST_CANCER_MCF7)
    
    # compare serms at same dose
    ligands = ['estradiol', 'tamoxifen', 'raloxifene']
    results = tc_sim.compare_ligands(ligands, dose=10.0, duration_hours=24.0)
    
    print("\npeak cyclin d1 expression:")
    for ligand, df in results.items():
        peak = df['gene_cyclin_d1_protein'].max()
        print(f"  {ligand:12s}: {peak:.3f}")
    
    return results


def example_3_dose_response():
    """dose-response curve for tamoxifen"""
    print("\n=== example 3: dose-response analysis ===\n")
    
    dr_sim = DoseResponseSimulation(CellType.BREAST_CANCER_MCF7)
    
    # test dose range (logarithmic spacing)
    doses = np.logspace(-2, 2, 10)  # 0.01 to 100 mg
    
    print(f"testing {len(doses)} doses from {doses[0]:.3f} to {doses[-1]:.1f} mg...")
    
    dr_df = dr_sim.run_dose_response(
        ligand_name='tamoxifen',
        doses=doses,
        duration_hours=24.0
    )
    
    # calculate ec50 for proliferation inhibition
    ec50 = dr_sim.calculate_ec50(dr_df, 'proliferation_rate_mean')
    print(f"\nec50 (proliferation): {ec50:.3f} mg")
    
    # display dose-response table
    print("\ndose-response data:")
    print(dr_df[['dose', 'proliferation_rate_mean', 'survival_signal_mean']].to_string(index=False))
    
    return dr_df


def example_4_tissue_specificity():
    """demonstrate tissue-specific effects of raloxifene"""
    print("\n=== example 4: tissue specificity ===\n")
    
    tissue_sim = TissueSpecificSimulation()
    
    # compare raloxifene effects across tissues
    tissues = [
        CellType.BREAST_CANCER_MCF7,
        CellType.OSTEOBLAST,
        CellType.ENDOMETRIAL
    ]
    
    print("simulating raloxifene (20 mg) across tissues...")
    results = tissue_sim.run_tissue_comparison(
        ligand_name='raloxifene',
        dose=20.0,
        tissues=tissues
    )
    
    # calculate selectivity indices
    print("\ntissue selectivity (relative to breast):")
    
    metrics = [
        ('er_alpha_transcriptional_activity', 'transcriptional activity'),
        ('gene_cyclin_d1_protein', 'cyclin d1 (proliferation)'),
        ('nongenomic_mapk_active', 'mapk activation')
    ]
    
    for metric, label in metrics:
        selectivity = tissue_sim.calculate_selectivity_index(
            results, metric, CellType.BREAST_CANCER_MCF7
        )
        
        print(f"\n  {label}:")
        for tissue, index in selectivity.items():
            activity = "agonist" if index > 0.5 else "antagonist"
            print(f"    {tissue:20s}: {index:.2f}x ({activity})")
    
    return results


def example_5_ligand_kinetics():
    """compare kinetics of different ligands"""
    print("\n=== example 5: ligand kinetics ===\n")
    
    # fast-acting vs slow-acting ligands
    ligands_info = [
        ('estradiol', 'rapid onset, short half-life'),
        ('tamoxifen', 'moderate onset, long half-life'),
        ('fulvestrant', 'slow onset, very long half-life')
    ]
    
    print("comparing ligand kinetics at 10 mg dose:\n")
    
    for ligand, description in ligands_info:
        sim = ERSimulation(CellType.BREAST_CANCER_MCF7, duration_hours=12.0)
        sim.add_treatment(ligand, dose_mg=10.0)
        df = sim.run()
        
        # find time to peak transcription
        peak_idx = df['er_alpha_transcriptional_activity'].idxmax()
        time_to_peak = df.loc[peak_idx, 'time_hours']
        peak_value = df.loc[peak_idx, 'er_alpha_transcriptional_activity']
        
        print(f"{ligand:12s} ({description})")
        print(f"  time to peak: {time_to_peak:.1f} hours")
        print(f"  peak activity: {peak_value:.3f}\n")


def example_6_phosphorylation_effects():
    """examine kinase-mediated er phosphorylation"""
    print("\n=== example 6: phosphorylation and crosstalk ===\n")
    
    # simulate with estradiol
    sim = ERSimulation(CellType.BREAST_CANCER_MCF7, duration_hours=24.0)
    sim.add_treatment('estradiol', dose_mg=1.0)
    df = sim.run()
    
    # analyze crosstalk between pathways
    print("pathway crosstalk analysis:")
    print(f"  peak mapk activation: {df['nongenomic_mapk_active'].max():.3f}")
    print(f"  peak akt activation: {df['nongenomic_akt_active'].max():.3f}")
    
    if 'er_alpha_phosphorylation' in df.columns:
        print(f"  peak er phosphorylation: {df['er_alpha_phosphorylation'].max():.3f}")
    
    # correlate kinase activity with transcription
    correlation = df['nongenomic_mapk_active'].corr(df['er_alpha_transcriptional_activity'])
    print(f"\ncorrelation (mapk vs transcription): {correlation:.3f}")
    
    return df


def example_7_target_gene_dynamics():
    """detailed analysis of target gene expression"""
    print("\n=== example 7: target gene dynamics ===\n")
    
    sim = ERSimulation(CellType.BREAST_CANCER_MCF7, duration_hours=48.0)
    sim.add_treatment('estradiol', dose_mg=10.0)
    df = sim.run()
    
    genes = ['pr', 'cyclin_d1', 'bcl2', 'tff1', 'vegf', 'c_myc']
    
    print("target gene expression profiles:")
    print("\ngene          peak mrna    peak protein    time to peak (h)")
    print("-" * 60)
    
    for gene in genes:
        mrna_col = f'gene_{gene}_mrna'
        protein_col = f'gene_{gene}_protein'
        
        if mrna_col in df.columns:
            peak_mrna = df[mrna_col].max()
            peak_protein = df[protein_col].max()
            time_to_peak = df.loc[df[protein_col].idxmax(), 'time_hours']
            
            print(f"{gene:12s}  {peak_mrna:10.3f}  {peak_protein:12.3f}  {time_to_peak:12.1f}")
    
    return df


def example_8_er_alpha_beta_balance():
    """examine er-α and er-β balance"""
    print("\n=== example 8: er-α/er-β balance ===\n")
    
    # compare cell types with different er ratios
    cell_types = [
        (CellType.BREAST_CANCER_MCF7, "mcf-7 (α-dominant)"),
        (CellType.OSTEOBLAST, "osteoblast (β-dominant)")
    ]
    
    print("effect of er-α/er-β ratio on estradiol response:\n")
    
    for cell_type, label in cell_types:
        sim = ERSimulation(cell_type, duration_hours=24.0)
        sim.add_treatment('estradiol', dose_mg=1.0)
        df = sim.run()
        
        alpha_activity = df['er_alpha_transcriptional_activity'].max()
        proliferation = df['cell_proliferation_rate'].mean()
        
        print(f"{label}:")
        print(f"  er-α activity: {alpha_activity:.3f}")
        print(f"  proliferation rate: {proliferation:.4f}\n")


def main():
    """run all examples"""
    print("\n" + "="*70)
    print("estrogen receptor signaling simulation examples")
    print("="*70)
    
    # run examples
    example_1_basic_timecourse()
    example_2_serm_comparison()
    example_3_dose_response()
    example_4_tissue_specificity()
    example_5_ligand_kinetics()
    example_6_phosphorylation_effects()
    example_7_target_gene_dynamics()
    example_8_er_alpha_beta_balance()
    
    print("\n" + "="*70)
    print("all examples completed!")
    print("="*70 + "\n")


if __name__ == "__main__":
    main()
