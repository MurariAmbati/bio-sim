"""
PPAR (Peroxisome Proliferator-Activated Receptor) Signaling Model
Comprehensive simulation of PPAR-α, PPAR-γ, and PPAR-δ signaling pathways
"""

import numpy as np
from scipy.integrate import odeint
from dataclasses import dataclass
from typing import Dict, List, Tuple
import pandas as pd


@dataclass
class PPARParameters:
    """Parameters for PPAR signaling model"""
    
    # PPAR isoform expression levels
    ppar_alpha_base: float = 1.0
    ppar_gamma_base: float = 1.0
    ppar_delta_base: float = 1.0
    
    # Ligand binding rates
    k_ligand_bind_alpha: float = 0.5
    k_ligand_bind_gamma: float = 0.6
    k_ligand_bind_delta: float = 0.4
    k_ligand_unbind: float = 0.1
    
    # RXR heterodimerization
    k_rxr_bind: float = 0.8
    k_rxr_unbind: float = 0.2
    rxr_total: float = 2.0
    
    # DNA binding and transcription
    k_dna_bind: float = 0.7
    k_dna_unbind: float = 0.15
    k_transcription: float = 1.0
    
    # Lipid metabolism genes (PPAR-α primary)
    k_fatty_acid_oxidation: float = 0.9  # CPT1, ACOX1
    k_ketogenesis: float = 0.7  # HMGCS2
    k_lipoprotein_metabolism: float = 0.8  # APOA1, LPL
    
    # Insulin sensitivity genes (PPAR-γ primary)
    k_glucose_uptake: float = 0.85  # GLUT4
    k_adipogenesis: float = 0.75  # CEBPA, FABP4
    k_lipid_storage: float = 0.7  # PLIN1
    
    # Inflammation genes (all isoforms)
    k_nfkb_inhibition: float = 0.6  # Anti-inflammatory
    k_cytokine_suppression: float = 0.5  # IL-6, TNF-α suppression
    k_m2_polarization: float = 0.65  # Macrophage M2 polarization (PPAR-δ)
    
    # mRNA and protein dynamics
    k_mrna_degradation: float = 0.3
    k_protein_synthesis: float = 0.5
    k_protein_degradation: float = 0.2
    
    # Feedback mechanisms
    k_feedback_inhibition: float = 0.3
    k_crosstalk_alpha_gamma: float = 0.15
    k_crosstalk_gamma_delta: float = 0.12
    
    # Metabolic effects
    free_fatty_acids: float = 1.0
    glucose_concentration: float = 5.0  # mM
    insulin_level: float = 1.0
    inflammatory_cytokines: float = 0.5


class PPARSignalingModel:
    """
    Comprehensive PPAR signaling pathway model
    
    State variables:
    0-2: Free PPAR (α, γ, δ)
    3-5: Ligand-bound PPAR (α, γ, δ)
    6-8: PPAR-RXR heterodimers (α, γ, δ)
    9-11: DNA-bound PPAR-RXR (α, γ, δ)
    12-14: Target gene mRNA - lipid oxidation genes (α-driven)
    15-17: Target gene mRNA - insulin sensitivity genes (γ-driven)
    18-20: Target gene mRNA - anti-inflammatory genes (all isoforms)
    21-23: Metabolic proteins (fatty acid oxidation, glucose uptake, anti-inflammatory)
    24: Free RXR
    25: Lipid accumulation
    26: Insulin sensitivity index
    27: Inflammatory state
    """
    
    def __init__(self, params: PPARParameters = None):
        self.params = params if params is not None else PPARParameters()
        self.state_names = self._get_state_names()
        
    def _get_state_names(self) -> List[str]:
        return [
            'PPAR-α', 'PPAR-γ', 'PPAR-δ',
            'PPAR-α*', 'PPAR-γ*', 'PPAR-δ*',
            'PPAR-α-RXR', 'PPAR-γ-RXR', 'PPAR-δ-RXR',
            'PPAR-α-RXR-DNA', 'PPAR-γ-RXR-DNA', 'PPAR-δ-RXR-DNA',
            'mRNA_FAO_α', 'mRNA_FAO_γ', 'mRNA_FAO_δ',
            'mRNA_IS_α', 'mRNA_IS_γ', 'mRNA_IS_δ',
            'mRNA_AI_α', 'mRNA_AI_γ', 'mRNA_AI_δ',
            'Protein_FAO', 'Protein_IS', 'Protein_AI',
            'RXR_free',
            'Lipid_Accumulation',
            'Insulin_Sensitivity',
            'Inflammatory_State'
        ]
    
    def derivatives(self, state: np.ndarray, t: float, ligand_alpha: float, 
                   ligand_gamma: float, ligand_delta: float) -> np.ndarray:
        """Calculate derivatives for ODE system"""
        p = self.params
        ds = np.zeros(28)
        
        # Unpack state variables
        ppar_free = state[0:3]  # α, γ, δ
        ppar_ligand = state[3:6]
        ppar_rxr = state[6:9]
        ppar_rxr_dna = state[9:12]
        mrna_fao = state[12:15]  # Fatty acid oxidation
        mrna_is = state[15:18]   # Insulin sensitivity
        mrna_ai = state[18:21]   # Anti-inflammatory
        protein_fao = state[21]
        protein_is = state[22]
        protein_ai = state[23]
        rxr_free = state[24]
        lipid_accum = state[25]
        insulin_sens = state[26]
        inflam_state = state[27]
        
        ligands = np.array([ligand_alpha, ligand_gamma, ligand_delta])
        k_bind = np.array([p.k_ligand_bind_alpha, p.k_ligand_bind_gamma, p.k_ligand_bind_delta])
        
        # PPAR dynamics - ligand binding
        ds[0:3] = (
            np.array([p.ppar_alpha_base, p.ppar_gamma_base, p.ppar_delta_base]) 
            - k_bind * ligands * ppar_free
            + p.k_ligand_unbind * ppar_ligand
            - p.k_feedback_inhibition * ppar_free * ppar_rxr_dna
        )
        
        # Ligand-bound PPAR
        ds[3:6] = (
            k_bind * ligands * ppar_free
            - p.k_ligand_unbind * ppar_ligand
            - p.k_rxr_bind * ppar_ligand * rxr_free
            + p.k_rxr_unbind * ppar_rxr
        )
        
        # PPAR-RXR heterodimers
        ds[6:9] = (
            p.k_rxr_bind * ppar_ligand * rxr_free
            - p.k_rxr_unbind * ppar_rxr
            - p.k_dna_bind * ppar_rxr
            + p.k_dna_unbind * ppar_rxr_dna
        )
        
        # DNA-bound PPAR-RXR complexes
        ds[9:12] = (
            p.k_dna_bind * ppar_rxr
            - p.k_dna_unbind * ppar_rxr_dna
        )
        
        # Transcription factors with isoform specificity
        # PPAR-α: Primary driver of fatty acid oxidation
        transcription_fao = np.array([
            p.k_transcription * p.k_fatty_acid_oxidation * ppar_rxr_dna[0],
            p.k_transcription * p.k_fatty_acid_oxidation * 0.3 * ppar_rxr_dna[1],
            p.k_transcription * p.k_fatty_acid_oxidation * 0.5 * ppar_rxr_dna[2]
        ])
        
        # PPAR-γ: Primary driver of insulin sensitivity
        transcription_is = np.array([
            p.k_transcription * p.k_glucose_uptake * 0.3 * ppar_rxr_dna[0],
            p.k_transcription * p.k_glucose_uptake * ppar_rxr_dna[1],
            p.k_transcription * p.k_glucose_uptake * 0.4 * ppar_rxr_dna[2]
        ])
        
        # Anti-inflammatory: All isoforms contribute
        transcription_ai = np.array([
            p.k_transcription * p.k_nfkb_inhibition * 0.8 * ppar_rxr_dna[0],
            p.k_transcription * p.k_nfkb_inhibition * 0.7 * ppar_rxr_dna[1],
            p.k_transcription * p.k_nfkb_inhibition * ppar_rxr_dna[2]  # δ strong in inflammation
        ])
        
        # mRNA dynamics
        ds[12:15] = transcription_fao - p.k_mrna_degradation * mrna_fao
        ds[15:18] = transcription_is - p.k_mrna_degradation * mrna_is
        ds[18:21] = transcription_ai - p.k_mrna_degradation * mrna_ai
        
        # Protein synthesis and degradation
        total_mrna_fao = np.sum(mrna_fao)
        total_mrna_is = np.sum(mrna_is)
        total_mrna_ai = np.sum(mrna_ai)
        
        ds[21] = p.k_protein_synthesis * total_mrna_fao - p.k_protein_degradation * protein_fao
        ds[22] = p.k_protein_synthesis * total_mrna_is - p.k_protein_degradation * protein_is
        ds[23] = p.k_protein_synthesis * total_mrna_ai - p.k_protein_degradation * protein_ai
        
        # Free RXR dynamics
        total_rxr_bound = np.sum(ppar_rxr) + np.sum(ppar_rxr_dna)
        ds[24] = -p.k_rxr_bind * rxr_free * np.sum(ppar_ligand) + p.k_rxr_unbind * np.sum(ppar_rxr)
        
        # Metabolic outputs
        # Lipid accumulation (decreased by FAO, increased by FFA)
        ds[25] = p.free_fatty_acids * 0.5 - protein_fao * lipid_accum * 0.3
        
        # Insulin sensitivity (increased by PPAR-γ activity, decreased by inflammation)
        ds[26] = protein_is * 0.4 - inflam_state * insulin_sens * 0.3 - p.k_protein_degradation * insulin_sens
        
        # Inflammatory state (suppressed by all PPARs)
        ds[27] = p.inflammatory_cytokines * 0.3 - protein_ai * inflam_state * 0.5
        
        return ds
    
    def get_initial_conditions(self) -> np.ndarray:
        """Set initial conditions for the model"""
        state0 = np.zeros(28)
        
        # Initial PPAR levels
        state0[0:3] = [self.params.ppar_alpha_base, 
                       self.params.ppar_gamma_base, 
                       self.params.ppar_delta_base]
        
        # Initial RXR
        state0[24] = self.params.rxr_total
        
        # Initial metabolic states
        state0[25] = 1.0  # Lipid accumulation
        state0[26] = 1.0  # Insulin sensitivity
        state0[27] = self.params.inflammatory_cytokines  # Inflammatory state
        
        return state0
    
    def simulate(self, t_span: Tuple[float, float], n_points: int = 1000,
                ligand_alpha: float = 0.0, ligand_gamma: float = 0.0, 
                ligand_delta: float = 0.0) -> pd.DataFrame:
        """
        Simulate PPAR signaling dynamics
        
        Parameters:
        -----------
        t_span: Tuple of (start_time, end_time)
        n_points: Number of time points
        ligand_alpha: Concentration of PPAR-α ligand (e.g., fibrates)
        ligand_gamma: Concentration of PPAR-γ ligand (e.g., thiazolidinediones)
        ligand_delta: Concentration of PPAR-δ ligand (e.g., GW501516)
        
        Returns:
        --------
        DataFrame with simulation results
        """
        t = np.linspace(t_span[0], t_span[1], n_points)
        state0 = self.get_initial_conditions()
        
        solution = odeint(self.derivatives, state0, t, 
                         args=(ligand_alpha, ligand_gamma, ligand_delta))
        
        # Create DataFrame
        df = pd.DataFrame(solution, columns=self.state_names)
        df['Time'] = t
        
        # Add derived metrics
        df['Total_PPAR_Activity'] = (df['PPAR-α-RXR-DNA'] + 
                                      df['PPAR-γ-RXR-DNA'] + 
                                      df['PPAR-δ-RXR-DNA'])
        
        df['Metabolic_Health'] = (df['Insulin_Sensitivity'] * 2 - 
                                  df['Lipid_Accumulation'] - 
                                  df['Inflammatory_State']) / 4
        
        return df
    
    def simulate_drug_response(self, drug_type: str, dose_range: np.ndarray,
                              t_span: Tuple[float, float] = (0, 100),
                              n_points: int = 500) -> Dict[str, pd.DataFrame]:
        """
        Simulate response to different PPAR agonists
        
        Parameters:
        -----------
        drug_type: 'alpha' (fibrates), 'gamma' (TZDs), 'delta', or 'pan' (pan-agonist)
        dose_range: Array of drug concentrations
        
        Returns:
        --------
        Dictionary with DataFrames for each dose
        """
        results = {}
        
        for dose in dose_range:
            if drug_type == 'alpha':
                ligands = (dose, 0.0, 0.0)
            elif drug_type == 'gamma':
                ligands = (0.0, dose, 0.0)
            elif drug_type == 'delta':
                ligands = (0.0, 0.0, dose)
            elif drug_type == 'pan':
                ligands = (dose, dose, dose)
            else:
                raise ValueError(f"Unknown drug type: {drug_type}")
            
            df = self.simulate(t_span, n_points, *ligands)
            results[f'dose_{dose:.2f}'] = df
        
        return results
    
    def get_steady_state_metrics(self, df: pd.DataFrame) -> Dict[str, float]:
        """Extract steady-state metrics from simulation"""
        # Take last 10% of simulation as steady state
        steady_state = df.iloc[-len(df)//10:]
        
        metrics = {
            'PPAR-α Activity': steady_state['PPAR-α-RXR-DNA'].mean(),
            'PPAR-γ Activity': steady_state['PPAR-γ-RXR-DNA'].mean(),
            'PPAR-δ Activity': steady_state['PPAR-δ-RXR-DNA'].mean(),
            'Fatty Acid Oxidation': steady_state['Protein_FAO'].mean(),
            'Insulin Sensitivity': steady_state['Insulin_Sensitivity'].mean(),
            'Anti-inflammatory': steady_state['Protein_AI'].mean(),
            'Lipid Accumulation': steady_state['Lipid_Accumulation'].mean(),
            'Inflammatory State': steady_state['Inflammatory_State'].mean(),
            'Metabolic Health': steady_state['Metabolic_Health'].mean(),
        }
        
        return metrics


def compare_isoform_specificity() -> pd.DataFrame:
    """
    Compare the specificity of different PPAR isoforms for their primary functions
    """
    model = PPARSignalingModel()
    
    # Test individual isoform activation
    results = []
    
    for isoform, idx in [('α', 0), ('γ', 1), ('δ', 2)]:
        ligands = [0.0, 0.0, 0.0]
        ligands[idx] = 1.0
        
        df = model.simulate((0, 100), 500, *ligands)
        metrics = model.get_steady_state_metrics(df)
        metrics['Isoform'] = f'PPAR-{isoform}'
        results.append(metrics)
    
    return pd.DataFrame(results)
