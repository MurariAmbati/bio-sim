"""
parameter management for ar signaling model
"""

import numpy as np
from dataclasses import dataclass, field
from typing import Dict, Optional


@dataclass
class ParameterSet:
    """
    comprehensive parameter set for ar signaling pathway
    
    all rate constants are in units that make sense for the reaction order:
    - unimolecular: s^-1
    - bimolecular: M^-1 s^-1 or (molecules/cell)^-1 s^-1
    """
    
    # receptor properties
    ar_total: float = 1000.0  # total ar molecules per cell
    ar_cytoplasm_fraction: float = 0.7  # initial cytoplasmic fraction
    ar_nucleus_fraction: float = 0.3  # initial nuclear fraction
    
    # ligand concentrations (M)
    testosterone_external: float = 1.0e-9  # 1 nm
    dht_external: float = 1.0e-9  # 1 nm
    
    # ligand binding kinetics
    # testosterone binding
    k_on_testosterone: float = 1.0e6  # M^-1 s^-1
    k_off_testosterone: float = 0.1  # s^-1 (kd ~ 100 nm)
    
    # dht binding (stronger affinity)
    k_on_dht: float = 1.0e6  # M^-1 s^-1
    k_off_dht: float = 0.01  # s^-1 (kd ~ 10 nm)
    
    # ligand conversion
    k_5alpha_reductase: float = 0.01  # s^-1, testosterone → dht
    
    # heat shock protein dissociation
    k_hsp_dissociation: float = 1.0  # s^-1
    
    # dimerization kinetics
    k_dimerization: float = 0.1  # (molecules/cell)^-1 s^-1
    k_dimer_dissociation: float = 0.01  # s^-1
    
    # nuclear transport
    k_nuclear_import: float = 0.5  # s^-1
    k_nuclear_export: float = 0.1  # s^-1
    
    # dna binding (androgen response elements)
    num_are_sites: int = 100  # number of are binding sites
    k_dna_on: float = 0.01  # (molecules/cell)^-1 s^-1
    k_dna_off: float = 0.001  # s^-1
    
    # transcription rates
    k_transcription_basal: float = 0.001  # s^-1
    k_transcription_activated: float = 0.1  # s^-1
    
    # mrna and protein dynamics
    k_translation: float = 0.1  # s^-1
    k_mrna_degradation: float = 0.001  # s^-1 (half-life ~ 11 min)
    k_protein_degradation: float = 0.0001  # s^-1 (half-life ~ 2 hr)
    
    # ar degradation
    k_ar_degradation: float = 0.0001  # s^-1 (half-life ~ 2 hr)
    k_ar_ubiquitination: float = 0.00005  # s^-1
    
    # coregulator recruitment
    k_coactivator_on: float = 0.05  # (molecules/cell)^-1 s^-1
    k_coactivator_off: float = 0.01  # s^-1
    coactivator_total: float = 500.0  # molecules per cell
    
    k_corepressor_on: float = 0.02  # (molecules/cell)^-1 s^-1
    k_corepressor_off: float = 0.005  # s^-1
    corepressor_total: float = 300.0  # molecules per cell
    
    # chromatin remodeling
    k_chromatin_open: float = 0.01  # s^-1
    k_chromatin_close: float = 0.005  # s^-1
    
    # feedback mechanisms
    k_feedback_psa: float = 0.0001  # negative feedback strength
    k_feedback_fkbp5: float = 0.0002  # fkbp5-mediated feedback
    
    # compartment volumes (liters)
    volume_cell: float = 1.0e-12  # 1 pl (typical cell volume)
    volume_nucleus: float = 1.0e-13  # 100 fl (10% of cell)
    
    # physical constants
    avogadro: float = 6.022e23  # molecules/mol
    
    def __post_init__(self):
        """calculate derived parameters"""
        # calculate kd values
        self.kd_testosterone = self.k_off_testosterone / self.k_on_testosterone
        self.kd_dht = self.k_off_dht / self.k_on_dht
        
        # initial molecule numbers
        self.ar_cytoplasm_init = int(self.ar_total * self.ar_cytoplasm_fraction)
        self.ar_nucleus_init = int(self.ar_total * self.ar_nucleus_fraction)
    
    def to_dict(self) -> Dict[str, float]:
        """convert parameters to dictionary"""
        return {
            k: v for k, v in self.__dict__.items() 
            if not k.startswith('_')
        }
    
    def update(self, **kwargs):
        """update parameters with new values"""
        for key, value in kwargs.items():
            if hasattr(self, key):
                setattr(self, key, value)
            else:
                raise ValueError(f"unknown parameter: {key}")
    
    def get_half_life(self, k_degradation: float) -> float:
        """calculate half-life from degradation rate constant"""
        return np.log(2) / k_degradation
    
    def get_steady_state_ratio(self, k_forward: float, k_reverse: float) -> float:
        """calculate steady-state concentration ratio"""
        return k_forward / k_reverse


@dataclass
class DrugParameters:
    """parameters for drug/antagonist effects"""
    
    name: str = "control"
    
    # competitive antagonist parameters
    is_antagonist: bool = False
    ic50: Optional[float] = None  # M
    ki: Optional[float] = None  # M
    
    # effects on different steps
    binding_inhibition: float = 0.0  # 0-1, fraction of binding blocked
    translocation_inhibition: float = 0.0  # 0-1
    dna_binding_inhibition: float = 0.0  # 0-1
    transcription_inhibition: float = 0.0  # 0-1
    
    # ar degrader (protac)
    is_degrader: bool = False
    k_degradation_enhancement: float = 0.0  # fold increase in degradation
    
    # cyp17 inhibitor (blocks androgen synthesis)
    blocks_testosterone: bool = False
    testosterone_reduction: float = 0.0  # 0-1, fraction reduced
    
    # 5-alpha reductase inhibitor
    blocks_5ar: bool = False
    k_5ar_inhibition: float = 0.0  # 0-1, fraction of conversion blocked


@dataclass
class MutationParameters:
    """parameters for ar mutations"""
    
    name: str = "wildtype"
    
    # mutation types
    is_mutant: bool = False
    
    # gain-of-function mutations (promiscuous activation)
    promiscuous_activation: bool = False
    alternative_ligand_affinity: float = 1.0  # fold change
    
    # increased stability
    stability_increase: float = 1.0  # fold change in half-life
    
    # enhanced transcriptional activity
    transcription_enhancement: float = 1.0  # fold change
    
    # increased nuclear localization
    nuclear_import_enhancement: float = 1.0  # fold change
    
    # ar amplification
    ar_gene_copies: int = 1  # 1 = normal, 2+ = amplification


@dataclass
class CancerParameters:
    """parameters for prostate cancer phenotypes"""
    
    phenotype: str = "normal"  # normal, androgen_dependent, crpc, nepc
    
    # ar expression level
    ar_expression_fold: float = 1.0  # fold change vs normal
    
    # ar variant expression (ar-v7, etc.)
    ar_variant_fraction: float = 0.0  # 0-1, fraction of total ar
    ar_variant_constitutive: bool = False  # ligand-independent activity
    
    # intratumoral androgen synthesis
    intracrine_androgen: float = 0.0  # additional local androgen (M)
    cyp17a1_upregulation: float = 1.0  # fold change
    
    # coactivator overexpression
    coactivator_upregulation: float = 1.0  # fold change
    
    # feedback alteration
    feedback_loss: float = 0.0  # 0-1, loss of negative feedback
    
    # proliferation coupling
    proliferation_rate: float = 0.0  # s^-1, cell division rate
    psa_threshold_proliferation: float = 1000.0  # psa molecules for division


def get_default_parameters() -> ParameterSet:
    """return default parameter set"""
    return ParameterSet()


def get_cancer_parameters(stage: str = "crpc") -> tuple[ParameterSet, CancerParameters]:
    """
    get parameter sets for different cancer stages
    
    stages:
    - androgen_dependent: hormone-sensitive prostate cancer
    - crpc: castration-resistant prostate cancer
    - nepc: neuroendocrine prostate cancer
    """
    params = ParameterSet()
    cancer_params = CancerParameters(phenotype=stage)
    
    if stage == "androgen_dependent":
        # upregulated ar, responsive to adt
        params.ar_total = 2000.0  # 2x expression
        cancer_params.ar_expression_fold = 2.0
        cancer_params.proliferation_rate = 0.00001
        
    elif stage == "crpc":
        # ar amplification, intracrine synthesis, variant expression
        params.ar_total = 5000.0  # 5x expression
        params.k_5alpha_reductase = 0.05  # 5x 5ar activity
        cancer_params.ar_expression_fold = 5.0
        cancer_params.ar_variant_fraction = 0.3
        cancer_params.ar_variant_constitutive = True
        cancer_params.intracrine_androgen = 1.0e-9
        cancer_params.cyp17a1_upregulation = 10.0
        cancer_params.coactivator_upregulation = 3.0
        cancer_params.feedback_loss = 0.5
        cancer_params.proliferation_rate = 0.00002
        
    elif stage == "nepc":
        # ar-independent, neuroendocrine differentiation
        params.ar_total = 100.0  # loss of ar
        cancer_params.ar_expression_fold = 0.1
        cancer_params.feedback_loss = 1.0
        cancer_params.proliferation_rate = 0.00003
    
    return params, cancer_params


def get_drug_parameters(drug_name: str, concentration: float = 10.0e-6) -> DrugParameters:
    """
    get parameters for common ar-targeted drugs
    
    drugs:
    - bicalutamide: first-generation antiandrogen
    - enzalutamide: second-generation antiandrogen
    - apalutamide: second-generation antiandrogen
    - darolutamide: second-generation antiandrogen
    - abiraterone: cyp17 inhibitor
    - finasteride: 5-alpha reductase inhibitor
    - arv110: ar degrader (protac)
    """
    drug = DrugParameters(name=drug_name)
    
    if drug_name == "bicalutamide":
        drug.is_antagonist = True
        drug.ic50 = 0.19e-6  # 190 nm
        # partial antagonist
        drug.binding_inhibition = concentration / (concentration + drug.ic50) * 0.7
        drug.dna_binding_inhibition = 0.3
        
    elif drug_name == "enzalutamide":
        drug.is_antagonist = True
        drug.ic50 = 0.36e-6  # 360 nm
        # potent antagonist, blocks multiple steps
        inhibition = concentration / (concentration + drug.ic50)
        drug.binding_inhibition = inhibition * 0.8
        drug.translocation_inhibition = inhibition * 0.6
        drug.dna_binding_inhibition = inhibition * 0.9
        drug.transcription_inhibition = inhibition * 0.5
        
    elif drug_name == "apalutamide":
        drug.is_antagonist = True
        drug.ic50 = 0.016e-6  # 16 nm
        inhibition = concentration / (concentration + drug.ic50)
        drug.binding_inhibition = inhibition * 0.85
        drug.translocation_inhibition = inhibition * 0.7
        drug.dna_binding_inhibition = inhibition * 0.95
        
    elif drug_name == "darolutamide":
        drug.is_antagonist = True
        drug.ic50 = 0.011e-6  # 11 nm
        inhibition = concentration / (concentration + drug.ic50)
        drug.binding_inhibition = inhibition * 0.9
        drug.dna_binding_inhibition = inhibition * 0.85
        
    elif drug_name == "abiraterone":
        # blocks androgen synthesis
        drug.blocks_testosterone = True
        drug.testosterone_reduction = 0.95  # 95% reduction
        
    elif drug_name == "finasteride":
        # 5-alpha reductase inhibitor
        drug.blocks_5ar = True
        drug.k_5ar_inhibition = 0.9  # 90% inhibition
        
    elif drug_name == "arv110":
        # ar degrader (protac)
        drug.is_degrader = True
        drug.k_degradation_enhancement = 10.0  # 10x faster degradation
    
    return drug
