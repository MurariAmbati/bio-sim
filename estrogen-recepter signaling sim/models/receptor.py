"""
estrogen receptor (er) models
implements molecular dynamics of er-alpha and er-beta
"""

import numpy as np
from dataclasses import dataclass
from typing import Dict, List, Optional
from enum import Enum


class ReceptorType(Enum):
    """estrogen receptor subtypes"""
    ER_ALPHA = "er_alpha"
    ER_BETA = "er_beta"


class ReceptorState(Enum):
    """receptor conformational states"""
    UNBOUND = "unbound"
    LIGAND_BOUND = "ligand_bound"
    DIMERIZED = "dimerized"
    DNA_BOUND = "dna_bound"
    DEGRADED = "degraded"


@dataclass
class ReceptorParameters:
    """parameters for er dynamics"""
    # synthesis rates (molecules/cell/min)
    synthesis_rate: float = 0.5
    
    # degradation rates (1/min)
    degradation_rate: float = 0.01
    
    # binding affinities (nM)
    estradiol_kd: float = 0.1  # high affinity
    tamoxifen_kd: float = 10.0  # moderate affinity
    raloxifene_kd: float = 5.0
    
    # transition rates (1/min)
    dimerization_rate: float = 0.5
    dna_binding_rate: float = 0.3
    nuclear_import_rate: float = 0.2
    
    # conformational parameters
    agonist_activation_prob: float = 0.9
    antagonist_activation_prob: float = 0.1
    serm_tissue_specific_factor: float = 0.5


class EstrogenReceptor:
    """
    models estrogen receptor signaling dynamics
    includes genomic and non-genomic pathways
    """
    
    def __init__(self, receptor_type: ReceptorType, params: Optional[ReceptorParameters] = None):
        self.receptor_type = receptor_type
        self.params = params or ReceptorParameters()
        
        # state variables
        self.state = ReceptorState.UNBOUND
        self.concentration = 0.0  # nM
        self.nuclear_concentration = 0.0
        self.cytoplasmic_concentration = 0.0
        
        # ligand binding
        self.bound_ligand = None
        self.ligand_affinity_modifier = 1.0
        
        # post-translational modifications
        self.phosphorylation_sites = {
            'ser118': 0.0,  # mapk target
            'ser167': 0.0,  # akt target
            'ser305': 0.0,  # pkc target
            'tyr537': 0.0   # src target
        }
        
        # interaction partners
        self.coactivators = {}
        self.corepressors = {}
        
    def bind_ligand(self, ligand_name: str, concentration: float) -> float:
        """
        calculate ligand binding based on concentration and kd
        returns fraction bound
        """
        kd_map = {
            'estradiol': self.params.estradiol_kd,
            'tamoxifen': self.params.tamoxifen_kd,
            'raloxifene': self.params.raloxifene_kd,
        }
        
        kd = kd_map.get(ligand_name, 10.0) * self.ligand_affinity_modifier
        
        # hill equation (n=1)
        fraction_bound = concentration / (concentration + kd)
        
        if fraction_bound > 0.5:
            self.bound_ligand = ligand_name
            self.state = ReceptorState.LIGAND_BOUND
            
        return fraction_bound
    
    def dimerize(self, dt: float) -> bool:
        """receptor dimerization following ligand binding"""
        if self.state == ReceptorState.LIGAND_BOUND:
            # cooperative dimerization
            dimerization_prob = 1 - np.exp(-self.params.dimerization_rate * dt)
            
            if np.random.random() < dimerization_prob:
                self.state = ReceptorState.DIMERIZED
                return True
                
        return False
    
    def bind_dna(self, dt: float, ere_availability: float = 1.0) -> bool:
        """
        dna binding at estrogen response elements (eres)
        """
        if self.state == ReceptorState.DIMERIZED:
            binding_prob = (1 - np.exp(-self.params.dna_binding_rate * dt)) * ere_availability
            
            if np.random.random() < binding_prob:
                self.state = ReceptorState.DNA_BOUND
                return True
                
        return False
    
    def phosphorylate(self, site: str, kinase_activity: float):
        """
        ligand-independent activation via phosphorylation
        """
        if site in self.phosphorylation_sites:
            # saturating kinetics
            self.phosphorylation_sites[site] = min(1.0, 
                self.phosphorylation_sites[site] + kinase_activity * 0.1)
            
            # phosphorylation can activate receptor
            if sum(self.phosphorylation_sites.values()) > 1.5:
                self.ligand_affinity_modifier *= 1.2
    
    def recruit_coregulator(self, regulator_name: str, is_activator: bool, 
                           strength: float):
        """recruit coactivators or corepressors"""
        if is_activator:
            self.coactivators[regulator_name] = strength
        else:
            self.corepressors[regulator_name] = strength
    
    def calculate_transcriptional_activity(self) -> float:
        """
        compute transcriptional output based on receptor state
        """
        if self.state != ReceptorState.DNA_BOUND:
            return 0.0
        
        # base activity from ligand type
        if self.bound_ligand == 'estradiol':
            base_activity = self.params.agonist_activation_prob
        elif self.bound_ligand in ['tamoxifen', 'raloxifene']:
            base_activity = (self.params.antagonist_activation_prob * 
                           self.params.serm_tissue_specific_factor)
        else:
            base_activity = 0.1
        
        # modulation by coregulators
        coactivator_effect = sum(self.coactivators.values())
        corepressor_effect = sum(self.corepressors.values())
        
        # modulation by phosphorylation
        phospho_effect = sum(self.phosphorylation_sites.values()) * 0.1
        
        total_activity = base_activity * (1 + coactivator_effect + phospho_effect) / (1 + corepressor_effect)
        
        return np.clip(total_activity, 0, 2.0)
    
    def update(self, dt: float) -> Dict:
        """
        update receptor state over timestep dt
        returns state summary
        """
        # synthesis
        synthesis = self.params.synthesis_rate * dt
        
        # degradation (ubiquitin-proteasome pathway)
        degradation = self.concentration * self.params.degradation_rate * dt
        
        # ligand-induced degradation
        if self.state == ReceptorState.LIGAND_BOUND:
            degradation *= 2.0  # accelerated degradation
        
        # update concentration
        self.concentration += synthesis - degradation
        self.concentration = max(0, self.concentration)
        
        return {
            'concentration': self.concentration,
            'state': self.state.value,
            'bound_ligand': self.bound_ligand,
            'transcriptional_activity': self.calculate_transcriptional_activity(),
            'phosphorylation': sum(self.phosphorylation_sites.values()),
        }
