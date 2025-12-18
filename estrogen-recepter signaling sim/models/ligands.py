"""
estrogen receptor ligands
agonists, antagonists, and selective estrogen receptor modulators (serms)
"""

import numpy as np
from dataclasses import dataclass
from typing import Dict, List, Optional
from enum import Enum


class LigandType(Enum):
    """classification of er ligands"""
    AGONIST = "agonist"
    ANTAGONIST = "antagonist"
    SERM = "serm"  # selective estrogen receptor modulator
    SERD = "serd"  # selective estrogen receptor degrader


@dataclass
class LigandProperties:
    """pharmacological properties of er ligands"""
    name: str
    ligand_type: LigandType
    
    # binding affinity
    kd_er_alpha: float  # nM
    kd_er_beta: float   # nM
    
    # efficacy
    efficacy_er_alpha: float  # 0-1 scale
    efficacy_er_beta: float
    
    # tissue selectivity (relative to uterus)
    bone_selectivity: float = 1.0
    breast_selectivity: float = 1.0
    cardiovascular_selectivity: float = 1.0
    brain_selectivity: float = 1.0
    
    # pharmacokinetics
    absorption_rate: float = 0.5  # 1/min
    distribution_volume: float = 100.0  # L
    clearance_rate: float = 0.05  # 1/min
    half_life: float = 14.0  # hours
    
    # receptor degradation
    induces_degradation: bool = False
    degradation_rate_modifier: float = 1.0


class LigandLibrary:
    """collection of er ligands with their properties"""
    
    @staticmethod
    def get_ligand(name: str) -> Optional[LigandProperties]:
        """retrieve ligand properties by name"""
        library = {
            'estradiol': LigandProperties(
                name='17β-estradiol',
                ligand_type=LigandType.AGONIST,
                kd_er_alpha=0.1,
                kd_er_beta=0.5,
                efficacy_er_alpha=1.0,
                efficacy_er_beta=1.0,
                bone_selectivity=1.0,
                breast_selectivity=1.0,
                cardiovascular_selectivity=1.0,
                brain_selectivity=1.0,
                half_life=1.0,
            ),
            
            'estrone': LigandProperties(
                name='estrone',
                ligand_type=LigandType.AGONIST,
                kd_er_alpha=0.5,
                kd_er_beta=1.0,
                efficacy_er_alpha=0.6,
                efficacy_er_beta=0.5,
                half_life=2.0,
            ),
            
            'estriol': LigandProperties(
                name='estriol',
                ligand_type=LigandType.AGONIST,
                kd_er_alpha=1.0,
                kd_er_beta=2.0,
                efficacy_er_alpha=0.4,
                efficacy_er_beta=0.3,
                half_life=0.5,
            ),
            
            'tamoxifen': LigandProperties(
                name='tamoxifen',
                ligand_type=LigandType.SERM,
                kd_er_alpha=10.0,
                kd_er_beta=20.0,
                efficacy_er_alpha=0.2,
                efficacy_er_beta=0.1,
                bone_selectivity=1.5,  # agonist in bone
                breast_selectivity=0.1,  # antagonist in breast
                cardiovascular_selectivity=1.2,
                brain_selectivity=0.3,
                half_life=120.0,  # 5 days
            ),
            
            'raloxifene': LigandProperties(
                name='raloxifene',
                ligand_type=LigandType.SERM,
                kd_er_alpha=5.0,
                kd_er_beta=10.0,
                efficacy_er_alpha=0.15,
                efficacy_er_beta=0.1,
                bone_selectivity=1.8,  # strong agonist in bone
                breast_selectivity=0.05,  # strong antagonist in breast
                cardiovascular_selectivity=0.8,
                brain_selectivity=0.2,
                half_life=27.0,
            ),
            
            'fulvestrant': LigandProperties(
                name='fulvestrant',
                ligand_type=LigandType.SERD,
                kd_er_alpha=0.5,
                kd_er_beta=1.0,
                efficacy_er_alpha=0.0,
                efficacy_er_beta=0.0,
                bone_selectivity=0.0,
                breast_selectivity=0.0,
                cardiovascular_selectivity=0.0,
                brain_selectivity=0.0,
                induces_degradation=True,
                degradation_rate_modifier=5.0,
                half_life=960.0,  # 40 days
            ),
            
            'ici182780': LigandProperties(
                name='ici 182,780',
                ligand_type=LigandType.ANTAGONIST,
                kd_er_alpha=0.4,
                kd_er_beta=0.8,
                efficacy_er_alpha=0.0,
                efficacy_er_beta=0.0,
                induces_degradation=True,
                degradation_rate_modifier=4.0,
                half_life=480.0,
            ),
            
            'genistein': LigandProperties(
                name='genistein',
                ligand_type=LigandType.SERM,
                kd_er_alpha=100.0,
                kd_er_beta=20.0,  # beta selective
                efficacy_er_alpha=0.3,
                efficacy_er_beta=0.7,
                bone_selectivity=1.3,
                breast_selectivity=0.6,
                cardiovascular_selectivity=1.5,
                brain_selectivity=1.2,
                half_life=8.0,
            ),
            
            'daidzein': LigandProperties(
                name='daidzein',
                ligand_type=LigandType.SERM,
                kd_er_alpha=200.0,
                kd_er_beta=50.0,
                efficacy_er_alpha=0.2,
                efficacy_er_beta=0.5,
                half_life=7.0,
            ),
        }
        
        return library.get(name.lower())
    
    @staticmethod
    def list_ligands() -> List[str]:
        """list all available ligands"""
        return [
            'estradiol', 'estrone', 'estriol',
            'tamoxifen', 'raloxifene', 'fulvestrant',
            'ici182780', 'genistein', 'daidzein'
        ]


class LigandPharmacokinetics:
    """
    adme (absorption, distribution, metabolism, excretion) model
    """
    
    def __init__(self, ligand: LigandProperties):
        self.ligand = ligand
        
        # compartments
        self.gut_concentration = 0.0
        self.plasma_concentration = 0.0
        self.tissue_concentration = 0.0
        
        # metabolites
        self.metabolite_concentration = 0.0
    
    def oral_dose(self, dose_mg: float):
        """administer oral dose"""
        # convert to nM (approximate molecular weight 300 g/mol)
        dose_nmol = (dose_mg / 300.0) * 1e6
        
        # initial gut concentration
        self.gut_concentration += dose_nmol / 1.0  # gut volume ~1L
    
    def update(self, dt: float, tissue_type: str = 'breast') -> Dict:
        """
        update pharmacokinetic model
        dt in minutes
        """
        # absorption from gut to plasma
        absorption = (self.gut_concentration * 
                     self.ligand.absorption_rate * dt)
        self.gut_concentration -= absorption
        self.gut_concentration = max(0, self.gut_concentration)
        
        # distribution to tissues
        distribution = (self.plasma_concentration * 0.1 * dt)
        
        # tissue-specific distribution
        tissue_modifier = {
            'breast': self.ligand.breast_selectivity,
            'bone': self.ligand.bone_selectivity,
            'cardiovascular': self.ligand.cardiovascular_selectivity,
            'brain': self.ligand.brain_selectivity,
        }.get(tissue_type, 1.0)
        
        distribution *= tissue_modifier
        
        # metabolism and clearance
        clearance = (self.plasma_concentration * 
                    self.ligand.clearance_rate * dt)
        
        # update compartments
        self.plasma_concentration += absorption - distribution - clearance
        self.plasma_concentration = max(0, self.plasma_concentration)
        
        self.tissue_concentration += distribution - (self.tissue_concentration * 
                                                     self.ligand.clearance_rate * 
                                                     0.5 * dt)
        self.tissue_concentration = max(0, self.tissue_concentration)
        
        # metabolites
        self.metabolite_concentration += clearance * 0.8
        
        return {
            'gut': self.gut_concentration,
            'plasma': self.plasma_concentration,
            'tissue': self.tissue_concentration,
            'metabolites': self.metabolite_concentration,
        }
    
    def get_effective_concentration(self, tissue_type: str = 'breast') -> float:
        """get concentration at target tissue"""
        return self.tissue_concentration
