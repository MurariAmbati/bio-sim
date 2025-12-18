"""
estrogen receptor signaling pathways
genomic and non-genomic mechanisms
"""

import numpy as np
from typing import Dict, List, Tuple
from dataclasses import dataclass


@dataclass
class PathwayParameters:
    """parameters for er signaling cascades"""
    # mapk/erk pathway
    mapk_activation_rate: float = 0.3
    mapk_deactivation_rate: float = 0.5
    
    # pi3k/akt pathway
    pi3k_activation_rate: float = 0.4
    akt_phosphorylation_rate: float = 0.6
    
    # src kinase pathway
    src_activation_rate: float = 0.5
    
    # gene expression
    mrna_synthesis_rate: float = 1.0
    mrna_degradation_rate: float = 0.1
    protein_synthesis_rate: float = 0.5
    protein_degradation_rate: float = 0.05
    
    # calcium signaling
    calcium_release_rate: float = 0.8
    calcium_reuptake_rate: float = 1.2


class GenomicPathway:
    """
    classical genomic pathway (hours timescale)
    er -> nucleus -> ere binding -> transcription
    """
    
    def __init__(self, params: PathwayParameters = None):
        self.params = params or PathwayParameters()
        
        # target genes regulated by er
        self.target_genes = {
            'pr': {'mrna': 0.0, 'protein': 0.0},  # progesterone receptor
            'cyclin_d1': {'mrna': 0.0, 'protein': 0.0},  # cell cycle
            'bcl2': {'mrna': 0.0, 'protein': 0.0},  # anti-apoptotic
            'tff1': {'mrna': 0.0, 'protein': 0.0},  # ps2 protein
            'vegf': {'mrna': 0.0, 'protein': 0.0},  # angiogenesis
            'c_myc': {'mrna': 0.0, 'protein': 0.0},  # proliferation
        }
        
        # epigenetic state
        self.chromatin_accessibility = {gene: 0.5 for gene in self.target_genes}
        
    def transcribe(self, gene: str, er_activity: float, dt: float):
        """transcription of er target gene"""
        if gene not in self.target_genes:
            return
        
        accessibility = self.chromatin_accessibility[gene]
        
        # mrna synthesis (depends on er activity and chromatin state)
        synthesis = (self.params.mrna_synthesis_rate * er_activity * 
                    accessibility * dt)
        
        # mrna degradation
        degradation = (self.target_genes[gene]['mrna'] * 
                      self.params.mrna_degradation_rate * dt)
        
        self.target_genes[gene]['mrna'] += synthesis - degradation
        self.target_genes[gene]['mrna'] = max(0, self.target_genes[gene]['mrna'])
    
    def translate(self, gene: str, dt: float):
        """translation of mrna to protein"""
        if gene not in self.target_genes:
            return
        
        # protein synthesis
        synthesis = (self.target_genes[gene]['mrna'] * 
                    self.params.protein_synthesis_rate * dt)
        
        # protein degradation
        degradation = (self.target_genes[gene]['protein'] * 
                      self.params.protein_degradation_rate * dt)
        
        self.target_genes[gene]['protein'] += synthesis - degradation
        self.target_genes[gene]['protein'] = max(0, self.target_genes[gene]['protein'])
    
    def update_chromatin(self, gene: str, histone_acetylation: float):
        """epigenetic regulation via histone modifications"""
        if gene in self.chromatin_accessibility:
            # acetylation increases accessibility
            self.chromatin_accessibility[gene] = 0.5 + 0.5 * histone_acetylation
    
    def update(self, er_activity: float, dt: float) -> Dict:
        """update all target genes"""
        for gene in self.target_genes:
            self.transcribe(gene, er_activity, dt)
            self.translate(gene, dt)
        
        return {gene: data.copy() for gene, data in self.target_genes.items()}


class NonGenomicPathway:
    """
    rapid non-genomic pathways (minutes timescale)
    membrane er -> kinase cascades -> cellular effects
    """
    
    def __init__(self, params: PathwayParameters = None):
        self.params = params or PathwayParameters()
        
        # signaling molecules
        self.mapk_active = 0.0
        self.erk_active = 0.0
        self.pi3k_active = 0.0
        self.akt_active = 0.0
        self.akt_phospho = 0.0
        self.src_active = 0.0
        
        # calcium dynamics
        self.calcium_cytosol = 0.1  # μM
        self.calcium_er_store = 10.0  # μM
        
        # nitric oxide
        self.enos_active = 0.0
        self.no_concentration = 0.0
    
    def activate_mapk(self, membrane_er_activity: float, dt: float):
        """mapk/erk cascade activation"""
        # activation via membrane er
        activation = (self.params.mapk_activation_rate * 
                     membrane_er_activity * dt)
        
        # deactivation by phosphatases
        deactivation = self.mapk_active * self.params.mapk_deactivation_rate * dt
        
        self.mapk_active += activation - deactivation
        self.mapk_active = np.clip(self.mapk_active, 0, 1)
        
        # downstream erk activation
        self.erk_active = self.mapk_active * 0.8
    
    def activate_pi3k_akt(self, membrane_er_activity: float, dt: float):
        """pi3k/akt survival pathway"""
        # pi3k activation
        activation = (self.params.pi3k_activation_rate * 
                     membrane_er_activity * dt)
        
        deactivation = self.pi3k_active * 0.3 * dt
        
        self.pi3k_active += activation - deactivation
        self.pi3k_active = np.clip(self.pi3k_active, 0, 1)
        
        # akt phosphorylation by pdk1
        phospho_rate = (self.params.akt_phosphorylation_rate * 
                       self.pi3k_active * dt)
        
        dephospho_rate = self.akt_phospho * 0.4 * dt
        
        self.akt_phospho += phospho_rate - dephospho_rate
        self.akt_phospho = np.clip(self.akt_phospho, 0, 1)
        
        self.akt_active = self.akt_phospho
    
    def activate_src(self, membrane_er_activity: float, dt: float):
        """src kinase activation"""
        activation = self.params.src_activation_rate * membrane_er_activity * dt
        deactivation = self.src_active * 0.6 * dt
        
        self.src_active += activation - deactivation
        self.src_active = np.clip(self.src_active, 0, 1)
    
    def modulate_calcium(self, er_activity: float, dt: float):
        """rapid calcium signaling"""
        # er triggers calcium release from er stores
        release = (self.params.calcium_release_rate * er_activity * 
                  self.calcium_er_store * dt)
        
        # reuptake by serca pumps
        reuptake = (self.params.calcium_reuptake_rate * 
                   self.calcium_cytosol * dt)
        
        self.calcium_cytosol += release - reuptake
        self.calcium_er_store -= release - reuptake
        
        # maintain physiological bounds
        self.calcium_cytosol = np.clip(self.calcium_cytosol, 0.05, 2.0)
        self.calcium_er_store = np.clip(self.calcium_er_store, 0, 20.0)
    
    def activate_enos(self, er_activity: float, akt_activity: float, dt: float):
        """endothelial nitric oxide synthase activation"""
        # enos activated by er and akt phosphorylation
        self.enos_active = 0.5 * er_activity + 0.5 * akt_activity
        
        # no production
        no_production = self.enos_active * 0.5 * dt
        no_degradation = self.no_concentration * 2.0 * dt  # rapid degradation
        
        self.no_concentration += no_production - no_degradation
        self.no_concentration = max(0, self.no_concentration)
    
    def update(self, membrane_er_activity: float, dt: float) -> Dict:
        """update all non-genomic pathways"""
        self.activate_mapk(membrane_er_activity, dt)
        self.activate_pi3k_akt(membrane_er_activity, dt)
        self.activate_src(membrane_er_activity, dt)
        self.modulate_calcium(membrane_er_activity, dt)
        self.activate_enos(membrane_er_activity, self.akt_active, dt)
        
        return {
            'mapk_active': self.mapk_active,
            'erk_active': self.erk_active,
            'pi3k_active': self.pi3k_active,
            'akt_active': self.akt_active,
            'src_active': self.src_active,
            'calcium_cytosol': self.calcium_cytosol,
            'enos_active': self.enos_active,
            'no_concentration': self.no_concentration,
        }


class CrossTalk:
    """
    pathway crosstalk and feedback loops
    """
    
    @staticmethod
    def mapk_to_er_phosphorylation(mapk_activity: float) -> float:
        """mapk phosphorylates er at ser118"""
        return mapk_activity * 0.6
    
    @staticmethod
    def akt_to_er_phosphorylation(akt_activity: float) -> float:
        """akt phosphorylates er at ser167"""
        return akt_activity * 0.5
    
    @staticmethod
    def src_to_er_phosphorylation(src_activity: float) -> float:
        """src phosphorylates er at tyr537"""
        return src_activity * 0.7
    
    @staticmethod
    def calcium_to_calmodulin(calcium: float) -> float:
        """calcium activates calmodulin"""
        # hill equation n=4 (cooperative binding)
        kd = 0.5
        return calcium**4 / (calcium**4 + kd**4)
    
    @staticmethod
    def cyclin_d1_feedback(cyclin_d1_protein: float) -> float:
        """cyclin d1 promotes cell cycle progression"""
        return min(1.0, cyclin_d1_protein / 5.0)
    
    @staticmethod
    def bcl2_survival(bcl2_protein: float) -> float:
        """bcl-2 inhibits apoptosis"""
        return min(1.0, bcl2_protein / 3.0)
