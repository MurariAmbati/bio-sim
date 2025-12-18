"""
cellular context for er signaling
integrates receptor, pathway, and cellular responses
"""

import numpy as np
from typing import Dict, List, Optional
from .receptor import EstrogenReceptor, ReceptorType, ReceptorParameters
from .pathways import GenomicPathway, NonGenomicPathway, PathwayParameters, CrossTalk
from .ligands import LigandPharmacokinetics, LigandLibrary


class CellType:
    """cell type definitions"""
    BREAST_EPITHELIAL = "breast_epithelial"
    BREAST_CANCER_MCF7 = "mcf7"
    BREAST_CANCER_T47D = "t47d"
    ENDOMETRIAL = "endometrial"
    OSTEOBLAST = "osteoblast"
    ENDOTHELIAL = "endothelial"
    NEURON = "neuron"


class CellularState:
    """
    state variables for a single cell
    """
    
    def __init__(self):
        # cell cycle
        self.cell_cycle_phase = "g1"
        self.cell_cycle_progress = 0.0
        
        # survival
        self.apoptosis_index = 0.0
        self.survival_signal = 1.0
        
        # proliferation
        self.proliferation_rate = 0.0
        self.division_count = 0
        
        # metabolism
        self.atp_level = 1.0
        self.glucose_uptake = 1.0
        
        # stress
        self.oxidative_stress = 0.0
        self.dna_damage = 0.0


class EstrogenResponsiveCell:
    """
    complete cell model with er signaling
    """
    
    def __init__(self, cell_type: str, 
                 er_alpha_params: Optional[ReceptorParameters] = None,
                 er_beta_params: Optional[ReceptorParameters] = None,
                 pathway_params: Optional[PathwayParameters] = None):
        
        self.cell_type = cell_type
        self.state = CellularState()
        
        # receptors
        self.er_alpha = EstrogenReceptor(ReceptorType.ER_ALPHA, er_alpha_params)
        self.er_beta = EstrogenReceptor(ReceptorType.ER_BETA, er_beta_params)
        
        # set initial receptor levels based on cell type
        self._initialize_receptors()
        
        # pathways
        self.genomic = GenomicPathway(pathway_params)
        self.non_genomic = NonGenomicPathway(pathway_params)
        
        # ligand pharmacokinetics
        self.ligand_pk = {}
        
        # time tracking
        self.time = 0.0
    
    def _initialize_receptors(self):
        """set cell-type specific receptor levels"""
        receptor_levels = {
            CellType.BREAST_CANCER_MCF7: {'alpha': 10.0, 'beta': 2.0},
            CellType.BREAST_CANCER_T47D: {'alpha': 8.0, 'beta': 3.0},
            CellType.BREAST_EPITHELIAL: {'alpha': 5.0, 'beta': 5.0},
            CellType.ENDOMETRIAL: {'alpha': 12.0, 'beta': 4.0},
            CellType.OSTEOBLAST: {'alpha': 3.0, 'beta': 7.0},
            CellType.ENDOTHELIAL: {'alpha': 4.0, 'beta': 6.0},
            CellType.NEURON: {'alpha': 2.0, 'beta': 8.0},
        }
        
        levels = receptor_levels.get(self.cell_type, {'alpha': 5.0, 'beta': 5.0})
        self.er_alpha.concentration = levels['alpha']
        self.er_beta.concentration = levels['beta']
    
    def add_ligand(self, ligand_name: str, dose_mg: float):
        """add ligand to cell environment"""
        ligand = LigandLibrary.get_ligand(ligand_name)
        if ligand:
            pk = LigandPharmacokinetics(ligand)
            pk.oral_dose(dose_mg)
            self.ligand_pk[ligand_name] = pk
    
    def update_ligand_binding(self, dt: float):
        """update all ligand-receptor interactions"""
        tissue_map = {
            CellType.BREAST_CANCER_MCF7: 'breast',
            CellType.BREAST_CANCER_T47D: 'breast',
            CellType.BREAST_EPITHELIAL: 'breast',
            CellType.OSTEOBLAST: 'bone',
            CellType.ENDOTHELIAL: 'cardiovascular',
            CellType.NEURON: 'brain',
        }
        tissue = tissue_map.get(self.cell_type, 'breast')
        
        for ligand_name, pk in self.ligand_pk.items():
            # update pharmacokinetics
            pk.update(dt, tissue)
            
            # get tissue concentration
            conc = pk.get_effective_concentration(tissue)
            
            # bind to receptors
            self.er_alpha.bind_ligand(ligand_name, conc)
            self.er_beta.bind_ligand(ligand_name, conc)
    
    def update_crosstalk(self, dt: float):
        """implement pathway crosstalk"""
        # kinases phosphorylate er
        mapk_phospho = CrossTalk.mapk_to_er_phosphorylation(
            self.non_genomic.mapk_active)
        akt_phospho = CrossTalk.akt_to_er_phosphorylation(
            self.non_genomic.akt_active)
        src_phospho = CrossTalk.src_to_er_phosphorylation(
            self.non_genomic.src_active)
        
        self.er_alpha.phosphorylate('ser118', mapk_phospho)
        self.er_alpha.phosphorylate('ser167', akt_phospho)
        self.er_alpha.phosphorylate('tyr537', src_phospho)
        
        # calcium signaling
        calmodulin_active = CrossTalk.calcium_to_calmodulin(
            self.non_genomic.calcium_cytosol)
        
        # feedback from target genes
        cyclin_d1 = self.genomic.target_genes['cyclin_d1']['protein']
        cell_cycle_signal = CrossTalk.cyclin_d1_feedback(cyclin_d1)
        
        bcl2 = self.genomic.target_genes['bcl2']['protein']
        survival_signal = CrossTalk.bcl2_survival(bcl2)
        
        # update cellular state
        self.state.cell_cycle_progress += cell_cycle_signal * dt * 0.01
        self.state.survival_signal = survival_signal
        self.state.apoptosis_index = max(0, 1.0 - survival_signal)
    
    def update_cellular_responses(self, dt: float):
        """update cell behavior based on er signaling"""
        # proliferation driven by target genes
        cyclin_d1 = self.genomic.target_genes['cyclin_d1']['protein']
        c_myc = self.genomic.target_genes['c_myc']['protein']
        
        self.state.proliferation_rate = (cyclin_d1 * 0.3 + c_myc * 0.3) * dt
        
        # cell cycle progression
        if self.state.cell_cycle_progress >= 1.0:
            self.state.division_count += 1
            self.state.cell_cycle_progress = 0.0
        
        # metabolism
        vegf = self.genomic.target_genes['vegf']['protein']
        self.state.glucose_uptake = 1.0 + vegf * 0.5
    
    def simulate_step(self, dt: float) -> Dict:
        """
        single simulation timestep
        dt in minutes
        """
        # update ligand binding and pharmacokinetics
        self.update_ligand_binding(dt)
        
        # receptor dynamics
        self.er_alpha.dimerize(dt)
        self.er_beta.dimerize(dt)
        
        ere_availability = 0.8
        self.er_alpha.bind_dna(dt, ere_availability)
        self.er_beta.bind_dna(dt, ere_availability)
        
        # calculate combined er activity
        er_alpha_activity = self.er_alpha.calculate_transcriptional_activity()
        er_beta_activity = self.er_beta.calculate_transcriptional_activity()
        total_er_activity = er_alpha_activity + er_beta_activity * 0.7
        
        # membrane er activity (fraction at membrane)
        membrane_er_fraction = 0.05
        membrane_er_activity = (self.er_alpha.concentration + 
                               self.er_beta.concentration) * membrane_er_fraction * 0.1
        
        # update pathways
        genomic_state = self.genomic.update(total_er_activity, dt)
        non_genomic_state = self.non_genomic.update(membrane_er_activity, dt)
        
        # crosstalk
        self.update_crosstalk(dt)
        
        # cellular responses
        self.update_cellular_responses(dt)
        
        # receptor state updates
        er_alpha_state = self.er_alpha.update(dt)
        er_beta_state = self.er_beta.update(dt)
        
        # increment time
        self.time += dt
        
        return {
            'time': self.time,
            'er_alpha': er_alpha_state,
            'er_beta': er_beta_state,
            'genomic': genomic_state,
            'non_genomic': non_genomic_state,
            'cellular_state': {
                'proliferation_rate': self.state.proliferation_rate,
                'apoptosis_index': self.state.apoptosis_index,
                'survival_signal': self.state.survival_signal,
                'division_count': self.state.division_count,
                'cell_cycle_progress': self.state.cell_cycle_progress,
            }
        }
