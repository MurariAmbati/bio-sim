"""
core androgen receptor pathway model
"""

import numpy as np
from typing import Dict, List, Optional, Tuple
from dataclasses import dataclass

from .parameters import ParameterSet, DrugParameters, MutationParameters, CancerParameters


@dataclass
class PathwayState:
    """represents the state of all molecular species in the pathway"""
    
    # free androgen receptor
    ar_free_cytoplasm: float = 0.0
    ar_free_nucleus: float = 0.0
    
    # ligands (free, intracellular)
    testosterone: float = 0.0
    dht: float = 0.0
    
    # ligand-bound ar (monomers)
    ar_testosterone_cytoplasm: float = 0.0
    ar_testosterone_nucleus: float = 0.0
    ar_dht_cytoplasm: float = 0.0
    ar_dht_nucleus: float = 0.0
    
    # ar dimers
    ar_testosterone_dimer_nucleus: float = 0.0
    ar_dht_dimer_nucleus: float = 0.0
    
    # dna-bound complexes
    ar_testosterone_dna: float = 0.0
    ar_dht_dna: float = 0.0
    
    # coregulators
    coactivator_free: float = 0.0
    corepressor_free: float = 0.0
    ar_dna_coactivator: float = 0.0
    ar_dna_corepressor: float = 0.0
    
    # gene expression
    mrna_psa: float = 0.0
    protein_psa: float = 0.0
    mrna_klk2: float = 0.0
    mrna_tmprss2: float = 0.0
    mrna_nkx31: float = 0.0
    mrna_fkbp5: float = 0.0
    
    # available dna binding sites
    are_sites_free: float = 0.0
    
    def to_array(self) -> np.ndarray:
        """convert state to numpy array"""
        return np.array([
            self.ar_free_cytoplasm,
            self.ar_free_nucleus,
            self.testosterone,
            self.dht,
            self.ar_testosterone_cytoplasm,
            self.ar_testosterone_nucleus,
            self.ar_dht_cytoplasm,
            self.ar_dht_nucleus,
            self.ar_testosterone_dimer_nucleus,
            self.ar_dht_dimer_nucleus,
            self.ar_testosterone_dna,
            self.ar_dht_dna,
            self.coactivator_free,
            self.corepressor_free,
            self.ar_dna_coactivator,
            self.ar_dna_corepressor,
            self.mrna_psa,
            self.protein_psa,
            self.mrna_klk2,
            self.mrna_tmprss2,
            self.mrna_nkx31,
            self.mrna_fkbp5,
            self.are_sites_free
        ])
    
    @classmethod
    def from_array(cls, arr: np.ndarray) -> 'PathwayState':
        """create state from numpy array"""
        return cls(
            ar_free_cytoplasm=arr[0],
            ar_free_nucleus=arr[1],
            testosterone=arr[2],
            dht=arr[3],
            ar_testosterone_cytoplasm=arr[4],
            ar_testosterone_nucleus=arr[5],
            ar_dht_cytoplasm=arr[6],
            ar_dht_nucleus=arr[7],
            ar_testosterone_dimer_nucleus=arr[8],
            ar_dht_dimer_nucleus=arr[9],
            ar_testosterone_dna=arr[10],
            ar_dht_dna=arr[11],
            coactivator_free=arr[12],
            corepressor_free=arr[13],
            ar_dna_coactivator=arr[14],
            ar_dna_corepressor=arr[15],
            mrna_psa=arr[16],
            protein_psa=arr[17],
            mrna_klk2=arr[18],
            mrna_tmprss2=arr[19],
            mrna_nkx31=arr[20],
            mrna_fkbp5=arr[21],
            are_sites_free=arr[22]
        )
    
    @property
    def species_names(self) -> List[str]:
        """get list of species names"""
        return [
            'ar_free_cytoplasm', 'ar_free_nucleus',
            'testosterone', 'dht',
            'ar_testosterone_cytoplasm', 'ar_testosterone_nucleus',
            'ar_dht_cytoplasm', 'ar_dht_nucleus',
            'ar_testosterone_dimer_nucleus', 'ar_dht_dimer_nucleus',
            'ar_testosterone_dna', 'ar_dht_dna',
            'coactivator_free', 'corepressor_free',
            'ar_dna_coactivator', 'ar_dna_corepressor',
            'mrna_psa', 'protein_psa',
            'mrna_klk2', 'mrna_tmprss2', 'mrna_nkx31', 'mrna_fkbp5',
            'are_sites_free'
        ]


class ArPathwayModel:
    """
    comprehensive androgen receptor signaling pathway model
    
    implements the full reaction network for ar signaling including:
    - ligand binding (testosterone, dht)
    - nuclear translocation
    - dimerization
    - dna binding
    - coregulator recruitment
    - transcriptional activation
    - feedback regulation
    """
    
    def __init__(self, 
                 parameters: Optional[ParameterSet] = None,
                 drug: Optional[DrugParameters] = None,
                 mutation: Optional[MutationParameters] = None,
                 cancer: Optional[CancerParameters] = None):
        """
        initialize ar pathway model
        
        args:
            parameters: parameter set for kinetic constants
            drug: drug/antagonist parameters
            mutation: ar mutation parameters
            cancer: cancer phenotype parameters
        """
        self.params = parameters or ParameterSet()
        self.drug = drug or DrugParameters()
        self.mutation = mutation or MutationParameters()
        self.cancer = cancer or CancerParameters()
        
        # apply cancer-related parameter modifications
        self._apply_cancer_modifications()
        
        # apply mutation-related modifications
        self._apply_mutation_modifications()
    
    def _apply_cancer_modifications(self):
        """modify parameters based on cancer phenotype"""
        if self.cancer.phenotype != "normal":
            self.params.ar_total *= self.cancer.ar_expression_fold
            self.params.coactivator_total *= self.cancer.coactivator_upregulation
    
    def _apply_mutation_modifications(self):
        """modify parameters based on ar mutations"""
        if self.mutation.is_mutant:
            if self.mutation.stability_increase > 1.0:
                self.params.k_ar_degradation /= self.mutation.stability_increase
            if self.mutation.transcription_enhancement > 1.0:
                self.params.k_transcription_activated *= self.mutation.transcription_enhancement
            if self.mutation.nuclear_import_enhancement > 1.0:
                self.params.k_nuclear_import *= self.mutation.nuclear_import_enhancement
    
    def get_initial_state(self) -> PathwayState:
        """
        get initial state with reasonable starting conditions
        
        returns:
            initial pathway state
        """
        state = PathwayState()
        
        # initialize ar in cytoplasm and nucleus
        state.ar_free_cytoplasm = self.params.ar_cytoplasm_init
        state.ar_free_nucleus = self.params.ar_nucleus_init
        
        # initialize ligands
        state.testosterone = self.params.testosterone_external * 1e9  # convert to molecules
        state.dht = self.params.dht_external * 1e9
        
        # apply drug effects on testosterone levels
        if self.drug.blocks_testosterone:
            state.testosterone *= (1 - self.drug.testosterone_reduction)
        
        # add intracrine androgen in cancer
        if self.cancer.intracrine_androgen > 0:
            state.testosterone += self.cancer.intracrine_androgen * 1e9
            state.dht += self.cancer.intracrine_androgen * 1e9 * 0.1  # some conversion
        
        # initialize coregulators
        state.coactivator_free = self.params.coactivator_total
        state.corepressor_free = self.params.corepressor_total
        
        # initialize available dna sites
        state.are_sites_free = self.params.num_are_sites
        
        return state
    
    def compute_derivatives(self, state: PathwayState, t: float) -> PathwayState:
        """
        compute time derivatives for all species (ode system)
        
        args:
            state: current pathway state
            t: current time
            
        returns:
            derivatives (dstate/dt)
        """
        d = PathwayState()  # derivatives
        
        # apply drug effects
        f_binding = 1.0 - self.drug.binding_inhibition
        f_translocation = 1.0 - self.drug.translocation_inhibition
        f_dna_binding = 1.0 - self.drug.dna_binding_inhibition
        f_transcription = 1.0 - self.drug.transcription_inhibition
        f_5ar = 1.0 - self.drug.k_5ar_inhibition if self.drug.blocks_5ar else 1.0
        
        # degradation enhancement for protac
        k_deg_ar = self.params.k_ar_degradation
        if self.drug.is_degrader:
            k_deg_ar *= self.drug.k_degradation_enhancement
        
        # --- ligand dynamics ---
        
        # testosterone to dht conversion (5-alpha reductase)
        r_5ar = self.params.k_5alpha_reductase * f_5ar * state.testosterone
        
        d.testosterone = -r_5ar
        d.dht = r_5ar
        
        # --- ligand binding in cytoplasm ---
        
        # testosterone binding
        r_t_bind_cyto = (self.params.k_on_testosterone * f_binding * 
                         state.ar_free_cytoplasm * state.testosterone)
        r_t_unbind_cyto = self.params.k_off_testosterone * state.ar_testosterone_cytoplasm
        
        # dht binding
        r_dht_bind_cyto = (self.params.k_on_dht * f_binding * 
                          state.ar_free_cytoplasm * state.dht)
        r_dht_unbind_cyto = self.params.k_off_dht * state.ar_dht_cytoplasm
        
        d.ar_free_cytoplasm = (r_t_unbind_cyto + r_dht_unbind_cyto - 
                               r_t_bind_cyto - r_dht_bind_cyto -
                               k_deg_ar * state.ar_free_cytoplasm)
        
        d.testosterone += (r_t_unbind_cyto - r_t_bind_cyto)
        d.dht += (r_dht_unbind_cyto - r_dht_bind_cyto)
        
        # --- nuclear translocation ---
        
        # testosterone-bound ar
        r_t_import = self.params.k_nuclear_import * f_translocation * state.ar_testosterone_cytoplasm
        r_t_export = self.params.k_nuclear_export * state.ar_testosterone_nucleus
        
        d.ar_testosterone_cytoplasm = (r_t_bind_cyto - r_t_unbind_cyto - 
                                       r_t_import + r_t_export -
                                       k_deg_ar * state.ar_testosterone_cytoplasm)
        
        # dht-bound ar
        r_dht_import = self.params.k_nuclear_import * f_translocation * state.ar_dht_cytoplasm
        r_dht_export = self.params.k_nuclear_export * state.ar_dht_nucleus
        
        d.ar_dht_cytoplasm = (r_dht_bind_cyto - r_dht_unbind_cyto - 
                              r_dht_import + r_dht_export -
                              k_deg_ar * state.ar_dht_cytoplasm)
        
        # --- nuclear ar ---
        
        # ligand binding in nucleus
        r_t_bind_nuc = (self.params.k_on_testosterone * f_binding * 
                       state.ar_free_nucleus * state.testosterone)
        r_t_unbind_nuc = self.params.k_off_testosterone * state.ar_testosterone_nucleus
        
        r_dht_bind_nuc = (self.params.k_on_dht * f_binding * 
                         state.ar_free_nucleus * state.dht)
        r_dht_unbind_nuc = self.params.k_off_dht * state.ar_dht_nucleus
        
        d.ar_free_nucleus = (r_t_unbind_nuc + r_dht_unbind_nuc - 
                            r_t_bind_nuc - r_dht_bind_nuc -
                            k_deg_ar * state.ar_free_nucleus)
        
        # --- dimerization ---
        
        # testosterone-bound dimer
        r_t_dimer = (self.params.k_dimerization * 
                     state.ar_testosterone_nucleus * state.ar_testosterone_nucleus)
        r_t_dimer_dissoc = self.params.k_dimer_dissociation * state.ar_testosterone_dimer_nucleus
        
        # dht-bound dimer
        r_dht_dimer = (self.params.k_dimerization * 
                       state.ar_dht_nucleus * state.ar_dht_nucleus)
        r_dht_dimer_dissoc = self.params.k_dimer_dissociation * state.ar_dht_dimer_nucleus
        
        d.ar_testosterone_nucleus = (r_t_bind_nuc - r_t_unbind_nuc + 
                                     r_t_import - r_t_export -
                                     2 * r_t_dimer + 2 * r_t_dimer_dissoc -
                                     k_deg_ar * state.ar_testosterone_nucleus)
        
        d.ar_dht_nucleus = (r_dht_bind_nuc - r_dht_unbind_nuc + 
                           r_dht_import - r_dht_export -
                           2 * r_dht_dimer + 2 * r_dht_dimer_dissoc -
                           k_deg_ar * state.ar_dht_nucleus)
        
        d.ar_testosterone_dimer_nucleus = (r_t_dimer - r_t_dimer_dissoc -
                                           k_deg_ar * state.ar_testosterone_dimer_nucleus)
        
        d.ar_dht_dimer_nucleus = (r_dht_dimer - r_dht_dimer_dissoc -
                                  k_deg_ar * state.ar_dht_dimer_nucleus)
        
        # --- dna binding ---
        
        # testosterone dimer to dna
        r_t_dna_bind = (self.params.k_dna_on * f_dna_binding * 
                       state.ar_testosterone_dimer_nucleus * state.are_sites_free)
        r_t_dna_unbind = self.params.k_dna_off * state.ar_testosterone_dna
        
        # dht dimer to dna
        r_dht_dna_bind = (self.params.k_dna_on * f_dna_binding * 
                         state.ar_dht_dimer_nucleus * state.are_sites_free)
        r_dht_dna_unbind = self.params.k_dna_off * state.ar_dht_dna
        
        d.ar_testosterone_dimer_nucleus += (-r_t_dna_bind + r_t_dna_unbind)
        d.ar_dht_dimer_nucleus += (-r_dht_dna_bind + r_dht_dna_unbind)
        
        d.ar_testosterone_dna = (r_t_dna_bind - r_t_dna_unbind -
                                k_deg_ar * state.ar_testosterone_dna)
        
        d.ar_dht_dna = (r_dht_dna_bind - r_dht_dna_unbind -
                       k_deg_ar * state.ar_dht_dna)
        
        d.are_sites_free = (-r_t_dna_bind - r_dht_dna_bind + 
                            r_t_dna_unbind + r_dht_dna_unbind)
        
        # --- coregulator dynamics ---
        
        total_ar_dna = state.ar_testosterone_dna + state.ar_dht_dna
        
        # coactivator recruitment
        r_coact_on = self.params.k_coactivator_on * total_ar_dna * state.coactivator_free
        r_coact_off = self.params.k_coactivator_off * state.ar_dna_coactivator
        
        d.coactivator_free = -r_coact_on + r_coact_off
        d.ar_dna_coactivator = r_coact_on - r_coact_off
        
        # corepressor recruitment (competes with coactivator)
        r_corep_on = self.params.k_corepressor_on * total_ar_dna * state.corepressor_free
        r_corep_off = self.params.k_corepressor_off * state.ar_dna_corepressor
        
        d.corepressor_free = -r_corep_on + r_corep_off
        d.ar_dna_corepressor = r_corep_on - r_corep_off
        
        # --- transcription ---
        
        # basal transcription
        k_tx_basal = self.params.k_transcription_basal
        
        # activated transcription (with coactivator)
        k_tx_activated = self.params.k_transcription_activated * f_transcription
        
        # feedback from fkbp5
        feedback_fkbp5 = 1.0 / (1.0 + self.params.k_feedback_fkbp5 * state.mrna_fkbp5)
        
        # psa transcription
        r_tx_psa = (k_tx_basal * state.are_sites_free + 
                   k_tx_activated * state.ar_dna_coactivator * feedback_fkbp5)
        
        d.mrna_psa = r_tx_psa - self.params.k_mrna_degradation * state.mrna_psa
        
        # klk2 transcription
        r_tx_klk2 = (k_tx_basal * 0.5 * state.are_sites_free + 
                    k_tx_activated * 0.8 * state.ar_dna_coactivator)
        
        d.mrna_klk2 = r_tx_klk2 - self.params.k_mrna_degradation * state.mrna_klk2
        
        # tmprss2 transcription
        r_tx_tmprss2 = (k_tx_basal * 0.3 * state.are_sites_free + 
                       k_tx_activated * 0.9 * state.ar_dna_coactivator)
        
        d.mrna_tmprss2 = r_tx_tmprss2 - self.params.k_mrna_degradation * state.mrna_tmprss2
        
        # nkx3.1 transcription
        r_tx_nkx31 = (k_tx_basal * 0.4 * state.are_sites_free + 
                     k_tx_activated * 0.7 * state.ar_dna_coactivator)
        
        d.mrna_nkx31 = r_tx_nkx31 - self.params.k_mrna_degradation * state.mrna_nkx31
        
        # fkbp5 transcription (negative feedback)
        r_tx_fkbp5 = (k_tx_basal * 0.2 * state.are_sites_free + 
                     k_tx_activated * 1.2 * state.ar_dna_coactivator)
        
        d.mrna_fkbp5 = r_tx_fkbp5 - self.params.k_mrna_degradation * state.mrna_fkbp5
        
        # --- translation (psa protein as example) ---
        
        r_tl_psa = self.params.k_translation * state.mrna_psa
        
        d.protein_psa = r_tl_psa - self.params.k_protein_degradation * state.protein_psa
        
        return d
    
    def get_reaction_list(self) -> List[Dict]:
        """
        get list of all reactions in the model
        
        returns:
            list of reaction dictionaries with reactants, products, rate
        """
        reactions = [
            {
                'name': 'testosterone_to_dht',
                'reactants': ['testosterone'],
                'products': ['dht'],
                'rate': 'k_5alpha_reductase',
                'type': 'conversion'
            },
            {
                'name': 'ar_testosterone_binding_cytoplasm',
                'reactants': ['ar_free_cytoplasm', 'testosterone'],
                'products': ['ar_testosterone_cytoplasm'],
                'rate': 'k_on_testosterone',
                'type': 'binding'
            },
            {
                'name': 'ar_dht_binding_cytoplasm',
                'reactants': ['ar_free_cytoplasm', 'dht'],
                'products': ['ar_dht_cytoplasm'],
                'rate': 'k_on_dht',
                'type': 'binding'
            },
            {
                'name': 'ar_testosterone_nuclear_import',
                'reactants': ['ar_testosterone_cytoplasm'],
                'products': ['ar_testosterone_nucleus'],
                'rate': 'k_nuclear_import',
                'type': 'transport'
            },
            {
                'name': 'ar_dht_nuclear_import',
                'reactants': ['ar_dht_cytoplasm'],
                'products': ['ar_dht_nucleus'],
                'rate': 'k_nuclear_import',
                'type': 'transport'
            },
            {
                'name': 'ar_testosterone_dimerization',
                'reactants': ['ar_testosterone_nucleus', 'ar_testosterone_nucleus'],
                'products': ['ar_testosterone_dimer_nucleus'],
                'rate': 'k_dimerization',
                'type': 'dimerization'
            },
            {
                'name': 'ar_dht_dimerization',
                'reactants': ['ar_dht_nucleus', 'ar_dht_nucleus'],
                'products': ['ar_dht_dimer_nucleus'],
                'rate': 'k_dimerization',
                'type': 'dimerization'
            },
            {
                'name': 'ar_testosterone_dna_binding',
                'reactants': ['ar_testosterone_dimer_nucleus', 'are_sites_free'],
                'products': ['ar_testosterone_dna'],
                'rate': 'k_dna_on',
                'type': 'dna_binding'
            },
            {
                'name': 'ar_dht_dna_binding',
                'reactants': ['ar_dht_dimer_nucleus', 'are_sites_free'],
                'products': ['ar_dht_dna'],
                'rate': 'k_dna_on',
                'type': 'dna_binding'
            },
            {
                'name': 'psa_transcription',
                'reactants': ['ar_testosterone_dna'],
                'products': ['ar_testosterone_dna', 'mrna_psa'],
                'rate': 'k_transcription_activated',
                'type': 'transcription'
            },
            {
                'name': 'psa_translation',
                'reactants': ['mrna_psa'],
                'products': ['mrna_psa', 'protein_psa'],
                'rate': 'k_translation',
                'type': 'translation'
            }
        ]
        
        return reactions
    
    def get_conserved_quantities(self, state: PathwayState) -> Dict[str, float]:
        """
        calculate conserved quantities (useful for validation)
        
        args:
            state: current pathway state
            
        returns:
            dictionary of conserved quantities
        """
        # total ar (all forms)
        total_ar = (state.ar_free_cytoplasm + state.ar_free_nucleus +
                   state.ar_testosterone_cytoplasm + state.ar_testosterone_nucleus +
                   state.ar_dht_cytoplasm + state.ar_dht_nucleus +
                   2 * state.ar_testosterone_dimer_nucleus +
                   2 * state.ar_dht_dimer_nucleus +
                   2 * state.ar_testosterone_dna +
                   2 * state.ar_dht_dna)
        
        # total coactivator
        total_coactivator = state.coactivator_free + state.ar_dna_coactivator
        
        # total corepressor
        total_corepressor = state.corepressor_free + state.ar_dna_corepressor
        
        # total are sites
        total_are = (state.are_sites_free + state.ar_testosterone_dna + 
                    state.ar_dht_dna)
        
        return {
            'total_ar': total_ar,
            'total_coactivator': total_coactivator,
            'total_corepressor': total_corepressor,
            'total_are_sites': total_are
        }
