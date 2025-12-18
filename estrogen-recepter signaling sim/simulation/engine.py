"""
simulation engine for er signaling dynamics
"""

import numpy as np
import pandas as pd
from typing import Dict, List, Optional, Tuple
from models.cell import EstrogenResponsiveCell, CellType
from models.receptor import ReceptorParameters
from models.pathways import PathwayParameters


class ERSimulation:
    """
    main simulation engine for estrogen receptor signaling
    """
    
    def __init__(self, 
                 cell_type: str = CellType.BREAST_CANCER_MCF7,
                 duration_hours: float = 48.0,
                 dt_minutes: float = 1.0):
        
        self.cell_type = cell_type
        self.duration_hours = duration_hours
        self.dt_minutes = dt_minutes
        
        # initialize cell
        self.cell = EstrogenResponsiveCell(cell_type)
        
        # simulation results storage
        self.results = []
        self.time_points = []
        
    def add_treatment(self, ligand_name: str, dose_mg: float, 
                     administration_time_hours: float = 0.0):
        """
        add ligand treatment at specified time
        """
        self.cell.add_ligand(ligand_name, dose_mg)
        
    def run(self, verbose: bool = False) -> pd.DataFrame:
        """
        execute simulation
        returns dataframe of results
        """
        n_steps = int((self.duration_hours * 60) / self.dt_minutes)
        
        self.results = []
        self.time_points = []
        
        for step in range(n_steps):
            # simulate one timestep
            state = self.cell.simulate_step(self.dt_minutes)
            
            # flatten nested state for dataframe
            flat_state = self._flatten_state(state)
            self.results.append(flat_state)
            self.time_points.append(state['time'])
            
            if verbose and step % 100 == 0:
                print(f"step {step}/{n_steps}, time={state['time']:.1f} min")
        
        # convert to dataframe
        df = pd.DataFrame(self.results)
        df['time_hours'] = df['time'] / 60.0
        
        return df
    
    def _flatten_state(self, state: Dict) -> Dict:
        """flatten nested state dictionary"""
        flat = {'time': state['time']}
        
        # er receptors
        for receptor_type in ['er_alpha', 'er_beta']:
            if receptor_type in state:
                for key, val in state[receptor_type].items():
                    flat[f'{receptor_type}_{key}'] = val
        
        # genomic pathway (target genes)
        if 'genomic' in state:
            for gene, data in state['genomic'].items():
                flat[f'gene_{gene}_mrna'] = data['mrna']
                flat[f'gene_{gene}_protein'] = data['protein']
        
        # non-genomic pathway
        if 'non_genomic' in state:
            for key, val in state['non_genomic'].items():
                flat[f'nongenomic_{key}'] = val
        
        # cellular state
        if 'cellular_state' in state:
            for key, val in state['cellular_state'].items():
                flat[f'cell_{key}'] = val
        
        return flat
    
    def get_summary_statistics(self, df: pd.DataFrame) -> Dict:
        """calculate summary statistics from simulation"""
        stats = {}
        
        # receptor levels
        stats['er_alpha_mean'] = df['er_alpha_concentration'].mean()
        stats['er_alpha_max'] = df['er_alpha_concentration'].max()
        stats['er_beta_mean'] = df['er_beta_concentration'].mean()
        
        # transcriptional activity
        stats['transcription_auc'] = np.trapz(
            df['er_alpha_transcriptional_activity'], 
            df['time_hours']
        )
        
        # target gene expression
        gene_cols = [col for col in df.columns if 'gene_' in col and '_protein' in col]
        for col in gene_cols:
            gene_name = col.replace('gene_', '').replace('_protein', '')
            stats[f'{gene_name}_peak'] = df[col].max()
            stats[f'{gene_name}_auc'] = np.trapz(df[col], df['time_hours'])
        
        # kinase activation
        stats['mapk_peak'] = df['nongenomic_mapk_active'].max()
        stats['akt_peak'] = df['nongenomic_akt_active'].max()
        
        # cellular outcomes
        stats['final_divisions'] = df['cell_division_count'].iloc[-1]
        stats['proliferation_rate_mean'] = df['cell_proliferation_rate'].mean()
        stats['survival_signal_mean'] = df['cell_survival_signal'].mean()
        
        return stats


class DoseResponseSimulation:
    """
    simulate dose-response relationships
    """
    
    def __init__(self, cell_type: str = CellType.BREAST_CANCER_MCF7):
        self.cell_type = cell_type
        
    def run_dose_response(self, 
                         ligand_name: str,
                         doses: List[float],
                         duration_hours: float = 24.0) -> pd.DataFrame:
        """
        run simulations across dose range
        """
        results = []
        
        for dose in doses:
            sim = ERSimulation(self.cell_type, duration_hours)
            sim.add_treatment(ligand_name, dose)
            df = sim.run()
            
            # get summary for this dose
            stats = sim.get_summary_statistics(df)
            stats['dose'] = dose
            stats['ligand'] = ligand_name
            results.append(stats)
        
        return pd.DataFrame(results)
    
    def calculate_ec50(self, dose_response_df: pd.DataFrame, 
                      response_column: str) -> float:
        """
        calculate ec50 from dose-response data
        """
        doses = dose_response_df['dose'].values
        responses = dose_response_df[response_column].values
        
        # normalize to 0-1
        responses_norm = (responses - responses.min()) / (responses.max() - responses.min())
        
        # find dose where response = 0.5
        idx = np.argmin(np.abs(responses_norm - 0.5))
        ec50 = doses[idx]
        
        return ec50


class TimeCourseSimulation:
    """
    detailed time-course analysis
    """
    
    def __init__(self, cell_type: str = CellType.BREAST_CANCER_MCF7):
        self.cell_type = cell_type
    
    def run_timecourse(self,
                      ligand_name: str,
                      dose: float,
                      duration_hours: float = 48.0,
                      sampling_interval_minutes: float = 5.0) -> pd.DataFrame:
        """
        run simulation with specified sampling
        """
        sim = ERSimulation(self.cell_type, duration_hours, dt_minutes=1.0)
        sim.add_treatment(ligand_name, dose)
        df = sim.run()
        
        # downsample to sampling interval
        sampling_points = np.arange(0, duration_hours, sampling_interval_minutes / 60.0)
        df_sampled = pd.DataFrame()
        
        for t in sampling_points:
            idx = (df['time_hours'] - t).abs().argmin()
            df_sampled = pd.concat([df_sampled, df.iloc[[idx]]], ignore_index=True)
        
        return df_sampled
    
    def compare_ligands(self,
                       ligands: List[str],
                       dose: float,
                       duration_hours: float = 24.0) -> Dict[str, pd.DataFrame]:
        """
        compare multiple ligands at same dose
        """
        results = {}
        
        for ligand in ligands:
            df = self.run_timecourse(ligand, dose, duration_hours)
            results[ligand] = df
        
        return results


class TissueSpecificSimulation:
    """
    simulate er signaling across different tissue types
    """
    
    def run_tissue_comparison(self,
                            ligand_name: str,
                            dose: float,
                            tissues: Optional[List[str]] = None) -> Dict[str, pd.DataFrame]:
        """
        compare ligand effects across tissues
        """
        if tissues is None:
            tissues = [
                CellType.BREAST_CANCER_MCF7,
                CellType.OSTEOBLAST,
                CellType.ENDOMETRIAL,
                CellType.ENDOTHELIAL,
            ]
        
        results = {}
        
        for tissue in tissues:
            sim = ERSimulation(tissue, duration_hours=24.0)
            sim.add_treatment(ligand_name, dose)
            df = sim.run()
            results[tissue] = df
        
        return results
    
    def calculate_selectivity_index(self, 
                                   tissue_results: Dict[str, pd.DataFrame],
                                   response_column: str,
                                   reference_tissue: str = CellType.BREAST_CANCER_MCF7) -> Dict:
        """
        calculate tissue selectivity indices
        """
        reference_response = tissue_results[reference_tissue][response_column].max()
        
        selectivity = {}
        for tissue, df in tissue_results.items():
            if tissue != reference_tissue:
                tissue_response = df[response_column].max()
                selectivity[tissue] = tissue_response / reference_response
        
        return selectivity
