import numpy as np
from scipy.integrate import odeint
from dataclasses import dataclass
from typing import Tuple, Dict, List


@dataclass
class ThyroidReceptorParameters:
    """parameters for thyroid hormone receptor signaling model"""
    
    # hormone synthesis and secretion
    k_syn_t4: float = 0.8           # t4 synthesis rate
    k_syn_t3: float = 0.2           # t3 synthesis rate
    k_sec: float = 0.5              # secretion rate
    
    # deiodinase activity
    k_d1: float = 0.4               # type 1 deiodinase (t4->t3)
    k_d2: float = 0.6               # type 2 deiodinase (t4->t3)
    k_d3: float = 0.3               # type 3 deiodinase (inactivation)
    
    # receptor dynamics
    k_bind_t3: float = 1.2          # t3-tr binding rate
    k_unbind: float = 0.2           # unbinding rate
    k_syn_tr: float = 0.5           # tr synthesis
    k_deg_tr: float = 0.1           # tr degradation
    
    # coregulator recruitment
    k_coreg_on: float = 0.8         # coactivator binding
    k_coreg_off: float = 0.3        # coactivator unbinding
    k_corep_on: float = 0.5         # corepressor binding
    k_corep_off: float = 0.4        # corepressor unbinding
    
    # gene transcription
    k_basal: float = 0.1            # basal transcription
    k_trans_act: float = 2.0        # activated transcription
    k_trans_rep: float = 0.05       # repressed transcription
    
    # metabolic effects
    k_met_act: float = 1.5          # metabolic activation
    k_met_deg: float = 0.3          # metabolic enzyme degradation
    
    # developmental effects
    k_dev_act: float = 0.8          # developmental program activation
    k_dev_prog: float = 0.2         # developmental progression
    
    # feedback regulation
    k_fb_tsh: float = 0.7           # tsh feedback strength
    k_fb_hormone: float = 0.5       # hormone feedback
    
    # degradation rates
    k_deg_t4: float = 0.05          # t4 degradation
    k_deg_t3: float = 0.15          # t3 degradation
    k_deg_mrna: float = 0.4         # mrna degradation
    k_deg_protein: float = 0.2      # protein degradation


class ThyroidReceptorModel:
    """comprehensive thyroid hormone receptor signaling model
    
    models:
    - thyroid hormone synthesis and metabolism
    - receptor binding and activation
    - coregulator dynamics
    - target gene transcription
    - metabolic and developmental effects
    - feedback regulation
    """
    
    def __init__(self, params: ThyroidReceptorParameters = None):
        self.params = params or ThyroidReceptorParameters()
        
    def derivatives(self, state: np.ndarray, t: float, 
                   tsh_signal: float = 1.0, stress: float = 0.0) -> np.ndarray:
        """calculate time derivatives for all state variables
        
        state variables:
        0: tsh (thyroid stimulating hormone)
        1: t4 (thyroxine)
        2: t3 (triiodothyronine)
        3: rt3 (reverse t3)
        4: tr_free (free thyroid receptor)
        5: tr_t3 (t3-bound receptor)
        6: tr_t3_coreg (receptor-coactivator complex)
        7: tr_corep (receptor-corepressor complex)
        8: target_mrna (target gene mrna)
        9: metabolic_enzyme (metabolic proteins)
        10: developmental_factor (developmental proteins)
        11: metabolic_rate (cellular metabolism)
        12: d2_enzyme (type 2 deiodinase)
        13: d3_enzyme (type 3 deiodinase)
        """
        
        p = self.params
        
        # unpack state
        tsh, t4, t3, rt3, tr_free, tr_t3, tr_t3_coreg, tr_corep = state[:8]
        target_mrna, met_enz, dev_factor, met_rate = state[8:12]
        d2, d3 = state[12:14]
        
        # feedback regulation
        hormone_feedback = 1.0 / (1.0 + p.k_fb_hormone * (t4 + 2.0 * t3))
        tsh_production = tsh_signal * hormone_feedback
        
        # thyroid hormone synthesis (tsh-dependent)
        t4_synthesis = p.k_syn_t4 * tsh * p.k_sec
        t3_synthesis = p.k_syn_t3 * tsh * p.k_sec
        
        # deiodinase-mediated conversion
        t4_to_t3 = (p.k_d1 + d2) * t4
        t4_to_rt3 = d3 * t4
        
        # receptor dynamics
        tr_binding = p.k_bind_t3 * t3 * tr_free
        tr_unbinding = p.k_unbind * tr_t3
        tr_synthesis = p.k_syn_tr
        tr_degradation = p.k_deg_tr * (tr_free + tr_t3)
        
        # coregulator recruitment
        coactivator_binding = p.k_coreg_on * tr_t3
        coactivator_release = p.k_coreg_off * tr_t3_coreg
        corepressor_binding = p.k_corep_on * tr_free
        corepressor_release = p.k_corep_off * tr_corep
        
        # transcriptional regulation
        basal_transcription = p.k_basal
        activated_transcription = p.k_trans_act * tr_t3_coreg
        repressed_transcription = -p.k_trans_rep * tr_corep
        total_transcription = (basal_transcription + activated_transcription + 
                              repressed_transcription)
        
        # metabolic effects
        metabolic_activation = p.k_met_act * met_enz
        
        # developmental effects
        developmental_activation = p.k_dev_act * tr_t3_coreg
        
        # deiodinase regulation (feedback)
        d2_regulation = p.k_d2 * (1.0 / (1.0 + 0.5 * t3))
        d3_regulation = p.k_d3 * (1.0 + 0.3 * t3)
        
        # derivatives
        d_state = np.zeros(14)
        
        # tsh dynamics
        d_state[0] = tsh_production - p.k_fb_tsh * tsh * (t4 + t3)
        
        # hormone levels
        d_state[1] = (t4_synthesis - t4_to_t3 - t4_to_rt3 - 
                     p.k_deg_t4 * t4)
        d_state[2] = (t3_synthesis + t4_to_t3 - tr_binding + tr_unbinding - 
                     p.k_deg_t3 * t3)
        d_state[3] = t4_to_rt3 - p.k_deg_t3 * rt3
        
        # receptor states
        d_state[4] = (tr_synthesis - tr_binding + tr_unbinding - 
                     corepressor_binding + corepressor_release - 
                     p.k_deg_tr * tr_free)
        d_state[5] = (tr_binding - tr_unbinding - coactivator_binding + 
                     coactivator_release - p.k_deg_tr * tr_t3)
        d_state[6] = coactivator_binding - coactivator_release
        d_state[7] = corepressor_binding - corepressor_release
        
        # target gene expression
        d_state[8] = total_transcription - p.k_deg_mrna * target_mrna
        d_state[9] = p.k_met_act * target_mrna - p.k_met_deg * met_enz
        d_state[10] = developmental_activation - p.k_deg_protein * dev_factor
        
        # metabolic rate
        d_state[11] = metabolic_activation - p.k_met_deg * met_rate
        
        # deiodinase levels
        d_state[12] = d2_regulation - 0.2 * d2
        d_state[13] = d3_regulation - 0.2 * d3
        
        return d_state
    
    def simulate(self, t_span: Tuple[float, float], 
                initial_state: np.ndarray = None,
                tsh_signal: float = 1.0,
                stress: float = 0.0,
                n_points: int = 1000) -> Tuple[np.ndarray, np.ndarray]:
        """run simulation
        
        returns:
            time points and state trajectories
        """
        
        if initial_state is None:
            # physiological steady state
            initial_state = np.array([
                1.0,    # tsh
                100.0,  # t4
                2.0,    # t3
                0.5,    # rt3
                5.0,    # tr_free
                0.5,    # tr_t3
                0.2,    # tr_t3_coreg
                1.0,    # tr_corep
                1.0,    # target_mrna
                1.0,    # metabolic_enzyme
                0.5,    # developmental_factor
                1.0,    # metabolic_rate
                0.6,    # d2
                0.3,    # d3
            ])
        
        t = np.linspace(t_span[0], t_span[1], n_points)
        solution = odeint(self.derivatives, initial_state, t, 
                         args=(tsh_signal, stress))
        
        return t, solution
    
    def get_state_labels(self) -> List[str]:
        """get labels for state variables"""
        return [
            'tsh', 't4', 't3', 'rt3',
            'tr_free', 'tr_t3', 'tr_t3_coreg', 'tr_corep',
            'target_mrna', 'metabolic_enzyme', 'developmental_factor',
            'metabolic_rate', 'd2', 'd3'
        ]
    
    def calculate_metrics(self, solution: np.ndarray) -> Dict[str, float]:
        """calculate physiological metrics from simulation"""
        
        metrics = {
            'mean_t3': np.mean(solution[:, 2]),
            'mean_t4': np.mean(solution[:, 1]),
            't4_t3_ratio': np.mean(solution[:, 1]) / np.mean(solution[:, 2]),
            'receptor_occupancy': np.mean(solution[:, 5]) / 
                                 (np.mean(solution[:, 4]) + np.mean(solution[:, 5])),
            'mean_metabolic_rate': np.mean(solution[:, 11]),
            'transcriptional_activity': np.mean(solution[:, 6]),
            'developmental_index': np.mean(solution[:, 10]),
            'mean_tsh': np.mean(solution[:, 0]),
        }
        
        return metrics


def simulate_hypothyroidism(model: ThyroidReceptorModel, 
                           severity: float = 0.5) -> Tuple[np.ndarray, np.ndarray]:
    """simulate hypothyroid state"""
    params = model.params
    params.k_syn_t4 *= (1.0 - severity)
    params.k_syn_t3 *= (1.0 - severity)
    
    return model.simulate((0, 200), tsh_signal=1.5)


def simulate_hyperthyroidism(model: ThyroidReceptorModel,
                            severity: float = 0.5) -> Tuple[np.ndarray, np.ndarray]:
    """simulate hyperthyroid state"""
    params = model.params
    params.k_syn_t4 *= (1.0 + severity)
    params.k_syn_t3 *= (1.0 + severity)
    
    return model.simulate((0, 200), tsh_signal=0.3)


def simulate_tr_mutation(model: ThyroidReceptorModel,
                        binding_defect: float = 0.7) -> Tuple[np.ndarray, np.ndarray]:
    """simulate thyroid receptor mutation"""
    params = model.params
    params.k_bind_t3 *= (1.0 - binding_defect)
    
    return model.simulate((0, 200))


def simulate_development(model: ThyroidReceptorModel) -> Tuple[np.ndarray, np.ndarray]:
    """simulate developmental trajectory with changing thyroid hormones"""
    
    # developmental stages with different hormone requirements
    t_early = np.linspace(0, 50, 300)
    t_mid = np.linspace(50, 100, 300)
    t_late = np.linspace(100, 150, 300)
    
    # early development - low hormones
    _, sol_early = model.simulate((0, 50), tsh_signal=0.5, n_points=300)
    
    # mid development - increasing hormones
    initial_mid = sol_early[-1, :]
    _, sol_mid = model.simulate((50, 100), initial_state=initial_mid,
                               tsh_signal=1.2, n_points=300)
    
    # late development - mature levels
    initial_late = sol_mid[-1, :]
    _, sol_late = model.simulate((100, 150), initial_state=initial_late,
                                tsh_signal=1.0, n_points=300)
    
    t_total = np.concatenate([t_early, t_mid, t_late])
    sol_total = np.vstack([sol_early, sol_mid, sol_late])
    
    return t_total, sol_total
