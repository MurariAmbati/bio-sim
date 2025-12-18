"""
simulation engine for ar pathway model
"""

import numpy as np
from scipy.integrate import odeint, solve_ivp
from typing import Dict, List, Optional, Tuple
import pandas as pd

from ..models.ar_pathway import ArPathwayModel, PathwayState
from ..models.parameters import ParameterSet


class SimulationResult:
    """container for simulation results"""
    
    def __init__(self, 
                 time: np.ndarray,
                 states: np.ndarray,
                 species_names: List[str],
                 parameters: ParameterSet):
        """
        initialize simulation result
        
        args:
            time: time points
            states: state values at each time point (shape: [n_times, n_species])
            species_names: names of species
            parameters: parameters used in simulation
        """
        self.time = time
        self.states = states
        self.species_names = species_names
        self.parameters = parameters
    
    def get_species(self, species_name: str) -> np.ndarray:
        """get time series for specific species"""
        idx = self.species_names.index(species_name)
        return self.states[:, idx]
    
    def to_dataframe(self) -> pd.DataFrame:
        """convert results to pandas dataframe"""
        df = pd.DataFrame(self.states, columns=self.species_names)
        df.insert(0, 'time', self.time)
        return df
    
    def get_steady_state(self, threshold: float = 1e-6) -> Optional[PathwayState]:
        """
        check if system reached steady state and return it
        
        args:
            threshold: threshold for considering system at steady state
            
        returns:
            steady state if reached, none otherwise
        """
        # check last 10% of simulation
        n_check = max(10, int(0.1 * len(self.time)))
        recent_states = self.states[-n_check:]
        
        # calculate coefficient of variation for each species
        cv = np.std(recent_states, axis=0) / (np.mean(recent_states, axis=0) + 1e-10)
        
        if np.all(cv < threshold):
            return PathwayState.from_array(self.states[-1])
        return None


class Simulator:
    """
    main simulation engine for ar pathway
    
    supports both ode (deterministic) and stochastic simulations
    """
    
    def __init__(self, 
                 model: ArPathwayModel,
                 time_end: float = 10000.0,
                 dt: float = 1.0,
                 method: str = 'ode'):
        """
        initialize simulator
        
        args:
            model: ar pathway model
            time_end: simulation end time (seconds)
            dt: time step for output (seconds)
            method: 'ode' for deterministic, 'gillespie' for stochastic
        """
        self.model = model
        self.time_end = time_end
        self.dt = dt
        self.method = method
    
    def run(self, 
            initial_state: Optional[PathwayState] = None,
            **kwargs) -> SimulationResult:
        """
        run simulation
        
        args:
            initial_state: initial state (if none, use model default)
            **kwargs: additional parameters to update before simulation
            
        returns:
            simulation results
        """
        # update model parameters if provided
        if kwargs:
            self.model.params.update(**kwargs)
        
        # get initial state
        if initial_state is None:
            initial_state = self.model.get_initial_state()
        
        # run appropriate simulation method
        if self.method == 'ode':
            return self._run_ode(initial_state)
        elif self.method == 'gillespie':
            return self._run_gillespie(initial_state)
        else:
            raise ValueError(f"unknown method: {self.method}")
    
    def _run_ode(self, initial_state: PathwayState) -> SimulationResult:
        """
        run deterministic ode simulation
        
        args:
            initial_state: initial pathway state
            
        returns:
            simulation results
        """
        # time points
        t_span = (0, self.time_end)
        t_eval = np.arange(0, self.time_end + self.dt, self.dt)
        
        # initial conditions
        y0 = initial_state.to_array()
        
        # define ode function
        def ode_func(t, y):
            state = PathwayState.from_array(y)
            derivs = self.model.compute_derivatives(state, t)
            return derivs.to_array()
        
        # solve ode system
        solution = solve_ivp(
            ode_func,
            t_span,
            y0,
            method='LSODA',
            t_eval=t_eval,
            rtol=1e-6,
            atol=1e-9
        )
        
        return SimulationResult(
            time=solution.t,
            states=solution.y.T,
            species_names=initial_state.species_names,
            parameters=self.model.params
        )
    
    def _run_gillespie(self, initial_state: PathwayState) -> SimulationResult:
        """
        run stochastic gillespie simulation
        
        args:
            initial_state: initial pathway state
            
        returns:
            simulation results
        """
        # initialize
        current_time = 0.0
        current_state = initial_state.to_array().copy()
        
        # storage for trajectory
        times = [0.0]
        states = [current_state.copy()]
        
        # output time points
        next_output_time = self.dt
        
        while current_time < self.time_end:
            # compute propensities for all reactions
            state_obj = PathwayState.from_array(current_state)
            propensities = self._compute_propensities(state_obj)
            
            # total propensity
            a0 = np.sum(propensities)
            
            if a0 == 0:
                # no more reactions possible
                break
            
            # time to next reaction
            tau = np.random.exponential(1.0 / a0)
            current_time += tau
            
            # which reaction fires?
            reaction_idx = np.random.choice(len(propensities), p=propensities / a0)
            
            # update state
            current_state = self._apply_reaction(current_state, reaction_idx)
            
            # store at output times
            while next_output_time <= current_time and next_output_time <= self.time_end:
                times.append(next_output_time)
                states.append(current_state.copy())
                next_output_time += self.dt
        
        # fill remaining time points if simulation ended early
        while next_output_time <= self.time_end:
            times.append(next_output_time)
            states.append(current_state.copy())
            next_output_time += self.dt
        
        return SimulationResult(
            time=np.array(times),
            states=np.array(states),
            species_names=initial_state.species_names,
            parameters=self.model.params
        )
    
    def _compute_propensities(self, state: PathwayState) -> np.ndarray:
        """
        compute propensities for all reactions (gillespie)
        
        args:
            state: current pathway state
            
        returns:
            array of propensities
        """
        p = self.model.params
        propensities = []
        
        # testosterone to dht
        propensities.append(p.k_5alpha_reductase * state.testosterone)
        
        # ar-testosterone binding (cytoplasm)
        propensities.append(p.k_on_testosterone * state.ar_free_cytoplasm * state.testosterone)
        propensities.append(p.k_off_testosterone * state.ar_testosterone_cytoplasm)
        
        # ar-dht binding (cytoplasm)
        propensities.append(p.k_on_dht * state.ar_free_cytoplasm * state.dht)
        propensities.append(p.k_off_dht * state.ar_dht_cytoplasm)
        
        # nuclear import
        propensities.append(p.k_nuclear_import * state.ar_testosterone_cytoplasm)
        propensities.append(p.k_nuclear_import * state.ar_dht_cytoplasm)
        
        # dimerization
        propensities.append(p.k_dimerization * state.ar_testosterone_nucleus * 
                          (state.ar_testosterone_nucleus - 1) / 2)
        propensities.append(p.k_dimerization * state.ar_dht_nucleus * 
                          (state.ar_dht_nucleus - 1) / 2)
        
        # dna binding
        propensities.append(p.k_dna_on * state.ar_testosterone_dimer_nucleus * state.are_sites_free)
        propensities.append(p.k_dna_on * state.ar_dht_dimer_nucleus * state.are_sites_free)
        
        # transcription
        propensities.append(p.k_transcription_activated * 
                          (state.ar_testosterone_dna + state.ar_dht_dna))
        
        # mrna degradation
        propensities.append(p.k_mrna_degradation * state.mrna_psa)
        
        return np.maximum(propensities, 0)  # ensure non-negative
    
    def _apply_reaction(self, state: np.ndarray, reaction_idx: int) -> np.ndarray:
        """
        apply reaction to state (gillespie)
        
        args:
            state: current state array
            reaction_idx: index of reaction to apply
            
        returns:
            updated state
        """
        new_state = state.copy()
        
        # apply stoichiometry for each reaction
        # (simplified version - full implementation would have complete stoichiometry matrix)
        
        if reaction_idx == 0:  # testosterone → dht
            new_state[2] -= 1  # testosterone
            new_state[3] += 1  # dht
        
        # add more reaction stoichiometries as needed...
        
        return new_state
    
    def run_multiple(self, 
                    n_runs: int = 100,
                    initial_state: Optional[PathwayState] = None,
                    **kwargs) -> List[SimulationResult]:
        """
        run multiple simulations (for stochastic averaging or parameter sampling)
        
        args:
            n_runs: number of simulation runs
            initial_state: initial state
            **kwargs: additional parameters
            
        returns:
            list of simulation results
        """
        results = []
        for i in range(n_runs):
            result = self.run(initial_state=initial_state, **kwargs)
            results.append(result)
        
        return results
    
    def run_time_course(self,
                       ligand_timecourse: Dict[str, np.ndarray],
                       time_points: np.ndarray) -> SimulationResult:
        """
        run simulation with time-varying ligand concentrations
        
        args:
            ligand_timecourse: dict with 'testosterone' and/or 'dht' as keys,
                             arrays of concentrations matching time_points
            time_points: time points for ligand changes
            
        returns:
            simulation result
        """
        # split simulation into segments with different ligand levels
        all_times = []
        all_states = []
        
        initial_state = self.model.get_initial_state()
        current_state = initial_state
        
        for i in range(len(time_points) - 1):
            t_start = time_points[i]
            t_end = time_points[i + 1]
            
            # update ligand concentrations
            if 'testosterone' in ligand_timecourse:
                current_state.testosterone = ligand_timecourse['testosterone'][i]
            if 'dht' in ligand_timecourse:
                current_state.dht = ligand_timecourse['dht'][i]
            
            # run segment
            segment_sim = Simulator(self.model, time_end=t_end - t_start, dt=self.dt)
            result = segment_sim.run(initial_state=current_state)
            
            # append results
            if i == 0:
                all_times.extend(result.time + t_start)
                all_states.extend(result.states)
            else:
                all_times.extend(result.time[1:] + t_start)
                all_states.extend(result.states[1:])
            
            # update state for next segment
            current_state = PathwayState.from_array(result.states[-1])
        
        return SimulationResult(
            time=np.array(all_times),
            states=np.array(all_states),
            species_names=initial_state.species_names,
            parameters=self.model.params
        )
