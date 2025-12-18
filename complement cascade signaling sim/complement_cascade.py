"""
Complement Cascade Signaling Simulation
Comprehensive model of classical, alternative, and lectin pathways
"""

import numpy as np
from dataclasses import dataclass, field
from typing import Dict, List, Tuple
from enum import Enum


class ComplementPathway(Enum):
    CLASSICAL = "classical"
    ALTERNATIVE = "alternative"
    LECTIN = "lectin"


@dataclass
class ComplementComponent:
    """Represents a complement protein component"""
    name: str
    concentration: float
    active: bool = False
    location: str = "plasma"  # plasma, membrane, surface
    
    
@dataclass
class PathogenSurface:
    """Represents a pathogen surface with complement binding sites"""
    pathogen_id: str
    surface_area: float
    c3b_deposited: int = 0
    mac_complexes: int = 0
    opsonized: bool = False
    lysed: bool = False
    

@dataclass
class SimulationState:
    """Complete state of the complement cascade simulation"""
    time: float = 0.0
    components: Dict[str, ComplementComponent] = field(default_factory=dict)
    pathogens: List[PathogenSurface] = field(default_factory=list)
    anaphylatoxins: Dict[str, float] = field(default_factory=lambda: {"C3a": 0.0, "C5a": 0.0})
    inflammatory_response: float = 0.0
    

class ComplementCascade:
    """
    Comprehensive complement cascade simulation
    Models all three activation pathways and downstream effects
    """
    
    def __init__(self):
        self.state = SimulationState()
        self.history: List[SimulationState] = []
        self.dt = 0.01  # time step - reduced for stability
        
        # Initialize complement components with physiological concentrations (μg/mL)
        self._initialize_components()
        
        # Rate constants
        self.rates = {
            # Classical pathway
            "c1_activation": 0.8,
            "c4_cleavage": 1.2,
            "c2_cleavage": 1.0,
            "c3_conv_classical": 2.5,
            
            # Alternative pathway
            "c3_tickover": 0.05,
            "factor_b_binding": 0.3,
            "factor_d_cleavage": 0.9,
            "c3_conv_alternative": 2.0,
            "properdin_stabilization": 0.4,
            
            # Lectin pathway
            "mbl_activation": 0.7,
            "masp_activation": 1.1,
            "c3_conv_lectin": 2.2,
            
            # Terminal pathway (MAC formation)
            "c5_cleavage": 1.5,
            "c6_binding": 2.0,
            "c7_binding": 1.8,
            "c8_binding": 1.6,
            "c9_polymerization": 2.5,
            
            # Regulatory
            "c3b_decay": 0.02,
            "anaphylatoxin_clearance": 0.1,
            "mac_lysis": 0.3,
        }
        
    def _initialize_components(self):
        """Initialize all complement components"""
        components = {
            # Classical pathway
            "C1q": ComplementComponent("C1q", 70.0),
            "C1r": ComplementComponent("C1r", 34.0),
            "C1s": ComplementComponent("C1s", 31.0),
            "C4": ComplementComponent("C4", 600.0),
            "C2": ComplementComponent("C2", 25.0),
            "C4b": ComplementComponent("C4b", 0.0),
            "C2a": ComplementComponent("C2a", 0.0),
            
            # Alternative pathway
            "C3": ComplementComponent("C3", 1300.0),
            "Factor B": ComplementComponent("Factor B", 200.0),
            "Factor D": ComplementComponent("Factor D", 1.8),
            "Properdin": ComplementComponent("Properdin", 25.0),
            "C3b": ComplementComponent("C3b", 0.0),
            "Bb": ComplementComponent("Bb", 0.0),
            
            # Lectin pathway
            "MBL": ComplementComponent("MBL", 1.0),
            "MASP-1": ComplementComponent("MASP-1", 12.0),
            "MASP-2": ComplementComponent("MASP-2", 0.5),
            
            # Terminal pathway
            "C5": ComplementComponent("C5", 85.0),
            "C6": ComplementComponent("C6", 45.0),
            "C7": ComplementComponent("C7", 55.0),
            "C8": ComplementComponent("C8", 55.0),
            "C9": ComplementComponent("C9", 60.0),
            "C5b": ComplementComponent("C5b", 0.0),
            "MAC": ComplementComponent("MAC", 0.0, location="membrane"),
            
            # Regulators
            "Factor H": ComplementComponent("Factor H", 500.0),
            "Factor I": ComplementComponent("Factor I", 35.0),
            "C4BP": ComplementComponent("C4BP", 250.0),
            "DAF": ComplementComponent("DAF", 0.0, location="membrane"),
        }
        self.state.components = components
        
    def add_pathogen(self, pathogen_id: str, surface_area: float = 10.0):
        """Add a pathogen surface to the simulation"""
        pathogen = PathogenSurface(pathogen_id, surface_area)
        self.state.pathogens.append(pathogen)
        
    def activate_classical_pathway(self, intensity: float = 1.0):
        """Activate the classical pathway (antibody-antigen complex)"""
        c1q = self.state.components["C1q"]
        c1q.active = True
        c1q.concentration *= intensity
        
    def activate_lectin_pathway(self, intensity: float = 1.0):
        """Activate the lectin pathway (pathogen surface carbohydrates)"""
        mbl = self.state.components["MBL"]
        mbl.active = True
        mbl.concentration *= intensity
        
    def activate_alternative_pathway(self):
        """Activate the alternative pathway (spontaneous C3 hydrolysis)"""
        # Alternative pathway is always active through C3 tickover
        pass
        
    def _classical_pathway_step(self):
        """Simulate one time step of the classical pathway"""
        c1q = self.state.components["C1q"]
        c4 = self.state.components["C4"]
        c2 = self.state.components["C2"]
        c4b = self.state.components["C4b"]
        c2a = self.state.components["C2a"]
        c3 = self.state.components["C3"]
        c3b = self.state.components["C3b"]
        
        if c1q.active:
            # C1 complex cleaves C4 → C4a + C4b
            c4_cleavage = self.rates["c4_cleavage"] * c1q.concentration * c4.concentration * self.dt
            c4.concentration -= c4_cleavage
            c4b.concentration += c4_cleavage
            
            # C4b binds C2, C1s cleaves C2 → C2a + C2b
            if c4b.concentration > 0:
                c2_cleavage = self.rates["c2_cleavage"] * c4b.concentration * c2.concentration * self.dt
                c2.concentration -= c2_cleavage
                c2a.concentration += c2_cleavage
                
                # C4b2a (C3 convertase) cleaves C3 → C3a + C3b
                if c2a.concentration > 0:
                    c3_convertase_activity = c4b.concentration * c2a.concentration
                    c3_cleavage = self.rates["c3_conv_classical"] * c3_convertase_activity * c3.concentration * self.dt
                    c3_cleavage = min(c3_cleavage, c3.concentration)  # can't cleave more than available
                    c3.concentration -= c3_cleavage
                    c3b.concentration += c3_cleavage
                    self.state.anaphylatoxins["C3a"] += c3_cleavage
                    
    def _alternative_pathway_step(self):
        """Simulate one time step of the alternative pathway"""
        c3 = self.state.components["C3"]
        c3b = self.state.components["C3b"]
        factor_b = self.state.components["Factor B"]
        factor_d = self.state.components["Factor D"]
        bb = self.state.components["Bb"]
        properdin = self.state.components["Properdin"]
        
        # Spontaneous C3 tickover (hydrolysis)
        c3_tickover = self.rates["c3_tickover"] * c3.concentration * self.dt
        c3.concentration -= c3_tickover
        c3b.concentration += c3_tickover
        self.state.anaphylatoxins["C3a"] += c3_tickover
        
        # C3b binds Factor B
        if c3b.concentration > 0:
            fb_binding = self.rates["factor_b_binding"] * c3b.concentration * factor_b.concentration * self.dt
            factor_b.concentration -= fb_binding
            
            # Factor D cleaves Factor B → Ba + Bb (C3bBb = C3 convertase)
            bb_formation = self.rates["factor_d_cleavage"] * fb_binding * factor_d.concentration * self.dt
            bb.concentration += bb_formation
            
            # C3bBb convertase cleaves more C3 (amplification loop)
            if bb.concentration > 0:
                c3_convertase_activity = c3b.concentration * bb.concentration
                
                # Properdin stabilizes C3bBb
                stabilization = 1.0 + self.rates["properdin_stabilization"] * properdin.concentration
                c3_convertase_activity *= stabilization
                
                c3_cleavage = self.rates["c3_conv_alternative"] * c3_convertase_activity * c3.concentration * self.dt
                c3_cleavage = min(c3_cleavage, c3.concentration)  # can't cleave more than available
                c3.concentration -= c3_cleavage
                c3b.concentration += c3_cleavage
                self.state.anaphylatoxins["C3a"] += c3_cleavage
                
    def _lectin_pathway_step(self):
        """Simulate one time step of the lectin pathway"""
        mbl = self.state.components["MBL"]
        masp2 = self.state.components["MASP-2"]
        c4 = self.state.components["C4"]
        c2 = self.state.components["C2"]
        c4b = self.state.components["C4b"]
        c2a = self.state.components["C2a"]
        c3 = self.state.components["C3"]
        c3b = self.state.components["C3b"]
        
        if mbl.active:
            # MBL-MASP complex cleaves C4 and C2 (similar to C1)
            c4_cleavage = self.rates["masp_activation"] * mbl.concentration * masp2.concentration * c4.concentration * self.dt
            c4.concentration -= c4_cleavage
            c4b.concentration += c4_cleavage
            
            if c4b.concentration > 0:
                c2_cleavage = self.rates["masp_activation"] * c4b.concentration * c2.concentration * self.dt
                c2.concentration -= c2_cleavage
                c2a.concentration += c2_cleavage
                
                # C4b2a cleaves C3
                if c2a.concentration > 0:
                    c3_convertase_activity = c4b.concentration * c2a.concentration
                    c3_cleavage = self.rates["c3_conv_lectin"] * c3_convertase_activity * c3.concentration * self.dt
                    c3_cleavage = min(c3_cleavage, c3.concentration)  # can't cleave more than available
                    c3.concentration -= c3_cleavage
                    c3b.concentration += c3_cleavage
                    self.state.anaphylatoxins["C3a"] += c3_cleavage
                    
    def _terminal_pathway_step(self):
        """Simulate MAC formation and membrane lysis"""
        c3b = self.state.components["C3b"]
        c5 = self.state.components["C5"]
        c5b = self.state.components["C5b"]
        c6 = self.state.components["C6"]
        c7 = self.state.components["C7"]
        c8 = self.state.components["C8"]
        c9 = self.state.components["C9"]
        mac = self.state.components["MAC"]
        
        # C3 convertase + C3b forms C5 convertase, cleaves C5 → C5a + C5b
        if c3b.concentration > 10.0:  # threshold for C5 convertase formation
            c5_cleavage = self.rates["c5_cleavage"] * c3b.concentration * c5.concentration * self.dt
            c5.concentration -= c5_cleavage
            c5b.concentration += c5_cleavage
            self.state.anaphylatoxins["C5a"] += c5_cleavage
            
            # C5b initiates MAC assembly
            if c5b.concentration > 0:
                # C5b-6 complex
                c6_binding = self.rates["c6_binding"] * c5b.concentration * c6.concentration * self.dt
                
                # C5b-6-7 complex
                c7_binding = self.rates["c7_binding"] * c6_binding * c7.concentration * self.dt
                
                # C5b-6-7-8 complex
                c8_binding = self.rates["c8_binding"] * c7_binding * c8.concentration * self.dt
                
                # C5b-6-7-8-9n (multiple C9 polymerize to form pore)
                c9_poly = self.rates["c9_polymerization"] * c8_binding * c9.concentration * self.dt
                
                mac.concentration += c9_poly
                c6.concentration -= c6_binding
                c7.concentration -= c7_binding
                c8.concentration -= c8_binding
                c9.concentration -= c9_poly * 10  # ~10-18 C9 molecules per MAC
                
    def _opsonization_step(self):
        """Handle C3b-mediated opsonization of pathogens"""
        c3b = self.state.components["C3b"]
        
        for pathogen in self.state.pathogens:
            if not pathogen.opsonized and c3b.concentration > 5.0:
                # C3b deposits on pathogen surface
                deposition_rate = 0.5 * c3b.concentration * pathogen.surface_area * self.dt
                c3b_deposited = int(deposition_rate)
                
                if c3b_deposited > 0:
                    pathogen.c3b_deposited += c3b_deposited
                    c3b.concentration -= c3b_deposited * 0.1  # consume C3b
                    
                    # Mark as opsonized when sufficient C3b coating
                    if pathogen.c3b_deposited > 50:
                        pathogen.opsonized = True
                        
    def _lysis_step(self):
        """Handle MAC-mediated cell lysis"""
        mac = self.state.components["MAC"]
        
        for pathogen in self.state.pathogens:
            if not pathogen.lysed and mac.concentration > 0:
                # MAC inserts into pathogen membrane
                insertion_rate = self.rates["mac_lysis"] * mac.concentration * pathogen.surface_area * self.dt
                mac_inserted = int(insertion_rate)
                
                if mac_inserted > 0:
                    pathogen.mac_complexes += mac_inserted
                    mac.concentration -= mac_inserted * 0.2
                    
                    # Lysis occurs with sufficient MAC pores
                    if pathogen.mac_complexes > 10:
                        pathogen.lysed = True
                        
    def _inflammation_step(self):
        """Calculate inflammatory response from anaphylatoxins"""
        # C5a is ~100x more potent than C3a
        c3a_effect = self.state.anaphylatoxins["C3a"] * 0.01
        c5a_effect = self.state.anaphylatoxins["C5a"] * 1.0
        
        # Inflammation increases with anaphylatoxin concentration
        inflammation_increase = (c3a_effect + c5a_effect) * self.dt
        self.state.inflammatory_response += inflammation_increase
        
        # Anaphylatoxins are cleared over time
        self.state.anaphylatoxins["C3a"] *= (1 - self.rates["anaphylatoxin_clearance"] * self.dt)
        self.state.anaphylatoxins["C5a"] *= (1 - self.rates["anaphylatoxin_clearance"] * self.dt)
        
    def _regulatory_step(self):
        """Apply regulatory mechanisms"""
        c3b = self.state.components["C3b"]
        factor_h = self.state.components["Factor H"]
        factor_i = self.state.components["Factor I"]
        
        # Factor H + Factor I degrade C3b → iC3b (inactive)
        c3b_decay = self.rates["c3b_decay"] * c3b.concentration * factor_h.concentration * factor_i.concentration * self.dt
        c3b.concentration -= c3b_decay
    
    def _clamp_concentrations(self):
        """Clamp all concentrations to reasonable bounds to prevent overflow"""
        MAX_CONCENTRATION = 1e6  # Maximum reasonable concentration
        
        for component in self.state.components.values():
            # Ensure non-negative
            component.concentration = max(0.0, component.concentration)
            # Cap at maximum
            component.concentration = min(component.concentration, MAX_CONCENTRATION)
        
        # Clamp anaphylatoxins
        for key in self.state.anaphylatoxins:
            self.state.anaphylatoxins[key] = max(0.0, min(self.state.anaphylatoxins[key], MAX_CONCENTRATION))
        
        # Clamp inflammation
        self.state.inflammatory_response = max(0.0, min(self.state.inflammatory_response, MAX_CONCENTRATION))
        
    def step(self):
        """Execute one simulation time step"""
        # Execute all pathway steps
        self._classical_pathway_step()
        self._alternative_pathway_step()
        self._lectin_pathway_step()
        self._terminal_pathway_step()
        
        # Execute downstream effects
        self._opsonization_step()
        self._lysis_step()
        self._inflammation_step()
        
        # Apply regulation
        self._regulatory_step()
        
        # Clamp all concentrations to prevent overflow
        self._clamp_concentrations()
        
        # Update time
        self.state.time += self.dt
        
        # Store history (deep copy key metrics)
        self.history.append(self._snapshot_state())
        
    def _snapshot_state(self) -> Dict:
        """Create a snapshot of current state for history"""
        return {
            "time": self.state.time,
            "C3": self.state.components["C3"].concentration,
            "C3b": self.state.components["C3b"].concentration,
            "C5": self.state.components["C5"].concentration,
            "MAC": self.state.components["MAC"].concentration,
            "C3a": self.state.anaphylatoxins["C3a"],
            "C5a": self.state.anaphylatoxins["C5a"],
            "inflammation": self.state.inflammatory_response,
            "opsonized_count": sum(1 for p in self.state.pathogens if p.opsonized),
            "lysed_count": sum(1 for p in self.state.pathogens if p.lysed),
        }
        
    def run(self, duration: float):
        """Run simulation for specified duration"""
        steps = int(duration / self.dt)
        for _ in range(steps):
            self.step()
            
    def get_history_dataframe(self):
        """Convert history to pandas DataFrame for analysis"""
        import pandas as pd
        return pd.DataFrame(self.history)
        
    def reset(self):
        """Reset simulation to initial state"""
        self.__init__()


def run_standard_simulation() -> ComplementCascade:
    """Run a standard complement cascade activation scenario"""
    sim = ComplementCascade()
    
    # Add pathogens
    for i in range(5):
        sim.add_pathogen(f"pathogen_{i}", surface_area=10.0)
    
    # Activate classical pathway (antibody recognition)
    sim.activate_classical_pathway(intensity=1.5)
    
    # Run simulation
    sim.run(duration=50.0)
    
    return sim


if __name__ == "__main__":
    sim = run_standard_simulation()
    df = sim.get_history_dataframe()
    
    print(f"\nComplement Cascade Simulation Results")
    print(f"Duration: {sim.state.time:.1f} time units")
    print(f"\nFinal State:")
    print(f"  C3b deposited: {sim.state.components['C3b'].concentration:.2f}")
    print(f"  MAC formed: {sim.state.components['MAC'].concentration:.2f}")
    print(f"  Opsonized pathogens: {df['opsonized_count'].iloc[-1]}")
    print(f"  Lysed pathogens: {df['lysed_count'].iloc[-1]}")
    print(f"  Inflammatory response: {sim.state.inflammatory_response:.2f}")
    print(f"  C5a concentration: {sim.state.anaphylatoxins['C5a']:.2f}")
