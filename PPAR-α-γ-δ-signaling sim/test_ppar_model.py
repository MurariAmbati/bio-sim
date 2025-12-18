"""
Unit tests for PPAR signaling model
"""

import numpy as np
import pytest
from ppar_model import PPARSignalingModel, PPARParameters, compare_isoform_specificity


class TestPPARModel:
    """Test suite for PPAR signaling model"""
    
    def test_initialization(self):
        """Test model initialization"""
        model = PPARSignalingModel()
        assert model is not None
        assert len(model.state_names) == 28
        
    def test_initial_conditions(self):
        """Test initial conditions are valid"""
        model = PPARSignalingModel()
        state0 = model.get_initial_conditions()
        
        assert len(state0) == 28
        assert np.all(state0 >= 0)  # All concentrations should be non-negative
        assert state0[0] > 0  # PPAR-α should be present
        assert state0[24] > 0  # RXR should be present
        
    def test_simulation_runs(self):
        """Test that simulation runs without errors"""
        model = PPARSignalingModel()
        df = model.simulate((0, 10), 100, 0.0, 0.0, 0.0)
        
        assert df is not None
        assert len(df) == 100
        assert 'Time' in df.columns
        assert 'PPAR-α' in df.columns
        
    def test_ligand_binding(self):
        """Test that ligands activate PPARs"""
        model = PPARSignalingModel()
        
        # Baseline
        df_baseline = model.simulate((0, 100), 500, 0.0, 0.0, 0.0)
        baseline_activity = df_baseline['PPAR-α-RXR-DNA'].iloc[-10:].mean()
        
        # With ligand
        df_ligand = model.simulate((0, 100), 500, 1.0, 0.0, 0.0)
        ligand_activity = df_ligand['PPAR-α-RXR-DNA'].iloc[-10:].mean()
        
        assert ligand_activity > baseline_activity
        
    def test_isoform_specificity(self):
        """Test that isoforms have specific effects"""
        model = PPARSignalingModel()
        
        # PPAR-α should increase fatty acid oxidation
        df_alpha = model.simulate((0, 100), 500, 1.5, 0.0, 0.0)
        fao_alpha = df_alpha['Protein_FAO'].iloc[-10:].mean()
        
        # PPAR-γ should increase insulin sensitivity
        df_gamma = model.simulate((0, 100), 500, 0.0, 1.5, 0.0)
        is_gamma = df_gamma['Protein_IS'].iloc[-10:].mean()
        
        # PPAR-α effect on FAO should be stronger
        fao_gamma = df_gamma['Protein_FAO'].iloc[-10:].mean()
        assert fao_alpha > fao_gamma
        
        # PPAR-γ effect on IS should be stronger
        is_alpha = df_alpha['Protein_IS'].iloc[-10:].mean()
        assert is_gamma > is_alpha
        
    def test_dose_response(self):
        """Test dose-response curves"""
        model = PPARSignalingModel()
        doses = np.array([0.0, 0.5, 1.0, 2.0])
        results = model.simulate_drug_response('alpha', doses, (0, 50), 200)
        
        assert len(results) == len(doses)
        
        # Activity should increase with dose
        activities = []
        for key in results.keys():
            df = results[key]
            activity = df['PPAR-α-RXR-DNA'].iloc[-10:].mean()
            activities.append(activity)
        
        # Check monotonic increase
        for i in range(len(activities)-1):
            assert activities[i+1] >= activities[i]
    
    def test_steady_state_metrics(self):
        """Test steady-state metric extraction"""
        model = PPARSignalingModel()
        df = model.simulate((0, 100), 500, 1.0, 1.0, 0.0)
        metrics = model.get_steady_state_metrics(df)
        
        assert 'PPAR-α Activity' in metrics
        assert 'Insulin Sensitivity' in metrics
        assert 'Metabolic Health' in metrics
        assert all(isinstance(v, float) for v in metrics.values())
        
    def test_parameter_modification(self):
        """Test custom parameter initialization"""
        params = PPARParameters(
            ppar_alpha_base=2.0,
            k_transcription=2.0
        )
        model = PPARSignalingModel(params)
        
        assert model.params.ppar_alpha_base == 2.0
        assert model.params.k_transcription == 2.0
        
    def test_comparison_function(self):
        """Test isoform comparison"""
        df = compare_isoform_specificity()
        
        assert df is not None
        assert len(df) == 3  # Three isoforms
        assert 'Isoform' in df.columns
        assert 'PPAR-α Activity' in df.columns
        
    def test_conservation_laws(self):
        """Test that key conservation laws hold"""
        model = PPARSignalingModel()
        df = model.simulate((0, 100), 500, 0.5, 0.5, 0.5)
        
        # RXR conservation (approximately)
        rxr_total = (df['RXR_free'] + 
                    df['PPAR-α-RXR'] + df['PPAR-γ-RXR'] + df['PPAR-δ-RXR'] +
                    df['PPAR-α-RXR-DNA'] + df['PPAR-γ-RXR-DNA'] + df['PPAR-δ-RXR-DNA'])
        
        # Should be relatively constant
        rxr_std = rxr_total.std()
        rxr_mean = rxr_total.mean()
        assert rxr_std / rxr_mean < 0.3  # Less than 30% variation
        
    def test_non_negative_states(self):
        """Test that all states remain non-negative"""
        model = PPARSignalingModel()
        df = model.simulate((0, 100), 500, 2.0, 2.0, 2.0)
        
        # Check all state variables
        for col in model.state_names:
            if col in df.columns:
                assert np.all(df[col] >= -1e-6)  # Allow small numerical errors


class TestPPARParameters:
    """Test parameter class"""
    
    def test_default_parameters(self):
        """Test default parameter values"""
        params = PPARParameters()
        
        assert params.ppar_alpha_base == 1.0
        assert params.ppar_gamma_base == 1.0
        assert params.ppar_delta_base == 1.0
        assert params.rxr_total == 2.0
        
    def test_custom_parameters(self):
        """Test custom parameter setting"""
        params = PPARParameters(
            ppar_alpha_base=0.5,
            k_transcription=1.5,
            free_fatty_acids=2.0
        )
        
        assert params.ppar_alpha_base == 0.5
        assert params.k_transcription == 1.5
        assert params.free_fatty_acids == 2.0


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
