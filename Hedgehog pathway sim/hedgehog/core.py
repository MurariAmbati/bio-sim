from __future__ import annotations

import numpy as np
import pandas as pd

from .math_utils import clamp01, hill
from .params import CoreParams


class CoreMapping:
    """Ligand→receptor→gli mapping shared by all simulators."""

    def __init__(self, core: CoreParams) -> None:
        self.core = core

    def transduction(self, hh: np.ndarray | float, ptch_level: np.ndarray | float) -> tuple[np.ndarray, np.ndarray]:
        """Map ligand + ptch to (smo, gli_ss) in 0..1.

        - hh reduces effective ptch repression.
        - ptch represses smo.
        - smo activates gli via Hill.
        """
        hh = np.asarray(hh)
        ptch_level = np.asarray(ptch_level)
        hh = np.clip(hh, 0.0, 1.0)
        ptch_level = np.maximum(ptch_level, 0.0)

        ptch_repression = ptch_level / (1.0 + self.core.k_hh * hh)
        smo = 1.0 / (1.0 + self.core.alpha_ptch * ptch_repression)
        gli_ss = hill(smo, n=self.core.hill_n, k=self.core.hill_k)
        return np.clip(smo, 0.0, 1.0), np.clip(gli_ss, 0.0, 1.0)

    def steady_state(self, hh: float, ptch_level: float) -> dict[str, float]:
        hh = clamp01(hh)
        ptch_level = max(0.0, float(ptch_level))
        smo, gli_ss = self.transduction(hh, ptch_level)
        ptch_repression = ptch_level / (1.0 + self.core.k_hh * hh)
        return {
            "hh": hh,
            "ptch": ptch_level,
            "ptch_repression": float(ptch_repression),
            "smo": float(np.asarray(smo)),
            "gli_ss": float(np.asarray(gli_ss)),
        }

    def dose_response(self, ptch_level: float, points: int = 201) -> pd.DataFrame:
        ptch_level = max(0.0, float(ptch_level))
        points = int(max(25, points))
        hh_grid = np.linspace(0.0, 1.0, points)
        smo_grid, gli_grid = self.transduction(hh_grid, ptch_level)
        return pd.DataFrame({"hh": hh_grid, "smo": np.asarray(smo_grid), "gli_ss": np.asarray(gli_grid)})
