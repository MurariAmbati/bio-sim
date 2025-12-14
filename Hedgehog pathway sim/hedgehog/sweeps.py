from __future__ import annotations

import numpy as np
import pandas as pd

from .core import CoreMapping


def grid_surface(mapping: CoreMapping, ptch_level: float, hh_points: int = 101) -> pd.DataFrame:
    """Return dose-response curve as a dataframe (hh, smo, gli_ss)."""
    return mapping.dose_response(ptch_level=ptch_level, points=hh_points)


def gli_surface(mapping: CoreMapping, ptch_min: float, ptch_max: float, hh_points: int = 81, ptch_points: int = 61) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Compute a 2D surface gli_ss(hh, ptch).

    Returns (hh_grid, ptch_grid, gli_matrix) where gli_matrix has shape (ptch_points, hh_points).
    """
    hh_points = int(max(25, hh_points))
    ptch_points = int(max(25, ptch_points))
    ptch_min = max(0.0, float(ptch_min))
    ptch_max = max(ptch_min + 1e-6, float(ptch_max))

    hh = np.linspace(0.0, 1.0, hh_points)
    ptch = np.linspace(ptch_min, ptch_max, ptch_points)

    gli = np.zeros((ptch_points, hh_points), dtype=float)
    for i, p in enumerate(ptch):
        _smo, gli_row = mapping.transduction(hh, p)
        gli[i, :] = np.asarray(gli_row, dtype=float)

    return hh, ptch, gli
