from __future__ import annotations

import numpy as np


def clamp01(x: float) -> float:
    return float(max(0.0, min(1.0, x)))


def hill(x: np.ndarray | float, n: float, k: float) -> np.ndarray | float:
    x = np.asarray(x)
    xn = np.power(np.maximum(x, 0.0), n)
    kn = k**n
    return xn / (xn + kn)
