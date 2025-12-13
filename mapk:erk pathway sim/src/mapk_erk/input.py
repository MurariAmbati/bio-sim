from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class Stimulus:
    """Time-varying ligand input L(t) (e.g., EGF).

    Concentration units are arbitrary but must be consistent with receptor activation rates.
    """

    kind: str = "step"  # step | pulse | exp
    amplitude: float = 1.0
    t_on: float = 0.0
    t_off: float = 30.0
    tau: float = 10.0  # for exp decay

    def value(self, t: float) -> float:
        if self.kind == "step":
            return float(self.amplitude) if t >= self.t_on else 0.0
        if self.kind == "pulse":
            return float(self.amplitude) if (self.t_on <= t <= self.t_off) else 0.0
        if self.kind == "exp":
            if t < self.t_on:
                return 0.0
            return float(self.amplitude) * float(np.exp(-(t - self.t_on) / max(self.tau, 1e-9)))
        raise ValueError(f"Unknown stimulus kind: {self.kind}")

    def vectorized(self, t: np.ndarray) -> np.ndarray:
        t = np.asarray(t, dtype=float)
        if self.kind == "step":
            return self.amplitude * (t >= self.t_on)
        if self.kind == "pulse":
            return self.amplitude * ((t >= self.t_on) & (t <= self.t_off))
        if self.kind == "exp":
            out = np.zeros_like(t, dtype=float)
            mask = t >= self.t_on
            out[mask] = self.amplitude * np.exp(-(t[mask] - self.t_on) / max(self.tau, 1e-9))
            return out
        raise ValueError(f"Unknown stimulus kind: {self.kind}")
