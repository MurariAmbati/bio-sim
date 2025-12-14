from __future__ import annotations

from dataclasses import dataclass
from typing import Protocol

import numpy as np


class OdeSystem(Protocol):
    """Vector-field interface for ODE integration."""

    def f(self, t: float, y: np.ndarray) -> np.ndarray: ...


@dataclass(frozen=True)
class IntegratorConfig:
    dt: float
    t_end: float
    t0: float = 0.0


class EulerIntegrator:
    """Explicit Euler integrator (fast, simple, less stable)."""

    def __init__(self, config: IntegratorConfig) -> None:
        self.config = config

    def solve(self, system: OdeSystem, y0: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        dt = max(1e-6, float(self.config.dt))
        t_end = max(dt, float(self.config.t_end))
        t0 = float(self.config.t0)

        steps = int(np.ceil((t_end - t0) / dt))
        t = t0 + dt * np.arange(steps + 1)

        y = np.zeros((steps + 1, int(y0.size)), dtype=float)
        y[0] = np.asarray(y0, dtype=float)

        for i in range(steps):
            y[i + 1] = y[i] + dt * system.f(float(t[i]), y[i])

        return t, y


class RK4Integrator:
    """Classic 4th-order Runge–Kutta integrator (more stable than Euler)."""

    def __init__(self, config: IntegratorConfig) -> None:
        self.config = config

    def solve(self, system: OdeSystem, y0: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        dt = max(1e-6, float(self.config.dt))
        t_end = max(dt, float(self.config.t_end))
        t0 = float(self.config.t0)

        steps = int(np.ceil((t_end - t0) / dt))
        t = t0 + dt * np.arange(steps + 1)

        y = np.zeros((steps + 1, int(y0.size)), dtype=float)
        y[0] = np.asarray(y0, dtype=float)

        for i in range(steps):
            ti = float(t[i])
            yi = y[i]
            k1 = system.f(ti, yi)
            k2 = system.f(ti + 0.5 * dt, yi + 0.5 * dt * k1)
            k3 = system.f(ti + 0.5 * dt, yi + 0.5 * dt * k2)
            k4 = system.f(ti + dt, yi + dt * k3)
            y[i + 1] = yi + (dt / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)

        return t, y
