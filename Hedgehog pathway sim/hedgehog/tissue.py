from __future__ import annotations

import numpy as np
import pandas as pd

from .core import CoreMapping
from .math_utils import hill
from .params import SpatialParams


def _stable_dt(dx: float, D: float) -> float:
    # explicit diffusion stability (rough): dt <= dx^2 / (2D)
    if D <= 0:
        return float("inf")
    return (dx * dx) / (2.0 * D)


class Tissue1DModel:
    """1D diffusion + signaling tissue model."""

    def __init__(self, params: SpatialParams) -> None:
        self.p = params
        self.mapping = CoreMapping(params.core)

    def simulate(self, nx: int, t_end: float, dt: float, ptch_level: float) -> tuple[pd.DataFrame, pd.DataFrame]:
        nx = int(max(25, min(401, nx)))
        dt = max(1e-4, float(dt))
        t_end = max(dt, float(t_end))
        ptch_level = max(0.0, float(ptch_level))

        x = np.linspace(0.0, 1.0, nx)
        dx = float(x[1] - x[0])

        dt_max = _stable_dt(dx, self.p.D)
        if dt > dt_max:
            dt = 0.95 * dt_max

        steps = int(np.ceil(t_end / dt))
        t = np.linspace(0.0, steps * dt, steps + 1)

        hh = np.zeros(nx)
        gli = np.zeros(nx)
        hhip = np.zeros(nx)

        src = x <= float(self.p.source_width)

        keep = int(max(1, steps // 240))
        rows: list[dict[str, float]] = []

        for i in range(steps + 1):
            if i % keep == 0 or i == steps:
                for j in range(nx):
                    rows.append(
                        {
                            "t": float(t[i]),
                            "x": float(x[j]),
                            "hh": float(hh[j]),
                            "gli": float(gli[j]),
                            "hhip": float(hhip[j]),
                        }
                    )

            if i == steps:
                break

            hh[src] = float(np.clip(self.p.source_hh, 0.0, 1.0))

            lap = np.zeros_like(hh)
            lap[1:-1] = (hh[2:] - 2.0 * hh[1:-1] + hh[:-2]) / (dx * dx)
            lap[0] = (hh[1] - hh[0]) / (dx * dx)
            lap[-1] = (hh[-2] - hh[-1]) / (dx * dx)

            sink = self.p.sink_strength * hhip * hh
            dhh = self.p.D * lap - self.p.k_decay * hh - sink
            hh = np.clip(hh + dt * dhh, 0.0, 1.0)
            hh[src] = float(np.clip(self.p.source_hh, 0.0, 1.0))

            _, gli_ss = self.mapping.transduction(hh, ptch_level)
            dgli = (gli_ss - gli) / max(1e-4, float(self.p.tau_gli))
            gli = np.clip(gli + dt * dgli, 0.0, 1.0)

            h_act = hill(gli, n=self.p.hhip_n, k=self.p.hhip_k)
            dh = self.p.k_hhip_prod * h_act - self.p.k_hhip_deg * hhip
            hhip = np.clip(hhip + dt * dh, 0.0, 3.0)

        space_time = pd.DataFrame(rows)
        final = pd.DataFrame({"x": x, "hh": hh, "gli": gli, "hhip": hhip})
        return space_time, final
