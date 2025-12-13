from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from mapk_erk import Stimulus, default_initial_conditions, default_parameters, simulate_ssa


def main() -> None:
    p = default_parameters()
    x0 = default_initial_conditions(p)

    stim = Stimulus(kind="step", amplitude=1.0, t_on=0.0)

    res = simulate_ssa(p=p, x0=x0, stim=stim, t_span=(0.0, 120.0), volume_fL=1.0, seed=1)
    df = res.to_frame()

    print(df[["t", "ERK_PP", "DUSP"]].tail(10))


if __name__ == "__main__":
    main()
