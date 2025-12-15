from __future__ import annotations

import argparse

import numpy as np

from .model import ModelParams, simulate


def main() -> None:
    ap = argparse.ArgumentParser(description="simulate simplified tcr signaling")
    ap.add_argument("--t-end", type=float, default=30.0)
    ap.add_argument("--points", type=int, default=601)
    args = ap.parse_args()

    t = np.linspace(0.0, float(args.t_end), int(args.points))
    df = simulate(t, ModelParams())
    print(df.to_csv(index=False))


if __name__ == "__main__":
    main()
