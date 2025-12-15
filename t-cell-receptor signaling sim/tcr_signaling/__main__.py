from __future__ import annotations

from .model import ModelParams, simulate


def main() -> None:
    import numpy as np

    t = np.linspace(0.0, 30.0, 601)
    df = simulate(t, ModelParams())
    print(df.tail())


if __name__ == "__main__":
    main()
