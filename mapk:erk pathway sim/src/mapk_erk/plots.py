from __future__ import annotations

from typing import Iterable, Optional

import matplotlib.pyplot as plt
import pandas as pd


def plot_timecourses(
    df: pd.DataFrame,
    species: Iterable[str],
    title: str = "",
    ax: Optional[plt.Axes] = None,
) -> plt.Axes:
    if ax is None:
        _, ax = plt.subplots(figsize=(8, 4))

    for s in species:
        ax.plot(df["t"], df[s], label=s)

    ax.set_xlabel("time")
    ax.set_ylabel("concentration")
    if title:
        ax.set_title(title)
    ax.legend(frameon=False)
    return ax
