from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Tuple

import matplotlib.pyplot as plt
import numpy as np
import yaml

from model import STATE_NAMES, build_ssa_reactions, ligand_pulse


def load_params(path: Path) -> Dict:
    with path.open("r", encoding="utf-8") as f:
        return yaml.safe_load(f)


def build_x0(p: Dict) -> np.ndarray:
    # SSA uses integer molecule counts; here we just scale the ODE inits.
    init = p["init"]
    # scale factor maps "uM" → counts in an effective volume
    scale = float(p.get("ssa", {}).get("scale", 300.0))
    x0 = np.array([int(round(float(init[name]) * scale)) for name in STATE_NAMES], dtype=int)
    return x0


@dataclass
class Trajectory:
    t: np.ndarray
    x: np.ndarray


def gillespie(
    t0: float,
    t_end: float,
    x0: np.ndarray,
    p: Dict,
    seed: int,
) -> Trajectory:
    rng = np.random.default_rng(seed)

    lig_cfg = p["ligand"]
    lig_fn = lambda t: ligand_pulse(t, lig_cfg)

    names, rxns = build_ssa_reactions(p, lig_fn)

    t = float(t0)
    x = x0.copy()

    ts = [t]
    xs = [x.copy()]

    while t < t_end:
        a = np.array([max(0.0, r.propensity(t, x, p)) for r in rxns], dtype=float)
        a0 = float(a.sum())
        if a0 <= 0.0:
            break

        r1 = rng.random()
        r2 = rng.random()
        tau = -np.log(max(r1, 1e-12)) / a0
        t_next = t + tau
        if t_next > t_end:
            break

        # choose reaction
        threshold = r2 * a0
        j = int(np.searchsorted(np.cumsum(a), threshold, side="right"))
        j = min(j, len(rxns) - 1)

        # apply stoichiometry, guard against negative counts
        for (k, delta) in rxns[j].updates:
            x[k] = x[k] + int(delta)
            if x[k] < 0:
                x[k] = 0

        t = t_next
        ts.append(t)
        xs.append(x.copy())

    return Trajectory(t=np.array(ts), x=np.stack(xs, axis=0))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--params", type=str, default="params.yaml")
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--out", type=str, default="out_ssa.png")
    args = ap.parse_args()

    p = load_params(Path(args.params))

    t0 = float(p["sim"]["t0"])
    t_end = float(p["sim"]["t_end"])

    x0 = build_x0(p)

    traj = gillespie(t0=t0, t_end=t_end, x0=x0, p=p, seed=args.seed)

    i = {n: k for k, n in enumerate(STATE_NAMES)}

    t = traj.t
    X = traj.x

    # plot a few key species counts
    fig, ax = plt.subplots(1, 1, figsize=(10, 4))
    ax.step(t, X[:, i["cAMP_mem"]], where="post", label="cAMP_mem (counts)")
    ax.step(t, X[:, i["cAMP_cyt"]], where="post", label="cAMP_cyt (counts)")
    ax.step(t, X[:, i["PKAstar"]], where="post", label="PKA* (counts)")
    ax.set_xlabel("time (s)")
    ax.set_ylabel("molecule counts")
    ax.legend(loc="best")
    fig.tight_layout()

    out_path = Path(args.out)
    fig.savefig(out_path, dpi=180)
    print(f"Saved {out_path.resolve()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
