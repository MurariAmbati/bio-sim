from __future__ import annotations

import argparse
from pathlib import Path
from typing import Dict

import matplotlib.pyplot as plt
import numpy as np
import yaml
from scipy.integrate import solve_ivp

from model import STATE_NAMES, ligand_pulse, ode_rhs


def load_params(path: Path) -> Dict:
    with path.open("r", encoding="utf-8") as f:
        return yaml.safe_load(f)


def build_y0(p: Dict) -> np.ndarray:
    init = p["init"]
    y0 = np.array([float(init[name]) for name in STATE_NAMES], dtype=float)
    return y0


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--params", type=str, default="params.yaml")
    ap.add_argument("--out", type=str, default="out_ode.png")
    args = ap.parse_args()

    params_path = Path(args.params)
    p = load_params(params_path)

    lig_cfg = p["ligand"]
    lig_fn = lambda t: ligand_pulse(t, lig_cfg)

    y0 = build_y0(p)

    t0 = float(p["sim"]["t0"])
    t_end = float(p["sim"]["t_end"])
    max_step = float(p["sim"].get("max_step", 0.5))

    rhs = lambda t, y: ode_rhs(t, y, p, lig_fn)

    sol = solve_ivp(
        rhs,
        t_span=(t0, t_end),
        y0=y0,
        method="LSODA",
        max_step=max_step,
        rtol=1e-7,
        atol=1e-9,
        dense_output=False,
    )

    if not sol.success:
        raise RuntimeError(sol.message)

    t = sol.t
    Y = sol.y
    i = {n: k for k, n in enumerate(STATE_NAMES)}

    # Observables
    L = np.array([lig_fn(tt) for tt in t])
    R_active = Y[i["Rstar"], :]
    cAMP_mem = Y[i["cAMP_mem"], :]
    cAMP_cyt = Y[i["cAMP_cyt"], :]
    PKAstar = Y[i["PKAstar"], :]
    pCREB = Y[i["pCREB"], :]
    RArr = Y[i["RArr"], :]
    Rint = Y[i["Rint"], :]

    fig, axes = plt.subplots(4, 1, figsize=(10, 10), sharex=True)

    axes[0].plot(t, L, label="Ligand L(t)")
    axes[0].set_ylabel("L")
    axes[0].legend(loc="best")

    axes[1].plot(t, R_active, label="R*")
    axes[1].plot(t, RArr, label="R:Arr")
    axes[1].plot(t, Rint, label="R_int")
    axes[1].set_ylabel("Receptor states")
    axes[1].legend(loc="best")

    axes[2].plot(t, cAMP_mem, label="cAMP_mem")
    axes[2].plot(t, cAMP_cyt, label="cAMP_cyt")
    axes[2].set_ylabel("cAMP")
    axes[2].legend(loc="best")

    axes[3].plot(t, PKAstar, label="PKA*")
    axes[3].plot(t, pCREB, label="pCREB")
    axes[3].set_ylabel("PKA / CREB")
    axes[3].set_xlabel("time (s)")
    axes[3].legend(loc="best")

    fig.tight_layout()
    out_path = Path(args.out)
    fig.savefig(out_path, dpi=180)
    print(f"Saved {out_path.resolve()}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
