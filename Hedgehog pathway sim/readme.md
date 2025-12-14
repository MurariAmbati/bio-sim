# hedgehog

slug: hedgehog

role: developmental patterning, stem cell maintenance, cancer.

this repo is a streamlit sandbox for qualitative hedgehog signaling with more depth:

hh → ptch ⟂ smo → gli

it includes:

- steady-state dose response (hh → smo, gli*)
- feedback dynamics (gli → ptch negative feedback)
- extended dynamics (gli activator/repressor + hhip ligand sink)
- 1d tissue simulation (hh diffusion + decay + hhip sink → spatial gli patterns)

## run

1. create env (optional)
2. install deps
3. start streamlit

```bash
pip install -r requirements.txt
streamlit run app.py
```

## notes

- outputs are dimensionless (0..1)
- this is not a calibrated biochemical model
- model code lives in ./hedgehog (core mapping, ode, tissue, integrators, network, sweeps, viz)
