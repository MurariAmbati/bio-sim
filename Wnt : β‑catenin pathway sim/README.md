# wnt_beta_catenin

minimal streamlit simulator for a wnt / β‑catenin pathway ode with parameter scanning for multistability.

## run

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
streamlit run app.py
```

## notes

- model is a compact ode capturing wnt-dependent inhibition of the destruction complex and a nonlinear feedback from β‑catenin onto destruction-complex activity.
- the steady-state scan integrates from multiple initial conditions per wnt input and clusters distinct endpoints to reveal possible multistability.
