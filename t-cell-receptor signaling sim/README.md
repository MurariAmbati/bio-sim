# tcr signaling (t cell receptor)

minimal streamlit simulator for a simplified t cell receptor signaling network with three readouts:
- nf-κb (ikt→nfkb)
- mapk (ras→raf→mek→erk)
- ca²⁺ (plcγ→ip3→calcium)

## quickstart

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
streamlit run app.py
```

## run

start the app:

```bash
source .venv/bin/activate
streamlit run app.py
```

open:
- http://localhost:8501

stop:
- ctrl+c in the terminal

## model scope

this is a compact, interpretable, not cell-type-specific model intended for educational exploration.
inputs:
- antigen (tcr stimulus)
- cd28 costimulation
- lck activity
- shp1 inhibition (negative feedback)

outputs:
- nfkb activity
- erk activity
- cytosolic calcium

## files

- app.py: streamlit ui
- tcr_signaling/model.py: ode model + simulation helper
- tcr_signaling/pathway.py: pathway graph definition
- tcr_signaling/plots.py: plotly helpers
