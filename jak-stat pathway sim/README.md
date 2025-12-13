# jak_stat

jak–stat pathway simulator (ode + discrete/logical).

role: cytokine and growth factor signaling; immune responses.

## quickstart

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
streamlit run app/streamlit_app.py
```

## models

- ode: cytokine→receptor complex→jak activity→pstat→dimer→nuclear stat; socs negative feedback.
- logical: boolean network with synchronous updates; same wiring and feedback.

## scenarios

- baseline: moderate cytokine input.
- pulse: cytokine pulse for a chosen duration.
- socs knockout: disables negative feedback.

## layout

- jak_stat/model_ode.py: ode rhs + simulation helper
- jak_stat/model_logical.py: logical update rules + simulation helper
- jak_stat/viz.py: plotting helpers
- app/streamlit_app.py: interactive ui
