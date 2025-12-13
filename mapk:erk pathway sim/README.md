# mapk/erk pathway simulator

deterministic ode model + optional stochastic (ssa) mode + streamlit ui.

## run

requirements: python 3.10+

note: this folder name contains `:`. on macos, python `venv` refuses to create environments inside paths containing `:`. either rename the folder or create the venv elsewhere.

### 1) create env (recommended)

```bash
python3 -m venv /users/murari/biosim/.venvs/mapk_erk
source /users/murari/biosim/.venvs/mapk_erk/bin/activate
```

### 2) install deps

```bash
python -m pip install -r requirements.txt
```

### 3) run the app

```bash
streamlit run app/streamlit_app.py
```

### 4) run examples

```bash
python examples/run_ode.py
python examples/run_ssa.py
```
