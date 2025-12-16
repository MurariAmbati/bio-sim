# hippo / yap–taz pathway

slug: hippo_yap_taz

role: organ size control, proliferation vs apoptosis, mechano-sensing

## what this is

this repo is a mechanistic (but still intentionally simplified) ode model of the hippo pathway controlling yap/taz localization and downstream transcription.

the model is designed to be:

- explicit about state variables (mst, lats, yap_c, yap_p, yap_n, gene, actin)
- explicit about couplings (stiffness → actin, density → mst, mst/lats → yap phosphorylation, actin → nuclear import)
- interactive via streamlit for exploration and hypothesis testing

## quickstart

1) create an environment and install deps:

- `python -m venv .venv`
- `source .venv/bin/activate`
- `pip install -r requirements.txt`

note: this workspace folder name contains a `:` character, and python’s `venv` refuses to create a venv inside paths containing `:`. if you hit that error, create the venv in a sibling directory (example below) and run the app with that interpreter:

- `cd ..`
- `python -m venv hippo_yap_taz_venv`
- `./hippo_yap_taz_venv/bin/python -m pip install -r "hippo : yap-taz pathway sim/requirements.txt"`

2) run the app:

- `streamlit run app.py`

or, if using the sibling venv:

- `../hippo_yap_taz_venv/bin/python -m streamlit run app.py`

## model (state + equations)

state vector:

- mst: upstream hippo kinase activity (mst1/2 proxy)
- lats: lats1/2 activity
- yap_c: unphosphorylated cytoplasmic yap/taz
- yap_p: phosphorylated cytoplasmic yap/taz
- yap_n: nuclear yap/taz (active)
- gene: representative yap/tead target gene output
- actin: mechanosensing / actin tension proxy

inputs:

- stiffness ∈ [0, 1] (0 soft, 1 stiff)
- density ∈ [0, 1] (contact inhibition / crowding)

key design choices (what to look for in the code):

- contact inhibition promotes mst activation (density → mst)
- mechanical tension suppresses mst/lats activation (actin ┤ mst, actin ┤ lats)
- lats phosphorylates yap/taz, lowering nuclear entry (lats → yap_p, yap_p ┤ yap_n)
- actin promotes nuclear import of yap/taz (actin → yap_n)
- yap_n drives gene output through a hill nonlinearity

implementation details:

- ode right-hand-side: hippo_yap_taz/sim/model.py
- simulation + steady state helper: hippo_yap_taz/sim/simulate.py
- visualizations: hippo_yap_taz/viz/plots.py

## streamlit views

the app provides:

- time series of all state variables
- a phase portrait (lats vs nuclear yap)
- a pathway network view (nodes colored by steady-state values)
- a local sensitivity plot (one-at-a-time) for steady nuclear yap fraction
- a steady-state heatmap over (stiffness, density)

## notes / limitations

- this is not a curated, isoform-accurate biochemical reconstruction; it is a compact control-oriented hippo/yap–taz model.
- parameters are dimensionless / normalized to keep exploration simple.
- if you want a closer-to-literature mapping (e.g., merlin/nf2, sav1, amot, 14-3-3 sequestration, explicit taz), we can extend the state graph and refit parameters.
