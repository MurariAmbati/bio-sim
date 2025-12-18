# c-kit / stem cell factor signaling simulation

comprehensive mathematical model of c-kit receptor tyrosine kinase and stem cell factor (scf) signaling pathway for hematopoietic and mast cell regulation.

## overview

this simulation models the complex molecular interactions and downstream signaling cascades initiated by scf binding to the c-kit receptor, including:

- **receptor dynamics**: binding, dimerization, and phosphorylation
- **pi3k/akt pathway**: cell survival and growth signaling
- **mapk/erk pathway**: proliferation and differentiation control  
- **jak/stat pathway**: transcriptional regulation
- **cell fate outcomes**: proliferation, differentiation, survival, apoptosis

## features

### biological accuracy
- based on established kinetic parameters from literature
- models receptor internalization and degradation
- captures pathway crosstalk and feedback loops
- physiologically relevant time scales (minutes to hours)

### simulation capabilities
- multiple scf stimulation profiles (constant, pulse, oscillatory, ramp, decay)
- adjustable kinetic parameters for sensitivity analysis
- high-resolution temporal dynamics (100-2000 time points)
- comprehensive state variable tracking (16 components)

### visualizations
- **receptor dynamics**: real-time tracking of receptor states
- **pathway activation**: temporal profiles of signaling molecules
- **cell fate outcomes**: integrated readouts of cellular decisions
- **heatmaps**: spatiotemporal patterns of pathway activation
- **phase portraits**: dynamic relationships between variables
- **correlation analysis**: pathway interaction networks

## installation

### prerequisites
- python 3.8+
- pip package manager

### setup

```bash
# clone or download the repository
cd "c-KIT -stem cell factor signaling sim"

# create virtual environment (recommended)
python -m venv venv
source venv/bin/activate  # on windows: venv\Scripts\activate

# install dependencies
pip install -r requirements.txt
```

## usage

### basic execution

```bash
streamlit run app.py
```

the application will launch in your default web browser at `http://localhost:8501`

### configuration

#### sidebar controls

**scf stimulation**
- profile type: select stimulation pattern
- amplitude: scf concentration (0-5 au)
- duration/period: temporal characteristic (10-200 min)

**time course**
- simulation duration: total time span (50-500 min)
- time steps: resolution (100-2000 points)

**advanced parameters**
- binding rate: scf-c-kit association kinetics
- phosphorylation rate: receptor activation speed
- pathway activations: akt, erk, stat kinetics

#### simulation workflow

1. configure scf stimulation profile
2. adjust temporal parameters
3. (optional) modify advanced kinetic parameters
4. click "run simulation"
5. explore results across multiple visualization tabs

## model architecture

### mathematical framework

the model employs a system of ordinary differential equations (odes) describing mass-action kinetics and michaelis-menten dynamics:

```
d[X]/dt = production - degradation - modification
```

### state variables (16 total)

**receptor layer**
- scf (free ligand)
- c-kit (free receptor)
- scf-c-kit complex
- c-kit dimer
- phospho-c-kit (activated receptor)

**signaling layer**
- pi3k active
- akt active
- ras active
- mek active
- erk active
- jak active
- stat active

**output layer**
- proliferation signal
- differentiation signal
- survival signal
- apoptosis signal

### kinetic parameters

#### receptor dynamics
- `k_bind`: scf-c-kit binding rate (0.1 min⁻¹)
- `k_unbind`: dissociation rate (0.05 min⁻¹)
- `k_dim`: dimerization rate (0.2 min⁻¹)
- `k_phos`: phosphorylation rate (0.15 min⁻¹)
- `k_dephos`: dephosphorylation rate (0.08 min⁻¹)

#### pathway activation
- `k_pi3k_act`: pi3k activation (0.12 min⁻¹)
- `k_akt_act`: akt activation (0.18 min⁻¹)
- `k_ras_act`: ras activation (0.14 min⁻¹)
- `k_mek_act`: mek activation (0.16 min⁻¹)
- `k_erk_act`: erk activation (0.20 min⁻¹)
- `k_jak_act`: jak activation (0.13 min⁻¹)
- `k_stat_act`: stat activation (0.17 min⁻¹)

#### cell outcomes
- `k_prolif`: proliferation rate (0.05 min⁻¹)
- `k_diff`: differentiation rate (0.03 min⁻¹)
- `k_survival`: survival signal (0.04 min⁻¹)
- `k_apop`: apoptosis rate (0.02 min⁻¹)

## output interpretation

### summary metrics
- **peak erk activation**: maximum mapk pathway response
- **peak akt activation**: maximum survival signaling
- **final proliferation**: integrated growth signal
- **final differentiation**: integrated maturation signal

### visualization tabs

#### receptor dynamics
tracks the molecular events at the receptor level:
- free c-kit receptor availability
- scf-c-kit complex formation
- receptor dimerization
- phosphorylation state

#### signaling pathways
monitors activation of three major cascades:
- pi3k → akt (survival/metabolism)
- ras → mek → erk (proliferation)
- jak → stat (transcription)

includes peak timing and amplitude analysis

#### cell fate
integrates pathway outputs into biological outcomes:
- proliferation (cell division)
- differentiation (maturation)
- survival (anti-apoptotic)
- apoptosis (cell death)

pie chart shows relative contribution of each fate

#### heatmap
spatiotemporal visualization of pathway activation:
- color intensity represents activation level
- time progression along x-axis
- pathway components along y-axis
- correlation matrix shows pathway interactions

#### phase space
dynamic trajectory analysis:
- 2d projections of system state
- temporal evolution shown by color gradient
- start (green star) and end (red x) markers
- reveals limit cycles, attractors, and bifurcations

## data export

results can be downloaded in two formats:
- **csv**: tabular data for external analysis
- **json**: structured data with metadata

both formats include all 16 state variables at each time point.

## biological context

### c-kit receptor
- type iii receptor tyrosine kinase
- cd117 surface marker
- 976 amino acids, ~145 kda
- expressed on hematopoietic stem cells, mast cells, germ cells

### stem cell factor (scf)
- also known as kit ligand, steel factor, mast cell growth factor
- exists in membrane-bound and soluble forms
- critical for hematopoiesis and mast cell development
- chromosome 12q22-q24 (human)

### clinical relevance
- **gain-of-function mutations**: gastrointestinal stromal tumors (gist), mastocytosis
- **loss-of-function mutations**: piebaldism, anemia
- **therapeutic targeting**: imatinib, sunitinib (c-kit inhibitors)

### hematopoietic role
- maintains hematopoietic stem cell (hsc) pool
- promotes myeloid and erythroid progenitor proliferation
- synergizes with other cytokines (il-3, gm-csf, epo)
- regulates hsc mobilization and homing

### mast cell function
- essential for mast cell development and survival
- enhances ige-mediated degranulation
- promotes chemotaxis and cytokine production
- modulates inflammatory responses

## technical details

### numerical integration
- method: `scipy.integrate.odeint` (lsoda algorithm)
- adaptive step size for efficiency
- handles stiff ode systems
- absolute and relative tolerance: 1e-6

### performance optimization
- vectorized numpy operations
- efficient dataframe manipulation with pandas
- gpu-accelerated plotting with plotly
- session state caching for responsiveness

### code structure
```
app.py
├── signaling_parameters (dataclass)
├── ckit_signaling_model (class)
│   ├── __init__
│   ├── model_equations (ode system)
│   └── simulate (integration wrapper)
├── create_scf_profile (stimulation generator)
├── plot_receptor_dynamics
├── plot_signaling_pathways
├── plot_cell_fate
├── plot_pathway_heatmap
├── plot_phase_portrait
└── main (streamlit app)
```

## customization

### adding new pathways
1. extend `model_equations` with new differential equations
2. add parameters to `signaling_parameters`
3. update initial conditions array
4. create visualization function
5. add to appropriate tab

### modifying kinetics
adjust rate constants in sidebar or `signaling_parameters` class to explore:
- sensitivity to parameter changes
- robustness of pathway activation
- threshold effects
- oscillatory behavior

### alternative models
replace ode system in `model_equations` to implement:
- stochastic simulations (gillespie algorithm)
- spatial models (partial differential equations)
- multi-compartment models (cytoplasm, nucleus, membrane)
- agent-based models (cellular populations)

## troubleshooting

### simulation fails to converge
- reduce time step size
- adjust initial conditions
- check parameter ranges
- increase solver tolerance

### unexpected dynamics
- verify parameter units (min⁻¹ vs sec⁻¹)
- check scf profile configuration
- ensure positive rate constants
- review biological constraints

### performance issues
- reduce number of time steps
- shorten simulation duration
- close other browser tabs
- increase available ram

## references

### key publications
1. lennartsson j, rönnstrand l. (2012) stem cell factor receptor/c-kit: from basic science to clinical implications. physiol rev. 92:1619-1649.

2. roskoski r. (2005) signaling by kit protein-tyrosine kinase—the stem cell factor receptor. biochem biophys res commun. 337:1-13.

3. broudy vc. (1997) stem cell factor and hematopoiesis. blood. 90:1345-1364.

4. ashman lk. (1999) the biology of stem cell factor and its receptor c-kit. int j biochem cell biol. 31:1037-1051.

5. rönnstrand l. (2004) signal transduction via the stem cell factor receptor/c-kit. cell mol life sci. 61:2535-2548.

### pathway databases
- kegg pathway: hsa04640 (hematopoietic cell lineage)
- reactome: r-hsa-1433557 (signaling by scf-kit)
- wikipathways: wp481 (kit receptor signaling)

## license

this simulation is provided for educational and research purposes. 

## author

biosim project - systems biology modeling suite

## version

1.0.0 (december 2025)

## acknowledgments

model development based on extensive literature review of c-kit/scf signaling mechanisms and hematopoietic regulation. kinetic parameters derived from experimental measurements in primary cells and cell lines.
