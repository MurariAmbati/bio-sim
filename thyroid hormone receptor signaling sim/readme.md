# thyroid hormone receptor signaling

comprehensive mechanistic model of thyroid hormone receptor signaling pathways including hormone synthesis, metabolism, receptor dynamics, and physiological effects on basal metabolism and development.

## features

- detailed molecular model of thyroid hormone receptor (tr) signaling
- t4/t3 synthesis, conversion, and degradation dynamics
- receptor binding and coregulator recruitment
- transcriptional regulation of target genes
- metabolic and developmental effects
- feedback regulation via tsh axis
- disease state simulations (hypothyroidism, hyperthyroidism, receptor mutations)
- interactive streamlit visualization interface

## installation

```bash
pip install -r requirements.txt
```

## usage

### run streamlit app

```bash
streamlit run app.py
```

### programmatic usage

```python
from thyroid_receptor import ThyroidReceptorModel, ThyroidReceptorParameters

# create model
model = ThyroidReceptorModel()

# run simulation
t, solution = model.simulate((0, 200))

# calculate metrics
metrics = model.calculate_metrics(solution)
print(f"mean t3: {metrics['mean_t3']:.2f}")
print(f"metabolic rate: {metrics['mean_metabolic_rate']:.2f}")
```

## model components

### state variables

- tsh: thyroid stimulating hormone
- t4: thyroxine (prohormone)
- t3: triiodothyronine (active hormone)
- rt3: reverse t3 (inactive metabolite)
- tr_free: unbound thyroid receptor
- tr_t3: hormone-bound receptor
- tr_t3_coreg: activated receptor-coactivator complex
- tr_corep: receptor-corepressor complex
- target_mrna: target gene transcripts
- metabolic_enzyme: metabolic proteins
- developmental_factor: developmental proteins
- metabolic_rate: cellular metabolism
- d2: type 2 deiodinase
- d3: type 3 deiodinase

### key processes

- hormone synthesis and secretion
- peripheral deiodination (t4→t3 conversion)
- receptor binding and activation
- coregulator recruitment
- gene transcription
- metabolic regulation
- developmental programs
- negative feedback via hpt axis

## simulation modes

- normal physiology
- hypothyroidism
- hyperthyroidism
- receptor mutations
- development
- parameter exploration

## references

based on established models of thyroid hormone signaling and receptor dynamics in mammalian systems.
