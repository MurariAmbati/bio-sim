# complement cascade signaling

comprehensive simulation of the complement system: classical, alternative, and lectin pathways. models opsonization, membrane attack complex (MAC) formation, and inflammatory responses.

## features

- **three activation pathways**: classical (antibody-mediated), alternative (spontaneous), lectin (carbohydrate-recognition)
- **opsonization**: C3b deposition on pathogen surfaces
- **lysis**: MAC-mediated membrane perforation
- **inflammation**: anaphylatoxin (C3a, C5a) generation
- **regulation**: factor H, factor I, decay mechanisms
- **interactive visualization**: streamlit dashboard with real-time dynamics

## installation

```bash
pip install -r requirements.txt
```

## usage

### run simulation

```bash
streamlit run app.py
```

### programmatic use

```python
from complement_cascade import ComplementCascade

sim = ComplementCascade()
sim.add_pathogen("bacteria_1", surface_area=10.0)
sim.activate_classical_pathway(intensity=1.5)
sim.run(duration=50.0)

df = sim.get_history_dataframe()
print(f"opsonized: {df['opsonized_count'].iloc[-1]}")
print(f"lysed: {df['lysed_count'].iloc[-1]}")
```

## architecture

```
classical pathway       alternative pathway      lectin pathway
     (C1qrs)               (C3 tickover)             (MBL-MASP)
        ↓                       ↓                        ↓
      C4b2a                   C3bBb                    C4b2a
        └──────────────────────┴──────────────────────┘
                         C3 convertase
                               ↓
                         C3 → C3b + C3a
                               ↓
                    ┌──────────┴──────────┐
              opsonization          C5 convertase
              (phagocytosis)              ↓
                                    C5b-6-7-8-9n
                                         MAC
                                         ↓
                                    cell lysis
```

## model details

- **physiological concentrations**: μg/mL plasma levels
- **kinetic rates**: based on known biochemical constants
- **regulatory control**: factor H/I degradation, DAF, C4BP
- **spatial effects**: surface area-dependent deposition
- **amplification**: alternative pathway feedback loop

## outputs

- time-series concentration profiles
- pathogen opsonization/lysis status
- inflammatory response quantification
- 3D phase space trajectories
- component activity heatmaps

## slug

`complement_cascade`
