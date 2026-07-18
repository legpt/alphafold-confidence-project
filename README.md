# AlphaFold Confidence Score Analysis — Human p53 (TP53)

Fetches AlphaFold's precomputed prediction for human p53 (UniProt `P04637`) from
the [AlphaFold Protein Structure Database](https://alphafold.ebi.ac.uk/) and
analyzes/visualizes its confidence metrics:

- **pLDDT** — per-residue confidence (0-100)
- **PAE** — predicted aligned error between residue pairs

## Contents

- `alphafold_confidence_analysis.ipynb` — main notebook: fetches data via the
  AlphaFold DB API, computes summary statistics and confidence intervals on
  pLDDT (overall and per structural domain), and visualizes:
  - per-residue pLDDT along the sequence with AlphaFold's confidence bands
  - a PAE heatmap
  - the 3D structure colored by pLDDT

## Setup

```
pip install -r requirements.txt
jupyter notebook alphafold_confidence_analysis.ipynb
```

## Findings

p53's folded DNA-binding domain (residues 94-312) scores high and tight on
pLDDT, while its intrinsically disordered N-/C-terminal regions score low and
noisy — which is the *correct* signal, not a modeling failure. The PAE heatmap
shows a block structure: AlphaFold is confident about each domain's internal
fold but not about how domains are positioned relative to one another, since
they're joined by flexible linkers.
