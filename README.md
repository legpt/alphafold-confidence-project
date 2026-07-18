# AlphaFold Confidence Score Analysis — Human p53 (TP53)

Fetches AlphaFold's precomputed prediction for human p53 (UniProt `P04637`) from
the [AlphaFold Protein Structure Database](https://alphafold.ebi.ac.uk/) and
analyzes/visualizes its confidence metrics:

- **pLDDT** — per-residue confidence (0-100)
- **PAE** — predicted aligned error between residue pairs

## Research question

**Does AlphaFold's per-residue confidence score (pLDDT) track real, experimentally-verified
protein disorder — and does that relationship hold generally, beyond a single protein?**

- **H1 (single-protein):** within p53, residues in experimentally-characterized disordered
  regions have significantly lower pLDDT than residues in the folded domains.
- **H2 (generalization):** across a panel of proteins with independently-verified disorder
  (from [DisProt](https://disprot.org/)) vs. proteins with no documented disorder, low pLDDT
  predicts disorder status better than chance.

Both hypotheses are tested against DisProt as an independent ground truth, not just inferred
from visual inspection of one protein.

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
