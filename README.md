# Analysis on Acne Enzymes

In this report, I have conducted research using Alphafold- an AI protein structure tool- to determine how acne is related to inflammasomes in the skin. The specifics of the analysis I completed can be found in the acne folder. I became interested in this topic after I realized I was experiencing acne breakouts at 15 and have been curious about the underlying biological mechanisms ever since. 

## Goal of the Study

The goal of this report is to take a simple approach with alphafold to determine how certain enzymes relate to a target enzyme (Caspase-1). The reason for choosing Caspase-1 is because this enzyme plays a key role of inflammation in acne. The enzymes I chose to compare Caspase-1 to are MMP1, ELANE, PTGS2, RIPK1, and TNFAIP3. Along with Caspase-1, these five enzymes were selected because they represent each key stage of acne, from initial inflammation and tissue breakdown to the body's natural healing response. The goal is to use simple statistical methods to visualize and analyze how these enzymes relate to the target enzyme Caspase-1 structurally and biochemically.

## Data Extraction

The data fetched for these six enzymes comes from Alphafold database. This database contains the following features for each enzyme:
- pLDDT: scores how confident AlphaFold is in the 3D shape of each amino acid
- Mean PAE: Estimates how accurately two different parts of the protein are positioned together
- C-alpha Coordinates: The 3D coordinates marking where each amino acid sits in physical space
- Distance to Centroid: How far an amino acid is from the center of the protein.
- SASA: Measures how much of an amino acid touches the surrounding fluid
- Hydropathy: Shows whether an amino acid prefers water or hides inside oily pockets
- Charge: The electrical charge carried by an amino acid inside human skin cells

To parse and extract data from the database, I used python and python libraries such as Pandas, NumPy, and Biopython. 

## Summary Statistics

### 1. Biophysical Summary Statistics Across Target Enzymes
```text
Comprehensive Biophysical Summary Statistics:
         residue_count  plddt_mean  plddt_std  plddt_median  plddt_iqr  plddt_skew  sasa_mean  sasa_std  sasa_median  sasa_iqr  pae_mean  pae_std
protein                                                                                                                                          
CASP1              404       81.68      22.61         91.38      17.16       -1.48      61.02     52.79        53.09     88.51     17.28     6.41
ELANE              267       88.20      17.96         97.50       8.34       -1.72      57.79     58.89        39.27     80.85      9.87     7.90
MMP1               469       91.22      14.96         97.00       4.57       -2.53      54.10     52.41        37.57     74.12      7.72     6.76
PTGS2              604       93.02      15.28         98.00       1.89       -3.14      48.14     48.11        32.94     68.58      6.47     6.66
RIPK1              671       69.74      25.29         80.38      51.03       -0.44      81.97     60.28        79.64     95.70     23.26     4.96
TNFAIP3            790       73.73      23.49         81.97      41.40       -0.66      77.74     58.35        77.30     93.84     23.42     5.56
```
**Significance:** Catalytic effector enzymes (PTGS2, MMP1, and ELANE) show high structural confidence (mean pLDDT > 88) and compact cores, whereas Caspase-1 exhibits intermediate stability, and signaling regulators (RIPK1 and TNFAIP3) display higher solvent exposure and flexibility.

---

### 2. Structural Confidence Tier Distribution
```text
Structural Confidence Tier Distribution (% per enzyme):
confidence_tier  Very High (>=90)  Confident (70-90)  Low (50-70)  Very Low (<50)
protein                                                                          
CASP1                        56.7               21.3          8.2            13.9
ELANE                        75.7                7.5         10.9             6.0
MMP1                         85.1                4.3          5.8             4.9
PTGS2                        90.7                0.8          3.1             5.3
RIPK1                        32.2               26.5         10.4            30.8
TNFAIP3                      35.4               30.1         10.1            24.3
```
**Significance:** Over 75% to 90% of residues in PTGS2, MMP1, and ELANE are modeled in the highest confidence tier, whereas RIPK1 and TNFAIP3 contain roughly 25% to 30% low-confidence/disordered regions, and Caspase-1 retains a solidly confident core (78% confident or higher) alongside a flexible pro-domain.

---

### 3. Residue Chemical Property Distribution
```text
Residue Chemical Property Distribution (% per enzyme):
chemical_class  Acidic (-)  Basic (+)  Polar Neutral  Nonpolar Hydrophobic
protein                                                                   
CASP1                 13.6       13.6           27.5                  45.3
ELANE                  4.9       11.6           24.7                  58.8
MMP1                  12.4       14.9           24.7                  48.0
PTGS2                 10.3       13.2           28.0                  48.5
RIPK1                 13.0       14.2           29.2                  43.7
TNFAIP3               11.0       16.7           30.0                  42.3
```
**Significance:** All six enzymes maintain a conserved hydrophobic structural core (42%–59% nonpolar residues), with balanced acidic and basic surface residues enabling specific protein-protein and catalytic interactions in inflammatory tissue.

---

### Visualizations

![Biophysical Summary Heatmap](images/acne_plot_1.png)
*Figure 1: Heatmap showing standardized (Z-score) biophysical metrics (pLDDT, SASA, PAE, hydropathy, and charge) across each of the six enzymes.*

![pLDDT and SASA Density Distributions](images/acne_plot_2.png)
*Figure 2: Kernel density curves illustrating the contrast in AlphaFold structural confidence (pLDDT) and solvent accessibility (SASA) between rigid enzymes and flexible signaling regulators.*

![Residue Chemical Property Bar Chart](images/acne_plot_4.png)
*Figure 3: Grouped bar chart comparing the percentage breakdown of acidic, basic, polar neutral, and nonpolar hydrophobic amino acids across the six targets.*

![3D Multi-Dimensional Residue Folding Space](images/acne_plot_6.png)
*Figure 4: 3D scatter visualization showing how residues organize in physical space by centroid distance, solvent accessibility, and AlphaFold confidence colored by hydrophobicity.*

### Key Observations
- **High Confidence in Rigid Enzymes**: PTGS2 (COX-2), MMP1, and ELANE show very high structural confidence (mean pLDDT > 88), reflecting stable, well-folded catalytic domains with low predicted error (PAE < 10 Å).
- **Target Enzyme (Caspase-1)**: Caspase-1 sits in an intermediate-to-high confidence tier (mean pLDDT 81.68, median 91.38). While its core catalytic domain is modeled with high confidence, its flexible pro-domain and regulatory loop regions show lower local confidence.
- **Flexible Regulatory Proteins**: Signaling and regulatory enzymes (RIPK1 and TNFAIP3) have lower average confidence (mean pLDDT ~70–74) and higher surface accessibility (SASA > 77 Å²), indicating extensive flexible linkers and disordered regulatory regions.