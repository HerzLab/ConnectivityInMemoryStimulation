# Structural and Functional Connectivity Predict Memory Enhancement by Direct Brain Stimulation

## Overview

This repository contains the analysis pipeline used in the study **“Structural and Functional Connectivity Predict the Effects of Direct Brain Stimulation on Memory.”** The project investigates how the *network embedding* of stimulation targets—quantified using normative **structural connectivity**, **functional connectivity**, and their **structure–function congruence**—determines inter-individual variability in stimulation-induced memory enhancement.


---

## Key Concepts and Definitions

* **Closed-loop stimulation**: Electrical stimulation triggered in real time by a classifier detecting low encoding states.
* **Random stimulation**: Stimulation delivered at randomly selected encoding events, independent of brain state.
* **Stimulation-related memory change (Δ recall)**: Difference in recall performance between stimulation and non-stimulation lists.
* **Structural connectivity (SC)**: Normative diffusion tractography–based connectivity profiles estimated using the [NeMo framework](https://github.com/kjamison/nemo).
* **Functional connectivity (FC)**: Normative task-based fMRI connectivity during verbal encoding.
* **Structure–function congruence**: Spatial overlap between a stimulation site’s structural connectivity profile and a normative verbal-encoding activation network, quantified using the Dice coefficient.

---

## Analysis Pipeline

All analyses are implemented in **MATLAB** and organized into three modules that directly correspond to the main analytical components of the paper.

### S1. Structural Connectivity of Stimulation Sites

**Script**: `S1_StructuralConnectivity.m`

**Goal**
To determine whether the structural embedding of stimulation sites predicts stimulation-related memory enhancement, and whether this relationship depends on stimulation mode.

**Methods**

* Input: Normative probabilistic tractography matrices generated with the **NeMo Tool** (HCP-based diffusion data).
* Each stimulation site is modeled as a 12 mm sphere in MNI space.
* Pairwise ChaCo values are aggregated to derive a parcel-wise structural connectivity strength profile.
* Correlation analysis (Kendall’s τ) relates structural connectivity strength to Δ recall separately for:

  * Closed-loop stimulation
  * Random stimulation
* Family-wise error is controlled using permutation-based **Tmax correction**.

**Key Outputs**

* Group-level structural connectivity maps
* Statistical maps identifying regions whose connectivity predicts memory enhancement in the closed-loop condition

---

### S2. Structure–Function Congruence with the Verbal-Encoding Network

**Script**: `S2_NetworkCongruence.m`

**Goal**
To test whether stimulation sites whose structural projections align with a normative verbal-encoding network yield greater memory enhancement.

**Methods**

* A normative verbal-encoding network is derived from task-based fMRI (word > rest contrast) in healthy participants.
* Structural connectivity maps from each stimulation site are thresholded across a range of values.
* Functional activation maps are thresholded across a range of top-percentile levels.
* Binary overlap between structural and functional maps is quantified using the **Dice similarity coefficient**.
* Spearman correlations relate Dice values to Δ recall for closed-loop and random stimulation.
* Multiple comparisons across threshold grids are controlled using **FDR**.

**Key Outputs**

* Dice–memory correlation matrices
* Visualizations demonstrating robust structure–function–behavior coupling in closed-loop stimulation only

---

### S3. Functional Connectivity Profiles

**Script**: `S3_FunctionalConnectivity.m`

**Goal**
To assess whether normative functional connectivity alone predicts stimulation-related memory enhancement.

**Methods**

* Normative task-based functional connectivity is computed from healthy participants performing a verbal encoding task.
* Pearson correlations relate functional connectivity strength to Δ recall.
* Permutation-based correction is applied across all parcels.

**Key Result (Consistent with the Paper)**
Although functional connectivity shows qualitative overlap with structural effects, **no functional connectivity effects survive multiple-comparison correction**, and FC does not independently predict memory enhancement in multivariate models.

---

## Directory Structure

```
├── src/
│   └── pipeline/
│       ├── S1_StructuralConnectivity.m
│       ├── S2_NetworkCongruence.m
│       └── S3_FunctionalConnectivity.m
│
├── data/
│   ├── NeMo/                  # Normative structural connectivity matrices
│   ├── Atlases/               # Parcellations and templates (e.g., 438-region atlas)
│   ├── fMRI/                  # Normative verbal-encoding activation and FC data
│   └── subjects.xlsx          # Behavioral and stimulation metadata
│
└── outputs/
    ├── StructuralConnectivity/
    ├── NetworkCongruence/
    └── FunctionalConnectivity/
```

---

## Interpretation and Scope

This codebase is designed to reproduce the **network-level analyses** reported in the manuscript. It is not intended to provide a clinical closed-loop stimulation system, but rather to:

* Quantify how stimulation targets are embedded within large-scale memory networks
* Disentangle the roles of stimulation timing, structural connectivity, and functional connectivity
* Provide a reproducible framework for **network-guided neuromodulation research**

Normative connectivity datasets are used throughout; patient-specific connectomes are not required to run the pipeline.

---

## License

This project is released under the **MIT License**. The code may be reused and adapted for academic research with appropriate citation.

---

## Citation

If you use this code, please cite the corresponding paper:

> **[Citation to be added after publication]**

A BibTeX entry and DOI will be inserted here once the article is formally published.
