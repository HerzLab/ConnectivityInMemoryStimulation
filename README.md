# Structural and Functional Connectivity Predict Memory Enhancement by Direct Brain Stimulation

## Overview

This repository contains the analysis pipeline used in the study **“Structural and Functional Connectivity Predict the Effects of Direct Brain Stimulation on Memory.”** The project examines how the *network embedding* of stimulation targets—quantified using normative **structural connectivity**, **functional connectivity**, and their **structure–function congruence**—accounts for inter-individual variability in stimulation-induced memory enhancement.

The analyses focus on **left lateral temporal cortex (LTC)** stimulation during a verbal free-recall task, contrasting **closed-loop (state-dependent)** and **random (open-loop)** stimulation paradigms.

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

All analyses are implemented in **MATLAB** and organized into three modules that correspond directly to the major analytical components of the manuscript.

### S1. Structural Connectivity of Stimulation Sites

**Script**: `S1_StructuralConnectivity.m`

**Goal**
To determine whether the structural embedding of stimulation sites predicts stimulation-related memory enhancement, and whether this relationship depends on stimulation mode.

**Methods**

* **Input**: Normative probabilistic tractography matrices generated with the [NeMo Tool](https://kuceyeski-wcm-web.s3.amazonaws.com/upload.html) (HCP-based diffusion data).
* Each stimulation site is modeled as a 12 mm sphere in MNI space.
* Pairwise ChaCo (Change in Connectivity) values are aggregated to derive parcel-wise structural connectivity strength profiles.
* Correlation analyses (Kendall’s τ) relate structural connectivity strength to Δ recall separately for:

  * Closed-loop stimulation
  * Random stimulation
* Family-wise error is controlled using permutation-based [Tmax correction](https://github.com/mickcrosse/PERMUTOOLS).

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
* Binary overlap between structural and functional maps is quantified using the Dice similarity coefficient.
* Spearman correlations relate Dice values to Δ recall for closed-loop and random stimulation.
* Multiple comparisons across threshold grids are controlled using false discovery rate (FDR).

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
* Family-wise error is controlled using permutation-based [Tmax correction](https://github.com/mickcrosse/PERMUTOOLS).

**Key Result (Consistent with the Manuscript)**
Although functional connectivity shows qualitative overlap with structural effects, **no functional connectivity effects survive multiple-comparison correction**, and FC does not independently predict memory enhancement in multivariate models.

---

## Directory Structure

```
├── src/
│   └── pipeline/
│       ├── S1_StructuralConnectivity.m    # Structural connectivity–memory associations
│       ├── S2_NetworkCongruence.m         # Structure–function congruence analysis
│       └── S3_FunctionalConnectivity.m    # Functional connectivity–memory associations
│
├── data/
│   ├── NeMo/                       # Normative structural connectivity matrices (ChaCo outputs)
│   ├── cocommp438*                 # Parcellations and templates (438-region atlas variants)
│   ├── HC_Parcels_ROI_signal.mat   # Pre-extracted fMRI timeseries for atlas parcels (healthy controls)
│   ├── HC_electrode_signal.mat     # Pre-extracted fMRI timeseries for stimulation-site ROIs
│   ├── WL_Active-rest_HC_1mm.nii   # Normative functional activation map (Word List > Rest)
│   └── subjects.xlsx               # Behavioral, stimulation, and session-level metadata
│
└── outputs/
    ├── StructuralConnectivity/
    │   ├── Mean_StructuralConnectivity.nii        # Group-average structural connectivity strength map
    │   ├── CV_StructuralConnectivity.nii          # Inter-subject coefficient of variation (CV) map
    │   ├── StructuralConnectivity_[group]_r.nii   # SC–Δrecall correlation maps (e.g., Close, Random)
    │   └── StructuralConnectivity_[group]_r_permu.nii # Permutation-corrected significant clusters (Tmax)
    │
    ├── NetworkCongruence/
    │   └── Dice_r_[group].tiff                    # Heatmaps of Dice–Δrecall correlations
    │       # X-axis: functional network thresholds (wlLevels)
    │       # Y-axis: structural connectivity thresholds (connLevels)
    │       # Markers: * p < 0.05 uncorrected; ** FDR q < 0.05
    │
    └── FunctionalConnectivity/
        ├── Mean_Rest_FC.tiff                      # Normative resting-state FC map
        ├── Mean_WL_FC.tiff                        # Normative task-based FC (Word List)
        ├── [Condition]_[Group]_r_p0p05_uncorr.tiff# FC–Δrecall correlation maps
        └── Func_feature.mat                       # Top weighted parcels/edges from closed-loop analysis
```

---

## Interpretation and Scope

This codebase is designed to reproduce the **network-level analyses** reported in the manuscript. It is not intended to provide a clinical closed-loop stimulation system, but rather to:

* Quantify how stimulation targets are embedded within large-scale memory networks
* Disentangle the roles of stimulation timing, structural connectivity, and functional connectivity
* Provide a reproducible framework for **network-guided neuromodulation research**

---

## License

This project is released under the **MIT License**. The code may be reused and adapted for academic research with appropriate citation.

---

## Citation

If you use this code, please cite the corresponding paper:

> **[Citation to be added after publication]**
