# Structural network embedding predicts memory outcomes of closed-loop brain stimulation

## Overview

This repository contains the analysis pipeline used in the study **“Structural network embedding predicts memory outcomes of closed-loop brain stimulation.”** The project examines how the *network embedding* of stimulation targets, quantified using normative **structural connectivity**, **functional connectivity**, and **structure–function congruence**, accounts for inter-individual variability in stimulation-induced memory enhancement.

The analyses focus on **left lateral temporal cortex (LTC)** stimulation during a verbal free-recall task, contrasting **closed-loop (state-dependent)** and **random (open-loop)** stimulation paradigms.

---

## Key Concepts and Definitions

* **Closed-loop stimulation**: Electrical stimulation triggered in real time by a classifier detecting low encoding states.
* **Random stimulation**: Stimulation delivered at randomly selected encoding events, independent of brain state.
* **Stimulation-related memory change (Δ recall)**: Difference in recall performance between stimulation and non-stimulation lists.
* **Structural connectivity (SC)**: Normative diffusion tractography–based connectivity profiles estimated using the [NeMo framework](https://github.com/kjamison/nemo).
* **Functional connectivity (FC)**: Normative task-based and resting-state fMRI connectivity.
* **Structure–function congruence**: Spatial overlap between a stimulation site’s structural connectivity profile and functional/anatomical network definitions (either derived from normative verbal-encoding activation or using Cole-Anticevic CAB-NP networks), quantified using the Dice similarity coefficient.

---

## Analysis Pipeline

All analyses are implemented in **MATLAB** and organized into four modules that correspond directly to the analytical components of the study.

### S1. Structural Connectivity of Stimulation Sites

**Script**: [S1_StructuralConnectivity.m](src/pipeline/S1_StructuralConnectivity.m)

**Goal**
To determine whether the structural embedding of stimulation sites predicts stimulation-related memory enhancement, and whether this relationship depends on stimulation mode.

**Methods**
* Normative probabilistic tractography matrices generated with the NeMo Tool (HCP-based diffusion data).
* Each stimulation site is modeled as a 12 mm sphere in MNI space.
* Pairwise ChaCo (Change in Connectivity) values are aggregated to derive parcel-wise structural connectivity strength (degree) profiles.
* Correlation analyses (Pearson, Spearman, or Kendall's τ) relate structural connectivity strength to Δ recall separately for:
  * Closed-loop stimulation (Close)
  * Random stimulation (Random)
  * All stimulation groups (All)
* Family-wise error is controlled using permutation-based Tmax correction (via `permucorr2`).

**Key Outputs**
* Group-average structural connectivity maps (`Mean_StructuralConnectivity.nii`, `CV_StructuralConnectivity.nii`, `Mean_StructuralConnectivity_log.nii`).
* Correlation maps (`StructuralConnectivity_[group]_r.nii`/`tiff`).
* Permutation-corrected statistical maps and surface plots (`StructuralConnectivity_[group]_r_permu.tiff`).

---

### S2. Functional Connectivity Profiles

**Script**: [S2_FunctionalConnectivity.m](src/pipeline/S2_FunctionalConnectivity.m)

**Goal**
To assess whether normative functional connectivity alone predicts stimulation-related memory enhancement, and to identify predictive functional features.

**Methods**
* Normative task-based functional connectivity (Word List) and resting-state functional connectivity are computed from healthy participants performing the task.
* Pearson correlations relate individual patient functional connectivity profiles (from the stimulation-site ROI to all other parcels) to Δ recall.
* Feature selection is performed to identify the top predictive edges/nodes (e.g., top 10 weighted features) in the closed-loop stimulation group.

**Key Outputs**
* Normative FC maps (`Mean_Rest_FC.tiff`, `Mean_WL_FC.tiff`).
* Condition/group-specific correlation maps (`[Condition]_[group]_r.tiff`).
* Top weighted parcels/edges from closed-loop analysis saved in `Func_feature.mat`.
* Statistical analysis matrices (`FunctionalConnectivity_Stats.mat`).

---

### S3. Structure–Function Congruence with the Verbal-Encoding Network

**Script**: [S3_NetworkCongruence.m](src/pipeline/S3_NetworkCongruence.m)

**Goal**
To test whether stimulation sites whose structural projections align with a normative verbal-encoding network yield greater memory enhancement.

**Methods**
* A normative verbal-encoding network is derived from task-based fMRI (word > rest contrast) in healthy participants.
* Structural connectivity maps from each stimulation site are thresholded across a range of values.
* Functional activation maps are thresholded across a range of top-percentile levels.
* Binary overlap between structural and functional maps is quantified using the Dice similarity coefficient.
* Spearman correlations relate Dice values to Δ recall for Close, Random, and All stimulation groups.
* Multiple comparisons across threshold grids are controlled using false discovery rate (FDR).

**Key Outputs**
* Dice–memory correlation matrix heatmaps (`Dice_r_[group].tiff`) across threshold grids.
* Group-level and individual stimulation site structural/functional overlap visualizations (`WL_node_map.tiff`, `02R1195E.tiff`).

---

### S4. Structure–Function Congruence with CAB-NP Networks

**Script**: [S4_NetworkCongruence_CABNP.m](src/pipeline/S4_NetworkCongruence_CABNP.m)

**Goal**
To test if structure-function congruence based on Cole-Anticevic (CAB-NP) network definitions predicts memory outcomes of closed-loop brain stimulation.

**Methods**
* Network definitions are loaded from `Nemo2HCP360.xlsx` (Column 1: ROI ID, Column 5: Network Name).
* A verbal memory reference network is generated (top 30% of verbal encoding activation map).
* For each Cole-Anticevic network, calculates the Dice overlap between the network mask and the subjects' structural connectivity (SC) networks.
* Correlates individual subject Dice overlap scores with memory performance (`delta_recall`) using Spearman correlation for Close, Random, and All stimulation groups.
* Applies Benjamini-Hochberg FDR correction on the Close group p-values.
* Saves rankings of networks and generates scatter plots and surface visualization maps for each network.

**Key Outputs**
* Rankings of CAB-NP networks by correlation with recall (`CABNP_Network_Rankings.csv`).
* Horizontal bar plot ranking all networks with significance markers (`CABNP_R_rankings.tiff`).
* Network surface maps and scatter plots inside `surfaces/` subdirectory.

---

## Directory Structure

```
├── src/
│   ├── pipeline/
│   │   ├── S1_StructuralConnectivity.m    # Structural connectivity–memory associations
│   │   ├── S2_FunctionalConnectivity.m    # Functional connectivity–memory associations
│   │   ├── S3_NetworkCongruence.m         # Structure–function congruence analysis (with functional network)
│   │   └── S4_NetworkCongruence_CABNP.m    # Structure–function congruence analysis (with Cole-Anticevic networks)
│   │
│   └── utils/
│       ├── mean_by_atlas.m               # Computes mean activation values per atlas parcel
│       ├── mean_by_atlas_4d.m            # Computes mean activation values for 4D datasets
│       ├── parcel_to_HCPsurface.m        # Projects parcel-level values onto HCP cortical surface
│       ├── parcel_to_volume.m            # Maps parcel-level values back into MNI volume space
│       └── permucorr2.m                  # Permutation-based correlation and correction
│
├── data/
│   ├── NeMo/                       # Normative structural connectivity matrices (NeMo outputs)
│   ├── MNI_ROI/                    # Contains ROI definitions of stimulation site (12 mm sphere) in MNI space
│   ├── cocommp438*                 # Parcellations and templates (438-region atlas variants)
│   ├── HC_Parcels_ROI_signal.mat   # Pre-extracted fMRI timeseries for atlas parcels (in healthy controls)
│   ├── HC_electrode_signal.mat     # Pre-extracted fMRI timeseries for stimulation-site ROIs (in healthy controls)
│   ├── WL_Active-rest_HC_1mm.nii   # Normative functional activation map (Word List > Rest)
│   ├── Nemo2HCP360.xlsx            # ROI mapping table to Glasser 360 & CAB-NP network names
│   ├── Region_values.xlsx          # Region heatmap values used for exclusions
│   └── subjects.xlsx               # Behavioral, stimulation, and session-level metadata
│
└── outputs/
    ├── StructuralConnectivity/
    │   ├── Mean_StructuralConnectivity.nii         # Group-average structural connectivity strength map
    │   ├── CV_StructuralConnectivity.nii           # Inter-subject coefficient of variation (CV) map
    │   ├── StructuralConnectivity_[group]_r.nii    # SC–Δrecall correlation maps (e.g., Close, Random, All)
    │   ├── StructuralConnectivity_[group]_r.tiff   # Surface visualization of correlation maps
    │   └── StructuralConnectivity_[group]_r_permu.tiff # Permutation-corrected significant surface clusters
    │
    ├── FunctionalConnectivity/
    │   ├── Mean_Rest_FC.tiff                      # Normative resting-state FC map
    │   ├── Mean_WL_FC.tiff                        # Normative task-based FC (Word List)
    │   ├── [Condition]_[Group]_r.tiff             # FC–Δrecall correlation surface plots
    │   └── Func_feature.mat                       # Top weighted parcels/edges from closed-loop analysis
    │
    ├── NetworkCongruence/
    │   ├── Dice_r_[group].tiff                    # Heatmaps of Dice–Δrecall correlations
    │   │   # X-axis: Normative verbal-encoding network thresholds 
    │   │   # Y-axis: Structural connectivity thresholds
    │   └── WL_node_map*                           # Functional network maps (NIfTI & Surface plots)
    │
    └── NetworkCongruence_CABNP/
        ├── CABNP_Network_Rankings.csv             # Rankings of CAB-NP networks by correlation with recall
        ├── CABNP_R_rankings.tiff                  # Horizontal bar plot ranking all networks
        └── surfaces/                              # Individual CAB-NP network surface and scatter plots
```

---

## License

This project is released under the **MIT License**. The code may be reused and adapted for academic research with appropriate citation.

---

## Citation

If you use this code, please cite the corresponding paper:

> **Zhang Q, Ezzyat Y, Cao R, Javidi SS, Sperling MR, Kahana MJ, Tracy JI, Herz N. Structural and Functional Connectivity Predict the Effects of Direct Brain Stimulation on Memory. bioRxiv [Preprint]. 2026 Mar 19:2026.03.17.712408. doi: 10.64898/2026.03.17.712408. PMID: 41890041; PMCID: PMC13015471.**
