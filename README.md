# Connectivity Analysis for Memory Recall

## Project Overview

This project investigates the relationship between brain connectivity (Structural and Functional) and memory recall performance. Guided by the draft paper *Structural and Functional Connectivity Predict the Effects of Direct Brain Stimulation on Memory*, the analysis focuses on how the congruence between structural networks and functional activation maps relates to memory outcomes, specifically under different stimulation conditions (Close vs. Random).

## Pipeline Description

The analysis pipeline is implemented in MATLAB and organized into three main stages. The scripts are located in `src/pipeline/`.

### 1. Structural Connectivity (S1)

**Script**: `S1_StructuralConnectivity.m`

* **Purpose**: To characterize the underlying structural scaffold of the brain and correlate nodal structural strength with memory performance.
* **Methods**:
  * **Input**: Individual probabilistic tractography matrices (from `data/NeMo/`).
  * **Analysis**: Calculates the Mean Structural Connectivity (SC) and Coefficient of Variation (CV) across subjects. It then performs a permutation-based correlation (Kendall's tau) between nodal degree and `delta_recall` scores for 'All', 'Close', and 'Random' groups.
  * **Key Output**: `Mean_StructuralConnectivity.nii` (Average SC map), `StructuralConnectivity_close_r.nii` (Correlation map for Close group).

### 2. Network Congruence (S2)

**Script**: `S2_NetworkCongruence.m`

* **Purpose**: To quantify the spatial "match" (congruence) between the functional task activation network and the structural connectivity network. The hypothesis is that higher congruence predicts better recall.
* **Methods**:
  * **Input**: Word List (WL) activation map (`WL_Active-rest_HC_1mm.nii`).
  * **Analysis**:
    * Define **WL Network**: Top $K$ activated regions.
    * Define **SC Network**: Regions structurally connected to the WL Network above a strength threshold.
    * **Metric**: Computes the **Dice Coefficient** between WL and SC networks across a range of thresholds.
    * **Statistics**: Correlates Dice coefficients with `delta_recall` using Spearman correlation, applying FDR correction for multiple comparisons.
  * **Key Output**: `Dice_r_close.tiff` (Matrix plot of correlations between Congruence and Recall).

### 3. Functional Connectivity (S3)

**Script**: `S3_FunctionalConnectivity.m`

* **Purpose**: To likely functional features that drive memory performance using Resting-state and Task-based connectivity.
* **Methods**:
  * **Input**: Pre-computed electrode and atlas ROI signals.
  * **Analysis**:
    * Computes Functional Connectivity (FC) matrices (Pearson correlation) for 'Rest' and 'WL' (Task) conditions.
    * Averages FC maps across healthy controls to create normative templates.
    * Correlates individual patient FC patterns with `delta_recall`.
    * **Feature Selection**: Identifies top predictive features (edges/nodes) in the 'Close' group for further analysis.
  * **Key Output**: `Mean_Rest_FC.tiff`, `Func_feature.mat` (Top 10 features).

---

## Data Directory Structure (`data/`)

The `data/` folder contains all necessary inputs for the pipeline.

* **NeMo/**: Contains pre-processed structural connectivity matrices for each subject (`*ifod2act_chacoconn_cocommp438subj_mean.csv`).
* **MNI_ROI/**: Contains ROI definitions and mask files used for signal extraction.
* **Atlases/Templates**:
  * `cocommp438_dil*.nii`: Parcellation atlases (dilated versions) used for mapping data to brain regions.
  * `WL_Active-rest_HC_1mm.nii`: Average functional activation map for the Word List task.
* **Clinical Data**:
  * `subjects.xlsx`: Master spreadsheet containing patient IDs, group assignments (`randomstim`: 0=Close, 1=Random), and behavioral scores (`delta_recall`).
  * `Nemo2HCP360.xlsx`: Mapping file for visualization on HCP surfaces.
* **Signals**:
  * `HC_Parcels_ROI_signal.mat`: Pre-extracted timeseries for atlas parcels.
  * `HC_electrode_signal.mat`: Pre-extracted timeseries for electrodes.

## Outputs and Results (`outputs/`)

Results are organized by analysis module.

### `outputs/StructuralConnectivity/`

* `Mean_StructuralConnectivity.nii`: The group-average structural connectivity strength map.
* `CV_StructuralConnectivity.nii`: Coefficient of Variation map showing inter-subject variability.
* `StructuralConnectivity_[group]_r.nii`: Map of correlations between SC degree and recall for a specific group.
* `Node_[group]_r_permu.nii`: Significant clusters after permutation testing.

### `outputs/NetworkCongruence/`

* `WL_node_map.nii`: The defined functional network derived from the activation map.
* `Dice_r_[group].tiff`: Heatmaps showing the correlation (R-value) between Network Congruence (Dice) and Recall.
  * **X-axis**: Threshold for defining the Functional Network (`wlLevels`).
  * **Y-axis**: Threshold for defining the Structural Network (`connLevels`).
  * **Markers**: `*` (p<0.05 uncorrected), `**` (FDR significant).

### `outputs/FunctionalConnectivity/`

* `Mean_Rest_FC.tiff` / `Mean_WL_FC.tiff`: Average functional connectivity maps for Rest and Task conditions.
* `[Condition]_[Group]_r.tiff`: Correlation surface maps showing relationships between FC and Recall.
* `Func_feature.mat`: A MATLAB data file containing the top 10 weighted features (parcels/edges) identified from the 'Close' group analysis, used for downstream processing.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
Both the source code (`src/`) and data outputs (`outputs/`, `data/`) are available for use under this license.
