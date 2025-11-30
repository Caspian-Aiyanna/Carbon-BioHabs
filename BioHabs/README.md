# BioHabs: Integrated Habitat Suitability and Biomass Carbon Dynamics
This repository contains the R scripts for the BioHabs pipeline, which integrates Step Selection Functions (SSF), Species Distribution Models (SSDM), and H2O AutoML to assess habitat suitability and carbon sequestration potential for elephants.

## Pipeline Overview
 **6 main stages**, executed sequentially:

1.  **`01_a_iSSA.R`**:
    *    Fits Integrated Step-Selection Analysis (iSSA) models (Signer et al. 2024), separating movement and habitat selection. Projects habitat suitability (RSF) onto the landscape.
    *    Cleaned telemetry CSVs, Environmental Raster Stack.
    *    RSF Rasters (`_SSF_rsf.tif`), Model Diagnostics (VIF, CV, Partial Plots), Coefficient Tables.

2.  **`02_dbscan_thin_degrees.R`**:
    *    Performs spatial thinning of occurrence data using DBSCAN clustering to reduce spatial autocorrelation. Generates replicates for uncertainty analysis.
    *    Cleaned telemetry CSVs.
    *    Thinned occurrence CSVs (Main + 3 Replicates).

3.  **`03_h2o_replicates.R`**:
    *    Trains H2O AutoML models on the thinned occurrence replicates.
    *    Thinned occurrence CSVs, Environmental Raster Stack.
    *    Prediction Rasters (`prediction_*.tif`), Model Metrics.

4.  **`04_ssdm_replicates.R`**:
    *    Trains Stacked Species Distribution Models (SSDM) on the thinned occurrence replicates using an ensemble of algorithms (GLM, GBM, MARS, SVM, CTA).
    *    Thinned occurrence CSVs, Environmental Raster Stack.
    *    Ensemble Prediction Rasters (`ESDM_*.tif`), Variable Importance.

5.  **`05_methods_comparison.R`**:
    *    Compares the predictions from H2O, SSDM, and SSF. Calculates correlation metrics (Pearson, Spearman), RMSE, and overlap indices.
    *    Prediction Rasters from all three methods.
    *    Comparison Tables (`per_rep_metrics.csv`), Difference Maps.

6.  **`06_uncertainty.R`**:
    *    Performs the final uncertainty decomposition and Carbon Sequestration Potential (CSP) calculation.
    *    Prediction Rasters (H2O, SSDM, SSF), Biomass Baseline Raster.
    *    Hybrid Ensemble Suitability, Total Uncertainty Map, CSP Map.

## How to Run

### Prerequisites
*   R (>= 4.0)
*   Required R packages: `terra`, `sf`, `dplyr`, `h2o`, `SSDM`, `amt`, `optparse`, `ggplot2`.
*   Java (for H2O).

### Execution

run the entire pipeline using the master script. There are two execution modes:

1.  **FAST Mode** (Default):
    *   Optimized for speed and testing.
    *   Uses fewer model iterations, fewer algorithms, and fewer cores for debugging or validation.
```bash
    Rscript BioHabs/run_pipeline.R --mode FAST
```

2.  **REPRO Mode** (Reproducible):
    *   Final results.
    *   Uses all algorithms, more iterations, and deterministic seeds.
    *   maximum robustness and reproducibility.
```bash
    Rscript BioHabs/run_pipeline.R --mode REPRO
```

### Running Individual Stages
can also run individual stages if needed:

```bash
# Example: Run Stage 1 for Run A in FAST mode
Rscript BioHabs/01_a_iSSA.R --run A --mode FAST

# Example: Run Stage 6 (Uncertainty)
Rscript BioHabs/06_uncertainty.R
```

## Directory Structure
*   `BioHabs/`: Contains all R scripts.
*   `data/`: Input data (telemetry, environmental rasters).
*   `results/`: Output results for each method (H2O, SSDM, SSF).
*   `paper_results/`: Final figures and maps for publication.
*   `logs/`: Execution logs.

## Configuration
The pipeline behavior (e.g., FAST vs REPRO mode) is controlled by `config.yml` in the project root.
