# BioHabs: Integrated Habitat Suitability and Biomass Carbon Dynamics

This repository contains the R scripts for the BioHabs pipeline, which integrates Step Selection Functions (SSF), Species Distribution Models (SSDM), and H2O AutoML to assess habitat suitability and carbon sequestration potential for elephants.

## Pipeline Overview

The pipeline consists of **6 main stages**, executed sequentially:

1.  **`01_ssf_to_rasters.R`**:
    *   **Purpose**: Fits Step Selection Functions (SSF) to telemetry data and projects relative selection strength (RSF) onto the landscape.
    *   **Input**: Cleaned telemetry CSVs, Environmental Raster Stack.
    *   **Output**: RSF Rasters (`_SSF_rsf.tif`).

2.  **`02_dbscan_thin_degrees.R`**:
    *   **Purpose**: Performs spatial thinning of occurrence data using DBSCAN clustering to reduce spatial autocorrelation. Generates replicates for uncertainty analysis.
    *   **Input**: Cleaned telemetry CSVs.
    *   **Output**: Thinned occurrence CSVs (Main + 3 Replicates).

3.  **`03_h2o_replicates.R`**:
    *   **Purpose**: Trains H2O AutoML models on the thinned occurrence replicates.
    *   **Input**: Thinned occurrence CSVs, Environmental Raster Stack.
    *   **Output**: Prediction Rasters (`prediction_*.tif`), Model Metrics.

4.  **`04_ssdm_replicates.R`**:
    *   **Purpose**: Trains Stacked Species Distribution Models (SSDM) on the thinned occurrence replicates using an ensemble of algorithms (GLM, GBM, MARS, SVM, CTA).
    *   **Input**: Thinned occurrence CSVs, Environmental Raster Stack.
    *   **Output**: Ensemble Prediction Rasters (`ESDM_*.tif`), Variable Importance.

5.  **`05_methods_comparison.R`**:
    *   **Purpose**: Compares the predictions from H2O, SSDM, and SSF. Calculates correlation metrics (Pearson, Spearman), RMSE, and overlap indices.
    *   **Input**: Prediction Rasters from all three methods.
    *   **Output**: Comparison Tables (`per_rep_metrics.csv`), Difference Maps.

6.  **`06_uncertainity.R`**:
    *   **Purpose**: Performs the final uncertainty decomposition and Carbon Sequestration Potential (CSP) calculation.
    *   **Input**: Prediction Rasters (H2O, SSDM, SSF), Biomass Baseline Raster.
    *   **Output**: Hybrid Ensemble Suitability, Total Uncertainty Map, CSP Map.

## How to Run

### Prerequisites
*   R (>= 4.0)
*   Required R packages: `terra`, `sf`, `dplyr`, `h2o`, `SSDM`, `amt`, `optparse`, `ggplot2`.
*   Java (for H2O).

### Execution

You can run the entire pipeline using the master script. There are two execution modes:

1.  **FAST Mode** (Default):
    *   Optimized for speed and testing.
    *   Uses fewer model iterations, fewer algorithms, and fewer cores.
    *   Ideal for debugging or quick validation.
```bash
    Rscript BioHabs/run_pipeline.R --mode FAST
```

2.  **REPRO Mode** (Reproducible):
    *   Full production run for publication.
    *   Uses all algorithms, more iterations, and deterministic seeds.
    *   Ensures maximum robustness and reproducibility.
```bash
    Rscript BioHabs/run_pipeline.R --mode REPRO
```

### Running Individual Stages
You can also run individual stages if needed:

```bash
# Example: Run Stage 1 for Run A in FAST mode
Rscript BioHabs/01_ssf_to_rasters.R --run A --mode FAST

# Example: Run Stage 6 (Uncertainty)
Rscript BioHabs/06_uncertainity.R
```

## Directory Structure

*   `BioHabs/`: Contains all R scripts.
*   `data/`: Input data (telemetry, environmental rasters).
*   `results/`: Output results for each method (H2O, SSDM, SSF).
*   `paper_results/`: Final figures and maps for publication.
*   `logs/`: Execution logs.

## Configuration
The pipeline behavior (e.g., FAST vs REPRO mode) is controlled by `config.yml` in the project root.
