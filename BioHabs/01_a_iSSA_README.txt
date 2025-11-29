# 01_a_iSSA.R - SSF Analysis
# Step-Selection Function (SSF) habitat suitability maps.

### Quick Start

```bash
# Run for all elephants in Run B (FAST mode)
Rscript BioHabs/01_a_iSSA.R --run B --mode FAST

# Run for specific elephant
Rscript BioHabs/01_a_iSSA.R --run B --mode FAST --species E5
# Run for Run A
Rscript BioHabs/01_a_iSSA.R --run A --mode FAST
```

### Options
- `--run`: A or B (default: B)
- `--mode`: FAST or REPRO (default: FAST)
- `--species`: Specific elephant (optional, e.g., E5B)

###Function
1. Loads GPS tracking data (already at 30-minute intervals)
2. Respects original data structure (NO resampling)
3. Generates steps from consecutive locations
4. Fits Step-Selection Function models
5. Projects habitat suitability to raster maps
6. Saves results in `results/SSF/{run}/{species}/`

### Outputs
For each elephant:
- `{Species}_SSF_rsf.tif` - Raw suitability map
- `{Species}_SSF_rsf_0to1.tif` - Normalized (0-1) suitability map
- Diagnostic plots
- Metadata CSV

### Notes
- Only removes GPS errors (< 0.1m steps)
- Preserves all real movement data


