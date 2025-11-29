# BioHabs Project Structure - Path Configuration

## ✅ VERIFIED: All scripts now use correct paths

### Project Root Structure
```
BTEH/BTEH/                          # Project root (.root)
├── BioHabs/                        # Scripts directory
│   ├── data/                       # ⚠️ DATA IS HERE (not at root)
│   │   ├── clean/                  # Cleaned GPS data
│   │   │   ├── A/                  # After period
│   │   │   ├── B/                  # Before period
│   │   │   └── OG/                 # Original raw data
│   │   ├── envi/                   # Environmental layers
│   │   │   ├── A/                  # After period
│   │   │   │   ├── stack/          # Stacked rasters
│   │   │   │   │   └── stack.tif
│   │   │   │   └── *.tif           # Individual layers
│   │   │   └── B/                  # Before period
│   │   │       ├── stack/
│   │   │       │   └── stack.tif
│   │   │       └── *.tif
│   │   ├── occ/                    # Occurrence data
│   │   └── shp/                    # Shapefiles
│   ├── 01_ssf_to_rasters.R        # ✅ FIXED
│   ├── 01_a_iSSA.R                # ⚠️ NEEDS FIXING (syntax error)
│   ├── 02_dbscan_thin_degrees.R   # ✅ OK (uses config)
│   ├── 03_h2o_replicates.R        # ✅ OK (uses config)
│   ├── 04_ssdm_replicates.R       # ✅ OK (uses config)
│   ├── 05_methods_comparison.R    # ✅ OK (uses config)
│   └── 06_uncertainity.R          # ✅ FIXED
├── config.yml                      # ✅ FIXED - all paths now include BioHabs/
├── results/                        # Output directory
├── logs/                           # Log files
└── replicates_DBSCAN/             # DBSCAN replicates

```

### Path Configuration (config.yml)
```yaml
paths:
  root: "."
  envi_before: "BioHabs/data/envi/B"     # ✅ FIXED
  envi_after:  "BioHabs/data/envi/A"     # ✅ FIXED
  occ: "BioHabs/data/occ/thinned_DBSCAN" # ✅ FIXED
  clean: "BioHabs/data/clean"            # ✅ FIXED
  raw: "BioHabs/data/raw"                # ✅ FIXED
  results_h2o: "results/H2O"
  results_ssdm: "results/SSDM"
  plans: "plans"
  logs: "logs"
```

### Script Status

#### ✅ Working Scripts
1. **01_ssf_to_rasters.R** - Fixed paths, successfully ran for Run A (E3A, E4A, E5A)
2. **02_dbscan_thin_degrees.R** - Uses config paths correctly
3. **03_h2o_replicates.R** - Uses config paths correctly
4. **04_ssdm_replicates.R** - Uses config paths correctly
5. **05_methods_comparison.R** - Just completed successfully for Run A
6. **06_uncertainity.R** - Fixed to handle SSF replicates

#### ⚠️ Needs Attention
1. **01_a_iSSA.R** - Has syntax error (unexpected end of input at line 402)
   - This is the NEW script for generating iSSA replicates with proper variance
   - Needs to be completely rewritten/fixed

### How Scripts Resolve Paths

All scripts follow this pattern:
```R
# 1. Detect script location
.script <- dirname(.this_file())

# 2. Set project root (parent of BioHabs/)
.root <- normalizePath(file.path(.script, ".."), winslash = "/")
# .root = "D:/PhD/Projects/SDM_projects/BTEH/SDM_ele/BTEH/BTEH"

# 3. Load config
cfg <- read_config(file.path(.root, "config.yml"))

# 4. Use paths from config (which now include BioHabs/ prefix)
env_dir <- cfg$paths$envi_after  # "BioHabs/data/envi/A"
full_path <- file.path(.root, env_dir)
# = "D:/PhD/Projects/SDM_projects/BTEH/SDM_ele/BTEH/BTEH/BioHabs/data/envi/A"
```

### ⚠️ CRITICAL RULE
**NEVER change the project structure!**
- Data MUST stay in `BioHabs/data/`
- Scripts MUST stay in `BioHabs/`
- All paths in `config.yml` MUST include `BioHabs/` prefix for data paths
- Results go in `results/` at project root
- Logs go in `logs/` at project root

### Next Steps for iSSA Maps

To generate 6 elephant-specific iSSA maps, you need to:

1. **Fix 01_a_iSSA.R** (currently has syntax error)
   - OR use the working 01_ssf_to_rasters.R with bootstrapping

2. **Run for all 6 elephants** (E1B, E2B, E3B, E4B, E5B, E6B):
   ```bash
   Rscript BioHabs/01_ssf_to_rasters.R --run B --mode FAST
   ```

3. **Generate replicates** using 02_dbscan_thin_degrees.R if needed

### Current Status
- ✅ Run A SSF completed (E3A, E4A, E5A)
- ✅ Methods comparison completed for Run A
- ⏳ Run B SSF pending (needs stack file check)
- ❌ 01_a_iSSA.R broken (syntax error)
