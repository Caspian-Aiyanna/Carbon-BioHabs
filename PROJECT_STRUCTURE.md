# BioHabs Project Structure - Path Configuration

### Project Root Structure
```
BTEH/BTEH/                          # Project root (.root)
├── BioHabs/                        # Scripts directory
│   ├── data/                       # INPUT DATA (not at root)
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
│   ├── 01_ssf_to_rasters.R        # FIXED
│   ├── 01_a_iSSA.R                
│   ├── 02_dbscan_thin_degrees.R   
│   ├── 03_h2o_replicates.R        
│   ├── 04_ssdm_replicates.R       
│   ├── 05_methods_comparison.R    
│   └── 06_uncertainity.R          
├── config.yml                      
├── results/                        
├── logs/                           
└── replicates_DBSCAN/             

```

### Path Configuration (config.yml)
```yaml
paths:
  root: "."
  envi_before: "BioHabs/data/envi/B"     
  envi_after:  "BioHabs/data/envi/A"     
  occ: "BioHabs/data/occ/thinned_DBSCAN" 
  clean: "BioHabs/data/clean"            
  raw: "BioHabs/data/raw"                
  results_h2o: "results/H2O"
  results_ssdm: "results/SSDM"
  plans: "plans"
  logs: "logs"
```

### Script Status

#### Working Scripts
1. **01_ssf_to_rasters.R** - Fixed paths, successfully ran for Run A (E3A, E4A, E5A)
2. **02_dbscan_thin_degrees.R** - Uses config paths correctly
3. **03_h2o_replicates.R** - Uses config paths correctly
4. **04_ssdm_replicates.R** - Uses config paths correctly
5. **05_methods_comparison.R** - Just completed successfully for Run A
6. **06_uncertainity.R** - Fixed to handle SSF replicates



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
