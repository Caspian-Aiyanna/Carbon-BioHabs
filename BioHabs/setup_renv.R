# ==============================================================================
# BioHabs Environment Setup (renv)
# ==============================================================================
# Run this script ONCE to initialize the reproducible environment.
# It will:
# 1. Install 'renv' if missing
# 2. Initialize the project
# 3. Discover dependencies
# 4. Create 'renv.lock'
# ==============================================================================

if (!requireNamespace("renv", quietly = TRUE)) {
    message("Installing renv...")
    install.packages("renv")
}

# Initialize renv
message("Initializing renv project...")
# We use 'bare = TRUE' to avoid auto-installing everything immediately,
# then hydrate to find what's actually used.
renv::init(bare = TRUE, restart = FALSE)

# Hydrate: Discover packages used in .R files and install them into the library
message("Hydrating dependencies (scanning scripts)...")
renv::hydrate()

# Snapshot: Record the state into renv.lock
message("Creating lockfile (renv.lock)...")
renv::snapshot(prompt = FALSE)

message("\n==============================================================================")
message("Setup Complete! 'renv.lock' has been created.")
message("To reproduce this environment on HPC, run: Rscript -e 'renv::restore()'")
message("==============================================================================")
