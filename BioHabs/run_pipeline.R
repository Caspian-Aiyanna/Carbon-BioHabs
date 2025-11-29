#!/usr/bin/env Rscript
# ============================================================================
# BTEH Pipeline - Master Execution Script
# ============================================================================
# This script runs the entire BioHabs pipeline end-to-end.
# It orchestrates the execution of the 6 main analysis scripts.
#
# Stages:
#   1. SSF to Rasters (Movement-based habitat selection)
#   2. DBSCAN Thinning (Spatial thinning with replicates)
#   3. H2O Replicates (AutoML with uncertainty)
#   4. SSDM Replicates (Ensemble SDM with uncertainty)
#   5. Methods Comparison (H2O vs SSDM vs SSF)
#   6. Uncertainty Analysis & Carbon Calculation
# ============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(optparse)
})

# Command-line options
option_list <- list(
  make_option(c("--mode"),
    type = "character", default = "FAST",
    help = "Execution mode: FAST (default) or REPRO [for publication]"
  )
)
opt <- parse_args(OptionParser(option_list = option_list))

# Configuration
ROOT <- normalizePath(".", winslash = "/")
SCRIPT_DIR <- "BioHabs" # Scripts are in this subdirectory
LOG_DIR <- file.path(ROOT, "logs")
dir.create(LOG_DIR, recursive = TRUE, showWarnings = FALSE)

MASTER_LOG <- file.path(LOG_DIR, "master_pipeline.log")
MODE <- toupper(opt$mode) # FAST or REPRO
RUNS <- c("A", "B")

# ============================================================================
# Logging Helper
# ============================================================================
log_msg <- function(..., level = "INFO") {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  msg <- paste0("[", timestamp, "] [", level, "] ", paste0(..., collapse = ""))
  cat(msg, "\n")
  cat(msg, "\n", file = MASTER_LOG, append = TRUE)
}

log_section <- function(title) {
  sep <- paste(rep("=", 80), collapse = "")
  log_msg(sep)
  log_msg(title)
  log_msg(sep)
}

# ============================================================================
# Execute Script Helper
# ============================================================================
run_script <- function(script_name, args = list(), required = TRUE) {
  log_section(paste("RUNNING:", script_name))

  # Look for script in SCRIPT_DIR
  script_path <- file.path(ROOT, SCRIPT_DIR, script_name)
  if (!file.exists(script_path)) {
    msg <- paste("Script not found:", script_path)
    if (required) {
      log_msg(msg, level = "ERROR")
      stop(msg)
    } else {
      log_msg(msg, level = "WARN")
      return(list(success = FALSE, skipped = TRUE))
    }
  }

  # Build command
  cmd_args <- c(script_path)
  for (arg_name in names(args)) {
    cmd_args <- c(cmd_args, paste0("--", arg_name, "=", args[[arg_name]]))
  }

  cmd <- paste("Rscript", paste(shQuote(cmd_args), collapse = " "))
  log_msg("Command: ", cmd)

  # Run
  start_time <- Sys.time()
  result <- system(cmd, intern = FALSE, wait = TRUE)
  end_time <- Sys.time()
  duration <- as.numeric(difftime(end_time, start_time, units = "secs"))

  if (result == 0) {
    log_msg("SUCCESS (", round(duration, 2), " seconds)", level = "SUCCESS")
    return(list(success = TRUE, duration = duration))
  } else {
    msg <- paste("FAILED with exit code:", result)
    if (required) {
      log_msg(msg, level = "ERROR")
      stop(msg)
    } else {
      log_msg(msg, level = "WARN")
      return(list(success = FALSE, duration = duration))
    }
  }
}

# ============================================================================
# Main Pipeline
# ============================================================================
main <- function() {
  log_section("BTEH PIPELINE - MASTER EXECUTION")
  log_msg("Root directory: ", ROOT)
  log_msg("Script directory: ", file.path(ROOT, SCRIPT_DIR))
  log_msg("Mode: ", MODE)
  log_msg("Runs: ", paste(RUNS, collapse = ", "))
  log_msg("Start time: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))

  pipeline_start <- Sys.time()

  # Run pipeline for each run (A and B)
  for (run in RUNS) {
    log_section(paste("PROCESSING RUN:", run))

    # Stage 1: iSSA to Rasters
    log_section(paste("STAGE 1: iSSA TO RASTERS - Run", run))
    run_script("01_a_iSSA.R",
      args = list(run = run, mode = MODE),
      required = TRUE
    )

    # Stage 2: DBSCAN Thinning
    log_section(paste("STAGE 2: DBSCAN THINNING - Run", run))
    run_script("02_dbscan_thin_degrees.R",
      args = list(run = run),
      required = TRUE
    )

    # Stage 3: H2O Replicates
    log_section(paste("STAGE 3: H2O AUTOML REPLICATES - Run", run))
    run_script("03_h2o_replicates.R",
      args = list(run = run, mode = MODE),
      required = TRUE
    )

    # Stage 4: SSDM Replicates
    log_section(paste("STAGE 4: SSDM REPLICATES - Run", run))
    run_script("04_ssdm_replicates.R",
      args = list(run = run, mode = MODE),
      required = TRUE
    )

    # Stage 5: Methods Comparison
    log_section(paste("STAGE 5: METHODS COMPARISON - Run", run))
    run_script("05_methods_comparison.R",
      args = list(run = run, mode = MODE),
      required = FALSE
    )
  }

  # Stage 6: Uncertainty Analysis (combines both runs)
  log_section("STAGE 6: UNCERTAINTY ANALYSIS & CARBON CALCULATION")
  run_script("06_uncertainty.R",
    args = list(
      tagA = "A",
      tagB = "B",
      root = "results" # Relative to project root
    ),
    required = FALSE
  )

  # ============================================================================
  # Summary
  # ============================================================================
  pipeline_end <- Sys.time()
  total_duration <- as.numeric(difftime(pipeline_end, pipeline_start, units = "mins"))

  log_section("PIPELINE COMPLETE")
  log_msg("Total duration: ", round(total_duration, 2), " minutes")
  log_msg("End time: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
  log_msg("Check logs in: ", LOG_DIR)
}

# Run the pipeline
tryCatch(
  {
    main()
  },
  error = function(e) {
    log_msg("PIPELINE FAILED: ", conditionMessage(e), level = "ERROR")
    quit(status = 1)
  }
)
