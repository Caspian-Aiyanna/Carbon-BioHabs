#!/usr/bin/env Rscript
# Uncertainty Decomposition & Carbon Analysis (Phases 4 & 5)
# Implements the Hybrid Multi-Species Carbon Framework as per study flowchart
# 1. Loads H2O/SSDM replicates + SSF single
# 2. Calculates Hybrid Ensemble (Inverse Variance Weighted)
# 3. Decomposes Uncertainty
# 4. Calculates Carbon Sequestration Potential (CSP)

suppressPackageStartupMessages({
  library(optparse)
  library(terra)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(readr)
  library(ggplot2)
  library(scales)
  library(viridisLite)
})

# Command-line options
option_list <- list(
  make_option(c("-e", "--elephant"), type = "character", default = NULL, help = "Elephants (e.g. E3,E4) or ALL"),
  make_option(c("-a", "--tagA"), type = "character", default = "A", help = "Before tag [default %default]"),
  make_option(c("-b", "--tagB"), type = "character", default = "B", help = "After tag [default %default]"),
  make_option(c("-r", "--root"), type = "character", default = "results", help = "Results root [default %default]"),
  make_option(c("--env_dir"), type = "character", default = "data/envi", help = "Environment dir for Biomass [default %default]"),
  make_option(c("-o", "--outdir"), type = "character", default = "paper_results/uncertainty", help = "Output root [default %default]")
)
opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$elephant)) opt$elephant <- "E3,E4,E5,E1,E2,E6"
elephants <- toupper(unlist(strsplit(opt$elephant, ",")))
if ("ALL" %in% elephants) elephants <- c("E1", "E2", "E3", "E4", "E5", "E6")

# Helpers
first_existing <- function(...) {
  cands <- c(...)
  for (p in cands) {
    if (nzchar(p) && file.exists(p)) {
      return(p)
    }
  }
  NULL
}

# Path finders
get_reps <- function(method, ele, run, root) {
  # Returns list of raster paths
  sp <- paste0(ele, run)
  if (method == "H2O") {
    d <- file.path(root, "H2O", run, sp, "replicates")
    if (!dir.exists(d)) {
      return(character(0))
    }
    # Look for prediction_*.tif in rep folders
    reps <- list.dirs(d, recursive = FALSE)
    f <- list.files(reps, pattern = paste0("prediction_", sp, ".*\\.tif$"), full.names = TRUE)
    return(f)
  } else if (method == "SSDM") {
    d <- file.path(root, "SSDM", run, sp, "replicates")
    if (!dir.exists(d)) {
      return(character(0))
    }
    reps <- list.dirs(d, recursive = FALSE)
    f <- list.files(reps, pattern = paste0("ESDM_", sp, ".*\\.tif$"), full.names = TRUE)
    return(f)
  } else if (method == "SSF") {
    # New logic for iSSA replicates
    d <- file.path(root, "iSSA", run, sp, "replicates")
    if (!dir.exists(d)) {
      return(character(0))
    }
    # Look for SSF_Species_repX.tif
    f <- list.files(d, pattern = paste0("SSF_", sp, "_rep.*\\.tif$"), full.names = TRUE)
    return(f)
  }
  character(0)
}

get_single <- function(method, ele, run, root) {
  sp <- paste0(ele, run)
  if (method == "SSF") {
    d <- file.path(root, "iSSA", run, sp)
    # Try the new iSSA single output if replicates fail, or the old one
    f <- first_existing(
      file.path(d, paste0(sp, "_SSF_rsf_0to1.tif")),
      file.path(d, paste0(sp, "_SSF_rsf.tif"))
    )
    return(f)
  }
  NULL
}

# Main Analysis
for (ele in elephants) {
  for (run in c(opt$tagA, opt$tagB)) {
    message(sprintf("\nProcessing %s Run %s...", ele, run))

    # 1. Load Models
    # --------------
    # H2O Replicates
    h2o_files <- get_reps("H2O", ele, run, opt$root)
    if (length(h2o_files) == 0) {
      message("  No H2O replicates found. Skipping.")
      next
    }
    r_h2o_reps <- rast(h2o_files)

    # Create template from first H2O raster for alignment
    template <- r_h2o_reps[[1]]

    # Calculate H2O mean and variance
    M_h2o <- mean(r_h2o_reps, na.rm = TRUE)
    V_h2o <- app(r_h2o_reps, fun = var, na.rm = TRUE)

    # SSDM Replicates
    ssdm_files <- get_reps("SSDM", ele, run, opt$root)
    if (length(ssdm_files) == 0) {
      message("  No SSDM replicates found. Skipping.")
      next
    }
    r_ssdm_reps <- rast(ssdm_files)
    r_ssdm_reps <- resample(r_ssdm_reps, template)

    # Calculate SSDM mean and variance
    M_ssdm <- mean(r_ssdm_reps, na.rm = TRUE)
    V_ssdm <- app(r_ssdm_reps, fun = var, na.rm = TRUE)

    # SSF (Replicates OR Single)
    ssf_files <- get_reps("SSF", ele, run, opt$root)

    if (length(ssf_files) > 1) {
      # Case A: We have replicates (from 01_a_iSSA.R)
      message(sprintf("  Found %d SSF replicates. Calculating true variance.", length(ssf_files)))
      r_ssf_reps <- rast(ssf_files)
      r_ssf_reps <- resample(r_ssf_reps, template)

      M_ssf <- mean(r_ssf_reps, na.rm = TRUE)
      V_ssf <- app(r_ssf_reps, fun = var, na.rm = TRUE)
    } else {
      # Case B: Fallback to single file (Old 01_ssf_to_rasters.R)
      ssf_file <- get_single("SSF", ele, run, opt$root)
      if (is.null(ssf_file)) {
        message("  No SSF file found. Skipping.")
        next
      }
      r_ssf <- rast(ssf_file)
      r_ssf <- resample(r_ssf, template)

      message("  Using single SSF file with variance proxy.")
      M_ssf <- r_ssf
      # Proxy variance: mean of H2O and SSDM variance
      V_ssf <- (V_h2o + V_ssdm) / 2
    }

    # 3. Hybrid Ensemble (Inverse Variance Weighted)
    # ----------------------------------------------
    # w_k = 1 / (sigma_k^2 + epsilon)
    eps <- 1e-6
    W_h2o <- 1 / (V_h2o + eps)
    W_ssdm <- 1 / (V_ssdm + eps)
    W_ssf <- 1 / (V_ssf + eps)

    W_sum <- W_h2o + W_ssdm + W_ssf

    R_hybrid <- (W_h2o * M_h2o + W_ssdm * M_ssdm + W_ssf * M_ssf) / W_sum
    names(R_hybrid) <- "Hybrid_Ensemble"

    # 4. Uncertainty Decomposition
    # ----------------------------
    # Total Uncertainty: sigma_total = sqrt(1 / sum(1/sigma_k^2)) = sqrt(1/W_sum)
    Sigma_total <- sqrt(1 / W_sum)
    names(Sigma_total) <- "Uncertainty_Total"

    # 5. Carbon Sequestration Potential (CSP)
    # ---------------------------------------
    # CSP = Base + (1.5 * R_hybrid)
    # Try to load Biomass (AGB)
    agb_file <- first_existing(
      "BioHabs/data/envi/raw/ORNL_Total_Carbon_2010_30m.tif",
      "BioHabs/data/envi/raw/ORNL_AGB_Carbon_2010_30m.tif",
      file.path(opt$env_dir, "raw", "ORNL_Total_Carbon_2010_30m.tif"),
      file.path(opt$env_dir, run, "AGB.tif"),
      file.path(opt$env_dir, run, "Biomass.tif"),
      file.path(opt$env_dir, "AGB.tif")
    )

    if (!is.null(agb_file)) {
      r_agb <- rast(agb_file)
      r_agb <- resample(r_agb, template)
      Base_Carbon <- r_agb # Assuming AGB is the baseline
    } else {
      message("  [WARN] AGB/Biomass layer not found. Using constant baseline (100).")
      Base_Carbon <- template * 0 + 100
    }

    Elephant_Effect <- 1.5 # MgC/yr
    CSP <- Base_Carbon + (Elephant_Effect * R_hybrid)
    names(CSP) <- "CSP"

    # 6. Save Outputs
    # ---------------
    out_dir <- file.path(opt$outdir, ele, "combined")
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

    writeRaster(R_hybrid, file.path(out_dir, sprintf("%s_%s_Hybrid_Ensemble.tif", ele, run)), overwrite = TRUE)
    writeRaster(Sigma_total, file.path(out_dir, sprintf("%s_%s_Uncertainty_Total.tif", ele, run)), overwrite = TRUE)
    writeRaster(CSP, file.path(out_dir, sprintf("%s_%s_CSP.tif", ele, run)), overwrite = TRUE)

    # Save plots
    fig_dir <- file.path(opt$outdir, ele, "figures")
    dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

    # Raster Maps
    png(file.path(fig_dir, sprintf("%s_%s_Hybrid.png", ele, run)), width = 1200, height = 1000, res = 150)
    plot(R_hybrid, main = paste(ele, run, "Hybrid Ensemble Suitability"), col = viridis(100), axes = FALSE)
    dev.off()

    png(file.path(fig_dir, sprintf("%s_%s_Uncertainty.png", ele, run)), width = 1200, height = 1000, res = 150)
    plot(Sigma_total, main = paste(ele, run, "Total Uncertainty"), col = inferno(100), axes = FALSE)
    dev.off()

    png(file.path(fig_dir, sprintf("%s_%s_CSP.png", ele, run)), width = 1200, height = 1000, res = 150)
    plot(CSP, main = paste(ele, run, "Carbon Sequestration Potential (MgC)"), col = viridis(100), axes = FALSE)
    dev.off()

    # Scatter Plots
    # Create a data frame for plotting (sample if too large)
    df_plot <- as.data.frame(c(R_hybrid, Sigma_total), na.rm = TRUE)
    colnames(df_plot) <- c("Suitability", "Uncertainty")

    if (nrow(df_plot) > 50000) df_plot <- df_plot[sample(nrow(df_plot), 50000), ]

    p1 <- ggplot(df_plot, aes(x = Suitability, y = Uncertainty)) +
      geom_point(alpha = 0.1, color = "#2c3e50") +
      geom_smooth(method = "gam", color = "red", se = FALSE) +
      theme_minimal() +
      labs(
        title = paste(ele, run, "- Uncertainty vs Suitability"),
        x = "Hybrid Suitability", y = "Total Uncertainty"
      )

    ggsave(file.path(fig_dir, sprintf("%s_%s_Scatter_Uncertainty_vs_Suitability.png", ele, run)), p1, width = 6, height = 5)

    message("  Done.")
  }
}
message("\nAll Phase 4/5 processing complete.")
