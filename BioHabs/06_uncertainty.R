#!/usr/bin/env Rscript
# Uncertainty Decomposition & Carbon Analysis (Phases 4 & 5)
# Implements the Hybrid Multi-Species Carbon Framework as per study flowchart
# 1. Loads H2O/SSDM replicates + SSF single
# 2. Calculates Hybrid Ensemble (Inverse Variance Weighted)
# 3. Decomposes Uncertainty


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
  make_option(c("-o", "--outdir"), type = "character", default = "results/uncertainty", help = "Output root [default %default]")
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
  # Initialize list to store plot data for this elephant across runs
  ele_plot_data <- list()

  for (run in c(opt$tagA, opt$tagB)) {
    # Explicitly clear variables to prevent leakage
    suppressWarnings(rm(M_h2o, V_h2o, M_ssdm, V_ssdm, M_ssf, V_ssf, W_h2o, W_ssdm, W_ssf))

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
    # Calculate Normalized Weights
    W_norm_h2o <- W_h2o / W_sum
    W_norm_ssdm <- W_ssdm / W_sum
    W_norm_ssf <- W_ssf / W_sum

    # Within-Model Variance (Weighted Average of Variances)
    Var_within <- (W_norm_h2o * V_h2o) + (W_norm_ssdm * V_ssdm) + (W_norm_ssf * V_ssf)

    # Between-Model Variance (H2O vs SSDM ONLY)
    W_sum_sdm <- W_h2o + W_ssdm
    R_sdm <- (W_h2o * M_h2o + W_ssdm * M_ssdm) / W_sum_sdm

    W_n_h2o_sdm <- W_h2o / W_sum_sdm
    W_n_ssdm_sdm <- W_ssdm / W_sum_sdm

    Var_between_sdm <- (W_n_h2o_sdm * (M_h2o - R_sdm)^2) +
      (W_n_ssdm_sdm * (M_ssdm - R_sdm)^2)

    Share_SDM <- W_sum_sdm / W_sum
    Var_between_contribution <- Var_between_sdm * Share_SDM

    # Total Uncertainty
    Sigma_total <- sqrt(Var_within + Var_between_contribution)
    names(Sigma_total) <- "Uncertainty_Total"

    # 5. Carbon Sequestration Potential (CSP) - REMOVED
    # ---------------------------------------
    # Concept removed as per user request.


    # 6. Save Outputs
    # ---------------
    out_dir <- file.path(opt$outdir, ele, "combined")
    dir.create(out_dir, recursive = TRUE, showWarnings = TRUE)

    writeRaster(R_hybrid, file.path(out_dir, sprintf("%s_%s_Hybrid_Ensemble.tif", ele, run)), overwrite = TRUE)
    writeRaster(Sigma_total, file.path(out_dir, sprintf("%s_%s_Uncertainty_Total.tif", ele, run)), overwrite = TRUE)

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

    # Combined Panel Plot (Hybrid | Uncertainty)
    # ------------------------------------------
    # Create a summary folder for all runs
    summary_dir <- file.path(opt$outdir, "all_run_panels")
    dir.create(summary_dir, recursive = TRUE, showWarnings = FALSE)

    panel_path <- file.path(fig_dir, sprintf("%s_%s_Panel.png", ele, run))
    summary_path <- file.path(summary_dir, sprintf("%s_%s_Panel.png", ele, run))

    png(panel_path, width = 2000, height = 1000, res = 150)
    par(mfrow = c(1, 2), mar = c(1, 1, 3, 4))
    plot(R_hybrid, main = paste(ele, run, "Hybrid Suitability"), col = viridis(100), axes = FALSE)
    plot(Sigma_total, main = paste(ele, run, "Uncertainty"), col = inferno(100), axes = FALSE)
    par(mfrow = c(1, 1))
    dev.off()

    # Copy to summary folder
    file.copy(panel_path, summary_path, overwrite = TRUE)


    # Collect data for advanced plotting
    # Stack: Hybrid, Uncertainty, H2O, SSDM, SSF
    stack_plot <- c(R_hybrid, Sigma_total, M_h2o, M_ssdm, M_ssf)
    names(stack_plot) <- c("Hybrid", "Uncertainty", "H2O", "SSDM", "SSF")

    df_plot <- as.data.frame(stack_plot, na.rm = TRUE)
    if (nrow(df_plot) > 50000) df_plot <- df_plot[sample(nrow(df_plot), 50000), ]
    df_plot$Run <- run

    # Store in list for this elephant
    ele_plot_data[[run]] <- df_plot

    message("  Run ", run, " processed.")
  } # End Run Loop

  # ====================================================
  # Generate Comparative Plots (A vs B) for this Elephant
  # ====================================================
  if (length(ele_plot_data) > 0) {
    all_data <- do.call(rbind, ele_plot_data)

    # Ensure Run is a factor for proper ordering/coloring
    all_data$Run <- factor(all_data$Run, levels = c(opt$tagA, opt$tagB))

    fig_dir <- file.path(opt$outdir, ele, "figures")
    dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

    # ---------------------------------------------------------
    # 1. Pairwise Method Comparisons (Scatter Plots)
    # ---------------------------------------------------------
    # H2O vs SSDM
    p_h2o_ssdm <- ggplot(all_data, aes(x = H2O, y = SSDM, color = Run)) +
      geom_point(alpha = 0.1, size = 0.5) +
      geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray") +
      scale_color_manual(values = c("A" = "blue", "B" = "#FF7F7F")) +
      theme_minimal() +
      labs(title = paste(ele, "- H2O vs SSDM Suitability"), x = "H2O", y = "SSDM")
    ggsave(file.path(fig_dir, sprintf("%s_Compare_Scatter_H2O_SSDM.png", ele)), p_h2o_ssdm, width = 6, height = 6)

    # H2O vs SSF
    p_h2o_ssf <- ggplot(all_data, aes(x = H2O, y = SSF, color = Run)) +
      geom_point(alpha = 0.1, size = 0.5) +
      geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray") +
      scale_color_manual(values = c("A" = "blue", "B" = "#FF7F7F")) +
      theme_minimal() +
      labs(title = paste(ele, "- H2O vs SSF Suitability"), x = "H2O", y = "SSF")
    ggsave(file.path(fig_dir, sprintf("%s_Compare_Scatter_H2O_SSF.png", ele)), p_h2o_ssf, width = 6, height = 6)

    # SSDM vs SSF
    p_ssdm_ssf <- ggplot(all_data, aes(x = SSDM, y = SSF, color = Run)) +
      geom_point(alpha = 0.1, size = 0.5) +
      geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray") +
      scale_color_manual(values = c("A" = "blue", "B" = "#FF7F7F")) +
      theme_minimal() +
      labs(title = paste(ele, "- SSDM vs SSF Suitability"), x = "SSDM", y = "SSF")
    ggsave(file.path(fig_dir, sprintf("%s_Compare_Scatter_SSDM_SSF.png", ele)), p_ssdm_ssf, width = 6, height = 6)

    # ---------------------------------------------------------
    # 2. Suitability vs Uncertainty (Scatter & Hexbin)
    # ---------------------------------------------------------
    p_suit_unc <- ggplot(all_data, aes(x = Hybrid, y = Uncertainty, color = Run)) +
      geom_point(alpha = 0.1, size = 0.5) +
      geom_smooth(method = "gam", se = FALSE) +
      scale_color_manual(values = c("A" = "blue", "B" = "#FF7F7F")) +
      theme_minimal() +
      labs(title = paste(ele, "- Hybrid Suitability vs Uncertainty"), x = "Hybrid Suitability", y = "Uncertainty")
    ggsave(file.path(fig_dir, sprintf("%s_Compare_Scatter_Suitability_Uncertainty.png", ele)), p_suit_unc, width = 7, height = 5)

    # Hexbin version for density
    p_hex <- ggplot(all_data, aes(x = Hybrid, y = Uncertainty)) +
      geom_hex(bins = 50) +
      scale_fill_viridis_c() +
      facet_wrap(~Run) +
      theme_minimal() +
      labs(title = paste(ele, "- Suitability vs Uncertainty Density"), x = "Hybrid Suitability", y = "Uncertainty")
    ggsave(file.path(fig_dir, sprintf("%s_Compare_Hexbin_Suitability_Uncertainty.png", ele)), p_hex, width = 10, height = 5)

    # ---------------------------------------------------------
    # 3. Jaccard Threshold Analysis (0.75 Threshold)
    # ---------------------------------------------------------
    # Calculate Jaccard Index at threshold 0.75 for each pair in each run
    threshold <- 0.75

    calc_jaccard <- function(a, b, thresh) {
      bin_a <- a > thresh
      bin_b <- b > thresh
      intersection <- sum(bin_a & bin_b, na.rm = TRUE)
      union <- sum(bin_a | bin_b, na.rm = TRUE)
      if (union == 0) {
        return(0)
      }
      return(intersection / union)
    }

    jaccard_results <- data.frame()

    for (r in unique(all_data$Run)) {
      sub_data <- all_data[all_data$Run == r, ]

      j_h2o_ssdm <- calc_jaccard(sub_data$H2O, sub_data$SSDM, threshold)
      j_h2o_ssf <- calc_jaccard(sub_data$H2O, sub_data$SSF, threshold)
      j_ssdm_ssf <- calc_jaccard(sub_data$SSDM, sub_data$SSF, threshold)

      jaccard_results <- rbind(jaccard_results, data.frame(
        Run = r,
        Pair = c("H2O-SSDM", "H2O-SSF", "SSDM-SSF"),
        Jaccard = c(j_h2o_ssdm, j_h2o_ssf, j_ssdm_ssf)
      ))
    }

    # Bar plot of Jaccard Indices
    p_jaccard <- ggplot(jaccard_results, aes(x = Pair, y = Jaccard, fill = Run)) +
      geom_bar(stat = "identity", position = "dodge") +
      scale_fill_manual(values = c("A" = "blue", "B" = "#FF7F7F")) +
      theme_minimal() +
      labs(title = paste(ele, "- Jaccard Similarity (Threshold > 0.75)"), y = "Jaccard Index")
    ggsave(file.path(fig_dir, sprintf("%s_Compare_Bar_Jaccard_075.png", ele)), p_jaccard, width = 6, height = 5)

    # ---------------------------------------------------------
    # 4. Spindle Graphs (Violin Plots of Suitability Distributions)
    # ---------------------------------------------------------
    # Reshape for violin plot
    long_data <- all_data %>%
      select(Run, H2O, SSDM, SSF, Hybrid) %>%
      pivot_longer(cols = c(H2O, SSDM, SSF, Hybrid), names_to = "Method", values_to = "Suitability")

    p_violin <- ggplot(long_data, aes(x = Method, y = Suitability, fill = Run)) +
      geom_violin(alpha = 0.7, trim = TRUE, scale = "width") +
      geom_boxplot(width = 0.1, position = position_dodge(0.9), fill = "white", outlier.shape = NA) +
      scale_fill_manual(values = c("A" = "blue", "B" = "#FF7F7F")) +
      theme_minimal() +
      labs(title = paste(ele, "- Suitability Distributions (Spindle/Violin)"), y = "Suitability")
    ggsave(file.path(fig_dir, sprintf("%s_Compare_Violin_Suitability.png", ele)), p_violin, width = 8, height = 6)

    # Clean up for next elephant
    rm(ele_plot_data)
  }
} # End Elephant Loop
message("\nAll Phase 4/5 processing complete.")
