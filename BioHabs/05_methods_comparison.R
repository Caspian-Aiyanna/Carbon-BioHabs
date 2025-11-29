# Three-Method Comparison: H2O vs SSDM vs SSF
# Compares habitat suitability predictions across three modeling frameworks
# Handles both single predictions and replicates with flexible file discovery
# Output: Correlation metrics, difference maps, and hotspot overlap statistics

suppressPackageStartupMessages({
  library(optparse)
  library(terra)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(readr)
  library(ggplot2)
  library(scales)
  library(tools)
  library(grid)
  library(purrr)
  library(stringr)
})

# Detect script location and project root
.this_file <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  filearg <- grep("^--file=", args, value = TRUE)
  if (length(filearg)) return(normalizePath(sub("^--file=", "", filearg)))
  if (!is.null(sys.frames()) && length(sys.frames())) {
    fi <- tryCatch(normalizePath(sys.frames()[[1]]$ofile), error = function(e) NA_character_)
    if (!is.na(fi)) return(fi)
  }
  stop("Cannot determine script path. Run via Rscript or set working directory to project root.")
}
.script <- dirname(.this_file())
.root   <- normalizePath(file.path(.script, ".."), winslash = "/", mustWork = TRUE)

source(file.path(.root, "R", "utils_io.R"))
source(file.path(.root, "R", "utils_repro.R"))

`%||%` <- function(x, y) if (is.null(x)) y else x

# Command-line options
opt <- list(
  make_option(c("--run"),     type = "character", default = "B",    help = "Run tag: A or B [default A]"),
  make_option(c("--mode"),    type = "character", default = "FAST", help = "REPRO or FAST [default FAST]"),
  make_option(c("--species"), type = "character", default = NULL,   help = "Optional single species (e.g., E6B)"),
  make_option(c("--q"),       type = "double",    default = 0.75,   help = "Hotspot quantile [0.75]")
)
opts <- parse_args(OptionParser(option_list = opt))

cfg <- read_config()
if (!is.null(opts$mode)) cfg$mode <- toupper(opts$mode)
mode <- set_mode(cfg)  # will set threads/cores etc.

# --- paths & logs ------------------------------------------------------------
h2o_root     <- cfg$paths$results_h2o  %||% file.path("results","H2O")
ssdm_root    <- cfg$paths$results_ssdm %||% file.path("results","SSDM")
ssf_root     <- cfg$paths$results_ssf  %||% file.path("results","SSF")
compare_root <- file.path("results","compare", toupper(opts$run))
logs_dir     <- cfg$paths$logs %||% "logs"
dir_ensure(compare_root); dir_ensure(logs_dir)
dir_ensure(file.path(compare_root, "01_between_methods/rasters"))
dir_ensure(file.path(compare_root, "01_between_methods/plots"))
dir_ensure(file.path(compare_root, "02_tables"))
dir_ensure(file.path(compare_root, "03_maps"))

logf <- file.path(logs_dir, sprintf("05_compare_methods_%s.log", toupper(opts$run)))
log_line(sprintf("Starting 05_h2o_vs_ssdm_vs_ssf_results.R (mode=%s, run=%s, q=%.2f)",
                 cfg$mode, toupper(opts$run), opts$q), logf)

# --- helpers -----------------------------------------------------------------
align_to <- function(r_ref, r_move, categorical = FALSE) {
  if (!compareGeom(r_ref, r_move, stopOnError = FALSE)) {
    r_move <- project(r_move, crs(r_ref), method = if (categorical) "near" else "bilinear")
    r_move <- resample(r_move, r_ref,  method = if (categorical) "near" else "bilinear")
    r_move <- crop(r_move, r_ref)
  }
  m <- !is.na(r_ref) & !is.na(r_move)
  r_ref  <- mask(r_ref,  m, maskvalues = 0)
  r_move <- mask(r_move, m, maskvalues = 0)
  list(r1 = r_ref, r2 = r_move)
}

jaccard_binary <- function(b1, b2) {
  inter <- suppressWarnings(global(b1 & b2, "sum", na.rm = TRUE)[[1]])
  union <- suppressWarnings(global(b1 | b2, "sum", na.rm = TRUE)[[1]])
  ifelse(union == 0, NA_real_, inter / union)
}

plot_raster_continuous <- function(
  r, title, out_png,
  center0 = FALSE,
  legend_title = "Suitability",
  barheight_cm = 12,
  barwidth_cm  = 0.7,
  base_size    = 12
) {
  df <- as.data.frame(r, xy = TRUE, na.rm = TRUE)
  names(df) <- c("x","y","val")
  p <- ggplot(df, aes(x = x, y = y, fill = val)) +
    geom_raster() + coord_equal() +
    labs(title = title, x = NULL, y = NULL, fill = legend_title) +
    theme_minimal(base_size = base_size) +
    theme(
      axis.text = element_blank(),
      panel.grid = element_blank(),
      plot.title = element_text(hjust = 0, face = "bold"),
      legend.position = "right",
      legend.title = element_text(size = base_size + 1, face = "bold"),
      legend.text  = element_text(size = base_size)
    )
  if (center0) {
    lim <- max(abs(range(df$val, na.rm = TRUE)))
    p <- p + scale_fill_gradientn(
      colors = c("#0c97f3","#ffffff","#e41b1e"),
      limits = c(-lim, lim),
      na.value = "#ffffff",
      guide = guide_colorbar(
        barheight = unit(barheight_cm, "cm"),
        barwidth  = unit(barwidth_cm,  "cm"),
        title.position = "top", title.hjust = 0
      )
    )
  } else {
    p <- p + scale_fill_viridis_c(
      na.value = NA,
      guide = guide_colorbar(
        barheight = unit(barheight_cm, "cm"),
        barwidth  = unit(barwidth_cm,  "cm"),
        title.position = "top", title.hjust = 0
      )
    )
  }
  ggsave(out_png, p, width = 8.2, height = 6.2, dpi = 300)
}

compute_pair_metrics <- function(rA, rB, q) {
  al <- align_to(rA, rB, categorical = FALSE)
  rA <- al$r1; rB <- al$r2
  vA <- values(rA, mat=FALSE); vB <- values(rB, mat=FALSE)
  idx <- !is.na(vA) & !is.na(vB)
  pearson  <- suppressWarnings(cor(vA[idx], vB[idx], method="pearson"))
  spearman <- suppressWarnings(cor(vA[idx], vB[idx], method="spearman"))
  rmse     <- sqrt(mean((vA[idx] - vB[idx])^2))
  mae      <- mean(abs(vA[idx] - vB[idx]))
  q1 <- quantile(vA[idx], q, na.rm=TRUE)
  q2 <- quantile(vB[idx], q, na.rm=TRUE)
  jac <- jaccard_binary(rA > q1, rB > q2)
  list(pearson=pearson, spearman=spearman, rmse=rmse, mae=mae, jaccard_hotspot_q=jac,
       rA=rA, rB=rB)
}

# --- discover species (present in at least one of H2O/SSDM/SSF) --------------
run_tag <- toupper(opts$run)
h2o_run_dir  <- file.path(h2o_root,  run_tag)
ssdm_run_dir <- file.path(ssdm_root, run_tag)
ssf_run_dir  <- file.path(ssf_root,  run_tag)

species_h2o  <- if (dir.exists(h2o_run_dir))  list.dirs(h2o_run_dir,  full.names = FALSE, recursive = FALSE) else character()
species_ssdm <- if (dir.exists(ssdm_run_dir)) list.dirs(ssdm_run_dir, full.names = FALSE, recursive = FALSE) else character()
species_ssf  <- if (dir.exists(ssf_run_dir))  list.dirs(ssf_run_dir,  full.names = FALSE, recursive = FALSE) else character()

if (!is.null(opts$species)) {
  species_h2o  <- intersect(species_h2o,  opts$species)
  species_ssdm <- intersect(species_ssdm, opts$species)
  species_ssf  <- intersect(species_ssf,  opts$species)
}
species_all <- sort(unique(c(species_h2o, species_ssdm, species_ssf)))
stopifnot("No species folders found in results/{H2O,SSDM,SSF}/<RUN>/." = length(species_all) > 0)
log_line(sprintf("Species discovered: %s", paste(species_all, collapse = ", ")), logf)

# ---- replicate listers (flexible) -------------------------------------------
list_h2o_rep_preds <- function(sp) {
  rep_dir <- file.path(h2o_run_dir, sp, "replicates")
  if (!dir.exists(rep_dir)) return(tibble(rep = integer(), path = character()))
  reps <- list.dirs(rep_dir, full.names = FALSE, recursive = FALSE)
  if (!length(reps)) return(tibble(rep = integer(), path = character()))
  rows <- lapply(reps, function(rn) {
    rp <- file.path(rep_dir, rn)
    # prefer prediction_<SP>.tif, then pred.tif, else first .tif
    cand <- file.path(rp, sprintf("prediction_%s.tif", sp))
    if (!file.exists(cand)) {
      cand2 <- file.path(rp, "pred.tif")
      if (file.exists(cand2)) cand <- cand2
    }
    if (!file.exists(cand)) {
      tifs <- list.files(rp, pattern="\\.tif$", full.names=TRUE)
      if (length(tifs)) cand <- tifs[1]
    }
    if (file.exists(cand)) tibble(rep = as.integer(sub("^rep","", rn)), path = cand) else NULL
  })
  bind_rows(rows)
}

list_ssdm_rep_preds <- function(sp) {
  rep_dir <- file.path(ssdm_run_dir, sp, "replicates")
  if (!dir.exists(rep_dir)) return(tibble(rep = integer(), path = character()))
  reps <- list.dirs(rep_dir, full.names = FALSE, recursive = FALSE)
  if (!length(reps)) return(tibble(rep = integer(), path = character()))
  rows <- lapply(reps, function(rn) {
    rp <- file.path(rep_dir, rn)
    cand <- file.path(rp, sprintf("ESDM_%s_rep%d.tif", sp, as.integer(sub("^rep","", rn))))
    if (!file.exists(cand)) {
      tifs <- list.files(rp, pattern="\\.tif$", full.names=TRUE)
      if (length(tifs)) cand <- tifs[1]
    }
    if (file.exists(cand)) tibble(rep = as.integer(sub("^rep","", rn)), path = cand) else NULL
  })
  bind_rows(rows)
}

# ---- singles readers (strict to your names, with final safe fallback) -------
read_h2o_single <- function(sp) {
  d <- file.path(h2o_run_dir, sp)
  cand <- file.path(d, sprintf("prediction_%s.tif", sp))
  if (file.exists(cand)) return(rast(cand))
  tifs <- if (dir.exists(d)) list.files(d, pattern="\\.tif$", full.names=TRUE) else character()
  if (length(tifs)) return(rast(tifs[1]))
  NULL
}

read_ssdm_single <- function(sp) {
  d <- file.path(ssdm_run_dir, sp)
  cand <- file.path(d, sprintf("ESDM_%s.tif", sp))
  if (file.exists(cand)) return(rast(cand))
  tifs <- if (dir.exists(d)) list.files(d, pattern="\\.tif$", full.names=TRUE) else character()
  if (length(tifs)) return(rast(tifs[1]))
  NULL
}

read_ssf_single <- function(sp) {
  d <- file.path(ssf_run_dir, sp)
  cand <- file.path(d, sprintf("%s_SSF_rsf_0to1.tif", sp))
  if (file.exists(cand)) return(rast(cand))
  # tolerate missing extension variant
  cand2 <- file.path(d, sprintf("%s_SSF_rsf_0to1", sp))
  if (file.exists(cand2)) return(rast(cand2))
  tifs <- if (dir.exists(d)) list.files(d, pattern="\\.tif$", full.names=TRUE) else character()
  if (length(tifs)) return(rast(tifs[1]))
  NULL
}

# --- main --------------------------------------------------------------------
per_rep_rows    <- list()
per_single_rows <- list()

for (sp in species_all) {
  # try replicates for H2O and SSDM, SSF is single
  h2o_tbl  <- list_h2o_rep_preds(sp)
  ssdm_tbl <- list_ssdm_rep_preds(sp)
  r_ssf    <- read_ssf_single(sp)

  have_reps <- nrow(h2o_tbl) > 0 && nrow(ssdm_tbl) > 0
  if (have_reps) {
    reps_both <- intersect(h2o_tbl$rep, ssdm_tbl$rep)
    if (!length(reps_both)) have_reps <- FALSE
  }

  if (have_reps) {
    log_line(sprintf("[replicates] %s | shared: %s",
                     sp, paste(sort(reps_both), collapse = ", ")), logf)

    predsH <- list(); predsS <- list()
    diffs_HS <- list(); diffs_HF <- list(); diffs_SF <- list()

    for (rk in sort(reps_both)) {
      f_h2o  <- h2o_tbl$path[h2o_tbl$rep == rk][1]
      f_ssdm <- ssdm_tbl$path[ssdm_tbl$rep == rk][1]

      r_h2o  <- rast(f_h2o)
      r_ssdm <- rast(f_ssdm)

      # Align all to H2O grid
      al_s   <- align_to(r_h2o, r_ssdm, categorical = FALSE)
      r_h2o  <- al_s$r1
      r_ssdm <- al_s$r2

      # If SSF exists, align; else skip pairs with SSF
      if (!is.null(r_ssf)) {
        al_f <- align_to(r_h2o, r_ssf, categorical = FALSE)
        r_h2o_al <- al_f$r1
        r_ssf_al <- al_f$r2
      } else {
        r_h2o_al <- r_h2o
        r_ssf_al <- NULL
      }

      # H2O vs SSDM
      m_hs <- compute_pair_metrics(r_h2o, r_ssdm, opts$q)

      # H2O vs SSF (if SSF available)
      if (!is.null(r_ssf_al)) {
        m_hf <- compute_pair_metrics(r_h2o_al, r_ssf_al, opts$q)
      } else {
        m_hf <- list(pearson=NA_real_, spearman=NA_real_, rmse=NA_real_, mae=NA_real_, jaccard_hotspot_q=NA_real_)
      }

      # SSDM vs SSF (if SSF available)
      if (!is.null(r_ssf_al)) {
        m_sf <- compute_pair_metrics(r_ssdm, r_ssf_al, opts$q)
      } else {
        m_sf <- list(pearson=NA_real_, spearman=NA_real_, rmse=NA_real_, mae=NA_real_, jaccard_hotspot_q=NA_real_)
      }

      per_rep_rows[[length(per_rep_rows) + 1L]] <- tibble(
        run = run_tag, dataset = sp, replicate = rk,
        h2o_vs_ssdm_pearson  = m_hs$pearson,
        h2o_vs_ssdm_spearman = m_hs$spearman,
        h2o_vs_ssdm_rmse     = m_hs$rmse,
        h2o_vs_ssdm_mae      = m_hs$mae,
        h2o_vs_ssdm_jaccard  = m_hs$jaccard_hotspot_q,
        h2o_vs_ssf_pearson   = m_hf$pearson,
        h2o_vs_ssf_spearman  = m_hf$spearman,
        h2o_vs_ssf_rmse      = m_hf$rmse,
        h2o_vs_ssf_mae       = m_hf$mae,
        h2o_vs_ssf_jaccard   = m_hf$jaccard_hotspot_q,
        ssdm_vs_ssf_pearson  = m_sf$pearson,
        ssdm_vs_ssf_spearman = m_sf$spearman,
        ssdm_vs_ssf_rmse     = m_sf$rmse,
        ssdm_vs_ssf_mae      = m_sf$mae,
        ssdm_vs_ssf_jaccard  = m_sf$jaccard_hotspot_q
      )

      predsH[[length(predsH) + 1L]] <- m_hs$rA
      predsS[[length(predsS) + 1L]] <- m_hs$rB
      diffs_HS[[length(diffs_HS) + 1L]] <- (m_hs$rA - m_hs$rB)
      if (!is.null(r_ssf_al)) {
        diffs_HF[[length(diffs_HF) + 1L]] <- (r_h2o_al - r_ssf_al)
        diffs_SF[[length(diffs_SF) + 1L]] <- (r_ssdm - r_ssf_al)
      }
    }

    # Mean maps (replicate means)
    rH_mean <- do.call(c, predsH) |> mean(na.rm = TRUE)
    rS_mean <- do.call(c, predsS) |> mean(na.rm = TRUE)
    writeRaster(rH_mean, file.path(compare_root,"01_between_methods/rasters", sprintf("%s_H2O_mean.tif", sp)),
                overwrite=TRUE, datatype="FLT4S", gdal=c("COMPRESS=LZW","BIGTIFF=IF_SAFER"))
    writeRaster(rS_mean, file.path(compare_root,"01_between_methods/rasters", sprintf("%s_SSDM_mean.tif", sp)),
                overwrite=TRUE, datatype="FLT4S", gdal=c("COMPRESS=LZW","BIGTIFF=IF_SAFER"))
    plot_raster_continuous(rH_mean, paste0(sp," — H2O mean suitability"),
                           file.path(compare_root,"03_maps", sprintf("%s_H2O_mean.png", sp)))
    plot_raster_continuous(rS_mean, paste0(sp," — SSDM mean suitability"),
                           file.path(compare_root,"03_maps", sprintf("%s_SSDM_mean.png", sp)))

    if (length(diffs_HS)) {
      d_hs_mean <- do.call(c, diffs_HS) |> mean(na.rm = TRUE)
      writeRaster(d_hs_mean, file.path(compare_root,"01_between_methods/rasters",
                                       sprintf("diff_mean_H2O_minus_SSDM_%s.tif", sp)),
                  overwrite=TRUE, datatype="FLT4S", gdal=c("COMPRESS=LZW","BIGTIFF=IF_SAFER"))
      plot_raster_continuous(d_hs_mean, paste0("Δ (H2O − SSDM): ", sp),
                             file.path(compare_root,"01_between_methods/plots",
                                       sprintf("diff_mean_H2O_minus_SSDM_%s.png", sp)),
                             center0=TRUE, legend_title="Δ suitability")
    }

    if (!is.null(r_ssf) && length(diffs_HF)) {
      rF_mean <- mean(do.call(c, lapply(diffs_HF, function(x) x + r_ssf - r_ssf)), na.rm = TRUE) # dummy to get same grid
      # Better: align SSF to H2O mean grid and save map
      alF <- align_to(rH_mean, r_ssf, categorical = FALSE)
      rF_mean <- alF$r2
      writeRaster(rF_mean, file.path(compare_root,"01_between_methods/rasters", sprintf("%s_SSF.tif", sp)),
                  overwrite=TRUE, datatype="FLT4S", gdal=c("COMPRESS=LZW","BIGTIFF=IF_SAFER"))
      plot_raster_continuous(rF_mean, paste0(sp," — SSF suitability"),
                             file.path(compare_root,"03_maps", sprintf("%s_SSF.png", sp)))
      d_hf_mean <- do.call(c, diffs_HF) |> mean(na.rm = TRUE)
      d_sf_mean <- do.call(c, diffs_SF) |> mean(na.rm = TRUE)
      writeRaster(d_hf_mean, file.path(compare_root,"01_between_methods/rasters",
                                       sprintf("diff_mean_H2O_minus_SSF_%s.tif", sp)),
                  overwrite=TRUE, datatype="FLT4S", gdal=c("COMPRESS=LZW","BIGTIFF=IF_SAFER"))
      writeRaster(d_sf_mean, file.path(compare_root,"01_between_methods/rasters",
                                       sprintf("diff_mean_SSDM_minus_SSF_%s.tif", sp)),
                  overwrite=TRUE, datatype="FLT4S", gdal=c("COMPRESS=LZW","BIGTIFF=IF_SAFER"))
      plot_raster_continuous(d_hf_mean, paste0("Δ (H2O − SSF): ", sp),
                             file.path(compare_root,"01_between_methods/plots",
                                       sprintf("diff_mean_H2O_minus_SSF_%s.png", sp)),
                             center0=TRUE, legend_title="Δ suitability")
      plot_raster_continuous(d_sf_mean, paste0("Δ (SSDM − SSF): ", sp),
                             file.path(compare_root,"01_between_methods/plots",
                                       sprintf("diff_mean_SSDM_minus_SSF_%s.png", sp)),
                             center0=TRUE, legend_title="Δ suitability")
    }

  } else {
    # --- SINGLE predictions fallback -----------------------------------------
    rH <- read_h2o_single(sp)
    rS <- read_ssdm_single(sp)
    rF <- read_ssf_single(sp)

    if (is.null(rH) && is.null(rS) && is.null(rF)) {
      log_line(sprintf("Skip %s — no single rasters found in H2O/SSDM/SSF.", sp), logf)
      next
    }
    log_line(sprintf("[single] %s | H2O=%s, SSDM=%s, SSF=%s",
                     sp, !is.null(rH), !is.null(rS), !is.null(rF)), logf)

    # Available pairs among the three
    if (!is.null(rH) && !is.null(rS)) {
      m <- compute_pair_metrics(rH, rS, opts$q)
      per_single_rows[[length(per_single_rows)+1L]] <- tibble(
        run=run_tag, dataset=sp, pair="H2O-SSDM",
        pearson=m$pearson, spearman=m$spearman, rmse=m$rmse, mae=m$mae,
        jaccard_hotspot_q=m$jaccard_hotspot_q
      )
      writeRaster(m$rA - m$rB, file.path(compare_root,"01_between_methods/rasters",
                 sprintf("diff_single_H2O_minus_SSDM_%s.tif", sp)),
                 overwrite=TRUE, datatype="FLT4S", gdal=c("COMPRESS=LZW","BIGTIFF=IF_SAFER"))
      plot_raster_continuous(m$rA - m$rB, paste0("Δ (H2O − SSDM, single): ", sp),
                             file.path(compare_root,"01_between_methods/plots",
                                       sprintf("diff_single_H2O_minus_SSDM_%s.png", sp)),
                             center0=TRUE, legend_title="Δ suitability")
      plot_raster_continuous(m$rA, paste0(sp," — H2O (single)"),
                             file.path(compare_root,"03_maps", sprintf("%s_H2O_single.png", sp)))
      plot_raster_continuous(m$rB, paste0(sp," — SSDM (single)"),
                             file.path(compare_root,"03_maps", sprintf("%s_SSDM_single.png", sp)))
    }
    if (!is.null(rH) && !is.null(rF)) {
      m <- compute_pair_metrics(rH, rF, opts$q)
      per_single_rows[[length(per_single_rows)+1L]] <- tibble(
        run=run_tag, dataset=sp, pair="H2O-SSF",
        pearson=m$pearson, spearman=m$spearman, rmse=m$rmse, mae=m$mae,
        jaccard_hotspot_q=m$jaccard_hotspot_q
      )
      writeRaster(m$rA - m$rB, file.path(compare_root,"01_between_methods/rasters",
                 sprintf("diff_single_H2O_minus_SSF_%s.tif", sp)),
                 overwrite=TRUE, datatype="FLT4S", gdal=c("COMPRESS=LZW","BIGTIFF=IF_SAFER"))
      plot_raster_continuous(m$rA - m$rB, paste0("Δ (H2O − SSF, single): ", sp),
                             file.path(compare_root,"01_between_methods/plots",
                                       sprintf("diff_single_H2O_minus_SSF_%s.png", sp)),
                             center0=TRUE, legend_title="Δ suitability")
      plot_raster_continuous(m$rB, paste0(sp," — SSF"),
                             file.path(compare_root,"03_maps", sprintf("%s_SSF.png", sp)))
    }
    if (!is.null(rS) && !is.null(rF)) {
      m <- compute_pair_metrics(rS, rF, opts$q)
      per_single_rows[[length(per_single_rows)+1L]] <- tibble(
        run=run_tag, dataset=sp, pair="SSDM-SSF",
        pearson=m$pearson, spearman=m$spearman, rmse=m$rmse, mae=m$mae,
        jaccard_hotspot_q=m$jaccard_hotspot_q
      )
      writeRaster(m$rA - m$rB, file.path(compare_root,"01_between_methods/rasters",
                 sprintf("diff_single_SSDM_minus_SSF_%s.tif", sp)),
                 overwrite=TRUE, datatype="FLT4S", gdal=c("COMPRESS=LZW","BIGTIFF=IF_SAFER"))
      plot_raster_continuous(m$rA - m$rB, paste0("Δ (SSDM − SSF, single): ", sp),
                             file.path(compare_root,"01_between_methods/plots",
                                       sprintf("diff_single_SSDM_minus_SSF_%s.png", sp)),
                             center0=TRUE, legend_title="Δ suitability")
    }
  }
}

# ----- Tables: replicates (may be empty) -------------------------------------
per_rep_df <- if (length(per_rep_rows)) bind_rows(per_rep_rows) else
  tibble(run=character(), dataset=character(), replicate=integer(),
         h2o_vs_ssdm_pearson=double(),  h2o_vs_ssdm_spearman=double(),
         h2o_vs_ssdm_rmse=double(),     h2o_vs_ssdm_mae=double(),
         h2o_vs_ssdm_jaccard=double(),
         h2o_vs_ssf_pearson=double(),   h2o_vs_ssf_spearman=double(),
         h2o_vs_ssf_rmse=double(),      h2o_vs_ssf_mae=double(),
         h2o_vs_ssf_jaccard=double(),
         ssdm_vs_ssf_pearson=double(),  ssdm_vs_ssf_spearman=double(),
         ssdm_vs_ssf_rmse=double(),     ssdm_vs_ssf_mae=double(),
         ssdm_vs_ssf_jaccard=double())

out_per_rep <- file.path(compare_root, "02_tables", "per_rep_metrics_long.csv")
readr::write_csv(per_rep_df, out_per_rep)
log_line(sprintf("Wrote per-replicate metrics (long): %s [n=%d]", out_per_rep, nrow(per_rep_df)), logf)

if (nrow(per_rep_df) > 0) {
  summary_long <- per_rep_df %>%
    mutate(pair = "replicate_means") %>%  # tag
    group_by(run, dataset) %>%
    summarise(
      n_reps = dplyr::n(),
      across(starts_with("h2o_vs_ssdm_") | starts_with("h2o_vs_ssf_") | starts_with("ssdm_vs_ssf_"),
             list(mean = ~mean(.x, na.rm = TRUE),
                  sd   = ~sd(.x,   na.rm = TRUE)),
             .names = "{.col}_{.fn}")
    ) %>% ungroup()
  out_summary <- file.path(compare_root, "02_tables", "per_species_summary_long.csv")
  readr::write_csv(summary_long, out_summary)
  log_line(sprintf("Wrote per-species summary (replicates): %s", out_summary), logf)
} else {
  log_line("No replicate rows -> skipping replicate summary.", logf)
}

# ----- Tables: singles (may be empty) ----------------------------------------
per_single_df <- if (length(per_single_rows)) bind_rows(per_single_rows) else
  tibble(run=character(), dataset=character(), pair=character(),
         pearson=double(), spearman=double(), rmse=double(), mae=double(), jaccard_hotspot_q=double())
out_per_single <- file.path(compare_root, "02_tables", "per_single_metrics_long.csv")
readr::write_csv(per_single_df, out_per_single)
log_line(sprintf("Wrote single-pred metrics (long): %s [n=%d]", out_per_single, nrow(per_single_df)), logf)

if (nrow(per_single_df) > 0) {
  summary_single <- per_single_df %>%
    group_by(run, dataset, pair) %>%
    summarise(
      n_pairs = dplyr::n(),
      across(c(pearson, spearman, rmse, mae, jaccard_hotspot_q),
             list(mean=~mean(.x, na.rm=TRUE), sd=~sd(.x, na.rm=TRUE)),
             .names="{.col}_{.fn}")
    ) %>% ungroup()
  out_summary_single <- file.path(compare_root, "02_tables", "per_species_summary_single.csv")
  readr::write_csv(summary_single, out_summary_single)
  log_line(sprintf("Wrote per-species summary (single): %s", out_summary_single), logf)
} else {
  log_line("No single rows -> skipping single summary.", logf)
}

log_line(sprintf("Done. Outputs in: %s", compare_root), logf)
