# H2O Utility Functions
# Helpers for initializing H2O clusters and making predictions on raster data

suppressPackageStartupMessages({
  library(h2o)
  library(terra)
})

# Calculate number of threads to use based on configuration
# Negative values mean "reserve N cores" (e.g., -2 = use all but 2 cores)
resolve_nthreads <- function(nthreads_cfg) {
  if (is.null(nthreads_cfg)) {
    return(1L)
  }
  if (nthreads_cfg > 0) {
    return(as.integer(nthreads_cfg))
  }

  # Negative value: reserve that many cores
  avail <- as.integer(parallel::detectCores(logical = TRUE))
  use <- max(1L, avail - abs(as.integer(nthreads_cfg)))
  return(use)
}

# Initialize H2O cluster with appropriate threading based on mode
init_h2o <- function(cfg) {
  mode <- toupper(cfg$mode %||% "REPRO")
  if (mode == "REPRO") {
    nthreads <- resolve_nthreads(cfg$h2o$nthreads_repro %||% 1L)
  } else {
    nthreads <- resolve_nthreads(cfg$h2o$nthreads_fast %||% -2L)
  }

  # Shut down any existing H2O cluster to avoid conflicts
  suppressWarnings(try(h2o::h2o.shutdown(prompt = FALSE), silent = TRUE))
  Sys.sleep(1)

  # Start fresh H2O cluster
  h2o::h2o.init(nthreads = nthreads, strict_version_check = FALSE)
  h2o::h2o.removeAll()

  invisible(nthreads)
}

# Make predictions on a raster using an H2O model
# Processes data in chunks to avoid memory issues
# Returns a SpatRaster with predicted probabilities
predict_raster_h2o <- function(r, model, block_rows = 50000) {
  stopifnot(inherits(r, "SpatRaster"))

  # Create empty output raster with same dimensions
  out <- rast(r[[1]])
  names(out) <- "p1"

  ncell_total <- ncell(r)
  idx <- seq_len(ncell_total)
  vals <- rep(NA_real_, ncell_total)

  # Process raster in chunks to manage memory
  chunks <- split(idx, ceiling(seq_along(idx) / block_rows))
  for (ch in chunks) {
    df <- as.data.frame(r[ch], na.rm = FALSE)
    hf <- as.h2o(df)
    pr <- h2o.predict(model, hf)
    vals[ch] <- as.vector(pr[["p1"]]) # Extract probability for class 1
  }

  out <- setValues(out, vals)
  return(out)
}

# Check if output files already exist (to skip redundant processing)
skip_if_done <- function(dir, files) {
  all(file.exists(file.path(dir, files)))
}

# Default value operator: return b if a is NULL
`%||%` <- function(a, b) if (is.null(a)) b else a
