#!/usr/bin/env Rscript
# =============================================================================
# FULL SCRIPT: Copernicus CDS → AOI-cropped GeoTIFFs with Robust NetCDF Handling
# - Reads AOI vector (SHP/GPKG/GeoJSON)
# - User-configurable dataset preset (ERA5-Land monthly OR CMIP6 monthly template)
# - Select variables, years, months, output CRS & resolution
# - Downloads via Python cdsapi (reticulate), with licence-aware error handling
# - Robust NetCDF opener: handles ZIP deliveries; falls back to stars+ncdf4 if terra/GDAL lacks NetCDF
# - Crops/masks, reprojects, resamples
# - Writes one GeoTIFF per time slice into tidy folders, with clean names
#
# REQUIREMENTS:
#   1) CDS credentials file:
#      Windows: C:/Users/<YOU>/.cdsapirc
#      Linux/Mac: ~/.cdsapirc
#      Contents:
#        url: https://cds.climate.copernicus.eu/api/v2
#        key: <UID>:<API-KEY>
#   2) Python available. If reticulate picks the wrong Python, set CONFIG$python_bin below.
#
# USAGE:
#   - Edit CONFIG block below.
#   - Run:  Rscript cds_to_tifs_full.R
# =============================================================================

# --------------------- PACKAGES ---------------------------------------------
options(repos = c(CRAN = "https://cloud.r-project.org"))
req <- c("sf", "terra", "reticulate", "jsonlite")
for (p in req) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
suppressPackageStartupMessages({
  library(sf)
  library(terra)
  library(reticulate)
  library(jsonlite)
})

# --------------------- CONFIG (EDIT THIS) -----------------------------------
CONFIG <- list(
  # AOI vector path (SHP/GPKG/GeoJSON). Multi-features will be dissolved.
  aoi_path = "BioHabs/data/shp/HV20233.shp", # <-- change to your AOI

  # Dataset preset: "era5_land_monthly" (ready) or "cmip6_monthly" (template)
  preset = "era5_land_monthly",

  # Time selection (Overridden by A/B Logic below)
  years = 2019:2021,
  months = sprintf("%02d", 1:12),

  # ERA5-Land monthly variables (names must match CDS dataset page)
  era5_vars = c(
    "2m_temperature",
    "total_precipitation",
    "leaf_area_index_high_vegetation",
    "leaf_area_index_low_vegetation",
    "volumetric_soil_water_layer_1",
    "potential_evaporation"
  ),

  # CMIP6 template (use if preset == "cmip6_monthly")
  cmip6 = list(
    temporal_resolution = "monthly",
    experiment          = "historical", # e.g., "ssp245"
    variable            = "near_surface_air_temperature",
    model               = "CNRM-CM6-1"
  ),

  # Output grid
  target_epsg = 4326, # 4326 (degrees) or a projected CRS (meters)
  target_res = 0.000269, # ~30m resolution (0.000269 degrees)

  # Output & behavior
  out_dir = "outputs/cds_tifs", # Overridden by A/B Logic
  overwrite = TRUE,
  keep_nc = TRUE,

  # Force a specific Python (optional; else reticulate auto-detects)
  python_bin = "C:/Users/lenovo/AppData/Local/Programs/Python/Python311/python.exe"
)

# Define Runs for A/B Periods
RUNS <- list(
  A = list(
    years = 2024:2025,
    out_dir = "BioHabs/data/envi/A/monthly"
  ),
  B = list(
    years = 2018:2023,
    out_dir = "BioHabs/data/envi/B/monthly"
  )
)

# --------------------- UTILITIES --------------------------------------------
ensure_dir <- function(x) if (!dir.exists(x)) dir.create(x, recursive = TRUE, showWarnings = FALSE)

aoi_to_cds_area <- function(aoi_sf) {
  aoi_wgs <- st_transform(aoi_sf, 4326)
  bb <- st_bbox(aoi_wgs)
  c(
    north = as.numeric(bb["ymax"]),
    west = as.numeric(bb["xmin"]),
    south = as.numeric(bb["ymin"]),
    east = as.numeric(bb["xmax"])
  )
}

write_time_tifs <- function(r, base_outdir, dataset_name, varname, overwrite = TRUE) {
  time_vals <- terra::time(r)
  if (is.null(time_vals)) {
    outdir <- file.path(base_outdir, dataset_name, varname)
    ensure_dir(outdir)
    outpath <- file.path(outdir, sprintf("%s_%s.tif", dataset_name, varname))
    terra::writeRaster(r, outpath, overwrite = overwrite)
    return(invisible(outpath))
  }
  for (i in seq_along(time_vals)) {
    dt <- time_vals[i]
    stamp <- tryCatch(format(as.Date(dt), "%Y%m"), error = function(e) sprintf("L%03d", i))
    outdir <- file.path(base_outdir, dataset_name, varname, substr(stamp, 1, 4))
    ensure_dir(outdir)
    outpath <- file.path(outdir, sprintf("%s_%s_%s.tif", dataset_name, varname, stamp))
    if (!file.exists(outpath) || overwrite) {
      terra::writeRaster(r[[i]], outpath, overwrite = overwrite)
    }
  }
}

cds_retrieve_safe <- function(client, dataset, request, target = NULL) {
  res <- try(
    {
      if (!is.null(target)) {
        client$retrieve(dataset, request, target = target)
      } else {
        client$retrieve(dataset, request)
      }
    },
    silent = TRUE
  )

  if (inherits(res, "try-error")) {
    err <- reticulate::py_last_error()
    msg <- paste(capture.output(str(err)), collapse = "\n")
    if (grepl("403 Client Error: Forbidden", msg) && grepl("licen", msg, ignore.case = TRUE)) {
      message("\n🚫 CDS licence not accepted for: ", dataset)
      message("➡ Please accept the licence on the CDS website for 'reanalysis-era5-land-monthly-means'.")
      stop("CDS licence not accepted.", call. = FALSE)
    } else {
      stop("CDS retrieve failed.\n", msg, call. = FALSE)
    }
  }
  res
}

ensure_pkg <- function(pkgs) {
  for (p in pkgs) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
}

maybe_unzip_nc <- function(path) {
  zlist <- try(utils::unzip(path, list = TRUE), silent = TRUE)
  if (inherits(zlist, "try-error")) {
    return(path)
  } # not a zip
  nc_rows <- zlist[grepl("\\.nc$", tolower(zlist$Name)), , drop = FALSE]
  if (nrow(nc_rows) == 0) stop("ZIP file does not contain a .nc")
  nc_inside <- nc_rows$Name[1]
  exdir <- file.path(dirname(path), "unzipped_nc")
  dir.create(exdir, showWarnings = FALSE, recursive = TRUE)
  utils::unzip(path, files = nc_inside, exdir = exdir, overwrite = TRUE)
  file.path(exdir, nc_inside)
}

open_nc_robust <- function(nc_path) {
  # 1. Try terra first (fastest, simplest)
  r <- try(terra::rast(nc_path), silent = TRUE)

  # Check if terra succeeded AND has valid extent/CRS
  if (!inherits(r, "try-error")) {
    e <- as.vector(terra::ext(r))
    is_generic <- all(e == c(0, ncol(r), 0, nrow(r)))
    is_no_crs <- is.na(terra::crs(r)) || terra::crs(r) == ""

    if (!is_generic && !is_no_crs) {
      return(r)
    }
    message("terra opened NetCDF but found no georeference. Switching to stars fallback...")
  } else {
    message("terra could not open NetCDF. Switching to stars fallback...")
  }

  # 2. Stars Fallback
  ensure_pkg(c("stars", "ncdf4", "ncmeta"))
  library(stars)
  library(ncmeta)

  # Get all variables
  meta <- ncmeta::nc_meta(nc_path)
  dinfo <- meta$variable_dimensions
  dnames <- meta$dimensions

  # Handle ncmeta structure variations
  if (!all(c("id", "name") %in% names(dnames))) {
    dnames <- data.frame(id = dnames$id, name = dnames$name)
  } else {
    dnames <- dnames[, c("id", "name")]
  }

  if (!"dimension_id" %in% names(dinfo)) {
    cand <- data.frame(variable = unique(dinfo$variable))
  } else {
    cand <- merge(dinfo, dnames, by.x = "dimension_id", by.y = "id", all.x = TRUE)
  }

  # Identify spatial variables
  # We want ALL variables that have spatial dimensions
  # Simplify logic: take all variables that are NOT dimensions and NOT 'number'/'expver'
  all_vars <- unique(cand$variable)
  spatial_vars <- setdiff(all_vars, c("time", "latitude", "longitude", "number", "expver", "valid_time"))

  # If that leaves nothing, try the previous heuristic
  if (length(spatial_vars) == 0) {
    spatial_vars <- names(Filter(function(z) {
      dn <- tolower(z$name)
      has_x <- any(grepl("lon|x", dn))
      has_y <- any(grepl("lat|y", dn))
      has_x && has_y
    }, split(cand, cand$variable)))
  }

  if (length(spatial_vars) == 0) {
    # Last resort: take everything that isn't a coordinate
    message("Warning: Could not identify spatial variables by dimensions. Using exclusion list.")
    nc_obj <- ncdf4::nc_open(nc_path)
    all_v <- names(nc_obj$var)
    ncdf4::nc_close(nc_obj)
    spatial_vars <- setdiff(all_v, c("time", "latitude", "longitude", "number", "expver", "valid_time"))
  }

  if (length(spatial_vars) == 0) spatial_vars <- unique(cand$variable)

  message("Found variables: ", paste(spatial_vars, collapse = ", "))

  r_list <- list()

  for (vname in spatial_vars) {
    message("  Reading variable: ", vname)
    s <- try(stars::read_ncdf(nc_path, var = vname, proxy = FALSE), silent = TRUE)
    if (inherits(s, "try-error")) next

    # Convert to terra
    rt <- try(terra::rast(s), silent = TRUE)

    # Fallback via temp TIF if direct conversion fails
    if (inherits(rt, "try-error")) {
      tmp_tif <- tempfile(pattern = paste0(vname, "_"), fileext = ".tif")
      if (is.na(sf::st_crs(s))) sf::st_crs(s) <- 4326

      # Handle degenerate dimensions for extent
      dims <- stars::st_dimensions(s)
      if (!is.null(dims$latitude) && is.na(dims$latitude$delta)) {
        lat_val <- dims$latitude$values[1]
        res_guess <- 0.1
        # We can't easily fix extent in stars object before writing without hacking attributes
        # So we write, read back, and fix extent in terra
      }

      stars::write_stars(s, tmp_tif)
      rt <- terra::rast(tmp_tif)

      # Fix extent if needed
      if (any(is.na(as.vector(terra::ext(rt))))) {
        lons <- stars::st_get_dimension_values(s, "longitude")
        lats <- stars::st_get_dimension_values(s, "latitude")
        res_x <- if (length(lons) > 1) diff(lons)[1] else 0.1
        res_y <- if (length(lats) > 1) diff(lats)[1] else 0.1
        xmin <- min(lons) - res_x / 2
        xmax <- max(lons) + res_x / 2
        ymin <- min(lats) - abs(res_y) / 2
        ymax <- max(lats) + abs(res_y) / 2
        terra::ext(rt) <- c(xmin, xmax, ymin, ymax)
        terra::crs(rt) <- "EPSG:4326"
      }
    }

    # Fix Time
    tm <- try(stars::st_get_dimension_values(s, "time"), silent = TRUE)
    if (!inherits(tm, "try-error") && !is.null(tm)) {
      suppressWarnings(try(terra::time(rt) <- as.POSIXct(tm), silent = TRUE))
    }

    # Fix Names: Use variable name + index/time
    # If multiple layers (time), terra names them file_1, file_2...
    # We want vname_1, vname_2...
    if (terra::nlyr(rt) > 1) {
      names(rt) <- paste0(vname, "_", 1:terra::nlyr(rt))
    } else {
      names(rt) <- vname
    }

    r_list[[vname]] <- rt
  }

  if (length(r_list) == 0) stop("Could not read any variables from NetCDF.")

  # Combine all variables into one raster
  # Note: They must have same extent/res/crs. Usually true for one NC file.
  r_final <- try(terra::rast(r_list), silent = TRUE)
  if (inherits(r_final, "try-error")) {
    message("Could not stack variables (different grids?). Returning list of rasters is not supported by main loop.")
    # Fallback: just return the first one or error?
    # For now, let's assume they stack. If not, we might need to resample them to the first one.
    r1 <- r_list[[1]]
    for (i in 2:length(r_list)) {
      terra::resample(r_list[[i]], r1)
    }
    r_final <- terra::rast(r_list)
  }

  return(r_final)
}

# --------------------- PREP: AOI + PYTHON + CDSAPI --------------------------
if (!file.exists(CONFIG$aoi_path)) stop("AOI not found: ", CONFIG$aoi_path)
aoi <- st_read(CONFIG$aoi_path, quiet = TRUE) |> st_make_valid()
if (nrow(aoi) > 1) aoi <- st_as_sf(st_union(aoi))
area_vec <- aoi_to_cds_area(aoi)

if (!is.null(CONFIG$python_bin)) {
  reticulate::use_python(CONFIG$python_bin, required = TRUE)
}
py_ok <- try(reticulate::py_config(), silent = TRUE)
if (inherits(py_ok, "try-error")) stop("Python not found. Set CONFIG$python_bin to your python.exe path.")

if (!reticulate::py_module_available("cdsapi")) {
  message("Installing Python 'cdsapi' via pip …")
  reticulate::py_install("cdsapi", pip = TRUE)
}

# CDS Credentials
CDS_URL <- "https://cds.climate.copernicus.eu/api"
CDS_KEY <- "f2837b9a-306e-4461-b2e7-769f2e4aec23"

cdsapi <- reticulate::import("cdsapi")
client <- cdsapi$Client(url = CDS_URL, key = CDS_KEY)

# --------------------- REQUEST BUILDERS -------------------------------------
build_req_era5 <- function(vars, years, months, area) {
  list(
    product_type = "monthly_averaged_reanalysis",
    variable     = as.list(vars),
    year         = as.list(as.character(years)),
    month        = as.list(months),
    time         = "00:00",
    area         = as.numeric(area), # [N, W, S, E] degrees
    format       = "netcdf"
  )
}

build_req_cmip6 <- function(cmip6, years, months, area) {
  period <- sprintf("%s-%s", min(as.integer(years)), max(as.integer(years)))
  list(
    temporal_resolution = cmip6$temporal_resolution,
    experiment          = cmip6$experiment,
    variable            = cmip6$variable,
    model               = cmip6$model,
    period              = period,
    month               = as.list(months),
    area                = as.numeric(area),
    format              = "netcdf"
  )
}

# --------------------- MAIN LOOP (A/B) --------------------------------------
dataset_name <- "reanalysis-era5-land-monthly-means"

for (run_id in names(RUNS)) {
  run_cfg <- RUNS[[run_id]]
  years <- as.character(run_cfg$years)
  out_dir <- run_cfg$out_dir

  message(sprintf("\n=== Processing Run %s (%s-%s) ===", run_id, min(years), max(years)))
  ensure_dir(out_dir)

  # Build Request
  if (CONFIG$preset == "era5_land_monthly") {
    request <- build_req_era5(CONFIG$era5_vars, years, CONFIG$months, area_vec)
  } else {
    stop("Only ERA5 supported for this loop currently.")
  }

  # Download
  tmp_dir <- file.path(tempdir(), "cds_download")
  ensure_dir(tmp_dir)
  nc_target <- file.path(tmp_dir, sprintf("ERA5_%s_%s.nc", run_id, format(Sys.time(), "%Y%m%d%H%M")))

  message("  Requesting data from CDS...")
  cds_retrieve_safe(client, dataset_name, request, target = nc_target)

  if (!file.exists(nc_target)) {
    warning(sprintf("  Download failed for Run %s", run_id))
    next
  }
  message("  Download complete.")

  # Unzip & Open
  nc_path <- maybe_unzip_nc(nc_target)
  r <- open_nc_robust(nc_path)

  # Process (Crop -> Project -> Resample)
  message("  Processing raster...")

  # 1. Fix extent/rotation if needed (0-360 to -180/180)
  if (any(ext(r)[1:2] > 180)) {
    message("    Detected 0-360 longitude, rotating to -180/180...")
    r <- terra::rotate(r)
  }

  # 2. Ensure CRS is set (ERA5 is usually WGS84)
  if (is.na(crs(r)) || crs(r) == "") {
    crs(r) <- "EPSG:4326"
  }

  # 3. Transform AOI to match Raster CRS for cropping
  aoi_r <- st_transform(aoi, crs(r))

  # Debug: Print extents
  message("    Raster Extent: ", paste(as.vector(ext(r)), collapse = ", "))
  message("    AOI Extent:    ", paste(as.vector(ext(aoi_r)), collapse = ", "))

  # 4. Crop and Mask
  # Use tryCatch to handle potential non-overlapping errors gracefully
  r <- tryCatch(
    {
      r_crop <- terra::crop(r, vect(aoi_r))
      terra::mask(r_crop, vect(aoi_r))
    },
    error = function(e) {
      stop(
        "Crop failed. The downloaded raster does not overlap with the AOI.\n",
        "Raster Extent: ", paste(as.vector(ext(r)), collapse = ", "), "\n",
        "AOI Extent: ", paste(as.vector(ext(aoi_r)), collapse = ", ")
      )
    }
  )

  target_crs_str <- sprintf("EPSG:%s", CONFIG$target_epsg)
  r <- terra::project(r, target_crs_str, method = "bilinear")

  aoi_t <- st_transform(aoi, CONFIG$target_epsg)
  tmpl <- terra::rast(
    ext = terra::ext(terra::vect(aoi_t)),
    crs = target_crs_str, resolution = CONFIG$target_res
  )
  r <- terra::resample(r, tmpl, method = "bilinear")

  # Write TIFs
  message("  Writing GeoTIFFs...")
  bn <- names(r)
  var_prefix <- sub("(_\\d+)$", "", bn)
  vars_found <- unique(var_prefix)

  for (v in vars_found) {
    idx <- which(var_prefix == v)
    r_v <- r[[idx]]
    write_time_tifs(r_v, out_dir, dataset_name, v, overwrite = CONFIG$overwrite)
  }

  # Optionally keep raw NetCDF
  if (isTRUE(CONFIG$keep_nc)) {
    final_nc_dir <- file.path(out_dir, dataset_name)
    ensure_dir(final_nc_dir)
    final_nc <- file.path(final_nc_dir, paste0(dataset_name, "_raw.nc"))
    file.copy(nc_path, final_nc, overwrite = TRUE)
  }

  message(sprintf("  Run %s complete. Output: %s", run_id, out_dir))
}

message("\n✅ All CDS downloads finished.")
# =============================================================================
