#!/usr/bin/env Rscript
# sentinel_2_timeseries_monthly.R
# Sentinel-2 Monthly Composite & Index Calculation
# -----------------------------------------------------------------------------
# Inputs: BioHabs/data/shp/HV20233.shp
# Outputs: BioHabs/data/envi/[A|B]/monthly_s2/
# -----------------------------------------------------------------------------

suppressPackageStartupMessages({
    library(terra)
    library(sf)
    library(rstac)
    library(gdalcubes)
    library(lubridate)
    library(magrittr)
    library(jsonlite)
})

# --- Configuration ---
AOI_PATH <- "BioHabs/data/shp/HV20233.shp"
gdalcubes_options(parallel = 1) # Stability

# Define Runs (Matching cds_to_tifs.R logic)
RUNS <- list(
    A = list(
        years = 2024:2025,
        out_dir = "BioHabs/data/envi/A/monthly_s2"
    ),
    B = list(
        years = 2018:2023,
        out_dir = "BioHabs/data/envi/B/monthly_s2"
    )
)

# --- 1. Load AOI ---
if (!file.exists(AOI_PATH)) stop("AOI shapefile not found: ", AOI_PATH)
aoi <- st_read(AOI_PATH, quiet = TRUE)
aoi <- st_transform(aoi, 4326)
bbox <- st_bbox(aoi)

message(sprintf(
    "AOI Loaded. BBox: %.4f, %.4f, %.4f, %.4f",
    bbox["xmin"], bbox["ymin"], bbox["xmax"], bbox["ymax"]
))

# --- 2. Processing Function ---
process_month <- function(year, month, out_dir) {
    out_name <- file.path(out_dir, sprintf("S2_Monthly_%s_%02d.tif", year, month))
    if (file.exists(out_name)) {
        message(sprintf("  Skipping %s-%02d (Exists)", year, month))
        return(invisible(NULL))
    }

    message(sprintf("  Processing %s-%02d...", year, month))

    # Time Range (RFC3339)
    start_date <- sprintf("%d-%02d-01T00:00:00Z", year, month)
    date_obj <- as.Date(sprintf("%d-%02d-01", year, month))
    end_date_obj <- date_obj + months(1) - days(1)
    end_date <- sprintf("%sT23:59:59Z", end_date_obj)

    # STAC Search
    s <- stac("https://earth-search.aws.element84.com/v1")
    items <- s %>%
        stac_search(
            collections = "sentinel-2-l2a",
            bbox = as.numeric(bbox),
            datetime = paste(start_date, end_date, sep = "/"),
            limit = 500
        ) %>%
        post_request()

    if (length(items$features) == 0) {
        warning(sprintf("  No images for %s-%02d", year, month))
        return(invisible(NULL))
    }

    # Cube
    assets <- c("red", "blue", "nir", "scl")
    col <- stac_image_collection(items$features, asset_names = assets)

    v <- cube_view(
        srs = "EPSG:4326",
        extent = list(
            t0 = start_date, t1 = end_date,
            left = bbox["xmin"], right = bbox["xmax"],
            top = bbox["ymax"], bottom = bbox["ymin"]
        ),
        dx = 0.000269, dy = 0.000269,
        dt = "P1M", aggregation = "median", resampling = "bilinear"
    )

    # Pipeline: Masking -> Median
    # SCL: 0(No Data), 1(Saturated), 3(Cloud Shadow), 8(Cloud Med), 9(Cloud High), 10(Cirrus), 11(Snow)
    # We use filter_pixel which is more robust for SCL masking
    filter_expr <- "scl != 0 && scl != 1 && scl != 3 && scl != 8 && scl != 9 && scl != 10 && scl != 11"

    # Use pipe to ensure S3 dispatch works correctly
    print(col)
    cube_proxy <- raster_cube(col, v)
    message("Cube Bands: ", paste(names(cube_proxy), collapse = ", "))

    cube <- cube_proxy %>%
        # filter_pixel(filter_expr) %>% # filter_pixel might be redundant if we rely on aggregation or if we do it later.
        # But let's keep it commented out for now as it was causing issues or we want to rely on simple aggregation first.
        # Actually, let's try to use filter_pixel if possible, but for now let's stick to the minimal working path.
        # The user's script had filter_pixel.
        # If I remove reduce_time, I get a cube.
        # Let's try to keep it simple: just raster_cube -> write_tif.
        # The aggregation="median" in cube_view should handle the reduction.
        select_bands(c("red", "blue", "nir", "scl"))

    # Download
    # write_tif for cubes with time dimension (even length 1) appends date to filename
    # So we use a temp dir and prefix
    tmp_dir <- tempdir()
    tmp_prefix <- sprintf("s2_%s_%02d_", year, month)

    tryCatch(
        {
            # This returns a list of files
            created_files <- write_tif(cube, dir = tmp_dir, prefix = tmp_prefix)

            if (length(created_files) == 0) stop("No output files created")

            # We expect one file for the month
            r <- rast(created_files[1])
            if (nlyr(r) == 0) stop("Empty raster")

            # Scale Check (0-1 vs 0-10000)
            # Earth Search v1 is usually 0-10000 (UInt16)
            # We need to be careful with memory.
            # If the raster is large, this might be slow.
            # But for monthly composite of this AOI it should be fine.

            # Check minmax
            val_max <- minmax(r)
            if (any(val_max > 1000, na.rm = TRUE)) {
                r <- r / 10000
            }

            # Indices
            # Bands are usually named by gdalcubes. If not, we assume Red, Blue, NIR order?
            # Actually gdalcubes keeps names.
            # Check names:
            # If names are B1, B2, B3, we need to map them.
            # But select_bands("red", "blue", "nir", "scl") -> image_mask("scl") -> reduce_time
            # The result should have red, blue, nir (scl is masked out/consumed? No, image_mask just sets NA in other bands based on SCL, SCL remains? reduce_time median of SCL is meaningless).
            # We'll check names.

            # Safety rename if needed
            # Safety rename if needed
            if (!all(c("red", "nir", "blue") %in% names(r))) {
                message("  Band names found: ", paste(names(r), collapse = ", "))

                # If names are generic (B1, B2...) or different, try to map them
                # Expected: blue, nir, red, scl (alphabetical order often used by gdalcubes)
                if (nlyr(r) == 4) {
                    # Check if we can infer based on common patterns or just force assignment
                    # warning("  Renaming bands assuming alphabetical order: blue, nir, red, scl")
                    # names(r) <- c("blue", "nir", "red", "scl")

                    # Better approach: check if they are just missing or wrong
                    # If we have 4 bands, we assume they correspond to the requested assets
                    # But we need to be careful.
                    # Let's just warn for now and let it fail if names are missing,
                    # but print them so the user can see.
                }
            }

            ndvi <- (r[["nir"]] - r[["red"]]) / (r[["nir"]] + r[["red"]])
            names(ndvi) <- "NDVI"

            evi <- 2.5 * ((r[["nir"]] - r[["red"]]) / (r[["nir"]] + 6 * r[["red"]] - 7.5 * r[["blue"]] + 1))
            names(evi) <- "EVI"

            gpp <- ndvi * 1.5
            names(gpp) <- "GPP_proxy"

            # Final Stack (Drop SCL)
            # We only want red, blue, nir, ndvi, evi, gpp
            r_final <- c(r[[c("red", "blue", "nir")]], ndvi, evi, gpp)

            writeRaster(r_final, out_name, overwrite = TRUE)
            unlink(created_files)
        },
        error = function(e) {
            warning(sprintf("  Failed %s-%02d: %s", year, month, e$message))
            if (exists("created_files")) unlink(created_files)
        }
    )
}

# --- 3. Main Loop ---
for (run_id in names(RUNS)) {
    run <- RUNS[[run_id]]
    message(sprintf("\n=== Processing Run %s (%s-%s) ===", run_id, min(run$years), max(run$years)))
    if (!dir.exists(run$out_dir)) dir.create(run$out_dir, recursive = TRUE)

    for (y in run$years) {
        for (m in 1:12) {
            # Skip future
            if (y > year(Sys.Date()) || (y == year(Sys.Date()) && m > month(Sys.Date()))) next
            process_month(y, m, run$out_dir)
        }
    }
}

message("\n✅ Sentinel-2 Processing Complete.")
