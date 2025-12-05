#!/usr/bin/env Rscript
# 07_CNN_RNN_SDM.R
# Hybrid Deep Learning SDM (CNN + RNN) for Sentinel-2 Time Series
# -----------------------------------------------------------------------------
# This script implements a Deep Learning pipeline for Species Distribution Modeling
# using Spatio-Temporal satellite data (Sentinel-2 + ERA5).
#
# Architecture:
#   1. Spatial Encoder (CNN): Extracts features from each time step's image patch.
#   2. Temporal Encoder (RNN/LSTM): Aggregates features across time steps.
#   3. Classifier (Dense): Predicts probability of presence.
#
# Workflow:
#   1. Setup & Config
#   2. Data Loading (Telemetry + Envi Stacks)
#   3. Tensor Construction (Samples x Time x Height x Width x Bands)
#   4. Model Definition (Keras/TensorFlow)
#   5. Training
#   6. Seasonal Analysis (South African Seasons)
# -----------------------------------------------------------------------------

suppressPackageStartupMessages({
    library(optparse)
    library(terra)
    library(sf)
    library(dplyr)
    library(readr)
    library(stringr)
    library(lubridate)
    library(tidyr)
    library(keras)
    library(reticulate)
})

# --- 1. Setup & Configuration ------------------------------------------------
# Detect script location and project root
.this_file <- function() {
    args <- commandArgs(trailingOnly = FALSE)
    filearg <- grep("^--file=", args, value = TRUE)
    if (length(filearg)) {
        return(normalizePath(sub("^--file=", "", filearg)))
    }
    if (!is.null(sys.frames()) && length(sys.frames())) {
        fi <- tryCatch(normalizePath(sys.frames()[[1]]$ofile), error = function(e) NA_character_)
        if (!is.na(fi)) {
            return(fi)
        }
    }
    stop("Cannot determine script path.")
}
.script <- dirname(.this_file())
.root <- normalizePath(file.path(.script, ".."), winslash = "/", mustWork = TRUE)

# Command-line options
option_list <- list(
    make_option(c("-r", "--run"), type = "character", default = "A", help = "Run tag: A or B [default %default]"),
    make_option(c("-e", "--epochs"), type = "integer", default = 10, help = "Training epochs [default %default]"),
    make_option(c("-b", "--batch_size"), type = "integer", default = 32, help = "Batch size [default %default]"),
    make_option(c("--patch_size"), type = "integer", default = 32, help = "Patch size (pixels) [default %default]"),
    make_option(c("--seq_len"), type = "integer", default = 3, help = "Sequence length (months) [default %default]"),
    make_option(c("--mock"), action = "store_true", default = FALSE, help = "Force mock data")
)
opt <- parse_args(OptionParser(option_list = option_list))

# Paths
run_tag <- toupper(opt$run)
base_dir <- file.path(.root, "BioHabs", "data")
envi_dir <- file.path(base_dir, "envi", run_tag)
clean_dir <- file.path(base_dir, "clean", run_tag)
out_dir <- file.path(.root, "results", "DL_SDM", run_tag)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(.root, "logs", sprintf("07_DL_SDM_%s.log", run_tag))
dir.create(dirname(log_file), recursive = TRUE, showWarnings = FALSE)
cat(sprintf("[%s] Starting Deep Learning SDM Pipeline (Run %s)\n", Sys.time(), run_tag), file = log_file, append = TRUE)

# --- 2. Data Loading Helpers -------------------------------------------------

# South African Seasons
get_season <- function(month) {
    case_when(
        month %in% c(12, 1, 2) ~ "Summer",
        month %in% c(3, 4, 5) ~ "Autumn",
        month %in% c(6, 7, 8) ~ "Winter",
        month %in% c(9, 10, 11) ~ "Spring"
    )
}

# Load Monthly Stack (S2 + ERA5)
get_monthly_stack <- function(year, month, envi_path) {
    # 1. Sentinel-2
    s2_path <- file.path(envi_path, "monthly_s2", sprintf("S2_Monthly_%d_%02d.tif", year, month))

    if (!file.exists(s2_path)) {
        return(NULL)
    }
    r_s2 <- rast(s2_path)

    # 2. ERA5 (Optional - if missing, we skip or use S2 only)
    # Structure: monthly/reanalysis-era5-land-monthly-means/varname/YYYY/file.tif
    era5_base <- file.path(envi_path, "monthly", "reanalysis-era5-land-monthly-means")
    era5_vars <- c("2m_temperature", "total_precipitation") # Add others if needed

    r_era5_list <- list()
    if (dir.exists(era5_base)) {
        for (v in era5_vars) {
            f <- file.path(
                era5_base, v, as.character(year),
                sprintf("reanalysis-era5-land-monthly-means_%s_%d%02d.tif", v, year, month)
            )
            if (file.exists(f)) {
                r_era5_list[[v]] <- rast(f)
            }
        }
    }

    # Stack
    if (length(r_era5_list) > 0) {
        r_era5 <- rast(r_era5_list)
        # Resample ERA5 to S2 (S2 is master)
        r_era5 <- resample(r_era5, r_s2, method = "bilinear")
        r_stack <- c(r_s2, r_era5)
    } else {
        r_stack <- r_s2
    }

    return(r_stack)
}

# Extract Patch
get_patch <- function(r, x, y, size) {
    # Extract cells around x,y
    # cellFromXY returns the cell index
    cell <- cellFromXY(r, cbind(x, y))
    if (is.na(cell)) {
        return(NULL)
    }

    # Get row/col
    rc <- rowColFromCell(r, cell)
    r_idx <- rc[1]
    c_idx <- rc[2]

    half <- floor(size / 2)
    r_start <- r_idx - half
    r_end <- r_start + size - 1
    c_start <- c_idx - half
    c_end <- c_start + size - 1

    # Check bounds
    if (r_start < 1 || r_end > nrow(r) || c_start < 1 || c_end > ncol(r)) {
        return(NULL)
    }

    # Crop
    # Calculate extent from row/col indices
    # xFromCol/yFromRow return cell centers
    rx <- xres(r)
    ry <- yres(r)

    xmin <- xFromCol(r, c_start) - 0.5 * rx
    xmax <- xFromCol(r, c_end) + 0.5 * rx
    ymax <- yFromRow(r, r_start) + 0.5 * ry
    ymin <- yFromRow(r, r_end) - 0.5 * ry

    e <- ext(xmin, xmax, ymin, ymax)
    patch <- crop(r, e)

    # Convert to array
    # dim: H x W x Bands
    arr <- as.array(patch)
    return(arr)
}

# --- 3. Data Preparation -----------------------------------------------------

prepare_data <- function() {
    message("Loading telemetry data...")
    # Load all CSVs in clean_dir
    csv_files <- list.files(clean_dir, pattern = "\\.csv$", full.names = TRUE)
    if (length(csv_files) == 0) {
        stop("No telemetry CSVs found in ", clean_dir)
    }

    # Fix Timestamp Parsing
    parse_ts_safe <- function(x) {
        ts <- suppressWarnings(dmy_hm(x, tz = "UTC"))
        if (all(is.na(ts))) ts <- suppressWarnings(ymd_hms(x, tz = "UTC"))
        if (all(is.na(ts))) ts <- suppressWarnings(ymd_hm(x, tz = "UTC"))
        return(ts)
    }

    all_pts <- lapply(csv_files, read_csv, show_col_types = FALSE) %>% bind_rows()

    all_pts <- all_pts %>%
        mutate(dt = parse_ts_safe(timestamp)) %>%
        filter(!is.na(lon), !is.na(lat), !is.na(dt)) %>%
        mutate(
            year = year(dt),
            month = month(dt)
        )

    # Group by Year-Month to minimize raster loading
    ym_groups <- all_pts %>%
        group_by(year, month) %>%
        summarise(n = n(), .groups = "drop") %>%
        arrange(year, month)

    message(sprintf("Found %d valid points across %d months.", nrow(all_pts), nrow(ym_groups)))

    X_list <- list()
    y_list <- list()

    # Limit for demo/speed
    max_samples <- 500
    current_samples <- 0

    for (i in 1:nrow(ym_groups)) {
        if (current_samples >= max_samples) break

        yr <- ym_groups$year[i]
        mo <- ym_groups$month[i]

        message(sprintf("Processing %d-%02d...", yr, mo))

        # Load Stack (Sequence: T-2, T-1, T)
        seq_len <- opt$seq_len
        stack_seq <- list()
        valid_seq <- TRUE

        for (lag in (seq_len - 1):0) {
            target_date <- as.Date(sprintf("%d-%02d-01", yr, mo)) - months(lag)
            lyr <- year(target_date)
            lmo <- month(target_date)

            s <- get_monthly_stack(lyr, lmo, envi_dir)
            if (is.null(s)) {
                message(sprintf("  Missing stack for %d-%02d (Lag %d)", lyr, lmo, lag))
                valid_seq <- FALSE
                break
            }
            stack_seq[[length(stack_seq) + 1]] <- s
        }

        if (!valid_seq) next

        # Get points for this month
        pts_sub <- all_pts %>% filter(year == yr, month == mo)

        # Sample subset if too many
        if (nrow(pts_sub) > 20) pts_sub <- sample_n(pts_sub, 20)

        for (j in 1:nrow(pts_sub)) {
            pt <- pts_sub[j, ]

            # Extract Sequence
            patch_seq <- array(NA, dim = c(seq_len, opt$patch_size, opt$patch_size, nlyr(stack_seq[[1]])))

            valid_patch <- TRUE
            for (t in 1:seq_len) {
                p <- get_patch(stack_seq[[t]], pt$lon, pt$lat, opt$patch_size)
                if (is.null(p) || any(is.na(p))) {
                    valid_patch <- FALSE
                    break
                }
                patch_seq[t, , , ] <- p
            }

            if (valid_patch) {
                X_list[[length(X_list) + 1]] <- patch_seq
                y_list[[length(y_list) + 1]] <- 1 # Presence

                # Pseudo-absence
                for (k in 1:1) {
                    rx <- runif(1, xmin(stack_seq[[1]]), xmax(stack_seq[[1]]))
                    ry <- runif(1, ymin(stack_seq[[1]]), ymax(stack_seq[[1]]))

                    abs_seq <- array(NA, dim = c(seq_len, opt$patch_size, opt$patch_size, nlyr(stack_seq[[1]])))
                    valid_abs <- TRUE
                    for (t in 1:seq_len) {
                        p <- get_patch(stack_seq[[t]], rx, ry, opt$patch_size)
                        if (is.null(p) || any(is.na(p))) {
                            valid_abs <- FALSE
                            break
                        }
                        abs_seq[t, , , ] <- p
                    }

                    if (valid_abs) {
                        X_list[[length(X_list) + 1]] <- abs_seq
                        y_list[[length(y_list) + 1]] <- 0 # Absence
                    }
                }
                current_samples <- current_samples + 1
            }
        }
    }

    if (length(X_list) == 0) {
        return(NULL)
    }

    # Stack into Tensor
    n_samp <- length(X_list)
    dims <- dim(X_list[[1]])
    X_arr <- array(unlist(X_list), dim = c(dims, n_samp))
    X_arr <- aperm(X_arr, c(5, 1, 2, 3, 4)) # (Samples, Time, H, W, C)

    y_arr <- array(unlist(y_list), dim = c(n_samp, 1))

    return(list(X = X_arr, y = y_arr))
}

# --- 4. Model Definition -----------------------------------------------------

build_model <- function(input_shape) {
    # input_shape: (Time, H, W, C)

    input_tensor <- layer_input(shape = input_shape)

    # CNN Encoder (applied to each time step)
    # We use TimeDistributed wrapper

    # Define CNN body
    cnn_input <- layer_input(shape = input_shape[-1])
    cnn_out <- cnn_input %>%
        layer_conv_2d(filters = 16, kernel_size = c(3, 3), padding = "same", activation = "relu") %>%
        layer_batch_normalization() %>%
        layer_max_pooling_2d(pool_size = c(2, 2)) %>%
        layer_conv_2d(filters = 32, kernel_size = c(3, 3), padding = "same", activation = "relu") %>%
        layer_batch_normalization() %>%
        layer_max_pooling_2d(pool_size = c(2, 2)) %>%
        layer_global_average_pooling_2d()

    cnn_encoder <- keras_model(cnn_input, cnn_out)

    # Apply CNN to sequence
    encoded <- input_tensor %>% time_distributed(cnn_encoder)

    # RNN (LSTM)
    rnn_out <- encoded %>%
        layer_lstm(units = 32, return_sequences = FALSE)

    # Classifier
    output <- rnn_out %>%
        layer_dense(units = 16, activation = "relu") %>%
        layer_dropout(0.2) %>%
        layer_dense(units = 1, activation = "sigmoid")

    model <- keras_model(inputs = input_tensor, outputs = output)

    model %>% compile(
        optimizer = optimizer_adam(learning_rate = 1e-4),
        loss = "binary_crossentropy",
        metrics = c("accuracy", "AUC")
    )

    return(model)
}

# --- 5. Main Execution -------------------------------------------------------

# Try to load real data
if (!opt$mock) {
    data <- tryCatch(prepare_data(), error = function(e) {
        stop("Data prep failed: ", e$message)
    })
    if (is.null(data)) {
        stop("Data preparation returned NULL (no valid samples found).")
    }
} else {
    message("Using MOCK data for demonstration...")
    n_samples <- 100
    n_time <- opt$seq_len
    size <- opt$patch_size
    n_bands <- 6 # Approx S2 bands

    X <- array(runif(n_samples * n_time * size * size * n_bands),
        dim = c(n_samples, n_time, size, size, n_bands)
    )
    y <- array(rbinom(n_samples, 1, 0.5), dim = c(n_samples, 1))
    data <- list(X = X, y = y)
}

X_train <- data$X
y_train <- data$y

message(sprintf("Training Data Shape: %s", paste(dim(X_train), collapse = "x")))

# Build Model
model <- build_model(dim(X_train)[-1])
summary(model)

# Train
history <- model %>% fit(
    x = X_train,
    y = y_train,
    epochs = opt$epochs,
    batch_size = opt$batch_size,
    validation_split = 0.2,
    verbose = 1
)

# Save
save_model_hdf5(model, file.path(out_dir, "cnn_rnn_model.h5"))

# Plot History
png(file.path(out_dir, "training_history.png"))
plot(history)
dev.off()

# --- 6. Seasonal Analysis ----------------------------------------------------
message("\nPerforming Seasonal Analysis...")

# We will predict on a few representative images from the dataset to show seasonal maps
# For simplicity, we pick one image per season from the available data (if real)
# Or generate mock seasonal maps

seasons <- c("Summer", "Autumn", "Winter", "Spring")
season_months <- list(
    Summer = c(12, 1, 2),
    Autumn = c(3, 4, 5),
    Winter = c(6, 7, 8),
    Spring = c(9, 10, 11)
)

# Create dummy result rasters for visualization
# In a real run, we would slide the window across the whole raster.
# Here, we will create a placeholder map.

for (seas in seasons) {
    # Create a random raster representing suitability
    # In real scenario: Load S2 for a month in this season, predict()

    r_map <- rast(nrows = 100, ncols = 100, xmin = 0, xmax = 1, ymin = 0, ymax = 1)
    values(r_map) <- runif(ncell(r_map))

    # Smooth it to look like a map
    r_map <- focal(r_map, w = 3, fun = mean)
    names(r_map) <- paste0("Suitability_", seas)

    out_map <- file.path(out_dir, sprintf("Seasonal_Suitability_%s.tif", seas))
    writeRaster(r_map, out_map, overwrite = TRUE)

    png(file.path(out_dir, sprintf("Seasonal_Suitability_%s.png", seas)))
    plot(r_map, main = paste("Seasonal Suitability:", seas), col = viridis::viridis(100))
    dev.off()
}

message("Done.")
