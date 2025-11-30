#!/usr/bin/env Rscript
suppressPackageStartupMessages({
    library(optparse)
    library(dplyr)
    library(readr)
    library(lubridate)
    library(amt)
    library(sf)
    library(terra)
    library(survival)
    library(stringr)
    library(tidyr)
    library(tibble)
    library(ggplot2)
    library(scales)
})

# script location and project root
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
    stop("Cannot determine script path. Run via Rscript or set working directory to project root.")
}
.script <- dirname(.this_file())
.root <- normalizePath(file.path(.script, ".."), winslash = "/", mustWork = TRUE)

setwd(.root)

source(file.path(.root, "R", "utils_io.R"))
source(file.path(.root, "R", "utils_repro.R"))

opt <- list(
    make_option(c("--run"), type = "character", default = "A", help = "Run tag: A or B [default: %default]"),
    make_option(c("--mode"), type = "character", default = "REPRO", help = "REPRO or FAST (overrides config.yml)"),
    make_option(c("--species"), type = "character", default = NULL, help = "Optional single species (e.g., E3A)")
)
opts <- parse_args(OptionParser(option_list = opt))

# Validate run parameter
opts$run <- toupper({
    x <- opts$run
    if (is.null(x) || is.na(x) || !nzchar(x)) "A" else x
})
if (!opts$run %in% c("A", "B")) stop("`--run` must be A or B (got: ", opts$run, ")")

cfg <- read_config(file.path(.root, "config.yml"))
if (!is.null(opts$mode)) cfg$mode <- toupper(opts$mode)
mode <- set_mode(cfg)

# ===== Configuration =====
# Analysis parameters
tz_data <- "UTC"
seed_val <- 20161113
n_controls <- 10 # Number of random control steps per used step
min_burst_n <- 4 # Minimum GPS fixes per burst
min_steps_for_fit <- 20 # Minimum steps needed to estimate step-length distribution
drop_na_prop <- 0.3 # Drop variables with >30% missing data
near_zero_sd <- 1e-8 # Threshold for near-constant variables

# Plot settings
max_vars_partial <- 6 # Number of top variables for partial response plots
fig_width <- 9
fig_height <- 6
dpi <- 300

set.seed(seed_val)

# Define directories based on run
env_dir <- if (opts$run == "A") {
    cfg$paths$envi_after %||% file.path("BioHabs", "data", "envi", "A")
} else {
    cfg$paths$envi_before %||% file.path("BioHabs", "data", "envi", "B")
}
# Note: SSF uses 'stack' subdirectory in env_dir
env_stack_dir <- file.path(env_dir, "stack")

sp_dir <- if (opts$run == "A") {
    file.path("BioHabs", "data", "clean", "A")
} else {
    file.path("BioHabs", "data", "clean", "B")
}
# Fallback if specific A/B clean folders don't exist, try generic clean
if (!dir.exists(file.path(.root, sp_dir))) {
    sp_dir <- cfg$paths$clean %||% file.path("BioHabs", "data", "clean")
}

res_dir <- file.path("results", "iSSA", opts$run)

# Ensure absolute paths
env_stack_dir <- file.path(.root, env_stack_dir)
sp_dir <- file.path(.root, sp_dir)
res_dir <- file.path(.root, res_dir)
logs_dir <- cfg$paths$logs %||% file.path(.root, "logs")

dir_ensure(res_dir)
dir_ensure(logs_dir)

logf <- file.path(logs_dir, sprintf("01_ssf_to_rasters_%s.log", opts$run))
log_line(sprintf("Starting 01_ssf_to_rasters.R (mode=%s, run=%s)", mode, opts$run), logf)
log_line(sprintf("Input species dir: %s", sp_dir), logf)
log_line(sprintf("Env stack dir: %s", env_stack_dir), logf)
log_line(sprintf("Results dir: %s", res_dir), logf)

# ===== Helper Functions =====
req_cols <- c("collar_id", "species", "timestamp", "lon", "lat")

parse_ts <- function(x, tz = "UTC") {
    ts <- suppressWarnings(lubridate::dmy_hm(x, tz = tz))
    if (all(is.na(ts))) ts <- suppressWarnings(lubridate::ymd_hms(x, tz = tz))
    if (all(is.na(ts))) ts <- suppressWarnings(lubridate::ymd_hm(x, tz = tz))
    if (all(is.na(ts))) stop("Timestamp parsing failed. Check 'timestamp' format.")
    ts
}

# Generate control endpoints (uniform angle; gamma step length)
make_controls <- function(steps_used, n_controls, gamma_shape, gamma_scale) {
    need <- c("x1_", "y1_", "step_id_")
    if (!all(need %in% names(steps_used))) {
        stop(
            "Expected columns x1_, y1_, step_id_ in steps table. Have: ",
            paste(names(steps_used), collapse = ", ")
        )
    }
    N <- nrow(steps_used)
    if (N == 0L) {
        return(dplyr::slice_head(steps_used, n = 0))
    }

    base <- steps_used[rep(seq_len(N), each = n_controls), c("x1_", "y1_", "step_id_")]
    L <- rgamma(n = N * n_controls, shape = gamma_shape, scale = gamma_scale)
    ang <- runif(n = N * n_controls, min = 0, max = 2 * pi)

    tibble(
        step_id_ = base$step_id_,
        x2_ = base$x1_ + L * cos(ang),
        y2_ = base$y1_ + L * sin(ang),
        case_ = 0L
    )
}



pretty_var <- function(v) gsub("_", " ", v)

theme_pub <- function() {
    theme_minimal(base_size = 12) +
        theme(
            panel.grid.minor = element_blank(),
            plot.title = element_text(face = "bold"),
            axis.title = element_text(face = "bold")
        )
}

# ============================ DIAGNOSTIC HELPERS ==============================
# Calculate VIF manually (avoids 'car' dependency)
calc_vif <- function(model) {
    tryCatch(
        {
            mm <- model.matrix(model)
            # Remove intercept if present (clogit doesn't have one)
            if ("(Intercept)" %in% colnames(mm)) mm <- mm[, -which(colnames(mm) == "(Intercept)"), drop = FALSE]

            n_vars <- ncol(mm)
            if (n_vars < 2) {
                return(setNames(rep(NA, n_vars), colnames(mm)))
            }

            vifs <- numeric(n_vars)
            names(vifs) <- colnames(mm)

            for (i in 1:n_vars) {
                y <- mm[, i]
                x <- mm[, -i, drop = FALSE]
                r2 <- summary(lm(y ~ x))$r.squared
                vifs[i] <- 1 / (1 - r2)
            }
            vifs
        },
        error = function(e) {
            warning("VIF calculation failed: ", e$message)
            NULL
        }
    )
}

# k-fold Cross Validation for iSSA
run_cv <- function(data, formula, k = 5) {
    tryCatch(
        {
            # Split by strata (step_id_) to keep case/control groups together
            strata_ids <- unique(data$step_id_)
            n_strata <- length(strata_ids)
            if (n_strata < k) {
                return(NULL)
            }

            folds <- sample(rep(1:k, length.out = n_strata))
            cv_res <- numeric(k)

            for (i in 1:k) {
                test_strata <- strata_ids[folds == i]
                train_data <- data %>% filter(!step_id_ %in% test_strata)
                test_data <- data %>% filter(step_id_ %in% test_strata)

                # Fit model on training
                m <- survival::clogit(formula, data = train_data)

                # Predict on test (risk score)
                # Note: predict(..., type="risk") gives exp(Xb)
                # But we need to handle new levels in strata if any (though strata shouldn't matter for Xb)
                pred <- predict(m, newdata = test_data, type = "risk")

                # Calculate conditional log-likelihood for test set
                # For each stratum: P(case) = exp(Xb_case) / sum(exp(Xb_all))
                test_data$risk <- pred

                loglik <- test_data %>%
                    group_by(step_id_) %>%
                    summarize(
                        prob = risk[case_ == 1] / sum(risk),
                        .groups = "drop"
                    ) %>%
                    summarize(ll = sum(log(prob))) %>%
                    pull(ll)

                cv_res[i] <- loglik
            }

            list(
                mean_ll = mean(cv_res),
                se_ll = sd(cv_res) / sqrt(k),
                folds = k
            )
        },
        error = function(e) {
            warning("CV failed: ", e$message)
            NULL
        }
    )
}

# ============================ CORE PER-FILE PROCESS ===========================
process_one_file <- function(run_tag, stack_path, csv_file, out_base, i_seed = 0L) {
    tag <- tools::file_path_sans_ext(basename(csv_file))
    log_line(sprintf("Processing: %s", tag), logf)

    out_dir <- file.path(out_base, tag)
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

    # Local log for the run
    run_log_file <- file.path(out_dir, "run_log.txt")
    log_msg <- function(...) {
        txt <- paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", paste0(..., collapse = ""))
        cat(txt, "\n", file = run_log_file, append = TRUE)
    }

    tryCatch(
        {
            # ---------- Load & prepare stack ----------
            stopifnot(file.exists(stack_path))
            r_stack <- terra::rast(stack_path)
            stack_crs <- sf::st_crs(r_stack)

            # Save original names, then assign safe names based on original
            orig_names <- names(r_stack)

            # Rename numeric names 1-19 to Bio1-Bio19 (WorldClim convention)
            # This addresses the user's request to have correct variable names
            new_names <- orig_names
            for (i in seq_along(new_names)) {
                if (grepl("^[0-9]+$", new_names[i])) {
                    val <- as.integer(new_names[i])
                    if (val >= 1 && val <= 19) {
                        new_names[i] <- paste0("Bio", val)
                    }
                }
            }

            # Use make.names to ensure they are valid R variable names
            # This preserves "Bio1" as "Bio1", but changes "dist road" to "dist.road"
            # It also handles "0" -> "X0" if not renamed
            new_names <- make.names(new_names, unique = TRUE)
            names(r_stack) <- new_names
            safe_names <- names(r_stack)

            # Save name mapping (original -> safe) for traceability
            readr::write_csv(
                tibble(original = orig_names, safe = safe_names),
                file.path(out_dir, "stack_name_mapping.csv")
            )

            # ---------- Load data ----------
            stopifnot(file.exists(csv_file))
            dat <- readr::read_csv(csv_file, show_col_types = FALSE)

            miss <- setdiff(req_cols, names(dat))
            if (length(miss) > 0) stop("Missing required columns: ", paste(miss, collapse = ", "))

            dat$timestamp <- parse_ts(dat$timestamp, tz = tz_data)
            if (any(is.na(dat$lon) | is.na(dat$lat))) stop("Missing lon/lat values.")
            if (any(abs(dat$lon) > 180 | abs(dat$lat) > 90)) stop("Invalid lon/lat values.")

            dat <- dat %>%
                arrange(collar_id, timestamp) %>%
                group_by(collar_id) %>%
                filter(!(duplicated(timestamp) & duplicated(lon) & duplicated(lat))) %>%
                ungroup()

            pts_wgs <- sf::st_as_sf(dat, coords = c("lon", "lat"), crs = 4326)
            pts_m <- sf::st_transform(pts_wgs, crs = stack_crs)
            crd <- sf::st_coordinates(pts_m)
            pts_m$X <- crd[, 1]
            pts_m$Y <- crd[, 2]

            trk <- amt::make_track(pts_m, .x = X, .y = Y, .t = timestamp, id = collar_id)

            # ---------- resample 30±10 -> fallback 60±20 ----------
            trk30 <- trk %>%
                amt::track_resample(rate = lubridate::minutes(30), tolerance = lubridate::minutes(10)) %>%
                amt::filter_min_n_burst(min_n = min_burst_n)

            trk_use <- if (nrow(trk30) >= 2 * min_burst_n) {
                log_msg("Using cadence 30±10 minutes.")
                trk30
            } else {
                log_msg("Insufficient fixes for 30±10; trying 60±20.")
                trk60 <- trk %>%
                    amt::track_resample(rate = lubridate::minutes(60), tolerance = lubridate::minutes(20)) %>%
                    amt::filter_min_n_burst(min_n = min_burst_n)
                if (nrow(trk60) < 2 * min_burst_n) stop("Insufficient fixes for both cadences; skipping.")
                log_msg("Using cadence 60±20 minutes.")
                trk60
            }

            # ---------- used steps ----------
            steps_raw <- amt::steps_by_burst(trk_use)

            if (!"step_id_" %in% names(steps_raw)) {
                if ("burst_" %in% names(steps_raw)) {
                    steps_raw <- steps_raw %>%
                        group_by(burst_) %>%
                        mutate(step_id_ = paste0(as.character(burst_), "_", row_number())) %>%
                        ungroup()
                } else {
                    steps_raw$step_id_ <- sprintf("s%07d", seq_len(nrow(steps_raw)))
                }
                log_msg("Added missing step_id_ column.")
            }

            eps <- 1e-5
            if ("sl_" %in% names(steps_raw)) {
                zero_n <- sum(is.finite(steps_raw$sl_) & steps_raw$sl_ == 0)
                if (zero_n > 0) {
                    log_msg("Found ", zero_n, " zero-length steps -> replacing with epsilon.")
                    steps_raw$sl_[is.finite(steps_raw$sl_) & steps_raw$sl_ == 0] <- eps
                }
            }

            needed_cols <- c("x1_", "y1_", "x2_", "y2_")
            miss_cols <- setdiff(needed_cols, names(steps_raw))
            if (length(miss_cols) > 0) stop("Expected columns missing in steps_by_burst output: ", paste(miss_cols, collapse = ", "))

            # ---------- Gamma for controls ----------
            sl_vec <- steps_raw$sl_[is.finite(steps_raw$sl_) & steps_raw$sl_ > 0]
            if (length(sl_vec) < min_steps_for_fit) stop("Too few valid step lengths for param estimation.")
            sl_mean <- mean(sl_vec)
            sl_sd <- stats::sd(sl_vec)
            if (!is.finite(sl_sd) || sl_sd == 0) {
                log_msg("Zero SD for step lengths; inflating SD slightly.")
                sl_sd <- max(eps, sl_mean * 0.05)
            }
            gamma_shape <- (sl_mean / sl_sd)^2
            gamma_scale <- (sl_sd^2) / sl_mean

            used <- steps_raw %>%
                select(step_id_, x2_, y2_) %>%
                mutate(case_ = 1L)

            set.seed(seed_val + i_seed)
            controls <- make_controls(
                steps_used   = steps_raw,
                n_controls   = n_controls,
                gamma_shape  = gamma_shape,
                gamma_scale  = gamma_scale
            )

            ssf_data <- bind_rows(used, controls) %>%
                filter(is.finite(x2_), is.finite(y2_))
            if (nrow(ssf_data) < 50) stop("Too few steps for modeling (n<50).")

            # ---------- Extract ALL covariates at endpoints ----------
            xy_mat <- as.matrix(ssf_data[, c("x2_", "y2_")])
            ext <- terra::extract(r_stack, xy_mat, method = "bilinear")
            if ("ID" %in% names(ext)) ext <- ext[, -1, drop = FALSE]
            colnames(ext) <- names(r_stack)
            covar_names <- colnames(ext)
            ssf_data <- bind_cols(ssf_data, as_tibble(ext))

            # ---------- Filter covariates by NA and variance ----------
            na_rate <- sapply(ssf_data[, covar_names, drop = FALSE], function(z) mean(!is.finite(z)))
            keep1 <- names(na_rate)[na_rate <= drop_na_prop]
            sds <- sapply(ssf_data[, keep1, drop = FALSE], stats::sd, na.rm = TRUE)
            keep2 <- names(sds)[is.finite(sds) & sds > near_zero_sd]
            cov_keep <- intersect(keep1, keep2)
            if (length(cov_keep) < 2) stop("Too few usable covariates after NA/variance filtering.")

            cc <- stats::complete.cases(ssf_data[, cov_keep, drop = FALSE])
            ssf_data <- ssf_data[cc, , drop = FALSE]
            if (nrow(ssf_data) < 50) stop("Too few steps after covariate filtering.")

            # ---------- Z-score ----------
            stats_df <- tibble(
                variable = cov_keep,
                mean = sapply(ssf_data[, cov_keep, drop = FALSE], function(z) mean(z, na.rm = TRUE)),
                sd = sapply(ssf_data[, cov_keep, drop = FALSE], function(z) stats::sd(z, na.rm = TRUE))
            )
            stats_df$sd[!is.finite(stats_df$sd) | stats_df$sd == 0] <- 1

            for (v in stats_df$variable) {
                zname <- paste0(v, "_z")
                mu <- stats_df$mean[stats_df$variable == v]
                sdv <- stats_df$sd[stats_df$variable == v]
                ssf_data[[zname]] <- (ssf_data[[v]] - mu) / sdv
            }

            # ========== ADD MOVEMENT COVARIATES (iSSA) ==========
            # Step length and turn angle from steps_raw
            ssf_data <- ssf_data %>%
                left_join(
                    steps_raw %>% select(step_id_, sl_, ta_),
                    by = "step_id_"
                )

            # Create movement covariates (Signer et al. 2024)
            ssf_data <- ssf_data %>%
                mutate(
                    log_sl_ = log(sl_ + 1e-6), # Log step length (avoid log(0))
                    cos_ta_ = cos(ta_) # Cosine turn angle (directional persistence)
                ) %>%
                filter(is.finite(log_sl_), is.finite(cos_ta_))

            if (nrow(ssf_data) < 50) stop("Too few steps after adding movement covariates.")

            log_msg(sprintf("Added movement covariates: sl_, log_sl_, cos_ta_ (n=%d)", nrow(ssf_data)))

            # ---------- Model (TRUE iSSA: Movement + Habitat) ----------
            z_vars <- paste0(cov_keep, "_z")

            # iSSA formula: Movement + Habitat + Strata
            movement_terms <- c("sl_", "log_sl_", "cos_ta_")
            all_terms <- c(movement_terms, z_vars)

            base_fml <- stats::reformulate(termlabels = all_terms, response = "case_")
            fml <- update(base_fml, . ~ . + strata(step_id_))

            log_msg(sprintf(
                "Fitting iSSA with %d habitat + %d movement covariates",
                length(z_vars), length(movement_terms)
            ))
            model_clogit <- survival::clogit(fml, data = ssf_data)

            # ====== Hardened coefficient/CI assembly (name-aligned; no out-of-bounds) ======
            sm <- summary(model_clogit)
            betas <- coef(model_clogit)
            nm <- names(betas)

            # Coef table as data.frame, row-aligned to nm
            coef_mat <- if (is.null(dim(sm$coefficients))) {
                cm <- as.data.frame(t(sm$coefficients))
                rownames(cm) <- nm
                cm
            } else {
                as.data.frame(sm$coefficients)
            }
            # Ensure rows are in the same order as nm
            if (!is.null(rownames(coef_mat))) {
                coef_mat <- coef_mat[nm, , drop = FALSE]
            } else {
                rownames(coef_mat) <- nm
            }

            # Confidence intervals aligned to nm
            ci_raw <- suppressWarnings(confint(model_clogit))
            if (is.null(dim(ci_raw))) {
                ci_mat <- matrix(ci_raw, nrow = 1, dimnames = list(nm, c("2.5 %", "97.5 %")))
            } else {
                ci_mat <- as.matrix(ci_raw)
                if (!is.null(rownames(ci_mat))) {
                    ci_mat <- ci_mat[nm, , drop = FALSE]
                } else {
                    rownames(ci_mat) <- nm
                }
            }

            se <- if ("SE(coef)" %in% colnames(coef_mat)) coef_mat[, "SE(coef)"] else rep(NA_real_, length(nm))
            zval <- if ("z" %in% colnames(coef_mat)) coef_mat[, "z"] else rep(NA_real_, length(nm))
            pval <- if ("Pr(>|z|)" %in% colnames(coef_mat)) coef_mat[, "Pr(>|z|)"] else rep(NA_real_, length(nm))

            beta_tbl <- tibble(
                variable = sub("_z$", "", nm),
                beta     = as.numeric(betas),
                se       = as.numeric(se),
                z        = as.numeric(zval),
                p        = as.numeric(pval),
                ci_lo    = as.numeric(ci_mat[, 1]),
                ci_hi    = as.numeric(ci_mat[, 2])
            ) %>%
                left_join(stats_df, by = "variable") %>%
                select(variable, beta, se, z, p, ci_lo, ci_hi, mean, sd) %>%
                arrange(desc(abs(beta)))

            readr::write_csv(beta_tbl, file.path(out_dir, "model_coefficients_and_scaling.csv"))
            capture.output(sm, file = file.path(out_dir, "clogit_summary.txt"))

            # ---------- Enhanced Diagnostics ----------
            # 1. VIF
            vifs <- calc_vif(model_clogit)
            if (!is.null(vifs)) {
                vif_df <- tibble(variable = names(vifs), vif = vifs)
                readr::write_csv(vif_df, file.path(out_dir, "diagnostics_vif.csv"))
                log_msg("Calculated VIFs.")
            }

            # 2. Cross-Validation (5-fold)
            cv_res <- run_cv(ssf_data, fml, k = 5)
            if (!is.null(cv_res)) {
                cv_txt <- sprintf("5-fold CV Log-Likelihood: %.2f (SE: %.2f)", cv_res$mean_ll, cv_res$se_ll)
                writeLines(cv_txt, file.path(out_dir, "diagnostics_cv.txt"))
                log_msg(paste("CV Result:", cv_txt))
            }

            # ---------- Partial Response Plots (Habitat + Movement) ----------
            # Select top variables by absolute z-score
            top_vars <- beta_tbl %>%
                arrange(desc(abs(z))) %>%
                slice_head(n = max_vars_partial) %>%
                pull(variable)

            # Ensure movement vars are considered for plotting if significant
            # (They are already in beta_tbl, so they will be picked up if z-score is high)

            if (length(top_vars) > 0) {
                pdf(file.path(out_dir, "partial_response_plots.pdf"), width = fig_width, height = fig_height)

                for (v in top_vars) {
                    # Handle movement vs habitat variables differently for range
                    is_move <- v %in% c("sl_", "log_sl_", "cos_ta_")

                    # Construct range for prediction
                    if (is_move) {
                        # Use observed range from data
                        vals <- seq(min(ssf_data[[v]], na.rm = TRUE), max(ssf_data[[v]], na.rm = TRUE), length.out = 100)
                    } else {
                        # Habitat vars are z-scaled, so -3 to 3 is reasonable
                        vals <- seq(-3, 3, length.out = 100)
                    }

                    # Create prediction data frame
                    # Set all other vars to their mean (0 for z-scores, mean for movement)
                    pred_df <- ssf_data[rep(1, 100), ]
                    # Reset all predictors to mean/0
                    for (vn in all_terms) {
                        if (vn %in% movement_terms) {
                            pred_df[[vn]] <- mean(ssf_data[[vn]], na.rm = TRUE)
                        } else {
                            pred_df[[vn]] <- 0
                        }
                    }
                    # Vary the target variable
                    pred_df[[v]] <- vals

                    # Predict relative selection strength (RSS)
                    # exp(beta * x) ignoring intercept/strata
                    beta_val <- beta_tbl$beta[beta_tbl$variable == v]
                    # Note: This is a simplified partial plot (marginal effect)
                    # RSS = exp(beta * (x - mean))
                    # Since we centered x around mean/0, it's just exp(beta * x) relative to x=0/mean

                    # Actually, for correct RSS, we compare x to a reference x0.
                    # Let's use the mean as reference.
                    ref_val <- if (is_move) mean(ssf_data[[v]], na.rm = TRUE) else 0
                    rss <- exp(beta_val * (vals - ref_val))

                    p <- ggplot(data.frame(x = vals, y = rss), aes(x, y)) +
                        geom_line(color = "blue", size = 1) +
                        labs(
                            title = paste("Partial Response:", pretty_var(v)),
                            subtitle = sprintf("Beta = %.3f (z = %.1f)", beta_val, beta_tbl$z[beta_tbl$variable == v]),
                            x = if (is_move) v else paste(v, "(z-score)"),
                            y = "Relative Selection Strength (RSS)"
                        ) +
                        theme_pub()
                    print(p)
                }
                dev.off()
                log_msg("Generated partial response plots.")
            }

            # ---------- RSF projection (guard missing bands) ----------
            # Filter for projection ONLY, preserving beta_tbl for plots
            avail <- intersect(beta_tbl$variable, names(r_stack))
            if (length(avail) == 0) stop("No model covariates found in raster stack for projection.")

            # Use a subset for projection calculations
            beta_proj <- beta_tbl %>% filter(variable %in% avail)
            sel <- r_stack[[avail]]

            mu_vec <- beta_proj$mean
            sd_vec <- beta_proj$sd
            sd_vec[!is.finite(sd_vec) | sd_vec == 0] <- 1
            b_vec <- beta_proj$beta

            rsf <- terra::app(sel, fun = function(v) {
                z <- (v - mu_vec) / sd_vec
                exp(sum(z * b_vec))
            })

            g <- terra::global(rsf, fun = c("min", "max"), na.rm = TRUE)
            rsf01 <- if (isTRUE(all.equal(g[1, 1], g[1, 2]))) rsf * 0 else (rsf - g[1, 1]) / (g[1, 2] - g[1, 1])

            out_tif <- file.path(out_dir, paste0(tag, "_SSF_rsf.tif"))
            out_tif01 <- file.path(out_dir, paste0(tag, "_SSF_rsf_0to1.tif"))
            terra::writeRaster(rsf, out_tif, overwrite = TRUE, wopt = list(datatype = "FLT4S", gdal = "COMPRESS=LZW"))
            terra::writeRaster(rsf01, out_tif01, overwrite = TRUE, wopt = list(datatype = "FLT4S", gdal = "COMPRESS=LZW"))

            # ===================== DIAGNOSTIC PLOTS ===================
            if (nrow(beta_tbl) > 0) {
                coef_df <- beta_tbl %>%
                    mutate(variable_label = pretty_var(variable)) %>%
                    arrange(beta) %>%
                    mutate(variable_label = factor(variable_label, levels = variable_label))

                p_coef <- ggplot(coef_df, aes(x = beta, y = variable_label)) +
                    geom_vline(xintercept = 0, linetype = 2) +
                    geom_errorbarh(aes(xmin = ci_lo, xmax = ci_hi), width = 0) +
                    geom_point(size = 2) +
                    labs(
                        title = paste0("SSF coefficients (", run_tag, " — ", tag, ")"),
                        x = "Coefficient (log-relative selection)", y = NULL
                    ) +
                    theme_pub()
                ggsave(file.path(out_dir, "coef_forest.png"), p_coef, width = fig_width, height = fig_height, dpi = dpi)

                p_vi <- ggplot(coef_df, aes(x = reorder(variable_label, abs(beta)), y = abs(beta))) +
                    geom_col() +
                    coord_flip() +
                    labs(
                        title = paste0("Variable importance (|β|) — ", run_tag, " — ", tag),
                        x = NULL, y = "|β|"
                    ) +
                    theme_pub()
                ggsave(file.path(out_dir, "var_importance.png"), p_vi, width = fig_width, height = fig_height, dpi = dpi)

                top_vars <- head(coef_df$variable, min(10, nrow(coef_df)))
                # Only plot densities for habitat variables (z-scores)
                # Movement vars don't have z-scores in the same way in ssf_data (they are raw sl_, log_sl_, cos_ta_)
                # But we named them with _z suffix in the model? No, movement vars are raw.
                # Let's filter top_vars to only those present in ssf_data with _z suffix or exact match

                # Check which top vars are habitat (have _z suffix in data)
                hab_vars <- top_vars[paste0(top_vars, "_z") %in% names(ssf_data)]

                if (length(hab_vars) > 0) {
                    dens_df <- ssf_data %>%
                        select(case_, all_of(paste0(hab_vars, "_z"))) %>%
                        mutate(case_ = factor(case_, levels = c(0, 1), labels = c("available", "used"))) %>%
                        pivot_longer(cols = -case_, names_to = "variable", values_to = "value") %>%
                        mutate(
                            variable = sub("_z$", "", variable),
                            variable = factor(pretty_var(variable), levels = pretty_var(hab_vars))
                        )

                    p_dens <- ggplot(dens_df, aes(x = value, fill = case_)) +
                        geom_density(alpha = 0.5) +
                        facet_wrap(~variable, scales = "free", ncol = 3) +
                        labs(
                            title = paste0("Endpoint covariate densities (z) — ", run_tag, " — ", tag),
                            x = "z-score", y = "Density", fill = NULL
                        ) +
                        theme_pub()
                    ggsave(file.path(out_dir, "covariate_densities.png"), p_dens, width = fig_width, height = fig_height, dpi = dpi)
                }

                # (Old partial response plot removed - superseded by new PDF)

                # Win rate
                # Construct predictor matrix (handle movement vs habitat var names)
                vars <- coef_df$variable
                col_names <- sapply(vars, function(v) {
                    if (v %in% c("sl_", "log_sl_", "cos_ta_")) v else paste0(v, "_z")
                })

                # Check if all columns exist
                missing_cols <- setdiff(col_names, names(ssf_data))
                if (length(missing_cols) > 0) stop("Missing columns for Win Rate: ", paste(missing_cols, collapse = ", "))

                Xz <- as.matrix(ssf_data[, col_names, drop = FALSE])
                lp <- drop(Xz %*% matrix(coef_df$beta, ncol = 1))
                pred_rsf <- exp(lp)

                wr_df <- tibble(step_id_ = ssf_data$step_id_, case_ = ssf_data$case_, rsf = pred_rsf) %>%
                    group_by(step_id_) %>%
                    summarise(
                        used_rsf = rsf[case_ == 1][1],
                        max_rsf = max(rsf),
                        win = as.integer(used_rsf >= max_rsf),
                        .groups = "drop"
                    )
                win_rate <- mean(wr_df$win, na.rm = TRUE)
                readr::write_csv(wr_df, file.path(out_dir, "winrate_per_step.csv"))

                p_wr <- ggplot(wr_df, aes(x = win)) +
                    geom_bar() +
                    scale_x_continuous(breaks = c(0, 1), labels = c("lost", "won")) +
                    labs(
                        title = paste0("Within-stratum win rate = ", percent(win_rate), " — ", run_tag, " — ", tag),
                        x = NULL, y = "Count of strata"
                    ) +
                    theme_pub()
                ggsave(file.path(out_dir, "winrate_bar.png"), p_wr, width = 6, height = 4, dpi = dpi)
            }

            # ---------- Meta ----------
            meta <- tibble(
                run = run_tag,
                tag, input_csv = csv_file, n_points = nrow(dat),
                n_steps_used = sum(ssf_data$case_ == 1L),
                n_controls = n_controls,
                cadence = if (exists("trk30") && identical(trk_use, trk30)) "30±10" else "60±20",
                timestamp_tz = tz_data, seed = seed_val + i_seed,
                raster_stack = stack_path, output_rsf = out_tif, output_rsf01 = out_tif01,
                n_covariates = nrow(beta_tbl),
                covariates = paste(beta_tbl$variable, collapse = ", "),
                run_time = format(Sys.time(), "%Y-%m-%d %H:%M:%S")
            )
            if (exists("wr_df")) meta$win_rate <- mean(wr_df$win, na.rm = TRUE)

            readr::write_csv(meta, file.path(out_dir, "run_meta.csv"))
            log_line(sprintf("SUCCESS → %s (+ 0–1, plots, meta)", out_tif), logf)
        },
        error = function(e) {
            log_msg("ERROR: ", conditionMessage(e))
            log_line(sprintf("ERROR in %s: %s", tag, conditionMessage(e)), logf)
        }
    )
}

# ============================ MAIN EXECUTION ===============================
stack_path <- file.path(env_stack_dir, "stack.tif")
if (!file.exists(stack_path)) {
    log_line(sprintf("Missing stack: %s — skipping run.", stack_path), logf)
    quit(status = 0)
}

csv_files <- list.files(sp_dir, pattern = "\\.csv$", full.names = TRUE)
if (length(csv_files) == 0) {
    log_line(sprintf("No CSV files in: %s — skipping run.", sp_dir), logf)
    quit(status = 0)
}

# Filter by species if requested
if (!is.null(opts$species)) {
    csv_files <- csv_files[grepl(opts$species, basename(csv_files))]
    if (length(csv_files) == 0) {
        log_line(sprintf("Species %s not found in %s", opts$species, sp_dir), logf)
        quit(status = 0)
    }
}

for (i in seq_along(csv_files)) {
    process_one_file(opts$run, stack_path, csv_files[[i]], res_dir, i_seed = i)
}

log_line("All species finished.", logf)
