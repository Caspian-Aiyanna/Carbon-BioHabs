suppressPackageStartupMessages({
  library(yaml)
})

# Read the project configuration file
read_config <- function(path = "config.yml") {
  yaml::read_yaml(path)
}

# Set execution mode (REPRO or FAST) and configure threading
# REPRO mode: single-threaded for deterministic results
# FAST mode: multi-threaded for faster execution
set_mode <- function(cfg) {
  mode <- toupper(cfg$mode %||% "REPRO")
  if (mode == "REPRO") {
    # Force single-threaded execution for reproducibility
    Sys.setenv(
      "OMP_NUM_THREADS" = "1",
      "MKL_NUM_THREADS" = "1",
      "OPENBLAS_NUM_THREADS" = "1"
    )
  }
  invisible(mode)
}

# Default value operator: return b if a is NULL
`%||%` <- function(a, b) if (is.null(a)) b else a

# Generate a deterministic seed from a text tag
# This ensures the same tag always produces the same seed
seed_for <- function(tag, base = 1L) {
  h <- sum(utf8ToInt(as.character(tag))) + as.integer(base)
  h <- as.integer(abs(h) %% .Machine$integer.max)
  set.seed(h)
  invisible(h)
}

save_or_load <- function(path, value_expr) {
  if (file.exists(path)) {
    return(readr::read_rds(path))
  } else {
    val <- eval(substitute(value_expr))
    dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
    readr::write_rds(val, path)
    return(val)
  }
}

skip_if_done <- function(dir, files = character()) {
  all(file.exists(file.path(dir, files)))
}
