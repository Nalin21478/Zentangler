#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

args <- commandArgs(trailingOnly = TRUE)
results_rdata <- if (length(args) >= 1L) args[[1L]] else "Parallel_CLR.RData"
out_file <- if (length(args) >= 2L) args[[2L]] else "blueberry_parallel_clr_nc_grid.csv"

env <- new.env(parent = emptyenv())
load(results_rdata, envir = env)

if (!exists("best_per_outcome", envir = env, inherits = FALSE)) {
  stop("Expected object 'best_per_outcome' in ", results_rdata, call. = FALSE)
}

grid <- as.data.frame(get("best_per_outcome", envir = env, inherits = FALSE), stringsAsFactors = FALSE)

needed <- c("mae", "fusion_mode", "sis_n", "glmnet_alpha", "outcome", "n_active_mediators")
missing_needed <- setdiff(needed, colnames(grid))
if (length(missing_needed) > 0L) {
  stop("best_per_outcome is missing columns: ", paste(missing_needed, collapse = ", "), call. = FALSE)
}

grid <- grid[, needed, drop = FALSE]
grid$job_id <- seq_len(nrow(grid))

utils::write.csv(grid, out_file, row.names = FALSE)
message("Wrote CLR negative-control grid to: ", out_file)
message("Rows: ", nrow(grid))
