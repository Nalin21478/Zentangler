#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

default_blueberry_negative_control_grid <- function() {
  data.frame(
    mae = c(
      "pathways",
      "ecs",
      "ecs",
      "pathways",
      "ecs",
      "all_omics",
      "all_omics",
      "ecs",
      "pathways"
    ),
    fusion_mode = c(
      "early",
      "late",
      "late",
      "late",
      "late",
      "late",
      "late",
      "late",
      "late"
    ),
    sis_n = c(25L, 100L, 10L, 25L, 25L, 25L, 25L, 25L, 10L),
    glmnet_alpha = c(1, 1, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1),
    outcome = c(
      "cGMP_change",
      "FMD_change",
      "Trigl_change",
      "HDLC_change",
      "DBP_change",
      "CRP_ug_ml20_change",
      "SBP_change",
      "IL6_pg_ml20_change",
      "LDLC_change"
    ),
    n_active_mediators = c(4L, 6L, 5L, 4L, 3L, 3L, 6L, 1L, 2L),
    stringsAsFactors = FALSE
  )
}

args <- commandArgs(trailingOnly = TRUE)
out_file <- if (length(args) >= 1L) args[[1L]] else "blueberry_parallel_negative_control_grid.csv"

grid <- default_blueberry_negative_control_grid()
grid$job_id <- seq_len(nrow(grid))

utils::write.csv(grid, out_file, row.names = FALSE)
message("Wrote negative-control grid to: ", out_file)
