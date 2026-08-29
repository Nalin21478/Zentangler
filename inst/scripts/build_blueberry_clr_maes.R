#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(MultiAssayExperiment)
  library(S4Vectors)
  library(SummarizedExperiment)
})

subset_blueberry_mae_views <- function(mae, view_names) {
  all_views <- names(MultiAssayExperiment::experiments(mae))
  missing_views <- setdiff(view_names, all_views)
  if (length(missing_views) > 0) {
    stop("Views not found in MAE: ", paste(missing_views, collapse = ", "), call. = FALSE)
  }

  exps <- MultiAssayExperiment::experiments(mae)[view_names]
  sm <- as.data.frame(MultiAssayExperiment::sampleMap(mae), stringsAsFactors = FALSE)
  sm <- sm[sm$assay %in% view_names, , drop = FALSE]
  sm$assay <- factor(sm$assay, levels = view_names)

  MultiAssayExperiment::MultiAssayExperiment(
    experiments = MultiAssayExperiment::ExperimentList(exps),
    colData = S4Vectors::DataFrame(MultiAssayExperiment::colData(mae)),
    sampleMap = S4Vectors::DataFrame(
      assay = sm$assay,
      primary = sm$primary,
      colname = sm$colname
    )
  )
}

save_mae_pair <- function(mae, base_path_no_ext, object_name = "mae") {
  saveRDS(mae, paste0(base_path_no_ext, ".rds"))
  assign(object_name, mae, envir = environment())
  save(list = object_name, file = paste0(base_path_no_ext, ".RData"), envir = environment())
  invisible(TRUE)
}

compute_default_pseudocount <- function(mat) {
  positive_vals <- mat[is.finite(mat) & mat > 0]
  if (length(positive_vals) == 0L) return(1e-06)
  min(positive_vals) / 2
}

clr_transform_matrix <- function(mat, pseudocount = NULL) {
  mat <- as.matrix(mat)
  storage.mode(mat) <- "numeric"
  if (any(!is.finite(mat), na.rm = TRUE)) {
    stop("CLR transform requires finite values in the assay matrix.", call. = FALSE)
  }
  if (any(mat < 0, na.rm = TRUE)) {
    stop("CLR transform requires non-negative assay values.", call. = FALSE)
  }

  if (is.null(pseudocount)) {
    pseudocount <- compute_default_pseudocount(mat)
  }
  if (!is.finite(pseudocount) || pseudocount <= 0) {
    stop("Pseudocount must be a positive finite number.", call. = FALSE)
  }

  mat_pc <- mat + pseudocount
  log_mat <- log(mat_pc)
  clr_mat <- sweep(log_mat, 2, colMeans(log_mat, na.rm = TRUE), FUN = "-")
  attr(clr_mat, "pseudocount") <- pseudocount
  clr_mat
}

clr_transform_se <- function(se, assay_name = 1L, pseudocount = NULL) {
  assay_name <- SummarizedExperiment::assayNames(se)[[assay_name]]
  mat <- SummarizedExperiment::assay(se, assay_name)
  clr_mat <- clr_transform_matrix(mat, pseudocount = pseudocount)
  SummarizedExperiment::assay(se, assay_name) <- clr_mat
  S4Vectors::metadata(se)$transform <- "CLR"
  S4Vectors::metadata(se)$clr_pseudocount <- attr(clr_mat, "pseudocount")
  se
}

build_blueberry_clr_maes <- function(
  input_mae_path = file.path("data", "Blueberry_MAE.rds")
) {
  if (!file.exists(input_mae_path)) {
    stop("Input MAE not found: ", input_mae_path, call. = FALSE)
  }

  mae_raw <- readRDS(input_mae_path)
  exp_names <- names(MultiAssayExperiment::experiments(mae_raw))

  transformed_experiments <- lapply(exp_names, function(view_name) {
    se <- mae_raw[[view_name]]
    clr_transform_se(se)
  })
  names(transformed_experiments) <- exp_names

  mae_clr <- MultiAssayExperiment::MultiAssayExperiment(
    experiments = MultiAssayExperiment::ExperimentList(transformed_experiments),
    colData = S4Vectors::DataFrame(MultiAssayExperiment::colData(mae_raw)),
    sampleMap = S4Vectors::DataFrame(MultiAssayExperiment::sampleMap(mae_raw))
  )
  S4Vectors::metadata(mae_clr)$transform <- "CLR"

  maes <- list(
    full = mae_clr,
    pathways = subset_blueberry_mae_views(
      mae_clr,
      c("species", "dna_pathways", "rna_pathways")
    ),
    ecs = subset_blueberry_mae_views(
      mae_clr,
      c("species", "dna_ecs", "rna_ecs")
    ),
    all_omics = subset_blueberry_mae_views(
      mae_clr,
      c("species", "dna_pathways", "rna_pathways", "dna_ecs", "rna_ecs")
    )
  )

  out_bases <- c(
    full = file.path("data", "Blueberry_MAE_CLR"),
    pathways = file.path("data", "Blueberry_MAE_pathways_CLR"),
    ecs = file.path("data", "Blueberry_MAE_ecs_CLR"),
    all_omics = file.path("data", "Blueberry_MAE_all_omics_CLR")
  )

  for (nm in names(maes)) {
    save_mae_pair(maes[[nm]], out_bases[[nm]], object_name = "mae")
  }

  maes
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  input_mae_path <- if (length(args) >= 1L) args[[1L]] else file.path("data", "Blueberry_MAE.rds")
  maes <- build_blueberry_clr_maes(input_mae_path = input_mae_path)

  for (nm in names(maes)) {
    exp_names <- names(MultiAssayExperiment::experiments(maes[[nm]]))
    message(
      nm,
      ": ",
      paste(exp_names, collapse = ", "),
      " | samples in colData = ",
      nrow(as.data.frame(MultiAssayExperiment::colData(maes[[nm]])))
    )
  }
}
