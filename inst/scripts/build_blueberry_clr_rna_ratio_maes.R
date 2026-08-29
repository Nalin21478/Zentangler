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

compute_ratio_pseudocount <- function(rna_mat, dna_mat) {
  pooled <- c(as.matrix(rna_mat), as.matrix(dna_mat))
  positive_vals <- pooled[is.finite(pooled) & pooled > 0]
  if (length(positive_vals) == 0L) return(1e-06)
  min(positive_vals) / 2
}

compute_log_rna_dna_ratio_se <- function(rna_se, dna_se, pseudocount = NULL) {
  assay_name_rna <- SummarizedExperiment::assayNames(rna_se)[[1L]]
  assay_name_dna <- SummarizedExperiment::assayNames(dna_se)[[1L]]
  rna_mat <- SummarizedExperiment::assay(rna_se, assay_name_rna)
  dna_mat <- SummarizedExperiment::assay(dna_se, assay_name_dna)

  common_features <- intersect(rownames(dna_mat), rownames(rna_mat))
  common_samples <- intersect(colnames(dna_mat), colnames(rna_mat))

  if (length(common_features) == 0L) {
    stop("No shared features found between RNA and DNA assays.", call. = FALSE)
  }
  if (length(common_samples) == 0L) {
    stop("No shared samples found between RNA and DNA assays.", call. = FALSE)
  }

  dna_mat <- dna_mat[common_features, common_samples, drop = FALSE]
  rna_mat <- rna_mat[common_features, common_samples, drop = FALSE]

  if (any(!is.finite(dna_mat), na.rm = TRUE) || any(!is.finite(rna_mat), na.rm = TRUE)) {
    stop("RNA/DNA ratio requires finite values in both assays.", call. = FALSE)
  }
  if (any(dna_mat < 0, na.rm = TRUE) || any(rna_mat < 0, na.rm = TRUE)) {
    stop("RNA/DNA ratio requires non-negative values in both assays.", call. = FALSE)
  }

  if (is.null(pseudocount)) {
    pseudocount <- compute_ratio_pseudocount(rna_mat, dna_mat)
  }
  if (!is.finite(pseudocount) || pseudocount <= 0) {
    stop("RNA/DNA ratio pseudocount must be a positive finite number.", call. = FALSE)
  }

  ratio_mat <- log((rna_mat + pseudocount) / (dna_mat + pseudocount))

  ratio_se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(abundance = ratio_mat),
    rowData = SummarizedExperiment::rowData(rna_se)[common_features, , drop = FALSE]
  )
  S4Vectors::metadata(ratio_se)$transform <- "log_RNA_over_DNA"
  S4Vectors::metadata(ratio_se)$ratio_pseudocount <- pseudocount
  ratio_se
}

build_blueberry_clr_rna_ratio_maes <- function(
  input_mae_path = file.path("data", "Blueberry_MAE.rds")
) {
  if (!file.exists(input_mae_path)) {
    stop("Input MAE not found: ", input_mae_path, call. = FALSE)
  }

  mae_raw <- readRDS(input_mae_path)

  experiments_transformed <- list(
    species = clr_transform_se(mae_raw[["species"]]),
    dna_pathways = clr_transform_se(mae_raw[["dna_pathways"]]),
    rna_pathways = compute_log_rna_dna_ratio_se(mae_raw[["rna_pathways"]], mae_raw[["dna_pathways"]]),
    dna_ecs = clr_transform_se(mae_raw[["dna_ecs"]]),
    rna_ecs = compute_log_rna_dna_ratio_se(mae_raw[["rna_ecs"]], mae_raw[["dna_ecs"]])
  )

  mae_transformed <- MultiAssayExperiment::MultiAssayExperiment(
    experiments = MultiAssayExperiment::ExperimentList(experiments_transformed),
    colData = S4Vectors::DataFrame(MultiAssayExperiment::colData(mae_raw)),
    sampleMap = S4Vectors::DataFrame(MultiAssayExperiment::sampleMap(mae_raw))
  )
  S4Vectors::metadata(mae_transformed)$transform <- c(
    species = "CLR",
    dna_pathways = "CLR",
    rna_pathways = "log_RNA_over_DNA",
    dna_ecs = "CLR",
    rna_ecs = "log_RNA_over_DNA"
  )

  maes <- list(
    full = mae_transformed,
    pathways = subset_blueberry_mae_views(
      mae_transformed,
      c("species", "dna_pathways", "rna_pathways")
    ),
    ecs = subset_blueberry_mae_views(
      mae_transformed,
      c("species", "dna_ecs", "rna_ecs")
    ),
    all_omics = subset_blueberry_mae_views(
      mae_transformed,
      c("species", "dna_pathways", "rna_pathways", "dna_ecs", "rna_ecs")
    )
  )

  out_bases <- c(
    full = file.path("data", "Blueberry_MAE_CLR_RNARatio"),
    pathways = file.path("data", "Blueberry_MAE_pathways_CLR_RNARatio"),
    ecs = file.path("data", "Blueberry_MAE_ecs_CLR_RNARatio"),
    all_omics = file.path("data", "Blueberry_MAE_all_omics_CLR_RNARatio")
  )

  for (nm in names(maes)) {
    save_mae_pair(maes[[nm]], out_bases[[nm]], object_name = "mae")
  }

  maes
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  input_mae_path <- if (length(args) >= 1L) args[[1L]] else file.path("data", "Blueberry_MAE.rds")
  maes <- build_blueberry_clr_rna_ratio_maes(input_mae_path = input_mae_path)

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
