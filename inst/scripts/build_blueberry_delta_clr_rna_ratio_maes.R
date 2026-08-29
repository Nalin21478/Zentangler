#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(MultiAssayExperiment)
  library(S4Vectors)
  library(SummarizedExperiment)
  library(readxl)
})

source(file.path("inst", "scripts", "build_blueberry_mae_from_excel.R"))

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

compute_ratio_pseudocount <- function(rna_mat, dna_mat) {
  pooled <- c(as.matrix(rna_mat), as.matrix(dna_mat))
  positive_vals <- pooled[is.finite(pooled) & pooled > 0]
  if (length(positive_vals) == 0L) return(1e-06)
  min(positive_vals) / 2
}

build_all_timepoint_experiments <- function(xlsx_path, metadata_samples) {
  list(
    species = build_species_experiment(xlsx_path, metadata_samples),
    dna_pathways = build_feature_by_sample_experiment(
      xlsx_path = xlsx_path,
      sheet_name = "Table 5 - pathways DNA",
      metadata_samples = metadata_samples,
      keep_unstratified_only = TRUE
    ),
    rna_pathways = build_feature_by_sample_experiment(
      xlsx_path = xlsx_path,
      sheet_name = "Table 6 - pathways RNA",
      metadata_samples = metadata_samples,
      keep_unstratified_only = TRUE
    ),
    dna_ecs = build_feature_by_sample_experiment(
      xlsx_path = xlsx_path,
      sheet_name = "Table 7 - ECs DNA",
      metadata_samples = metadata_samples,
      keep_unstratified_only = TRUE
    ),
    rna_ecs = build_feature_by_sample_experiment(
      xlsx_path = xlsx_path,
      sheet_name = "Table 8 - ECs RNA",
      metadata_samples = metadata_samples,
      keep_unstratified_only = TRUE
    )
  )
}

compute_delta_clr_experiment <- function(se, baseline_samples, month6_samples) {
  assay_name <- SummarizedExperiment::assayNames(se)[[1L]]
  mat <- SummarizedExperiment::assay(se, assay_name)
  clr_mat <- clr_transform_matrix(mat)

  delta_mat <- clr_mat[, month6_samples, drop = FALSE] - clr_mat[, baseline_samples, drop = FALSE]
  colnames(delta_mat) <- month6_samples

  delta_se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(abundance = delta_mat),
    rowData = SummarizedExperiment::rowData(se)
  )
  S4Vectors::metadata(delta_se)$transform <- "CLR_then_delta"
  S4Vectors::metadata(delta_se)$clr_pseudocount <- attr(clr_mat, "pseudocount")
  delta_se
}

compute_delta_log_ratio_experiment <- function(rna_se, dna_se, baseline_samples, month6_samples) {
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

  pseudocount <- compute_ratio_pseudocount(rna_mat, dna_mat)
  ratio_mat <- log((rna_mat + pseudocount) / (dna_mat + pseudocount))
  delta_mat <- ratio_mat[, month6_samples, drop = FALSE] - ratio_mat[, baseline_samples, drop = FALSE]
  colnames(delta_mat) <- month6_samples

  delta_se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(abundance = delta_mat),
    rowData = SummarizedExperiment::rowData(rna_se)[common_features, , drop = FALSE]
  )
  S4Vectors::metadata(delta_se)$transform <- "log_RNA_over_DNA_then_delta"
  S4Vectors::metadata(delta_se)$ratio_pseudocount <- pseudocount
  delta_se
}

build_blueberry_delta_clr_rna_ratio_maes <- function(
  xlsx_path = file.path("data", "Blueberry_Multiomics.xlsx"),
  metadata_path = file.path("data", "metadata_emma_longitudinal_only.csv")
) {
  if (!file.exists(xlsx_path)) {
    stop("Workbook not found: ", xlsx_path, call. = FALSE)
  }
  if (!file.exists(metadata_path)) {
    stop("Metadata not found: ", metadata_path, call. = FALSE)
  }

  metadata_df <- utils::read.csv(
    metadata_path,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  metadata_all <- normalize_blueberry_metadata(metadata_df, sample_col = "Sample")
  metadata_all$timepoint_num <- suppressWarnings(as.numeric(as.character(metadata_all$timepoint)))
  month6_metadata <- build_month6_outcome_metadata(metadata_df)

  month6_rows <- metadata_all[metadata_all$timepoint_num == 6, , drop = FALSE]
  baseline_rows <- metadata_all[metadata_all$timepoint_num == 0, , drop = FALSE]
  baseline_match <- match(month6_rows$subject, baseline_rows$subject)

  if (anyNA(baseline_match)) {
    stop("Some month-6 subjects do not have baseline samples.", call. = FALSE)
  }

  baseline_samples <- baseline_rows$sample[baseline_match]
  month6_samples <- month6_rows$sample

  experiments_all <- build_all_timepoint_experiments(
    xlsx_path = xlsx_path,
    metadata_samples = rownames(metadata_all)
  )

  keep_subject <- Reduce(
    `&`,
    lapply(experiments_all, function(se) {
      baseline_samples %in% colnames(se) & month6_samples %in% colnames(se)
    })
  )

  if (!all(keep_subject)) {
    baseline_samples <- baseline_samples[keep_subject]
    month6_samples <- month6_samples[keep_subject]
    month6_rows <- month6_rows[keep_subject, , drop = FALSE]
    month6_metadata <- month6_metadata[month6_samples, , drop = FALSE]
  } else {
    month6_metadata <- month6_metadata[month6_samples, , drop = FALSE]
  }

  experiments_delta <- list(
    species = compute_delta_clr_experiment(experiments_all[["species"]], baseline_samples, month6_samples),
    dna_pathways = compute_delta_clr_experiment(experiments_all[["dna_pathways"]], baseline_samples, month6_samples),
    rna_pathways = compute_delta_log_ratio_experiment(
      experiments_all[["rna_pathways"]],
      experiments_all[["dna_pathways"]],
      baseline_samples,
      month6_samples
    ),
    dna_ecs = compute_delta_clr_experiment(experiments_all[["dna_ecs"]], baseline_samples, month6_samples),
    rna_ecs = compute_delta_log_ratio_experiment(
      experiments_all[["rna_ecs"]],
      experiments_all[["dna_ecs"]],
      baseline_samples,
      month6_samples
    )
  )

  sample_map_rows <- lapply(names(experiments_delta), function(view) {
    sample_ids <- colnames(experiments_delta[[view]])
    data.frame(
      assay = view,
      primary = sample_ids,
      colname = sample_ids,
      stringsAsFactors = FALSE
    )
  })
  sample_map <- do.call(rbind, sample_map_rows)
  sample_map <- S4Vectors::DataFrame(
    assay = factor(sample_map$assay, levels = names(experiments_delta)),
    primary = sample_map$primary,
    colname = sample_map$colname
  )

  mae_delta <- MultiAssayExperiment::MultiAssayExperiment(
    experiments = MultiAssayExperiment::ExperimentList(experiments_delta),
    colData = S4Vectors::DataFrame(
      month6_metadata,
      row.names = rownames(month6_metadata)
    ),
    sampleMap = sample_map
  )
  S4Vectors::metadata(mae_delta)$transform <- c(
    species = "CLR_then_delta",
    dna_pathways = "CLR_then_delta",
    rna_pathways = "log_RNA_over_DNA_then_delta",
    dna_ecs = "CLR_then_delta",
    rna_ecs = "log_RNA_over_DNA_then_delta"
  )
  S4Vectors::metadata(mae_delta)$delta_definition <- "month6_minus_baseline"

  maes <- list(
    full = mae_delta,
    pathways = subset_blueberry_mae_views(
      mae_delta,
      c("species", "dna_pathways", "rna_pathways")
    ),
    ecs = subset_blueberry_mae_views(
      mae_delta,
      c("species", "dna_ecs", "rna_ecs")
    ),
    all_omics = subset_blueberry_mae_views(
      mae_delta,
      c("species", "dna_pathways", "rna_pathways", "dna_ecs", "rna_ecs")
    )
  )

  out_bases <- c(
    full = file.path("data", "Blueberry_MAE_delta_CLR_RNARatio"),
    pathways = file.path("data", "Blueberry_MAE_pathways_delta_CLR_RNARatio"),
    ecs = file.path("data", "Blueberry_MAE_ecs_delta_CLR_RNARatio"),
    all_omics = file.path("data", "Blueberry_MAE_all_omics_delta_CLR_RNARatio")
  )

  for (nm in names(maes)) {
    save_mae_pair(maes[[nm]], out_bases[[nm]], object_name = "mae")
  }

  maes
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  xlsx_path <- if (length(args) >= 1L) args[[1L]] else file.path("data", "Blueberry_Multiomics.xlsx")
  metadata_path <- if (length(args) >= 2L) args[[2L]] else file.path("data", "metadata_emma_longitudinal_only.csv")
  maes <- build_blueberry_delta_clr_rna_ratio_maes(
    xlsx_path = xlsx_path,
    metadata_path = metadata_path
  )

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
