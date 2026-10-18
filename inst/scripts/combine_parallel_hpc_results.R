#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

parse_args <- function(args) {
  out <- list()
  i <- 1L
  while (i <= length(args)) {
    key <- args[[i]]
    if (!grepl("^--", key)) {
      i <- i + 1L
      next
    }
    nm <- sub("^--", "", key)
    val <- if (i + 1L <= length(args) && !grepl("^--", args[[i + 1L]])) args[[i + 1L]] else TRUE
    out[[nm]] <- val
    i <- i + if (isTRUE(val)) 1L else 2L
  }
  out
}

split_arg <- function(x, default = character(0)) {
  if (is.null(x) || length(x) == 0L || all(is.na(x)) || !nzchar(paste(x, collapse = ""))) return(default)
  out <- trimws(unlist(strsplit(as.character(x), ",", fixed = TRUE), use.names = FALSE))
  out[nzchar(out)]
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
if (is.null(args$`results-dir`)) stop("Missing --results-dir", call. = FALSE)

results_dir <- args$`results-dir`
summary_out <- if (!is.null(args$`summary-out`)) args$`summary-out` else file.path(results_dir, "all_observed_summaries.csv")
best_out <- if (!is.null(args$`best-out`)) args$`best-out` else file.path(results_dir, "best_per_outcome.csv")
fusion_pref <- split_arg(args$`fusion-pref`, default = c("multiview", "late", "early"))
sis_pref <- suppressWarnings(as.numeric(split_arg(args$`sis-pref`, default = c("10", "25", "50", "100"))))
alpha_pref <- suppressWarnings(as.numeric(split_arg(args$`alpha-pref`, default = c("0.5", "1", "0"))))

files <- list.files(results_dir, pattern = "_observed_summary\\.csv$", full.names = TRUE)
if (length(files) == 0L) stop("No observed summary CSVs found in ", results_dir, call. = FALSE)

tabs <- lapply(files, function(f) utils::read.csv(f, stringsAsFactors = FALSE, check.names = FALSE))
amc <- do.call(rbind, tabs)
rownames(amc) <- NULL
if (!("exposure_mode" %in% colnames(amc))) amc$exposure_mode <- "existing"
if (!("exposure_label" %in% colnames(amc))) amc$exposure_label <- amc$exposure_mode

# fusion_mode column values: multiview (legacy "intermediate"), late, early.
amc$fusion_rank <- match(ifelse(amc$fusion_mode %in% "intermediate", "multiview", amc$fusion_mode), fusion_pref)
amc$sis_rank <- match(amc$sis_n, sis_pref)
amc$alpha_rank <- match(amc$glmnet_alpha, alpha_pref)
amc$fusion_rank[is.na(amc$fusion_rank)] <- length(fusion_pref) + 1L
amc$sis_rank[is.na(amc$sis_rank)] <- length(sis_pref) + 1L
amc$alpha_rank[is.na(amc$alpha_rank)] <- length(alpha_pref) + 1L

ord <- do.call(order, list(
  amc$exposure_mode,
  amc$outcome,
  -amc$n_active_mediators,
  amc$fusion_rank,
  amc$sis_rank,
  amc$alpha_rank,
  amc$mae
))
amc_sorted <- amc[ord, , drop = FALSE]
best_idx <- !duplicated(paste(amc_sorted$exposure_mode, amc_sorted$outcome, sep = "||"))
best <- amc_sorted[best_idx, , drop = FALSE]
best <- best[, setdiff(colnames(best), c("fusion_rank", "sis_rank", "alpha_rank")), drop = FALSE]
rownames(best) <- NULL

utils::write.csv(amc[, setdiff(colnames(amc), c("fusion_rank", "sis_rank", "alpha_rank")), drop = FALSE], summary_out, row.names = FALSE)
utils::write.csv(best, best_out, row.names = FALSE)
message("Wrote observed summary table: ", summary_out)
message("Wrote best-per-outcome table: ", best_out)
