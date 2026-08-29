# -----------------------------------------------------------------------------
# Negative-control permutation workflows
# -----------------------------------------------------------------------------

zentangler_validate_negative_control_type <- function(null_type) {
  match.arg(null_type, choices = c("permute_y"))
}

zentangler_validate_n_perm <- function(n_perm) {
  n_perm <- as.integer(n_perm)
  if (!is.finite(n_perm) || length(n_perm) != 1L || is.na(n_perm) || n_perm < 1L) {
    stop("n_perm must be a single integer >= 1.", call. = FALSE)
  }
  n_perm
}

zentangler_perm_seed <- function(seed, perm_id, offset = 100000L) {
  as.integer(seed) + as.integer(offset) + as.integer(perm_id)
}

zentangler_fit_seed <- function(seed, perm_id, offset = 0L) {
  as.integer(seed) + as.integer(offset) + as.integer(perm_id)
}

zentangler_permute_mae_coldata <- function(mae, var, seed = NULL) {
  zentangler_require_mae()
  if (!is.null(seed)) set.seed(seed)
  mae_null <- mae
  cd <- MultiAssayExperiment::colData(mae_null)
  if (!(var %in% colnames(cd))) {
    stop("Outcome column not found in MAE colData: ", var, call. = FALSE)
  }
  cd[[var]] <- sample(as.vector(cd[[var]]))
  MultiAssayExperiment::colData(mae_null) <- cd
  mae_null
}

zentangler_make_negative_control_mae <- function(mae, y_var, null_type = c("permute_y"), seed = NULL) {
  null_type <- zentangler_validate_negative_control_type(null_type)
  if (identical(null_type, "permute_y")) {
    return(zentangler_permute_mae_coldata(mae = mae, var = y_var, seed = seed))
  }
  stop("Unsupported null_type: ", null_type, call. = FALSE)
}

zentangler_get_primary_q_values <- function(tab) {
  if (!is.data.frame(tab) || nrow(tab) == 0L) return(numeric(0))
  if ("q_primary" %in% names(tab)) return(suppressWarnings(as.numeric(tab$q_primary)))
  if ("q_primary_global" %in% names(tab)) return(suppressWarnings(as.numeric(tab$q_primary_global)))
  numeric(0)
}

zentangler_min_finite_or_na <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x)]
  if (length(x) == 0L) return(NA_real_)
  min(x)
}

zentangler_max_finite_or_na <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x)]
  if (length(x) == 0L) return(NA_real_)
  max(x)
}

zentangler_mean_or_na <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x)]
  if (length(x) == 0L) return(NA_real_)
  mean(x)
}

zentangler_empirical_ge <- function(null_vals, obs_val) {
  null_vals <- suppressWarnings(as.numeric(null_vals))
  if (!is.finite(obs_val)) return(NA_real_)
  ok <- is.finite(null_vals)
  if (!any(ok)) return(NA_real_)
  mean(null_vals[ok] >= obs_val)
}

zentangler_empirical_le <- function(null_vals, obs_val) {
  null_vals <- suppressWarnings(as.numeric(null_vals))
  if (!is.finite(obs_val)) return(NA_real_)
  ok <- is.finite(null_vals)
  if (!any(ok)) return(NA_real_)
  mean(null_vals[ok] <= obs_val)
}

zentangler_bind_nonempty_rows <- function(rows) {
  keep <- vapply(rows, function(x) is.data.frame(x) && nrow(x) > 0L, logical(1))
  if (!any(keep)) return(data.frame())
  out <- do.call(rbind, rows[keep])
  rownames(out) <- NULL
  out
}

zentangler_make_parallel_mediator_stat_table <- function(
  fit,
  q_threshold = 0.25,
  perm_id = 0L,
  null_type = "observed",
  fit_seed = NA_integer_,
  permutation_seed = NA_integer_
) {
  tab <- zentangler_all_mediators(fit)
  if (!is.data.frame(tab) || nrow(tab) == 0L) return(data.frame())
  q <- zentangler_get_primary_q_values(tab)
  if (length(q) == 0L) q <- rep(NA_real_, nrow(tab))
  out <- data.frame(
    perm_id = as.integer(perm_id),
    null_type = as.character(null_type),
    fit_seed = as.integer(fit_seed),
    permutation_seed = as.integer(permutation_seed),
    omics = as.character(tab$omics %||% rep(NA_character_, nrow(tab))),
    mediator = as.character(tab$mediator %||% rep(NA_character_, nrow(tab))),
    mediator_id = if (all(c("omics", "mediator") %in% names(tab))) paste(tab$omics, tab$mediator, sep = "::") else as.character(seq_len(nrow(tab))),
    score = if ("score" %in% names(tab)) suppressWarnings(as.numeric(tab$score)) else NA_real_,
    abs_score = if ("abs_score" %in% names(tab)) suppressWarnings(as.numeric(tab$abs_score)) else if ("score" %in% names(tab)) abs(suppressWarnings(as.numeric(tab$score))) else NA_real_,
    p_primary = if ("p_primary" %in% names(tab)) suppressWarnings(as.numeric(tab$p_primary)) else NA_real_,
    q_primary = q,
    is_active = is.finite(q) & q <= q_threshold & is.finite(tab$b) & tab$b != 0,
    stringsAsFactors = FALSE
  )
  rownames(out) <- NULL
  out
}

zentangler_make_sequential_path_stat_table <- function(
  fit,
  q_threshold = 0.25,
  perm_id = 0L,
  null_type = "observed",
  fit_seed = NA_integer_,
  permutation_seed = NA_integer_
) {
  tab <- zentangler_sequential_paths(fit)
  if (!is.data.frame(tab) || nrow(tab) == 0L) return(data.frame())
  q <- zentangler_get_primary_q_values(tab)
  if (length(q) == 0L) q <- rep(NA_real_, nrow(tab))
  path_id <- if ("key_path" %in% names(tab)) {
    as.character(tab$key_path)
  } else if ("mediator_path" %in% names(tab)) {
    as.character(tab$mediator_path)
  } else {
    as.character(seq_len(nrow(tab)))
  }
  out <- data.frame(
    perm_id = as.integer(perm_id),
    null_type = as.character(null_type),
    fit_seed = as.integer(fit_seed),
    permutation_seed = as.integer(permutation_seed),
    path_id = path_id,
    sequential_score = if ("sequential_score" %in% names(tab)) suppressWarnings(as.numeric(tab$sequential_score)) else NA_real_,
    abs_sequential_score = if ("abs_sequential_score" %in% names(tab)) suppressWarnings(as.numeric(tab$abs_sequential_score)) else if ("sequential_score" %in% names(tab)) abs(suppressWarnings(as.numeric(tab$sequential_score))) else NA_real_,
    p_primary = if ("p_primary" %in% names(tab)) suppressWarnings(as.numeric(tab$p_primary)) else NA_real_,
    q_primary = q,
    is_active = is.finite(q) & q <= q_threshold & is.finite(tab$b) & tab$b != 0,
    stringsAsFactors = FALSE
  )
  rownames(out) <- NULL
  out
}

zentangler_parallel_empirical_mediator_summary <- function(observed_stats, null_stats) {
  if (!is.data.frame(observed_stats) || nrow(observed_stats) == 0L) return(data.frame())
  if (!is.data.frame(null_stats) || nrow(null_stats) == 0L) {
    out <- observed_stats
    out$n_null <- 0L
    out$n_null_active <- 0L
    out$null_mean_abs_score <- NA_real_
    out$null_max_abs_score <- NA_real_
    out$null_selection_frequency <- NA_real_
    out$empirical_p_is_active <- NA_real_
    out$empirical_p_abs_score <- NA_real_
    rownames(out) <- NULL
    return(out)
  }
  null_split <- split(null_stats, null_stats$mediator_id)
  rows <- lapply(seq_len(nrow(observed_stats)), function(ii) {
    obs <- observed_stats[ii, , drop = FALSE]
    nn <- null_split[[as.character(obs$mediator_id[[1L]])]]
    if (is.null(nn) || nrow(nn) == 0L) {
      obs$n_null <- 0L
      obs$n_null_active <- 0L
      obs$null_mean_abs_score <- NA_real_
      obs$null_max_abs_score <- NA_real_
      obs$null_selection_frequency <- NA_real_
      obs$empirical_p_is_active <- NA_real_
      obs$empirical_p_abs_score <- NA_real_
    } else {
      obs$n_null <- nrow(nn)
      obs$n_null_active <- sum(as.logical(nn$is_active), na.rm = TRUE)
      obs$null_mean_abs_score <- zentangler_mean_or_na(nn$abs_score)
      obs$null_max_abs_score <- zentangler_max_finite_or_na(nn$abs_score)
      obs$null_selection_frequency <- zentangler_mean_or_na(as.numeric(nn$is_active))
      obs$empirical_p_is_active <- zentangler_empirical_ge(
        null_vals = as.numeric(nn$is_active),
        obs_val = as.numeric(obs$is_active[[1L]])
      )
      obs$empirical_p_abs_score <- zentangler_empirical_ge(
        null_vals = nn$abs_score,
        obs_val = obs$abs_score[[1L]]
      )
    }
    obs
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

zentangler_sequential_empirical_path_summary <- function(observed_stats, null_stats) {
  if (!is.data.frame(observed_stats) || nrow(observed_stats) == 0L) return(data.frame())
  if (!is.data.frame(null_stats) || nrow(null_stats) == 0L) {
    out <- observed_stats
    out$n_null <- 0L
    out$n_null_active <- 0L
    out$null_mean_abs_sequential_score <- NA_real_
    out$null_max_abs_sequential_score <- NA_real_
    out$null_selection_frequency <- NA_real_
    out$empirical_p_is_active <- NA_real_
    out$empirical_p_abs_sequential_score <- NA_real_
    rownames(out) <- NULL
    return(out)
  }
  null_split <- split(null_stats, null_stats$path_id)
  rows <- lapply(seq_len(nrow(observed_stats)), function(ii) {
    obs <- observed_stats[ii, , drop = FALSE]
    nn <- null_split[[as.character(obs$path_id[[1L]])]]
    if (is.null(nn) || nrow(nn) == 0L) {
      obs$n_null <- 0L
      obs$n_null_active <- 0L
      obs$null_mean_abs_sequential_score <- NA_real_
      obs$null_max_abs_sequential_score <- NA_real_
      obs$null_selection_frequency <- NA_real_
      obs$empirical_p_is_active <- NA_real_
      obs$empirical_p_abs_sequential_score <- NA_real_
    } else {
      obs$n_null <- nrow(nn)
      obs$n_null_active <- sum(as.logical(nn$is_active), na.rm = TRUE)
      obs$null_mean_abs_sequential_score <- zentangler_mean_or_na(nn$abs_sequential_score)
      obs$null_max_abs_sequential_score <- zentangler_max_finite_or_na(nn$abs_sequential_score)
      obs$null_selection_frequency <- zentangler_mean_or_na(as.numeric(nn$is_active))
      obs$empirical_p_is_active <- zentangler_empirical_ge(
        null_vals = as.numeric(nn$is_active),
        obs_val = as.numeric(obs$is_active[[1L]])
      )
      obs$empirical_p_abs_sequential_score <- zentangler_empirical_ge(
        null_vals = nn$abs_sequential_score,
        obs_val = obs$abs_sequential_score[[1L]]
      )
    }
    obs
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

zentangler_summarize_parallel_negative_control <- function(
  fit,
  q_threshold = 0.25,
  perm_id = 0L,
  null_type = "observed",
  fit_seed = NA_integer_,
  permutation_seed = NA_integer_
) {
  tab <- zentangler_all_mediators(fit)
  q <- zentangler_get_primary_q_values(tab)
  score <- if (is.data.frame(tab) && "score" %in% names(tab)) tab$score else numeric(0)
  data.frame(
    perm_id = as.integer(perm_id),
    null_type = as.character(null_type),
    fit_seed = as.integer(fit_seed),
    permutation_seed = as.integer(permutation_seed),
    n_mediators = if (is.data.frame(tab)) nrow(tab) else 0L,
    n_active_mediators = sum(is.finite(q) & q <= q_threshold, na.rm = TRUE),
    best_q = zentangler_min_finite_or_na(q),
    max_abs_score = zentangler_max_finite_or_na(abs(score)),
    stringsAsFactors = FALSE
  )
}

zentangler_summarize_sequential_negative_control <- function(
  fit,
  q_threshold = 0.25,
  perm_id = 0L,
  null_type = "observed",
  fit_seed = NA_integer_,
  permutation_seed = NA_integer_
) {
  paths <- zentangler_sequential_paths(fit)
  q <- zentangler_get_primary_q_values(paths)
  score <- if (is.data.frame(paths) && "abs_sequential_score" %in% names(paths)) {
    paths$abs_sequential_score
  } else {
    numeric(0)
  }
  data.frame(
    perm_id = as.integer(perm_id),
    null_type = as.character(null_type),
    fit_seed = as.integer(fit_seed),
    permutation_seed = as.integer(permutation_seed),
    n_paths = if (is.data.frame(paths)) nrow(paths) else 0L,
    n_significant_paths = sum(is.finite(q) & q <= q_threshold, na.rm = TRUE),
    best_q = zentangler_min_finite_or_na(q),
    max_abs_sequential_score = zentangler_max_finite_or_na(score),
    stringsAsFactors = FALSE
  )
}

zentangler_parallel_empirical_summary <- function(observed_summary, null_summaries) {
  data.frame(
    null_type = "permute_y",
    n_null = nrow(null_summaries),
    obs_n_active_mediators = observed_summary$n_active_mediators[[1L]],
    obs_best_q = observed_summary$best_q[[1L]],
    obs_max_abs_score = observed_summary$max_abs_score[[1L]],
    null_mean_n_active_mediators = zentangler_mean_or_na(null_summaries$n_active_mediators),
    null_max_n_active_mediators = zentangler_max_finite_or_na(null_summaries$n_active_mediators),
    null_mean_max_abs_score = zentangler_mean_or_na(null_summaries$max_abs_score),
    empirical_p_n_active_mediators = zentangler_empirical_ge(
      null_vals = null_summaries$n_active_mediators,
      obs_val = observed_summary$n_active_mediators[[1L]]
    ),
    empirical_p_best_q = zentangler_empirical_le(
      null_vals = null_summaries$best_q,
      obs_val = observed_summary$best_q[[1L]]
    ),
    empirical_p_max_abs_score = zentangler_empirical_ge(
      null_vals = null_summaries$max_abs_score,
      obs_val = observed_summary$max_abs_score[[1L]]
    ),
    stringsAsFactors = FALSE
  )
}

zentangler_sequential_empirical_summary <- function(observed_summary, null_summaries) {
  data.frame(
    null_type = "permute_y",
    n_null = nrow(null_summaries),
    obs_n_significant_paths = observed_summary$n_significant_paths[[1L]],
    obs_best_q = observed_summary$best_q[[1L]],
    obs_max_abs_sequential_score = observed_summary$max_abs_sequential_score[[1L]],
    null_mean_n_significant_paths = zentangler_mean_or_na(null_summaries$n_significant_paths),
    null_max_n_significant_paths = zentangler_max_finite_or_na(null_summaries$n_significant_paths),
    null_mean_max_abs_sequential_score = zentangler_mean_or_na(null_summaries$max_abs_sequential_score),
    empirical_p_n_significant_paths = zentangler_empirical_ge(
      null_vals = null_summaries$n_significant_paths,
      obs_val = observed_summary$n_significant_paths[[1L]]
    ),
    empirical_p_best_q = zentangler_empirical_le(
      null_vals = null_summaries$best_q,
      obs_val = observed_summary$best_q[[1L]]
    ),
    empirical_p_max_abs_sequential_score = zentangler_empirical_ge(
      null_vals = null_summaries$max_abs_sequential_score,
      obs_val = observed_summary$max_abs_sequential_score[[1L]]
    ),
    stringsAsFactors = FALSE
  )
}

#' Run Y-permutation negative controls for parallel Zentangler
#'
#' Fits the observed parallel Zentangler model once, then permutes the outcome
#' column in `colData(mae)` `n_perm` times and refits the same model on each
#' permuted dataset. This produces an empirical negative-control reference for
#' mediator discoveries under a broken outcome relationship.
#'
#' Phase 1 supports `null_type = "permute_y"` for Gaussian and binomial
#' outcomes. Survival outcomes are not yet supported because they do not map to
#' a single `y_var` column.
#'
#' @param mae A `MultiAssayExperiment`.
#' @param x_var Exposure column in `colData(mae)`.
#' @param y_var Outcome column in `colData(mae)`.
#' @param n_perm Number of outcome permutations.
#' @param null_type Negative-control type. Currently only `"permute_y"`.
#' @param q_threshold Q-value cutoff used for active mediator counts.
#' @param seed Base seed used for the observed fit and to derive permutation
#'   seeds.
#' @param return_null_fits Logical. If `TRUE`, return the full fitted objects
#'   for all null runs; otherwise only return their summaries.
#' @param ... Additional arguments passed to `fit_multiview_parallel_zentangler()`.
#'
#' @return A list containing the observed fit, observed summary,
#'   observed mediator-level summary, null summaries, null mediator-level
#'   summaries, empirical summary, empirical mediator-level summary, and
#'   optionally the null fits.
#' @export
run_parallel_negative_controls <- function(
  mae,
  x_var = NULL,
  y_var = NULL,
  n_perm = 100L,
  null_type = c("permute_y"),
  q_threshold = 0.25,
  seed = 1,
  return_null_fits = FALSE,
  ...
) {
  null_type <- zentangler_validate_negative_control_type(null_type)
  n_perm <- zentangler_validate_n_perm(n_perm)
  q_threshold <- max(zentangler_clean_q_thresholds(q_threshold))

  dots <- list(...)
  y_family <- dots$y_family %||% "gaussian"
  if (identical(y_family, "survival")) {
    stop("run_parallel_negative_controls() currently supports only non-survival outcomes with a single y_var.", call. = FALSE)
  }

  observed_seed <- zentangler_fit_seed(seed = seed, perm_id = 0L)
  observed_fit <- fit_multiview_parallel_zentangler(
    mae = mae,
    x_var = x_var,
    y_var = y_var,
    seed = observed_seed,
    ...
  )
  observed_summary <- zentangler_summarize_parallel_negative_control(
    fit = observed_fit,
    q_threshold = q_threshold,
    perm_id = 0L,
    null_type = "observed",
    fit_seed = observed_seed,
    permutation_seed = NA_integer_
  )
  observed_mediator_stats <- zentangler_make_parallel_mediator_stat_table(
    fit = observed_fit,
    q_threshold = q_threshold,
    perm_id = 0L,
    null_type = "observed",
    fit_seed = observed_seed,
    permutation_seed = NA_integer_
  )

  null_rows <- vector("list", n_perm)
  null_mediator_rows <- vector("list", n_perm)
  null_fits <- if (isTRUE(return_null_fits)) vector("list", n_perm) else NULL
  for (perm_id in seq_len(n_perm)) {
    perm_seed <- zentangler_perm_seed(seed = seed, perm_id = perm_id)
    fit_seed <- zentangler_fit_seed(seed = seed, perm_id = perm_id)
    mae_null <- zentangler_make_negative_control_mae(
      mae = mae,
      y_var = y_var,
      null_type = null_type,
      seed = perm_seed
    )
    fit_null <- fit_multiview_parallel_zentangler(
      mae = mae_null,
      x_var = x_var,
      y_var = y_var,
      seed = fit_seed,
      ...
    )
    null_rows[[perm_id]] <- zentangler_summarize_parallel_negative_control(
      fit = fit_null,
      q_threshold = q_threshold,
      perm_id = perm_id,
      null_type = null_type,
      fit_seed = fit_seed,
      permutation_seed = perm_seed
    )
    null_mediator_rows[[perm_id]] <- zentangler_make_parallel_mediator_stat_table(
      fit = fit_null,
      q_threshold = q_threshold,
      perm_id = perm_id,
      null_type = null_type,
      fit_seed = fit_seed,
      permutation_seed = perm_seed
    )
    if (isTRUE(return_null_fits)) null_fits[[perm_id]] <- fit_null
  }

  null_summaries <- zentangler_bind_nonempty_rows(null_rows)
  null_mediator_summaries <- zentangler_bind_nonempty_rows(null_mediator_rows)
  empirical_summary <- zentangler_parallel_empirical_summary(
    observed_summary = observed_summary,
    null_summaries = null_summaries
  )
  empirical_mediator_summary <- zentangler_parallel_empirical_mediator_summary(
    observed_stats = observed_mediator_stats,
    null_stats = null_mediator_summaries
  )

  out <- list(
    observed_fit = observed_fit,
    observed_summary = observed_summary,
    observed_mediator_summary = observed_mediator_stats,
    null_summaries = null_summaries,
    null_mediator_summaries = null_mediator_summaries,
    empirical_summary = empirical_summary,
    empirical_mediator_summary = empirical_mediator_summary,
    settings = list(
      model_type = "parallel",
      null_type = null_type,
      n_perm = n_perm,
      q_threshold = q_threshold,
      seed = as.integer(seed),
      return_null_fits = isTRUE(return_null_fits)
    )
  )
  if (isTRUE(return_null_fits)) out$null_fits <- null_fits
  class(out) <- c("zentangler_parallel_negative_control", "list")
  out
}

#' Run Y-permutation negative controls for sequential Zentangler
#'
#' Fits the observed sequential Zentangler model once, then permutes the
#' outcome column in `colData(mae)` `n_perm` times and refits the same route on
#' each permuted dataset. This provides an empirical negative-control reference
#' for retained sequential paths under a broken outcome relationship.
#'
#' Phase 1 supports `null_type = "permute_y"` for Gaussian and binomial
#' outcomes. Survival outcomes are not yet supported because they do not map to
#' a single `y_var` column.
#'
#' @param mae A `MultiAssayExperiment`.
#' @param x_var Exposure column in `colData(mae)`.
#' @param y_var Outcome column in `colData(mae)`.
#' @param path_templates Route definition passed to `fit_sequential_zentangler()`.
#' @param n_perm Number of outcome permutations.
#' @param null_type Negative-control type. Currently only `"permute_y"`.
#' @param q_threshold Q-value cutoff used for significant path counts.
#' @param seed Base seed used for the observed fit and to derive permutation
#'   seeds.
#' @param return_null_fits Logical. If `TRUE`, return the full fitted objects
#'   for all null runs; otherwise only return their summaries.
#' @param ... Additional arguments passed to `fit_sequential_zentangler()`.
#'
#' @return A list containing the observed fit, observed summary,
#'   observed path-level summary, null summaries, null path-level summaries,
#'   empirical summary, empirical path-level summary, and optionally the null
#'   fits.
#' @export
run_sequential_negative_controls <- function(
  mae,
  x_var = NULL,
  y_var = NULL,
  path_templates = NULL,
  n_perm = 100L,
  null_type = c("permute_y"),
  q_threshold = 0.25,
  seed = 1,
  return_null_fits = FALSE,
  ...
) {
  null_type <- zentangler_validate_negative_control_type(null_type)
  n_perm <- zentangler_validate_n_perm(n_perm)
  q_threshold <- max(zentangler_clean_q_thresholds(q_threshold))

  dots <- list(...)
  y_family <- dots$y_family %||% "gaussian"
  if (identical(y_family, "survival")) {
    stop("run_sequential_negative_controls() currently supports only non-survival outcomes with a single y_var.", call. = FALSE)
  }

  observed_seed <- zentangler_fit_seed(seed = seed, perm_id = 0L)
  observed_fit <- fit_sequential_zentangler(
    mae = mae,
    x_var = x_var,
    y_var = y_var,
    path_templates = path_templates,
    seed = observed_seed,
    ...
  )
  observed_summary <- zentangler_summarize_sequential_negative_control(
    fit = observed_fit,
    q_threshold = q_threshold,
    perm_id = 0L,
    null_type = "observed",
    fit_seed = observed_seed,
    permutation_seed = NA_integer_
  )
  observed_path_stats <- zentangler_make_sequential_path_stat_table(
    fit = observed_fit,
    q_threshold = q_threshold,
    perm_id = 0L,
    null_type = "observed",
    fit_seed = observed_seed,
    permutation_seed = NA_integer_
  )

  null_rows <- vector("list", n_perm)
  null_path_rows <- vector("list", n_perm)
  null_fits <- if (isTRUE(return_null_fits)) vector("list", n_perm) else NULL
  for (perm_id in seq_len(n_perm)) {
    perm_seed <- zentangler_perm_seed(seed = seed, perm_id = perm_id)
    fit_seed <- zentangler_fit_seed(seed = seed, perm_id = perm_id)
    mae_null <- zentangler_make_negative_control_mae(
      mae = mae,
      y_var = y_var,
      null_type = null_type,
      seed = perm_seed
    )
    fit_null <- fit_sequential_zentangler(
      mae = mae_null,
      x_var = x_var,
      y_var = y_var,
      path_templates = path_templates,
      seed = fit_seed,
      ...
    )
    null_rows[[perm_id]] <- zentangler_summarize_sequential_negative_control(
      fit = fit_null,
      q_threshold = q_threshold,
      perm_id = perm_id,
      null_type = null_type,
      fit_seed = fit_seed,
      permutation_seed = perm_seed
    )
    null_path_rows[[perm_id]] <- zentangler_make_sequential_path_stat_table(
      fit = fit_null,
      q_threshold = q_threshold,
      perm_id = perm_id,
      null_type = null_type,
      fit_seed = fit_seed,
      permutation_seed = perm_seed
    )
    if (isTRUE(return_null_fits)) null_fits[[perm_id]] <- fit_null
  }

  null_summaries <- zentangler_bind_nonempty_rows(null_rows)
  null_path_summaries <- zentangler_bind_nonempty_rows(null_path_rows)
  empirical_summary <- zentangler_sequential_empirical_summary(
    observed_summary = observed_summary,
    null_summaries = null_summaries
  )
  empirical_path_summary <- zentangler_sequential_empirical_path_summary(
    observed_stats = observed_path_stats,
    null_stats = null_path_summaries
  )

  out <- list(
    observed_fit = observed_fit,
    observed_summary = observed_summary,
    observed_path_summary = observed_path_stats,
    null_summaries = null_summaries,
    null_path_summaries = null_path_summaries,
    empirical_summary = empirical_summary,
    empirical_path_summary = empirical_path_summary,
    settings = list(
      model_type = "sequential",
      null_type = null_type,
      n_perm = n_perm,
      q_threshold = q_threshold,
      seed = as.integer(seed),
      return_null_fits = isTRUE(return_null_fits)
    )
  )
  if (isTRUE(return_null_fits)) out$null_fits <- null_fits
  class(out) <- c("zentangler_sequential_negative_control", "list")
  out
}
