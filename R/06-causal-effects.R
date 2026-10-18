infer_effect_contrast <- function(X, x0 = NULL, x1 = NULL) {
  x_fin <- as.numeric(X[is.finite(X)])
  if (length(x_fin) == 0) stop("Cannot infer effect contrast from non-finite X.")

  if (is.null(x0) || !is.finite(as.numeric(x0))) {
    ux <- sort(unique(x_fin))
    x0 <- if (length(ux) == 2L) ux[1] else as.numeric(stats::quantile(x_fin, probs = 0.25, names = FALSE))
  }
  if (is.null(x1) || !is.finite(as.numeric(x1))) {
    ux <- sort(unique(x_fin))
    x1 <- if (length(ux) == 2L) ux[2] else as.numeric(stats::quantile(x_fin, probs = 0.75, names = FALSE))
  }
  x0 <- as.numeric(x0)
  x1 <- as.numeric(x1)
  if (!is.finite(x0) || !is.finite(x1)) stop("effect_x0 and effect_x1 must be finite.")
  if (isTRUE(all.equal(x0, x1))) stop("effect_x0 and effect_x1 must be different.")
  c(x0 = x0, x1 = x1, delta = x1 - x0)
}

estimate_reduced_x_summary <- function(Y, X, C = NULL, y_family = c("gaussian", "binomial", "survival"), x0 = NULL, x1 = NULL) {
  y_family <- match.arg(y_family)
  contrast <- infer_effect_contrast(X, x0 = x0, x1 = x1)
  dat <- data.frame(Y = I(Y), X = as.numeric(X), check.names = TRUE)
  if (!is.null(C) && ncol(C) > 0) {
    C_df <- as.data.frame(C, check.names = TRUE)
    colnames(C_df) <- make.names(colnames(C_df), unique = TRUE)
    dat <- cbind(dat, C_df)
  }

  fit <- try(
    if (identical(y_family, "survival")) {
      survival::coxph(Y ~ ., data = dat)
    } else if (identical(y_family, "gaussian")) {
      stats::lm(Y ~ ., data = dat)
    } else {
      suppressWarnings(stats::glm(Y ~ ., data = dat, family = stats::binomial()))
    },
    silent = TRUE
  )
  if (inherits(fit, "try-error") || is.null(stats::coef(fit))) {
    return(data.frame(
      reduced_total_effect = NA_real_,
      reduced_total_effect_coefficient = NA_real_,
      p_reduced_total_effect = NA_real_,
      stringsAsFactors = FALSE
    ))
  }

  co <- stats::coef(fit)
  x_coef <- suppressWarnings(as.numeric(co["X"]))
  if (!is.finite(x_coef)) x_coef <- NA_real_

  p_value <- NA_real_
  if (identical(y_family, "survival")) {
    sm <- try(summary(fit), silent = TRUE)
    if (!inherits(sm, "try-error") && !is.null(sm$coefficients) && "X" %in% rownames(sm$coefficients)) {
      p_value <- suppressWarnings(as.numeric(sm$coefficients["X", "Pr(>|z|)"]))
    }
    effect <- if (is.finite(x_coef)) x_coef * contrast[["delta"]] else NA_real_
  } else if (identical(y_family, "gaussian")) {
    sm <- try(summary(fit), silent = TRUE)
    if (!inherits(sm, "try-error") && !is.null(sm$coefficients) && "X" %in% rownames(sm$coefficients)) {
      p_value <- suppressWarnings(as.numeric(sm$coefficients["X", "Pr(>|t|)"]))
    }
    effect <- if (is.finite(x_coef)) x_coef * contrast[["delta"]] else NA_real_
  } else {
    sm <- try(summary(fit), silent = TRUE)
    if (!inherits(sm, "try-error") && !is.null(sm$coefficients) && "X" %in% rownames(sm$coefficients)) {
      p_value <- suppressWarnings(as.numeric(sm$coefficients["X", "Pr(>|z|)"]))
    }
    co[!is.finite(co)] <- 0
    W <- stats::model.matrix(fit)
    W0 <- W
    W1 <- W
    x_col <- which(colnames(W) == "X")
    if (length(x_col) == 1L) {
      W0[, x_col] <- contrast[["x0"]]
      W1[, x_col] <- contrast[["x1"]]
      effect <- mean(stats::plogis(as.numeric(W1 %*% co)) - stats::plogis(as.numeric(W0 %*% co)), na.rm = TRUE)
    } else {
      effect <- NA_real_
    }
  }

  data.frame(
    reduced_total_effect = effect,
    reduced_total_effect_coefficient = x_coef,
    p_reduced_total_effect = p_value,
    stringsAsFactors = FALSE
  )
}

# -----------------------------------------------------------------------------
# Post-selection causal effect stage (parallel)
#
# Effects are estimated for a fixed effect set of mediators (by default the
# mediators with q_primary <= effect_q_threshold). For that set the stage fits
#   mediator models   M_j ~ X + C                (linear, one per mediator)
#   outcome model     Y ~ X + C + M_set          (lm / logistic glm / Cox)
# without penalization, and computes
#   NDE = E[Y(x1, M(x0))] - E[Y(x0, M(x0))]
#   NIE = E[Y(x1, M(x1))] - E[Y(x1, M(x0))]
#   TE  = NDE + NIE
# by parametric g-computation. The outcome scale is the mean difference
# (gaussian), the risk difference (binomial) or the conditional log hazard
# ratio (survival, computed on the Cox linear predictor). Per-mediator indirect
# effects shift one mediator at a time and are on the same scale as NIE. The
# model assumes no exposure-mediator interaction.
# -----------------------------------------------------------------------------

zentangler_effect_set <- function(combined_mediators, q_threshold = 0.25, q_col = "q_primary") {
  if (is.null(combined_mediators) || nrow(combined_mediators) == 0L) {
    return(data.frame(omics = character(0), mediator = character(0), stringsAsFactors = FALSE))
  }
  q <- combined_mediators[[q_col]]
  keep <- is.finite(q) & q <= q_threshold & is.finite(combined_mediators$b) & combined_mediators$b != 0
  data.frame(
    omics = as.character(combined_mediators$omics[keep]),
    mediator = as.character(combined_mediators$mediator[keep]),
    stringsAsFactors = FALSE
  )
}

make_effect_mediator_matrix <- function(effect_set, blocks, n) {
  if (is.null(effect_set) || nrow(effect_set) == 0L) {
    return(list(M = matrix(nrow = n, ncol = 0), map = data.frame(omics = character(0), mediator = character(0), design_name = character(0), stringsAsFactors = FALSE)))
  }
  mats <- list()
  keep <- logical(nrow(effect_set))
  for (ii in seq_len(nrow(effect_set))) {
    view <- effect_set$omics[ii]
    med <- effect_set$mediator[ii]
    if (!(view %in% names(blocks)) || !(med %in% colnames(blocks[[view]]))) next
    keep[ii] <- TRUE
    mats[[length(mats) + 1L]] <- as.numeric(blocks[[view]][, med])
  }
  map <- effect_set[keep, c("omics", "mediator"), drop = FALSE]
  map$design_name <- if (nrow(map) > 0L) make.names(paste0("M_", map$omics, "__", map$mediator), unique = TRUE) else character(0)
  if (length(mats) == 0L) return(list(M = matrix(nrow = n, ncol = 0), map = map))
  M <- do.call(cbind, mats)
  colnames(M) <- map$design_name
  rownames(map) <- NULL
  list(M = M, map = map)
}

fit_active_outcome_effect_model <- function(Y, X, M_active, C = NULL, y_family = c("gaussian", "binomial", "survival")) {
  # Unpenalized outcome model used by the causal effect stage. Returns
  # coefficients for the design matrix (including a zero intercept for Cox,
  # whose linear predictor has no intercept), or NULL coefficients with a
  # warning when the fit fails.
  y_family <- match.arg(y_family)
  Z <- cbind(`(Intercept)` = 1, X = as.numeric(X))
  if (!is.null(C) && ncol(C) > 0) Z <- cbind(Z, as.matrix(C))
  if (!is.null(M_active) && ncol(M_active) > 0) Z <- cbind(Z, as.matrix(M_active))
  Z <- as.matrix(Z)

  if (ncol(Z) >= nrow(Z) - 2L) {
    warning("Effect model has as many parameters as samples; causal effects are NA. Use a stricter effect_q_threshold.", call. = FALSE)
    return(list(coefficients = NULL, design = Z))
  }
  fit <- try(
    if (identical(y_family, "survival")) {
      survival::coxph(Y ~ Z[, -1, drop = FALSE])
    } else if (identical(y_family, "gaussian")) {
      stats::lm.fit(x = Z, y = as.numeric(Y))
    } else {
      suppressWarnings(stats::glm.fit(x = Z, y = as.numeric(Y), family = stats::binomial(), intercept = TRUE))
    },
    silent = TRUE
  )
  if (inherits(fit, "try-error") || is.null(fit$coefficients)) {
    warning("Causal effect outcome model failed; causal effects are NA.", call. = FALSE)
    return(list(coefficients = NULL, design = Z))
  }
  co <- as.numeric(fit$coefficients)
  if (identical(y_family, "survival")) co <- c(0, co)
  names(co) <- colnames(Z)
  if (any(!is.finite(co))) {
    warning("Causal effect outcome model is rank deficient or separated; causal effects are NA.", call. = FALSE)
    return(list(coefficients = NULL, design = Z))
  }
  if (identical(y_family, "binomial") && isFALSE(fit$converged)) {
    warning("Causal effect logistic model did not converge (possible separation); causal effects are NA.", call. = FALSE)
    return(list(coefficients = NULL, design = Z))
  }
  list(coefficients = co, design = Z)
}

predict_effect_mean <- function(co, Z_template, X_value, M_value = NULL, y_family = c("gaussian", "binomial", "survival")) {
  y_family <- match.arg(y_family)
  Z <- Z_template
  Z[, "X"] <- as.numeric(X_value)
  if (!is.null(M_value) && ncol(M_value) > 0) {
    hit <- intersect(colnames(M_value), colnames(Z))
    Z[, hit] <- M_value[, hit, drop = FALSE]
  }
  eta <- as.numeric(Z %*% co[colnames(Z)])
  if (identical(y_family, "binomial")) stats::plogis(eta) else eta
}

zentangler_effect_labels <- function(y_family, sequential = FALSE) {
  suffix <- if (sequential) "_sequential" else ""
  list(
    method = paste0(switch(y_family,
      gaussian = "parametric_g_computation_linear",
      binomial = "parametric_g_computation_logistic",
      survival = "parametric_g_computation_cox_linear_predictor"
    ), suffix),
    scale = switch(y_family,
      gaussian = "mean_difference",
      binomial = "risk_difference",
      survival = "log_hazard_ratio"
    )
  )
}

zentangler_na_effect_row <- function(y_family, contrast, reduced_summary, n_active, view_names = character(0), sequential = FALSE) {
  lab <- zentangler_effect_labels(y_family, sequential)
  out <- data.frame(
    effect_method = lab$method,
    effect_scale = lab$scale,
    x0 = contrast[["x0"]], x1 = contrast[["x1"]], delta_x = contrast[["delta"]],
    nde = NA_real_, nie_total = NA_real_, te = NA_real_, prop_mediated = NA_real_,
    reduced_total_effect = reduced_summary$reduced_total_effect[[1L]],
    reduced_total_effect_coefficient = reduced_summary$reduced_total_effect_coefficient[[1L]],
    p_reduced_total_effect = reduced_summary$p_reduced_total_effect[[1L]],
    direct_coefficient_active_model = NA_real_,
    n_active_total = as.integer(n_active),
    stringsAsFactors = FALSE
  )
  for (view in view_names) {
    safe <- make.names(view)
    out[[paste0("nie_", safe)]] <- NA_real_
    out[[paste0("n_active_", safe)]] <- 0L
  }
  out
}

compute_effect_decomposition_multiview <- function(
  effect_set,
  Y,
  X,
  blocks,
  C = NULL,
  y_family = c("gaussian", "binomial", "survival"),
  x0,
  x1
) {
  y_family <- match.arg(y_family)
  contrast <- c(x0 = as.numeric(x0), x1 = as.numeric(x1), delta = as.numeric(x1) - as.numeric(x0))
  n <- length(as.numeric(X))
  reduced_summary <- estimate_reduced_x_summary(Y = Y, X = X, C = C, y_family = y_family, x0 = x0, x1 = x1)
  info <- make_effect_mediator_matrix(effect_set, blocks, n)
  M_obs <- info$M
  map <- info$map
  view_names <- names(blocks)
  lab <- zentangler_effect_labels(y_family)

  effect_fit <- fit_active_outcome_effect_model(Y = Y, X = X, M_active = M_obs, C = C, y_family = y_family)
  if (is.null(effect_fit$coefficients)) {
    return(list(
      effects = zentangler_na_effect_row(y_family, contrast, reduced_summary, nrow(map), view_names),
      mediator_effects = data.frame(omics = map$omics, mediator = map$mediator, a_effect_model = rep(NA_real_, nrow(map)), b_effect_model = rep(NA_real_, nrow(map)), indirect_effect = rep(NA_real_, nrow(map)), stringsAsFactors = FALSE)
    ))
  }
  co <- effect_fit$coefficients
  Z_template <- effect_fit$design

  # Linear mediator models M_j ~ X + C on the same samples.
  a_eff <- setNames(rep(NA_real_, nrow(map)), map$design_name)
  if (ncol(M_obs) > 0L) {
    a_tab <- zentangler_a_stage_lm(M_obs, X, C)
    a_eff[] <- a_tab$a_hat
  }
  M_x0 <- M_obs
  M_x1 <- M_obs
  for (dn in colnames(M_obs)) {
    M_x0[, dn] <- M_obs[, dn] + a_eff[[dn]] * (x0 - as.numeric(X))
    M_x1[, dn] <- M_obs[, dn] + a_eff[[dn]] * (x1 - as.numeric(X))
  }

  mu_00 <- predict_effect_mean(co, Z_template, X_value = x0, M_value = M_x0, y_family = y_family)
  mu_10 <- predict_effect_mean(co, Z_template, X_value = x1, M_value = M_x0, y_family = y_family)
  mu_11 <- predict_effect_mean(co, Z_template, X_value = x1, M_value = M_x1, y_family = y_family)

  nde <- mean(mu_10 - mu_00)
  nie <- mean(mu_11 - mu_10)
  te <- nde + nie

  out <- data.frame(
    effect_method = lab$method,
    effect_scale = lab$scale,
    x0 = x0,
    x1 = x1,
    delta_x = contrast[["delta"]],
    nde = nde,
    nie_total = nie,
    te = te,
    prop_mediated = if (is.finite(te) && abs(te) > 1e-08) nie / te else NA_real_,
    reduced_total_effect = reduced_summary$reduced_total_effect[[1L]],
    reduced_total_effect_coefficient = reduced_summary$reduced_total_effect_coefficient[[1L]],
    p_reduced_total_effect = reduced_summary$p_reduced_total_effect[[1L]],
    direct_coefficient_active_model = as.numeric(co[["X"]]),
    n_active_total = nrow(map),
    stringsAsFactors = FALSE
  )

  # Per-view and per-mediator indirect effects, always reported for every view
  # so bootstrap replicates have identical columns.
  for (view in view_names) {
    safe <- make.names(view)
    idx <- map$omics == view
    if (any(idx)) {
      M_view <- M_x0
      cols <- map$design_name[idx]
      M_view[, cols] <- M_x1[, cols, drop = FALSE]
      mu_view <- predict_effect_mean(co, Z_template, X_value = x1, M_value = M_view, y_family = y_family)
      out[[paste0("nie_", safe)]] <- mean(mu_view - mu_10)
    } else {
      out[[paste0("nie_", safe)]] <- 0
    }
    out[[paste0("n_active_", safe)]] <- sum(idx)
  }

  med_eff <- data.frame(
    omics = map$omics,
    mediator = map$mediator,
    a_effect_model = as.numeric(a_eff),
    b_effect_model = as.numeric(co[map$design_name]),
    indirect_effect = rep(NA_real_, nrow(map)),
    stringsAsFactors = FALSE
  )
  for (ii in seq_len(nrow(map))) {
    dn <- map$design_name[ii]
    M_one <- M_x0
    M_one[, dn] <- M_x1[, dn]
    mu_one <- predict_effect_mean(co, Z_template, X_value = x1, M_value = M_one, y_family = y_family)
    med_eff$indirect_effect[ii] <- mean(mu_one - mu_10)
  }

  list(effects = out, mediator_effects = med_eff)
}

# -----------------------------------------------------------------------------
# Post-selection causal effect stage (sequential)
#
# For the fixed set of significant paths, every mediator on those paths gets a
# linear mediator model conditional on X, C, its parents and all upstream
# ancestors on the effect paths; the outcome model conditions on X, C and every
# mediator in the set. NDE / NIE / TE then follow by g-computation, propagating
# the exposure through the mediator models layer by layer. Survival effects
# are on the Cox linear-predictor (log hazard ratio) scale.
# -----------------------------------------------------------------------------

zentangler_sequential_effect_paths <- function(paths, q_threshold = 0.25) {
  if (is.null(paths) || nrow(paths) == 0L) return(data.frame(key_path = character(0), stringsAsFactors = FALSE))
  keep <- is.finite(paths$q_primary) & paths$q_primary <= q_threshold & is.finite(paths$b) & paths$b != 0
  paths[keep, , drop = FALSE]
}

zentangler_sequential_key_structure <- function(paths) {
  keys_list <- strsplit(as.character(paths$key_path), " -> ", fixed = TRUE)
  stage <- list()
  parents <- list()
  ancestors <- list()
  for (keys in keys_list) {
    for (j in seq_along(keys)) {
      k <- keys[[j]]
      stage[[k]] <- j
      if (j > 1L) {
        parents[[k]] <- unique(c(parents[[k]], keys[[j - 1L]]))
        ancestors[[k]] <- unique(c(ancestors[[k]], keys[seq_len(j - 1L)]))
      }
    }
  }
  list(stage = unlist(stage), parents = parents, ancestors = ancestors)
}

zentangler_fit_sequential_mediator_models <- function(paths, data_by_key, X, C = NULL) {
  ks <- zentangler_sequential_key_structure(paths)
  keys <- names(sort(ks$stage))
  models <- list()
  for (key in keys) {
    regs <- ks$ancestors[[key]] %||% character(0)
    D <- cbind(`(Intercept)` = 1, X = as.numeric(X))
    if (!is.null(C) && ncol(C) > 0L) D <- cbind(D, as.matrix(C))
    if (length(regs) > 0L) {
      R <- do.call(cbind, lapply(regs, function(k) as.numeric(data_by_key[[k]])))
      colnames(R) <- regs
      D <- cbind(D, R)
    }
    fit <- stats::lm.fit(D, as.numeric(data_by_key[[key]]))
    co <- fit$coefficients
    names(co) <- colnames(D)
    if (any(!is.finite(co))) stop("Sequential mediator model for ", key, " is rank deficient.", call. = FALSE)
    models[[key]] <- list(stage = ks$stage[[key]], regressors = regs, coefficients = co)
  }
  list(models = models, keys = keys)
}

zentangler_predict_sequential_mediators <- function(model_info, X_value, n, C = NULL) {
  keys <- model_info$keys
  out <- matrix(NA_real_, nrow = n, ncol = length(keys), dimnames = list(NULL, keys))
  for (key in keys) {
    m <- model_info$models[[key]]
    co <- m$coefficients
    val <- rep(co[["(Intercept)"]], n) + co[["X"]] * as.numeric(X_value)
    if (!is.null(C) && ncol(C) > 0L) val <- val + as.numeric(as.matrix(C) %*% co[colnames(C)])
    for (r in m$regressors) val <- val + co[[r]] * out[, r]
    out[, key] <- val
  }
  out
}

compute_effect_decomposition_sequential <- function(
  sequential_paths,
  Y,
  X,
  data_by_key,
  C = NULL,
  y_family = c("gaussian", "binomial", "survival"),
  x0,
  x1,
  n_layers = 2L
) {
  y_family <- match.arg(y_family)
  contrast <- c(x0 = as.numeric(x0), x1 = as.numeric(x1), delta = as.numeric(x1) - as.numeric(x0))
  n <- length(as.numeric(X))
  lab <- zentangler_effect_labels(y_family, sequential = TRUE)
  reduced_summary <- estimate_reduced_x_summary(Y = Y, X = X, C = C, y_family = y_family, x0 = x0, x1 = x1)
  stage_names <- paste0("stage", seq_len(n_layers))
  na_row <- function(n_paths, n_keys) {
    out <- zentangler_na_effect_row(y_family, contrast, reduced_summary, n_keys, sequential = TRUE)
    out$n_active_paths <- n_paths
    out$n_active_mediators <- n_keys
    for (st in stage_names) out[[paste0("nie_", st)]] <- NA_real_
    out
  }

  has_paths <- !is.null(sequential_paths) && nrow(sequential_paths) > 0L
  if (has_paths) {
    model_info <- try(zentangler_fit_sequential_mediator_models(sequential_paths, data_by_key, X, C), silent = TRUE)
    if (inherits(model_info, "try-error")) {
      warning("Sequential mediator models failed; causal effects are NA.", call. = FALSE)
      return(na_row(nrow(sequential_paths), NA_integer_))
    }
    keys <- model_info$keys
    M_obs <- do.call(cbind, lapply(keys, function(k) as.numeric(data_by_key[[k]])))
    colnames(M_obs) <- keys
  } else {
    model_info <- list(models = list(), keys = character(0))
    keys <- character(0)
    M_obs <- matrix(nrow = n, ncol = 0)
  }
  design_names <- if (length(keys) > 0L) make.names(paste0("M_", keys), unique = TRUE) else character(0)
  M_named <- M_obs
  colnames(M_named) <- design_names
  effect_fit <- fit_active_outcome_effect_model(Y = Y, X = X, M_active = M_named, C = C, y_family = y_family)
  if (is.null(effect_fit$coefficients)) return(na_row(if (has_paths) nrow(sequential_paths) else 0L, length(keys)))
  co <- effect_fit$coefficients
  Z_template <- effect_fit$design

  rename <- function(M) {
    colnames(M) <- design_names
    M
  }
  M_x0 <- rename(zentangler_predict_sequential_mediators(model_info, x0, n, C))
  M_x1 <- rename(zentangler_predict_sequential_mediators(model_info, x1, n, C))
  # Natural-course mediators under x0/x1 are model predictions; keep the
  # observed residual so the g-computation is centred on the observed data.
  if (length(keys) > 0L) {
    M_fit_obs <- rename(zentangler_predict_sequential_mediators(model_info, as.numeric(X), n, C))
    resid <- M_named - M_fit_obs
    M_x0 <- M_x0 + resid
    M_x1 <- M_x1 + resid
  }
  mu_00 <- predict_effect_mean(co, Z_template, X_value = x0, M_value = M_x0, y_family = y_family)
  mu_10 <- predict_effect_mean(co, Z_template, X_value = x1, M_value = M_x0, y_family = y_family)
  mu_11 <- predict_effect_mean(co, Z_template, X_value = x1, M_value = M_x1, y_family = y_family)
  nde <- mean(mu_10 - mu_00)
  nie <- mean(mu_11 - mu_10)
  te <- nde + nie

  out <- data.frame(
    effect_method = lab$method,
    effect_scale = lab$scale,
    x0 = x0, x1 = x1, delta_x = contrast[["delta"]],
    nde = nde, nie_total = nie, te = te,
    prop_mediated = if (is.finite(te) && abs(te) > 1e-08) nie / te else NA_real_,
    reduced_total_effect = reduced_summary$reduced_total_effect[[1L]],
    reduced_total_effect_coefficient = reduced_summary$reduced_total_effect_coefficient[[1L]],
    p_reduced_total_effect = reduced_summary$p_reduced_total_effect[[1L]],
    direct_coefficient_active_model = as.numeric(co[["X"]]),
    n_active_total = length(keys),
    n_active_paths = if (has_paths) nrow(sequential_paths) else 0L,
    n_active_mediators = length(keys),
    indirect_product_sum = if (has_paths) sum(sequential_paths$sequential_score, na.rm = TRUE) * contrast[["delta"]] else 0,
    stringsAsFactors = FALSE
  )
  stage_of <- if (length(keys) > 0L) model_info$models[keys] else list()
  for (st in seq_len(n_layers)) {
    nm <- paste0("nie_stage", st)
    cols <- design_names[vapply(stage_of, function(m) identical(as.integer(m$stage), st), logical(1))]
    if (length(cols) == 0L) {
      out[[nm]] <- 0
    } else {
      # Shift only stage-st mediators to their x1 values (holding the rest at x0).
      M_st <- M_x0
      M_st[, cols] <- M_x1[, cols, drop = FALSE]
      out[[nm]] <- mean(predict_effect_mean(co, Z_template, X_value = x1, M_value = M_st, y_family = y_family) - mu_10)
    }
  }
  out
}
