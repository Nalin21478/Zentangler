# -----------------------------------------------------------------------------
# Fusion / Y-stage helpers
# -----------------------------------------------------------------------------

validate_glmnet_alpha <- function(glmnet_alpha) {
  if (length(glmnet_alpha) != 1L || !is.numeric(glmnet_alpha) || !is.finite(glmnet_alpha)) {
    stop("glmnet_alpha must be one finite numeric value between 0 and 1.")
  }
  if (glmnet_alpha < 0 || glmnet_alpha > 1) {
    stop("glmnet_alpha must be between 0 and 1. Use 1 for lasso, 0.5 for elastic net, and 0 for ridge.")
  }
  as.numeric(glmnet_alpha)
}

glmnet_alpha_label <- function(glmnet_alpha) {
  if (isTRUE(all.equal(glmnet_alpha, 1))) return("lasso")
  if (isTRUE(all.equal(glmnet_alpha, 0))) return("ridge")
  paste0("elastic net alpha=", format(glmnet_alpha, trim = TRUE, scientific = FALSE))
}

zentangler_glmnet_family <- function(y_family) {
  if (identical(y_family, "survival")) "cox" else y_family
}

zentangler_y_design <- function(X, C = NULL) {
  W <- cbind(X = as.numeric(X))
  if (!is.null(C) && ncol(C) > 0) W <- cbind(W, as.matrix(C))
  W
}

zentangler_make_foldid <- function(n, nfolds = 10L) {
  # Inside a bootstrap replicate, duplicated rows of the same original sample
  # (or cluster) must stay in the same fold; 07-bootstrap.R sets the grouping
  # through options(zentangler.fold_groups).
  groups <- getOption("zentangler.fold_groups")
  if (!is.null(groups) && length(groups) == n) {
    ug <- unique(groups)
    nfolds <- max(3L, min(as.integer(nfolds), length(ug)))
    g_fold <- sample(rep(seq_len(nfolds), length.out = length(ug)))
    return(g_fold[match(groups, ug)])
  }
  nfolds <- max(3L, min(as.integer(nfolds), n))
  sample(rep(seq_len(nfolds), length.out = n))
}

zentangler_subset_y <- function(Y, idx) {
  if (inherits(Y, "Surv")) Y[idx, , drop = FALSE] else Y[idx]
}

# -----------------------------------------------------------------------------
# Supplementary: IntegratedLearner-style early fusion (one glmnet model)
# -----------------------------------------------------------------------------

fit_y_multiview_glmnet_lasso <- function(
  Y,
  X,
  blocks,
  C = NULL,
  y_family = c("gaussian", "binomial", "survival"),
  lambda_choice = c("lambda.1se", "lambda.min"),
  glmnet_alpha = 1
) {
  y_family <- match.arg(y_family)
  lambda_choice <- match.arg(lambda_choice)
  glmnet_alpha <- validate_glmnet_alpha(glmnet_alpha)

  pref <- prefix_multiview_blocks(blocks)
  W <- zentangler_y_design(X, C)
  Z <- cbind(W, pref$X)

  pf <- c(rep(0, ncol(W)), rep(1, ncol(pref$X)))

  cvfit <- glmnet::cv.glmnet(
    x = as.matrix(Z),
    y = Y,
    family = zentangler_glmnet_family(y_family),
    alpha = glmnet_alpha,
    penalty.factor = pf,
    standardize = TRUE,
    foldid = zentangler_make_foldid(nrow(Z))
  )

  beta <- as.matrix(coef(cvfit, s = lambda_choice))[, 1]
  names(beta) <- gsub("`", "", names(beta), fixed = TRUE)

  coef_by_view <- lapply(names(blocks), function(view) {
    dn <- paste0(view, "::", colnames(blocks[[view]]))
    out <- setNames(rep(0, length(dn)), colnames(blocks[[view]]))
    hit <- intersect(dn, names(beta))
    if (length(hit) > 0) out[sub(paste0("^", view, "::"), "", hit)] <- beta[hit]
    out
  })
  names(coef_by_view) <- names(blocks)

  x_coef <- beta["X"]
  if (!is.finite(x_coef)) x_coef <- 0

  list(
    method = paste0("IntegratedLearner-style early fusion glmnet (", glmnet_alpha_label(glmnet_alpha), ")"),
    coef_by_view = coef_by_view,
    x_coef = as.numeric(x_coef),
    w_coef = setNames(as.numeric(beta[colnames(W)]), colnames(W)),
    fit = cvfit,
    feature_map = pref$map
  )
}

# -----------------------------------------------------------------------------
# Main engine: cooperative learning via the `multiview` package
# (Ding, Li, Narasimhan & Tibshirani, PNAS 2022)
#
#   min  loss(Y, W g + sum_k Mt_k b_k) + (rho/2) sum_{k<l} ||Mt_k b_k - Mt_l b_l||^2
#        + lambda * ||b||_1
#
# W = [X, C] is unpenalized. Mt_k are the view-k mediators residualized on
# [1, W]; this is a reparametrization (Frisch-Waugh), so b is still the
# coefficient of M_k in Y ~ X + C + M, while the agreement penalty acts on the
# covariate-adjusted view contributions. W enters through an offset that is
# refreshed by backfitting after tuning.
#
# rho and lambda are tuned jointly by K-fold cross-validation over
# `rho_grid` (shared folds; deviance for binomial/Cox, MSE for Gaussian).
# rho runs from 0 to 1: rho = 0 is early fusion (one lasso over all views),
# rho = 1 is late-fusion-like (views fitted to agree), and values in between
# trade the two off. The data choose rho, so no fusion type is fixed a priori.
# -----------------------------------------------------------------------------

zentangler_default_rho_grid <- function() c(0, 0.2, 0.4, 0.6, 0.8, 1)

zentangler_multiview_family <- function(y_family) {
  switch(y_family, gaussian = stats::gaussian(), binomial = stats::binomial(), survival = "cox")
}

zentangler_unpenalized_fit <- function(Y, W, y_family, offset = NULL) {
  # Unpenalized Y ~ W (+ offset); returns coefficients for W and the intercept.
  W <- as.matrix(W)
  off <- if (is.null(offset)) rep(0, nrow(W)) else offset
  if (identical(y_family, "survival")) {
    fit <- survival::coxph(Y ~ W + offset(off), ties = "breslow")
    co <- stats::coef(fit)
    return(list(intercept = 0, coef = setNames(as.numeric(co), colnames(W))))
  }
  fam <- if (identical(y_family, "binomial")) stats::binomial() else stats::gaussian()
  fit <- suppressWarnings(stats::glm.fit(cbind(1, W), as.numeric(Y), family = fam, offset = off))
  co <- fit$coefficients
  co[!is.finite(co)] <- 0
  list(intercept = co[[1L]], coef = setNames(as.numeric(co[-1L]), colnames(W)))
}

zentangler_cv_loss <- function(Y, eta, y_family) {
  if (identical(y_family, "gaussian")) return(mean((as.numeric(Y) - eta)^2))
  if (identical(y_family, "binomial")) {
    mu <- pmin(pmax(stats::plogis(eta), 1e-10), 1 - 1e-10)
    y <- as.numeric(Y)
    return(-2 * mean(y * log(mu) + (1 - y) * log(1 - mu)))
  }
  glmnet::coxnet.deviance(pred = eta, y = Y) / length(eta)
}

zentangler_multiview_path <- function(x_list, Y, rho, y_family, offset, lambda = NULL, alpha = 1) {
  args <- list(
    x_list = x_list, y = Y, rho = rho, family = zentangler_multiview_family(y_family),
    offset = offset, alpha = alpha, standardize = TRUE
  )
  if (!is.null(lambda)) args$lambda <- lambda
  if (identical(y_family, "survival")) args$intercept <- NULL
  nonconv <- FALSE
  fit <- withCallingHandlers(
    do.call(multiview::multiview, args),
    warning = function(w) {
      msg <- conditionMessage(w)
      if (grepl("did not converge", msg, fixed = TRUE)) nonconv <<- TRUE
      # glmnet >= 5.1 Cox tie-handling notice; multiview pins Breslow internally.
      if (grepl("converge", msg, fixed = TRUE) || grepl("cox.ties", msg, fixed = TRUE)) invokeRestart("muffleWarning")
    }
  )
  attr(fit, "zentangler_nonconverged") <- nonconv
  fit
}

zentangler_multiview_coef <- function(fit, s, n_med) {
  # Coefficients by position: intercept first (absent for Cox), then mediators
  # in x_list order. multiview's coefficient row names are not stable across
  # families, so names are not used.
  co <- as.numeric(as.matrix(stats::coef(fit, s = s))[, 1])
  if (length(co) == n_med + 1L) list(a0 = co[[1L]], beta = co[-1L]) else list(a0 = 0, beta = co[seq_len(n_med)])
}

fit_y_multiview_cooperative <- function(
  Y,
  X,
  blocks,
  C = NULL,
  y_family = c("gaussian", "binomial", "survival"),
  lambda_choice = c("lambda.1se", "lambda.min"),
  glmnet_alpha = 1,
  rho = NULL,
  rho_grid = zentangler_default_rho_grid(),
  nfolds = 10L,
  backfit_iter = 3L,
  ...
) {
  y_family <- match.arg(y_family)
  lambda_choice <- match.arg(lambda_choice)
  if (!requireNamespace("multiview", quietly = TRUE)) {
    stop("Package 'multiview' is required for fusion = 'multiview'. Install it with install.packages('multiview').", call. = FALSE)
  }
  n <- length(as.numeric(X))
  W <- zentangler_y_design(X, C)
  W1 <- cbind(1, W)
  qrW <- qr(W1)
  views <- names(blocks)
  x_list <- lapply(views, function(v) {
    M <- as.matrix(blocks[[v]])
    Mt <- M - W1 %*% qr.coef(qrW, M)
    colnames(Mt) <- colnames(M)
    Mt
  })
  names(x_list) <- views
  n_med <- sum(vapply(x_list, ncol, integer(1)))

  if (!is.null(rho)) {
    rho_grid <- as.numeric(rho)
  } else if (length(views) < 2L) {
    rho_grid <- 0
  }
  rho_grid <- sort(unique(as.numeric(rho_grid)))
  if (length(rho_grid) == 0L || any(!is.finite(rho_grid)) || any(rho_grid < 0 | rho_grid > 1)) {
    stop("rho / rho_grid must lie in [0, 1] (0 = early fusion, 1 = late fusion).", call. = FALSE)
  }

  # Offset from the unpenalized exposure/covariate model. Because the
  # mediators are orthogonal to W, this offset is close to the joint solution.
  base <- zentangler_unpenalized_fit(Y, W, y_family)
  offset <- base$intercept + as.numeric(W %*% base$coef)

  foldid <- zentangler_make_foldid(n, nfolds)
  K <- max(foldid)
  cv_tab <- list()
  paths <- list()
  dropped_rho <- numeric(0)
  for (rho_i in rho_grid) {
    full <- zentangler_multiview_path(x_list, Y, rho_i, y_family, offset, alpha = glmnet_alpha)
    lam <- full$lambda
    # Guard against non-converged solution paths (multiview 1.0 Cox fits with
    # rho > 0 currently stop after a handful of lambdas without converging).
    if (rho_i > 0 && (isTRUE(attr(full, "zentangler_nonconverged")) || length(lam) < 15L)) {
      dropped_rho <- c(dropped_rho, rho_i)
      next
    }
    err <- matrix(NA_real_, nrow = K, ncol = length(lam))
    for (k in seq_len(K)) {
      te <- which(foldid == k)
      tr <- which(foldid != k)
      fk <- try(zentangler_multiview_path(
        lapply(x_list, function(M) M[tr, , drop = FALSE]), zentangler_subset_y(Y, tr), rho_i,
        y_family, offset[tr], lambda = lam, alpha = glmnet_alpha
      ), silent = TRUE)
      if (inherits(fk, "try-error")) next
      x_te <- lapply(x_list, function(M) M[te, , drop = FALSE])
      eta <- stats::predict(fk, newx = x_te, s = lam, newoffset = offset[te], type = "link")
      eta <- as.matrix(eta)
      for (j in seq_len(min(ncol(eta), length(lam)))) {
        err[k, j] <- zentangler_cv_loss(zentangler_subset_y(Y, te), eta[, j], y_family)
      }
    }
    cvm <- colMeans(err, na.rm = TRUE)
    cvsd <- apply(err, 2, stats::sd, na.rm = TRUE) / sqrt(colSums(is.finite(err)))
    ok <- is.finite(cvm)
    i_min <- which(ok)[which.min(cvm[ok])]
    i_1se <- which(ok & cvm <= cvm[i_min] + cvsd[i_min])
    i_1se <- if (length(i_1se) > 0L) min(i_1se) else i_min
    cv_tab[[length(cv_tab) + 1L]] <- data.frame(
      rho = rho_i, cv_min = cvm[i_min], lambda.min = lam[i_min], lambda.1se = lam[i_1se], stringsAsFactors = FALSE
    )
    paths[[as.character(rho_i)]] <- full
  }
  if (length(dropped_rho) > 0L) {
    warning(
      "multiview did not converge for rho = ", paste(dropped_rho, collapse = ", "),
      if (identical(y_family, "survival")) " (known multiview 1.0 issue for Cox models with rho > 0)" else "",
      "; these values were removed from the rho grid.",
      call. = FALSE
    )
  }
  if (length(cv_tab) == 0L) stop("No rho value produced a converged multiview fit.", call. = FALSE)
  cv_tab <- do.call(rbind, cv_tab)
  best <- which.min(cv_tab$cv_min)
  rho_hat <- cv_tab$rho[best]
  lam_hat <- cv_tab[[lambda_choice]][best]

  # Backfit the unpenalized W block at the chosen (rho, lambda).
  fit <- paths[[as.character(rho_hat)]]
  for (it in seq_len(max(1L, backfit_iter))) {
    cf <- zentangler_multiview_coef(fit, lam_hat, n_med)
    med_lp <- as.numeric(do.call(cbind, x_list) %*% cf$beta)
    base <- zentangler_unpenalized_fit(Y, W, y_family, offset = med_lp)
    offset <- base$intercept + as.numeric(W %*% base$coef)
    fit <- zentangler_multiview_path(x_list, Y, rho_hat, y_family, offset, lambda = paths[[as.character(rho_hat)]]$lambda, alpha = glmnet_alpha)
  }
  cf <- zentangler_multiview_coef(fit, lam_hat, n_med)
  beta <- cf$beta
  sizes <- vapply(x_list, ncol, integer(1))
  idx <- split(seq_len(n_med), rep(views, sizes))[views]
  coef_by_view <- lapply(views, function(v) setNames(beta[idx[[v]]], colnames(x_list[[v]])))
  names(coef_by_view) <- views

  # Back to the original parametrization: W coefficients of Y ~ W + M.
  M_all <- do.call(cbind, lapply(views, function(v) as.matrix(blocks[[v]])))
  med_lp_raw <- as.numeric(M_all %*% beta)
  shift <- qr.coef(qrW, med_lp_raw)
  shift[!is.finite(shift)] <- 0
  w_coef <- base$coef - shift[-1L]
  names(w_coef) <- colnames(W)

  list(
    method = paste0("multiview cooperative learning (rho = ", format(rho_hat), if (length(rho_grid) > 1L) " by CV" else "", "; ", glmnet_alpha_label(glmnet_alpha), ")"),
    coef_by_view = coef_by_view,
    x_coef = as.numeric(w_coef[["X"]]),
    w_coef = w_coef,
    rho = rho_hat,
    lambda = lam_hat,
    fit = list(multiview = fit, cv = cv_tab, rho = rho_hat, lambda = lam_hat)
  )
}

# -----------------------------------------------------------------------------
# Supplementary: IntegratedLearner-style late fusion (per-view models + stacking)
#
# Each view k gets its own sparse outcome model Y ~ X + C + M_k. Its
# mediator-only linear predictor s_k = M_k b_k is cross-fitted (out-of-fold),
# and an unpenalized meta model Y ~ X + C + s_1 + ... + s_K estimates
# non-negative stacking weights w_k (views with negative weights are dropped
# and the meta model is refit). Everything lives on the link scale (identity,
# logit or log-hazard), so the implied mediator coefficient w_k * b_kj is on the
# same scale as the early-fusion coefficients, and X enters the meta model
# once rather than once per view.
# -----------------------------------------------------------------------------

fit_block_lasso_for_late <- function(
  Y,
  X,
  M,
  C = NULL,
  y_family = c("gaussian", "binomial", "survival"),
  lambda_choice = c("lambda.1se", "lambda.min"),
  glmnet_alpha = 1,
  foldid = NULL
) {
  y_family <- match.arg(y_family)
  lambda_choice <- match.arg(lambda_choice)
  glmnet_alpha <- validate_glmnet_alpha(glmnet_alpha)

  M <- as.matrix(M)
  W <- zentangler_y_design(X, C)
  Z <- cbind(W, M)
  med_cols <- (ncol(W) + 1L):ncol(Z)
  pf <- c(rep(0, ncol(W)), rep(1, ncol(M)))
  fam <- zentangler_glmnet_family(y_family)
  n <- nrow(Z)
  if (is.null(foldid)) foldid <- zentangler_make_foldid(n)

  cvfit <- glmnet::cv.glmnet(
    x = Z, y = Y, family = fam, alpha = glmnet_alpha,
    penalty.factor = pf, standardize = TRUE, foldid = foldid
  )
  lam <- as.numeric(cvfit[[lambda_choice]])
  beta_full <- as.matrix(coef(cvfit, s = lam))[, 1]
  intercept_name <- "(Intercept)"
  beta_full <- beta_full[setdiff(names(beta_full), intercept_name)]
  med_coef <- setNames(as.numeric(beta_full[med_cols]), colnames(M))

  # Cross-fitted mediator-only linear predictor.
  s_cf <- rep(NA_real_, n)
  for (k in sort(unique(foldid))) {
    test <- which(foldid == k)
    train <- which(foldid != k)
    fit_k <- try(
      glmnet::glmnet(
        x = Z[train, , drop = FALSE], y = zentangler_subset_y(Y, train), family = fam,
        alpha = glmnet_alpha, penalty.factor = pf, standardize = TRUE
      ),
      silent = TRUE
    )
    if (inherits(fit_k, "try-error")) next
    b_k <- as.matrix(glmnet::coef.glmnet(fit_k, s = lam))[, 1]
    b_k <- b_k[setdiff(names(b_k), intercept_name)]
    s_cf[test] <- as.numeric(M[test, , drop = FALSE] %*% b_k[med_cols])
  }
  s_cf[!is.finite(s_cf)] <- 0

  x_coef <- as.numeric(beta_full[[1L]])
  list(cvfit = cvfit, lambda = lam, score_cf = s_cf, x_coef = x_coef, med_coef = med_coef)
}

zentangler_fit_meta_model <- function(Y, W, S, y_family) {
  dat <- as.data.frame(cbind(W, S), check.names = FALSE)
  colnames(dat) <- make.names(colnames(dat), unique = TRUE)
  if (identical(y_family, "survival")) {
    fit <- survival::coxph(Y ~ ., data = dat)
  } else if (identical(y_family, "gaussian")) {
    fit <- stats::lm(Y ~ ., data = cbind(Y = as.numeric(Y), dat))
  } else {
    fit <- suppressWarnings(stats::glm(Y ~ ., data = cbind(Y = as.numeric(Y), dat), family = stats::binomial()))
  }
  co <- stats::coef(fit)
  names(co) <- gsub("`", "", names(co), fixed = TRUE)
  map <- setNames(colnames(cbind(W, S)), colnames(dat))
  co <- co[names(co) %in% names(map)]
  names(co) <- map[names(co)]
  list(fit = fit, coef = co)
}

fit_y_multiview_late_fusion_lasso <- function(
  Y,
  X,
  blocks,
  C = NULL,
  y_family = c("gaussian", "binomial", "survival"),
  lambda_choice = c("lambda.1se", "lambda.min"),
  glmnet_alpha = 1
) {
  y_family <- match.arg(y_family)
  lambda_choice <- match.arg(lambda_choice)
  glmnet_alpha <- validate_glmnet_alpha(glmnet_alpha)

  n <- length(as.numeric(X))
  foldid <- zentangler_make_foldid(n)
  block_fits <- lapply(names(blocks), function(view) {
    fit_block_lasso_for_late(
      Y = Y, X = X, M = blocks[[view]], C = C, y_family = y_family,
      lambda_choice = lambda_choice, glmnet_alpha = glmnet_alpha, foldid = foldid
    )
  })
  names(block_fits) <- names(blocks)

  W <- zentangler_y_design(X, C)
  S_all <- do.call(cbind, lapply(block_fits, `[[`, "score_cf"))
  S_all <- matrix(S_all, nrow = n)
  colnames(S_all) <- paste0("score_", names(blocks))

  weights <- setNames(rep(0, length(blocks)), names(blocks))
  keep <- names(blocks)[apply(S_all, 2, stats::sd) > 1e-10]
  meta <- NULL
  while (length(keep) > 0L) {
    meta <- zentangler_fit_meta_model(Y, W, S_all[, paste0("score_", keep), drop = FALSE], y_family)
    w_keep <- meta$coef[paste0("score_", keep)]
    w_keep[!is.finite(w_keep)] <- -1
    if (all(w_keep >= 0)) {
      weights[keep] <- as.numeric(w_keep)
      break
    }
    # Non-negative stacking: drop the most negative view and refit.
    keep <- keep[-which.min(w_keep)]
  }
  if (length(keep) == 0L) {
    meta <- zentangler_fit_meta_model(Y, W, NULL, y_family)
  }

  coef_by_view <- lapply(names(blocks), function(view) block_fits[[view]]$med_coef * weights[[view]])
  names(coef_by_view) <- names(blocks)
  x_coef <- as.numeric(meta$coef[["X"]] %||% 0)
  if (!is.finite(x_coef)) x_coef <- 0

  list(
    method = paste0("IntegratedLearner-style late fusion, cross-fitted stacking (", glmnet_alpha_label(glmnet_alpha), ")"),
    coef_by_view = coef_by_view,
    x_coef = x_coef,
    w_coef = setNames(as.numeric(meta$coef[colnames(W)]), colnames(W)),
    fit = list(blocks = lapply(block_fits, `[[`, "cvfit"), meta = meta$fit),
    meta_weights = weights
  )
}

zentangler_resolve_fusion <- function(fusion, integrated_learner_fusion = "late") {
  # User-facing choice -> internal outcome-model engine.
  #   fusion = "multiview"            cooperative learning, rho tuned in [0, 1] (main)
  #   fusion = "integrated_learner"   supplementary IntegratedLearner-style fusion,
  #                                   integrated_learner_fusion = "late" or "early"
  fusion <- match.arg(as.character(fusion)[1], c("multiview", "integrated_learner"))
  if (identical(fusion, "multiview")) return("multiview")
  il <- match.arg(as.character(integrated_learner_fusion)[1], c("late", "early"))
  paste0("il_", il)
}

fit_y_multiview_stage <- function(
  Y,
  X,
  blocks,
  C = NULL,
  engine = c("multiview", "il_late", "il_early"),
  y_family = c("gaussian", "binomial", "survival"),
  lambda_choice = c("lambda.1se", "lambda.min"),
  glmnet_alpha = 1,
  rho = NULL,
  rho_grid = zentangler_default_rho_grid()
) {
  engine <- match.arg(engine)
  y_family <- match.arg(y_family)
  lambda_choice <- match.arg(lambda_choice)
  glmnet_alpha <- validate_glmnet_alpha(glmnet_alpha)

  n_med_total <- sum(vapply(blocks, function(M) ncol(as.matrix(M)), integer(1)))
  if (identical(engine, "multiview") && n_med_total < 2L) {
    # multiview needs at least two penalized columns; with a single mediator
    # cooperative learning reduces to the single-model (early) fit.
    engine <- "il_early"
  }
  if (identical(engine, "multiview")) {
    return(fit_y_multiview_cooperative(
      Y = Y,
      X = X,
      blocks = blocks,
      C = C,
      y_family = y_family,
      lambda_choice = lambda_choice,
      glmnet_alpha = glmnet_alpha,
      rho = rho,
      rho_grid = rho_grid
    ))
  }

  if (identical(engine, "il_late")) {
    return(fit_y_multiview_late_fusion_lasso(
      Y = Y,
      X = X,
      blocks = blocks,
      C = C,
      y_family = y_family,
      lambda_choice = lambda_choice,
      glmnet_alpha = glmnet_alpha
    ))
  }

  fit_y_multiview_glmnet_lasso(
    Y = Y,
    X = X,
    blocks = blocks,
    C = C,
    y_family = y_family,
    lambda_choice = lambda_choice,
    glmnet_alpha = glmnet_alpha
  )
}
