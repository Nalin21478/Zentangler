# -----------------------------------------------------------------------------
# B-stage inference
#
# Two inference engines are supported everywhere in the package:
#
#   b_inference = "debiased"   De-biased (de-sparsified) lasso on the joint
#                              outcome model Y ~ X + C + all screened mediators.
#                              Every screened mediator is a target, not only the
#                              lasso-active ones, so p_b exists for the whole
#                              screened set (HIMA2-style). Gaussian and logistic
#                              outcomes use nodewise-lasso de-biasing (van de
#                              Geer et al. 2014); Cox outcomes use the
#                              de-sparsified partial-likelihood estimator
#                              (Fang, Ning & Liu 2017; Kong et al. 2021) with a
#                              nodewise lasso on the partial-likelihood
#                              information matrix.
#
#   b_inference = "bootstrap"  Full-pipeline bootstrap. Each replicate re-runs
#                              screening and the outcome model; p_b is the
#                              percentile-interval inversion p-value of the
#                              fusion coefficient (see 07-bootstrap.R).
#
# Naive active-set refits (OLS / GLM / Cox Wald on the lasso-selected set) have
# been removed: they ignore selection and are anti-conservative.
# -----------------------------------------------------------------------------

zentangler_validate_b_inference <- function(method) {
  method <- as.character(method)[1]
  if (identical(method, "refit")) {
    stop("b_inference = 'refit' was removed because active-set refits ignore lasso selection. Use 'debiased' or 'bootstrap'.", call. = FALSE)
  }
  if (method %in% c("debiased_lasso", "debiased_logistic_lasso", "debiased_cox_lasso")) {
    stop("b_inference = '", method, "' is now spelled 'debiased' (the family is taken from y_family).", call. = FALSE)
  }
  # "none" is internal: bootstrap replicates need point estimates only.
  match.arg(method, c("debiased", "bootstrap", "none"))
}

pick_lambda_from_cv <- function(cvfit, lambda_choice = c("lambda.1se", "lambda.min")) {
  lambda_choice <- match.arg(lambda_choice)
  lam <- suppressWarnings(as.numeric(cvfit[[lambda_choice]]))
  if (!is.finite(lam) || lam <= 0) {
    lam_grid <- suppressWarnings(as.numeric(cvfit$lambda))
    lam <- if (length(lam_grid) > 0) max(min(lam_grid), 1e-06) else 1e-04
  }
  lam
}

zentangler_standardize_columns <- function(Z, weights = NULL) {
  Z <- as.matrix(Z)
  if (is.null(weights)) {
    ctr <- colMeans(Z)
  } else {
    ctr <- colSums(Z * weights) / sum(weights)
  }
  sc <- apply(Z, 2, stats::sd)
  sc[!is.finite(sc) | sc < 1e-10] <- 1
  list(Z = sweep(sweep(Z, 2, ctr, "-"), 2, sc, "/"), center = ctr, scale = sc)
}

# Nodewise regression of column j on the remaining columns (possibly
# row-weighted through sqrt_w). Returns the residual r_j on the unweighted scale.
zentangler_nodewise_residual <- function(Zw, j, Z = Zw, lambda_choice = "lambda.min", nfolds = 5L) {
  z_j <- Zw[, j]
  z_m <- Zw[, -j, drop = FALSE]
  if (ncol(z_m) == 0L) return(Z[, j])
  if (ncol(z_m) == 1L) {
    denom <- sum(z_m^2)
    gamma <- if (denom > 0) sum(z_m * z_j) / denom else 0
  } else {
    cv_j <- glmnet::cv.glmnet(
      x = z_m, y = z_j, family = "gaussian", alpha = 1,
      standardize = FALSE, intercept = FALSE, nfolds = nfolds
    )
    lam_j <- pick_lambda_from_cv(cv_j, lambda_choice = lambda_choice)
    gamma <- as.numeric(glmnet::coef.glmnet(cv_j, s = lam_j))[-1]
  }
  Z[, j] - as.numeric(Z[, -j, drop = FALSE] %*% gamma)
}

zentangler_order_targets <- function(targets, beta_hat, max_targets) {
  max_targets <- max(1L, as.integer(max_targets))
  if (length(targets) <= max_targets) return(list(keep = targets, dropped = character(0)))
  ord <- order(abs(beta_hat[targets]), decreasing = TRUE, na.last = TRUE)
  list(keep = targets[ord[seq_len(max_targets)]], dropped = targets[ord[-seq_len(max_targets)]])
}

# Coordinate-descent lasso on a Gram matrix:
#   min_g (1/2) g' S g - c' g + lambda ||g||_1
zentangler_gram_lasso <- function(S, cvec, lambda, maxit = 500L, tol = 1e-07) {
  p <- length(cvec)
  g <- rep(0, p)
  if (p == 0L) return(g)
  d <- diag(S)
  d[d <= 0] <- 1e-10
  for (it in seq_len(maxit)) {
    g_old <- g
    for (k in seq_len(p)) {
      rho_k <- cvec[k] - sum(S[k, -k] * g[-k])
      g[k] <- sign(rho_k) * max(abs(rho_k) - lambda, 0) / d[k]
    }
    if (max(abs(g - g_old)) < tol) break
  }
  g
}

# Partial-likelihood gradient (of the average log partial likelihood) and
# information matrix at eta = Z beta, Breslow ties.
zentangler_cox_grad_info <- function(time, status, Z, eta) {
  n <- nrow(Z)
  p <- ncol(Z)
  ord <- order(time, decreasing = TRUE)
  time_o <- time[ord]
  status_o <- status[ord]
  Z_o <- Z[ord, , drop = FALSE]
  r_o <- exp(eta[ord] - max(eta))
  S0 <- 0
  S1 <- rep(0, p)
  S2 <- matrix(0, p, p)
  grad <- rep(0, p)
  info <- matrix(0, p, p)
  ut <- unique(time_o)
  pos <- 1L
  for (t in ut) {
    idx <- integer(0)
    while (pos <= n && time_o[pos] == t) {
      idx <- c(idx, pos)
      pos <- pos + 1L
    }
    zi <- Z_o[idx, , drop = FALSE]
    ri <- r_o[idx]
    S0 <- S0 + sum(ri)
    S1 <- S1 + colSums(zi * ri)
    S2 <- S2 + crossprod(zi * sqrt(ri))
    ev <- idx[status_o[idx] == 1]
    d_t <- length(ev)
    if (d_t == 0L) next
    zbar <- S1 / S0
    grad <- grad + colSums(Z_o[ev, , drop = FALSE]) - d_t * zbar
    info <- info + d_t * (S2 / S0 - tcrossprod(zbar))
  }
  list(grad = grad / n, info = info / n)
}

# Generic de-biased lasso for the joint outcome model.
#   Y          outcome (numeric, 0/1, or Surv)
#   W          unpenalized design (X and covariates)
#   P          penalized design (mediators), colnames are target names
#   targets    which columns of P to de-bias
# Returns a data.frame(name, estimate, se, p) on the original column scale.
zentangler_debiased_lasso <- function(
  Y,
  W,
  P,
  targets = colnames(P),
  y_family = c("gaussian", "binomial", "survival"),
  lambda_choice = c("lambda.1se", "lambda.min"),
  max_targets = 200L,
  beta_init = NULL
) {
  # beta_init: optional initial estimator on the original column scale, named
  # by colnames(cbind(W, P)) (missing names are treated as 0). When supplied,
  # the one-step de-biasing correction is applied to it instead of to a fresh
  # lasso; the de-biased estimator targets the same coefficient of the joint
  # model, so the fused (e.g. cooperative-learning) fit can be used directly.
  y_family <- match.arg(y_family)
  lambda_choice <- match.arg(lambda_choice)
  empty <- data.frame(name = character(0), estimate = numeric(0), se = numeric(0), p = numeric(0), stringsAsFactors = FALSE)
  W <- as.matrix(W)
  P <- as.matrix(P)
  if (ncol(P) == 0L || length(targets) == 0L) return(list(table = empty, dropped = character(0), method = "no_targets"))

  Z_raw <- cbind(W, P)
  n <- nrow(Z_raw)
  p <- ncol(Z_raw)
  min_n <- if (identical(y_family, "gaussian")) 20L else 30L
  if (n < min_n) {
    return(list(table = empty, dropped = targets, method = paste0("debiased_", y_family, "_too_few_samples")))
  }
  pf <- c(rep(0, ncol(W)), rep(1, ncol(P)))
  st <- zentangler_standardize_columns(Z_raw)
  Z <- st$Z
  colnames(Z) <- colnames(Z_raw)
  out_rows <- list()
  init_std <- NULL
  if (!is.null(beta_init)) {
    init_std <- setNames(rep(0, ncol(Z)), colnames(Z))
    hit <- intersect(names(beta_init), colnames(Z))
    init_std[hit] <- as.numeric(beta_init[hit])
    init_std[!is.finite(init_std)] <- 0
    init_std <- init_std * st$scale
  }

  if (identical(y_family, "gaussian")) {
    y0 <- as.numeric(Y) - mean(as.numeric(Y))
    if (is.null(init_std)) {
      cv_main <- glmnet::cv.glmnet(x = Z, y = y0, family = "gaussian", alpha = 1, penalty.factor = pf, standardize = FALSE, intercept = FALSE)
      lam <- pick_lambda_from_cv(cv_main, lambda_choice)
      beta_hat <- as.numeric(glmnet::coef.glmnet(cv_main, s = lam))[-1]
      names(beta_hat) <- colnames(Z)
    } else {
      beta_hat <- init_std
    }
    resid <- y0 - as.numeric(Z %*% beta_hat)
    s_hat <- sum(abs(beta_hat) > 0)
    sigma2 <- max(sum(resid^2) / max(1, n - s_hat - 1), 1e-10)
    tg <- zentangler_order_targets(intersect(targets, colnames(Z)), beta_hat, max_targets)
    for (nm in tg$keep) {
      j <- match(nm, colnames(Z))
      r_j <- try(zentangler_nodewise_residual(Z, j), silent = TRUE)
      if (inherits(r_j, "try-error")) next
      denom <- sum(Z[, j] * r_j)
      if (!is.finite(denom) || abs(denom) < 1e-10) next
      est <- beta_hat[[nm]] + sum(r_j * resid) / denom
      se <- sqrt(sigma2 * sum(r_j^2)) / abs(denom)
      out_rows[[nm]] <- c(est, se)
    }
    method <- "debiased_lasso_gaussian"
  } else if (identical(y_family, "binomial")) {
    y_raw <- as.numeric(Y)
    yu <- sort(unique(y_raw[is.finite(y_raw)]))
    if (length(yu) != 2L || !all(yu %in% c(0, 1))) stop("Binomial outcomes must be coded 0/1.", call. = FALSE)
    y <- y_raw
    if (is.null(init_std)) {
      cv_main <- glmnet::cv.glmnet(x = Z, y = y, family = "binomial", alpha = 1, penalty.factor = pf, standardize = FALSE, intercept = TRUE)
      lam <- pick_lambda_from_cv(cv_main, lambda_choice)
      co <- as.numeric(glmnet::coef.glmnet(cv_main, s = lam))
      a0 <- co[1]
      beta_hat <- co[-1]
      names(beta_hat) <- colnames(Z)
    } else {
      beta_hat <- init_std
      lp <- as.numeric(Z %*% beta_hat)
      a0 <- suppressWarnings(stats::glm.fit(matrix(1, nrow = n), y, family = stats::binomial(), offset = lp)$coefficients[[1L]])
    }
    mu <- stats::plogis(a0 + as.numeric(Z %*% beta_hat))
    w <- pmax(mu * (1 - mu), 1e-05)
    score <- y - mu
    # Weighted centering projects out the intercept direction, so the nodewise
    # step can run without an intercept on the sqrt(w)-weighted design.
    wc <- sweep(Z, 2, colSums(Z * w) / sum(w), "-")
    Zw <- wc * sqrt(w)
    tg <- zentangler_order_targets(intersect(targets, colnames(Z)), beta_hat, max_targets)
    for (nm in tg$keep) {
      j <- match(nm, colnames(Z))
      r_j <- try(zentangler_nodewise_residual(Zw, j, Z = wc), silent = TRUE)
      if (inherits(r_j, "try-error")) next
      tau <- mean(w * wc[, j] * r_j)
      if (!is.finite(tau) || abs(tau) < 1e-10) next
      est <- beta_hat[[nm]] + mean(r_j * score) / tau
      se <- sqrt(mean(w * r_j^2) / (n * tau^2))
      out_rows[[nm]] <- c(est, se)
    }
    method <- "debiased_lasso_logistic"
  } else {
    if (!inherits(Y, "Surv")) stop("Survival de-biased inference requires a Surv outcome.", call. = FALSE)
    time <- as.numeric(Y[, 1])
    status <- as.numeric(Y[, 2])
    if (is.null(init_std)) {
      cv_main <- glmnet::cv.glmnet(x = Z, y = Y, family = "cox", alpha = 1, penalty.factor = pf, standardize = FALSE)
      lam <- pick_lambda_from_cv(cv_main, lambda_choice)
      beta_hat <- as.numeric(glmnet::coef.glmnet(cv_main, s = lam))
      names(beta_hat) <- colnames(Z)
    } else {
      beta_hat <- init_std
    }
    gi <- zentangler_cox_grad_info(time, status, Z, as.numeric(Z %*% beta_hat))
    Sigma <- gi$info
    lam_node <- sqrt(log(max(p, 2)) / n)
    tg <- zentangler_order_targets(intersect(targets, colnames(Z)), beta_hat, max_targets)
    for (nm in tg$keep) {
      j <- match(nm, colnames(Z))
      gam <- zentangler_gram_lasso(Sigma[-j, -j, drop = FALSE], Sigma[-j, j], lambda = lam_node * sqrt(max(Sigma[j, j], 1e-10)))
      tau2 <- Sigma[j, j] - sum(Sigma[j, -j] * gam)
      if (!is.finite(tau2) || tau2 <= 1e-10) next
      theta <- rep(0, p)
      theta[j] <- 1
      theta[-j] <- -gam
      theta <- theta / tau2
      est <- beta_hat[[nm]] + sum(theta * gi$grad)
      v <- as.numeric(t(theta) %*% Sigma %*% theta) / n
      if (!is.finite(v) || v <= 0) next
      out_rows[[nm]] <- c(est, sqrt(v))
    }
    method <- "debiased_lasso_cox"
  }

  if (length(out_rows) == 0L) return(list(table = empty, dropped = tg$dropped, method = method))
  mat <- do.call(rbind, out_rows)
  nm <- rownames(mat)
  sc <- st$scale[nm]
  est <- mat[, 1] / sc
  se <- mat[, 2] / sc
  tab <- data.frame(
    name = nm,
    estimate = as.numeric(est),
    se = as.numeric(se),
    p = as.numeric(2 * stats::pnorm(-abs(mat[, 1] / mat[, 2]))),
    stringsAsFactors = FALSE
  )
  list(table = tab, dropped = tg$dropped, method = method)
}

zentangler_fusion_init <- function(yfit, blocks) {
  # Initial estimator for the de-biasing step, taken from the fused outcome
  # model: unpenalized W coefficients plus view-prefixed mediator coefficients.
  if (is.null(yfit$w_coef) || is.null(yfit$coef_by_view)) return(NULL)
  med <- unlist(lapply(names(blocks), function(v) {
    b <- yfit$coef_by_view[[v]]
    setNames(as.numeric(b), paste0(v, "::", names(b)))
  }), use.names = TRUE)
  out <- c(yfit$w_coef, med)
  out[!is.finite(out)] <- 0
  out
}

empty_p_b_multiview <- function(blocks) {
  lapply(blocks, function(M) setNames(rep(NA_real_, ncol(M)), colnames(M)))
}

infer_p_b_multiview <- function(
  Y,
  X,
  blocks,
  C = NULL,
  coef_by_view = NULL,
  y_family = c("gaussian", "binomial", "survival"),
  method = c("debiased", "bootstrap", "none"),
  lambda_choice = c("lambda.1se", "lambda.min"),
  max_debias_targets = 200L,
  extra_penalized = NULL,
  beta_init = NULL
) {
  # Returns per-view lists of de-biased estimates (b), standard errors and
  # p-values for every screened mediator. `extra_penalized` adds penalized
  # adjustment columns (used by the sequential terminal model for upstream
  # mediators) that are part of the joint model but are not targets.
  y_family <- match.arg(y_family)
  method <- zentangler_validate_b_inference(method)
  lambda_choice <- match.arg(lambda_choice)
  out_p <- empty_p_b_multiview(blocks)
  out_b <- empty_p_b_multiview(blocks)
  out_se <- empty_p_b_multiview(blocks)

  if (!identical(method, "debiased")) {
    return(list(p_b = out_p, b = out_b, se = out_se, method = if (identical(method, "bootstrap")) "bootstrap_b_path" else "none", n_untested = 0L))
  }

  pref <- prefix_multiview_blocks(blocks)
  P <- pref$X
  if (!is.null(extra_penalized) && ncol(extra_penalized) > 0L) P <- cbind(P, as.matrix(extra_penalized))
  W <- zentangler_y_design(X, C)
  db <- try(
    zentangler_debiased_lasso(
      Y = Y, W = W, P = P, targets = colnames(pref$X), y_family = y_family,
      lambda_choice = lambda_choice, max_targets = max_debias_targets,
      beta_init = beta_init
    ),
    silent = TRUE
  )
  if (inherits(db, "try-error")) {
    msg <- conditionMessage(attr(db, "condition"))
    warning("De-biased B-path inference failed; returning NA p-values. Reason: ", msg, call. = FALSE)
    return(list(p_b = out_p, b = out_b, se = out_se, method = paste0("debiased_failed: ", msg), n_untested = ncol(pref$X)))
  }
  if (length(db$dropped) > 0L) {
    warning(
      length(db$dropped), " screened mediators exceeded debias_max_targets and were not de-biased; ",
      "they enter the FDR step with p = 1. Increase debias_max_targets or reduce sis_n.",
      call. = FALSE
    )
  }
  tab <- db$table
  for (i in seq_len(nrow(tab))) {
    hit <- match(tab$name[i], pref$map$design_name)
    if (is.na(hit)) next
    view <- pref$map$view[hit]
    med <- pref$map$mediator[hit]
    out_p[[view]][med] <- tab$p[i]
    out_b[[view]][med] <- tab$estimate[i]
    out_se[[view]][med] <- tab$se[i]
  }
  list(p_b = out_p, b = out_b, se = out_se, method = db$method, n_untested = ncol(pref$X) - nrow(tab))
}

add_mediator_effect_columns <- function(combined, effect_decomposition, mediator_effects = NULL) {
  # Per-mediator indirect effects on the same scale as nie_total. They come
  # from the post-selection effect model (06-causal-effects.R) and exist only
  # for mediators in the effect set; other mediators get NA.
  if (is.null(combined) || nrow(combined) == 0) return(combined)
  combined$indirect_effect <- NA_real_
  combined$prop_nie <- NA_real_
  combined$in_effect_set <- FALSE
  if (!is.null(mediator_effects) && nrow(mediator_effects) > 0L) {
    key_tab <- paste(combined$omics, combined$mediator, sep = "::")
    key_eff <- paste(mediator_effects$omics, mediator_effects$mediator, sep = "::")
    hit <- match(key_tab, key_eff)
    ok <- !is.na(hit)
    combined$indirect_effect[ok] <- mediator_effects$indirect_effect[hit[ok]]
    combined$in_effect_set[ok] <- TRUE
    total_nie <- suppressWarnings(as.numeric(effect_decomposition$nie_total[1]))
    if (is.finite(total_nie) && abs(total_nie) > 1e-08) {
      combined$prop_nie[ok] <- combined$indirect_effect[ok] / total_nie
    }
  }
  combined$abs_indirect_effect <- abs(combined$indirect_effect)
  combined
}

assemble_mediator_table <- function(screen_tab, b_vec, b_fusion_vec, se_b_vec, p_b_vec, fdr_method = c("BH", "BY")) {
  fdr_method <- validate_fdr_method(fdr_method)
  med <- screen_tab$mediator
  selected <- as.logical(screen_tab$selected)

  out <- data.frame(
    mediator = med,
    a = screen_tab$a_hat,
    p_a = screen_tab$p_a,
    q_a = screen_tab$q_a,
    b = as.numeric(b_vec[med]),
    se_b = as.numeric(se_b_vec[med]),
    p_b = as.numeric(p_b_vec[med]),
    b_fusion = as.numeric(b_fusion_vec[med]),
    selected_by_screen = selected,
    selected_by_sis = selected,
    selected_by_y_model = is.finite(as.numeric(b_fusion_vec[med])) & as.numeric(b_fusion_vec[med]) != 0,
    stringsAsFactors = FALSE
  )
  out$b[!is.finite(out$b)] <- 0
  out$b_fusion[!is.finite(out$b_fusion)] <- 0
  out <- zentangler_finalize_joint_p(out)
  out$q_primary_within_view <- p.adjust(out$p_primary, method = fdr_method)
  out$q_primary <- out$q_primary_within_view

  out <- out[order(out$abs_score, decreasing = TRUE), ]
  rownames(out) <- NULL
  out
}

zentangler_finalize_joint_p <- function(tab) {
  # Joint-significance (max-P) test of H0: a = 0 or b = 0. Mediators that were
  # not tested in the B-stage (screened out, or beyond debias_max_targets) get
  # p_primary = 1 so that FDR is controlled over every tested feature rather
  # than over the selected subset.
  tab$score <- tab$a * tab$b
  tab$abs_score <- abs(tab$score)
  both <- is.finite(tab$p_a) & is.finite(tab$p_b)
  tab$joint_p_ab <- NA_real_
  tab$joint_p_ab[both] <- pmax(tab$p_a[both], tab$p_b[both])
  tab$p_primary <- tab$joint_p_ab
  tab$p_primary[!is.finite(tab$p_primary)] <- 1
  tab$tested_b_stage <- is.finite(tab$p_b)
  tab
}
