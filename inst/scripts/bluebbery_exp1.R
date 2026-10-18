library(Zentangler)
library(MultiAssayExperiment)

# -----------------------------
# Load MAEs
# -----------------------------
load("~/Desktop/Zentangler/data/Blueberry_MAE_pathways.RData")   # loads object: mae
mae_pathways <- mae

load("~/Desktop/Zentangler/data/Blueberry_MAE_ecs.RData")
mae_ecs <- mae

load("~/Desktop/Zentangler/data/Blueberry_MAE_all_omics.RData")
mae_all <- mae

# -----------------------------
# Covariates
# -----------------------------
covars <- c(
  "Age",
  "sex",
  "bmi",
  "bp_med",
  "statin",
  "metabolicsyndromecriteria"
)

# -----------------------------
# Outcomes you may want to run
# -----------------------------
primary_outcomes <- c(
  "cGMP_change",
  "cGMP_month6"
)

secondary_change_outcomes <- c(
  "FMD_change",
  "AIx75bpm_change",
  "PWV_change",
  "SBP_change",
  "DBP_change",
  "Trigl_change",
  "HDLC_change",
  "LDLC_change",
  "CRP_ug_ml20_change",
  "IL6_pg_ml20_change",
  "TNF.._pg_ml20_change"
)

all_outcomes <- c(secondary_change_outcomes)

# Map a fusion label to the zentangler outcome-model arguments:
#   "multiview" (default; cooperative learning, rho in [0, 1] tuned by CV)
#   "early" / "late" (supplementary IntegratedLearner-style fusion)
fusion_args <- function(label) {
  label <- as.character(label)[1]
  if (identical(label, "multiview")) {
    return(list(fusion = "multiview", integrated_learner_fusion = "late"))
  }
  if (label %in% c("early", "late")) {
    return(list(fusion = "integrated_learner", integrated_learner_fusion = label))
  }
  stop("Unknown fusion label: ", label, call. = FALSE)
}

# -----------------------------
# Parallel Zentangler helper
# -----------------------------
run_parallel_blueberry <- function(mae, outcome, fusion_label = "early") {
  fa <- fusion_args(fusion_label)
  outcome_covars <- covars
  if (identical(outcome, "cGMP_month6")) {
    outcome_covars <- c("cGMP_baseline", outcome_covars)
  }
  
  fit_multiview_parallel_zentangler(
    mae = mae,
    x_var = "treatment_any",
    y_var = outcome,
    y_family = "gaussian",
    covariates = outcome_covars,
    screen_method = "sis",
    sis_n = 25,
    glmnet_alpha = 0.5,
    lambda_choice = "lambda.min",
    fusion = fa$fusion,
    integrated_learner_fusion = fa$integrated_learner_fusion,
    b_inference = "debiased",
    fdr_method = "BH",
    fdr_scope = "global"
  )
}

# -----------------------------
# Run parallel models
# -----------------------------
mae_list <- list(
  pathways = mae_pathways,
  ecs = mae_ecs,
  all_omics = mae_all
)

fusion_labels <- c("early", "late")

parallel_fits <- list()

for (mae_name in names(mae_list)) {
  parallel_fits[[mae_name]] <- list()
  for (fusion in fusion_labels) {
    parallel_fits[[mae_name]][[fusion]] <- lapply(all_outcomes, function(outcome) {
      run_parallel_blueberry(
        mae = mae_list[[mae_name]],
        outcome = outcome,
        fusion_label = fusion
      )
    })
    names(parallel_fits[[mae_name]][[fusion]]) <- all_outcomes
  }
}

# Example:
# zentangler_mediators(parallel_fits[["all_omics"]][["early"]][["Trigl_change"]], q_threshold = 0.25)
# zentangler_mediators(parallel_fits[["pathways"]][["early"]][["cGMP_change"]], top = 20)
# summary(parallel_fits[["ecs"]][["late"]][["IL6_pg_ml20_change"]], q_threshold = 0.25)$views

# -----------------------------
# Sequential Zentangler helper
# -----------------------------
run_sequential_blueberry <- function(
    mae,
    outcome,
    path_template,
    fusion_label = "early"
) {
  fa <- fusion_args(fusion_label)
  outcome_covars <- covars
  if (identical(outcome, "cGMP_month6")) {
    outcome_covars <- c("cGMP_baseline", outcome_covars)
  }
  
  fit_sequential_zentangler(
    mae = mae,
    x_var = "treatment_any",
    y_var = outcome,
    route = path_template,
    covariates = outcome_covars,
    sis_n = 100,
    cor_method = "spearman",
    min_abs_cor = 0.2,
    cor_q_threshold = 0.25,
    y_family = "gaussian",
    fusion = fa$fusion,
    integrated_learner_fusion = fa$integrated_learner_fusion,
    lambda_choice = "lambda.min",
    glmnet_alpha = 0.5,
    b_inference = "debiased"
  )
}

# -----------------------------
# Sequential routes
# taxa -> function -> outcome
# -----------------------------
seq_routes <- list(
  species_dna_pathways = c("species", "dna_pathways"),
  species_rna_pathways = c("species", "rna_pathways"),
  species_dna_ecs = c("species", "dna_ecs"),
  species_rna_ecs = c("species", "rna_ecs")
)

# Use the full MAE so all route views are available.
seq_outcomes <- c(
  "cGMP_change",
  "FMD_change",
  "AIx75bpm_change",
  "PWV_change",
  "SBP_change",
  "DBP_change",
  "Trigl_change",
  "HDLC_change",
  "LDLC_change",
  "CRP_ug_ml20_change",
  "IL6_pg_ml20_change",
  "TNF.._pg_ml20_change"
)

seq_fits <- list()
# Sequential fits: supplementary early fusion and the default multiview model.
seq_fusion_labels <- c("early", "multiview")

for (route_name in names(seq_routes)) {
  seq_fits[[route_name]] <- list()
  for (fusion in seq_fusion_labels) {
    seq_fits[[route_name]][[fusion]] <- lapply(seq_outcomes, function(outcome) {
      run_sequential_blueberry(
        mae = mae_all,
        outcome = outcome,
        path_template = seq_routes[[route_name]],
        fusion_label = fusion
      )
    })
    names(seq_fits[[route_name]][[fusion]]) <- seq_outcomes
  }
}

# Example:
# summary(seq_fits[["species_dna_pathways"]][["early"]][["cGMP_change"]])
# zentangler_paths(seq_fits[["species_rna_ecs"]][["multiview"]][["Trigl_change"]])
# seq_fits[["species_dna_ecs"]][["multiview"]][["IL6_pg_ml20_change"]]$edges

# -----------------------------
# Optional: quick extractor
# -----------------------------
get_parallel_active <- function(fit) {
  zentangler_mediators(fit, q_threshold = 0.25)
}

get_sequential_paths <- function(fit) {
  zentangler_paths(fit)
}