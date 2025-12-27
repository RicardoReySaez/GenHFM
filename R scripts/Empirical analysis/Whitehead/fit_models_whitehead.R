# ╔════════════════════════════════════════════════════════════════════════════╗
# ║                             SCRIPT OVERVIEW                                ║
# ╠════════════════════════════════════════════════════════════════════════════╣
# ║ Script Name   : fit_models_whitehead.R                                     ║
# ║ Author        : Ricardo Rey-Sáez                                           ║
# ║ Role          : PhD Student in Psychology                                  ║
# ║ Institution   : Autonomous University of Madrid, Spain                     ║
# ║ Email         : ricardoreysaez95@gmail.es                                  ║
# ║ Date          : 17-05-2025                                                 ║
# ╚════════════════════════════════════════════════════════════════════════════╝

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 1: Description
# ─────────────────────────────────────────────────────────────────────────────
# Script destinated to:
#   1. Fit gaussian, ex-Gaussian and sh-lognormal GHFMs for each experiment
#   2. Estimate Mixture-IS models
#   3. Compute ELPDs per model
# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 2: Load Packages
# ─────────────────────────────────────────────────────────────────────────────

library(cmdstanr)
library(posterior)

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 3: Global settings
# ─────────────────────────────────────────────────────────────────────────────

# Extra-functions to analyze data
source("R scripts/R functions/empirical_functions.R")

# Save all BGHFM in one list
BGHFM_models <- list(
  gaussian          = cmdstan_model("Stan models/BGHFM/C_BGHFM_gaussian.stan"),
  exgaussian_tau    = cmdstan_model("Stan models/BGHFM/C_BGHFM_exgaussian_tau.stan"),
  shifted_lognormal = cmdstan_model("Stan models/BGHFM/C_BGHFM_shlognormal.stan")
)

# Global Stan arguments
stan_arguments <- list(
  data = NULL, 
  iter_warmup = 1500,
  iter_sampling = 2000,
  chains = 4, 
  parallel_chains = 4,
  adapt_delta = 0.80,
  max_treedepth = 10,
  seed = 2025,
  refresh = 50,
  init = NULL)

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 4: Set prior parameters for all models (diffuse priors)
# ─────────────────────────────────────────────────────────────────────────────

# Same priors as in the simulation study
priors_all <- list(
  gaussian = list(
    # Population means
    pr_mu_alpha = matrix(c(0, 1), ncol = 2, nrow = 3, byrow = TRUE),
    pr_mu_beta = matrix(c(0, 1), ncol = 2, nrow = 3, byrow = TRUE),
    pr_shift_sigma = matrix(c(0, 1), ncol = 2, nrow = 3, byrow = TRUE),
    # Population std.devs
    pr_sd_alpha = matrix(c(4, 0, .5), ncol = 3, nrow = 3, byrow = TRUE),
    pr_sd_beta = matrix(c(4, 0, .5), ncol = 3, nrow = 3, byrow = TRUE),
    pr_scale_sigma = matrix(c(4, 0, .5), ncol = 3, nrow = 3, byrow = TRUE),
    # LKJ priors
    pr_L_R_alpha = 1,
    pr_L_Phi = 1,
    # Factor loadings
    pr_L_raw = c(1,1)
    ),
  exgaussian_tau = list(
    # Population means
    pr_mu_alpha = matrix(c(0, 1), ncol = 2, nrow = 3, byrow = TRUE),
    pr_mu_beta = matrix(c(0, 1), ncol = 2, nrow = 3, byrow = TRUE),
    pr_shift_sigma = matrix(c(0, 1), ncol = 2, nrow = 3, byrow = TRUE),
    pr_shift_tau = matrix(c(0, 1), ncol = 2, nrow = 3, byrow = TRUE),
    # Population std.devs
    pr_sd_alpha = matrix(c(4, 0, .5), ncol = 3, nrow = 3, byrow = TRUE),
    pr_sd_beta = matrix(c(4, 0, .5), ncol = 3, nrow = 3, byrow = TRUE),
    pr_scale_sigma = matrix(c(4, 0, .5), ncol = 3, nrow = 3, byrow = TRUE),
    pr_scale_tau = matrix(c(4, 0, .5), ncol = 3, nrow = 3, byrow = TRUE),
    # LKJ priors
    pr_L_R_alpha = 1,
    pr_L_Phi = 1,
    # Factor loadings
    pr_L_raw = c(1,1)
  ),
  shlognormal = list(
    # Population means
    pr_mu_alpha = matrix(c(0, 5), ncol = 2, nrow = 3, byrow = TRUE),
    pr_mu_beta = matrix(c(0, 1), ncol = 2, nrow = 3, byrow = TRUE),
    pr_mu_delta = matrix(c(0, 1), ncol = 2, nrow = 3, byrow = TRUE),
    pr_shift_sigma = matrix(c(0, 5), ncol = 2, nrow = 3, byrow = TRUE),
    # Population std.devs
    pr_sd_alpha = matrix(c(4, 0, .5), ncol = 3, nrow = 3, byrow = TRUE),
    pr_sd_beta = matrix(c(4, 0, .5), ncol = 3, nrow = 3, byrow = TRUE),
    pr_sd_delta = matrix(c(4, 0, .5), ncol = 3, nrow = 3, byrow = TRUE),
    pr_scale_sigma = matrix(c(4, 0, .5), ncol = 3, nrow = 3, byrow = TRUE),
    # LKJ priors
    pr_L_R_alpha = 1,
    pr_L_Phi = 1,
    # Factor loadings
    pr_L_raw = c(1,1)
  )
)

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 5: prepare whitehead et al. datasets
# ─────────────────────────────────────────────────────────────────────────────

# Prepare whitehead data frame
load(file = "Data/Processed/whitehead_data.rdata")
whitehead.data$RT <- whitehead.data$RT/1000

# Modify first two column names
colnames(whitehead.data)[1:2] <- c("subject", "condition")

# Generate specific data frames per experiment
whitehead.data.E1 <- whitehead.data[whitehead.data$experiment == 1,
                                    c("subject", "task", "condition", "RT")]
whitehead.data.E2 <- whitehead.data[whitehead.data$experiment == 2,
                                    c("subject", "task", "condition", "RT")]
whitehead.data.E3 <- whitehead.data[whitehead.data$experiment == 3,
                                    c("subject", "task", "condition", "RT")]

# Prepare Stan data: Experiment 1
sdata.E1 <- Stan.list.data(
  data = whitehead.data.E1, 
  subject_var = "subject", 
  task_var = "task", 
  condition_var = "condition",
  congruent_val = 0, 
  incongruent_val = 1, 
  rt_var = "RT")

# Prepare Stan data: Experiment 2
sdata.E2 <- Stan.list.data(
  data = whitehead.data.E2, 
  subject_var = "subject", 
  task_var = "task", 
  condition_var = "condition",
  congruent_val = 0, 
  incongruent_val = 1, 
  rt_var = "RT")

# Prepare Stan data: Experiment 3
sdata.E3 <- Stan.list.data(
  data = whitehead.data.E3, 
  subject_var = "subject", 
  task_var = "task", 
  condition_var = "condition",
  congruent_val = 0, 
  incongruent_val = 1, 
  rt_var = "RT")

# Add latent factors and Lambda index
sdata.E1$M <- sdata.E2$M <- sdata.E3$M <- 1
sdata.E1$Lambda_index <- sdata.E2$Lambda_index <- sdata.E3$Lambda_index <- matrix(1, nrow = 3)

# We're not interested in in LOO-CV right now.
sdata.E1$loo_mix_IS <- sdata.E2$loo_mix_IS <- sdata.E3$loo_mix_IS <- 0
sdata.E1$save_log_lik <- sdata.E2$save_log_lik <- sdata.E3$save_log_lik <- 0

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 6: Gaussian Hierarchical Factor Model
# ─────────────────────────────────────────────────────────────────────────────

# Prepare starting values: Experiment 1
gaussian_inits.E1 <- replicate(
  n = 4, expr = make_inits(
    data           = sdata_to_longdf(sdata.E1), 
    subject_var    = "subject", 
    task_var       = "task", 
    condition_var  = "condition",
    rt_var         = "rt", 
    congruent_id   = 1, 
    incongruent_id = 2, 
    M              = 1, 
    model          = "gaussian"
  ), 
  simplify = FALSE
  )

# Prepare starting values: Experiment 2
gaussian_inits.E2 <- replicate(
  n = 4, expr = make_inits(
    data           = sdata_to_longdf(sdata.E2), 
    subject_var    = "subject", 
    task_var       = "task", 
    condition_var  = "condition",
    rt_var         = "rt", 
    congruent_id   = 1, 
    incongruent_id = 2, 
    M              = 1, 
    model          = "gaussian"), 
  simplify = FALSE
)

# Prepare starting values: Experiment 3
gaussian_inits.E3 <- replicate(
  n = 4, expr = make_inits(
    data           = sdata_to_longdf(sdata.E3), 
    subject_var    = "subject", 
    task_var       = "task", 
    condition_var  = "condition",
    rt_var         = "rt", 
    congruent_id   = 1, 
    incongruent_id = 2, 
    M              = 1, 
    model          = "gaussian"), 
  simplify = FALSE
)

# Fit Gaussian Hierarchical Factor model: Experiment 1
stan_arguments$data <- c(sdata.E1, priors_all$gaussian)
stan_arguments$init <- gaussian_inits.E1
gaussian_HFM_fit_E1 <- do.call(BGHFM_models$gaussian$sample, stan_arguments)

# Fit Gaussian Hierarchical Factor model: Experiment 2
stan_arguments$data <- c(sdata.E2, priors_all$gaussian)
stan_arguments$init <- gaussian_inits.E2
gaussian_HFM_fit_E2 <- do.call(BGHFM_models$gaussian$sample, stan_arguments)

# Fit Gaussian Hierarchical Factor model: Experiment 3
stan_arguments$data <- c(sdata.E3, priors_all$gaussian)
stan_arguments$init <- gaussian_inits.E3
gaussian_HFM_fit_E3 <- do.call(BGHFM_models$gaussian$sample, stan_arguments)

# Save model results
gaussian_HFM_fit_E1$save_object(file = "Results/Stan/Whitehead/Fitted models/whitehead_gaussian_HFM_E1.rds")
gaussian_HFM_fit_E2$save_object(file = "Results/Stan/Whitehead/Fitted models/whitehead_gaussian_HFM_E2.rds")
gaussian_HFM_fit_E3$save_object(file = "Results/Stan/Whitehead/Fitted models/whitehead_gaussian_HFM_E3.rds")

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 6: ex-Gaussian Hierarchical Factor Model
# ─────────────────────────────────────────────────────────────────────────────

# Prepare starting values: Experiment 1
exgaussian_inits.E1 <- replicate(n = 4, expr = make_inits(
  data           = sdata_to_longdf(sdata.E1), 
  subject_var    = "subject", 
  task_var       = "task", 
  condition_var  = "condition",
  rt_var         = "rt", 
  congruent_id   = 1, 
  incongruent_id = 2, 
  M              = 1, 
  model          = "exgaussian"), 
  simplify = FALSE)

# Prepare starting values: Experiment 2
exgaussian_inits.E2 <- replicate(n = 4, expr = make_inits(
  data           = sdata_to_longdf(sdata.E2), 
  subject_var    = "subject", 
  task_var       = "task", 
  condition_var  = "condition",
  rt_var         = "rt", 
  congruent_id   = 1, 
  incongruent_id = 2, 
  M              = 1, 
  model          = "exgaussian"), 
  simplify = FALSE)

# Prepare starting values: Experiment 3
exgaussian_inits.E3 <- replicate(n = 4, expr = make_inits(
  data           = sdata_to_longdf(sdata.E3), 
  subject_var    = "subject", 
  task_var       = "task", 
  condition_var  = "condition",
  rt_var         = "rt", 
  congruent_id   = 1, 
  incongruent_id = 2, 
  M              = 1, 
  model          = "exgaussian"), 
  simplify = FALSE)

# Fit ex-Gaussian Hierarchical Factor model: Experiment 1
stan_arguments$data <- c(sdata.E1, priors_all$exgaussian_tau)
stan_arguments$init <- exgaussian_inits.E1
exgaussian_GHFM_fit_E1 <- do.call(BGHFM_models$exgaussian_tau$sample, stan_arguments)

# Fit ex-Gaussian Hierarchical Factor model: Experiment 2
stan_arguments$data <- c(sdata.E2, priors_all$exgaussian_tau)
stan_arguments$init <- exgaussian_inits.E2
exgaussian_GHFM_fit_E2 <- do.call(BGHFM_models$exgaussian_tau$sample, stan_arguments)

# Fit ex-Gaussian Hierarchical Factor model: Experiment 3
stan_arguments$data <- c(sdata.E3, priors_all$exgaussian_tau)
stan_arguments$init <- exgaussian_inits.E3
exgaussian_GHFM_fit_E3 <- do.call(BGHFM_models$exgaussian_tau$sample, stan_arguments)

# Save model results
exgaussian_GHFM_fit_E1$save_object(file = "Results/Stan/Whitehead/Fitted models/whitehead_exgaussian_GHFM_E1.rds")
exgaussian_GHFM_fit_E2$save_object(file = "Results/Stan/Whitehead/Fitted models/whitehead_exgaussian_GHFM_E2.rds")
exgaussian_GHFM_fit_E3$save_object(file = "Results/Stan/Whitehead/Fitted models/whitehead_exgaussian_GHFM_E3.rds")

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 7: Shifted-lognormal Hierarchical Factor Model
# ─────────────────────────────────────────────────────────────────────────────

# Prepare starting values: Experiment 1
set.seed(2025)
shlognormal_inits.E1 <- replicate(n = 4, expr = make_inits(
  data           = sdata_to_longdf(sdata.E1), 
  subject_var    = "subject", 
  task_var       = "task", 
  condition_var  = "condition",
  rt_var         = "rt", 
  congruent_id   = 1, 
  incongruent_id = 2, 
  M              = 1, 
  model          = "shifted.lognormal"), 
  simplify = FALSE)

# Prepare starting values: Experiment 2
set.seed(2025)
shlognormal_inits.E2 <- replicate(n = 4, expr = make_inits(
  data           = sdata_to_longdf(sdata.E2), 
  subject_var    = "subject", 
  task_var       = "task", 
  condition_var  = "condition",
  rt_var         = "rt", 
  congruent_id   = 1, 
  incongruent_id = 2, 
  M              = 1, 
  model          = "shifted.lognormal"), 
  simplify = FALSE)

# Prepare starting values: Experiment 3
set.seed(2025)
shlognormal_inits.E3 <- replicate(n = 4, expr = make_inits(
  data           = sdata_to_longdf(sdata.E3), 
  subject_var    = "subject", 
  task_var       = "task", 
  condition_var  = "condition",
  rt_var         = "rt", 
  congruent_id   = 1, 
  incongruent_id = 2, 
  M              = 1, 
  model          = "shifted.lognormal"), 
  simplify = FALSE)

# Fit shifted-lognormal Hierarchical Factor model: Experiment 1
stan_arguments$data <- c(sdata.E1, priors_all$shlognormal)
stan_arguments$init <- shlognormal_inits.E1
shlognormal_GHFM_fit_E1 <- do.call(BGHFM_models$shifted_lognormal$sample, stan_arguments)

# Fit shifted-lognormal Hierarchical Factor model: Experiment 2
stan_arguments$data <- c(sdata.E2, priors_all$shlognormal)
stan_arguments$init <- shlognormal_inits.E2
shlognormal_GHFM_fit_E2 <- do.call(BGHFM_models$shifted_lognormal$sample, stan_arguments)

# Fit shifted-lognormal Hierarchical Factor model: Experiment 3
stan_arguments$data <- c(sdata.E3, priors_all$shlognormal)
stan_arguments$init <- shlognormal_inits.E3
shlognormal_GHFM_fit_E3 <- do.call(BGHFM_models$shifted_lognormal$sample, stan_arguments)

# Save model results
shlognormal_GHFM_fit_E1$save_object(file = "Results/Stan/Whitehead/Fitted models/whitehead_shlognormal_GHFM_E1.rds")
shlognormal_GHFM_fit_E2$save_object(file = "Results/Stan/Whitehead/Fitted models/whitehead_shlognormal_GHFM_E2.rds")
shlognormal_GHFM_fit_E3$save_object(file = "Results/Stan/Whitehead/Fitted models/whitehead_shlognormal_GHFM_E3.rds")

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 8: Mixture Importance Sampling in Gaussian HFM
# ─────────────────────────────────────────────────────────────────────────────

# Set LOO-Mixture-IS in all models
sdata.E1$loo_mix_IS <- sdata.E2$loo_mix_IS <- sdata.E3$loo_mix_IS <- 1

# Fit Mixture-IS gaussian Hierarchical Factor Model: Experiment 1
stan_arguments$data <- c(sdata.E1, priors_all$gaussian)
stan_arguments$init <- gaussian_inits.E1
gaussian_HFM_mixture_IS_E1 <- do.call(BGHFM_models$gaussian$sample, stan_arguments)

# Fit Mixture-IS gaussian Hierarchical Factor Model: Experiment 2
stan_arguments$data <- c(sdata.E2, priors_all$gaussian)
stan_arguments$init <- gaussian_inits.E2
gaussian_HFM_mixture_IS_E2 <- do.call(BGHFM_models$gaussian$sample, stan_arguments)

# Fit Mixture-IS gaussian Hierarchical Factor Model: Experiment 3
stan_arguments$data <- c(sdata.E3, priors_all$gaussian)
stan_arguments$init <- gaussian_inits.E3
gaussian_HFM_mixture_IS_E3 <- do.call(BGHFM_models$gaussian$sample, stan_arguments)

# Save Mixture-IS models
gaussian_HFM_mixture_IS_E1$save_object(file = "Results/Stan/Whitehead/Mixture IS models/whitehead_gaussian_HFM_mixtureIS_E1.rds")
gaussian_HFM_mixture_IS_E2$save_object(file = "Results/Stan/Whitehead/Mixture IS models/whitehead_gaussian_HFM_mixtureIS_E2.rds")
gaussian_HFM_mixture_IS_E3$save_object(file = "Results/Stan/Whitehead/Mixture IS models/whitehead_gaussian_HFM_mixtureIS_E3.rds")

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 9: Mixture Importance Sampling in ex-Gaussian HFM
# ─────────────────────────────────────────────────────────────────────────────

# Fit Mixture-IS ex-Gaussian Hierarchical Factor Model: Experiment 1
stan_arguments$data <- c(sdata.E1, priors_all$exgaussian_tau)
stan_arguments$init <- exgaussian_inits.E1
exgaussian_GHFM_mixture_IS_E1 <- do.call(BGHFM_models$exgaussian_tau$sample, stan_arguments)

# Fit Mixture-IS ex-Gaussian Hierarchical Factor Model: Experiment 2
stan_arguments$data <- c(sdata.E2, priors_all$exgaussian_tau)
stan_arguments$init <- exgaussian_inits.E2
exgaussian_GHFM_mixture_IS_E2 <- do.call(BGHFM_models$exgaussian_tau$sample, stan_arguments)

# Fit Mixture-IS ex-Gaussian Hierarchical Factor Model: Experiment 3
stan_arguments$data <- c(sdata.E3, priors_all$exgaussian_tau)
stan_arguments$init <- exgaussian_inits.E3
exgaussian_GHFM_mixture_IS_E3 <- do.call(BGHFM_models$exgaussian_tau$sample, stan_arguments)

# Save Mixture-IS models
exgaussian_GHFM_mixture_IS_E1$save_object(file = "Results/Stan/Whitehead/Mixture IS models/whitehead_exgaussian_GHFM_mixtureIS_E1.rds")
exgaussian_GHFM_mixture_IS_E2$save_object(file = "Results/Stan/Whitehead/Mixture IS models/whitehead_exgaussian_GHFM_mixtureIS_E2.rds")
exgaussian_GHFM_mixture_IS_E3$save_object(file = "Results/Stan/Whitehead/Mixture IS models/whitehead_exgaussian_GHFM_mixtureIS_E3.rds")

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 10: Mixture Importance Sampling in shifted-lognormal HFM
# ─────────────────────────────────────────────────────────────────────────────

# Fit Mixture-IS shifted-lognormal Hierarchical Factor Model: Experiment 1
stan_arguments$data <- c(sdata.E1, priors_all$shlognormal)
stan_arguments$init <- shlognormal_inits.E1
shlognormal_GHFM_mixture_IS_E1 <- do.call(BGHFM_models$shifted_lognormal$sample, stan_arguments)

# Fit Mixture-IS shifted-lognormal Hierarchical Factor Model: Experiment 2
stan_arguments$data <- c(sdata.E2, priors_all$shlognormal)
stan_arguments$init <- shlognormal_inits.E2
shlognormal_GHFM_mixture_IS_E2 <- do.call(BGHFM_models$shifted_lognormal$sample, stan_arguments)

# Fit Mixture-IS shifted-lognormal Hierarchical Factor Model: Experiment 3
stan_arguments$data <- c(sdata.E3, priors_all$shlognormal)
stan_arguments$init <- shlognormal_inits.E3
shlognormal_GHFM_mixture_IS_E3 <- do.call(BGHFM_models$shifted_lognormal$sample, stan_arguments)

# Save all Mixture-IS models
shlognormal_GHFM_mixture_IS_E1$save_object(file = "Results/Stan/Whitehead/Mixture IS models/whitehead_shlognormal_GHFM_mixtureIS_E1.rds")
shlognormal_GHFM_mixture_IS_E2$save_object(file = "Results/Stan/Whitehead/Mixture IS models/whitehead_shlognormal_GHFM_mixtureIS_E2.rds")
shlognormal_GHFM_mixture_IS_E3$save_object(file = "Results/Stan/Whitehead/Mixture IS models/whitehead_shlognormal_GHFM_mixtureIS_E3.rds")

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 11: Estimate ELPDs using Mixture-IS 
# ─────────────────────────────────────────────────────────────────────────────

# Estimate Gaussian HFM Mixture-IS ELPDs: Experiment 1
loo_mixis_gaussian_E1 <- loo_mixis_BGHFM(
  fit = gaussian_HFM_mixture_IS_E1, 
  model = "gaussian", 
  data = sdata_to_longdf(sdata.E1), 
  subject_var = "subject", 
  task_var = "task", 
  condition_var = "condition", 
  rt_var = "rt", 
  condition_center = 1.5, 
  progress_bar = TRUE)

# Estimate Gaussian HFM Mixture-IS ELPDs: Experiment 2
loo_mixis_gaussian_E2 <- loo_mixis_BGHFM(
  fit = gaussian_HFM_mixture_IS_E2, 
  model = "gaussian", 
  data = sdata_to_longdf(sdata.E2), 
  subject_var = "subject", 
  task_var = "task", 
  condition_var = "condition", 
  rt_var = "rt", 
  condition_center = 1.5, 
  progress_bar = TRUE)

# Estimate Gaussian HFM Mixture-IS ELPDs: Experiment 3
loo_mixis_gaussian_E3 <- loo_mixis_BGHFM(
  fit = gaussian_HFM_mixture_IS_E3, 
  model = "gaussian", 
  data = sdata_to_longdf(sdata.E3), 
  subject_var = "subject", 
  task_var = "task", 
  condition_var = "condition", 
  rt_var = "rt", 
  condition_center = 1.5, 
  progress_bar = TRUE)

# Estimate ex-Gaussian HFM Mixture-IS ELPDs: Experiment 1
loo_mixis_exgaussian_E1 <- loo_mixis_BGHFM(
  fit = exgaussian_GHFM_mixture_IS_E1, 
  model = "exgaussian", 
  data = sdata_to_longdf(sdata.E1), 
  subject_var = "subject", 
  task_var = "task", 
  condition_var = "condition", 
  rt_var = "rt", 
  condition_center = 1.5, 
  progress_bar = TRUE)

# Estimate ex-Gaussian HFM Mixture-IS ELPDs: Experiment 2
loo_mixis_exgaussian_E2 <- loo_mixis_BGHFM(
  fit = exgaussian_GHFM_mixture_IS_E2, 
  model = "exgaussian", 
  data = sdata_to_longdf(sdata.E2), 
  subject_var = "subject", 
  task_var = "task", 
  condition_var = "condition", 
  rt_var = "rt", 
  condition_center = 1.5, 
  progress_bar = TRUE)

# Estimate ex-Gaussian HFM Mixture-IS ELPDs: Experiment 3
loo_mixis_exgaussian_E3 <- loo_mixis_BGHFM(
  fit = exgaussian_GHFM_mixture_IS_E3, 
  model = "exgaussian", 
  data = sdata_to_longdf(sdata.E3), 
  subject_var = "subject", 
  task_var = "task", 
  condition_var = "condition", 
  rt_var = "rt", 
  condition_center = 1.5, 
  progress_bar = TRUE)

# Estimate shifted-lognormal HFM Mixture-IS ELPDs: Experiment 1
loo_mixis_shlognormal_E1 <- loo_mixis_BGHFM(
  fit = shlognormal_GHFM_mixture_IS_E1, 
  model = "shifted.lognormal", 
  data = sdata_to_longdf(sdata.E1), 
  subject_var = "subject", 
  task_var = "task", 
  condition_var = "condition", 
  rt_var = "rt", 
  condition_center = 1.5, 
  progress_bar = TRUE)

# Estimate shifted-lognormal HFM Mixture-IS ELPDs: Experiment 2
loo_mixis_shlognormal_E2 <- loo_mixis_BGHFM(
  fit = shlognormal_GHFM_mixture_IS_E2, 
  model = "shifted.lognormal", 
  data = sdata_to_longdf(sdata.E2), 
  subject_var = "subject", 
  task_var = "task", 
  condition_var = "condition", 
  rt_var = "rt", 
  condition_center = 1.5, 
  progress_bar = TRUE)

# Estimate shifted-lognormal HFM Mixture-IS ELPDs: Experiment 3
loo_mixis_shlognormal_E3 <- loo_mixis_BGHFM(
  fit = shlognormal_GHFM_mixture_IS_E3, 
  model = "shifted.lognormal", 
  data = sdata_to_longdf(sdata.E3), 
  subject_var = "subject", 
  task_var = "task", 
  condition_var = "condition", 
  rt_var = "rt", 
  condition_center = 1.5, 
  progress_bar = TRUE)

# Save all ELPD estimates
whitehead_experiments_ELPDs <- list(
  ELPD_gaussian_E1    = loo_mixis_gaussian_E1,
  ELPD_gaussian_E2    = loo_mixis_gaussian_E2,
  ELPD_gaussian_E3    = loo_mixis_gaussian_E3,
  ELPD_exgaussian_E1  = loo_mixis_exgaussian_E1,
  ELPD_exgaussian_E2  = loo_mixis_exgaussian_E2,
  ELPD_exgaussian_E3  = loo_mixis_exgaussian_E3,
  ELPD_shlognormal_E1 = loo_mixis_shlognormal_E1,
  ELPD_shlognormal_E2 = loo_mixis_shlognormal_E2,
  ELPD_shlognormal_E3 = loo_mixis_shlognormal_E3
)

# Save mixture ELPDs 
saveRDS(whitehead_experiments_ELPDs, file = "Results/Rdata/ELPDs/whitehead_experiments_ELPDs.rds")

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 9: Interpret model factor loadings
# ─────────────────────────────────────────────────────────────────────────────

# Posterior distribution: factor loadings
gaussian_GHFM_fit_E1$summary("Lambda_std", mean, ~quantile(.x, probs = c(.025, .975)))
exgaussian_GHFM_fit_E1$summary("Lambda_std", mean, ~quantile(.x, probs = c(.025, .975)))
shlognorm_GHFM_fit_E1$summary("Lambda_std", mean, ~quantile(.x, probs = c(.025, .975)))
gaussian_GHFM_fit_E2$summary("Lambda_std", mean, ~quantile(.x, probs = c(.025, .975)))
exgaussian_GHFM_fit_E2$summary("Lambda_std", mean, ~quantile(.x, probs = c(.025, .975)))
shlognorm_GHFM_fit_E2$summary("Lambda_std", mean, ~quantile(.x, probs = c(.025, .975)))
gaussian_GHFM_fit_E3$summary("Lambda_std", mean, ~quantile(.x, probs = c(.025, .975)))
exgaussian_GHFM_fit_E3$summary("Lambda_std", mean, ~quantile(.x, probs = c(.025, .975)))
shlognorm_GHFM_fit_E3$summary("Lambda_std", mean, ~quantile(.x, probs = c(.025, .975)))

# ─────────────────────────────────────────────────────────────────────────────

