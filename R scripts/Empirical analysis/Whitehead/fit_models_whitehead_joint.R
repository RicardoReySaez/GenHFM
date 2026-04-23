# ╔════════════════════════════════════════════════════════════════════════════╗
# ║                             SCRIPT OVERVIEW                                ║
# ╠════════════════════════════════════════════════════════════════════════════╣
# ║ Script Name   : fit_models_whitehead_joint.R                               ║
# ║ Author        : Ricardo Rey-Sáez                                           ║
# ║ Role          : PhD Student in Psychology                                  ║
# ║ Institution   : Autonomous University of Madrid, Spain                     ║
# ║ Email         : ricardoreysaez95@gmail.es                                  ║
# ║ Date          : 17-05-2025                                                 ║
# ╚════════════════════════════════════════════════════════════════════════════╝

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 1: Description
# ─────────────────────────────────────────────────────────────────────────────
# Fit GHFMs and Mixture-IS models for all experiments toguether
# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 2: Load Packages
# ─────────────────────────────────────────────────────────────────────────────

library(cmdstanr)
library(posterior)
library(bayesplot)

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
  lognormal         = cmdstan_model("Stan models/BGHFM/C_BGHFM_lognormal.stan"),
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
# SECTION 5: Prepare Whitehead et al. Datasets
# ─────────────────────────────────────────────────────────────────────────────

# Prepare whitehead data frame
load(file = "Data/Processed/whitehead_data.rdata")
whitehead.data$RT <- whitehead.data$RT/1000

# Modify first two column names
colnames(whitehead.data)[1:2] <- c("subject", "condition")

# Prepare Stan data
sdata <- Stan.list.data(
  data = whitehead.data, 
  subject_var = "subject", 
  task_var = "task", 
  condition_var = "condition",
  congruent_val = 0, 
  incongruent_val = 1, 
  rt_var = "RT")

# Add latent factors and Lambda index
sdata$M <-  1
sdata$Lambda_index <- matrix(1, nrow = 3)

# We're not interested in in LOO-CV right now.
sdata$loo_mix_IS <- 0
sdata$save_log_lik <- 0

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 6: Gaussian Hierarchical Factor Model
# ─────────────────────────────────────────────────────────────────────────────

# Prepare starting values
gaussian_inits <- replicate(
  n = 4, expr = make_inits(
    data           = sdata_to_longdf(sdata), 
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

# Fit Gaussian Hierarchical Factor model
stan_arguments$data <- c(sdata, priors_all$gaussian)
stan_arguments$init <- gaussian_inits
gaussian_HFM_fit <- do.call(BGHFM_models$gaussian$sample, stan_arguments)

# Save model results
gaussian_HFM_fit$save_object(file = "Results/Stan/Whitehead/Fitted models/whitehead_gaussian_HFM_joint.rds")

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 7: ex-Gaussian Hierarchical Factor Model
# ─────────────────────────────────────────────────────────────────────────────

# Prepare starting values
exgaussian_inits <- replicate(n = 4, expr = make_inits(
  data           = sdata_to_longdf(sdata), 
  subject_var    = "subject", 
  task_var       = "task", 
  condition_var  = "condition",
  rt_var         = "rt", 
  congruent_id   = 1, 
  incongruent_id = 2, 
  M              = 1, 
  model          = "exgaussian"), 
  simplify = FALSE)

# Fit ex-Gaussian Hierarchical Factor model
stan_arguments$data <- c(sdata, priors_all$exgaussian_tau)
stan_arguments$init <- exgaussian_inits
exgaussian_GHFM_fit <- do.call(BGHFM_models$exgaussian_tau$sample, stan_arguments)

# Save model results
exgaussian_GHFM_fit$save_object(file = "Results/Stan/Whitehead/Fitted models/whitehead_exgaussian_GHFM_joint.rds")
# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 8: Lognormal Hierarchical Factor Model
# ─────────────────────────────────────────────────────────────────────────────

# Prepare starting values
set.seed(2025)
lognormal_inits <- replicate(n = 4, expr = make_inits(
  data           = sdata_to_longdf(sdata), 
  subject_var    = "subject", 
  task_var       = "task", 
  condition_var  = "condition",
  rt_var         = "rt", 
  congruent_id   = 1, 
  incongruent_id = 2, 
  M              = 1, 
  model          = "lognormal"), 
  simplify = FALSE)

# Fit lognormal Hierarchical Factor model
stan_arguments$data <- c(sdata, priors_all$shlognormal)
stan_arguments$init <- lognormal_inits
lognormal_GHFM_fit <- do.call(BGHFM_models$lognormal$sample, stan_arguments)

# Save model results
lognormal_GHFM_fit$save_object(file = "Results/Stan/Whitehead/Fitted models/whitehead_lognormal_GHFM_joint.rds")

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 9: Shifted-lognormal Hierarchical Factor Model
# ─────────────────────────────────────────────────────────────────────────────

# Prepare starting values
set.seed(2025)
shlognormal_inits <- replicate(n = 4, expr = make_inits(
  data           = sdata_to_longdf(sdata), 
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
stan_arguments$data <- c(sdata, priors_all$shlognormal)
stan_arguments$init <- shlognormal_inits
shlognormal_GHFM_fit <- do.call(BGHFM_models$shifted_lognormal$sample, stan_arguments)

# Save model results
shlognormal_GHFM_fit$save_object(file = "Results/Stan/Whitehead/Fitted models/whitehead_shlognormal_GHFM_joint.rds")

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 10: Mixture Importance Sampling in all GenHFMs
# ─────────────────────────────────────────────────────────────────────────────

# Set LOO-Mixture-IS in all models
sdata$loo_mix_IS <- 1

# Fit Mixture-IS gaussian Hierarchical Factor Model
stan_arguments$data <- c(sdata, priors_all$gaussian)
stan_arguments$init <- gaussian_inits
gaussian_HFM_mixture_IS <- do.call(BGHFM_models$gaussian$sample, stan_arguments)
gaussian_HFM_mixture_IS$save_object(file = "Results/Stan/Whitehead/Fitted models/whitehead_gaussian_HFM_mixtureIS_joint.rds")

# Fit Mixture-IS ex-gaussian Hierarchical Factor Model
stan_arguments$data <- c(sdata, priors_all$exgaussian_tau)
stan_arguments$init <- exgaussian_inits
exgaussian_GHFM_mixture_IS <- do.call(BGHFM_models$exgaussian_tau$sample, stan_arguments)
exgaussian_GHFM_mixture_IS$save_object(file = "Results/Stan/Whitehead/Fitted models/whitehead_exgaussian_GHFM_mixtureIS_joint.rds")

# Fit Mixture-IS lognormal Hierarchical Factor Model
stan_arguments$data <- c(sdata, priors_all$shlognormal)
stan_arguments$init <- lognormal_inits
lognormal_GHFM_mixture_IS <- do.call(BGHFM_models$lognormal$sample, stan_arguments)
lognormal_GHFM_mixture_IS$save_object(file = "Results/Stan/Whitehead/Fitted models/whitehead_lognormal_GHFM_mixtureIS_joint.rds")

# Fit Mixture-IS ex-gaussian Hierarchical Factor Model
stan_arguments$data <- c(sdata, priors_all$shlognormal)
stan_arguments$init <- shlognormal_inits
shlognormal_GHFM_mixture_IS <- do.call(BGHFM_models$shifted_lognormal$sample, stan_arguments)
shlognormal_GHFM_mixture_IS$save_object(file = "Results/Stan/Whitehead/Fitted models/whitehead_shlognormal_GHFM_mixtureIS_joint.rds")

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 10: Estimate ELPDs using Mixture-IS 
# ─────────────────────────────────────────────────────────────────────────────

# Estimate Gaussian HFM Mixture-IS ELPDs
loo_mixis_gaussian <- loo_mixis_BGHFM(
  fit = gaussian_HFM_mixture_IS, 
  model = "gaussian", 
  data = sdata_to_longdf(sdata), 
  subject_var = "subject", 
  task_var = "task", 
  condition_var = "condition", 
  rt_var = "rt", 
  condition_center = 1.5, 
  progress_bar = TRUE)

# Estimate ex-Gaussian HFM Mixture-IS ELPDs
loo_mixis_exgaussian <- loo_mixis_BGHFM(
  fit = exgaussian_GHFM_mixture_IS, 
  model = "exgaussian", 
  data = sdata_to_longdf(sdata), 
  subject_var = "subject", 
  task_var = "task", 
  condition_var = "condition", 
  rt_var = "rt", 
  condition_center = 1.5, 
  progress_bar = TRUE)

# Estimate lognormal HFM Mixture-IS ELPDs
loo_mixis_lognormal <- loo_mixis_BGHFM(
  fit = lognormal_GHFM_mixture_IS, 
  model = "lognormal", 
  data = sdata_to_longdf(sdata), 
  subject_var = "subject", 
  task_var = "task", 
  condition_var = "condition", 
  rt_var = "rt", 
  condition_center = 1.5, 
  progress_bar = TRUE)

# Estimate shifted-lognormal HFM Mixture-IS ELPDs: Experiment 3
loo_mixis_shlognormal <- loo_mixis_BGHFM(
  fit = shlognormal_GHFM_mixture_IS, 
  model = "shifted.lognormal", 
  data = sdata_to_longdf(sdata), 
  subject_var = "subject", 
  task_var = "task", 
  condition_var = "condition", 
  rt_var = "rt", 
  condition_center = 1.5, 
  progress_bar = TRUE)

# Save all ELPD estimates
whitehead_joint_ELPDs <- list(
  ELPD_gaussian_joint    = loo_mixis_gaussian,
  ELPD_exgaussian_joint  = loo_mixis_exgaussian,
  ELPD_lognormal_joint   = loo_mixis_lognormal,
  ELPD_shlognormal_joint = loo_mixis_shlognormal
)

# Save mixture ELPDs 
saveRDS(whitehead_joint_ELPDs, file = "Results/Rdata/ELPDs/whitehead_joint_ELPDs.rds")

# ─────────────────────────────────────────────────────────────────────────────
