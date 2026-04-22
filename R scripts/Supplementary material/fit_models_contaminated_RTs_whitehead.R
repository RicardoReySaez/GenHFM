# ╔════════════════════════════════════════════════════════════════════════════╗
# ║                             SCRIPT OVERVIEW                                ║
# ╠════════════════════════════════════════════════════════════════════════════╣
# ║ Script Name   : fit_models_contaminated_RTs_whitehead.R                    ║
# ║ Author        : Ricardo Rey-Sáez                                           ║
# ║ Role          : PhD Student in Psychology                                  ║
# ║ Institution   : Autonomous University of Madrid, Spain                     ║
# ║ Email         : ricardoreysaez95@gmail.es                                  ║
# ║ Date          : 17-05-2025                                                 ║
# ╚════════════════════════════════════════════════════════════════════════════╝

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 1: Description
# ─────────────────────────────────────────────────────────────────────────────
# Fit Shifted-lognormal mixture GHFMs
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
  shlognormal_mix = cmdstan_model("Stan models/BGHFM/C_BGHFM_shlognormal_mix.stan")
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
    pr_L_raw = c(1,1),
    # Contamination process
    pr_mu_pi = matrix(c(-1.65, .1), ncol = 2, nrow = 3, byrow = TRUE),
    pr_sd_pi = matrix(c(4, 0, .5), ncol = 3, nrow = 3, byrow = TRUE)
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
# SECTION 7: Shifted-Lognormal Mixture-RT Hierarchical Factor Model
# ─────────────────────────────────────────────────────────────────────────────

# Prepare starting values: Experiment 1
set.seed(2025)
shlognormal_mix_inits.E1 <- replicate(n = 4, expr = make_inits(
  data           = sdata_to_longdf(sdata.E1), 
  subject_var    = "subject", 
  task_var       = "task", 
  condition_var  = "condition",
  rt_var         = "rt", 
  congruent_id   = 1, 
  incongruent_id = 2, 
  M              = 1, 
  model          = "shifted.lognormal.mix"), 
  simplify = FALSE)

# Prepare starting values: Experiment 2
set.seed(2025)
shlognormal_mix_inits.E2 <- replicate(n = 4, expr = make_inits(
  data           = sdata_to_longdf(sdata.E2), 
  subject_var    = "subject", 
  task_var       = "task", 
  condition_var  = "condition",
  rt_var         = "rt", 
  congruent_id   = 1, 
  incongruent_id = 2, 
  M              = 1, 
  model          = "shifted.lognormal.mix"), 
  simplify = FALSE)

# Prepare starting values: Experiment 3
set.seed(2025)
shlognormal_mix_inits.E3 <- replicate(n = 4, expr = make_inits(
  data           = sdata_to_longdf(sdata.E3), 
  subject_var    = "subject", 
  task_var       = "task", 
  condition_var  = "condition",
  rt_var         = "rt", 
  congruent_id   = 1, 
  incongruent_id = 2, 
  M              = 1, 
  model          = "shifted.lognormal.mix"), 
  simplify = FALSE)

# Fit shifted-lognormal Mixture-RT Hierarchical Factor model: Experiment 1
stan_arguments$data <- c(sdata.E1, priors_all$shlognormal)
stan_arguments$init <- shlognormal_mix_inits.E1
shlognormal_mix_GHFM_fit_E1 <- do.call(BGHFM_models$shlognormal_mix$sample, stan_arguments)

# Fit shifted-lognormal Mixture-RT Hierarchical Factor model: Experiment 2
stan_arguments$data <- c(sdata.E2, priors_all$shlognormal)
stan_arguments$init <- shlognormal_mix_inits.E2
shlognormal_mix_GHFM_fit_E2 <- do.call(BGHFM_models$shlognormal_mix$sample, stan_arguments)

# Fit shifted-lognormal Mixture-RT Hierarchical Factor model: Experiment 3
stan_arguments$data <- c(sdata.E3, priors_all$shlognormal)
stan_arguments$init <- shlognormal_mix_inits.E3
shlognormal_mix_GHFM_fit_E3 <- do.call(BGHFM_models$shlognormal_mix$sample, stan_arguments)

# Save model results
shlognormal_mix_GHFM_fit_E1$save_object(
  file = "Results/Stan/Whitehead/Fitted models/whitehead_shlognormal_mix_GHFM_E1.rds")
shlognormal_mix_GHFM_fit_E2$save_object(
  file = "Results/Stan/Whitehead/Fitted models/whitehead_shlognormal_mix_GHFM_E2.rds")
shlognormal_mix_GHFM_fit_E3$save_object(
  file = "Results/Stan/Whitehead/Fitted models/whitehead_shlognormal_mix_GHFM_E3.rds")

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 8: Joint Lognormal Mixture Model (all experiments pooled)
# ─────────────────────────────────────────────────────────────────────────────

# Prepare Stan data: all experiments jointly
sdata.joint <- Stan.list.data(
  data = whitehead.data, 
  subject_var = "subject", 
  task_var = "task", 
  condition_var = "condition",
  congruent_val = 0, 
  incongruent_val = 1, 
  rt_var = "RT")

# Add latent factors and Lambda index
sdata.joint$M <- 1
sdata.joint$Lambda_index <- matrix(1, nrow = 3)

# We're not interested in LOO-CV right now
sdata.joint$loo_mix_IS <- 0
sdata.joint$save_log_lik <- 0

# Prepare starting values: joint model
set.seed(2025)
shlognormal_mix_inits.joint <- replicate(n = 4, expr = make_inits(
  data           = sdata_to_longdf(sdata.joint), 
  subject_var    = "subject", 
  task_var       = "task", 
  condition_var  = "condition",
  rt_var         = "rt", 
  congruent_id   = 1, 
  incongruent_id = 2, 
  M              = 1, 
  model          = "shifted.lognormal.mix"), 
  simplify = FALSE)

# Fit shifted-lognormal Mixture Hierarchical Factor model: joint
stan_arguments$data <- c(sdata.joint, priors_all$shlognormal)
stan_arguments$init <- shlognormal_mix_inits.joint
shlognormal_mix_GHFM_fit_joint <- do.call(BGHFM_models$shlognormal_mix$sample, stan_arguments)

# Save model results
shlognormal_mix_GHFM_fit_joint$save_object(
  file = "Results/Stan/Whitehead/Fitted models/whitehead_shlognormal_mix_GHFM_joint.rds")

# ─────────────────────────────────────────────────────────────────────────────
