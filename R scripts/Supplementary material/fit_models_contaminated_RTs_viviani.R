# ╔════════════════════════════════════════════════════════════════════════════╗
# ║                             SCRIPT OVERVIEW                                ║
# ╠════════════════════════════════════════════════════════════════════════════╣
# ║ Script Name   : fit_models_contaminated_RTs_viviani.R                      ║
# ║ Author        : Ricardo Rey-Sáez                                           ║
# ║ Role          : PhD Student in Psychology                                  ║
# ║ Institution   : Autonomous University of Madrid, Spain                     ║
# ║ Email         : ricardoreysaez95@gmail.es                                  ║
# ║ Date          : 17-05-2025                                                 ║
# ╚════════════════════════════════════════════════════════════════════════════╝

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 1: Description
# ─────────────────────────────────────────────────────────────────────────────
# Fit mixture models to account for contaminated RTs in Viviani dataset
# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 2: Load Packages
# ─────────────────────────────────────────────────────────────────────────────

library(cmdstanr)
library(posterior)
library(openxlsx) 
library(dplyr)    
library(tidyr)    
library(lme4)

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 3: Global settings
# ─────────────────────────────────────────────────────────────────────────────

# Extra-functions to analyze data
source("R scripts/R functions/empirical_functions.R")

# Load the mixture-shifted-lognormal contamination RT model
# Based on Haines et al., (2025)
# https://github.com/Nathaniel-Haines/Reliability_2020/blob/master/Code/Stan/joint_RT_shiftlnorm_mix.stan
BGHFM_models <- list(
  shifted_lognormal_mix = cmdstan_model("Stan models/BGHFM/E_BGHFM_shlognormal_mix.stan")
)

# Global Stan arguments
stan_arguments <- list(
  data = NULL, 
  iter_warmup = 1500,
  iter_sampling = 2000,
  chains = 4, 
  parallel_chains = 4,
  adapt_delta = 0.99,
  max_treedepth = 15,
  seed = 2025,
  refresh = 50,
  init = NULL)

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 4: Set prior parameters for the shifted-lognormal model
# ─────────────────────────────────────────────────────────────────────────────

# Meta-parameters
load("Results/Rdata/Meta-parameters/meta_parameters_full.rdata")

# Model names
name_map <- list(
  shlognormal      = "meta_parameters_full_shlognormal"
)

# Parameters in each model
mappings <- list(
  shlognormal = list(
    means = list(
      pr_mu_alpha    = "mu_alpha",
      pr_mu_beta     = "mu_theta",
      pr_mu_delta    = "mu_delta",
      pr_shift_sigma = "shift_sigma"
    ),
    sds   = list(
      pr_sd_alpha    = "sd_alpha",
      pr_sd_beta     = "sd_theta",
      pr_sd_delta    = "sd_delta",
      pr_scale_sigma = "scale_sigma"
    )
  )
)

# Helper function to compute SEs
parse_hw <- function(ci_str) {
  nums <- as.numeric(strsplit(gsub("\\[|\\]", "", ci_str), ",")[[1]])
  (nums[2] - nums[1]) / 2
}
J  <- 6L  
nu <- 3L # Student t degrees of freedom for standard deviation priors

# Set prior distribution parameters per model:
priors_all <- setNames(
  lapply(names(mappings), function(model_key) {
    # Results from the stroop meta-analysis
    df_stroop <- subset(
      meta_parameters_full[[ name_map[[model_key]] ]],
      task == "stroop"
    )
    
    map_spec <- mappings[[model_key]]
    model_list <- list()
    
    # Mapping population means
    for (stan_name in names(map_spec$means)) {
      col_data <- map_spec$means[[stan_name]]
      ci_col   <- paste0("CI_", col_data)
      v        <- c(df_stroop[[ col_data ]], parse_hw(df_stroop[[ ci_col ]]))
      model_list[[stan_name]] <- matrix(rep(v, J), nrow = J, byrow = TRUE)
    }
    
    # Mapping population standard deviations
    for (stan_name in names(map_spec$sds)) {
      col_data <- map_spec$sds[[stan_name]]
      ci_col   <- paste0("CI_", col_data)
      v        <- c(nu, df_stroop[[ col_data ]], parse_hw(df_stroop[[ ci_col ]]))
      model_list[[stan_name]] <- matrix(rep(v, J), nrow = J, byrow = TRUE)
    }
    
    # Exploratory Factor Analysis parameters
    model_list$pr_L_R_alpha <- 1
    model_list$pr_h2        <- c(1, 1)
    
    # Add contamination mixture prior parameters
    # Contamination process
    model_list$pr_mu_pi = matrix(c(-1.65, .1), ncol = 2, nrow = J, byrow = TRUE)
    model_list$pr_sd_pi = matrix(c(3, 0, .2), ncol = 3, nrow = J, byrow = TRUE)
    
    model_list
  }),
  names(mappings)
)

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 5: Prepare Viviani et al. Datasets
# ─────────────────────────────────────────────────────────────────────────────

# Load long data frame with level 1 observed response times
load(file = "Data/Processed/viviani_data.rdata")

# Condition ID
# tapply(vivianni_2024_df$RT, vivianni_2024_df$condition, mean)

# Prepare Stan data
sdata <- Stan.list.data(
  data = vivianni_2024_df, 
  subject_var = "subject", 
  task_var = "task", 
  condition_var = "condition",
  congruent_val = 1, 
  incongruent_val = 0,
  rt_var = "RT")

# Add latent factos and repulsion parameter
sdata$M <- 2
sdata$xi <- 100

# We're not interested in in LOO-CV right now.
sdata$loo_mix_IS <- 0
sdata$save_log_lik <- 0

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 6: Shifted-lognormal Mixture-RTs Hierarchical Factor Model
# ─────────────────────────────────────────────────────────────────────────────

# Prepare initial values
set.seed(2025)
shlognormal_mix_inits <- replicate(n = 4, expr = make_inits(
  data           = sdata_to_longdf(sdata), 
  subject_var    = "subject", 
  task_var       = "task", 
  condition_var  = "condition",
  rt_var         = "rt", 
  congruent_id   = 1, 
  incongruent_id = 2, 
  M              = 2, 
  model          = "shifted.lognormal.mix"), 
  simplify = FALSE)

# Fit shifted-lognormal Hierarchical Factor model
stan_arguments$init <- shlognormal_mix_inits
stan_arguments$data <- c(sdata, priors_all$shlognormal)
shlognormal_mix_GHFM_fit <- do.call(BGHFM_models$shifted_lognormal_mix$sample, stan_arguments)

# Save model results
shlognormal_mix_GHFM_fit$save_object(
  file = "Results/Stan/Viviani/Fitted models/viviani_shlognormal_mix_GHFM.rds")

# ─────────────────────────────────────────────────────────────────────────────
