# ╔════════════════════════════════════════════════════════════════════════════╗
# ║                             SCRIPT OVERVIEW                                ║
# ╠════════════════════════════════════════════════════════════════════════════╣
# ║ Script Name   : analyse_vivianni.R                                         ║
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
#   1. Fit gaussian, ex-Gaussian and shifted-lognormal GHFMs
#   2. Estimate Mixture-IS models
#   3. Compute ELPDs per model
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

# Save all BGHFM in one list
BGHFM_models <- list(
  gaussian          = cmdstan_model("Stan models/BGHFM/E_BGHFM_gaussian.stan"),
  exgaussian_rate   = cmdstan_model("Stan models/BGHFM/E_BGHFM_exgaussian_rate.stan"),
  shifted_lognormal = cmdstan_model("Stan models/BGHFM/E_BGHFM_shlognormal.stan")
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
# SECTION 4: Set prior parameters for all models (informative priors)
# ─────────────────────────────────────────────────────────────────────────────

# Meta-parameters
load("Results/Rdata/meta_parameters_full.rdata")

# Model names
name_map <- list(
  gaussian         = "meta_parameters_full_gauss",
  exgaussian_rate  = "meta_parameters_full_exgauss_rate",
  shlognormal      = "meta_parameters_full_shlognormal"
)

# Parameters in each model
mappings <- list(
  gaussian = list(
    means = list(
      pr_mu_alpha    = "mu_alpha",
      pr_mu_beta     = "mu_theta",
      pr_shift_sigma = "shift_sigma"
    ),
    sds   = list(
      pr_sd_alpha    = "sd_alpha",
      pr_sd_beta    = "sd_theta",
      pr_scale_sigma = "scale_sigma"
    )
  ),
  exgaussian_rate = list(
    means = list(
      pr_mu_alpha    = "mu_alpha",
      pr_mu_beta    = "mu_theta",
      pr_shift_sigma = "shift_sigma",
      pr_shift_rate  = "shift_lambda" 
    ),
    sds   = list(
      pr_sd_alpha    = "sd_alpha",
      pr_sd_beta    = "sd_theta",
      pr_scale_sigma = "scale_sigma",
      pr_scale_rate  = "scale_lambda" 
    )
  ),
  shlognormal = list(
    means = list(
      pr_mu_alpha    = "mu_alpha",
      pr_mu_beta    = "mu_theta",
      pr_mu_delta    = "mu_delta",
      pr_shift_sigma = "shift_sigma"
    ),
    sds   = list(
      pr_sd_alpha    = "sd_alpha",
      pr_sd_beta    = "sd_theta",
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
    
    model_list
  }),
  names(mappings)
)

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 5: prepare Viviani et al. datasets
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
# SECTION 6: Gaussian Hierarchical Factor Model
# ─────────────────────────────────────────────────────────────────────────────

# Prepare initial values
gaussian_inits <- replicate(n = 4, expr = make_inits(
  data           = sdata_to_longdf(sdata), 
  subject_var    = "subject", 
  task_var       = "task", 
  condition_var  = "condition",
  rt_var         = "rt", 
  congruent_id   = 1, 
  incongruent_id = 2, 
  M              = 2, 
  model          = "gaussian"), 
  simplify = FALSE)

# Fit Gaussian Hierarchical Factor model
stan_arguments$data <- c(sdata, priors_all$gaussian)
stan_arguments$init <- gaussian_inits
gaussian_HFM_fit <- do.call(BGHFM_models$gaussian$sample, stan_arguments)

# Save model results
gaussian_HFM_fit$save_object(file = "Results/Stan/Viviani/Fitted models/viviani_gaussian_HFM.rds")

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 7: ex-Gaussian Hierarchical Factor Model
# ─────────────────────────────────────────────────────────────────────────────

# Prepare initial values
exgaussian_inits <- replicate(n = 4, expr = make_inits(
  data           = sdata_to_longdf(sdata), 
  subject_var    = "subject", 
  task_var       = "task", 
  condition_var  = "condition",
  rt_var         = "rt", 
  congruent_id   = 1, 
  incongruent_id = 2, 
  M              = 2, 
  model          = "exgaussian"), 
  simplify = FALSE)

# Fit ex-Gaussian Hierarchical Factor model
stan_arguments$init <- exgaussian_inits
stan_arguments$data <- c(sdata, priors_all$exgaussian_rate)
exgaussian_GHFM_fit <- do.call(BGHFM_models$exgaussian_rate$sample, stan_arguments)

# Save model results
exgaussian_GHFM_fit$save_object(file = "Results/Stan/Viviani/Fitted models/viviani_exgaussian_GHFM.rds")

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 8: Shifted-lognormal Hierarchical Factor Model
# ─────────────────────────────────────────────────────────────────────────────

# Prepare initial values
set.seed(2026)
shlognormal_inits <- replicate(n = 4, expr = make_inits(
  data           = sdata_to_longdf(sdata), 
  subject_var    = "subject", 
  task_var       = "task", 
  condition_var  = "condition",
  rt_var         = "rt", 
  congruent_id   = 1, 
  incongruent_id = 2, 
  M              = 2, 
  model          = "shifted.lognormal"), 
  simplify = FALSE)

# Fit shifted-lognormal Hierarchical Factor model
stan_arguments$init <- shlognormal_inits
stan_arguments$data <- c(sdata, priors_all$shlognormal)
shlognormal_GHFM_fit <- do.call(BGHFM_models$shifted_lognormal$sample, stan_arguments)

# Save model results
shlognormal_GHFM_fit$save_object(file = "Results/Stan/Viviani/Fitted models/viviani_shlognormal_GHFM.rds")

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 9: Mixture Importance Sampling in all models
# ─────────────────────────────────────────────────────────────────────────────

# Set LOO-Mixture-IS in all models
sdata$loo_mix_IS <- 1

# Fit Mixture-IS gaussian Hierarchical Factor Model
stan_arguments$data <- c(sdata, priors_all$gaussian)
stan_arguments$init <- gaussian_inits
gaussian_HFM_mixture_IS <- do.call(BGHFM_models$gaussian$sample, stan_arguments)

# Fit Mixture-IS ex-Gaussian Hierarchical Factor Model
stan_arguments$data <- c(sdata, priors_all$exgaussian_rate)
stan_arguments$init <- exgaussian_inits
exgaussian_GHFM_mixture_IS <- do.call(BGHFM_models$exgaussian_rate$sample, stan_arguments)

# Fit Mixture-IS shifted-lognormal Hierarchical Factor Model
stan_arguments$data <- c(sdata, priors_all$shlognormal)
stan_arguments$init <- shlognormal_inits
shlognormal_GHFM_mixture_IS <- do.call(BGHFM_models$shifted_lognormal$sample, stan_arguments)

# Save Mixture-IS models
gaussian_HFM_mixture_IS$save_object(file = "Results/Stan/Viviani/Mixture IS models/viviani_gaussian_HFM_mixtureIS.rds")
exgaussian_GHFM_mixture_IS$save_object(file = "Results/Stan/Viviani/Mixture IS models/viviani_exgaussian_GHFM_mixtureIS.rds")
shlognormal_GHFM_mixture_IS$save_object(file = "Results/Stan/Viviani/Mixture IS models/viviani_shlognormal_GHFM_mixtureIS.rds")

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

# Estimate shifted-lognormal HFM Mixture-IS ELPDs
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
viviani_ELPDs <- list(
  ELPD_gaussian    = loo_mixis_gaussian,
  ELPD_exgaussian  = loo_mixis_exgaussian,
  ELPD_shlognormal = loo_mixis_shlognormal
)

# Save mixture ELPDs 
saveRDS(viviani_ELPDs, file = "Results/Rdata/ELPDs/viviani_HFM_ELPDs.rds")

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 11: Approximated ELPD in LMM (Zewotir & Galpin, 2005, p.160)
# ─────────────────────────────────────────────────────────────────────────────

# Download raw data from OSF
# download.file("https://osf.io/6qx98/download", mode = "wb",
#               destfile = "Data/Raw/Vivianni_2024_ManyStroopData.xlsx")

df <- read.xlsx(xlsxFile = "Data/Raw/Viviani_2024_ManyStroopData.xlsx") 

# Filter by experimental trial and remove block prefix
df <- df |> 
  filter(ExpTrial == 1) |> 
  mutate(Block = substr(Block, 5, nchar(Block)))

# Compute SSok flag per subkect
ss_ok <- df %>%
  group_by(SSID) %>%
  summarise(
    mean_iRT = mean(iRTok, na.rm = TRUE),
    mean_acc = mean(Correct, na.rm = TRUE)
  ) %>%
  mutate(
    zRT  = as.numeric(scale(mean_iRT)),
    ssOK = (abs(zRT) < 3) & (mean_acc > 0.75)
  ) %>%
  dplyr::select(SSID, ssOK)

# Join both data frames and filter by valid subjects
df_clean <- df %>%
  left_join(ss_ok, by = "SSID") %>%
  filter(
    !is.na(iRTok),
    !is.na(iRTpre),
    ssOK.x == 1,
    Correct == 1
  )

# Prepare variables for linear mixed models
df_clean <- df_clean %>%
  mutate(
    SubBlock_z = as.numeric(scale(SubBlockNum)),
    Trial_z    = as.numeric(scale(Trial)),
    iRTpre_z   = as.numeric(scale(iRTpre)),
    hRespC     = factor(RelevH, levels = c(-1, 1)),
    vRespC     = factor(RelevV, levels = c(-1, 1)),
    postERRC   = factor(PostERR, levels = c(0, 1)),
    CongC      = factor(CONG, levels = c(1, 0)),
    Block = factor(Block, levels = c("Peripheral", "Perifoveal",
                                     "Navon", "FigureGround",
                                     "Flanker", "Saliency"))
  )

# LMM with predictors
LMM_full <- lmer(formula = iRTok ~ Trial_z*SubBlock_z + iRTpre_z + postERRC
                 + hRespC + vRespC + Block*CongC + (Block*CongC|SSID), 
                 data = df_clean, 
                 REML = FALSE)

# LMM without predictors: we define a vraible with sum-to-zero contrast
df_clean$cond <- ifelse(df_clean$CONG == 1, -.5, .5)
LMM_basic <- lmer(formula = iRTok ~ 0 + Block + Block:cond + 
                    (0 + Block|SSID) + (0 + Block:cond|SSID), 
                  data = df_clean, 
                  control = lmerControl(optimizer = "bobyqa"),
                  REML = FALSE)

# Save both fitted models
saveRDS(list(LMM_full = LMM_full, LMM_basic = LMM_basic), 
        file = "Results/Rdata/LMM/viviani_LMM_models.rds")

# Compute ELPD LOO for each model 
ELPD_LMMs <- list(
  ELPD_LMM_full  = elpd_freq_LMM(fit = LMM_full, J_adj = TRUE),
  ELPD_LMM_basic = elpd_freq_LMM(fit = LMM_basic, J_adj = TRUE)
)

# Save LMM ELPDs
saveRDS(ELPD_LMMs, file = "Results/Rdata/ELPDs/viviani_LMM_ELPDs.rds")

# ─────────────────────────────────────────────────────────────────────────────
