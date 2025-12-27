# ╔════════════════════════════════════════════════════════════════════════════╗
# ║                             SCRIPT OVERVIEW                                ║
# ╠════════════════════════════════════════════════════════════════════════════╣
# ║ Script Name   : reanalysis_whitehead_exp3.R                                ║
# ║ Author        : Ricardo Rey-Sáez                                           ║
# ║ Role          : PhD Student in Psychology                                  ║
# ║ Institution   : Autonomous University of Madrid, Spain                     ║
# ║ Email         : ricardoreysaez95@gmail.es                                  ║
# ║ Date          : 15-11-2025                                                 ║
# ╚════════════════════════════════════════════════════════════════════════════╝

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 1: Description
# ─────────────────────────────────────────────────────────────────────────────
# Estimate GenHFMs using Rey-Mermet et al., data-inclusion criteria (Exp. 3)
# Original datasets, fitted models and scripts from Rey-Mermet et al. (2025)
# avalible here: https://osf.io/tyv9g/files/osfstorage 
# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 2: Load Packages
# ─────────────────────────────────────────────────────────────────────────────
# Libraries necessary for the script to function
library(cmdstanr)
library(posterior)
library(tidyverse)
library(loo)

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 3: User-Defined Functions
# ─────────────────────────────────────────────────────────────────────────────

# Extra-functions to analyze data
source("R scripts/R functions/empirical_functions.R")

# Save all BGHFM in one list
BGHFM_models <- list(
  gaussian          = cmdstan_model("Stan models/BGHFM/C_BGHFM_gaussian.stan"),
  exgaussian_tau    = cmdstan_model("Stan models/BGHFM/C_BGHFM_exgaussian_tau.stan"),
  shifted_lognormal = cmdstan_model("Stan models/BGHFM/C_BGHFM_shlognormal.stan")
)

# Diffuse prior distributions (same as simulation study)
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

# Global Stan arguments
stan_arguments <- list(
  data = NULL, 
  iter_warmup = 500,
  iter_sampling = 500,
  chains = 6, 
  parallel_chains = 6,
  adapt_delta = 0.80,
  max_treedepth = 10,
  seed = 2025,
  refresh = 50,
  init = NULL)

# Pointwise wiener PSIS-LOO: difference coding
loo_wiener_est_dc <- function(fit, data, n.cores, estimate_r_eff = FALSE) {
  
  # Wiener log-density function for data_i
  llfun_wiener <- function(data_i, draws, log = TRUE) {
    
    # Data from draw i
    subj_i <- data_i$subject
    task_i <- data_i$task
    rt_i   <- data_i$rt
    corr_i <- ifelse(data_i$corr == 1, "upper", "lower")
    
    # Dummy variables 
    fb <- as.numeric(data_i[[paste0(task_i, "_b")]])
    fi <- as.numeric(data_i[[paste0(task_i, "_i")]])
    
    # Posterior draws column ids
    taskname <- paste0("task",task_i)
    mu_b_ids     <- paste0("b_", c(task_i), c("_b", "_i"))
    mu_bs_ids    <- paste0("b_bs_", taskname)
    mu_ndt_ids   <- paste0("b_ndt_", taskname)
    r_mu_b_ids   <- paste0("r_subject[", subj_i,",",task_i,c("_b","_i"),"]")
    r_mu_bs_ids  <- paste0("r_subject__bs[", subj_i,",",taskname,"]")
    r_mu_ndt_ids <- paste0("r_subject__ndt[", subj_i,",",taskname,"]")
    
    # Boundary separation
    a  <- draws[, mu_bs_ids] + draws[, r_mu_bs_ids]
    # Drift rates
    v  <- fb * (draws[, mu_b_ids[1]] + draws[, r_mu_b_ids[1]]) + 
          fi * (draws[, mu_b_ids[2]] + draws[, r_mu_b_ids[2]])
    # Non-decision time
    t0 <- draws[, mu_ndt_ids] + draws[, r_mu_ndt_ids]
    
    # Density values
    dens <- rtdists::ddiffusion(
      rt = rep(rt_i, length(a)), 
      response = corr_i, 
      a = a, 
      v = v, 
      t0 = t0, 
      z = a * .5, 
      precision = 4)
    
    # Return log-density if log = TRUE, otherwise return the density
    if (log) {
      return(log(dens))
    } else {
      return(dens)
    }
  }
  
  # Posterior draws
  posterior_draws <- posterior::as_draws_matrix(fit)
  
  # Compute relative efficiency
  if(estimate_r_eff){
    r_eff <- relative_eff(llfun_wiener, 
                          log = FALSE, 
                          chain_id = rep(1:brms::nchains(fit), each = nrow(posterior_draws) / brms::nchains(fit)), 
                          data = data, 
                          draws = posterior_draws, 
                          cores = n.cores)
  } else {
    r_eff <- 1
  }
  
  # Compute PSIS-LOO (and save time)
  loo_model <- loo(llfun_wiener, draws = posterior_draws, 
                   r_eff = r_eff, data = data, cores = n.cores)
  
  # Return LOO and time
  return(loo_model)
}

# Pointwise wiener PSIS-LOO: condition coding
loo_wiener_est_cc <- function(fit, data, n.cores, estimate_r_eff = FALSE) {
  
  # Wiener log-density function for data_i
  llfun_wiener <- function(data_i, draws, log = TRUE) {
    
    # Data from draw i
    subj_i <- data_i$subject
    task_i <- data_i$task
    cond_i <- data_i$congruency
    rt_i   <- data_i$rt
    corr_i <- ifelse(data_i$corr == 1, "upper", "lower")
    
    # Task and condition names
    taskname <- paste0("task",task_i)
    condname <- paste0("congruency", cond_i)
    
    # Posterior draws column ids
    mu_b_ids     <- paste0("b_", condname, ":", taskname)
    mu_bs_ids    <- paste0("b_bs_", taskname)
    mu_ndt_ids   <- paste0("b_ndt_", taskname)
    r_mu_b_ids   <- paste0("r_subject[", subj_i,",",condname, ":", taskname,"]")
    r_mu_bs_ids  <- paste0("r_subject__bs[", subj_i,",", taskname,"]")
    r_mu_ndt_ids <- paste0("r_subject__ndt[", subj_i,",", taskname,"]")
    
    # Boundary separation
    a  <- draws[, mu_bs_ids] + draws[, r_mu_bs_ids]
    # Drift rates
    v  <- draws[, mu_b_ids] + draws[, r_mu_b_ids]
    # Non-decision time
    t0 <- draws[, mu_ndt_ids] + draws[, r_mu_ndt_ids]
    
    # Density values
    dens <- rtdists::ddiffusion(
      rt = rep(rt_i, length(a)), 
      response = corr_i, 
      a = a, 
      v = v, 
      t0 = t0, 
      z = a * .5, 
      precision = 4)
    
    # Return log-density if log = TRUE, otherwise return the density
    if (log) {
      return(log(dens))
    } else {
      return(dens)
    }
  }
  
  # Posterior draws
  posterior_draws <- posterior::as_draws_matrix(fit)
  
  # Compute relative efficiency
  if(estimate_r_eff){
    r_eff <- relative_eff(llfun_wiener, 
                          log = FALSE, 
                          chain_id = rep(1:brms::nchains(fit), each = nrow(posterior_draws) / brms::nchains(fit)), 
                          data = data, 
                          draws = posterior_draws, 
                          cores = n.cores)
  } else {
    r_eff <- 1
  }
  
  # Compute PSIS-LOO (and save time)
  loo_model <- loo(llfun_wiener, draws = posterior_draws, 
                   r_eff = r_eff, data = data, cores = n.cores)
  
  # Return LOO and time
  return(loo_model)
}

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 4: Rey-Mermet, Singmann and Oberauer data settings
# ─────────────────────────────────────────────────────────────────────────────

data.simon<- read.csv("Data/Raw/whitehead_2020_Exp3_simon.csv", header=TRUE,na.strings=c(""))
data.stroop<- read.csv("Data/Raw/whitehead_2020_Exp3_stroop.csv", header=TRUE,na.strings=c(""))
data.flanker<- read.csv("Data/Raw/whitehead_2020_Exp3_flanker.csv", header=TRUE,na.strings=c(""))


#subjects to include
includesubs <- c(105,108,109,110,111,112,113,114,115,116,117,118,119,120,121,
                 122,123,124,125,126,127,128,129,130,
                 131,132,133,134,135,136,137,138,139,140,141,142,143,144,145,
                 146,147,148,149,150,151,152,153,154,
                 155,156,157,158,159,160,161,162,163,165,166,167,168,169,170,
                 171,172,173,174,175,176,177,178,179,
                 180,181,182,183,184,185,187,188,189,190,191,192,193,194,195,
                 196,197,198,199,200,201,202,203,204,
                 205,206,207,208,209,210,211,212,213,214,215,216,217,218,219,
                 220,221,222,223,224,225,226,228,229,
                 231,232,233,234,235,236,237,238,240,241,242,243,244,246,247,
                 248,249,250,251,253,254,255,256,257,
                 259,260,261,262,263,264,265,266,267,268,269,270,271,273,274,
                 275,277,278,279,280,281,282,283,284,
                 285,286,287,288,289,290,291,293,294,295,296,297,299,301,303,
                 305,306,307,308,309,311,312,313,314,
                 316,317,318,319,320,321,322,323,324,325,326,327,328,330,331,
                 332,333,334)
#data$StimSlideSimon.RT <- as.numeric(as.character(data$StimSlideSimon.RT))

#Filter and clean data
df.simon <- data.simon %>% 
  mutate(
    StimSlideSimon.RT = as.integer(as.character(StimSlideSimon.RT)),
    StimSlideSimon.ACC = as.integer(as.character(StimSlideSimon.ACC))
  ) %>%
  filter(Subject %in% includesubs &
           #(StimSlideSimon.RT > 200 & StimSlideSimon.RT < 3000) &
           #StimSlideSimon.ACC == 1 &
           PracExp == "Exp") %>% 
  mutate(
    rt = StimSlideSimon.RT/1000,
    acc = StimSlideSimon.ACC,
    task = factor("simon")
  ) %>% 
  rename(subject = Subject,
         congruency = Congruency) %>% 
  dplyr::select(subject, task, congruency, acc, rt)

table(df.simon$acc)

str(df.simon)

df.flanker <- data.flanker %>% 
  mutate(
    StimSlideFlanker.RT = as.integer(as.character(StimSlideFlanker.RT)),
    StimSlideFlanker.ACC = as.integer(as.character(StimSlideFlanker.ACC))
  ) %>%
  filter(Subject %in% includesubs &
           #(StimSlideFlanker.RT > 200 & StimSlideFlanker.RT < 3000) &
           #StimSlideFlanker.ACC == 1 &
           PracExp == "Exp") %>% 
  mutate(
    rt = StimSlideFlanker.RT/1000,
    acc = StimSlideFlanker.ACC,
    task = factor("flanker")
  ) %>% 
  rename(subject = Subject,
         congruency = Congruency) %>% 
  dplyr::select(subject, task, congruency, acc, rt)

str(df.flanker)
table(df.flanker$acc)


df.stroop <- data.stroop %>% 
  mutate(StimSlideStroop.RT = as.integer(as.character(StimSlideStroop.RT)),
         StimSlideStroop.ACC = as.integer(as.character(StimSlideStroop.ACC)),
         Subject = as.integer(as.character(Subject))) %>%
  filter(Subject %in% includesubs & 
           #(StimSlideStroop.RT > 200 & StimSlideStroop.RT < 3000) &
           #StimSlideStroop.ACC == 1 &
           PracExp == "Exp") %>% 
  mutate(
    rt = StimSlideStroop.RT/1000,
    acc = StimSlideStroop.ACC,
    task = factor("stroop")
  ) %>% 
  rename(subject = Subject,
         congruency = Congruency) %>% 
  dplyr::select(subject, task, congruency, acc, rt)

e3_1 <- bind_rows(df.simon, df.stroop, df.flanker) %>% 
  as_tibble() %>% 
  mutate(congruency = factor(congruency, levels = c(0, 1), 
                             labels = c("incongruent", "congruent")))
e3 <- e3_1 %>%
  filter(rt > 0.200 & rt < 3)

1 - nrow(e3) / nrow(e3_1) 

# %>% 
#   group_by(task, congruency) %>% 
#   mutate(weight = 1/n()) %>% 
#   mutate(macc = abs((1 - acc) - mean(acc)))
#   ungroup()
e3


###### descriptives 
e3 %>% 
  group_by(task, congruency) %>% 
  summarise(macc = mean(acc), 
            mrt = mean(rt)) %>% 
  ungroup

e3 %>% 
  group_by(task, congruency) %>% 
  count()

ggplot(e3, aes(rt, color = congruency, fill = congruency,
               group = congruency)) +
  geom_density(alpha = 0.1, mapping = aes(y = after_stat(count))) +
  facet_wrap(~acc+task)

e3 %>% 
  summarise(n = n_distinct(subject))

# Save this data frame
save(e3, file = "R scripts/Empirical analysis/Rey-Mermet/datasets/e3.rdata")

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 5: Fit GenHFMs
# ─────────────────────────────────────────────────────────────────────────────

# Prepare Stan data
sdata.E3 <- Stan.list.data(
  data = e3, 
  subject_var = "subject", 
  task_var = "task", 
  condition_var = "congruency", 
  congruent_val = "congruent", 
  incongruent_val = "incongruent", 
  rt_var = "rt")

# Starting values: gaussian
gaussian_inits.E3 <- replicate(
  n = 6, expr = make_inits(
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

# Starting values: exgaussian
exgaussian_inits.E3 <- replicate(
  n = 6, expr = make_inits(
    data           = sdata_to_longdf(sdata.E3), 
    subject_var    = "subject", 
    task_var       = "task", 
    condition_var  = "condition",
    rt_var         = "rt", 
    congruent_id   = 1, 
    incongruent_id = 2, 
    M              = 1, 
    model          = "exgaussian"), 
  simplify = FALSE
)

# Starting values: shifted-lognormal
shlognormal_inits.E3 <- replicate(
  n = 6, expr = make_inits(
    data           = sdata_to_longdf(sdata.E3), 
    subject_var    = "subject", 
    task_var       = "task", 
    condition_var  = "condition",
    rt_var         = "rt", 
    congruent_id   = 1, 
    incongruent_id = 2, 
    M              = 1, 
    model          = "shifted.lognormal"), 
  simplify = FALSE
)

# Add latent factos and Lambda index
sdata.E3$M <- 1
sdata.E3$Lambda_index <- matrix(1, nrow = 3)

# We're not interested in in LOO-CV right now.
sdata.E3$loo_mix_IS <- 0
sdata.E3$save_log_lik <- 0

# Fit Gaussian Hierarchical Factor model: Experiment 1
stan_arguments$data <- c(sdata.E3, priors_all$gaussian)
stan_arguments$init <- gaussian_inits.E3
gaussian_HFM_fit_E3 <- do.call(BGHFM_models$gaussian$sample, stan_arguments)

# Fit ex-Gaussian Hierarchical Factor model: Experiment 1
stan_arguments$data <- c(sdata.E3, priors_all$exgaussian_tau)
stan_arguments$init <- exgaussian_inits.E3
exgaussian_HFM_fit_E3 <- do.call(BGHFM_models$exgaussian_tau$sample, stan_arguments)

# Fit shifted-lognormal Hierarchical Factor model: Experiment 1
stan_arguments$data <- c(sdata.E3, priors_all$shlognormal)
stan_arguments$init <- shlognormal_inits.E3
shlognormal_HFM_fit_E3 <- do.call(BGHFM_models$shifted_lognormal$sample, stan_arguments)

# Save new fitted models
gaussian_HFM_fit_E3$save_object("Results/Stan/Rey-Mermet/gaussian_HFM_fit_E3.rds")
exgaussian_HFM_fit_E3$save_object("Results/Stan/Rey-Mermet/exgaussian_HFM_fit_E3.rds")
shlognormal_HFM_fit_E3$save_object("Results/Stan/Rey-Mermet/shlognormal_HFM_fit_E3.rds")

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION extra: PSIS-LOO estimates
# ─────────────────────────────────────────────────────────────────────────────

# PSIS-LOO: gaussian HFM
loo_gaussian_HFM <- loo_BGHFM(
  fit = gaussian_HFM_fit_E3,
  data           = sdata_to_longdf(sdata.E3), 
  subject_var    = "subject", 
  task_var       = "task", 
  condition_var  = "condition",
  rt_var         = "rt", 
  n.cores = 8,  
  model = "gaussian",
  estimate_r_eff = TRUE)

# PSIS-LOO: exgaussian HFM
loo_exgaussian_HFM <- loo_BGHFM(
  fit = exgaussian_HFM_fit_E3,
  data           = sdata_to_longdf(sdata.E3), 
  subject_var    = "subject", 
  task_var       = "task", 
  condition_var  = "condition",
  rt_var         = "rt", 
  n.cores = 8,  
  model = "exgaussian",
  estimate_r_eff = TRUE)

# PSIS-LOO: exgaussian HFM
loo_shlognormal_HFM <- loo_BGHFM(
  fit = shlognormal_HFM_fit_E3,
  data           = sdata_to_longdf(sdata.E3), 
  subject_var    = "subject", 
  task_var       = "task", 
  condition_var  = "condition",
  rt_var         = "rt", 
  n.cores = 8,  
  model = "shifted.lognormal",
  estimate_r_eff = TRUE)

# PSIS-LOO: Wiener difussion model
load("R scripts/Empirical analysis/DDM_fitted_models/fit_wh_e3_1.rda")
rm(pred_fit_wh_e3_1)
load("R scripts/Empirical analysis/DDM_fitted_models/fit_wh_e3_2.rda")
rm(pred_fit_wh_e3_2)

# difference coding
loo_wiener_dc_exp3 <- loo_wiener_est_dc(fit = fit_wh_e3_1, data = fit_wh_e3_1$data, 
                                        n.cores = 8, estimate_r_eff = TRUE)
loo_wiener_cc_exp3 <- loo_wiener_est_cc(fit = fit_wh_e3_2, data = fit_wh_e3_2$data, 
                                        n.cores = 8, estimate_r_eff = TRUE)

# Save loo-results
loo_reanalysis_whitehead_exp3 <- list(
  gaussian_HFM_loo_exp3 = loo_gaussian_HFM,
  exgaussian_HFM_loo_exp3 = loo_exgaussian_HFM,
  shlognormal_HFM_loo_exp3 = loo_shlognormal_HFM,
  wiener_dc_loo_exp3 = loo_wiener_dc_exp3,
  wiener_cc_loo_exp3 = loo_wiener_cc_exp3
)

# PSIS-LOO estimates
saveRDS(object = loo_reanalysis_whitehead_exp3, 
        file = "Results/Rdata/Supplemental material/loo_reanalysis_whitehead_exp3.rds")

# ─────────────────────────────────────────────────────────────────────────────
