# ╔════════════════════════════════════════════════════════════════════════════╗
# ║                             SCRIPT OVERVIEW                                ║
# ╠════════════════════════════════════════════════════════════════════════════╣
# ║ Script Name   : reanalysis_whitehead_exp1.R                                ║
# ║ Author        : Ricardo Rey-Sáez                                           ║
# ║ Role          : PhD Student in Psychology                                  ║
# ║ Institution   : Autonomous University of Madrid, Spain                     ║
# ║ Email         : ricardoreysaez95@gmail.es                                  ║
# ║ Date          : 15-11-2025                                                 ║
# ╚════════════════════════════════════════════════════════════════════════════╝

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 1: Description
# ─────────────────────────────────────────────────────────────────────────────
# Estimate GenHFMs using Rey-Mermet et al., data-inclusion criteria (Exp. 1)
# Original datasets, fitted models and scripts from Rey-Mermet et al. (2025)
# avalible here: https://osf.io/tyv9g/files/osfstorage 
# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 2: Load Packages
# ─────────────────────────────────────────────────────────────────────────────
# Libraries necessary for the script to function
library(cmdstanr)
library(loo)
library(posterior)
library(tidyverse)

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

# Load experiment 1 data
data <- read.csv("Data/Raw/whitehead_2020_Exp1_all.csv", header=TRUE,na.strings=c(""))

#subjects to include
includesubs <- c(507,508,513,514,517,518,519,520,521,506,515,523,524,525,526,527,528,529,530,532,
                 533,534,535,537,539,540,541,542,544,545,546,547,548,549,551,553,555,559,560,561,
                 562,567,568,569,571,573,576,577,578,579,583,584,585,586,587,588,591,592,593,594,
                 595,596,597,599,601,602,603,605,606,607,608,609,610,611,613,614,615,616,617,619,
                 620,621,622,623,624,625,626,629,630,631,633,634,635,637,639,640,641,642,643,645,
                 647,648,649,650,651,652,654,655,656,657,658,659,660,661,662,663,664,665,666,667,
                 668,669,670,671,672,673,674,675,676,680,678,679,681,683,685,687,689,690,692,693,
                 694,696,697,698,699,700,701,702,704,705,706,708,710,711,712,713,715,716,718,719,
                 720,721,723,724,725,726,728,729,730,732,733,734,736,737,738,739,740,741)

#data$StimSlideSimon.RT <- as.numeric(as.character(data$StimSlideSimon.RT))

#Filter and clean data
df.simon <- data %>% mutate(prevcon = lag(Congruency)) %>%
  mutate(StimSlideSimon.RT = as.integer(as.character(StimSlideSimon.RT)),
         StimSlideSimon.ACC = as.integer(as.character(StimSlideSimon.ACC)),
         BlockNum = as.integer(as.character(BlockNum))) %>%
  filter(Subject %in% includesubs & StimSlideSimon.RT != "" &
           #(StimSlideSimon.RT > 200 & StimSlideSimon.RT < 3000) &
           prevcon != 'NA' &
           BlockNum > 2) %>% 
  mutate(
    rt = StimSlideSimon.RT/1000,
    acc = StimSlideSimon.ACC,
    task = factor("simon")
  ) %>% 
  rename(subject = Subject,
         congruency = Congruency) %>% 
  dplyr::select(subject, BlockNum, task, congruency, acc, rt)


df.flanker <- data %>% mutate(prevcon = lag(Congruency)) %>%
  mutate(StimSlideFlanker.RT = as.integer(as.character(StimSlideFlanker.RT)),
         StimSlideFlanker.ACC = as.integer(as.character(StimSlideFlanker.ACC)),
         BlockNum = as.integer(as.character(BlockNum))) %>%
  filter(Subject %in% includesubs & StimSlideFlanker.RT != "" &
           #(StimSlideFlanker.RT > 200 & StimSlideFlanker.RT < 3000) &
           prevcon != 'NA' &
           BlockNum > 2) %>% 
  mutate(
    rt = StimSlideFlanker.RT/1000,
    acc = StimSlideFlanker.ACC,
    task = factor("flanker")
  ) %>% 
  rename(subject = Subject,
         congruency = Congruency) %>% 
  dplyr::select(subject, BlockNum, task, congruency, acc, rt)


df.stroop <- data %>% mutate(prevcon = lag(Congruency)) %>%
  mutate(StimSlideStroop.RT = as.integer(as.character(StimSlideStroop.RT)),
         StimSlideStroop.ACC = as.integer(as.character(StimSlideStroop.ACC)),
         BlockNum = as.integer(as.character(BlockNum))) %>%
  filter(Subject %in% includesubs &  StimSlideStroop.RT != "" &
           #(StimSlideStroop.RT > 200 & StimSlideStroop.RT < 3000) &
           prevcon != 'NA' &
           BlockNum > 2) %>% 
  mutate(
    rt = StimSlideStroop.RT/1000,
    acc = StimSlideStroop.ACC,
    task = factor("stroop")
  ) %>% 
  rename(subject = Subject,
         congruency = Congruency) %>% 
  dplyr::select(subject, BlockNum, task, congruency, acc, rt)


e1_1 <- bind_rows(df.simon, df.stroop, df.flanker) %>% 
  as_tibble() %>% 
  mutate(congruency = factor(congruency, levels = c(0, 1), 
                             labels = c("incongruent", "congruent")))

e1 <- e1_1 %>%
  filter(rt > 0.200 & rt < 3)

1 - nrow(e1) / nrow(e1_1) 

tasks <- unique(e1$task)

add_cont_vars <- function(data, cur_task) {
  data %>% mutate(
    "{cur_task}_b" := ifelse(task == cur_task, 1, 0),
    "{cur_task}_i" := ifelse(task == cur_task & congruency == "incongruent", 1, 0)
    #"{cur_task}_n" := ifelse(task == cur_task & congruency == "neutral", 1, 0)
  )
}

e1fit <- e1 %>% 
  mutate(corr = acc)
for (t in tasks) {
  e1fit <- add_cont_vars(e1fit, t)
}

# Save this data frame
save(e1fit, file = "R scripts/Empirical analysis/Rey-Mermet/datasets/e1fit.rdata")

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 5: Fit GenHFMs
# ─────────────────────────────────────────────────────────────────────────────

# Prepare Stan data
sdata.E1 <- Stan.list.data(
  data = e1fit, 
  subject_var = "subject", 
  task_var = "task", 
  condition_var = "congruency", 
  congruent_val = "congruent", 
  incongruent_val = "incongruent", 
  rt_var = "rt")

# Starting values: gaussian
gaussian_inits.E1 <- replicate(
  n = 6, expr = make_inits(
    data           = sdata_to_longdf(sdata.E1), 
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
exgaussian_inits.E1 <- replicate(
  n = 6, expr = make_inits(
    data           = sdata_to_longdf(sdata.E1), 
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
shlognormal_inits.E1 <- replicate(
  n = 6, expr = make_inits(
    data           = sdata_to_longdf(sdata.E1), 
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
sdata.E1$M <- 1
sdata.E1$Lambda_index <- matrix(1, nrow = 3)

# We're not interested in in LOO-CV right now.
sdata.E1$loo_mix_IS <- 0
sdata.E1$save_log_lik <- 0

# Fit Gaussian Hierarchical Factor model: Experiment 1
stan_arguments$data <- c(sdata.E1, priors_all$gaussian)
stan_arguments$init <- gaussian_inits.E1
gaussian_HFM_fit_E1 <- do.call(BGHFM_models$gaussian$sample, stan_arguments)

# Fit ex-Gaussian Hierarchical Factor model: Experiment 1
stan_arguments$data <- c(sdata.E1, priors_all$exgaussian_tau)
stan_arguments$init <- exgaussian_inits.E1
exgaussian_HFM_fit_E1 <- do.call(BGHFM_models$exgaussian_tau$sample, stan_arguments)

# Fit shifted-lognormal Hierarchical Factor model: Experiment 1
stan_arguments$data <- c(sdata.E1, priors_all$shlognormal)
stan_arguments$init <- shlognormal_inits.E1
shlognormal_HFM_fit_E1 <- do.call(BGHFM_models$shifted_lognormal$sample, stan_arguments)

# Save new fitted models
gaussian_HFM_fit_E1$save_object("Results/Stan/Rey-Mermet/gaussian_HFM_fit_E1.rds")
exgaussian_HFM_fit_E1$save_object("Results/Stan/Rey-Mermet/exgaussian_HFM_fit_E1.rds")
shlognormal_HFM_fit_E1$save_object("Results/Stan/Rey-Mermet/shlognormal_HFM_fit_E1.rds")

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION extra: PSIS-LOO estimates
# ─────────────────────────────────────────────────────────────────────────────

# PSIS-LOO: gaussian HFM
loo_gaussian_HFM <- loo_BGHFM(
  fit = gaussian_HFM_fit_E1,
  data           = sdata_to_longdf(sdata.E1), 
  subject_var    = "subject", 
  task_var       = "task", 
  condition_var  = "condition",
  rt_var         = "rt", 
  n.cores = 8,  
  model = "gaussian",
  estimate_r_eff = TRUE)

# PSIS-LOO: exgaussian HFM
loo_exgaussian_HFM <- loo_BGHFM(
  fit = exgaussian_HFM_fit_E1,
  data           = sdata_to_longdf(sdata.E1), 
  subject_var    = "subject", 
  task_var       = "task", 
  condition_var  = "condition",
  rt_var         = "rt", 
  n.cores = 8,  
  model = "exgaussian",
  estimate_r_eff = TRUE)

# PSIS-LOO: exgaussian HFM
loo_shlognormal_HFM <- loo_BGHFM(
  fit = shlognormal_HFM_fit_E1,
  data           = sdata_to_longdf(sdata.E1), 
  subject_var    = "subject", 
  task_var       = "task", 
  condition_var  = "condition",
  rt_var         = "rt", 
  n.cores = 8,  
  model = "shifted.lognormal",
  estimate_r_eff = TRUE)

# PSIS-LOO: Wiener difussion model
load("R scripts/Empirical analysis/DDM_fitted_models/fit_wh_e1_1.rda")
rm(pred_fit_wh_e1_1)
load("R scripts/Empirical analysis/DDM_fitted_models/fit_wh_e1_2.rda")
rm(pred_fit_wh_e1_2)

# difference coding
loo_wiener_dc_exp1 <- loo_wiener_est_dc(fit = fit_wh_e1_1, data = fit_wh_e1_1$data, 
                                        n.cores = 8, estimate_r_eff = TRUE)
loo_wiener_cc_exp1 <- loo_wiener_est_cc(fit = fit_wh_e1_2, data = fit_wh_e1_2$data, 
                                        n.cores = 8, estimate_r_eff = TRUE)

# Save loo-results
loo_reanalysis_whitehead_exp1 <- list(
  gaussian_HFM_loo_exp1 = loo_gaussian_HFM,
  exgaussian_HFM_loo_exp1 = loo_exgaussian_HFM,
  shlognormal_HFM_loo_exp1 = loo_shlognormal_HFM,
  wiener_dc_loo_exp1 = loo_wiener_dc_exp1,
  wiener_cc_loo_exp1 = loo_wiener_cc_exp1
)

# PSIS-LOO estimates
saveRDS(object = loo_reanalysis_whitehead_exp1, 
        file = "Results/Rdata/Supplemental material/loo_reanalysis_whitehead_exp1.rds")

# ─────────────────────────────────────────────────────────────────────────────