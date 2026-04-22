# ╔════════════════════════════════════════════════════════════════════════════╗
# ║                             SCRIPT OVERVIEW                                ║
# ╠════════════════════════════════════════════════════════════════════════════╣
# ║ Script Name   : fit_BLMMs.R                                                ║
# ║ Author        : Ricardo Rey-Sáez                                           ║
# ║ Role          : PhD Student in Psychology                                  ║
# ║ Institution   : Autonomous University of Madrid, Spain                     ║
# ║ Email         : ricardoreysaez95@gmail.es                                  ║
# ║ Date          : 14-04-2025                                                 ║
# ╚════════════════════════════════════════════════════════════════════════════╝

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 1: Description
# ─────────────────────────────────────────────────────────────────────────────
# Fit Bayesian Linear Mixed Models for each experimental study
# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 2: Load Packages
# ─────────────────────────────────────────────────────────────────────────────
# Libraries necessary for the script to function
library(acdcquery) # Load SQL datasets
library(cmdstanr)  # Stan bayesian LMM
library(posterior) # Posterior summaries
library(doFuture)  # Parallelized estimation
library(progressr) # Progress bar in future framework

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 3: User-Defined Functions
# ─────────────────────────────────────────────────────────────────────────────

# Load meta-analysis user-defined functions
source("R scripts/R functions/meta_analysis_functions.R")

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 4: Load and prepare data
# Based on Haat ef al. (2025) code: 
# https://github.com/jstbcs/acdc-paper/blob/main/paper/p.Rmd
# ─────────────────────────────────────────────────────────────────────────────

# Download full SQL database
# download.file(url = "https://github.com/jstbcs/acdc-database/raw/main/acdc.db",
#               destfile = "Data/Raw/acdc.db", mode = "wb")

# Make SQL connection
conn <- connect_to_db("Data/Raw/acdc.db")

# Prepare arguments
arguments <-  list()
arguments <- add_argument(
  list = arguments,
  conn = conn,
  variable = "study_id",
  operator = "greater",
  values = c(-1)
)

# Data frame
query_results <- query_db(
  conn = conn,
  arguments = arguments,
  target_vars = c("default", "publication_id", "publication_code", "task_name")
)

# Select only Stroop, Flanker and Simon tasks
all.df <- query_results[which(query_results$task_name %in% c("stroop", "flanker", "simon")),]

# Check observed response distributions
cnt <- 1
for(i in unique(all.df$dataset_id)) {
  study <- unique(all.df$publication_code[which(all.df$dataset_id == i)])
  temp.df <- clean.data(all.df[which(all.df$dataset_id == i),])
  hist(temp.df$rt, main = paste("Dataset", cnt, ":", study), breaks = 100)
  cnt <- cnt + 1
}

# Remove dataset 1 (chetverikov et al., 2017): weird observer RT distribution
all.df <- all.df[-which(all.df$publication_id == 1),]

# Relabel dataset id
all.df$dataset_id <- as.integer(as.factor(all.df$dataset_id))

# Number of experiments
length(unique(all.df$dataset_id))

# Number of subjects
length(unique(all.df$subject))

# Number of experiments per task
rowSums(table(all.df$task_name, all.df$dataset_id)!=0)

# Save data for meta-analysis
save(all.df, file = "Data/Processed/acdc_data.rdata")

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 5: Fit all Bayesian Linear Mixed Models to each experimental study
# Method: Classic MCMC bayesian estimation
# ─────────────────────────────────────────────────────────────────────────────
# Fit Gaussian Linear Mixed Model
gaussian.blmm <- fit_LMM_models(data = all.df, model = "gaussian", 
                                method = "MCMC")
# Fit Exgaussian Linear Mixed Model: tau parameterization
exgaussian.blmm.tau <- fit_LMM_models(data = all.df, model = "exgaussian_tau", 
                                      method = "MCMC") 
# Fit Exgaussian Linear Mixed Model: rate parameterization
exgaussian.blmm.rate <- fit_LMM_models(data = all.df, model = "exgaussian_rate", 
                                       method = "MCMC")
# Fit Shifted-Lognormal Linear Mixed Model
shlognormal.blmm <- fit_LMM_models(data = all.df, model = "shifted.lognormal", 
                                   method = "MCMC")

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 6: Fit all Bayesian Linear Mixed Models to each experimental study
# Method: Laplace approximation
# ─────────────────────────────────────────────────────────────────────────────

# Fit Gaussian Linear Mixed Model
gaussian.laplace <- fit_LMM_models(data = all.df, model = "gaussian", 
                                   method = "Laplace")
# Fit Exgaussian Linear Mixed Model: tau parameterization
exgaussian.laplace.tau <- fit_LMM_models(data = all.df, model = "exgaussian_tau", 
                                         method = "Laplace") 
# Fit Exgaussian Linear Mixed Model: rate parameterization
exgaussian.laplace.rate <- fit_LMM_models(data = all.df, model = "exgaussian_rate", 
                                          method = "Laplace")
# Fit Shifted-Lognormal Linear Mixed Model
shlognormal.laplace <- fit_LMM_models(data = all.df, model = "shifted.lognormal", 
                                      method = "Laplace")

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 7: Save Parameter Estimates Results
# ─────────────────────────────────────────────────────────────────────────────

# Save author, task and other informative variables for each experiment
author <- task <- n_subj <- n_cong_trials <- n_incong_trials <- vector(length = length(unique(all.df$dataset_id))) 
for(i in unique(all.df$dataset_id)){
  author[i]   <- unique(all.df$publication_code[which(all.df$dataset_id==i)])
  task[i]     <- unique(all.df$task_name[which(all.df$dataset_id==i)])
  n_subj[i]   <- length(unique(all.df$subject[which(all.df$dataset_id==i)]))
  n_cong_trials[i] <- length(unique(all.df$trial[which(all.df$dataset_id==i & all.df$congruency== "congruent")]))
  n_incong_trials[i] <- length(unique(all.df$trial[which(all.df$dataset_id==i & all.df$congruency== "incongruent")]))
}

# Save this variables in a basic data frame
base_df <- data.frame(author = author, task = task, n_subj = n_subj,
                      n_cong_trials = n_cong_trials, 
                      n_incong_trials = n_incong_trials)

# Create data frames and save it in a list
parameter_estimates <- list(
  # Gaussian parameter estimates
  gaussian = cbind(base_df, as.data.frame(gaussian.blmm)), 
  # Exgaussian parameter estimates (tau)
  exgaussian_tau = cbind(base_df, as.data.frame(exgaussian.blmm.tau)),
  # Exgaussian parameter estimates (rate)
  exgaussian_rate = cbind(base_df, as.data.frame(exgaussian.blmm.rate)),
  # Shifted lognormal parameter estimates
  shifted.lognormal = cbind(base_df, as.data.frame(shlognormal.blmm))
  )

# Also Laplace approximation results
laplace_estimates <- list(
  # Gaussian parameter estimates
  gaussian = cbind(base_df, as.data.frame(gaussian.laplace)), 
  # Exgaussian parameter estimates (tau)
  exgaussian_tau = cbind(base_df, as.data.frame(exgaussian.laplace.tau)),
  # Exgaussian parameter estimates (rate)
  exgaussian_rate = cbind(base_df, as.data.frame(exgaussian.laplace.rate)),
  # Shifted lognormal parameter estimates
  shifted.lognormal = cbind(base_df, as.data.frame(shlognormal.laplace))
)

# Save results
saveRDS(object = parameter_estimates, file = "Results/Rdata/Meta-parameters/LMM_full_parestimates.rds")
saveRDS(object = laplace_estimates,   file = "Results/Rdata/Meta-parameters/LMM_laplace_parestimates.rds")

# Let's make checks
parameter_estimates <- readRDS(file = "Results/Rdata/LMM_full_parestimates.rds")
parameter_estimates <- readRDS(file = "Results/Rdata/LMM_laplace_parestimates.rds")

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 8: Minimal Checks on Model Estimates
# ─────────────────────────────────────────────────────────────────────────────

# Results by task: parameters posterior means
lapply(parameter_estimates, function(x) {
  x |> group_by(task) |> select(ends_with(match = "postMean")) |> summarise_all(mean)
})

# Results by task: parameters posterior standard deviation
lapply(parameter_estimates, function(x) {
  x |> group_by(task) |> select(ends_with(match = "postSd")) |> summarise_all(mean)
})

# Number of divergent chains
lapply(parameter_estimates, function(x) {
  x |> group_by(task) |> select(starts_with(match = "num_divergent")) |> summarise_all(sum)
})

# Maximum treedepth
lapply(parameter_estimates, function(x) {
  x |> group_by(task) |> select(starts_with(match = "num_max")) |> summarise_all(mean)
})

# Check that maximum treedepth values
parameter_estimates$shifted.lognormal |> 
  select("avg_rhat", "disp_rhat", starts_with("num_max")) |> 
  filter(num_max_treedepth1 > 0 | 
         num_max_treedepth2 > 0 |
         num_max_treedepth3 > 0) # Only two draws of 12.000, and Rhats are okey!

# EBFMI checks
lapply(parameter_estimates, function(x) {
  x |> group_by(task) |> select(starts_with(match = "ebfmi")) |> summarise_all(mean)
})

# Average R-hats
lapply(parameter_estimates, function(x) {
  x |> group_by(task) |> select(starts_with(match = "avg") | starts_with(match = "disp")) |> summarise_all(mean)
})


# ─────────────────────────────────────────────────────────────────────────────
