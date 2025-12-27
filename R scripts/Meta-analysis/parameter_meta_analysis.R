# ╔════════════════════════════════════════════════════════════════════════════╗
# ║                             SCRIPT OVERVIEW                                ║
# ╠════════════════════════════════════════════════════════════════════════════╣
# ║ Script Name   : paramete_meta_analysis.R                                   ║
# ║ Author        : Ricardo Rey-Sáez                                           ║
# ║ Role          : PhD Student in Psychology                                  ║
# ║ Institution   : Autonomous University of Madrid, Spain                     ║
# ║ Email         : ricardoreysaez95@gmail.es                                  ║
# ║ Date          : 22-04-2025                                                 ║
# ╚════════════════════════════════════════════════════════════════════════════╝

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 1: Description
# ─────────────────────────────────────────────────────────────────────────────
# Meta-analysis for each model parameter
# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 2: Load Packages
# ─────────────────────────────────────────────────────────────────────────────
# Libraries necessary for the script to function
library(metafor)   # Meta-analysis
library(ggpubr)    # Quick Q-Q plots
library(tidyverse) # Quick aggregate analysis
library(knitr)     # Meta-analytic results table
library(kableExtra)# Meta-analytic results table
library(stringr)   # string helpers                  
# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 3: User-Defined Functions
# ─────────────────────────────────────────────────────────────────────────────

# Function to create CIs
CIs <- function(parameter, digits = 3, log_param = FALSE) {
  lb <- parameter$ci.lb
  ub <- parameter$ci.ub
  if (log_param) {
    lb <- exp(lb)
    ub <- exp(ub)
  }
  # Number of digits:
  fmt <- paste0("%.", digits, "f")
  # Apply it for each element:
  ci <- mapply(function(l, u) {
    paste0("[", sprintf(fmt, l), ", ", sprintf(fmt, u), "]")
  }, l = lb, u = ub)
  return(ci)
}

# Number format for meta-analytic table
fmt_pair <- function(est, ci, digits = 3) {
  sprintf("%.*f %s", digits, est, ci)   # e.g. 0.123 [0.045, 0.202]
}

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 4: Load BLMM parameter estimates
# ─────────────────────────────────────────────────────────────────────────────

# Load LMM parameter estimates for each dataset
parameter_estimates <- readRDS(file = "Results/Rdata/Meta-parameters/LMM_full_parestimates.rds")

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 5: Meta-analysis of the gaussian model
# Note: Fixed effects are means and random effects are (log) variances
# ─────────────────────────────────────────────────────────────────────────────

# Selected cols
selected_cols <- c("author", "task", "n_subj", "n_cong_trials", "n_incong_trials",
                   # Fixed effects and uncertainty
                   "mu_alpha_postMean","mu_alpha_postSd", 
                   "mu_theta_postMean", "mu_theta_postSd", 
                   "logmu_sigma_postMean", "logmu_sigma_postSd",
                   "logshift_sigma_postMean", "logshift_sigma_postSd",
                   # Random effects and uncertainty
                   "logsd_alpha_postMean","logsd_alpha_postSd", 
                   "logsd_theta_postMean", "logsd_theta_postSd", 
                   "logsd_sigma_postMean", "logsd_sigma_postSd",
                   "logscale_sigma_postMean", "logscale_sigma_postSd")

# Final gaussian data frame
gaussian_estimates <- parameter_estimates$gaussian[,selected_cols]
gaussian_estimates$author <- as.factor(sapply(strsplit(parameter_estimates$gaussian$author, "_"), function(x){
  paste0(x[1], " (", x[2], ")")
}))

# Change posterior sds to posterior variances
gaussian_estimates$mu_alpha_postVar <- gaussian_estimates$mu_alpha_postSd^2
gaussian_estimates$mu_theta_postVar <- gaussian_estimates$mu_theta_postSd^2
gaussian_estimates$logmu_sigma_postVar <- gaussian_estimates$logmu_sigma_postSd^2
gaussian_estimates$logshift_sigma_postVar <- gaussian_estimates$logshift_sigma_postSd^2
gaussian_estimates$logsd_alpha_postVar <- gaussian_estimates$logsd_alpha_postSd^2
gaussian_estimates$logsd_theta_postVar <- gaussian_estimates$logsd_theta_postSd^2
gaussian_estimates$logsd_sigma_postVar <- gaussian_estimates$logsd_sigma_postSd^2
gaussian_estimates$logscale_sigma_postVar <- gaussian_estimates$logscale_sigma_postSd^2

# Check normality assumption
ggqqplot(data = gaussian_estimates, x = "mu_alpha_postMean")
ggqqplot(data = gaussian_estimates, x = "mu_theta_postMean")
ggqqplot(data = gaussian_estimates, x = "logsd_alpha_postMean")
ggqqplot(data = gaussian_estimates, x = "logsd_theta_postMean")
ggqqplot(data = gaussian_estimates, x = "logmu_sigma_postMean")
ggqqplot(data = gaussian_estimates, x = "logsd_sigma_postMean")
ggqqplot(data = gaussian_estimates, x = "logshift_sigma_postMean")
ggqqplot(data = gaussian_estimates, x = "logscale_sigma_postMean")

# Meta-analytic common arguments
common_args <- list(random = ~ factor(task) - 1 | author,
                    mods = ~ factor(task) - 1,
                    data = gaussian_estimates,
                    struct = "DIAG",
                    method = "REML",
                    test = "knha")

# Specify each model parameter and uncertainty
model_specs <- list(
  list(yi = "mu_alpha_postMean",       V = "mu_alpha_postVar"),       
  list(yi = "mu_theta_postMean",       V = "mu_theta_postVar"),       
  list(yi = "logmu_sigma_postMean",    V = "logmu_sigma_postVar"),
  list(yi = "logshift_sigma_postMean", V = "logshift_sigma_postVar"),
  list(yi = "logsd_alpha_postMean",    V = "logsd_alpha_postVar"),  
  list(yi = "logsd_theta_postMean",    V = "logsd_theta_postVar"), 
  list(yi = "logsd_sigma_postMean",    V = "logsd_sigma_postVar"),
  list(yi = "logscale_sigma_postMean", V = "logscale_sigma_postVar")
)

# Estimate all meta analysis at once
gauss_meta_results <-  lapply(model_specs, function(spec) {
  args <- modifyList(common_args, list(yi = gaussian_estimates[[spec$yi]], 
                                       V = gaussian_estimates[[spec$V]]))
  do.call(rma.mv, args)
})

# Rename meta-analytic results
names(gauss_meta_results) <- c("mu_alpha", "mu_theta", "logmu_sigma", "logshift_sigma",
                               "logsd_alpha", "logsd_theta", "logsd_sigma", "logscale_sigma")

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 5: Meta-analysis of the exgaussian model (tau parameterization)
# Note: Fixed effects are means and random effects are (log) variances
# ─────────────────────────────────────────────────────────────────────────────

# Selected cols
selected_cols <- c("author", "task", "n_subj", "n_cong_trials", "n_incong_trials",
                   # Fixed effects and uncertainty
                   "mu_alpha_postMean","mu_alpha_postSd", 
                   "mu_theta_postMean", "mu_theta_postSd", 
                   "logmu_sigma_postMean", "logmu_sigma_postSd",
                   "logshift_sigma_postMean", "logshift_sigma_postSd",
                   "logmu_tau_postMean", "logmu_tau_postSd",
                   "logshift_tau_postMean", "logshift_tau_postSd",
                   # Random effects and uncertainty
                   "logsd_alpha_postMean","logsd_alpha_postSd", 
                   "logsd_theta_postMean", "logsd_theta_postSd", 
                   "logsd_sigma_postMean", "logsd_sigma_postSd",
                   "logscale_sigma_postMean", "logscale_sigma_postSd",
                   "logsd_tau_postMean", "logsd_tau_postSd",
                   "logscale_tau_postMean", "logscale_tau_postSd")

# Final gaussian data frame
exgaussian_estimates <- parameter_estimates$exgaussian_tau[,selected_cols]
exgaussian_estimates$author <- as.factor(sapply(strsplit(parameter_estimates$exgaussian_tau$author, "_"), function(x){
  paste0(x[1], " (", x[2], ")")
}))


# Change posterior sds to posterior variances
exgaussian_estimates$mu_alpha_postVar       <- exgaussian_estimates$mu_alpha_postSd^2
exgaussian_estimates$mu_theta_postVar       <- exgaussian_estimates$mu_theta_postSd^2
exgaussian_estimates$logmu_sigma_postVar    <- exgaussian_estimates$logmu_sigma_postSd^2
exgaussian_estimates$logshift_sigma_postVar <- exgaussian_estimates$logshift_sigma_postSd^2
exgaussian_estimates$logshift_tau_postVar   <- exgaussian_estimates$logshift_tau_postSd^2
exgaussian_estimates$logmu_tau_postVar      <- exgaussian_estimates$logmu_tau_postSd^2
exgaussian_estimates$logsd_alpha_postVar    <- exgaussian_estimates$logsd_alpha_postSd^2
exgaussian_estimates$logsd_theta_postVar    <- exgaussian_estimates$logsd_theta_postSd^2
exgaussian_estimates$logsd_sigma_postVar    <- exgaussian_estimates$logsd_sigma_postSd^2
exgaussian_estimates$logsd_tau_postVar      <- exgaussian_estimates$logsd_tau_postSd^2
exgaussian_estimates$logscale_sigma_postVar <- exgaussian_estimates$logscale_sigma_postSd^2
exgaussian_estimates$logscale_tau_postVar   <- exgaussian_estimates$logscale_tau_postSd^2


# Check normality assumption
ggqqplot(data = exgaussian_estimates, x = "mu_alpha_postMean")
ggqqplot(data = exgaussian_estimates, x = "mu_theta_postMean")
ggqqplot(data = exgaussian_estimates, x = "logmu_sigma_postMean")
ggqqplot(data = exgaussian_estimates, x = "logmu_tau_postMean")
ggqqplot(data = exgaussian_estimates, x = "logshift_sigma_postMean")
ggqqplot(data = exgaussian_estimates, x = "logshift_tau_postMean")
ggqqplot(data = exgaussian_estimates, x = "logsd_alpha_postMean")
ggqqplot(data = exgaussian_estimates, x = "logsd_theta_postMean")
ggqqplot(data = exgaussian_estimates, x = "logsd_sigma_postMean")
ggqqplot(data = exgaussian_estimates, x = "logsd_tau_postMean")
ggqqplot(data = exgaussian_estimates, x = "logscale_sigma_postMean")
ggqqplot(data = exgaussian_estimates, x = "logscale_tau_postMean")

# Common arguments
common_args <- list(random = ~ factor(task) - 1 | author,
                    mods = ~ factor(task) - 1,
                    data = exgaussian_estimates,
                    struct = "DIAG",
                    method = "REML",
                    test = "knha")

# Specify each model parameter and uncertainty
model_specs <- list(
  list(yi = "mu_alpha_postMean",       V = "mu_alpha_postVar"),       
  list(yi = "mu_theta_postMean",       V = "mu_theta_postVar"),       
  list(yi = "logmu_sigma_postMean",    V = "logmu_sigma_postVar"),    
  list(yi = "logshift_sigma_postMean", V = "logshift_sigma_postVar"),
  list(yi = "logmu_tau_postMean",      V = "logmu_tau_postVar"),  
  list(yi = "logshift_tau_postMean",   V = "logshift_tau_postVar"),
  list(yi = "logsd_alpha_postMean",    V = "logsd_alpha_postVar"),  
  list(yi = "logsd_theta_postMean",    V = "logsd_theta_postVar"), 
  list(yi = "logsd_sigma_postMean",    V = "logsd_sigma_postVar"),
  list(yi = "logscale_sigma_postMean", V = "logscale_sigma_postVar"),
  list(yi = "logsd_tau_postMean",      V = "logsd_tau_postVar"),
  list(yi = "logscale_tau_postMean",   V = "logscale_tau_postVar")
)

# Estimate all meta analysis at once
exgauss_meta_results_tau <-  lapply(model_specs, function(spec) {
  args <- modifyList(common_args, list(yi = exgaussian_estimates[[spec$yi]], 
                                       V = exgaussian_estimates[[spec$V]]))
  do.call(rma.mv, args)
})

names(exgauss_meta_results_tau) <- c("mu_alpha", "mu_theta", "logmu_sigma", "logshift_sigma",
                                     "logmu_tau", "logshift_tau", "logsd_alpha", "logsd_theta", 
                                     "logsd_sigma", "logscale_sigma", "logsd_tau", "logscale_tau")

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 6: Meta-analysis of the exgaussian model (rate parameterization)
# Note: Fixed effects are means and random effects are (log) variances
# ─────────────────────────────────────────────────────────────────────────────

# Selected cols
selected_cols <- c("author", "task", "n_subj", "n_cong_trials", "n_incong_trials",
                   # Fixed effects and uncertainty
                   "mu_alpha_postMean","mu_alpha_postSd", 
                   "mu_theta_postMean", "mu_theta_postSd", 
                   "logmu_sigma_postMean", "logmu_sigma_postSd",
                   "logshift_sigma_postMean", "logshift_sigma_postSd",
                   "logmu_lambda_postMean", "logmu_lambda_postSd",
                   "logshift_lambda_postMean", "logshift_lambda_postSd",
                   # Random effects and uncertainty
                   "logsd_alpha_postMean","logsd_alpha_postSd", 
                   "logsd_theta_postMean", "logsd_theta_postSd", 
                   "logsd_sigma_postMean", "logsd_sigma_postSd",
                   "logscale_sigma_postMean", "logscale_sigma_postSd",
                   "logsd_lambda_postMean", "logsd_lambda_postSd",
                   "logscale_lambda_postMean", "logscale_lambda_postSd")

# Final gaussian data frame
exgaussian_estimates <- parameter_estimates$exgaussian_rate[,selected_cols]
exgaussian_estimates$author <- as.factor(sapply(strsplit(parameter_estimates$exgaussian_rate$author, "_"), function(x){
  paste0(x[1], " (", x[2], ")")
}))


# Change posterior sds to posterior variances
exgaussian_estimates$mu_alpha_postVar       <- exgaussian_estimates$mu_alpha_postSd^2
exgaussian_estimates$mu_theta_postVar       <- exgaussian_estimates$mu_theta_postSd^2
exgaussian_estimates$logmu_sigma_postVar    <- exgaussian_estimates$logmu_sigma_postSd^2
exgaussian_estimates$logshift_sigma_postVar <- exgaussian_estimates$logshift_sigma_postSd^2
exgaussian_estimates$logshift_lambda_postVar   <- exgaussian_estimates$logshift_lambda_postSd^2
exgaussian_estimates$logmu_lambda_postVar      <- exgaussian_estimates$logmu_lambda_postSd^2
exgaussian_estimates$logsd_alpha_postVar    <- exgaussian_estimates$logsd_alpha_postSd^2
exgaussian_estimates$logsd_theta_postVar    <- exgaussian_estimates$logsd_theta_postSd^2
exgaussian_estimates$logsd_sigma_postVar    <- exgaussian_estimates$logsd_sigma_postSd^2
exgaussian_estimates$logsd_lambda_postVar      <- exgaussian_estimates$logsd_lambda_postSd^2
exgaussian_estimates$logscale_sigma_postVar <- exgaussian_estimates$logscale_sigma_postSd^2
exgaussian_estimates$logscale_lambda_postVar   <- exgaussian_estimates$logscale_lambda_postSd^2


# Check normality assumption
ggqqplot(data = exgaussian_estimates, x = "mu_alpha_postMean")
ggqqplot(data = exgaussian_estimates, x = "mu_theta_postMean")
ggqqplot(data = exgaussian_estimates, x = "logmu_sigma_postMean")
ggqqplot(data = exgaussian_estimates, x = "logmu_lambda_postMean")
ggqqplot(data = exgaussian_estimates, x = "logshift_sigma_postMean")
ggqqplot(data = exgaussian_estimates, x = "logshift_lambda_postMean")
ggqqplot(data = exgaussian_estimates, x = "logsd_alpha_postMean")
ggqqplot(data = exgaussian_estimates, x = "logsd_theta_postMean")
ggqqplot(data = exgaussian_estimates, x = "logsd_sigma_postMean")
ggqqplot(data = exgaussian_estimates, x = "logsd_lambda_postMean")
ggqqplot(data = exgaussian_estimates, x = "logscale_sigma_postMean")
ggqqplot(data = exgaussian_estimates, x = "logscale_lambda_postMean")

# Common arguments
common_args <- list(random = ~ factor(task) - 1 | author,
                    mods = ~ factor(task) - 1,
                    data = exgaussian_estimates,
                    struct = "DIAG",
                    method = "REML",
                    test = "knha")

# Specify each model parameter and uncertainty
model_specs <- list(
  list(yi = "mu_alpha_postMean",       V = "mu_alpha_postVar"),       
  list(yi = "mu_theta_postMean",       V = "mu_theta_postVar"),       
  list(yi = "logmu_sigma_postMean",    V = "logmu_sigma_postVar"),    
  list(yi = "logshift_sigma_postMean", V = "logshift_sigma_postVar"),
  list(yi = "logmu_lambda_postMean",      V = "logmu_lambda_postVar"),  
  list(yi = "logshift_lambda_postMean",   V = "logshift_lambda_postVar"),
  list(yi = "logsd_alpha_postMean",    V = "logsd_alpha_postVar"),  
  list(yi = "logsd_theta_postMean",    V = "logsd_theta_postVar"), 
  list(yi = "logsd_sigma_postMean",    V = "logsd_sigma_postVar"),
  list(yi = "logscale_sigma_postMean", V = "logscale_sigma_postVar"),
  list(yi = "logsd_lambda_postMean",      V = "logsd_lambda_postVar"),
  list(yi = "logscale_lambda_postMean",   V = "logscale_lambda_postVar")
)

# Estimate all meta analysis at once
exgauss_meta_results_rate <-  lapply(model_specs, function(spec) {
  args <- modifyList(common_args, list(yi = exgaussian_estimates[[spec$yi]], 
                                       V = exgaussian_estimates[[spec$V]]))
  do.call(rma.mv, args)
})

names(exgauss_meta_results_rate) <- c("mu_alpha", "mu_theta", "logmu_sigma", "logshift_sigma",
                                     "logmu_lambda", "logshift_lambda", "logsd_alpha", "logsd_theta", 
                                     "logsd_sigma", "logscale_sigma", "logsd_lambda", "logscale_lambda")

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 7: Meta-analysis of the shifted-lognormal model
# Note: Fixed effects are means and random effects are (log) variances
# ─────────────────────────────────────────────────────────────────────────────

# Selected cols
selected_cols <- c("author", "task", "n_subj", "n_cong_trials", "n_incong_trials",
                   # Fixed effects and uncertainty
                   "mu_alpha_postMean","mu_alpha_postSd", 
                   "mu_theta_postMean", "mu_theta_postSd", 
                   "logmu_sigma_postMean", "logmu_sigma_postSd",
                   "logshift_sigma_postMean", "logshift_sigma_postSd",
                   "mu_delta_postMean", "mu_delta_postSd",
                   # Random effects and uncertainty
                   "logsd_alpha_postMean", "logsd_alpha_postSd", 
                   "logsd_theta_postMean", "logsd_theta_postSd", 
                   "logsd_sigma_postMean", "logsd_sigma_postSd",
                   "logscale_sigma_postMean", "logscale_sigma_postSd",
                   "logsd_delta_postMean", "logsd_delta_postSd")

# Final gaussian data frame
shlognormal_estimates <- parameter_estimates$shifted.lognormal[,selected_cols]
shlognormal_estimates$author <- as.factor(sapply(strsplit(parameter_estimates$shifted.lognormal$author, "_"), function(x){
  paste0(x[1], " (", x[2], ")")
}))


# Change posterior sds to posterior variances
shlognormal_estimates$mu_alpha_postVar        <- shlognormal_estimates$mu_alpha_postSd^2
shlognormal_estimates$mu_theta_postVar        <- shlognormal_estimates$mu_theta_postSd^2
shlognormal_estimates$logmu_sigma_postVar     <- shlognormal_estimates$logmu_sigma_postSd^2
shlognormal_estimates$logshift_sigma_postVar  <- shlognormal_estimates$logshift_sigma_postSd^2
shlognormal_estimates$mu_delta_postVar        <- shlognormal_estimates$mu_delta_postSd^2
shlognormal_estimates$logsd_alpha_postVar     <- shlognormal_estimates$logsd_alpha_postSd^2
shlognormal_estimates$logsd_theta_postVar     <- shlognormal_estimates$logsd_theta_postSd^2
shlognormal_estimates$logsd_sigma_postVar     <- shlognormal_estimates$logsd_sigma_postSd^2
shlognormal_estimates$logscale_sigma_postVar  <- shlognormal_estimates$logscale_sigma_postSd^2
shlognormal_estimates$logsd_delta_postVar     <- shlognormal_estimates$logsd_delta_postSd^2

# Check normality assumption
ggqqplot(data = shlognormal_estimates, x = "mu_alpha_postMean")
ggqqplot(data = shlognormal_estimates, x = "mu_theta_postMean")
ggqqplot(data = shlognormal_estimates, x = "logmu_sigma_postMean")
ggqqplot(data = shlognormal_estimates, x = "logshift_sigma_postMean")
ggqqplot(data = shlognormal_estimates, x = "mu_delta_postMean")
ggqqplot(data = shlognormal_estimates, x = "logsd_alpha_postMean")
ggqqplot(data = shlognormal_estimates, x = "logsd_theta_postMean")
ggqqplot(data = shlognormal_estimates, x = "logsd_sigma_postMean")
ggqqplot(data = shlognormal_estimates, x = "logscale_sigma_postMean")
ggqqplot(data = shlognormal_estimates, x = "logsd_delta_postMean")

# Common arguments
common_args <- list(random = ~ factor(task) - 1 | author,
                    mods = ~ factor(task) - 1,
                    data = shlognormal_estimates,
                    struct = "DIAG",
                    method = "REML",
                    test = "knha")

# Specify each model parameter and uncertainty
model_specs <- list(
  list(yi = "mu_alpha_postMean",       V = "mu_alpha_postVar"),       
  list(yi = "mu_theta_postMean",       V = "mu_theta_postVar"),       
  list(yi = "logmu_sigma_postMean",    V = "logmu_sigma_postVar"),  
  list(yi = "logshift_sigma_postMean", V = "logshift_sigma_postVar"),
  list(yi = "mu_delta_postMean",       V = "mu_delta_postVar"),  
  list(yi = "logsd_alpha_postMean",    V = "logsd_alpha_postVar"),  
  list(yi = "logsd_theta_postMean",    V = "logsd_theta_postVar"), 
  list(yi = "logsd_sigma_postMean",    V = "logsd_sigma_postVar"),
  list(yi = "logscale_sigma_postMean", V = "logscale_sigma_postVar"),
  list(yi = "logsd_delta_postMean",    V = "logsd_delta_postVar") 
)

# Estimate all meta analysis at once
shlognormal_meta_results <-  lapply(model_specs, function(spec) {
  args <- modifyList(common_args, list(yi = shlognormal_estimates[[spec$yi]], 
                                       V = shlognormal_estimates[[spec$V]]))
  do.call(rma.mv, args)
})

names(shlognormal_meta_results) <- c("mu_alpha", "mu_theta", "logmu_sigma", "logshift_sigma",
                                     "mu_delta", "logsd_alpha", "logsd_theta", "logsd_sigma", 
                                     "logscale_sigma", "logsd_delta")

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 8: summaries and save meta-analytic results
# ─────────────────────────────────────────────────────────────────────────────

# Final data frame results per model: gaussian model
meta_parameters_gauss <- data.frame(
  task = c("flanker", "simon", "stroop"),
  mu_alpha    = coef(gauss_meta_results$mu_alpha),
  mu_theta    = coef(gauss_meta_results$mu_theta),
  mu_sigma    = exp(coef(gauss_meta_results$logmu_sigma)),
  sd_alpha    = exp(coef(gauss_meta_results$logsd_alpha)),
  sd_theta    = exp(coef(gauss_meta_results$logsd_theta)),
  sd_sigma    = exp(coef(gauss_meta_results$logsd_sigma)),
  shift_sigma = exp(coef(gauss_meta_results$logshift_sigma)),
  scale_sigma = exp(coef(gauss_meta_results$logscale_sigma)))

# Final data frame results per model: exgaussian model
meta_parameters_exgauss_tau <- data.frame(
  task = c("flanker", "simon", "stroop"),
  mu_alpha    = coef(exgauss_meta_results_tau$mu_alpha),
  mu_theta    = coef(exgauss_meta_results_tau$mu_theta),
  mu_sigma    = exp(coef(exgauss_meta_results_tau$logmu_sigma)),
  mu_tau      = exp(coef(exgauss_meta_results_tau$logmu_tau)),
  sd_alpha    = exp(coef(exgauss_meta_results_tau$logsd_alpha)),
  sd_theta    = exp(coef(exgauss_meta_results_tau$logsd_theta)),
  sd_tau      = exp(coef(exgauss_meta_results_tau$logsd_tau)),
  sd_sigma    = exp(coef(exgauss_meta_results_tau$logsd_sigma)),
  shift_sigma = exp(coef(exgauss_meta_results_tau$logshift_sigma)),
  shift_tau   = exp(coef(exgauss_meta_results_tau$logshift_tau)),
  scale_sigma = exp(coef(exgauss_meta_results_tau$logscale_sigma)),
  scale_tau   = exp(coef(exgauss_meta_results_tau$logscale_tau)))

# Final data frame results per model: exgaussian model
meta_parameters_exgauss_rate <- data.frame(
  task = c("flanker", "simon", "stroop"),
  mu_alpha     = coef(exgauss_meta_results_rate$mu_alpha),
  mu_theta     = coef(exgauss_meta_results_rate$mu_theta),
  mu_sigma     = exp(coef(exgauss_meta_results_rate$logmu_sigma)),
  mu_lambda    = exp(coef(exgauss_meta_results_rate$logmu_lambda)),
  sd_alpha     = exp(coef(exgauss_meta_results_rate$logsd_alpha)),
  sd_theta     = exp(coef(exgauss_meta_results_rate$logsd_theta)),
  sd_lambda    = exp(coef(exgauss_meta_results_rate$logsd_lambda)),
  sd_sigma     = exp(coef(exgauss_meta_results_rate$logsd_sigma)),
  shift_sigma  = exp(coef(exgauss_meta_results_rate$logshift_sigma)),
  shift_lambda = exp(coef(exgauss_meta_results_rate$logshift_lambda)),
  scale_sigma  = exp(coef(exgauss_meta_results_rate$logscale_sigma)),
  scale_lambda = exp(coef(exgauss_meta_results_rate$logscale_lambda)))

# Final data frame results per model: shifted-lognormal model
meta_parameters_shlognormal <- data.frame(
  task = c("flanker", "simon", "stroop"),
  mu_alpha    = coef(shlognormal_meta_results$mu_alpha),
  mu_theta    = coef(shlognormal_meta_results$mu_theta),
  mu_delta    = coef(shlognormal_meta_results$mu_delta),
  mu_sigma    = exp(coef(shlognormal_meta_results$logmu_sigma)),
  sd_alpha    = exp(coef(shlognormal_meta_results$logsd_alpha)),
  sd_theta    = exp(coef(shlognormal_meta_results$logsd_theta)),
  sd_delta    = exp(coef(shlognormal_meta_results$logsd_delta)),
  sd_sigma    = exp(coef(shlognormal_meta_results$logsd_sigma)),
  shift_sigma = exp(coef(shlognormal_meta_results$logshift_sigma)),
  scale_sigma = exp(coef(shlognormal_meta_results$logscale_sigma)))

# Remove row names and save it
rownames(meta_parameters_gauss) <- rownames(meta_parameters_exgauss_tau) <- 
  rownames(meta_parameters_exgauss_rate) <- rownames(meta_parameters_shlognormal) <- NULL

# Save it
meta_parameters <- list(meta_gaussian = meta_parameters_gauss,
                        meta_exgaussian_tau = meta_parameters_exgauss_tau,
                        meta_exgaussian_rate = meta_parameters_exgauss_rate,
                        meta_shlognormal = meta_parameters_shlognormal)

save(meta_parameters, file = "Results/Rdata/Meta-parameters/meta_parameters.rdata")

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 9: Summaries with 95% confidence intervals
# ─────────────────────────────────────────────────────────────────────────────

# Final data frame results per model: gaussian model
meta_parameters_full_gauss <- data.frame(
  task = c("flanker", "simon", "stroop"),
  mu_alpha        = coef(gauss_meta_results$mu_alpha),
  CI_mu_alpha     = CIs(gauss_meta_results$mu_alpha), 
  mu_theta        = coef(gauss_meta_results$mu_theta),
  CI_mu_theta     = CIs(gauss_meta_results$mu_theta), 
  mu_sigma        = exp(coef(gauss_meta_results$logmu_sigma)),
  CI_mu_sigma     = CIs(gauss_meta_results$logmu_sigma, log_param = TRUE), 
  sd_alpha        = exp(coef(gauss_meta_results$logsd_alpha)),
  CI_sd_alpha     = CIs(gauss_meta_results$logsd_alpha, log_param = TRUE), 
  sd_theta        = exp(coef(gauss_meta_results$logsd_theta)),
  CI_sd_theta     = CIs(gauss_meta_results$logsd_theta, log_param = TRUE), 
  sd_sigma        = exp(coef(gauss_meta_results$logsd_sigma)),
  CI_sd_sigma     = CIs(gauss_meta_results$logsd_sigma, log_param = TRUE), 
  shift_sigma     = exp(coef(gauss_meta_results$logshift_sigma)),
  CI_shift_sigma  = CIs(gauss_meta_results$logshift_sigma, log_param = TRUE), 
  scale_sigma     = exp(coef(gauss_meta_results$logscale_sigma)),
  CI_scale_sigma  = CIs(gauss_meta_results$logscale_sigma, log_param = TRUE))

# Ex-Gaussian (τ-based) model
meta_parameters_full_exgauss_tau <- data.frame(
  task           = c("flanker", "simon", "stroop"),
  mu_alpha       = coef(exgauss_meta_results_tau$mu_alpha),
  CI_mu_alpha    = CIs(exgauss_meta_results_tau$mu_alpha),
  mu_theta       = coef(exgauss_meta_results_tau$mu_theta),
  CI_mu_theta    = CIs(exgauss_meta_results_tau$mu_theta),
  mu_sigma       = exp(coef(exgauss_meta_results_tau$logmu_sigma)),
  CI_mu_sigma    = CIs(exgauss_meta_results_tau$logmu_sigma,    log_param = TRUE),
  mu_tau         = exp(coef(exgauss_meta_results_tau$logmu_tau)),
  CI_mu_tau      = CIs(exgauss_meta_results_tau$logmu_tau,      log_param = TRUE),
  sd_alpha       = exp(coef(exgauss_meta_results_tau$logsd_alpha)),
  CI_sd_alpha    = CIs(exgauss_meta_results_tau$logsd_alpha,    log_param = TRUE),
  sd_theta       = exp(coef(exgauss_meta_results_tau$logsd_theta)),
  CI_sd_theta    = CIs(exgauss_meta_results_tau$logsd_theta,    log_param = TRUE),
  sd_tau         = exp(coef(exgauss_meta_results_tau$logsd_tau)),
  CI_sd_tau      = CIs(exgauss_meta_results_tau$logsd_tau,      log_param = TRUE),
  sd_sigma       = exp(coef(exgauss_meta_results_tau$logsd_sigma)),
  CI_sd_sigma    = CIs(exgauss_meta_results_tau$logsd_sigma,    log_param = TRUE),
  shift_sigma    = exp(coef(exgauss_meta_results_tau$logshift_sigma)),
  CI_shift_sigma = CIs(exgauss_meta_results_tau$logshift_sigma, log_param = TRUE),
  shift_tau      = exp(coef(exgauss_meta_results_tau$logshift_tau)),
  CI_shift_tau   = CIs(exgauss_meta_results_tau$logshift_tau,   log_param = TRUE),
  scale_sigma    = exp(coef(exgauss_meta_results_tau$logscale_sigma)),
  CI_scale_sigma = CIs(exgauss_meta_results_tau$logscale_sigma, log_param = TRUE),
  scale_tau      = exp(coef(exgauss_meta_results_tau$logscale_tau)),
  CI_scale_tau   = CIs(exgauss_meta_results_tau$logscale_tau,   log_param = TRUE)
)


# Ex-Gaussian (rate-based) model
meta_parameters_full_exgauss_rate <- data.frame(
  task             = c("flanker", "simon", "stroop"),
  mu_alpha         = coef(exgauss_meta_results_rate$mu_alpha),
  CI_mu_alpha      = CIs(exgauss_meta_results_rate$mu_alpha),
  mu_theta         = coef(exgauss_meta_results_rate$mu_theta),
  CI_mu_theta      = CIs(exgauss_meta_results_rate$mu_theta),
  mu_sigma         = exp(coef(exgauss_meta_results_rate$logmu_sigma)),
  CI_mu_sigma      = CIs(exgauss_meta_results_rate$logmu_sigma,    log_param = TRUE),
  mu_lambda        = exp(coef(exgauss_meta_results_rate$logmu_lambda)),
  CI_mu_lambda     = CIs(exgauss_meta_results_rate$logmu_lambda,   log_param = TRUE),
  sd_alpha         = exp(coef(exgauss_meta_results_rate$logsd_alpha)),
  CI_sd_alpha      = CIs(exgauss_meta_results_rate$logsd_alpha,    log_param = TRUE),
  sd_theta         = exp(coef(exgauss_meta_results_rate$logsd_theta)),
  CI_sd_theta      = CIs(exgauss_meta_results_rate$logsd_theta,    log_param = TRUE),
  sd_lambda        = exp(coef(exgauss_meta_results_rate$logsd_lambda)),
  CI_sd_lambda     = CIs(exgauss_meta_results_rate$logsd_lambda,   log_param = TRUE),
  sd_sigma         = exp(coef(exgauss_meta_results_rate$logsd_sigma)),
  CI_sd_sigma      = CIs(exgauss_meta_results_rate$logsd_sigma,    log_param = TRUE),
  shift_sigma      = exp(coef(exgauss_meta_results_rate$logshift_sigma)),
  CI_shift_sigma   = CIs(exgauss_meta_results_rate$logshift_sigma, log_param = TRUE),
  shift_lambda     = exp(coef(exgauss_meta_results_rate$logshift_lambda)),
  CI_shift_lambda  = CIs(exgauss_meta_results_rate$logshift_lambda,log_param = TRUE),
  scale_sigma      = exp(coef(exgauss_meta_results_rate$logscale_sigma)),
  CI_scale_sigma   = CIs(exgauss_meta_results_rate$logscale_sigma, log_param = TRUE),
  scale_lambda     = exp(coef(exgauss_meta_results_rate$logscale_lambda)),
  CI_scale_lambda  = CIs(exgauss_meta_results_rate$logscale_lambda,log_param = TRUE)
)


# Shifted-Lognormal model
meta_parameters_full_shlognormal <- data.frame(
  task             = c("flanker", "simon", "stroop"),
  mu_alpha         = coef(shlognormal_meta_results$mu_alpha),
  CI_mu_alpha      = CIs(shlognormal_meta_results$mu_alpha),
  mu_theta         = coef(shlognormal_meta_results$mu_theta),
  CI_mu_theta      = CIs(shlognormal_meta_results$mu_theta),
  mu_delta         = coef(shlognormal_meta_results$mu_delta),
  CI_mu_delta      = CIs(shlognormal_meta_results$mu_delta),
  mu_sigma         = exp(coef(shlognormal_meta_results$logmu_sigma)),
  CI_mu_sigma      = CIs(shlognormal_meta_results$logmu_sigma,    log_param = TRUE),
  sd_alpha         = exp(coef(shlognormal_meta_results$logsd_alpha)),
  CI_sd_alpha      = CIs(shlognormal_meta_results$logsd_alpha,    log_param = TRUE),
  sd_theta         = exp(coef(shlognormal_meta_results$logsd_theta)),
  CI_sd_theta      = CIs(shlognormal_meta_results$logsd_theta,    log_param = TRUE),
  sd_delta         = exp(coef(shlognormal_meta_results$logsd_delta)),
  CI_sd_delta      = CIs(shlognormal_meta_results$logsd_delta,    log_param = TRUE),
  sd_sigma         = exp(coef(shlognormal_meta_results$logsd_sigma)),
  CI_sd_sigma      = CIs(shlognormal_meta_results$logsd_sigma,    log_param = TRUE),
  shift_sigma      = exp(coef(shlognormal_meta_results$logshift_sigma)),
  CI_shift_sigma   = CIs(shlognormal_meta_results$logshift_sigma, log_param = TRUE),
  scale_sigma      = exp(coef(shlognormal_meta_results$logscale_sigma)),
  CI_scale_sigma   = CIs(shlognormal_meta_results$logscale_sigma, log_param = TRUE)
)

# Save full meta-parameters
meta_parameters_full <- list(
  meta_parameters_full_gauss = meta_parameters_full_gauss,
  meta_parameters_full_exgauss_tau = meta_parameters_full_exgauss_tau,
  meta_parameters_full_exgauss_rate = meta_parameters_full_exgauss_rate,
  meta_parameters_full_shlognormal = meta_parameters_full_shlognormal
)

save(meta_parameters_full, file = "Results/Rdata/meta-parameters/meta_parameters_full.rdata")

# ─────────────────────────────────────────────────────────────────────────────

