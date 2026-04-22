# ╔════════════════════════════════════════════════════════════════════════════╗
# ║                             SCRIPT OVERVIEW                                ║
# ╠════════════════════════════════════════════════════════════════════════════╣
# ║ Script Name   : simulation_sensitivity_script.R                            ║
# ║ Author        : Ricardo Rey-Sáez                                           ║
# ║ Role          : PhD Student in Psychology                                  ║
# ║ Institution   : Autonomous University of Madrid, Spain                     ║
# ║ Email         : ricardoreysaez95@gmail.es                                  ║
# ║ Date          : 18-03-2026                                                 ║
# ╚════════════════════════════════════════════════════════════════════════════╝

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 1: Description
# ─────────────────────────────────────────────────────────────────────────────
# Run Sensitivity Simulation Study with Gaussian data-generating process
# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 2: Load Packages
# ─────────────────────────────────────────────────────────────────────────────

library(SimDesign)
library(MASS)
library(cmdstanr)
library(posterior)
library(loo)
library(future)
library(progressr)
library(dplyr)
library(tidyr)
library(stringr)

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 3: User-Defined Functions
# ─────────────────────────────────────────────────────────────────────────────

# All user-defined functions to simulate and estimate different models
source("R scripts/R functions/simulation_functions.R")

# Meta-analytic parameters
load("Results/Rdata/Meta-parameters/meta_parameters.rdata")

# ─────────────────────────────────────────────────────────────────────────────


# ─────────────────────────────────────────────────────────────────────────────
# SECTION 4: Simulation study workflow: Design
# ─────────────────────────────────────────────────────────────────────────────

# To expand: I = c(100, 200), M = c("uni", "2_EFA", "2_CFA")
Design <- createDesign(I = 100,
                       M = "uni",
                       reliab = c(.3, .5, .7),
                       std_lambda_j = c(.3, .5, .7))

# Fixed objects: all Stan models (confirmatory + exploratory) for all 3 families
fixed_objects <- list(C_Gauss     = cmdstan_model("Stan models/BGHFM (sim study)/C_BGHFM_gaussian.stan"),
                      E_Gauss     = cmdstan_model("Stan models/BGHFM (sim study)/E_BGHFM_gaussian.stan"),
                      C_Exgauss   = cmdstan_model("Stan models/BGHFM (sim study)/C_BGHFM_exgaussian_tau.stan"),
                      C_shlognorm = cmdstan_model("Stan models/BGHFM (sim study)/C_BGHFM_shlognormal.stan"),
                      meta_parameters = meta_parameters)

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 5: Simulation study workflow: Data generation function
# ─────────────────────────────────────────────────────────────────────────────

Generate <- function(condition, fixed_objects) {
  
  # Attach condition to simulate data
  Attach(condition)
  
  # Fixed values
  I <- 100; J <- 6; K <- 2; L <- 100
  
  # ─────────────────────────────────────── #
  #    Specify latent structure matrices    #
  # ─────────────────────────────────────── #
  
  # Unidimensional model only
  lambda_std <- matrix(rep(std_lambda_j, J), ncol = 1)
  Phi_cor <- diag(1)
  
  # ───────────────────────────────────────────────── #
  #    Select Gaussian population parameters          #
  # ───────────────────────────────────────────────── #
  
  # Gaussian meta-analytic parameters (true DGP)
  meta_pars <- "meta_gaussian"
  mu_alpha_j    = rep(fixed_objects$meta_parameters[[meta_pars]]$mu_alpha, length.out = J)
  mu_theta_j    = rep(fixed_objects$meta_parameters[[meta_pars]]$mu_theta, length.out = J)
  sd_alpha_j    = rep(fixed_objects$meta_parameters[[meta_pars]]$sd_alpha, length.out = J)
  shift_sigma_j = rep(fixed_objects$meta_parameters[[meta_pars]]$shift_sigma, length.out = J)
  scale_sigma_j = rep(fixed_objects$meta_parameters[[meta_pars]]$scale_sigma, length.out = J)
  
  # Global parameters
  reliab_j = rep(reliab, J)
  
  # Compute total sample size
  N <- I * J * K * L
  
  # ────────────────────────────────── #
  #    Prepare trial-level settings    #
  # ────────────────────────────────── #
  
  # Level-1 residual variance expected value and variance
  mu_sigma_j <- shift_sigma_j + scale_sigma_j * sqrt(2 / pi)
  sd_sigma_j <- scale_sigma_j * sqrt(1 - 2 / pi)
  
  # Simulate level-1 residual variance
  sigma <- matrix(NA, nrow = I, ncol = J)
  for(j in 1:J){
    # Half-normal values
    sigma_tilde <- truncnorm::rtruncnorm(n = I, a = 0, b = Inf,
                                         mean = 0, sd = 1)
    # Scaled values with correct E[X] and V[X]
    sigma[,j] <- shift_sigma_j[j] + scale_sigma_j[j] * sigma_tilde
  }
  
  # ─────────────────────────────────────── #
  #    Prepare individual-level settings    #
  # ─────────────────────────────────────── #
  
  # Model communality
  communality <- lambda_std %*% Phi_cor %*% t(lambda_std)
  
  # Population model-implied correlation matrix
  R_theta_j <- communality + diag(1 - diag(communality))
  
  # Signal-to-Noise ratio
  gamma2_j <- (2 * reliab_j) / (L * (1 - reliab_j))
  
  # Expected value of the trial-level variance under the Gaussian model
  # (no tau or delta variance components)
  Ev <- mu_sigma_j^2 + sd_sigma_j^2
  
  # Population standard deviation
  sd_theta_j <- sqrt(gamma2_j * Ev)
  
  # Population model-implied covariance matrix
  Sigma_theta_j <- diag(sd_theta_j) %*% R_theta_j %*% diag(sd_theta_j)
  
  # Random intercepts and slopes
  alpha <- mvrnorm(n = I, mu = mu_alpha_j, Sigma = diag(sd_alpha_j^2))
  theta <- mvrnorm(n = I, mu = mu_theta_j, Sigma = Sigma_theta_j)
  
  # ──────────────────────────────── #
  #    Simulate hierarchical data    #
  # ──────────────────────────────── #
  
  # Empty matrix to store simulated data
  dat <- matrix(NA, nrow = N, ncol = 5)
  
  # Simulate data
  index = 0
  for(i in 1:I){ # For each subject
    for(j in 1:J){ # In each task
      for(k in 1:K){ # In each condition
        for(l in 1:L){ # In each trial
          index = index + 1
          dat[index,1] <- i # Subject ID
          dat[index,2] <- j # Task ID
          dat[index,3] <- k # Condition ID
          dat[index,4] <- l # Trial ID
          
          # Level 2 linear prediction
          mu = alpha[i,j] + theta[i,j] * (k - 1.5)
          # Level 1 observation: Gaussian (TRUE DGP, nly positive RTs)
          dat[index, 5] <- truncnorm::rtruncnorm(1, a = 0, b = Inf, mean = mu + .2, sd = sigma[i, j])        }
      }
    }
  }
  
  # Add variable names
  colnames(dat) <- c('subject','task','condition','trial','RT')
  dat <- as.data.frame(dat)
  # Return dataset
  dat
}

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 6: Simulation study workflow: Data analysis
# ─────────────────────────────────────────────────────────────────────────────

Analyse <- function(condition, dat, fixed_objects) {
  
  # ────────────────────────────────────── #
  #    Estimation via method of moments    #
  # ────────────────────────────────────── #
  source("R scripts/R functions/simulation_functions.R")
  
  # Method-of-moments effects, correlations and SEs
  effects_mm <- analytical.effect(data = dat)
  Rho_mm_est <- cor(effects_mm)
  Rho_mm_se <- (1 - Rho_mm_est^2)/sqrt(nrow(effects_mm) - 3)
  Rho_mm_LL <- Rho_mm_est + qnorm(0.025) * Rho_mm_se
  Rho_mm_UL <- Rho_mm_est + qnorm(0.975) * Rho_mm_se
  Rho_mm <- c(Rho_mm_est = Rho_mm_est[lower.tri(Rho_mm_est)],
              Rho_mm_se = Rho_mm_se[lower.tri(Rho_mm_se)],
              Rho_mm_LL = Rho_mm_LL[lower.tri(Rho_mm_LL)],
              Rho_mm_UL = Rho_mm_UL[lower.tri(Rho_mm_UL)])
  
  # Spearman adjustment and SEs via lavaan
  sp_adj <- spearman.adj(data = dat)
  Rho_sp_est <- sp_adj$rho_est
  Rho_sp_se <- Rho_mm_se / sqrt(tcrossprod(sp_adj$reliability))
  Rho_sp_LL <- Rho_sp_est + qnorm(0.025) * Rho_sp_se
  Rho_sp_UL <- Rho_sp_est + qnorm(0.975) * Rho_sp_se
  Rho_sp <- c(Rho_sp_est = Rho_sp_est[lower.tri(Rho_sp_est)],
              Rho_sp_se  = Rho_sp_se[lower.tri(Rho_sp_se)],
              Rho_sp_LL  = Rho_sp_LL[lower.tri(Rho_sp_LL)],
              Rho_sp_UL  = Rho_sp_UL[lower.tri(Rho_sp_UL)])
  if(any(is.na(Rho_sp))) {stop("NA in the spearman adjusted correlation estimation")}
  
  # ─────────────────────── #
  #    Prepare Stan data    #
  # ─────────────────────── #
  
  # Prepare Stan data
  stan.data <- Stan.list.data(dat)
  
  # Unidimensional settings (always)
  stan.data$M <- 1
  stan.data$Lambda_index <- matrix(rep(1,6), ncol=1)
  
  # Define starting values for each model
  gauss_inits <- list(
    make_inits(df = dat, M = 1, model = "gaussian"), 
    make_inits(df = dat, M = 1, model = "gaussian")
  )
  
  exgauss_inits <- list(
    make_inits(df = dat, M = 1, model = "exgaussian"), 
    make_inits(df = dat, M = 1, model = "exgaussian")
  )
  
  shlognorm_inits <- list(
    make_inits(df = dat, M = 1, model = "shifted.lognormal"), 
    make_inits(df = dat, M = 1, model = "shifted.lognormal")
  )
  
  # Stan model arguments
  stan_model_arguments <- list(
    data <- stan.data,
    iter_sampling = 1000,
    iter_warmup = 500,
    chains = 2,
    parallel_chains = 2,
    adapt_delta = .8,
    refresh = 50,
    init = NULL,
    show_messages = FALSE
  )
  
  # ────────────────────── #
  #    Model estimation    #
  # ────────────────────── #
  
  # Gaussian model
  stan_model_arguments$init <- gauss_inits
  fit_gaussian <- do.call(fixed_objects$C_Gauss$sample, stan_model_arguments)
  
  # Ex-Gaussian model
  stan_model_arguments$init <- exgauss_inits
  fit_exgaussian <- do.call(fixed_objects$C_Exgauss$sample, stan_model_arguments)
  
  # Shifted-Lognormal model
  stan_model_arguments$init <- shlognorm_inits
  fit_shlognormal <- do.call(fixed_objects$C_shlognorm$sample, stan_model_arguments)
  
  # ───────────────────────── #
  #    Parameter estimates    #
  # ───────────────────────── #
  
  # Group-level parameters from the Gaussian model (true DGP)
  mu_alpha <- posterior_summaries(fit_gaussian, "mu_alpha")
  mu_theta <- posterior_summaries(fit_gaussian, "mu_theta")
  mu_sigma <- posterior_summaries(fit_gaussian, "mu_sigma")
  sd_alpha <- posterior_summaries(fit_gaussian, "sd_alpha")
  sd_theta <- posterior_summaries(fit_gaussian, "sd_theta")
  sd_sigma <- posterior_summaries(fit_gaussian, "sd_sigma")
  
  # Model-implied correlations per model
  Rho_gaussian    <- posterior_summaries(fit_gaussian,    "Rho_theta", cormat = TRUE, parnames = "Rho_gaussian")
  Rho_exgaussian  <- posterior_summaries(fit_exgaussian,  "Rho_theta", cormat = TRUE, parnames = "Rho_exgaussian")
  Rho_shlognormal <- posterior_summaries(fit_shlognormal, "Rho_theta", cormat = TRUE, parnames = "Rho_shlognormal")
  
  # Reliability per model
  reliab_gaussian    <- posterior_summaries(fit_gaussian,    "reliability", parnames = "reliab_gaussian")
  reliab_exgaussian  <- posterior_summaries(fit_exgaussian,  "reliability", parnames = "reliab_exgaussian")
  reliab_shlognormal <- posterior_summaries(fit_shlognormal, "reliability", parnames = "reliab_shlognormal")
  
  # Lambda (factor loadings) per model — unidimensional
  lambda_gaussian    <- posterior_summaries(fit_gaussian,    parameter = "Lambda_std", parnames = "lambda_gaussian")
  lambda_exgaussian  <- posterior_summaries(fit_exgaussian,  parameter = "Lambda_std", parnames = "lambda_exgaussian")
  lambda_shlognormal <- posterior_summaries(fit_shlognormal, parameter = "Lambda_std", parnames = "lambda_shlognormal")
  
  # Tucker index of factor congruence per model
  true_lambda <- matrix(rep(condition$std_lambda_j, 6), ncol = 1)
  
  Tucker_phi_gaussian <- diag(
    psych::factor.congruence(
      x = true_lambda, 
      y = matrix(fit_gaussian$summary("Lambda_std")$mean, nrow = 6, ncol = 1)
    )
  )
  
  Tucker_phi_exgaussian <- diag(
    psych::factor.congruence(
      x = true_lambda, 
      y = matrix(fit_exgaussian$summary("Lambda_std")$mean, nrow = 6, ncol = 1)
    )
  )
  
  Tucker_phi_shlognormal <- diag(
    psych::factor.congruence(
      x = true_lambda, 
      y = matrix(fit_shlognormal$summary("Lambda_std")$mean, nrow = 6, ncol = 1)
    )
  )
  
  # ──────────────────────── #
  #    LOO model comparison  #
  # ──────────────────────── #
  
  # Compute PSIS-LOO for all three models
  loo_gaussian    <- loo_BHGFM(fit_gaussian,    model = "gaussian",           data = dat, n.cores = 2)
  loo_exgaussian  <- loo_BHGFM(fit_exgaussian,  model = "exgaussian",         data = dat, n.cores = 2)
  loo_shlognormal <- loo_BHGFM(fit_shlognormal, model = "shifted.lognormal",  data = dat, n.cores = 2)
  
  # LOO model comparison: Gaussian vs Ex-Gaussian
  loo_comp_exg <- loo_compare(list(exgaussian = loo_exgaussian, gaussian = loo_gaussian))
  gauss_first_exg <- rownames(loo_comp_exg)[1] == "gaussian"
  loo_null_weak_CI_exg <- loo_comp_exg[2,1] + c(qnorm(p = .025), qnorm(p = .975)) * loo_comp_exg[2,2]
  loo_null_strong_CI_exg <- loo_comp_exg[2,1] + c(-4, 4) * loo_comp_exg[2,2]
  loo_weak_diff_exg <- !(loo_null_weak_CI_exg[1] <= 0 & 0 <= loo_null_weak_CI_exg[2])
  loo_strong_diff_exg <- !(loo_null_strong_CI_exg[1] <= 0 & 0 <= loo_null_strong_CI_exg[2])
  gauss_weak_diff_exg <- all(gauss_first_exg, loo_weak_diff_exg)
  gauss_strong_diff_exg <- all(gauss_first_exg, loo_strong_diff_exg)
  
  # LOO model comparison: Gaussian vs Shifted-Lognormal
  loo_comp_shl <- loo_compare(list(shlognormal = loo_shlognormal, gaussian = loo_gaussian))
  gauss_first_shl <- rownames(loo_comp_shl)[1] == "gaussian"
  loo_null_weak_CI_shl <- loo_comp_shl[2,1] + c(qnorm(p = .025), qnorm(p = .975)) * loo_comp_shl[2,2]
  loo_null_strong_CI_shl <- loo_comp_shl[2,1] + c(-4, 4) * loo_comp_shl[2,2]
  loo_weak_diff_shl <- !(loo_null_weak_CI_shl[1] <= 0 & 0 <= loo_null_weak_CI_shl[2])
  loo_strong_diff_shl <- !(loo_null_strong_CI_shl[1] <= 0 & 0 <= loo_null_strong_CI_shl[2])
  gauss_weak_diff_shl <- all(gauss_first_shl, loo_weak_diff_shl)
  gauss_strong_diff_shl <- all(gauss_first_shl, loo_strong_diff_shl)
  
  # LOO 3-way comparison (all three models)
  loo_comp_all <- loo_compare(list(gaussian = loo_gaussian, 
                                   exgaussian = loo_exgaussian, 
                                   shlognormal = loo_shlognormal))
  gauss_best <- rownames(loo_comp_all)[1] == "gaussian"
  
  # ──────────────────────── #
  #    Simulation results    #    
  # ──────────────────────── #
  
  ret <- c(
    # Group-level parameters (from Gaussian = true DGP)
    unlist(mu_alpha), unlist(mu_theta), unlist(mu_sigma),
    unlist(sd_alpha), unlist(sd_theta), unlist(sd_sigma),
    # Factor loadings per model
    unlist(lambda_gaussian), unlist(lambda_exgaussian), unlist(lambda_shlognormal),
    # Correlations per model
    unlist(Rho_gaussian), unlist(Rho_exgaussian), unlist(Rho_shlognormal),
    # Reliability per model
    unlist(reliab_gaussian), unlist(reliab_exgaussian), unlist(reliab_shlognormal),
    # Method-of-moments correlations
    unlist(Rho_mm), unlist(Rho_sp),
    # Tucker congruence
    Tucker_gaussian_f1    = Tucker_phi_gaussian,
    Tucker_exgaussian_f1  = Tucker_phi_exgaussian,
    Tucker_shlognormal_f1 = Tucker_phi_shlognormal,
    # LOO: Gaussian vs Ex-Gaussian
    gauss_first_exg          = gauss_first_exg,
    loo_elpd_diff_exg        = loo_comp_exg[2,1],
    loo_diff_weak_exg        = loo_weak_diff_exg,
    loo_diff_strong_exg      = loo_strong_diff_exg,
    loo_diff_gauss_weak_exg  = gauss_weak_diff_exg,
    loo_diff_gauss_strong_exg = gauss_strong_diff_exg,
    # LOO: Gaussian vs Shifted-Lognormal
    gauss_first_shl          = gauss_first_shl,
    loo_elpd_diff_shl        = loo_comp_shl[2,1],
    loo_diff_weak_shl        = loo_weak_diff_shl,
    loo_diff_strong_shl      = loo_strong_diff_shl,
    loo_diff_gauss_weak_shl  = gauss_weak_diff_shl,
    loo_diff_gauss_strong_shl = gauss_strong_diff_shl,
    # LOO: 3-way best model
    gauss_best = gauss_best,
    # Computation time
    gaussian_time    = fit_gaussian$time()$total, 
    exgaussian_time  = fit_exgaussian$time()$total, 
    shlognormal_time = fit_shlognormal$time()$total,
    # Divergent transitions
    gaussian_divergence    = fit_gaussian$diagnostic_summary()$num_divergent,
    exgaussian_divergence  = fit_exgaussian$diagnostic_summary()$num_divergent,
    shlognormal_divergence = fit_shlognormal$diagnostic_summary()$num_divergent,
    # Max treedepth
    gaussian_treedepth    = fit_gaussian$diagnostic_summary()$num_max_treedepth,
    exgaussian_treedepth  = fit_exgaussian$diagnostic_summary()$num_max_treedepth,
    shlognormal_treedepth = fit_shlognormal$diagnostic_summary()$num_max_treedepth,
    # EBFMI
    gaussian_ebfmi    = fit_gaussian$diagnostic_summary()$ebfmi,
    exgaussian_ebfmi  = fit_exgaussian$diagnostic_summary()$ebfmi,
    shlognormal_ebfmi = fit_shlognormal$diagnostic_summary()$ebfmi
  )
  
  ret
}

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 7: Simulation study workflow: Final summaries
# ─────────────────────────────────────────────────────────────────────────────

Summarise <- function(condition, results, fixed_objects) {
  
  # Load utils functions
  source("R scripts/R functions/simulation_functions.R")
  
  # ───────────────────────────────────────── #
  #    Specify population parameter values    #
  # ───────────────────────────────────────── #
  
  # Number of trials, tasks and conditions
  I <- 100; J <- 6; K <- 2; L <- 100
  
  # Unidimensional latent structure (always)
  lambda_std <- matrix(rep(condition$std_lambda_j, J), ncol = 1)
  Phi_cor <- diag(1)
  
  # Gaussian meta-analytic parameters (true DGP)
  meta_pars <- "meta_gaussian"
  mu_alpha_j    = rep(fixed_objects$meta_parameters[[meta_pars]]$mu_alpha, length.out = J)
  mu_theta_j    = rep(fixed_objects$meta_parameters[[meta_pars]]$mu_theta, length.out = J)
  sd_alpha_j    = rep(fixed_objects$meta_parameters[[meta_pars]]$sd_alpha, length.out = J)
  shift_sigma_j = rep(fixed_objects$meta_parameters[[meta_pars]]$shift_sigma, length.out = J)
  scale_sigma_j = rep(fixed_objects$meta_parameters[[meta_pars]]$scale_sigma, length.out = J)
  
  # Global parameters
  reliab_j = rep(condition$reliab, J)
  
  # Level-1 residual variance expected value and variance
  mu_sigma_j <- shift_sigma_j + scale_sigma_j * sqrt(2 / pi)
  sd_sigma_j <- scale_sigma_j * sqrt(1 - 2 / pi)
  
  # Model communality
  communality <- lambda_std %*% Phi_cor %*% t(lambda_std)
  
  # Population model-implied correlation matrix
  R_theta_j <- communality + diag(1 - diag(communality))
  
  # Signal-to-Noise ratio
  gamma2_j <- (2 * reliab_j) / (L * (1 - reliab_j))
  
  # Expected value of the trial-level variance under the Gaussian model
  Ev <- mu_sigma_j^2 + sd_sigma_j^2
  
  # Population standard deviation
  sd_theta_j <- sqrt(gamma2_j * Ev)
  
  # ───────────────────────────────────── #
  #    Bias in all parameter estimates    #
  # ───────────────────────────────────── #
  
  # Global summaries for group-level parameters (from Gaussian/true model)
  mu_alpha <- global_summaries(results = results, parameter_name = "mu_alpha", parameter = mu_alpha_j)
  mu_theta <- global_summaries(results = results, parameter_name = "mu_theta", parameter = mu_theta_j)
  mu_sigma <- global_summaries(results = results, parameter_name = "mu_sigma", parameter = mu_sigma_j)
  sd_alpha <- global_summaries(results = results, parameter_name = "sd_alpha", parameter = sd_alpha_j)
  sd_theta <- global_summaries(results = results, parameter_name = "sd_theta", parameter = sd_theta_j)
  sd_sigma <- global_summaries(results = results, parameter_name = "sd_sigma", parameter = sd_sigma_j)
  
  # Correlations per model
  Rho_gaussian    <- global_summaries(results = results, parameter_name = "Rho_gaussian",    parameter = R_theta_j, Rho = TRUE)
  Rho_exgaussian  <- global_summaries(results = results, parameter_name = "Rho_exgaussian",  parameter = R_theta_j, Rho = TRUE)
  Rho_shlognormal <- global_summaries(results = results, parameter_name = "Rho_shlognormal", parameter = R_theta_j, Rho = TRUE)
  Rho_mm <- global_summaries(results = results, parameter_name = "Rho_mm", parameter = R_theta_j, method_of_moments = TRUE)
  Rho_sp <- global_summaries(results = results, parameter_name = "Rho_sp", parameter = R_theta_j, method_of_moments = TRUE)
  
  # Reliability per model
  reliab_gaussian    <- global_summaries(results = results, parameter_name = "reliab_gaussian",    parameter = reliab_j)
  reliab_exgaussian  <- global_summaries(results = results, parameter_name = "reliab_exgaussian",  parameter = reliab_j)
  reliab_shlognormal <- global_summaries(results = results, parameter_name = "reliab_shlognormal", parameter = reliab_j)
  
  # Factor loadings per model (unidimensional)
  lambda_gaussian    <- global_summaries(results = results, parameter_name = "lambda_gaussian",    parameter = lambda_std[,1])
  lambda_exgaussian  <- global_summaries(results = results, parameter_name = "lambda_exgaussian",  parameter = lambda_std[,1])
  lambda_shlognormal <- global_summaries(results = results, parameter_name = "lambda_shlognormal", parameter = lambda_std[,1])
  
  # ──────────────────────────────────────── #
  #    LOO model comparison summaries        #
  # ──────────────────────────────────────── #
  
  # LOO: Gaussian vs Ex-Gaussian
  loo_avg_elpd_exg  <- mean(results[, "loo_elpd_diff_exg"])
  loo_disp_elpd_exg <- sd(results[, "loo_elpd_diff_exg"])
  loo_avg_weak_exg  <- mean(results[, "loo_diff_weak_exg"])
  loo_avg_strong_exg <- mean(results[, "loo_diff_strong_exg"])
  loo_avg_gauss_weak_exg  <- mean(results[, "loo_diff_gauss_weak_exg"])
  loo_avg_gauss_strong_exg <- mean(results[, "loo_diff_gauss_strong_exg"])
  gauss_first_exg_avg <- mean(results[, "gauss_first_exg"])
  
  # LOO: Gaussian vs Shifted-Lognormal
  loo_avg_elpd_shl  <- mean(results[, "loo_elpd_diff_shl"])
  loo_disp_elpd_shl <- sd(results[, "loo_elpd_diff_shl"])
  loo_avg_weak_shl  <- mean(results[, "loo_diff_weak_shl"])
  loo_avg_strong_shl <- mean(results[, "loo_diff_strong_shl"])
  loo_avg_gauss_weak_shl  <- mean(results[, "loo_diff_gauss_weak_shl"])
  loo_avg_gauss_strong_shl <- mean(results[, "loo_diff_gauss_strong_shl"])
  gauss_first_shl_avg <- mean(results[, "gauss_first_shl"])
  
  # LOO: 3-way best model
  gauss_best_avg <- mean(results[, "gauss_best"])
  
  # ────────────────────────────────── #
  #    Diagnostic summary statistics   #
  # ────────────────────────────────── #
  
  # EBFMI
  ebfmi_gaussian_avg     <- mean(colMeans(results[, grep("^gaussian_ebfmi", colnames(results)), drop = FALSE]))
  ebfmi_gaussian_disp    <- mean(apply(results[, grep("^gaussian_ebfmi", colnames(results)), drop = FALSE], 2, sd))
  ebfmi_exgaussian_avg   <- mean(colMeans(results[, grep("^exgaussian_ebfmi", colnames(results)), drop = FALSE]))
  ebfmi_exgaussian_disp  <- mean(apply(results[, grep("^exgaussian_ebfmi", colnames(results)), drop = FALSE], 2, sd))
  ebfmi_shlognormal_avg  <- mean(colMeans(results[, grep("^shlognormal_ebfmi", colnames(results)), drop = FALSE]))
  ebfmi_shlognormal_disp <- mean(apply(results[, grep("^shlognormal_ebfmi", colnames(results)), drop = FALSE], 2, sd))
  
  # Time
  time_gaussian_avg     <- mean(results[, grep("^gaussian_time", colnames(results))])
  time_gaussian_disp    <- sd(results[, grep("^gaussian_time", colnames(results))])
  time_exgaussian_avg   <- mean(results[, grep("^exgaussian_time", colnames(results))])
  time_exgaussian_disp  <- sd(results[, grep("^exgaussian_time", colnames(results))])
  time_shlognormal_avg  <- mean(results[, grep("^shlognormal_time", colnames(results))])
  time_shlognormal_disp <- sd(results[, grep("^shlognormal_time", colnames(results))])
  
  # Divergent transitions
  gaussian_divergence_n    <- sum(rowSums(results[, grep("^gaussian_divergence", colnames(results)), drop = FALSE]) > 0)
  exgaussian_divergence_n  <- sum(rowSums(results[, grep("^exgaussian_divergence", colnames(results)), drop = FALSE]) > 0)
  shlognormal_divergence_n <- sum(rowSums(results[, grep("^shlognormal_divergence", colnames(results)), drop = FALSE]) > 0)
  gaussian_divergence_avg     <- mean(rowSums(results[, grep("^gaussian_divergence", colnames(results)), drop = FALSE]) / 3000 * 100)
  gaussian_divergence_disp    <- sd(rowSums(results[, grep("^gaussian_divergence", colnames(results)), drop = FALSE]) / 3000 * 100)
  exgaussian_divergence_avg   <- mean(rowSums(results[, grep("^exgaussian_divergence", colnames(results)), drop = FALSE]) / 3000 * 100)
  exgaussian_divergence_disp  <- sd(rowSums(results[, grep("^exgaussian_divergence", colnames(results)), drop = FALSE]) / 3000 * 100)
  shlognormal_divergence_avg  <- mean(rowSums(results[, grep("^shlognormal_divergence", colnames(results)), drop = FALSE]) / 3000 * 100)
  shlognormal_divergence_disp <- sd(rowSums(results[, grep("^shlognormal_divergence", colnames(results)), drop = FALSE]) / 3000 * 100)
  
  # Maximum treedepth
  gaussian_treedepth_n    <- sum(rowSums(results[, grep("^gaussian_treedepth", colnames(results)), drop = FALSE]) > 0)
  exgaussian_treedepth_n  <- sum(rowSums(results[, grep("^exgaussian_treedepth", colnames(results)), drop = FALSE]) > 0)
  shlognormal_treedepth_n <- sum(rowSums(results[, grep("^shlognormal_treedepth", colnames(results)), drop = FALSE]) > 0)
  gaussian_treedepth_avg     <- mean(rowSums(results[, grep("^gaussian_treedepth", colnames(results)), drop = FALSE]) / 3000 * 100)
  gaussian_treedepth_disp    <- sd(rowSums(results[, grep("^gaussian_treedepth", colnames(results)), drop = FALSE]) / 3000 * 100)
  exgaussian_treedepth_avg   <- mean(rowSums(results[, grep("^exgaussian_treedepth", colnames(results)), drop = FALSE]) / 3000 * 100)
  exgaussian_treedepth_disp  <- sd(rowSums(results[, grep("^exgaussian_treedepth", colnames(results)), drop = FALSE]) / 3000 * 100)
  shlognormal_treedepth_avg  <- mean(rowSums(results[, grep("^shlognormal_treedepth", colnames(results)), drop = FALSE]) / 3000 * 100)
  shlognormal_treedepth_disp <- sd(rowSums(results[, grep("^shlognormal_treedepth", colnames(results)), drop = FALSE]) / 3000 * 100)
  
  # Tucker congruence index
  Tucker_gaussian_f1_avg     <- mean(results[, "Tucker_gaussian_f1"])
  Tucker_exgaussian_f1_avg   <- mean(results[, "Tucker_exgaussian_f1"])
  Tucker_shlognormal_f1_avg  <- mean(results[, "Tucker_shlognormal_f1"])
  Tucker_gaussian_f1_disp    <- sd(results[, "Tucker_gaussian_f1"])
  Tucker_exgaussian_f1_disp  <- sd(results[, "Tucker_exgaussian_f1"])
  Tucker_shlognormal_f1_disp <- sd(results[, "Tucker_shlognormal_f1"])
  
  # ──────────────────────── #
  #    Final results         #
  # ──────────────────────── #
  
  res <- c(
    # Correlations per model
    Rho_gaussian, Rho_exgaussian, Rho_shlognormal,
    # Method-of-moments correlations
    Rho_mm, Rho_sp,
    # Factor loadings per model
    lambda_gaussian, lambda_exgaussian, lambda_shlognormal,
    # Reliability per model
    reliab_gaussian, reliab_exgaussian, reliab_shlognormal,
    # Group-level parameters (true DGP model)
    mu_alpha, mu_theta, mu_sigma,
    sd_alpha, sd_theta, sd_sigma,
    # Tucker congruence
    Tucker_gaussian_f1_avg       = Tucker_gaussian_f1_avg,
    Tucker_exgaussian_f1_avg     = Tucker_exgaussian_f1_avg,
    Tucker_shlognormal_f1_avg    = Tucker_shlognormal_f1_avg,
    Tucker_gaussian_f1_disp      = Tucker_gaussian_f1_disp,
    Tucker_exgaussian_f1_disp    = Tucker_exgaussian_f1_disp,
    Tucker_shlognormal_f1_disp   = Tucker_shlognormal_f1_disp,
    # LOO: Gaussian vs Ex-Gaussian
    gauss_first_exg_avg          = gauss_first_exg_avg,
    loo_avg_elpd_exg             = loo_avg_elpd_exg,
    loo_disp_elpd_exg            = loo_disp_elpd_exg,
    loo_avg_weak_exg             = loo_avg_weak_exg,
    loo_avg_strong_exg           = loo_avg_strong_exg,
    loo_avg_gauss_weak_exg       = loo_avg_gauss_weak_exg,
    loo_avg_gauss_strong_exg     = loo_avg_gauss_strong_exg,
    # LOO: Gaussian vs Shifted-Lognormal
    gauss_first_shl_avg          = gauss_first_shl_avg,
    loo_avg_elpd_shl             = loo_avg_elpd_shl,
    loo_disp_elpd_shl            = loo_disp_elpd_shl,
    loo_avg_weak_shl             = loo_avg_weak_shl,
    loo_avg_strong_shl           = loo_avg_strong_shl,
    loo_avg_gauss_weak_shl       = loo_avg_gauss_weak_shl,
    loo_avg_gauss_strong_shl     = loo_avg_gauss_strong_shl,
    # LOO: 3-way
    gauss_best_avg               = gauss_best_avg,
    # Diagnostics: Divergences
    divergence_gaussian_n        = gaussian_divergence_n,
    divergence_exgaussian_n      = exgaussian_divergence_n,
    divergence_shlognormal_n     = shlognormal_divergence_n,
    divergence_gaussian_avg      = gaussian_divergence_avg,
    divergence_gaussian_disp     = gaussian_divergence_disp,
    divergence_exgaussian_avg    = exgaussian_divergence_avg,
    divergence_exgaussian_disp   = exgaussian_divergence_disp,
    divergence_shlognormal_avg   = shlognormal_divergence_avg,
    divergence_shlognormal_disp  = shlognormal_divergence_disp,
    # Diagnostics: Treedepth
    treedepth_gaussian_n         = gaussian_treedepth_n,
    treedepth_exgaussian_n       = exgaussian_treedepth_n,
    treedepth_shlognormal_n      = shlognormal_treedepth_n,
    treedepth_gaussian_avg       = gaussian_treedepth_avg,
    treedepth_gaussian_disp      = gaussian_treedepth_disp,
    treedepth_exgaussian_avg     = exgaussian_treedepth_avg,
    treedepth_exgaussian_disp    = exgaussian_treedepth_disp,
    treedepth_shlognormal_avg    = shlognormal_treedepth_avg,
    treedepth_shlognormal_disp   = shlognormal_treedepth_disp,
    # Diagnostics: EBFMI
    ebfmi_gaussian_avg           = ebfmi_gaussian_avg,
    ebfmi_gaussian_disp          = ebfmi_gaussian_disp,
    ebfmi_exgaussian_avg         = ebfmi_exgaussian_avg,
    ebfmi_exgaussian_disp        = ebfmi_exgaussian_disp,
    ebfmi_shlognormal_avg        = ebfmi_shlognormal_avg,
    ebfmi_shlognormal_disp       = ebfmi_shlognormal_disp,
    # Diagnostics: Time
    time_gaussian_avg            = time_gaussian_avg,
    time_gaussian_disp           = time_gaussian_disp,
    time_exgaussian_avg          = time_exgaussian_avg,
    time_exgaussian_disp         = time_exgaussian_disp,
    time_shlognormal_avg         = time_shlognormal_avg,
    time_shlognormal_disp        = time_shlognormal_disp
  )
  
  res
}

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 8: Run parallelized simulation study
# ─────────────────────────────────────────────────────────────────────────────

# Prepare progress bar
handlers(list(
  handler_progress(
    format   = ":spin :current/:total [:bar] :percent in :elapsed ETA: :eta",
    complete = "+", 
    clear = FALSE
  )
))

# Set up our parallel plan
plan(multisession, workers = 15)

# Wrap runSimulation inside with_progress
with_progress({
  res <- runSimulation(
    design         = Design,
    replications   = 50,
    generate       = Generate, 
    fixed_objects  = fixed_objects,
    analyse        = Analyse,
    summarise      = Summarise,
    parallel       = "future",
    store_results  = TRUE,
    save           = TRUE,
    max_errors     = 5000L,
    ncores         = 2, 
    packages       = c("cmdstanr", "posterior", "MASS", "dplyr", 
                       "fungible", "infinitefactor", "loo"),
    control        = list(save_seeds = TRUE),
    filename       = "Results/Rdata/Supplementary material/sensitivity_simulation_results.rds"
  )
})

# Return to sequential plan
plan(sequential)

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 9: Save full and aggregated results
# ─────────────────────────────────────────────────────────────────────────────

# Full results
saveRDS(object = res, file = "Results/Rdata/Supplementary material/aggregated_sensitivity_results.rds")

# Aggregated results
raw_res <- SimExtract(res, what = "results")
saveRDS(object = raw_res, file = "Results/Rdata/Supplementary material/raw_sensitivity_results.rds")

# ─────────────────────────────────────────────────────────────────────────────
