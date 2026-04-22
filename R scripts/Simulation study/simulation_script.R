# ╔════════════════════════════════════════════════════════════════════════════╗
# ║                             SCRIPT OVERVIEW                                ║
# ╠════════════════════════════════════════════════════════════════════════════╣
# ║ Script Name   : simulation_script.R                                        ║
# ║ Author        : Ricardo Rey-Sáez                                           ║
# ║ Role          : PhD Student in Psychology                                  ║
# ║ Institution   : Autonomous University of Madrid, Spain                     ║
# ║ Email         : ricardoreysaez95@gmail.es                                  ║
# ║ Date          : 25-04-2025                                                 ║
# ╚════════════════════════════════════════════════════════════════════════════╝

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 1: Description
# ─────────────────────────────────────────────────────────────────────────────
# Main workflow for the BGHFM simulation study
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

# Simulation design
Design <- createDesign(Model = c("shifted.lognormal", "exgaussian"),
                       M = c("uni", "2_EFA","2_CFA"),
                       I = c(100, 200),
                       reliab = c(.3, .5, .7),
                       std_lambda_j = c(.3, .5, .7))


# Fixed objects accross all simulation conditions: stan models and parameters
fixed_objects <- list(C_Gauss     = cmdstan_model("Stan models/BGHFM (sim study)/C_BGHFM_gaussian.stan"),
                      E_Gauss     = cmdstan_model("Stan models/BGHFM (sim study)/E_BGHFM_gaussian.stan"),
                      C_Exgauss   = cmdstan_model("Stan models/BGHFM (sim study)/C_BGHFM_exgaussian_tau.stan"),
                      E_Exgauss   = cmdstan_model("Stan models/BGHFM (sim study)/E_BGHFM_exgaussian_tau.stan"),
                      C_shlognorm = cmdstan_model("Stan models/BGHFM (sim study)/C_BGHFM_shlognormal.stan"),
                      E_shlognorm = cmdstan_model("Stan models/BGHFM (sim study)/E_BGHFM_shlognormal.stan"),
                      meta_parameters = meta_parameters)

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 5: Simulation study workflow: Data generation function
# ─────────────────────────────────────────────────────────────────────────────

Generate <- function(condition, fixed_objects) {
  
  # Attach condition to simulate data
  Attach(condition)
  
  # Fixed values
  J <- 6; K <- 2; L <- 100
  
  # ─────────────────────────────────────── #
  #    Specify latent structure matrices    #
  # ─────────────────────────────────────── #
  
  # Specify latent matrices
  if(M == "uni"){
    lambda_std <- matrix(rep(std_lambda_j, J), ncol = 1)
    Phi_cor <- diag(1)
  } else if(M == "2_CFA"){
    lambda_std <- matrix(
      c(rep(std_lambda_j, 3),rep(0, 3),
        rep(0, 3),rep(std_lambda_j, 3)), 
      ncol = 2)
    Phi_cor <- diag(1 - 0.4, ncol = 2, nrow = 2) + 0.4
  } else if(M == "2_EFA"){
    lambda_std <- matrix(
      c(rep(std_lambda_j, 3),rep(std_lambda_j/3, 3),
        rep(std_lambda_j/3, 3),rep(std_lambda_j, 3)), 
      ncol = 2)
    Phi_cor <- diag(1, ncol = 2, nrow = 2)
  }
  
  # ───────────────────────────────────────────────── #
  #    Select model-specific population parameters    #
  # ───────────────────────────────────────────────── #
  
  # Common parameters in both models
  meta_pars <- ifelse(Model == "exgaussian", "meta_exgaussian_tau", "meta_shlognormal")
  mu_alpha_j    = rep(fixed_objects$meta_parameters[[meta_pars]]$mu_alpha, length.out = J)
  mu_theta_j    = rep(fixed_objects$meta_parameters[[meta_pars]]$mu_theta, length.out = J)
  sd_alpha_j    = rep(fixed_objects$meta_parameters[[meta_pars]]$sd_alpha, length.out = J)
  shift_sigma_j = rep(fixed_objects$meta_parameters[[meta_pars]]$shift_sigma, length.out = J)
  scale_sigma_j = rep(fixed_objects$meta_parameters[[meta_pars]]$scale_sigma, length.out = J)
  
  # Model-specific parameters
  if(Model == "exgaussian"){
    shift_tau_j  = rep(fixed_objects$meta_parameters[[meta_pars]]$shift_tau, length.out = J)
    scale_tau_j  = rep(fixed_objects$meta_parameters[[meta_pars]]$scale_tau, length.out = J)
  } else {
    mu_delta_j = rep(fixed_objects$meta_parameters[[meta_pars]]$mu_delta, length.out = J)
    sd_delta_j = rep(fixed_objects$meta_parameters[[meta_pars]]$sd_delta, length.out = J)
  }
  
  # Global parameters
  reliab_j   = rep(reliab, J)
  
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
    
  # Shift parameter 
  if(Model == "shifted.lognormal"){
    delta <- mvrnorm(n = I, mu = mu_delta_j, Sigma = diag(sd_delta_j^2))  
  }
  
  # Rate parameter (rate = 1 / tau)
  if(Model == "exgaussian"){
    mu_tau_j <- shift_tau_j + scale_tau_j * sqrt(2 / pi)
    sd_tau_j <- scale_tau_j * sqrt(1 - 2 / pi)
    
    # Simulate rate values
    rate <- matrix(NA, nrow = I, ncol = J)
    # Simulate level-1 residual values
    for(j in 1:J){
      # Half-normal values
      tau_tilde <- truncnorm::rtruncnorm(n = I, a = 0, b = Inf,
                                         mean = 0, sd = 1)
      # Scaled values with correct E[X] and V[X]
      rate[,j] <- 1 / (shift_tau_j[j] + scale_tau_j[j] * tau_tilde)
    }
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
  
  # Expected value of the trial-level variance under each model
  if(Model == "exgaussian"){
    Ev <- mu_sigma_j^2 + sd_sigma_j^2 + mu_tau_j^2 + sd_tau_j^2
  } else {
    Ev <- mu_sigma_j^2 + sd_sigma_j^2
  }
  
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
  for(i in 1:I){ # For each subjetct
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
          # Level 1 observation: shifted-lognormal
          if(Model == "shifted.lognormal"){
            dat[index,5] <- delta[i,j] + rlnorm(1, mu, sigma[i,j])  
          }
          # Level 1 observation: ex-gaussian
          if(Model == "exgaussian"){
            dat[index,5] <- rnorm(1, mu, sigma[i,j]) + rexp(1, rate[i,j])
          }
        }
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
  
  # Method-of-moments effects, correlations and SEs with Bonnet (2008, p. 174) equation (suggested by Gnambs, 2023)
  # https://online.ucpress.edu/collabra/article/9/1/87615/197169/A-Brief-Note-on-the-Standard-Error-of-the-Pearson
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
  
  # Unidimensional settings
  if(grepl(pattern = "uni", condition$M)){
    stan.data$M <- 1
    stan.data$Lambda_index <- matrix(rep(1,6),ncol=1)
    # Multidimensional settings
  } else {
    stan.data$M <- 2
    stan.data$xi <- 100
    stan.data$Lambda_index <- matrix(c(rep(1,3), rep(0,3),
                                       rep(0,3), rep(1,3)),
                                     ncol=2)
  }
  
  # Define starting values: gaussian model
  gauss_inits <- list(
    make_inits(df = dat, M = stan.data$M, model = "gaussian"), 
    make_inits(df = dat, M = stan.data$M, model = "gaussian")
  )
  
  # Define starting values: skewed model
  skew_inits <- list(
    make_inits(df = dat, M = stan.data$M, model = condition$Model), 
    make_inits(df = dat, M = stan.data$M, model = condition$Model)
  )
  
  # Stan model arguments
  stan_model_arguments <- list(
    data <- stan.data,
    iter_sampling = 1500,
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
  
  # Unidimensional estimation and 2 factors CFA
  if(condition$M == "uni" | condition$M == "2_CFA"){
    # Gaussian model
    stan_model_arguments$init <- gauss_inits
    fit_gauss  <- do.call(fixed_objects$C_Gauss$sample, stan_model_arguments)
    stan_model_arguments$init <- skew_inits
    if(condition$Model == "exgaussian"){
      # Exgaussian model
      fit_skewed <- do.call(fixed_objects$C_Exgauss$sample, stan_model_arguments)
    } else {
      # Shifted-lognormal model
      fit_skewed <- do.call(fixed_objects$C_shlognorm$sample, stan_model_arguments)
    }
    # 2 factors EFA
  } else {
    # Gaussian model
    stan_model_arguments$init <- gauss_inits
    fit_gauss  <- do.call(fixed_objects$E_Gauss$sample, stan_model_arguments)
    # Check if exists a unique rotation w.r.t the population factor loadings (if not, redraw)
    # True target matrix
    target <- matrix(c(rep(condition$std_lambda_j, 3),  rep(condition$std_lambda_j/3, 3),
                       rep(condition$std_lambda_j/3, 3),rep(condition$std_lambda_j, 3)), 
                     ncol = 2)
    # Rotate factor loadings matrix
    rotated_gauss  <- simulation_alignment(fit = fit_gauss, target = target)
    
    # Make initial values for skewed models
    stan_model_arguments$init <- skew_inits
    if(condition$Model == "exgaussian"){
      # Exgaussian model
      fit_skewed <- do.call(fixed_objects$E_Exgauss$sample, stan_model_arguments)
    } else {
      # Shifted-lognormal model
      fit_skewed <- do.call(fixed_objects$E_shlognorm$sample, stan_model_arguments)
    }
    
    # Rotate skewed factor loadings
    rotated_skewed <- simulation_alignment(fit = fit_skewed, target = target)
  }
  
  # ───────────────────────── #
  #    Parameter estimates    #
  # ───────────────────────── #
  
  # Posterior parameter summaries
  mu_alpha <- posterior_summaries(fit_skewed, "mu_alpha")
  mu_theta <- posterior_summaries(fit_skewed, "mu_theta")
  mu_sigma <- posterior_summaries(fit_skewed, "mu_sigma")
  sd_alpha <- posterior_summaries(fit_skewed, "sd_alpha")
  sd_theta <- posterior_summaries(fit_skewed, "sd_theta")
  sd_sigma <- posterior_summaries(fit_skewed, "sd_sigma")
  
  # Model-implied correlations and signal-to-noise ratios per model
  Rho_gauss <- posterior_summaries(fit_gauss, "Rho_theta", cormat = TRUE, parnames = "Rho_gauss")
  Rho_skew  <- posterior_summaries(fit_skewed, "Rho_theta", cormat = TRUE, parnames = "Rho_skew")
  reliab_gauss <- posterior_summaries(fit_gauss, "reliability", parnames = "reliab_gauss")
  reliab_skew  <- posterior_summaries(fit_skewed, "reliability", parnames = "reliab_skew")
  
  # Empty values (depends on latent structure)
  Phi_gauss <- posterior_summaries(fit_gauss, "Phi_gauss", empty = TRUE, times = 1)
  Phi_skew  <- posterior_summaries(fit_skewed, "Phi_skew", empty = TRUE, times = 1)
  crossload_gauss <- posterior_summaries(fit_gauss, "crossload_gauss", empty = TRUE, times = 6)
  crossload_skew <- posterior_summaries(fit_skewed, "crossload_skew", empty = TRUE, times = 6)
  
  # Empty values (depends on specific model)
  mu_delta <- posterior_summaries(fit_skewed, "mu_delta", empty = TRUE, times = 6)
  sd_delta <- posterior_summaries(fit_skewed, "sd_delta", empty = TRUE, times = 6)
  mu_tau <- posterior_summaries(fit_skewed, "mu_tau", empty = TRUE, times = 6)
  sd_tau <- posterior_summaries(fit_skewed, "sd_tau", empty = TRUE, times = 6)
  
  # Unidimensional settings
  if(condition$M == "uni"){
    # Gaussian and skewed lambdas
    lambda_gauss <- posterior_summaries(fit = fit_gauss, parameter = "Lambda_std", parnames = "lambda_gauss")
    lambda_skew <- posterior_summaries(fit = fit_skewed, parameter = "Lambda_std", parnames = "lambda_skew")
    
    # Tucker index of factor congruence: Gaussian model
    Tucker_phi_gauss <- diag(
      psych::factor.congruence(
        x = matrix(rep(condition$std_lambda_j, 6), ncol = 1), 
        y = matrix(fit_gauss$summary("Lambda_std")$mean, nrow = 6, ncol = 1)
      )
    )
    
    # Tucker index of factor congruence: Skew model
    Tucker_phi_skew <- diag(
      psych::factor.congruence(
        x = matrix(rep(condition$std_lambda_j, 6), ncol = 1), 
        y = matrix(fit_skewed$summary("Lambda_std")$mean, nrow = 6, ncol = 1)
      )
    )
    
    # Tucker mean and sd
    Tucker_gauss_f1 <- Tucker_phi_gauss
    Tucker_gauss_f2 <- 0
    Tucker_skew_f1 <- Tucker_phi_skew
    Tucker_skew_f2 <- 0
  }
  
  # Exploratory factor analysis
  if(condition$M == "2_EFA"){
    # Main factor loadings
    lambda_gauss <- posterior_summaries(rotated_gauss[c(1:3,10:12),],
                                        parameter = "Lambda_std", 
                                        rotated = TRUE, 
                                        parnames = "lambda_gauss")
    lambda_skew  <- posterior_summaries(rotated_skewed[c(1:3,10:12),], 
                                        parameter = "Lambda_std", 
                                        rotated = TRUE, 
                                        parnames = "lambda_skew")
    # Cross-loadings
    crossload_gauss <- posterior_summaries(rotated_gauss[4:9,], 
                                           parameter = "Lambda_std", 
                                           rotated = TRUE, 
                                           parnames = "crossload_gauss")
    crossload_skew  <- posterior_summaries(rotated_skewed[4:9,], 
                                           parameter = "Lambda_std", 
                                           rotated = TRUE, 
                                           parnames = "crossload_skew")
    
    # Tucker index of factor congruence: Gaussian model
    Tucker_phi_gauss <- diag(
      psych::factor.congruence(
        x = target, 
        y = matrix(rotated_gauss$mean, nrow = 6, ncol = 2)
      )
    )
    
    # Tucker index of factor congruence: Skew model
    Tucker_phi_skew <- diag(
      psych::factor.congruence(
        x = target, 
        y = matrix(rotated_skewed$mean, nrow = 6, ncol = 2)
      )
    )
    
    # Tucker mean and sd
    Tucker_gauss_f1 <- Tucker_phi_gauss[1]
    Tucker_gauss_f2 <- Tucker_phi_gauss[2]
    Tucker_skew_f1  <- Tucker_phi_skew[1]
    Tucker_skew_f2  <- Tucker_phi_skew[2]
  }
  
  # Exploratory factor analysis
  if(condition$M == "2_CFA"){
    # Factor loadings
    lambda_gauss <- posterior_summaries(fit_gauss, 
                                        "Lambda_std", 
                                        CFA_loadings = TRUE, 
                                        parnames = "lambda_gauss")
    lambda_skew <- posterior_summaries(fit_skewed, 
                                       "Lambda_std",
                                       CFA_loadings = TRUE, 
                                       parnames = "lambda_skew")
    # Common-factor correlation matrix
    Phi_gauss <- posterior_summaries(fit_gauss, "Phi_lv", 
                                     CFA_phi = TRUE, 
                                     parnames = "Phi_gauss")
    Phi_skew <- posterior_summaries(fit_skewed, "Phi_lv",
                                    CFA_phi = TRUE, 
                                    parnames = "Phi_skew")
    # Pop. Lambda matrix
    L_true <- matrix(
      c(rep(condition$std_lambda_j, 3),rep(0, 3),
        rep(0, 3),rep(condition$std_lambda_j, 3)), 
      ncol = 2)
    
    # Tucker index of factor congruence: Gaussian model
    Tucker_phi_gauss <- diag(
      psych::factor.congruence(
        x = L_true, 
        y = matrix(fit_gauss$summary("Lambda_std")$mean, nrow = 6, ncol = 2)
      )
    )
    
    # Tucker index of factor congruence: Skew model
    Tucker_phi_skew <- diag(
      psych::factor.congruence(
        x = L_true, 
        y = matrix(fit_skewed$summary("Lambda_std")$mean, nrow = 6, ncol = 2)
      )
    )
    
    # Tucker mean and sd
    Tucker_gauss_f1 <- Tucker_phi_gauss[1]
    Tucker_gauss_f2 <- Tucker_phi_gauss[2]
    Tucker_skew_f1  <- Tucker_phi_skew[1]
    Tucker_skew_f2  <- Tucker_phi_skew[2]
  }
  
  # Extra parameters
  if(condition$Model == "exgaussian"){
    # Save ex-gaussian parameters
    mu_tau <- posterior_summaries(fit_skewed, "mu_tau")
    sd_tau <- posterior_summaries(fit_skewed, "sd_tau")
  } else {
    # Save shifted-lognormal parameters
    mu_delta <- posterior_summaries(fit_skewed, "mu_delta")
    sd_delta <- posterior_summaries(fit_skewed, "sd_delta")
  }
  
  # Compute PSIS-LOO
  loo_gaussian <- loo_BHGFM(fit_gauss, model = "gaussian", data = dat, n.cores = 2)
  loo_skewed <- loo_BHGFM(fit_skewed, model = condition$Model, data = dat, n.cores = 2)
  
  # LOO model comparison
  loo_comparison <- loo_compare(list(skew = loo_skewed, gauss = loo_gaussian))
  
  # Best model
  skew_first <- rownames(loo_comparison)[1] == "skew"
  
  # Soft and hard null-CIs
  loo_null_weak_CI <- loo_comparison[2,1] + c(qnorm(p = .025), qnorm(p = .975)) * loo_comparison[2,2]
  loo_null_strong_CI <- loo_comparison[2,1] + c(-4, 4) * loo_comparison[2,2]
  
  loo_weak_diff <- !(loo_null_weak_CI[1] <= 0 & 0 <= loo_null_weak_CI[2])
  loo_strong_diff <- !(loo_null_strong_CI[1] <= 0 & 0 <= loo_null_strong_CI[2])
  
  # Selected as best model
  skew_weak_diff <- all(skew_first, loo_weak_diff)
  skew_strong_diff <- all(skew_first, loo_strong_diff)
  
  # Probability that the skewed model outperforms the null model
  if(skew_first) {
    prob_skew <- pnorm(0, mean = loo_comparison[2,1], sd = loo_comparison[2,2])
  } else {
    prob_skew <- 999 # Don't use
  }
  
  # ──────────────────────── #
  #    Simulation results    #    
  # ──────────────────────── #
  
  ret <- c(
    unlist(mu_alpha), unlist(mu_theta), unlist(mu_sigma), unlist(mu_tau), unlist(mu_delta),
    unlist(sd_alpha), unlist(sd_theta), unlist(sd_sigma), unlist(sd_tau), unlist(sd_delta),
    unlist(Phi_gauss), unlist(Phi_skew), unlist(lambda_gauss), unlist(lambda_skew),
    unlist(crossload_gauss), unlist(crossload_skew), unlist(Rho_gauss), unlist(Rho_skew),
    unlist(reliab_gauss), unlist(reliab_skew), 
    unlist(Rho_mm), unlist(Rho_sp), 
    skew_first = skew_first,
    Tucker_gauss_f1 = Tucker_gauss_f1,
    Tucker_gauss_f2 = Tucker_gauss_f2,
    Tucker_skew_f1  = Tucker_skew_f1,
    Tucker_skew_f2  = Tucker_skew_f2,
    prob_skew = prob_skew,
    loo_elpd_diff = loo_comparison[2,1],
    loo_diff_weak = loo_weak_diff,
    loo_diff_strong = loo_strong_diff,
    loo_diff_skew_weak = skew_weak_diff,
    loo_diff_skew_strong = skew_strong_diff,
    skew_time = fit_skewed$time()$total, gauss_time = fit_gauss$time()$total,
    gauss_divergence = fit_gauss$diagnostic_summary()$num_divergent,
    skew_divergence = fit_skewed$diagnostic_summary()$num_divergent,
    gauss_treedepth = fit_gauss$diagnostic_summary()$num_max_treedepth,
    skew_treedepth = fit_skewed$diagnostic_summary()$num_max_treedepth,
    gauss_ebfmi = fit_gauss$diagnostic_summary()$ebfmi,
    skew_ebfmi = fit_skewed$diagnostic_summary()$ebfmi
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
  I <- condition$I; J <- 6; K <- 2; L <- 100
  
  # Specify latent matrices
  if(condition$M == "uni"){
    lambda_std <- matrix(rep(condition$std_lambda_j, J), ncol = 1)
    Phi_cor <- diag(1)
  } else if(condition$M == "2_CFA"){
    lambda_std <- matrix(
      c(rep(condition$std_lambda_j, 3),rep(0, 3),
        rep(0, 3),rep(condition$std_lambda_j, 3)), 
      ncol = 2)
    Phi_cor <- diag(1 - 0.4, ncol = 2, nrow = 2) + 0.4
  } else if(condition$M == "2_EFA"){
    lambda_std <- matrix(
      c(rep(condition$std_lambda_j, 3),rep(condition$std_lambda_j/3, 3),
        rep(condition$std_lambda_j/3, 3),rep(condition$std_lambda_j, 3)), 
      ncol = 2)
    Phi_cor <- diag(1, ncol = 2, nrow = 2)
  }
  
  # Initial empty parameters
  mu_alpha_j <- mu_theta_j <- log_mu_sigma_j <- log_mu_tau_j <- mu_delta_j <- rep(NA, J)
  sd_alpha_j <- log_sd_sigma_j <- log_sd_tau_j <- sd_delta_j <- rep(NA, J)
  
  # Common parameters in both models
  meta_pars <- ifelse(condition$Model == "exgaussian", "meta_exgaussian_tau", "meta_shlognormal")
  mu_alpha_j    = rep(fixed_objects$meta_parameters[[meta_pars]]$mu_alpha, length.out = J)
  mu_theta_j    = rep(fixed_objects$meta_parameters[[meta_pars]]$mu_theta, length.out = J)
  sd_alpha_j    = rep(fixed_objects$meta_parameters[[meta_pars]]$sd_alpha, length.out = J)
  shift_sigma_j = rep(fixed_objects$meta_parameters[[meta_pars]]$shift_sigma, length.out = J)
  scale_sigma_j = rep(fixed_objects$meta_parameters[[meta_pars]]$scale_sigma, length.out = J)
  
  # Model-specific parameters
  if(condition$Model == "exgaussian"){
    shift_tau_j = rep(fixed_objects$meta_parameters[[meta_pars]]$shift_tau, length.out = J)
    scale_tau_j = rep(fixed_objects$meta_parameters[[meta_pars]]$scale_tau, length.out = J)
  } else {
    mu_delta_j = rep(fixed_objects$meta_parameters[[meta_pars]]$mu_delta, length.out = J)
    sd_delta_j = rep(fixed_objects$meta_parameters[[meta_pars]]$sd_delta, length.out = J)
  }
  
  # Global parameters
  reliab_j = rep(condition$reliab, J)
  
  # Level-1 residual variance expected value and variance
  mu_sigma_j <- shift_sigma_j + scale_sigma_j * sqrt(2 / pi)
  sd_sigma_j <- scale_sigma_j * sqrt(1 - 2 / pi)
  
  # Exponential rate expected value and variance
  if(condition$Model == "exgaussian"){
    mu_tau_j <- shift_tau_j + scale_tau_j * sqrt(2 / pi)
    sd_tau_j <- scale_tau_j * sqrt(1 - 2 / pi)
  }
  
  # Model communality
  communality <- lambda_std %*% Phi_cor %*% t(lambda_std)
  
  # Population model-implied correlation matrix
  R_theta_j <- communality + diag(1 - diag(communality))
  
  # Signal-to-Noise ratio
  gamma2_j <- (2 * reliab_j) / (L * (1 - reliab_j))
  
  # Expected value of the trial-level variance under each model
  if(condition$Model == "exgaussian"){
    Ev <- mu_sigma_j^2 + sd_sigma_j^2 + mu_tau_j^2 + sd_tau_j^2
  } else {
    Ev <- mu_sigma_j^2 + sd_sigma_j^2
  }
  
  # Population standard deviation
  sd_theta_j <- sqrt(gamma2_j * Ev)
  
  # ───────────────────────────────────── #
  #    bias in all parameter estimates    #
  # ───────────────────────────────────── #
  
  # Global summaries per parameter
  mu_alpha <- global_summaries(results = results, parameter_name = "mu_alpha", parameter = mu_alpha_j)
  mu_theta <- global_summaries(results = results, parameter_name = "mu_theta", parameter = mu_theta_j)
  mu_sigma <- global_summaries(results = results, parameter_name = "mu_sigma", parameter = mu_sigma_j)
  sd_alpha <- global_summaries(results = results, parameter_name = "sd_alpha", parameter = sd_alpha_j)
  sd_theta <- global_summaries(results = results, parameter_name = "sd_theta", parameter = sd_theta_j)
  sd_sigma <- global_summaries(results = results, parameter_name = "sd_sigma", parameter = sd_sigma_j)
  
  # Summaries for correlation matrices
  Rho_gauss <- global_summaries(results = results, parameter_name = "Rho_gauss", parameter = R_theta_j, Rho = TRUE)
  Rho_skew <- global_summaries(results = results, parameter_name = "Rho_skew", parameter = R_theta_j, Rho = TRUE)
  Rho_mm <- global_summaries(results = results, parameter_name = "Rho_mm", parameter = R_theta_j, method_of_moments = TRUE)
  Rho_sp <- global_summaries(results = results, parameter_name = "Rho_sp", parameter = R_theta_j, method_of_moments = TRUE)
  
  # Summaries for reliability
  reliab_gauss <- global_summaries(results = results, parameter_name = "reliab_gauss", parameter = reliab_j)
  reliab_skew  <- global_summaries(results = results, parameter_name = "reliab_skew", parameter = reliab_j)
  
  # Empty values (depends on latent structure)
  Phi_gauss <- global_summaries(parameter_name = "Phi_gauss", empty = TRUE, Phi.lv = TRUE)
  Phi_skew <- global_summaries(parameter_name = "Phi_skew", empty = TRUE, Phi.lv = TRUE)
  crossload_gauss <- global_summaries(parameter_name = "crossload_gauss", empty = TRUE)
  crossload_skew <- global_summaries(parameter_name = "crossload_skew", empty = TRUE)
  
  # Empty values (depends on specific model)
  mu_delta <- global_summaries(parameter_name = "mu_delta", empty = TRUE)
  sd_delta <- global_summaries(parameter_name = "sd_delta", empty = TRUE)
  mu_tau <- global_summaries(parameter_name = "mu_tau", empty = TRUE)
  sd_tau <- global_summaries(parameter_name = "sd_tau", empty = TRUE)
  
  # Unidimensional settings
  if(condition$M == "uni"){
    # Gaussian and skewed lambdas
    lambda_gauss <- global_summaries(results = results, 
                                     parameter_name = "lambda_gauss", 
                                     parameter = lambda_std[,1])
    lambda_skew <- global_summaries(results = results, 
                                    parameter_name = "lambda_skew", 
                                    parameter = lambda_std[,1])
  }
  
  # Exploratory factor analysis
  if(condition$M == "2_EFA"){
    # Main factor loadings
    lambda_gauss <- global_summaries(results = results, 
                                     parameter_name = "lambda_gauss", 
                                     parameter = c(lambda_std[1:3,1], 
                                                   lambda_std[4:6,2]))
    lambda_skew  <- global_summaries(results = results, 
                                     parameter_name = "lambda_skew", 
                                     parameter = c(lambda_std[1:3,1], 
                                                   lambda_std[4:6,2]))
    # Cross-loadings
    crossload_gauss <- global_summaries(results = results, 
                                        parameter_name = "crossload_gauss", 
                                        parameter = c(lambda_std[4:6,1], 
                                                      lambda_std[1:3,2]))
    crossload_skew  <- global_summaries(results = results, 
                                        parameter_name = "crossload_skew", 
                                        parameter = c(lambda_std[4:6,1], 
                                                      lambda_std[1:3,2]))
  }
  
  # Exploratory factor analysis
  if(condition$M == "2_CFA"){
    # Factor loadings
    lambda_gauss <- global_summaries(results = results, 
                                     parameter_name = "lambda_gauss", 
                                     parameter = c(lambda_std[1:3,1], 
                                                   lambda_std[4:6,2]))
    lambda_skew  <- global_summaries(results = results, 
                                     parameter_name = "lambda_skew", 
                                     parameter = c(lambda_std[1:3,1], 
                                                   lambda_std[4:6,2]))
    # Common-factor correlation matrix
    Phi_gauss <- global_summaries(results = results, 
                                  parameter_name = "Phi_gauss", 
                                  parameter = Phi_cor[1,2], 
                                  Phi.lv = TRUE)
    Phi_skew <- global_summaries(results = results, 
                                 parameter_name = "Phi_skew",
                                 parameter = Phi_cor[1,2], 
                                 Phi.lv = TRUE)
  }
  
  # Extra parameters
  if(condition$Model == "exgaussian"){
    # Save ex-gaussian parameters
    mu_tau <- global_summaries(results = results, parameter_name = "mu_tau", parameter = mu_tau_j)
    sd_tau <- global_summaries(results = results, parameter_name = "sd_tau", parameter = sd_tau_j)
  } else {
    # Save shifted-lognormal parameters
    mu_delta <- global_summaries(results = results, parameter_name = "mu_delta", parameter = mu_delta_j)
    sd_delta <- global_summaries(results = results, parameter_name = "sd_delta", parameter = sd_delta_j)
  }
  
  # PSIS-LOO expeted log-pointwise posterior density (ELPD)
  loo_avg_elpd <- mean(results[,grep("loo_elpd", colnames(results))])
  loo_disp_elpd <- sd(results[,grep("loo_elpd", colnames(results))])
  
  # Model comparison: PSIS-LOO differences
  loo_avg_weak <- mean(results[,grep("loo_diff_weak", colnames(results))])
  loo_disp_weak <- sd(results[,grep("loo_diff_weak", colnames(results))])
  loo_avg_strong <- mean(results[,grep("loo_diff_strong", colnames(results))])
  loo_disp_strong <- sd(results[,grep("loo_diff_strong", colnames(results))])
  
  # Model comparison: PSIS-LOO differences in favour of skewed mdel
  loo_avg_skew_weak <- mean(results[,grep("loo_diff_skew_weak", colnames(results))])
  loo_disp_skew_weak <- sd(results[,grep("loo_diff_skew_weak", colnames(results))])
  loo_avg_skew_strong <- mean(results[,grep("loo_diff_skew_strong", colnames(results))])
  loo_disp_skew_strong <- sd(results[,grep("loo_diff_skew_strong", colnames(results))])
  
  # Diagnostic summaries: average ebfmi and time
  ebfmi_gauss_avg  <- mean(colMeans(results[, grep("gauss_ebfmi", colnames(results))]))
  ebfmi_gauss_disp <- mean(apply(results[, grep("gauss_ebfmi", colnames(results))], 2, sd))
  ebfmi_skew_avg  <- mean(colMeans(results[, grep("skew_ebfmi", colnames(results))]))
  ebfmi_skew_disp <- mean(apply(results[, grep("skew_ebfmi", colnames(results))], 2, sd))
  time_gauss_avg  <- mean(results[, grep("gauss_time", colnames(results))])
  time_gauss_disp <- sd(results[, grep("gauss_time", colnames(results))])
  time_skew_avg   <- mean(results[, grep("skew_time", colnames(results))])
  time_skew_disp  <- sd(results[, grep("skew_time", colnames(results))])
  
  # Diagnostic summaries: number of divergent iterations
  gauss_divergence_n <- sum(rowSums(results[, grep("gauss_divergence", colnames(results))]) > 0)
  skew_divergence_n  <- sum(rowSums(results[, grep("skew_divergence", colnames(results))]) > 0)
  # Percentage of divergent iterations (avg and sd)
  gauss_divergence_avg  <- mean(rowSums(results[, grep("gauss_divergence", colnames(results))]) / 3000 * 100)
  gauss_divergence_disp <- sd(rowSums(results[, grep("gauss_divergence", colnames(results))]) / 3000 * 100)
  skew_divergence_avg   <- mean(rowSums(results[, grep("skew_divergence", colnames(results))]) / 3000 * 100)
  skew_divergence_disp  <- sd(rowSums(results[, grep("skew_divergence", colnames(results))]) / 3000 * 100)
  
  # Diagnostic summaries: number of maximum treedepth iterations
  gauss_treedepth_n <- sum(rowSums(results[, grep("gauss_treedepth", colnames(results))]) > 0)
  skew_treedepth_n  <- sum(rowSums(results[, grep("skew_treedepth", colnames(results))]) > 0)
  # Percentage of maximum treedepth iterations (avg and sd)
  gauss_treedepth_avg  <- mean(rowSums(results[, grep("gauss_treedepth", colnames(results))]) / 3000 * 100)
  gauss_treedepth_disp <- sd(rowSums(results[, grep("gauss_treedepth", colnames(results))]) / 3000 * 100)
  skew_treedepth_avg   <- mean(rowSums(results[, grep("skew_treedepth", colnames(results))]) / 3000 * 100)
  skew_treedepth_disp  <- sd(rowSums(results[, grep("skew_treedepth", colnames(results))]) / 3000 * 100)
  
  # Skew model always selected as the best?
  skew_first_avg  <- mean(results[,"skew_first"])
  skew_first_disp <- sd(results[,"skew_first"])
  
  # Probability that the skew model outperform the gaussian model in predicting new data. 
  prob_skew_avg  <- mean(results[,"prob_skew"][results[,"prob_skew"] != 999])
  prob_skew_disp <- sd(results[,"prob_skew"][results[,"prob_skew"] != 999])

  # Tucker congruence index
  Tucker_gauss_f1_avg <- mean(results[,"Tucker_gauss_f1"])
  Tucker_gauss_f2_avg <- mean(results[,"Tucker_gauss_f2"])
  Tucker_skew_f1_avg  <- mean(results[,"Tucker_skew_f1"])
  Tucker_skew_f2_avg  <- mean(results[,"Tucker_skew_f2"])
  Tucker_gauss_f1_disp <- sd(results[,"Tucker_gauss_f1"])
  Tucker_gauss_f2_disp <- sd(results[,"Tucker_gauss_f2"])
  Tucker_skew_f1_disp  <- sd(results[,"Tucker_skew_f1"])
  Tucker_skew_f2_disp  <- sd(results[,"Tucker_skew_f2"])

  # Final results
  res <- c(
    Rho_gauss, Rho_skew, 
    Rho_mm, Rho_sp,
    lambda_gauss, lambda_skew, reliab_gauss, reliab_skew, 
    mu_alpha, mu_theta, mu_sigma, mu_tau, mu_delta,
    sd_alpha, sd_theta, sd_sigma, sd_tau, sd_delta,
    crossload_gauss, crossload_skew, Phi_gauss, Phi_skew,
    Tucker_gauss_f1_avg     = Tucker_gauss_f1_avg,
    Tucker_gauss_f2_avg     = Tucker_gauss_f2_avg,
    Tucker_skew_f1_avg      = Tucker_skew_f1_avg,
    Tucker_skew_f2_avg      = Tucker_skew_f2_avg,
    Tucker_gauss_f1_disp    = Tucker_gauss_f1_disp,
    Tucker_gauss_f2_disp    = Tucker_gauss_f2_disp,
    Tucker_skew_f1_disp     = Tucker_skew_f1_disp,
    Tucker_skew_f2_disp     = Tucker_skew_f2_disp,
    skew_first_avg          = skew_first_avg,
    skew_first_disp         = skew_first_disp,
    prob_skew_avg           = prob_skew_avg,
    prob_skew_disp          = prob_skew_avg,
    loo_avg_elpd            = loo_avg_elpd,
    loo_disp_elpd           = loo_disp_elpd,
    loo_avg_both_weak       = loo_avg_weak,
    loo_disp_both_weak      = loo_disp_weak,
    loo_avg_both_strong     = loo_avg_strong,
    loo_disp_both_strong    = loo_disp_strong,
    loo_avg_skew_weak       = loo_avg_skew_weak,
    loo_disp_skew_weak      = loo_disp_skew_weak,
    loo_avg_skew_strong     = loo_avg_skew_strong,
    loo_disp_skew_strong    = loo_disp_skew_strong,
    divergence_gauss_n      = gauss_divergence_n, 
    divergence_skew_n       = skew_divergence_n,
    divergence_gauss_avg    = gauss_divergence_avg,
    divergence_skew_avg     = skew_divergence_avg,
    divergence_gauss_disp   = gauss_divergence_disp,
    divergence_skew_disp    = skew_divergence_disp,
    treedepth_gauss_n       = gauss_treedepth_n,
    treedepth_skew_n        = skew_treedepth_n,
    treedepth_gauss_avg     = gauss_treedepth_avg,
    treedepth_gauss_disp    = gauss_treedepth_disp,
    treedepth_skew_avg      = skew_treedepth_avg,
    treedepth_skew_disp     = skew_treedepth_disp,
    ebfmi_gauss_avg         = ebfmi_gauss_avg,
    ebfmi_gauss_disp        = ebfmi_gauss_disp,
    ebfmi_skew_avg          = ebfmi_skew_avg,
    ebfmi_skew_disp         = ebfmi_skew_disp,
    time_gauss_avg          = time_gauss_avg,
    time_gauss_disp         = time_gauss_disp,
    time_skew_avg           = time_skew_avg,
    time_skew_disp          = time_skew_disp
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
plan(multisession, workers = 8)

# Wrap runSimulation inside with_progress
with_progress({
  res <- runSimulation(
    design         = Design,
    replications   = 100,
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
    filename       = "Results/Rdata/Simulation study/BGHFM_simulation_results.rds"
  )
})

# Return to sequential plan
plan(sequential)

# Full results
res <- readRDS("Results/Rdata/Simulation study/BGHFM_simulation_results.rds")

# For some weird reason, the first condition didn't store the full results
SimExtract(res, what = "results") # Only two rows for condition 1.

# We just repeat that simulation condition.

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 9: Repeating the first condition
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
plan(multisession, workers = 8)

# Wrap runSimulation inside with_progress
with_progress({
  res_cond1 <- runSimulation(
    design         = Design[1,],
    replications   = 100,
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
    filename       = "Results/Rdata/Simulation study/BGHFM_simulation_results_cond1.rds"
  )
})

# Return to sequential plan
plan(sequential)

# Save aggregated and full results in unique objects
res_cond1 <- readRDS("Results/Rdata/Simulation study/BGHFM_simulation_results_cond1.rds")
res[1,] <- res_cond1[1,]
raw_res <- rbind(
  # New first condition
  SimExtract(res_cond1, what = "results"),
  # Remove the wrong first two observations
  SimExtract(res, what = "results")[-c(1:2),])

# Save both objects
saveRDS(object = res, file = "Results/Rdata/Simulation study/aggregated_simulation_results.rds")
saveRDS(object = raw_res, file = "Results/Rdata/Simulation study/raw_simulation_results.rds")
  
# ─────────────────────────────────────────────────────────────────────────────
