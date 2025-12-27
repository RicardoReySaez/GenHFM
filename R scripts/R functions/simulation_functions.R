# ╔════════════════════════════════════════════════════════════════════════════╗
# ║                             SCRIPT OVERVIEW                                ║
# ╠════════════════════════════════════════════════════════════════════════════╣
# ║ Script Name   : simulation_functions.R                                     ║
# ║ Author        : Ricardo Rey-Sáez                                           ║
# ║ Role          : PhD Student in Psychology                                  ║
# ║ Institution   : Autonomous University of Madrid, Spain                     ║
# ║ Email         : ricardoreysaez95@gmail.es                                  ║
# ║ Date          : 25-04-2025                                                 ║
# ╚════════════════════════════════════════════════════════════════════════════╝

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 1: Description
# ─────────────────────────────────────────────────────────────────────────────

# R functions necessary in the simulation study

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 2: Load Packages
# ─────────────────────────────────────────────────────────────────────────────

library(MASS)  # Multivariate gaussian distribution
library(dplyr) # Stan array data
library(posterior) # cmdstanr summaries
library(infinitefactor) # MatchAlign algorithm functions
library(fungible) # Align rotated matrix with population matrix

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 3: User-Defined Functions
# ─────────────────────────────────────────────────────────────────────────────

# Ensure packages inside each function; if not, install it.
ensure_packages <- function(packages) {
  # Iterar sobre la lista de paquetes
  for (pkg in packages) {
    # Verificar si el paquete está instalado
    if (!requireNamespace(pkg, quietly = TRUE)) {
      # Si no está instalado, instalarlo
      message(sprintf("Instalando el paquete: %s", pkg))
      install.packages(pkg, dependencies = TRUE)
    }
    # Cargar el paquete
    library(pkg, character.only = TRUE)
  }
}

# Long data frame to Stan data in list format
Stan.list.data <- function(x, subject_var = "subject", task_var = "task", 
                           condition_var = "condition", rt_var = "RT") {
  # Ensure dplyr
  ensure_packages(packages = "dplyr")
  
  # Rename columns, remove NAs
  data_renamed <- x %>%
    rename(
      subject   = all_of(subject_var),
      condition = all_of(condition_var),
      task      = all_of(task_var),
      RT        = all_of(rt_var)
    ) %>%
    na.omit()
  
  # Filter out subjects who do not have observations in all condition x task combos
  n_cond_total <- n_distinct(data_renamed$condition)
  n_task_total <- n_distinct(data_renamed$task)
  n_comb_total <- n_cond_total * n_task_total
  
  # Count number of distinct (condition, task) combos per subject
  subject_counts <- data_renamed %>%
    group_by(subject) %>%
    summarise(n_combos = n_distinct(interaction(condition, task)),
              .groups = "drop")
  
  # Identify who has all combos
  valid_subjects <- subject_counts %>%
    filter(n_combos == n_comb_total) %>%
    pull(subject)
  
  # Warn about excluded subjects
  excluded_subjects <- setdiff(unique(data_renamed$subject), valid_subjects)
  if (length(excluded_subjects) > 0) {
    warning(
      "The following subjects have been removed for not having observations in *all* task x condition combos:\n",
      paste(excluded_subjects, collapse = ", ")
    )
  }
  
  # Filter data to keep only valid subjects
  data_valid <- data_renamed %>%
    filter(subject %in% valid_subjects)
  
  # Re-factor (numeric) subject, condition, and task
  data_valid <- data_valid %>%
    mutate(
      subject   = as.numeric(factor(subject)),
      condition = as.numeric(factor(condition)),
      task      = as.numeric(factor(task))
    )
  
  # Number of subjects, conditions, tasks
  n_subj <- n_distinct(data_valid$subject)
  n_cond <- n_distinct(data_valid$condition)
  n_task <- n_distinct(data_valid$task)
  
  # Summaries: T_subj, RT_min, RT_max
  # First, summarise by (subject, condition, task)
  df_summ <- data_valid %>%
    group_by(subject, condition, task) %>%
    summarise(
      n_trials = n(),
      rt_min   = min(RT),
      rt_max   = max(RT),
      .groups  = "drop"
    )
  
  # Convert these summaries into 3D arrays with `xtabs()`
  T_subj <- xtabs(n_trials ~ subject + condition + task, data = df_summ)
  RT_min <- xtabs(rt_min   ~ subject + condition + task, data = df_summ)
  RT_max <- xtabs(rt_max   ~ subject + condition + task, data = df_summ)
  
  # Maximum number of trials across all (subject, condition, task)
  T_max <- max(df_summ$n_trials)
  
  # Build the 4D array of RT in a single pass (subject, condition, task, trial_id)
  # Create a trial_id from 1..n within each group
  data_with_id <- data_valid %>%
    group_by(subject, condition, task) %>%
    mutate(trial_id = row_number()) %>%
    ungroup()
  
  # Initialize the 4D array
  RT <- array(0, dim = c(n_subj, n_cond, n_task, T_max))
  
  # Build index matrix [row, col, slice, fourth_dim]
  # cbind() yields an N x 4 matrix for all rows
  idx_matrix <- cbind(
    data_with_id$subject,
    data_with_id$condition,
    data_with_id$task,
    data_with_id$trial_id
  )
  
  # Fill the RT array
  RT[idx_matrix] <- data_with_id$RT
  
  # Return the final list for Stan
  stan_dat <- list(
    I      = n_subj,
    K      = n_cond,
    J      = n_task,
    L_max  = T_max,
    T_subj = T_subj,
    RT     = RT,
    RT_min = RT_min,
    RT_max = RT_max
  )
  
  return(stan_dat)
}

# Informative starting values per model
make_inits <- function(df, M, model = NULL){
  
  # Ensure truncated gaussian distribution
  ensure_packages("truncnorm")
  
  # Gaussian Copula Function: standardized half-normal distribution
  scale_half_std_normal <- function(x){
    # Scale observed variable
    std_x <- (x - mean(x)) / sd(x)
    # Compute cumulative density values with std gaussian distribution
    p_x <- pnorm(std_x, mean = 0, sd = 1)
    # Generate standardized half-normal values
    y <- truncnorm::qtruncnorm(p = p_x, a = 0, b = Inf, mean = 0, sd = 1)
    return(y)
  }
  
  # Unit-vector parameters
  J <- length(unique(df$task))
  z_mat <- matrix(runif(length(unique(df$task)) * M, -1, 1), nrow = J, ncol=M)
  row_norms <- sqrt(rowSums(z_mat^2))
  z_normed  <- z_mat / row_norms
  
  if(model == "gaussian"){
    # Compute individual-level parameters
    alpha <- tapply(df$RT, list(df$subject, df$task), mean)
    pre_theta <- tapply(df$RT, list(df$subject, df$condition, df$task), mean)
    theta <- pre_theta[,2,] - pre_theta[,1,]
    sigma <- tapply(df$RT, list(df$subject, df$task), sd)
    # Exploratory factor analysis: varimax
    h2 <- psych::fa(theta, M, fm = "minrank", rotate = "varimax")$communality
    # Exploratory factor analysis: oblimin
    oblFA <- psych::fa(theta, M, fm = "minrank", rotate = "oblimin")
    if(M > 1){ 
      L_Phi <- t(chol(oblFA$Phi))
    } else {
      L_Phi <- diag(1)
    }
    
    # Mean and sd of level-1 residual variance
    mu_sigma = colMeans(sigma)
    sd_sigma = apply(sigma, 2, sd)
    
    # Compute offsets and scale values
    scale_sigma = sd_sigma / sqrt(1 - 2/pi)
    shift_sigma = mu_sigma - scale_sigma * sqrt(2/pi)
    
    # Ensure positive shifts
    if(any(shift_sigma<0)){ shift_sigma[shift_sigma<0] <- 1e-3 }
    
    # ---------------------------------- # 
    #    Gaussian initial values list    #
    # ---------------------------------- # 
    
    inital_values <- list(
      # Population means
      mu_alpha = colMeans(alpha),
      mu_theta = colMeans(theta),
      shift_sigma = shift_sigma,
      # Population standard deviations
      sd_alpha = apply(alpha, 2, sd),
      sd_theta = apply(theta, 2, sd),
      scale_sigma = scale_sigma,
      # Cholesky factor decomposition
      L_R_Alpha = t(chol(cor(alpha))),
      L_Phi     = L_Phi,
      # Individual-level latent scores
      Alpha_tilde = t(apply(alpha, 2, function(x) (x - mean(x)) / sd(x))),
      Theta_tilde = t(apply(theta, 2, function(x) (x - mean(x)) / sd(x))),
      sigma_tilde = t(apply(sigma, 2, scale_half_std_normal)),
      # Communalities, unit vector and factor loadings
      h2 = rep(0.5, J),
      L_raw = matrix(c(oblFA$loadings), ncol = M),
      Z = z_normed
    )
  }
  if(model == "exgaussian"){
    # Ex-gaussian method of moments function
    # Based on retimes: https://github.com/cran/retimes/blob/master/R/timefit.R
    mexgauss <- function(x,n=length(x)) {
      k <- start <- c(mu=NaN,sigma=NaN,tau=NaN) # Momenti
      k[1] <- mean(x)
      xdev <- x-k[1]
      k[2] <- sum(xdev^2)/(n-1)
      k[3] <- sum(xdev^3)/(n-1)
      if(k[3] > 0)
        start[3] <- (k[3]/2)^(1/3)
      else
        start[3] <- 0.8*sd(x)
      start[2] <- sqrt(abs(k[2]-start[3]^2))
      start[1] <- k[1]-start[3]
      return(start)
    }
    
    # Compute ex-gaussian parameters
    alpha <- tapply(df$RT, list(df$subject, df$task), function(x) mexgauss(x)[1])
    pre_theta <- tapply(df$RT, list(df$subject, df$condition, df$task), function(x) mexgauss(x)[1])
    theta <- pre_theta[,2,] - pre_theta[,1,]
    sigma <- tapply(df$RT, list(df$subject, df$task), function(x) mexgauss(x)[2])
    tau   <- tapply(df$RT, list(df$subject, df$task), function(x) mexgauss(x)[3])
    
    # Exploratory factor analysis: varimax
    h2 <- psych::fa(theta, M, fm = "minrank", rotate = "varimax")$communality
    # Exploratory factor analysis: oblimin
    oblFA <- psych::fa(theta, M, fm = "minrank", rotate = "oblimin")
    if(M > 1){ 
      L_Phi <- t(chol(oblFA$Phi))
    } else {
      L_Phi <- diag(1)
    }
    
    # Mean and sd of level-1 residual variance
    mu_sigma = colMeans(sigma)
    sd_sigma = apply(sigma, 2, sd)
    mu_tau   = colMeans(tau)
    sd_tau   = apply(tau, 2, sd)
    mu_rate  = 1/mu_tau
    sd_rate  = apply(1/tau, 2, sd)
    
    # Compute offsets and scale values
    scale_sigma  = sd_sigma / sqrt(1 - 2/pi)
    shift_sigma = mu_sigma - scale_sigma * sqrt(2/pi)
    scale_tau    = sd_tau / sqrt(1 - 2/pi)
    shift_tau   = mu_tau - scale_tau * sqrt(2/pi)
    scale_rate   = sd_rate / sqrt(1 - 2/pi)
    shift_rate  = mu_rate - scale_rate * sqrt(2/pi)
    
    # Ensure positive shifts
    if(any(shift_sigma<0)){ shift_sigma[shift_sigma<0] <- 1e-3 }
    if(any(shift_tau<0)){ shift_tau[shift_tau<0] <- 1e-3 }
    if(any(shift_rate<0)){ shift_rate[shift_rate<0] <- 1e-3 }
    
    # ------------------------------------- # 
    #    Ex-Gaussian initial values list    #
    # ------------------------------------- # 
    
    inital_values <- list(
      # Population means
      mu_alpha = colMeans(alpha),
      mu_theta = colMeans(theta),
      shift_sigma = shift_sigma,
      shift_rate  = shift_rate,
      shift_tau  = shift_tau,
      # Population standard deviations
      sd_alpha = apply(alpha, 2, sd),
      sd_theta = apply(theta, 2, sd),
      scale_sigma = scale_sigma,
      scale_rate  = scale_rate,
      scale_tau   = scale_tau,
      # Cholesky factor decomposition
      L_R_Alpha = t(chol(cor(alpha))),
      L_Phi     = L_Phi,
      # Communalities, unit vector and factor loadings
      h2 = rep(0.5, J),
      L_raw = matrix(c(oblFA$loadings), ncol = M),
      Z = z_normed,
      # Individual-level latent scores
      Alpha_tilde = t(apply(alpha, 2, function(x) (x - mean(x)) / sd(x))),
      Theta_tilde = t(apply(theta, 2, function(x) (x - mean(x)) / sd(x))),
      sigma_tilde = t(apply(sigma, 2, scale_half_std_normal)),
      rate_tilde  = t(apply(1/tau, 2, scale_half_std_normal)),
      tau_tilde   = t(apply(tau, 2, scale_half_std_normal))
    )
  }
  if(model == "shifted.lognormal"){
    # Compute individual-level parameters
    min_RT <- tapply(df$RT, list(df$subject, df$task), min)
    delta <- matrix(
      truncnorm::rtruncnorm(length(unique(df$subject)) * length(unique(df$task)), 
                            a = 0, b = c(min_RT), mean = c(min_RT/1.5), sd = .02), 
      nrow = length(unique(df$subject)), ncol = length(unique(df$task))
    )
    # Compute individual-level parameters
    alpha <- theta <- sigma <- matrix(NA, ncol = length(unique(df$task)), 
                                      nrow = length(unique(df$subject)))
    for(i in 1:length(unique(df$subject))){
      for(j in 1:length(unique(df$task))){
        alpha[i,j] <- mean(log(df$RT[which(df$subject==i & df$task==j)] - delta[i,j]))
        theta[i,j] <- mean(log(df$RT[which(df$subject==i & df$task==j & df$condition == 2)] - delta[i,j])) - 
          mean(log(df$RT[which(df$subject==i & df$task==j & df$condition == 1)] - delta[i,j]))
        sigma[i,j] <- sd(log(df$RT[which(df$subject==i & df$task==j)] - delta[i,j]))
      }
    }
    
    # Exploratory factor analysis: varimax
    h2 <- psych::fa(theta, M, fm = "minrank", rotate = "varimax")$communality
    # Exploratory factor analysis: oblimin
    oblFA <- psych::fa(theta, M, fm = "minrank", rotate = "oblimin")
    if(M > 1){ 
      L_Phi <- t(chol(oblFA$Phi))
    } else {
      L_Phi <- diag(1)
    }
    
    # Mean and sd of level-1 residual variance
    mu_sigma = colMeans(sigma)
    sd_sigma = apply(sigma, 2, sd)
    
    # Compute offsets and scale values
    scale_sigma = sd_sigma / sqrt(1 - 2/pi)
    shift_sigma = mu_sigma - scale_sigma * sqrt(2/pi)
    
    # Ensure positive shifts
    if(any(shift_sigma<0)){ shift_sigma[shift_sigma<0] <- 1e-3 }
    
    # ------------------------------------------- # 
    #    Shifted-lognormal initial values list    #
    # ------------------------------------------- # 
    
    inital_values <- list(
      # Population means
      mu_alpha = colMeans(alpha),
      mu_theta = colMeans(theta),
      shift_sigma = shift_sigma,
      mu_delta = colMeans(delta),
      # Population standard deviations
      sd_alpha = apply(alpha, 2, sd),
      sd_theta = apply(theta, 2, sd),
      scale_sigma = scale_sigma,
      sd_delta = apply(delta, 2, sd),
      # Cholesky factor decomposition
      L_R_Alpha = t(chol(cor(alpha))),
      L_Phi     = L_Phi,
      # Individual-level latent scores
      Alpha_tilde = t(apply(alpha, 2, FUN = function(x) (x - mean(x)) / sd(x))),
      Theta_tilde = t(apply(theta, 2, FUN = function(x) (x - mean(x)) / sd(x))),
      sigma_tilde = t(apply(sigma, 2, scale_half_std_normal)),
      Delta       = delta,
      # Communalities, unit vector and factor loadings
      h2 = rep(0.5, J),
      L_raw = matrix(c(oblFA$loadings), ncol = M),
      Z = z_normed
    )
  }
  
  return(inital_values)
}

# Posterior distribution summaries
posterior_summaries <- function(fit, parameter, empty = FALSE, times, rotated = FALSE, 
                                cormat = FALSE, CFA_loadings = FALSE, CFA_phi = FALSE, 
                                parnames = NULL){
  
  # Same parameter name as input parameter
  if(is.null(parnames)) parnames <- parameter
  
  # Empty vector
  if(empty){
    res <- c(
      postmean = rep(999, times),
      postsd   = rep(999, times),
      LL       = rep(999, times),
      UL       = rep(999, times),
      rhat     = 999,
      essb     = 999,
      esst     = 999
    )
    
    # Results for factor loadings in EFA settings
  } else if(rotated){
    res <- c(
      postmean = fit$mean,
      postsd   = fit$sd,
      LL       = fit$postLL,
      UL       = fit$postUL,
      rhat     = mean(fit$rhat),
      essb     = mean(fit$ess_bulk),
      esst     = mean(fit$ess_tail)
    )
    
  } else {
    # Cache the summary and quantile summary once.
    sumFit   <- fit$summary(parameter)
    sumFit_q <- fit$summary(parameter, quantiles = ~ quantile2(., probs = c(0.025, 0.975)))
    
    # Latent correlation matrix
    if(cormat){
      # Hard-coded dimension: 6x6 matrix assumed.
      mat_mean <- matrix(sumFit$mean, ncol = 6, nrow = 6)
      mat_sd   <- matrix(sumFit$sd, ncol = 6, nrow = 6)
      mat_LL   <- matrix(unlist(sumFit_q[,2]), ncol = 6, nrow = 6)
      mat_UL   <- matrix(unlist(sumFit_q[,3]), ncol = 6, nrow = 6)
      # Extract lower triangle values.
      res <- c(
        postmean = mat_mean[lower.tri(diag(6))],
        postsd   = mat_sd[lower.tri(diag(6))],
        LL       = mat_LL[lower.tri(diag(6))],
        UL       = mat_UL[lower.tri(diag(6))],
        rhat     = mean(matrix(sumFit$rhat, ncol = 6, nrow = 6)[lower.tri(diag(6))]),
        essb     = mean(matrix(sumFit$ess_bulk, ncol = 6, nrow = 6)[lower.tri(diag(6))]),
        esst     = mean(matrix(sumFit$ess_tail, ncol = 6, nrow = 6)[lower.tri(diag(6))])
      )
      
      # Results for factor loadings in CFA settings
    } else if(CFA_loadings){
      idx <- c(1:3, 10:12)
      res <- c(
        postmean = sumFit$mean[idx],
        postsd   = sumFit$sd[idx],
        LL       = unlist(sumFit_q[idx, 2], use.names = FALSE),
        UL       = unlist(sumFit_q[idx, 3], use.names = FALSE),
        rhat     = mean(sumFit$rhat[idx]),
        essb     = mean(sumFit$ess_bulk[idx]),
        esst     = mean(sumFit$ess_tail[idx])
      )
      
      # Results for latent common-factos correlation in CFA settings
    } else if(CFA_phi){
      idx <- 2
      res <- c(
        postmean = sumFit$mean[idx],
        postsd   = sumFit$sd[idx],
        LL       = unlist(sumFit_q[idx, 2], use.names = FALSE),
        UL       = unlist(sumFit_q[idx, 3], use.names = FALSE),
        rhat     = sumFit$rhat[idx],
        essb     = sumFit$ess_bulk[idx],
        esst     = sumFit$ess_tail[idx]
      )
      # Results for the rest of parameters
    } else {
      res <- c(
        postmean = sumFit$mean,
        postsd   = sumFit$sd,
        LL       = unlist(sumFit_q[, 2], use.names = FALSE),
        UL       = unlist(sumFit_q[, 3], use.names = FALSE),
        rhat     = mean(sumFit$rhat),
        essb     = mean(sumFit$ess_bulk),
        esst     = mean(sumFit$ess_tail)
      )
    }
  }
  # Adjust output names
  names(res) <- paste0(parnames, "_", names(res))
  res
}

# Global summaries for each condition
global_summaries <- function(results, parameter_name, parameter, Rho = FALSE, 
                             empty = FALSE, parnames = NULL, method_of_moments = FALSE, 
                             Phi.lv = FALSE) {
  
  # If empty, return a vector of NA values with the appropriate names.
  if (empty) {
    if(Phi.lv) {
      est_names <- c("abias_avg", "arbias_avg", "stdbias_avg", "rmse_avg", 
                     "postsd_avg", "postsd_disp", "empirical_disp", "rhat_avg", 
                     "rhat_disp", "essb_avg", "essb_disp", "esst_avg", "esst_disp",
                     "ECR_avg", "power_avg")
      res <- setNames(rep(NA, length(est_names)), paste0(est_names, "_", parameter_name))
      return(res)
    } else {
      est_names <- c("abias_avg", "abias_disp", "arbias_avg", "arbias_disp", 
                     "stdbias_avg", "stdbias_disp", "rmse_avg", "rmse_disp", 
                     "postsd_avg", "postsd_disp", "empirical_disp", "rhat_avg", 
                     "rhat_disp", "essb_avg", "essb_disp", "esst_avg", "esst_disp", 
                     "ECR_avg", "ECR_disp", "power_avg", "power_disp")
      res <- setNames(rep(NA, length(est_names)), paste0(est_names, "_", parameter_name))
      return(res)
    }
  }
  
  # Helper function to retrieve column indices based on a pattern.
  get_cols <- function(pattern) {
    grep(paste0(parameter_name, pattern), colnames(results))
  }
  
  # --- Branch 1: Method of Moments ---
  if (method_of_moments) {
    tmp.data   <- results[, get_cols("_est")]
    tmp.data.se <- results[, get_cols("_se")]
    param_lower <- parameter[lower.tri(parameter)]
    
    # Compute bias measures only once per type.
    bias_abs     <- bias(tmp.data, parameter = param_lower, abs = TRUE)
    bias_abs_rel <- bias(tmp.data, parameter = param_lower, type = "abs_relative", percent = TRUE, abs = TRUE)
    bias_std     <- bias(tmp.data, parameter = param_lower, type = "standardized", abs = TRUE)
    
    # Average and dispersion in all statistics
    abias_avg   <- mean(bias_abs)
    abias_disp  <- sd(bias_abs)
    arbias_avg  <- mean(bias_abs_rel)
    arbias_disp <- sd(bias_abs_rel)
    stdbias_avg <- mean(bias_std)
    stdbias_disp<- sd(bias_std)
    
    # Root Mean Square Error
    rmse_vals <- apply(tmp.data, 1, function(x) Rho_RMSE(param_lower, x))
    rmse_avg  <- mean(rmse_vals)
    rmse_disp <- sd(rmse_vals)
    
    # Standardized Root Mean Residuals
    srmr_vals <- apply(tmp.data, 1, function(x) SRMR(parameter, lower.tri2mat(lower_vals = x)))
    srmr_avg  <- mean(srmr_vals)
    srmr_disp <- sd(srmr_vals)
    
    # Posterior standard errors
    SEs_avg  <- mean(colMeans(tmp.data.se))
    SEs_disp <- mean(apply(tmp.data.se, 2, sd))
    
    # Efficiency estimate
    empirical_disp <- mean(apply(tmp.data, 2, sd))
    
    # Coverage and power estimates
    LL_cols <- get_cols("_LL")
    UL_cols <- get_cols("_UL")
    temp    <- abind::abind(results[, LL_cols], results[, UL_cols], along = 3)
    resECR  <- sapply(1:ncol(temp), function(j) SimDesign::ECR(cbind(temp[, j, 1], temp[, j, 2]), param_lower[j]))
    power   <- sapply(1:ncol(temp), function(j) 1 - SimDesign::ECR(cbind(temp[, j, 1], temp[, j, 2]), 0))
    
    # Final results
    res <- c(abias_avg, abias_disp, arbias_avg, arbias_disp, stdbias_avg, stdbias_disp, 
             rmse_avg, rmse_disp, srmr_avg, srmr_disp, SEs_avg, SEs_disp, 
             empirical_disp, mean(resECR), sd(resECR), mean(power), sd(power))
    
    # Modify names
    names(res) <- paste0(c("abias_avg", "abias_disp", "arbias_avg", "arbias_disp", 
                           "stdbias_avg", "stdbias_disp", "rmse_avg", "rmse_disp", 
                           "srmr_avg", "srmr_disp", "postsd_avg", "postsd_disp",
                           "empirical_disp", "ECR_avg", "ECR_disp", "power_avg", "power_disp"), "_", parameter_name)
    return(res)
  }
  
  # --- Common Data Extraction for Remaining Branches ---
  tmp.data      <- results[, get_cols("_postmean")]
  tmp.data.pstsd<- results[, get_cols("_postsd")]
  
  # --- Branch 2: When Rho is TRUE ---
  if (Rho) {
    # Select nonredundant correlation values
    param_lower <- parameter[lower.tri(parameter)]
    
    # Compute bias measures only once per type.
    bias_abs     <- bias(tmp.data, parameter = param_lower, abs = TRUE)
    bias_abs_rel <- bias(tmp.data, parameter = param_lower, type = "abs_relative", percent = TRUE, abs = TRUE)
    bias_std     <- bias(tmp.data, parameter = param_lower, type = "standardized", abs = TRUE)
    
    # Average and dispersion in all statistics
    abias_avg   <- mean(bias_abs)
    abias_disp  <- sd(bias_abs)
    arbias_avg  <- mean(bias_abs_rel)
    arbias_disp <- sd(bias_abs_rel)
    stdbias_avg <- mean(bias_std)
    stdbias_disp<- sd(bias_std)
    
    # Root Mean Square Error
    rmse_vals <- apply(tmp.data, 1, function(x) Rho_RMSE(param_lower, x))
    rmse_avg  <- mean(rmse_vals)
    rmse_disp <- sd(rmse_vals)
    
    # Standardized Root Mean Residuals
    srmr_vals <- apply(tmp.data, 1, function(x) SRMR(parameter, lower.tri2mat(lower_vals = x)))
    srmr_avg  <- mean(srmr_vals)
    srmr_disp <- sd(srmr_vals)
    
    # Posterior standard errors
    postsd_avg <- mean(colMeans(tmp.data.pstsd))
    postsd_disp<- mean(apply(tmp.data.pstsd, 2, sd))
    
    # Efficiency and convergence
    empirical_disp <- mean(apply(tmp.data, 2, sd))
    rhat   <- results[, get_cols("_rhat")]
    rhat_avg  <- mean(rhat)
    rhat_disp <- sd(rhat)
    essb   <- results[, get_cols("_essb")]
    essb_avg  <- mean(essb)
    essb_disp <- sd(essb)
    esst   <- results[, get_cols("_esst")]
    esst_avg  <- mean(esst)
    esst_disp <- sd(esst)
    
    # Coverage and power estimates
    LL_cols <- get_cols("_LL")
    UL_cols <- get_cols("_UL")
    temp    <- abind::abind(results[, LL_cols], results[, UL_cols], along = 3)
    resECR  <- sapply(1:ncol(temp), function(j) SimDesign::ECR(cbind(temp[, j, 1], temp[, j, 2]), param_lower[j]))
    power   <- sapply(1:ncol(temp), function(j) 1 - SimDesign::ECR(cbind(temp[, j, 1], temp[, j, 2]), 0))
    
    # Final results
    res <- c(abias_avg, abias_disp, arbias_avg, arbias_disp, stdbias_avg, stdbias_disp, 
             rmse_avg, rmse_disp, srmr_avg, srmr_disp, postsd_avg, postsd_disp, empirical_disp, 
             rhat_avg, rhat_disp, essb_avg, essb_disp, esst_avg, esst_disp,
             mean(resECR), sd(resECR), mean(power), sd(power))
    
    # Modify names
    names(res) <- paste0(c("abias_avg", "abias_disp", "arbias_avg", "arbias_disp", 
                           "stdbias_avg", "stdbias_disp", "rmse_avg", "rmse_disp", 
                           "srmr_avg", "srmr_disp", "postsd_avg", "postsd_disp",
                           "empirical_disp", "rhat_avg", "rhat_disp", "essb_avg", 
                           "essb_disp", "esst_avg", "esst_disp", "ECR_avg", "ECR_disp", 
                           "power_avg", "power_disp"), "_", parameter_name)
    return(res)
    
  } else if (Phi.lv){
    # Phi it's only one parameter
    abias_avg     <- bias(tmp.data, parameter = parameter, abs = TRUE)
    arbias_avg <- bias(tmp.data, parameter = parameter, abs = TRUE, type = "abs_relative", percent = TRUE)
    stdbias_avg     <- bias(tmp.data, parameter = parameter, type = "standardized", abs = TRUE)
    rmse_avg <- RMSE(tmp.data, parameter = parameter)
    # Posterior standard errors
    postsd_avg <- mean(tmp.data.pstsd)
    postsd_disp<- sd(tmp.data.pstsd)
    # Efficiency and convergence
    empirical_disp <- sd(tmp.data)
    rhat   <- results[, get_cols("_rhat")]
    rhat_avg  <- mean(rhat)
    rhat_disp <- sd(rhat)
    essb   <- results[, get_cols("_essb")]
    essb_avg  <- mean(essb)
    essb_disp <- sd(essb)
    esst   <- results[, get_cols("_esst")]
    esst_avg  <- mean(esst)
    esst_disp <- sd(esst)
    # Coverage and power estimates
    LL_cols <- get_cols("_LL")
    UL_cols <- get_cols("_UL")
    temp    <- cbind(results[, LL_cols], results[, UL_cols])
    resECR  <- ECR(CIs = temp, parameter)
    power   <- 1 - ECR(CIs = temp, parameter = 0)
    
    # Final results
    res <- c(abias_avg, arbias_avg, stdbias_avg, rmse_avg, postsd_avg, postsd_disp, empirical_disp, 
             rhat_avg, rhat_disp, essb_avg, essb_disp, esst_avg, esst_disp,
             resECR, power)
    
    # Modify names
    names(res) <- paste0(c("abias_avg", "arbias_avg", "stdbias_avg", "rmse_avg", 
                           "postsd_avg", "postsd_disp", "empirical_disp", "rhat_avg", 
                           "rhat_disp", "essb_avg", "essb_disp", "esst_avg", "esst_disp",
                           "ECR_avg", "power_avg"), "_", parameter_name)
    return(res)
  } else {
    # --- Branch 3: When Rho is FALSE ---
    # Compute bias measures only once per type.
    bias_abs     <- bias(tmp.data, parameter = parameter, abs = TRUE)
    bias_abs_rel <- bias(tmp.data, parameter = parameter, abs = TRUE, type = "abs_relative", percent = TRUE)
    bias_std     <- bias(tmp.data, parameter = parameter, type = "standardized", abs = TRUE)
    
    # Average and dispersion in all statistics
    abias_avg   <- mean(bias_abs)
    abias_disp  <- sd(bias_abs)
    arbias_avg  <- mean(bias_abs_rel)
    arbias_disp <- sd(bias_abs_rel)
    stdbias_avg <- mean(bias_std)
    stdbias_disp<- sd(bias_std)
    
    # Root Mean Square Error
    rmse_vals <- RMSE(tmp.data, parameter = parameter)
    rmse_avg  <- mean(rmse_vals)
    rmse_disp <- sd(rmse_vals)
    
    # Posterior standard errors
    postsd_avg <- mean(colMeans(tmp.data.pstsd))
    postsd_disp<- mean(apply(tmp.data.pstsd, 2, sd))
    
    # Efficiency and convergence
    empirical_disp <- mean(apply(tmp.data, 2, sd))
    rhat   <- results[, get_cols("_rhat")]
    rhat_avg  <- mean(rhat)
    rhat_disp <- sd(rhat)
    essb   <- results[, get_cols("_essb")]
    essb_avg  <- mean(essb)
    essb_disp <- sd(essb)
    esst   <- results[, get_cols("_esst")]
    esst_avg  <- mean(esst)
    esst_disp <- sd(esst)
    
    # Coverage and power estimates
    LL_cols <- get_cols("_LL")
    UL_cols <- get_cols("_UL")
    temp    <- abind::abind(results[, LL_cols], results[, UL_cols], along = 3)
    resECR  <- sapply(1:ncol(temp), function(j) SimDesign::ECR(cbind(temp[, j, 1], temp[, j, 2]), parameter[j]))
    power   <- sapply(1:ncol(temp), function(j) 1 - SimDesign::ECR(cbind(temp[, j, 1], temp[, j, 2]), 0))
    
    # Final results
    res <- c(abias_avg, abias_disp, arbias_avg, arbias_disp, stdbias_avg, stdbias_disp,
             rmse_avg, rmse_disp, postsd_avg, postsd_disp, empirical_disp, 
             rhat_avg, rhat_disp, essb_avg, essb_disp, esst_avg, esst_disp,
             mean(resECR), sd(resECR), mean(power), sd(power))
    
    # Modify names
    names(res) <- paste0(c("abias_avg", "abias_disp", "arbias_avg", "arbias_disp", 
                           "stdbias_avg", "stdbias_disp", "rmse_avg", "rmse_disp", 
                           "postsd_avg", "postsd_disp", "empirical_disp", "rhat_avg", 
                           "rhat_disp", "essb_avg", "essb_disp", "esst_avg", "esst_disp",
                           "ECR_avg", "ECR_disp", "power_avg", "power_disp"), "_", parameter_name)
    return(res)
  }
}

# MatchAlign Algorithm: orthogonal factor rotation across draws
# Paper: https://arxiv.org/abs/2107.13783
# Based on: https://github.com/poworoznek/infinitefactor/blob/master/R/jointRot.R
MatchAlign <- function(lambda_draws){
  # Ensure infinitefactor package
  ensure_packages("infinitefactor")
  # Save array dimensions
  dims <- dim(lambda_draws)
  # Apply varimax rotation
  vari_all <- apply(lambda_draws, 1, varimax)
  # Select rotated loadings
  loads <- lapply(vari_all, `[[`, 1)
  # Apply singular value decomposition and select largest singular value
  norms <- sapply(loads, norm, "2")
  # Use singular values to select pivot matrix
  piv <- loads[order(norms)][[round(dims[1]/2)]]
  # Apply transformation based on pivot matrix
  matches <- lapply(loads, infinitefactor::msfOUT, piv)
  rotated_lambdas <- mapply(infinitefactor::aplr, loads, matches, SIMPLIFY = FALSE)
  # Save in array
  array_lambdas <- array(unlist(rotated_lambdas), dim = c(dims[2], dims[3], dims[1]))
  # Return arrays
  return(array_lambdas)
}

# Rotated EFA factor loadings
rotated_EFA <- function(fit, lambda_name, probs = c(0.025, 0.975), 
                        only.summary = TRUE){
  # Ensure infinitefactor package
  ensure_packages(c("posterior", "tibble"))
  # Extract draws
  lambda_draws <- draws_of(as_draws_rvars(fit)[[lambda_name]])
  # Number of chains
  nchains <- fit$num_chains()
  # Estimate rotated factor loadings
  rotated <- MatchAlign(lambda_draws)
  niter <- dim(rotated)[3]/nchains
  # Compute sumamry table
  postMean   <- apply(rotated, 1:2, mean)
  postMedian <- apply(rotated, 1:2, median)
  postSD     <- apply(rotated, 1:2, sd)
  postMad    <- apply(rotated, 1:2, mad)
  postLL     <- apply(rotated, 1:2, quantile, probs = probs[1])
  postUL     <- apply(rotated, 1:2, quantile, probs = probs[2])
  
  # Estimate the intra-chain variability index
  temp <- matrix(0, niter, nchains)
  Rhat <- ESS_tail <- ESS_bulk <- vector(length = nrow(rotated) * ncol(rotated))
  ind <- 0
  for(i in 1:nrow(rotated)){
    for(j in 1:ncol(rotated)){
      ind <- ind + 1
      for(c in 1:(nchains)-1){
        temp[,(c + 1)] <- rotated[i, j, (niter * c + 1):(niter * (c + 1))]
      }
      Rhat[ind] <- posterior::rhat(temp)
      ESS_tail[ind] <- posterior::ess_tail(temp)
      ESS_bulk[ind] <- posterior::ess_bulk(temp)
    }
  }
  
  # Create a summary table
  summres <- dplyr::tibble(variable = paste0("Lambda[", 1:nrow(rotated), ",", rep(1:ncol(rotated), each = nrow(rotated)) , "]"),
                           mean = c(postMean), median = c(postMedian), sd = c(postSD), mad = c(postMad), postLL = c(postLL), 
                           postUL = c(postUL), rhat = Rhat, ess_bulk = ESS_bulk, ess_tail = ESS_tail)
  
  # Return summary results
  if(only.summary){ 
    return(summres)
  } else {
    return(list(summary = summres, rotated_draws = rotated))
  }
}

# Align posterior means with target population lambda matrix
simulation_alignment <- function(fit, target){
  # Ensure fungible and package
  ensure_packages("fungible")
  # Apply MathAlign method
  rotated_draws <- rotated_EFA(fit, lambda_name = "Lambda_std", 
                               probs = c(.025, .975), only.summary = TRUE)
  # Align loadings: Least Squares and 
  LS_Align <- faAlign(F1 = target, F2 = matrix(rotated_draws$mean, 6, 2))
  if(!LS_Align$UniqueMatch){
    # Try congruency coefficients
    CC_Align <- faAlign(F1 = target, F2 = matrix(rotated_draws$mean, 6, 2), MatchMethod = "CC")
    if(!CC_Align$UniqueMatch){
      stop("Non-unique rotation match with faAlign")
    } else {
      aligned_draws <- CC_Align$F2
    }
  } else {
    aligned_draws <- LS_Align$F2
  }
  
  # Align rotated output
  original_order <- rotated_draws$mean   # Original posterior mean values
  aligned_order <- c(aligned_draws)      # Aligned vector
  
  # Reorder both columns without the sign
  reorder_index <- match(abs(aligned_order), abs(original_order))
  
  # Identify Sign changes
  sign_change <- sign(aligned_order) != sign(original_order[reorder_index])
  
  # Reorder and adjust posterior outcome
  summary_aligned <- rotated_draws[reorder_index, ] %>%
    # Save original lower and upper limits
    mutate(
      oldLL = postLL,
      oldUL = postUL
    ) %>%
    mutate(
      mean    = ifelse(sign_change, -mean, mean),
      median  = ifelse(sign_change, -median, median),
      # If sign changes, postLL its -oldUL and postUL its -oldLL
      postLL  = ifelse(sign_change, -oldUL, oldLL),
      postUL  = ifelse(sign_change, -oldLL, oldUL)
    ) %>%
    # Remove temporal colums
    dplyr::select(-oldLL, -oldUL)
  
  # Return aligned summary
  return(summary_aligned)
}

# From lower triangular matrix to full matrix
lower.tri2mat <- function(lower_vals) {
  # Number of lower-tri off-diagonal elements
  m <- length(lower_vals)
  
  # Solve for n in the equation n*(n - 1)/2 = m
  # This is just the quadratic formula for n^2 - n - 2m = 0
  n <- (1 + sqrt(1 + 8 * m)) / 2
  
  # Check if n is an integer
  if (abs(n - round(n)) > 1e-10) {
    stop("Length of 'lower_vals' does not correspond to a valid (integer) n.")
  }
  n <- as.integer(round(n))
  
  # Initialize an n-by-n matrix with 1 on the diagonal
  R <- diag(1, n)
  
  # Fill in the lower-triangle off-diagonal entries
  R[lower.tri(R)] <- lower_vals
  
  # Mirror them to the upper triangle
  R[upper.tri(R)] <- t(R)[upper.tri(R)]
  
  return(R)
}

# Standardized Root Mean Residuals (SRMR, a normalized frobenius norm)
SRMR <- function(population, estimated) {
  
  # Extract only the unique elements (upper triangle including diagonal)
  idx <- upper.tri(population, diag = TRUE)
  
  # Compute the srmr statistic
  return(sqrt(mean((population[idx] - estimated[idx])^2)))
}

# Average RMSE for all correlation matrix
Rho_RMSE <- function(lower.tri.pop, lower.tri.est) { 
  sqrt(mean((lower.tri.pop - lower.tri.est)^2)) 
}

# Method-of-moments experimental effect
# Based on: https://github.com/PerceptionAndCognitionLab/ctx-inhibition/blob/public/papers/rev3/lib.R
analytical.effect <- function(data, subject_var = "subject", task_var = "task", 
                              condition_var = "condition", rt_var = "RT"){
  # Estimate response time mean per subject, taska and condition
  mean.rt.IJK <- tapply(data[,rt_var], INDEX = list(data[,subject_var], 
                                                    data[,task_var], 
                                                    data[,condition_var]), FUN = mean)
  # Compute experimental effect via method of moments
  effect <- mean.rt.IJK[,,2] - mean.rt.IJK[,,1]
  return(effect)
}

# Spearman correction with SEs via method of moments
# Based on: https://github.com/PerceptionAndCognitionLab/ctx-inhibition/blob/public/papers/rev3/lib.R
spearman.adj <- function(data, subject_var = "subject", task_var = "task", 
                         condition_var = "condition", rt_var = "RT"){
  # Estimate response time mean per subject, taska and condition
  mean.rt.IJK <- tapply(data[,rt_var], INDEX = list(data[,subject_var], 
                                                    data[,task_var], 
                                                    data[,condition_var]), FUN = mean)
  
  # Compute experimental effect via method of moments
  effects.mm <- mean.rt.IJK[,,2] - mean.rt.IJK[,,1]
  
  # Standard error of the mean (i.e., effects residual-variance)
  se <- function(x) sqrt(var(x)/length(x))
  
  # Residual variances of experimental effect per participant, task and condition
  effects.SEs <- tapply(data[,rt_var], INDEX = list(data[,subject_var], 
                                                    data[,task_var], 
                                                    data[,condition_var]), FUN = se)
  
  # Experimental effects residual variances (adding residual variances of each condition and task)
  effects.resvar <- colMeans(effects.SEs[,,1]^2 + effects.SEs[,,2]^2)
  
  return(list(rho_est = cov2cor(cov(effects.mm) - diag(effects.resvar)),
              reliability = (diag(cov(effects.mm)) - effects.resvar)/ diag(cov(effects.mm))))
}

# Pareto smoothed importance sampling leave-one-out cross-validation for all models
loo_BHGFM <- function(fit, model, data, n.cores = 2) {
  
  # Ensure {loo} package
  ensure_packages("loo")
  
  # Conditional log-likelihood functions in R: Gaussian model
  llfun_gaussian <- function(data_i, draws, log = TRUE){
    # Extract relevant data for the i-th observation
    RT_i <- data_i$RT
    condition <- data_i$condition
    task <- data_i$task
    subject <- data_i$subject
    
    # Select
    colpars <- paste0(c("Alpha", "Theta", "sigma"),"[",subject, ",", task, "]")
    
    # Compute the linear predictor
    mu_ijk <- draws[,colpars[1]] + (condition - 1.5) * draws[,colpars[2]]
    
    # Calculate the log-likelihood for each posterior draw
    log_lik_values <- dnorm(RT_i, mean = mu_ijk, sd = draws[,colpars[3]], log = TRUE)
    
    return(log_lik_values)
  }
  
  # Conditional log-likelihood functions in R: Exgaussian model
  llfun_exgaussian <- function(data_i, draws, log = TRUE){
    
    # Ex-gaussian log-probability density function
    dexgauss <- function(x, mu, sigma, lambda, log = FALSE) {
      # Input control
      if (any(sigma <= 0)) {
        stop("sigma must be positive.")
      }
      if (any(lambda <= 0)) {
        stop("lambda must be positive.")
      }
      
      # Complementary error function (erfc)
      erfc <- function(z) 2 * pnorm(-z * sqrt(2))
      
      # Log-probability density function (log-PDF)
      log_pdf <- log(lambda / 2) +
        0.5 * lambda * (2 * mu + lambda * sigma^2 - 2 * x) +
        log(erfc((mu + lambda * sigma^2 - x) / (sqrt(2) * sigma)))
      
      # If log = FALSE, return the density (PDF), otherwise return log-PDF
      if (log) {
        return(log_pdf)
      } else {
        return(exp(log_pdf))
      }
    }
    
    # Extract relevant data for the i-th observation
    RT_i <- data_i$RT
    condition <- data_i$condition
    task <- data_i$task
    subject <- data_i$subject
    
    # Select
    colpars <- paste0(c("Alpha", "Theta", "sigma", "rate"),"[",subject,",", task, "]")
    
    # Compute the linear predictor
    mu_ijk <- draws[,colpars[1]] + (condition - 1.5) * draws[,colpars[2]]
    
    # Calculate the log-likelihood for each posterior draw
    log_lik_values <- dexgauss(x = RT_i, mu = mu_ijk, sigma = draws[,colpars[3]], 
                               lambda = draws[,colpars[4]], log = TRUE)
    
    return(log_lik_values)
  }
  
  # Conditional log-likelihood functions in R: Shifted-lognormal model
  llfun_shiftlognormal <- function(data_i, draws, log = TRUE){
    
    # Ex-gaussian log-probability density function
    dshifted_lognorm <- function(x, mu, sigma, delta, log = FALSE) {
      # Input control
      if (any(sigma <= 0)) {
        stop("sigma must be positive.")
      }
      if (any(x < delta)) {
        stop("x must be greater than delta.")
      }
      
      # Log-probability density function (log-PDF)
      log_pdf <- -log(x - delta) - log(sigma * sqrt(2 * pi)) - 
                 ((log(x - delta) - mu)^2 / (2 * sigma^2))
      
      # Return log-PDF if log = TRUE, otherwise return the PDF
      if (log) {
        return(log_pdf)
      } else {
        return(exp(log_pdf))
      }
    }
    
    # Extract relevant data for the i-th observation
    RT_i <- data_i$RT
    condition <- data_i$condition
    task <- data_i$task
    subject <- data_i$subject
    
    # Select
    colpars <- paste0(c("Alpha", "Theta", "sigma", "Delta"),"[",subject, ",", task, "]")
    
    # Compute the linear predictor
    mu_ijk <- draws[,colpars[1]] + (condition - 1.5) * draws[,colpars[2]]
    
    # Calculate the log-likelihood for each posterior draw
    log_lik_values <- dshifted_lognorm(x = RT_i, mu = mu_ijk, sigma = draws[,colpars[3]], 
                                       delta = draws[,colpars[4]], log = TRUE)
    
    return(log_lik_values)
  }
  
  # All conditional log-likelihood functions
  ll_functions <- list(
    "gaussian" = llfun_gaussian,
    "exgaussian" = llfun_exgaussian,
    "shifted.lognormal" = llfun_shiftlognormal
  )
  
  # Parameters by model
  model_parameters <- list(
    "gaussian" = c("Alpha", "Theta", "sigma"),
    "exgaussian" = c("Alpha", "Theta", "sigma", "rate"),
    "shifted.lognormal" = c("Alpha", "Theta", "sigma", "Delta")
  )
  
  # Posterior draws per model
  posterior_draws <- posterior::subset_draws(fit$draws(model_parameters[[model]], format = "matrix"))
  
  # Compute PSIS-LOO 
  loo_model <- loo(ll_functions[[model]], draws = posterior_draws, data = data, cores = n.cores)
  
  # Return LOO object
  return(loo_model)
}

# Population parameters per condition
cond_parameters <- function(condition, fixed_objects) {
  # Experiment sizes
  I <- condition$I; J <- 6; K <- 2; L <- 100
  
  # Latent structure (lambda, Phi)
  if (condition$M == "uni") {
    lambda_std <- matrix(rep(condition$std_lambda_j, J), ncol = 1)
    Phi_cor    <- diag(1)
  } else if (condition$M == "2_CFA") {
    lambda_std <- matrix(
      c(rep(condition$std_lambda_j, 3), rep(0, 3),
        rep(0, 3), rep(condition$std_lambda_j, 3)),
      ncol = 2
    )
    Phi_cor <- diag(1 - 0.4, ncol = 2, nrow = 2) + 0.4
  } else if (condition$M == "2_EFA") {
    lambda_std <- matrix(
      c(rep(condition$std_lambda_j, 3), rep(condition$std_lambda_j/3, 3),
        rep(condition$std_lambda_j/3, 3), rep(condition$std_lambda_j, 3)),
      ncol = 2
    )
    Phi_cor <- diag(1, ncol = 2, nrow = 2)
  } else {
    stop("condition$M must be 'uni', '2_CFA', or '2_EFA'.")
  }
  
  # Initialize all vectors to NA (keys always exist)
  mu_alpha_j    <- rep(NA, J)
  sd_alpha_j    <- rep(NA, J)
  mu_theta_j    <- rep(NA, J)
  
  shift_sigma_j <- rep(NA, J)
  scale_sigma_j <- rep(NA, J)
  mu_sigma_j    <- rep(NA, J)
  sd_sigma_j    <- rep(NA, J)
  
  shift_tau_j   <- rep(NA, J)
  scale_tau_j   <- rep(NA, J)
  mu_tau_j      <- rep(NA, J)
  sd_tau_j      <- rep(NA, J)
  
  mu_delta_j    <- rep(NA, J)
  sd_delta_j    <- rep(NA, J)
  
  reliab_j      <- rep(NA, J)
  gamma2_j      <- rep(NA, J)
  sd_theta_j    <- rep(NA, J)
  
  # Select meta-parameter block
  if (identical(condition$Model, "exgaussian")) {
    meta_pars <- "meta_exgaussian_tau"
  } else if (identical(condition$Model, "shifted.lognormal")) {
    meta_pars <- "meta_shlognormal"
  } else {
    stop("condition$Model must be 'exgaussian' or 'shifted.lognormal'.")
  }
  
  # Common parameters across models
  mu_alpha_j    <- rep(fixed_objects$meta_parameters[[meta_pars]]$mu_alpha,   length.out = J)
  mu_theta_j    <- rep(fixed_objects$meta_parameters[[meta_pars]]$mu_theta,   length.out = J)
  sd_alpha_j    <- rep(fixed_objects$meta_parameters[[meta_pars]]$sd_alpha,   length.out = J)
  shift_sigma_j <- rep(fixed_objects$meta_parameters[[meta_pars]]$shift_sigma,length.out = J)
  scale_sigma_j <- rep(fixed_objects$meta_parameters[[meta_pars]]$scale_sigma,length.out = J)
  
  # Model-specific parameters
  if (identical(condition$Model, "exgaussian")) {
    shift_tau_j <- rep(fixed_objects$meta_parameters[[meta_pars]]$shift_tau, length.out = J)
    scale_tau_j <- rep(fixed_objects$meta_parameters[[meta_pars]]$scale_tau, length.out = J)
  } else { # shifted.lognormal
    mu_delta_j  <- rep(fixed_objects$meta_parameters[[meta_pars]]$mu_delta,  length.out = J)
    sd_delta_j  <- rep(fixed_objects$meta_parameters[[meta_pars]]$sd_delta,  length.out = J)
  }
  
  # Derived quantities
  reliab_j <- rep(condition$reliab, J)
  
  mu_sigma_j <- shift_sigma_j + scale_sigma_j * sqrt(2 / pi)
  sd_sigma_j <- scale_sigma_j * sqrt(1 - 2 / pi)
  
  if (identical(condition$Model, "exgaussian")) {
    mu_tau_j <- shift_tau_j + scale_tau_j * sqrt(2 / pi)
    sd_tau_j <- scale_tau_j * sqrt(1 - 2 / pi)
  }
  
  communality <- lambda_std %*% Phi_cor %*% t(lambda_std)     # J x J
  R_theta_j   <- communality + diag(1 - diag(communality))    # J x J
  
  gamma2_j <- (2 * reliab_j) / (L * (1 - reliab_j))
  # Expected value of the trial-level variance under each model
  if(identical(condition$Model, "exgaussian")) {
    Ev <- mu_sigma_j^2 + mu_tau_j^2
  } else {
    Ev <- mu_sigma_j^2
  }
  
  # Population standard deviation
  sd_theta_j <- sqrt(gamma2_j * Ev)
  
  # Output (flat, consistent list)
  list(
    lambda_std  = lambda_std,
    Phi_cor     = Phi_cor,
    communality = communality,
    R_theta_j   = R_theta_j,
    
    mu_alpha    = mu_alpha_j,
    sd_alpha    = sd_alpha_j,
    mu_theta    = mu_theta_j,
    sd_theta    = sd_theta_j,
    
    shift_sigma = shift_sigma_j,
    scale_sigma = scale_sigma_j,
    mu_sigma    = mu_sigma_j,
    sd_sigma    = sd_sigma_j,
    
    shift_tau   = shift_tau_j,
    scale_tau   = scale_tau_j,
    mu_tau      = mu_tau_j,
    sd_tau      = sd_tau_j,
    
    mu_delta    = mu_delta_j,
    sd_delta    = sd_delta_j,
    
    reliab      = reliab_j,
    gamma2      = gamma2_j
  )
}

# Non-aggregated summaries per condition
nonagg_summaries <- function(results, parameter_name, parameter,
                             Rho = FALSE, empty = FALSE, parnames = NULL,
                             method_of_moments = FALSE, Phi.lv = FALSE) {
  # Empty template
  empty_df <- data.frame(
    rep           = integer(0),
    parameter_set = character(0),
    parameter     = character(0),
    true_value    = numeric(0),
    est           = numeric(0),
    se            = numeric(0),
    bias          = numeric(0),
    abias         = numeric(0),
    rbias         = numeric(0),
    arbias_pct    = numeric(0),
    squared_error = numeric(0),
    stringsAsFactors = FALSE
  )
  if (empty) return(empty_df)
  
  # --- helpers ---
  get_cols <- function(pattern) grep(paste0(parameter_name, pattern), colnames(results))
  make_rho_names <- function(Phi) {
    idx <- which(lower.tri(Phi), arr.ind = TRUE)
    paste0("rho[", idx[,1], ",", idx[,2], "]")
  }
  get_parnames_standard <- function(cols) {
    base <- sub("_postmean$", "", cols)
    base <- sub(paste0("^", parameter_name, "_?"), "", base)
    base
  }
  to_num_mat <- function(x) {                # <- robust coercion
    out <- as.matrix(as.data.frame(x))
    storage.mode(out) <- "double"
    out
  }
  build_long_df <- function(est_mat, se_mat, true_vec, parnames_local) {
    est_mat <- to_num_mat(est_mat)
    if (!is.null(se_mat)) se_mat <- to_num_mat(se_mat)
    R <- nrow(est_mat); J <- ncol(est_mat)
    if (length(true_vec) != J) {
      stop("Length of 'parameter' must match number of selected estimate columns.")
    }
    if (is.null(parnames_local)) parnames_local <- paste0("par", seq_len(J))
    # replicate truths per replicate
    true_rep <- matrix(rep(as.numeric(true_vec), each = R), nrow = R, ncol = J)
    bias_mat <- est_mat - true_rep
    
    denom <- as.numeric(true_rep)
    denom[denom == 0] <- NA_real_  # avoid div-by-zero; yields NA for rbias/arbias_pct
    
    data.frame(
      rep           = rep(seq_len(R), times = J),
      parameter_set = parameter_name,
      parameter     = rep(parnames_local, each = R),
      true_value    = as.numeric(true_rep),
      est           = as.numeric(est_mat),
      se            = if (is.null(se_mat)) NA_real_ else as.numeric(se_mat),
      bias          = as.numeric(bias_mat),
      abias         = abs(as.numeric(bias_mat)),
      rbias         = as.numeric(bias_mat) / denom,
      arbias_pct    = 100 * abs(as.numeric(bias_mat)) / abs(denom),
      squared_error = as.numeric(bias_mat)^2,
      stringsAsFactors = FALSE
    )
  }
  # ---------------
  
  # Method of Moments
  if (method_of_moments) {
    est_mat <- results[, get_cols("_est"), drop = FALSE]
    se_mat  <- results[, get_cols("_se"),  drop = FALSE]
    if (ncol(est_mat) == 0) return(empty_df)
    
    if (is.matrix(parameter)) {
      true_vec <- parameter[lower.tri(parameter)]
      parnames_local <- if (is.null(parnames)) make_rho_names(parameter) else parnames
    } else {
      true_vec <- as.numeric(parameter)
      if (is.null(parnames)) {
        cols <- colnames(est_mat)
        parnames_local <- sub(paste0("^", parameter_name, "_?"), "", cols)
      } else parnames_local <- parnames
    }
    return(build_long_df(est_mat, se_mat, true_vec, parnames_local))
  }
  
  # Posterior means / sds
  postmean_mat <- results[, get_cols("_postmean"), drop = FALSE]
  postsd_mat   <- results[, get_cols("_postsd"),   drop = FALSE]
  if (ncol(postmean_mat) == 0) return(empty_df)
  
  if (Rho) {
    if (!is.matrix(parameter)) {
      stop("For Rho=TRUE, 'parameter' must be a correlation/covariance matrix.")
    }
    true_vec <- parameter[lower.tri(parameter)]
    parnames_local <- if (is.null(parnames)) make_rho_names(parameter) else parnames
    return(build_long_df(postmean_mat, postsd_mat, true_vec, parnames_local))
    
  } else if (Phi.lv) {
    true_vec <- if (length(parameter) == 1) rep(as.numeric(parameter), ncol(postmean_mat))
    else as.numeric(parameter)
    if (length(true_vec) != ncol(postmean_mat)) {
      stop("For Phi.lv=TRUE, length(parameter) must equal ncol(postmean_mat) or be a scalar.")
    }
    parnames_local <- if (!is.null(parnames)) parnames else sub("_postmean$", "", colnames(postmean_mat))
    return(build_long_df(postmean_mat, postsd_mat, true_vec, parnames_local))
    
  } else {
    true_vec <- as.numeric(parameter)
    if (length(true_vec) != ncol(postmean_mat)) {
      stop("Standard branch: length(parameter) must match number of selected '*_postmean' columns.")
    }
    parnames_local <- if (!is.null(parnames)) parnames else get_parnames_standard(colnames(postmean_mat))
    return(build_long_df(postmean_mat, postsd_mat, true_vec, parnames_local))
  }
}

# ─────────────────────────────────────────────────────────────────────────────