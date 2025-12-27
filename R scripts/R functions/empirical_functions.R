# ╔════════════════════════════════════════════════════════════════════════════╗
# ║                             SCRIPT OVERVIEW                                ║
# ╠════════════════════════════════════════════════════════════════════════════╣
# ║ Script Name   : empirical_functions.R                                      ║
# ║ Author        : Ricardo Rey-Sáez                                           ║
# ║ Role          : PhD Student in Psychology                                  ║
# ║ Institution   : Autonomous University of Madrid, Spain                     ║
# ║ Email         : ricardoreysaez95@gmail.es                                  ║
# ║ Date          : 08-04-2025                                                 ║
# ╚════════════════════════════════════════════════════════════════════════════╝

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 1: Description
# ─────────────────────────────────────────────────────────────────────────────
# Script with user-defined functions to use BHGFM in empirical settings
# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 2: Load Packages
# ─────────────────────────────────────────────────────────────────────────────

library(MASS)           # Multivariate gaussian distribution
library(dplyr)          # Stan array data
library(posterior)      # cmdstanr summaries
library(infinitefactor) # MatchAlign algorithm functions
library(lavaan)         # Latent variable analysis: Spearman adjustment

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 3: Data manipulation and starting values
# ─────────────────────────────────────────────────────────────────────────────

# Ensure packages inside each function; if not, install it.
ensure_packages <- function(packages) {
  for (pkg in packages) {
    # Check if package is installed
    if (!requireNamespace(pkg, quietly = TRUE)) {
      # If not, install it
      message(sprintf("Installing package: %s", pkg))
      install.packages(pkg, dependencies = TRUE)
    }
    # Load packages
    library(pkg, character.only = TRUE)
  }
}

# Transform long data frame to Stan input data
Stan.list.data <- function(data, 
                           # Column names in original data frame
                           subject_var = "subject", 
                           task_var = "task", 
                           condition_var = "condition", 
                           congruent_val = -0.5,
                           incongruent_val = 0.5,
                           rt_var = "RT") {
  # Ensure dplyr
  ensure_packages(packages = "dplyr")
  
  # Ensure condition values
  cond_vals <- unique(data[[condition_var]])
  if (!all(c(congruent_val, incongruent_val) %in% cond_vals)) {
    stop(
      "Column '", condition_var, 
      "' must containt both specified values: ",
      paste(c(congruent_val, incongruent_val), collapse = ", "),
      ". Unique values detected: ",
      paste(head(cond_vals, 10), collapse = ", "),
      if (length(cond_vals) > 10) " ..." else "",
      call. = FALSE
    )
  }
  
  ## 1) Rename columns, remove NAs
  data_renamed <- data %>%
    rename(
      subject   = all_of(subject_var),
      condition = all_of(condition_var),
      task      = all_of(task_var),
      RT        = all_of(rt_var)
    ) %>%
    na.omit()
  
  # Rename condition values 
  data_renamed <- data_renamed %>%
    dplyr::mutate(
      condition = ifelse(condition == congruent_val, -0.5, 0.5)
    )
  
  ## 2) Filter out subjects who do *not* have observations in *all* condition x task combos
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
  
  # Save task labels and factor correspondence
  task_levels  <- levels(factor(data_valid$task))
  task_link  <- data.frame(
    task_num   = seq_along(task_levels),
    task_label = task_levels,
    stringsAsFactors = FALSE
  )
  
  ## 3) Re-factor (numeric) subject, condition, and task
  data_valid <- data_valid %>%
    mutate(
      subject   = as.numeric(factor(subject)),
      condition = as.numeric(factor(condition)),
      task      = as.numeric(factor(task, levels = task_levels))
    )
  
  # Number of subjects, conditions, tasks
  n_subj <- n_distinct(data_valid$subject)
  n_cond <- n_distinct(data_valid$condition)
  n_task <- n_distinct(data_valid$task)
  
  ## 4) Summaries: T_subj, RT_min, RT_max
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
  
  ## 5) Build the 4D array of RT in a *single pass*
  #  (subject, condition, task, trial_id)
  #  Create a trial_id from 1..n within each group
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
  
  ## 6) Return the final list for Stan
  stan_dat <- list(
    I      = n_subj,
    K      = n_cond,
    J      = n_task,
    L_max  = T_max,
    T_subj = T_subj,
    RT     = RT,
    RT_min = RT_min,
    RT_max = RT_max,
    task_link = task_link
  )
  
  return(stan_dat)
}

# Transform data from Stan input to long data frame
sdata_to_longdf <- function(stan_dat, 
                            # Desired column names in final data frame
                            subject_col = "subject", 
                            condition_col = "condition", 
                            task_col = "task", 
                            trial_col = "trial", 
                            rt_col = "rt") {
  
  # Extraemos los componentes del stan_dat
  RT_array <- stan_dat$RT
  n_subj   <- stan_dat$I
  n_cond   <- stan_dat$K
  n_task   <- stan_dat$J
  T_subj   <- stan_dat$T_subj
  
  # Predefinir la lista para almacenar data frames
  temp_list <- vector("list", n_subj * n_cond * n_task)
  counter <- 1
  
  # Recorremos cada combinación de sujeto, condición y tarea para reconstruir el data frame
  for (i in 1:n_subj) {
    for (c in 1:n_cond) {
      for (t in 1:n_task) {
        n_trials <- T_subj[i, c, t]  # Número de pruebas para esta combinación
        if (n_trials > 0) {
          # Extraemos los tiempos de reacción para estas pruebas
          rt_values <- RT_array[i, c, t, 1:n_trials]
          
          # Creamos un data frame temporal con nombres provisionales
          temp_df <- data.frame(
            subject   = i,
            condition = c,
            task      = t,
            trial     = 1:n_trials,
            rt        = rt_values
          )
          
          # Renombramos las columnas usando los argumentos pasados a la función
          names(temp_df) <- c(subject_col, condition_col, task_col, trial_col, rt_col)
          
          # Almacenamos el data frame temporal en la lista
          temp_list[[counter]] <- temp_df
          counter <- counter + 1
        }
      }
    }
  }
  
  # Combina todos los data frames almacenados en la lista en un solo data frame
  reconstructed_data <- do.call(rbind, temp_list[1:(counter - 1)])
  
  return(reconstructed_data)
}

# Informative starting values per model
make_inits <- function(data, 
                       # Column names in original data frame
                       subject_var = "subject", 
                       task_var = "task", 
                       condition_var = "condition", 
                       rt_var = "RT",
                       # Identifier for congruent trials in df (i.e, 0, "congruent", etc...)
                       congruent_id = 1,    
                       # Identifier for congruent trials in df (i.e, 1, "incongruent", etc...)
                       incongruent_id = 2,
                       # Number of latent factors
                       M, 
                       # Desired model: gaussian, exgaussian or shifted.lognormal
                       model = NULL){
  
  # Ensure dplyr
  ensure_packages(c("dplyr", "truncnorm"))
  
  # Rename columns, remove NAs
  df <- data %>%
    rename(
      subject   = all_of(subject_var),
      condition = all_of(condition_var),
      task      = all_of(task_var),
      RT        = all_of(rt_var)
    ) %>%
    na.omit()
  
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
    pre_beta <- tapply(df$RT, list(df$subject, df$condition, df$task), mean)
    beta <- pre_beta[,2,] - pre_beta[,1,]
    sigma <- tapply(df$RT, list(df$subject, df$task), sd)
    # Exploratory factor analysis: varimax
    h2 <- psych::fa(beta, M, fm = "minrank", rotate = "varimax")$communality
    # Exploratory factor analysis: oblimin
    oblFA <- psych::fa(beta, M, fm = "minrank", rotate = "oblimin")
    if(M > 1){ 
      L_Phi <- t(chol(oblFA$Phi))
    } else {
      L_Phi <- diag(1)
    }
    
    # Mean and sd of level-1 residual variance
    mu_sigma = colMeans(sigma)
    sd_sigma = apply(sigma, 2, sd)
    
    # Compute shifts and scale values
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
      mu_beta = colMeans(beta),
      shift_sigma = shift_sigma,
      # Population standard deviations
      sd_alpha = apply(alpha, 2, sd),
      sd_beta = apply(beta, 2, sd),
      scale_sigma = scale_sigma,
      # Cholesky factor decomposition
      L_R_alpha = t(chol(cor(alpha))),
      L_Phi     = L_Phi,
      # Individual-level latent scores
      alpha_tilde = t(apply(alpha, 2, function(x) (x - mean(x)) / sd(x))),
      beta_tilde = t(apply(beta, 2, function(x) (x - mean(x)) / sd(x))),
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
    pre_beta <- tapply(df$RT, list(df$subject, df$condition, df$task), function(x) mexgauss(x)[1])
    beta <- pre_beta[,2,] - pre_beta[,1,]
    sigma <- tapply(df$RT, list(df$subject, df$task), function(x) mexgauss(x)[2])
    tau   <- tapply(df$RT, list(df$subject, df$task), function(x) mexgauss(x)[3])
    
    # Exploratory factor analysis: varimax
    h2 <- psych::fa(beta, M, fm = "minrank", rotate = "varimax")$communality
    # Exploratory factor analysis: oblimin
    oblFA <- psych::fa(beta, M, fm = "minrank", rotate = "oblimin")
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
    
    # Compute shifts and scale values
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
      mu_beta = colMeans(beta),
      shift_sigma = shift_sigma,
      shift_rate  = shift_rate,
      shift_tau  = shift_tau,
      # Population standard deviations
      sd_alpha = apply(alpha, 2, sd),
      sd_beta = apply(beta, 2, sd),
      scale_sigma = scale_sigma,
      scale_rate  = scale_rate,
      scale_tau   = scale_tau,
      # Cholesky factor decomposition
      L_R_alpha = t(chol(cor(alpha))),
      L_Phi     = L_Phi,
      # Communalities, unit vector and factor loadings
      h2 = rep(0.5, J),
      L_raw = matrix(c(oblFA$loadings), ncol = M),
      Z = z_normed,
      # Individual-level latent scores
      alpha_tilde = t(apply(alpha, 2, function(x) (x - mean(x)) / sd(x))),
      beta_tilde = t(apply(beta, 2, function(x) (x - mean(x)) / sd(x))),
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
    alpha <- beta <- sigma <- matrix(NA, ncol = length(unique(df$task)), 
                                      nrow = length(unique(df$subject)))
    for(i in 1:length(unique(df$subject))){
      for(j in 1:length(unique(df$task))){
        alpha[i,j] <- mean(log(df$RT[which(df$subject==i & df$task==j)] - delta[i,j]))
        beta[i,j] <- mean(log(df$RT[which(df$subject==i & df$task==j & df$condition == incongruent_id)] - delta[i,j])) - 
                      mean(log(df$RT[which(df$subject==i & df$task==j & df$condition == congruent_id)] - delta[i,j]))
        sigma[i,j] <- sd(log(df$RT[which(df$subject==i & df$task==j)] - delta[i,j]))
      }
    }
    
    # Exploratory factor analysis: varimax
    h2 <- psych::fa(beta, M, fm = "minrank", rotate = "varimax")$communality
    # Exploratory factor analysis: oblimin
    oblFA <- psych::fa(beta, M, fm = "minrank", rotate = "oblimin")
    if(M > 1){ 
      L_Phi <- t(chol(oblFA$Phi))
    } else {
      L_Phi <- diag(1)
    }
    
    # Mean and sd of level-1 residual variance
    mu_sigma = colMeans(sigma)
    sd_sigma = apply(sigma, 2, sd)
    
    # Compute shifts and scale values
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
      mu_beta = colMeans(beta),
      shift_sigma = shift_sigma,
      mu_delta = colMeans(delta),
      # Population standard deviations
      sd_alpha = apply(alpha, 2, sd),
      sd_beta = apply(beta, 2, sd),
      scale_sigma = scale_sigma,
      sd_delta = apply(delta, 2, sd),
      # Cholesky factor decomposition
      L_R_alpha = t(chol(cor(alpha))),
      L_Phi     = L_Phi,
      # Individual-level latent scores
      alpha_tilde = t(apply(alpha, 2, FUN = function(x) (x - mean(x)) / sd(x))),
      beta_tilde = t(apply(beta, 2, FUN = function(x) (x - mean(x)) / sd(x))),
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

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 4: Exploratory Factor Analysis rotation and label switching 
# ─────────────────────────────────────────────────────────────────────────────

# MatchAlign Algorithm: orthogonal factor rotation across draws
# Paper: https://arxiv.org/abs/2107.13783
# Based on: https://github.com/poworoznek/infinitefactor/blob/master/R/jointRot.R
MatchAlign <- function(fit, lambda_name, only_summary = TRUE) {
  # Unrotated posterior distribution
  unrotated_draws <- fit$draws(lambda_name, format = "matrix")
  # Number of items and latent common factors
  J <- max(as.integer( sub(",.*", "", sub(".*\\[", "", colnames(unrotated_draws)))))
  M <- max(as.integer( sub(".*,", "", sub("\\]",   "", colnames(unrotated_draws)))))
  # Apply varimax rotation criteria to all unrotated posterior draws
  vari_all <- apply(unrotated_draws, 1, function(x) varimax(matrix(x, nrow = J, ncol = M)))
  # Varimax rotated posterior distribution
  varimax_loadings <- lapply(vari_all, `[[`, 1)
  # Apply singular value decomposition and select largest singular value
  singular_values <- sapply(varimax_loadings, norm, "2")
  # Target varimax matrix: median singular value
  target <- varimax_loadings[order(singular_values)][[round(length(singular_values)/2)]]
  # Rotate posterior draws using selected target matrix
  rotated_draws <- t(sapply(varimax_loadings, function(x) infinitefactor::msf(x, target)))
  # Ensure positive reflection
  sign_factors <- rep(sign(colMeans(matrix(colMeans(rotated_draws), nrow = J, ncol = M))), each = J)
  rotated_draws <- t(apply(rotated_draws, 1, function(x) x * sign_factors))
  # Posterior distribution with posterior basic format
  colnames(rotated_draws) <- colnames(unrotated_draws)
  rotated_draws <- posterior::as_draws_matrix(rotated_draws)
  # Posterior summaries
  post_sum <- posterior::summarise_draws(rotated_draws)
  # Return results
  if(only_summary) { 
    # Posterior distribution summary
    return(post_sum) 
  } else {
    # Rotated posterior distribution draws and posterior summary
    return(list(post_sum = post_sum, rotated_draws = rotated_draws)) }
}

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 5: Posterior Predictive Model Checks
# ─────────────────────────────────────────────────────────────────────────────

# Simulate new level-1 observations
PPC_simdata <- function(fit,        # cmdstanr fitted model
                        sdata,      # Input data in the fitted model
                        n_draws,    # Number of desired level-1 simulated responses
                        model       # gaussian, exgaussian and shifted.lognormal
                        ) {
  # 1) Extract necessary values from the Stan data
  I      <- sdata$I
  J      <- sdata$J
  K      <- sdata$K
  L_max  <- sdata$L_max
  T_subj <- sdata$T_subj
  
  # 2) Extract the posterior draws in matrix form once
  post_samp <- fit$draws(format = "draws_matrix")
  cnames    <- colnames(post_samp)
  
  # 3) Identify columns for each parameter just once
  alpha_idx <- grep("^alpha\\[",  cnames)
  beta_idx  <- grep("^beta\\[",   cnames)
  Sigma_idx <- grep("^sigma\\[\\d+,\\d+\\]$", cnames)
  Delta_idx <- grep("^Delta\\[\\d+,\\d+\\]$", cnames)
  Rate_idx  <- grep("^rate\\[\\d+,\\d+\\]$",  cnames)
  
  # 4) Select random draws
  s_draws <- sample(seq_len(nrow(post_samp)), size = n_draws)
  
  # 5) Subset the posterior samples to those draws only
  post_samp_sub <- post_samp[s_draws, ]
  
  # 6) Reshape each parameter into a 3D array of dimensions (n_draws x I x J)
  #    so we can slice them easily in the nested loops
  alpha_array <- array(post_samp_sub[, alpha_idx], dim = c(n_draws, I, J))
  beta_array  <- array(post_samp_sub[, beta_idx], dim = c(n_draws, I, J))
  
  Sigma_array <- if (model %in% c("gaussian", "exgaussian", "shifted.lognormal")) {
    array(post_samp_sub[, Sigma_idx], dim = c(n_draws, I, J))
  } else {
    NULL
  }
  Delta_array <- if (model == "shifted.lognormal") {
    array(post_samp_sub[, Delta_idx], dim = c(n_draws, I, J))
  } else {
    NULL
  }
  Rate_array <- if (model == "exgaussian") {
    array(post_samp_sub[, Rate_idx], dim = c(n_draws, I, J))
  } else {
    NULL
  }
  
  # 7) We'll build up a list of data frames (one per (i,k,j) combo)
  output_list <- vector("list", I*K*J)
  idx_list    <- 1
  
  # 8) Nested loops: subject i, condition k, task j
  for (i in seq_len(I)) {
    for (k in seq_len(K)) {
      for (j in seq_len(J)) {
        n_trials <- T_subj[i, k, j]
        
        # Only build a data frame if n_trials > 0
        if (n_trials > 0) {
          
          # Observed data
          temp_df <- data.frame(
            subject   = i,
            condition = k,
            task      = j,
            trial     = seq_len(n_trials),
            RT_obs    = sdata$RT[i, k, j, seq_len(n_trials)]
          )
          
          # For each posterior draw, simulate RTs directly into a new column
          for (n in seq_len(n_draws)) {
            alpha_ij <- alpha_array[n, i, j]
            beta_ij  <- beta_array[n, i, j]
            sigma_ij <- if (!is.null(Sigma_array)) Sigma_array[n, i, j] else NA_real_
            delta_ij <- if (!is.null(Delta_array)) Delta_array[n, i, j] else NA_real_
            rate_ij  <- if (!is.null(Rate_array))  Rate_array[n, i, j]  else NA_real_
            
            # Mean structure
            mu_pred <- if (k == 1) {
              alpha_ij - 0.5 * beta_ij
            } else {
              alpha_ij + 0.5 * beta_ij
            }
            
            # Generate random data according to the model
            if (model == "gaussian") {
              sim_vals <- rnorm(n_trials, mean = mu_pred, sd = sigma_ij)
            } else if (model == "shifted.lognormal") {
              sim_vals <- delta_ij + rlnorm(n_trials, meanlog = mu_pred, sdlog = sigma_ij)
            } else if (model == "exgaussian") {
              sim_vals <- rnorm(n_trials, mean = mu_pred, sd = sigma_ij) +
                rexp(n_trials, rate = rate_ij)
            } else {
              # Fallback or error-handling for unknown model
              sim_vals <- rep(NA_real_, n_trials)
            }
            
            # Store in the data frame
            temp_df[[paste0("RT_draw_", n)]] <- sim_vals
          }
          
          # Store this in our output list
          output_list[[idx_list]] <- temp_df
          idx_list <- idx_list + 1
        }
      }
    }
  }
  
  # 9) Combine into one large data frame and return
  postpred_df <- do.call(rbind, output_list[seq_len(idx_list - 1)])
  return(postpred_df)
}

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 6: PSIS-LOO, Mixture-IS-LOO and moment_match estimation outside Stan
# ─────────────────────────────────────────────────────────────────────────────

# Basic PSIS-LOO-CV elpd estimation
loo_BGHFM <- function(fit, 
                      model, 
                      data, 
                      subject_var = "subject", 
                      task_var = "task", 
                      condition_var = "condition", 
                      rt_var = "RT",
                      n.cores = 2, 
                      estimate_r_eff = FALSE) {
  
  # Ensure packages
  ensure_packages(c("loo", "posterior", "dplyr"))
  
  # Prepare data
  data <- data %>%
    rename(
      subject   = all_of(subject_var),
      condition = all_of(condition_var),
      task      = all_of(task_var),
      RT        = all_of(rt_var)
    ) %>%
    na.omit()
  
  # Gaussian model conditional log-likelihood
  llfun_gaussian <- function(data_i, draws, log = TRUE){
    # Extract relevant data for the i-th observation
    RT_i <- data_i$RT
    condition <- data_i$condition
    task <- data_i$task
    subject <- data_i$subject
    
    # Select
    colpars <- paste0(c("alpha", "beta", "sigma"),"[",subject, ",", task, "]")
    
    # Compute the linear predictor
    mu_ijk <- draws[,colpars[1]] + (condition - 1.5) * draws[,colpars[2]]
    
    # Calculate the log-likelihood for each posterior draw
    log_lik_values <- dnorm(RT_i, mean = mu_ijk, sd = draws[,colpars[3]], log = log)
    
    return(log_lik_values)
  }
  
  # Exgaussian model conditional log-likelihood
  llfun_exgaussian <- function(data_i, draws, log = TRUE){
    
    # Ex-gaussian log-probability density function
    dexgauss <- function(x, mu, sigma, lambda, log = TRUE) {
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
    colpars <- paste0(c("alpha", "beta", "sigma", "rate"),"[",subject,",", task, "]")
    
    # Compute the linear predictor
    mu_ijk <- draws[,colpars[1]] + (condition - 1.5) * draws[,colpars[2]]
    
    # Calculate the log-likelihood for each posterior draw
    log_lik_values <- dexgauss(x = RT_i, mu = mu_ijk, sigma = draws[,colpars[3]], 
                               lambda = draws[,colpars[4]], log = log)
    
    return(log_lik_values)
  }
  
  # Shifted-lognormal model conditional log-likelihood
  llfun_shiftlognormal <- function(data_i, draws, log = TRUE){
    
    # Ex-gaussian log-probability density function
    dshifted_lognorm <- function(x, mu, sigma, delta, log = TRUE) {
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
    colpars <- paste0(c("alpha", "beta", "sigma", "Delta"),"[",subject, ",", task, "]")
    
    # Compute the linear predictor
    mu_ijk <- draws[,colpars[1]] + (condition - 1.5) * draws[,colpars[2]]
    
    # Calculate the log-likelihood for each posterior draw
    log_lik_values <- dshifted_lognorm(x = RT_i, mu = mu_ijk, sigma = draws[,colpars[3]], 
                                       delta = draws[,colpars[4]], log = log)
    
    return(log_lik_values)
  }
  
  # All conditional log-likelihood functions
  ll_functions <- list(
    "gaussian"          = llfun_gaussian,
    "exgaussian"        = llfun_exgaussian,
    "shifted.lognormal" = llfun_shiftlognormal
  )
  
  # Parameters by model
  model_parameters <- list(
    "gaussian"          = c("alpha", "beta", "sigma"),
    "exgaussian"        = c("alpha", "beta", "sigma", "rate"),
    "shifted.lognormal" = c("alpha", "beta", "sigma", "Delta")
  )
  
  # Posterior draws per model
  posterior_draws <- posterior::subset_draws(fit$draws(model_parameters[[model]], format = "matrix"))
  
  # Compute relative efficiency
  if(estimate_r_eff){
    r_eff <- relative_eff(ll_functions[[model]], 
                          log = FALSE, 
                          chain_id = rep(1:fit$num_chains(), each = nrow(posterior_draws) / fit$num_chains()), 
                          data = data, 
                          draws = posterior_draws, 
                          cores = n.cores)
  } else {
    r_eff <- 1
  }
  
  # Compute PSIS-LOO (and save time)
  loo_model <- loo(ll_functions[[model]], draws = posterior_draws, 
                   r_eff = r_eff, data = data, cores = n.cores)
  
  # Return LOO and time
  return(loo_model)
}

# Moment Matching for high pareto-k values
loo_mm_BGHFM <- function(fit,
                         model,
                         loo_object,
                         data,
                         subject_var = "subject", 
                         task_var = "task",
                         condition_var = "condition", 
                         rt_var = "RT",
                         n.cores = 1) {
  
  # Ensure packages
  ensure_packages(c("loo", "posterior", "dplyr"))
  
  # Prepare data
  data <- data %>%
    rename(
      subject   = all_of(subject_var),
      condition = all_of(condition_var),
      task      = all_of(task_var),
      RT        = all_of(rt_var)
    ) %>%
    na.omit()
  
  # Compile additiona model methods
  fit$init_model_methods()
  
  # ────────────────────────────── #
  #    Log-likelihood functions    #
  # ────────────────────────────── #
  
  llfun_gaussian <- function(x, i, data, ...){
    # Extract relevant data for the i-th observation
    RT_i <- data$RT[i]
    condition <- data$condition[i]
    task <- data$task[i]
    subject <- data$subject[i]
    
    # Parameter names for raw i
    colpars <- paste0(c("alpha", "beta", "sigma"),"[",subject, ",", task, "]")
    
    # Posterior draws for this parameters
    post_draws <- x$draws(variables = colpars, format = "draws_array")
    
    # Compute the linear predictor
    mu_ijk <- post_draws[,,colpars[1]] + (condition - 1.5) * post_draws[,,colpars[2]]
    
    # Calculate the log-likelihood for each posterior draw
    log_lik_values <- dnorm(RT_i, mean = mu_ijk, sd = post_draws[,,colpars[3]], log = TRUE)
    
    return(log_lik_values)
  }
  
  llfun_exgaussian <- function(x, i, data, ...){
    # Ex-gaussian log-probability density function
    dexgauss <- function(x, mu, sigma, lambda, log = TRUE) {
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
    RT_i <- data$RT[i]
    condition <- data$condition[i]
    task <- data$task[i]
    subject <- data$subject[i]
    
    # Parameter names for raw i
    colpars <- paste0(c("alpha", "beta", "sigma", "rate"),"[",subject, ",", task, "]")
    
    # Posterior draws for this parameters
    post_draws <- x$draws(variables = colpars, format = "draws_array")
    
    # Compute the linear predictor
    mu_ijk <- post_draws[,,colpars[1]] + (condition - 1.5) * post_draws[,,colpars[2]]
    
    # Calculate the log-likelihood for each posterior draw
    log_lik_values <- dexgauss(RT_i, mu = mu_ijk, sigma = post_draws[,,colpars[3]], 
                               lambda = post_draws[,,colpars[4]], log = TRUE)
    
    return(log_lik_values)
  }
  
  llfun_shlognormal <- function(x, i, data, ...){
    
    # Ex-gaussian log-probability density function
    dshifted_lognorm <- function(x, mu, sigma, delta, log = TRUE) {
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
    RT_i <- data$RT[i]
    condition <- data$condition[i]
    task <- data$task[i]
    subject <- data$subject[i]
    
    # Parameter names for raw i
    colpars <- paste0(c("alpha", "beta", "sigma", "Delta"),"[",subject, ",", task, "]")
    
    # Posterior draws for this parameters
    post_draws <- x$draws(variables = colpars, format = "draws_array")
    
    # Compute the linear predictor
    mu_ijk <- post_draws[,,colpars[1]] + (condition - 1.5) * post_draws[,,colpars[2]]
    
    # Calculate the log-likelihood for each posterior draw
    log_lik_values <- dshifted_lognorm(x = RT_i, mu = mu_ijk, sigma = post_draws[,,colpars[3]], 
                                       delta = post_draws[,,colpars[4]], log = TRUE)
    
    return(log_lik_values)
  }
  
  # ──────────────────────────────────────────────────────── #
  #    Log-likelihood functions on unconstrain parameters    #
  # ──────────────────────────────────────────────────────── #
  
  llfun_gaussian_i_upars <- function(x, upars, i, data, ...) {
    # Extract relevant data for the i-th observation
    RT_i <- data$RT[i]
    condition <- data$condition[i]
    task <- data$task[i]
    subject <- data$subject[i]
    
    # Number of draws
    S <- nrow(upars)
    
    # Empty vector of log-likelihoods
    ll_values <- numeric(S)
    
    # For each draw s
    for(s in 1:S) {
      # Constraint parameters for draw s
      cons_pars <- x$constrain_variables(upars[s,])
      # Parameters for subject i in task j
      alpha <- cons_pars[["alpha"]][subject, task]
      beta <- cons_pars[["beta"]][subject, task]
      sigma <- cons_pars[["sigma"]][subject, task]
      # Linear prediction
      mu_ijk <- alpha + beta * (condition - 1.5)
      # Log-likelihood in draw s
      ll_values[s] <- dnorm(x = RT_i, mean = mu_ijk, sd = sigma, log = TRUE)
    }
    # Return ll values
    return(ll_values)
  }
  
  llfun_exgaussian_i_upars <- function(x, upars, i, data, ...) {
    
    # Ex-gaussian log-probability density function
    dexgauss <- function(x, mu, sigma, lambda, log = TRUE) {
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
    RT_i <- data$RT[i]
    condition <- data$condition[i]
    task <- data$task[i]
    subject <- data$subject[i]
    
    # Number of draws
    S <- nrow(upars)
    
    # Empty vector of log-likelihoods
    ll_values <- numeric(S)
    
    # For each draw s
    for(s in 1:S) {
      # Constraint parameters for draw s
      cons_pars <- x$constrain_variables(upars[s,])
      # Parameters for subject i in task j
      alpha <- cons_pars[["alpha"]][subject, task]
      beta <- cons_pars[["beta"]][subject, task]
      sigma <- cons_pars[["sigma"]][subject, task]
      rate  <- cons_pars[["rate"]][subject, task]
      # Linear prediction
      mu_ijk <- alpha + beta * (condition - 1.5)
      # Log-likelihood in draw s
      ll_values[s] <- dexgauss(x = RT_i, mu = mu_ijk, sigma = sigma,
                               lambda = rate, log = TRUE)
    }
    # Return ll values
    return(ll_values)
  }
  
  llfun_shlognormal_i_upars <- function(x, upars, i, data, ...) {
    
    # Ex-gaussian log-probability density function
    dshifted_lognorm <- function(x, mu, sigma, delta, log = TRUE) {
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
    RT_i <- data$RT[i]
    condition <- data$condition[i]
    task <- data$task[i]
    subject <- data$subject[i]
    
    # Number of draws
    S <- nrow(upars)
    
    # Empty vector of log-likelihoods
    ll_values <- numeric(S)
    
    # For each draw s
    for(s in 1:S) {
      # Constraint parameters for draw s
      cons_pars <- x$constrain_variables(upars[s,])
      # Parameters for subject i in task j
      alpha <- cons_pars[["alpha"]][subject, task]
      beta <- cons_pars[["beta"]][subject, task]
      sigma <- cons_pars[["sigma"]][subject, task]
      Delta <- cons_pars[["Delta"]][subject, task]
      # Linear prediction
      mu_ijk <- alpha + beta * (condition - 1.5)
      # Log-likelihood in draw s
      ll_values[s] <- dshifted_lognorm(x = RT_i, mu = mu_ijk, sigma = sigma,
                                       delta = Delta, log = TRUE)
    }
    # Return ll values
    return(ll_values)
  }
  
  # ──────────────────────────────────────────────────── #
  #    Select log-likelihood functions for each model    #
  # ──────────────────────────────────────────────────── #
  
  # All conditional log-likelihood functions
  ll_functions <- list(
    "gaussian"          = llfun_gaussian,
    "exgaussian"        = llfun_exgaussian,
    "shifted.lognormal" = llfun_shlognormal
  )
  
  # All unconstrain log-likelihood functions
  ll_uncons_functions <- list(
    "gaussian"          = llfun_gaussian_i_upars,
    "exgaussian"        = llfun_exgaussian_i_upars,
    "shifted.lognormal" = llfun_shlognormal_i_upars
  )
  
  # Parameters by model
  model_parameters <- list(
    "gaussian"          = c("alpha", "beta", "sigma"),
    "exgaussian"        = c("alpha", "beta", "sigma", "rate"),
    "shifted.lognormal" = c("alpha", "beta", "sigma", "Delta")
  )
  
  # ──────────────────────────────────────────────────── #
  #    Pre-load posterior draws to improve efficiency    #
  # ──────────────────────────────────────────────────── #
  
  # Posterior draws per model
  posterior_draws <- fit$draws(model_parameters[[model]], format = "matrix")
  
  # Unconstrain parameters per model
  uncon_draws <- fit$unconstrain_draws()
  
  # ─────────────────────────────── #
  #    Apply moment-match method    #
  # ─────────────────────────────── #
  
  loo_mm_model <- loo_moment_match(
    x                 = fit,
    loo               = loo_object,
    post_draws        = function(x, ...) { x$draws(format = "matrix") },
    log_lik_i         = ll_functions[[model]],
    unconstrain_pars  = function(x, pars, ...) { x$unconstrain_draws(format = "matrix") },
    log_prob_upars    = function(x, upars, ...) { apply(upars, 1, x$log_prob) },
    log_lik_i_upars   = ll_uncons_functions[[model]], 
    cores             = n.cores, 
    data              = data)
  
  # Return LOO and time
  return(loo_mm_model)
}

# Mixture-IS elpd estimation
# note: it needs a model fitted with mixture-IS approach!
loo_mixis_BGHFM <- function(
    fit,
    model = c("gaussian", "exgaussian", "shifted.lognormal"),
    data,
    subject_var      = "subject",
    task_var         = "task",
    condition_var    = "condition",
    rt_var           = "RT",
    condition_center = 1.5,
    progress_bar     = TRUE,
    update_every     = 1L
) {
  if (!requireNamespace("posterior",   quietly = TRUE)) stop("Install 'posterior'.")
  if (!requireNamespace("matrixStats", quietly = TRUE)) stop("Install 'matrixStats'.")
  if (!requireNamespace("dplyr",       quietly = TRUE)) stop("Install 'dplyr'.")
  if (isTRUE(progress_bar) && !requireNamespace("progressr", quietly = TRUE)) {
    stop("Install 'progressr' or set progress_bar = FALSE.")
  }
  
  logaddexp_scalar <- function(a, b) {
    # Casos con NA/NaN
    if (is.na(a) && is.na(b)) return(NA_real_)
    if (is.na(a)) return(b)
    if (is.na(b)) return(a)
    
    # Casos con Inf
    if (a ==  Inf || b ==  Inf) return( Inf)
    if (a == -Inf && b == -Inf) return(-Inf)
    if (a == -Inf && is.finite(b)) return(b)
    if (b == -Inf && is.finite(a)) return(a)
    
    # Ambos finitos
    m <- max(a, b)
    m + log1p(exp(-abs(a - b)))
  }
  
  logaddexp_vec <- function(a, b) {
    stopifnot(length(a) == length(b))
    res <- rep(NA_real_, length(a))
    
    # where a is NA and b not NA -> res=b
    idx <- is.na(a) & !is.na(b)
    if (any(idx)) res[idx] <- b[idx]
    
    # where b is NA and a not NA -> res=a
    idx <- is.na(b) & !is.na(a)
    if (any(idx)) res[idx] <- a[idx]
    
    # both -Inf
    idx <- (a == -Inf) & (b == -Inf)
    idx[is.na(idx)] <- FALSE
    if (any(idx)) res[idx] <- -Inf
    
    # any +Inf
    idx <- (a ==  Inf) | (b ==  Inf)
    idx[is.na(idx)] <- FALSE
    if (any(idx)) res[idx] <- Inf
    
    # one -Inf and the other finite
    idx <- ((a == -Inf) & is.finite(b)) | ((b == -Inf) & is.finite(a))
    idx[is.na(idx)] <- FALSE
    if (any(idx)) res[idx] <- ifelse(a[idx] == -Inf, b[idx], a[idx])
    
    # both finite
    idx <- is.finite(a) & is.finite(b)
    idx[is.na(idx)] <- FALSE
    if (any(idx)) {
      aa <- a[idx]; bb <- b[idx]
      m  <- pmax(aa, bb)
      res[idx] <- m + log1p(exp(-abs(aa - bb)))
    }
    
    res
  }
  
  .erfc <- function(z) 2 * pnorm(-z * sqrt(2))
  
  dexgauss <- function(x, mu, sigma, lambda, log = TRUE) {
    if (any(sigma <= 0))  stop("sigma must be > 0")
    if (any(lambda <= 0)) stop("lambda (rate) must be > 0")
    lp <- log(lambda / 2) +
          0.5 * lambda * (2 * mu + lambda * sigma^2 - 2 * x) +
          log(.erfc((mu + lambda * sigma^2 - x) / (sqrt(2) * sigma)))
    if (log) lp else exp(lp)
  }
  dshifted_lognorm <- function(x, mu, sigma, delta, log = TRUE) {
    if (any(sigma <= 0)) stop("sigma must be > 0")
    if (any(x < delta))  stop("All x must be >= delta")
    lp <- -log(x - delta) - log(sigma * sqrt(2 * pi)) -
          ((log(x - delta) - mu)^2 / (2 * sigma^2))
    if (log) lp else exp(lp)
  }
  
  model <- match.arg(model)
  model_specs <- list(
    gaussian = list(
      params = c("alpha","beta","sigma"),
      ll_fun = function(x, pars, mu) dnorm(x, mean = mu, sd = pars$sigma, log = TRUE)
    ),
    exgaussian = list(
      params = c("alpha","beta","sigma","rate"),
      ll_fun = function(x, pars, mu) dexgauss(x, mu = mu, sigma = pars$sigma, lambda = pars$rate, log = TRUE)
    ),
    shifted.lognormal = list(
      params = c("alpha","beta","sigma","Delta"),
      ll_fun = function(x, pars, mu) dshifted_lognorm(x, mu = mu, sigma = pars$sigma, delta = pars$Delta, log = TRUE)
    )
  )
  spec <- model_specs[[model]]
  
  data <- dplyr::rename(
    data,
    subject   = dplyr::all_of(subject_var),
    condition = dplyr::all_of(condition_var),
    task      = dplyr::all_of(task_var),
    RT        = dplyr::all_of(rt_var)
  )
  data <- stats::na.omit(data)
  N      <- nrow(data)
  cond_c <- data$condition - condition_center
  
  posterior_draws <- fit$draws(variables = spec$params, format = "matrix")
  S  <- nrow(posterior_draws)
  cn <- colnames(posterior_draws)
  name_to_idx <- setNames(seq_along(cn), cn)
  
  ij_str  <- paste0("[", data$subject, ",", data$task, "]")
  get_idx <- function(p) unname(as.integer(name_to_idx[paste0(p, ij_str)]))
  
  alpha_idx <- get_idx("alpha")
  beta_idx  <- get_idx("beta")
  sigma_idx <- get_idx("sigma")
  extra_idx <- list()
  if ("rate"  %in% spec$params)  extra_idx$rate  <- get_idx("rate")
  if ("Delta" %in% spec$params) extra_idx$Delta <- get_idx("Delta")
  
  log_sum_invZtilde <- -Inf
  log_sum_pi_i      <- rep(-Inf, N)
  
  # --- progressr ---
  runner <- function() {
    p <- progressr::progressor(steps = S)
    for (s in seq_len(S)) {
      alpha_s <- posterior_draws[s, alpha_idx, drop = TRUE]
      beta_s  <- posterior_draws[s, beta_idx,  drop = TRUE]
      sigma_s <- posterior_draws[s, sigma_idx, drop = TRUE]
      mu_s    <- alpha_s + beta_s * cond_c
      
      pars_s <- list(sigma = sigma_s)
      if (length(extra_idx)) {
        if (!is.null(extra_idx$rate))  pars_s$rate  <- posterior_draws[s, extra_idx$rate,  drop = TRUE]
        if (!is.null(extra_idx$Delta)) pars_s$Delta <- posterior_draws[s, extra_idx$Delta, drop = TRUE]
      }
      
      ll_s <- spec$ll_fun(x = data$RT, pars = pars_s, mu = mu_s)
      
      # Small control for weird cases: i.e., Delta == x in shifted-lognormal
      if(any(is.na(ll_s))) { ll_s[is.na(ll_s)] <- Inf }
      
      log_Ztilde_s      <<- matrixStats::logSumExp(-ll_s)
      log_sum_invZtilde <<- logaddexp_scalar(log_sum_invZtilde, -log_Ztilde_s)
      log_sum_pi_i      <<- logaddexp_vec(log_sum_pi_i,        -log_Ztilde_s - ll_s)
      
      if (s %% update_every == 0L || s == S) {
        # Mensaje opcional (aparece en el handler si lo soporta)
        p(message = sprintf("Draw %d/%d", s, S))
      } else {
        p(amount = 1)  # avanza sin mensaje
      }
    }
  }
  
  if (isTRUE(progress_bar)) {
    progressr::handlers(progressr::handler_progress(format = ":spin :message [:bar] :percent in :elapsed ETA: :eta"))
    progressr::with_progress(runner())
  } else {
    runner()
  }
  
  elpd_i  <- unname(log_sum_invZtilde - log_sum_pi_i)
  elpd    <- sum(elpd_i)
  se_elpd <- sqrt(N * stats::var(elpd_i))
  
  list(
    elpd    = elpd,
    se_elpd = se_elpd,
    elpd_i  = elpd_i,
    draws   = S,
    N       = N,
    model   = model
  )
}

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 7: analytical effect and spearman adjustment formula
# ─────────────────────────────────────────────────────────────────────────────

# Experimental effects estimation via method of moments
analytical.effect <- function(data, 
                              subject_var = "subject", 
                              task_var = "task", 
                              condition_var = "condition", 
                              rt_var = "RT"){
  # Estimate response time mean per subject, tasks and condition
  mean.rt.IJK <- tapply(data[,rt_var], INDEX = list(data[,subject_var], 
                                                    data[,task_var], 
                                                    data[,condition_var]), FUN = mean)
  # Compute experimental effect via method of moments
  effect <- mean.rt.IJK[,,2] - mean.rt.IJK[,,1]
  
  # Transform int
  return(effect)
}

# Classic Spearman correction and reliability estimates
# based on: https://github.com/PerceptionAndCognitionLab/ctx-inhibition/blob/public/papers/rev3/lib.R
spearman.adj <- function(data, 
                         subject_var = "subject", 
                         task_var = "task", 
                         condition_var = "condition", 
                         rt_var = "RT"){
  
  # Estimate response time mean per subject, task and condition
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
  
  # Reliability estimates
  reliability <- (diag(cov(effects.mm)) - effects.resvar)/ diag(cov(effects.mm))
  
  # Spearman adjusted correlation
  rho_est <- cov2cor(cov(effects.mm) - diag(effects.resvar))
  
  # Return rho and reliability
  return(list(rho_est = rho_est, reliability = reliability))
}

# Spearman correction with SEs via lavaan
spearman.lavaan <- function(data, 
                         subject_var = "subject", 
                         task_var = "task", 
                         condition_var = "condition", 
                         rt_var = "RT"){
  
  # Ensure lavaan package
  ensure_packages("lavaan")
  # Estimate response time mean per subject, taska and condition
  mean.rt.IJK <- tapply(data[,rt_var], INDEX = list(data[,subject_var], 
                                                    data[,task_var], 
                                                    data[,condition_var]), FUN = mean)
  
  # Compute experimental effect via method of moments
  effects.mm <- as.data.frame(mean.rt.IJK[,,2] - mean.rt.IJK[,,1])
  colnames(effects.mm) <- paste0("V", 1:ncol(effects.mm))
    
    # Prepare lavaan syntax
    sint <- paste(
      # 1. Fix factor loadings as sqrt(reliability) * observed standard deviation
      paste0("F", 1:ncol(effects.mm), "=~(", sqrt(reliab), "*",apply(effects.mm, 2, sd),  ")*", 
             colnames(effects.mm), "\n", collapse = "\n"),
      # 2. Fix residual variances as (1 - reliability) * observed variance
      paste0("V", 1:ncol(effects.mm), "~~(", apply(effects.mm, 2, var), "*", 1 - reliab, ")*", 
             "V", 1:ncol(effects.mm), "\n", collapse = "\n"), 
      collapse = "\n")
    
    # Fit lavaan model 
    cfa(sint, effects.mm, estimator = "WLS", std.lv = FALSE)
    
    # De-attenuated correlation matrix
    rho_spearman <- matrix(c(lavInspect(lavaan.fit, what = "std")$psi), ncol = ncol(effects.mm))
    
    # WLS standard errors
    rho_SEs <- diag(ncol(effects.mm))
    rho_SEs[lower.tri(rho_SEs)] <- standardizedsolution(lavaan.fit)[19:33,5]
    rho_SEs[upper.tri(rho_SEs)] <- t(rho_SEs)[upper.tri(rho_SEs)]
    
    # Return values
    return(list(rho_est = rho_spearman,
                rho_SE = rho_SEs)) 
}

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 8: Approximated ELPD in LMM (based on Zewotir & Galpin, 2005)
# ─────────────────────────────────────────────────────────────────────────────

# ELPD LOO approximation (based on Zewotir & Galpin, 2005)
elpd_freq_LMM <- function(fit, J_adj = TRUE) {
  # Observed values
  y <- fit@resp$y
  
  # Fitted values
  y_fit <- fitted(fit)
  
  # Influential values
  H_vals <- hatvalues(fit)
  
  # Level-1 residual variance
  sigma_sq <- sigma(fit)^2 
  
  # e_(-i) = e_i / r_ii  (prediction/conditional residual, p.160).
  # Notice that r_ii = 1 - hat_values
  E_loo_resid <- (y - y_fit) / (1 - H_vals)
  
  # Since Var(e_(i)) = σ_e^2 * r_ii (p.156) and e_(-i) = e_i / r_ii (p.160),
  # Var(e_(-i)) = Var(e_(i)) / r_ii^2 = σ_e^2 * r_ii / r_ii^2 = σ_e^2 / r_ii^2
  V_loo_resid <- sigma_sq / (1 - H_vals)
  
  # Compute expected log-pointwise predictive density
  elpd_i <- dnorm(y,  mean = y - E_loo_resid, sd = sqrt(V_loo_resid), 
                  log = TRUE)
  
  # Add Jacobian adjustment to elpds
  # y = -1/x, so jacobian = 1/x^2, and log-jacobian is -2log(x)
  if(J_adj) { elpd_i <- elpd_i - 2 * log(-1/y) } 
  return(list(elpd = sum(elpd_i), se_elpd = sqrt(var(elpd_i) * length(y)), 
              elpd_i = elpd_i))
}


# ─────────────────────────────────────────────────────────────────────────────