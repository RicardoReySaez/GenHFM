# ╔════════════════════════════════════════════════════════════════════════════╗
# ║                             SCRIPT OVERVIEW                                ║
# ╠════════════════════════════════════════════════════════════════════════════╣
# ║ Script Name   : meta_analysis_functions.R                                  ║
# ║ Author        : Ricardo Rey-Sáez                                           ║
# ║ Role          : PhD Student in Psychology                                  ║
# ║ Institution   : Autonomous University of Madrid, Spain                     ║
# ║ Email         : ricardoreysaez95@gmail.es                                  ║
# ║ Date          : 08-04-2025                                                 ║
# ╚════════════════════════════════════════════════════════════════════════════╝

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 1: Description
# ─────────────────────────────────────────────────────────────────────────────

# Script with user-defined functions to run meta-analysis:
#   1. Clean empirical datasets
#   2. Fit experiment-level bayesian linear mixed model
#   3. Run meta-analysis with parameter estimates

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 2: Load Packages
# ─────────────────────────────────────────────────────────────────────────────

library(MASS)           # Multivariate gaussian distribution
library(dplyr)          # Stan array data
library(cmdstanr)       # cmdstanr
library(posterior)      # cmdstanr summaries
library(infinitefactor) # MatchAlign algorithm functions
library(lavaan)         # Latent variable analysis: Spearman adjustment

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 3: Data cleaning functions
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

# Solve issue with subjects
# retrieved from: https://github.com/jstbcs/acdc-paper/blob/main/share/helper_functions.R
# 20-12-2024
issue.sub <- function(dats){
  tmp <- table(dats$subject, dats$congruency)
  tmp <- cbind(tmp, as.numeric(rownames(tmp)))
  return(tmp[tmp[, 1] < 10 | tmp[, 2] < 10, 3])
}

# Clean data function
# retrieved from: https://github.com/jstbcs/acdc-paper/blob/main/paper/p.Rmd
# 20-12-2024
clean.data <- function(dat){
  dat_sub <- subset(dat, accuracy == 1)
  dat_sub <- subset(dat_sub, rt > .2)
  dat_sub <- subset(dat_sub, rt < 2.5)
  dat_sub <- subset(dat_sub, congruency != "neutral")
  dat_sub$congruency <- ifelse(dat_sub$congruency=="congruent", -.5, .5)
  dat_sub <- subset(dat_sub, block != -999)
  issue <- issue.sub(dat_sub)
  dat_sub <- subset(dat_sub, !(subject %in% issue))
  return(dat_sub[,c("subject", "congruency", "rt")])
}

# Stan Array data (one task)
Stan.array.data <- function(x) {
  library(dplyr)
  
  # Remove missing values
  tmp_data <- na.omit(x)
  
  # Number of subjects and conditions
  n_subj <- length(unique(tmp_data$subject))
  n_cond <- length(unique(tmp_data$condition))
  
  # Determine number of trials within subjects and conditions
  T_subj_all <- tmp_data %>%
    group_by(subject, condition) %>%
    summarize(n_trials = n(), .groups = 'drop')
  
  # Maximum number of trials across all subjects and conditions
  T_max <- max(T_subj_all$n_trials)
  
  # Initialize RT data array for Stan; dimensions = (subject, condition, trial)
  RT <- array(0, dim = c(n_subj, n_cond, T_max))
  T_subj <- matrix(0, nrow = n_subj, ncol = n_cond)
  
  # Initialize arrays for RT_min and RT_max
  RT_min <- matrix(NA, nrow = n_subj, ncol = n_cond)
  RT_max <- matrix(NA, nrow = n_subj, ncol = n_cond)
  
  # Create indices for subjects and conditions
  subject_indices <- sort(unique(tmp_data$subject))
  condition_indices <- sort(unique(tmp_data$condition))
  
  # Map subject and condition to indices
  tmp_data <- tmp_data %>%
    mutate(
      subj_idx = match(subject, subject_indices),
      cond_idx = match(condition, condition_indices)
    )
  
  # Fill T_subj, RT, RT_min, and RT_max arrays
  for (i in 1:n_subj) {
    for (c in 1:n_cond) {
      # Filter data for the current subject and condition
      data_subset <- tmp_data %>%
        filter(subj_idx == i, cond_idx == c)
      
      # Number of trials for the current subject and condition
      n_trials <- nrow(data_subset)
      T_subj[i, c] <- n_trials
      
      if (n_trials > 0) {
        # Extract RT values and store them in the RT array
        RT_values <- data_subset$rt
        RT[i, c, 1:n_trials] <- RT_values
        
        # Calculate and store the minimum and maximum RT values
        RT_min[i, c] <- min(RT_values)
        RT_max[i, c] <- max(RT_values)
      }
    }
  }
  
  # Stan-ready data list
  stan_dat <- list(
    I = n_subj,
    K = n_cond,
    L_max = T_max,
    T_subj = T_subj,
    RT = RT,
    RT_min = RT_min,
    RT_max = RT_max
  )
  
  return(stan_dat)
}

# Informative initial values per model (one task)
make.inits.1task <- function(df, model, congruent_id = -.5, incongruent_id = .5){
  
  # Ensure truncated gaussian distribution package
  ensure_packages("truncnorm")
  
  # Gaussian Copula Function: standardized half-normal distribution
  scale_half_std_normal <- function(x){
    # Scale observed variable
    std_x <- (x - mean(x)) / sd(x)
    # Compute cumulative density values with std gaussian distribution
    p_x <- pnorm(std_x, mean = 0, sd = 1)
    # Generate standardized half-normal values
    y <- truncnorm::qtruncnorm(p = p_x, a = 0, b = Inf, mean = 0, sd = 1)
    # Adjust infinity values to maximum observed values
    if(any(is.infinite(y))) {
      y[which(is.infinite(y))] <- max(y[-which(is.infinite(y))])
    }
    # Return transformed values
    return(y)
  }
  
  # ──────────────────── #
  #    Gaussian model    #
  # ──────────────────── #
  
  if(model == "gaussian"){
    # Subject-level initial values in raw scale
    alpha <- tapply(df$rt, df$subject, mean)
    condition_means <- tapply(df$rt, list(df$subject, df$condition), mean)
    theta <- condition_means[,2] - condition_means[,1]
    sigma <- tapply(df$rt, df$subject, sd)
    
    # Mean and sd of level-1 residual variance
    mu_sigma <- mean(sigma)
    sd_sigma <- sd(sigma)
    
    # Compute offsets and scale values
    scale_sigma <- sd_sigma / sqrt(1 - 2/pi)
    shift_sigma <- mu_sigma - scale_sigma * sqrt(2/pi)
    
    # Gaussian starting values
    initial_values <- list(
      # Population means
      mu_alpha    = mean(alpha),
      mu_theta    = mean(theta),
      shift_sigma = shift_sigma,
      # Population standard deviations
      sd_alpha    = sd(alpha),
      sd_theta    = sd(theta),
      scale_sigma = scale_sigma,
      # Individual-level latent scores
      Alpha_tilde = scale(alpha)[,1],
      Theta_tilde = scale(theta)[,1],
      Sigma_tilde = scale_half_std_normal(sigma)
    )
  }
  
  # ────────────────────── #
  #    Exgaussian model    #
  # ────────────────────── #
  
  if(model == "exgaussian_tau" | model == "exgaussian_rate") {
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
    
    # Subject-level initial values in raw scale
    exgauss_moments <- do.call(rbind, tapply(df$rt, df$subject, mexgauss))
    exgauss_cond_moments <- tapply(df$rt, list(df$subject, df$condition), function(x) mexgauss(x)["mu"])
    alpha <- exgauss_moments[,"mu"]
    theta <- exgauss_cond_moments[,2] - exgauss_cond_moments[,1]
    sigma <- exgauss_moments[,"sigma"]
    tau   <- exgauss_moments[,"tau"]
    lambda  <- 1/exgauss_moments[,"tau"]
    
    # Mean and sd of level-1 residual variance and tau
    mu_sigma   <- mean(sigma)
    sd_sigma   <- sd(sigma)
    mu_tau     <- mean(tau)
    sd_tau     <- sd(tau)
    mu_lambda  <- mean(lambda)
    sd_lambda  <- sd(lambda)
    
    # Compute offsets and scale values
    scale_sigma   <- sd_sigma / sqrt(1 - 2/pi)
    scale_tau     <- sd_tau / sqrt(1 - 2/pi)
    scale_lambda  <- sd_lambda / sqrt(1 - 2/pi)
    shift_sigma   <- mu_sigma - scale_sigma * sqrt(2/pi)
    shift_tau     <- mu_tau - scale_tau * sqrt(2/pi)
    shift_lambda  <- mu_lambda - scale_lambda * sqrt(2/pi)
    
    # Ensure positive shift values
    if(shift_sigma < 0)  { shift_sigma  <- 1e-3 }
    if(shift_tau < 0)    { shift_tau    <- 1e-3 } 
    if(shift_lambda < 0) { shift_lambda <- 1e-3 }
    
    # Exgaussian starting values
    initial_values <- list(
      # Population means
      mu_alpha     = mean(alpha),
      mu_theta     = mean(theta),
      shift_sigma  = shift_sigma,
      shift_tau    = shift_tau,
      shift_lambda = shift_lambda,
      # Population standard deviations
      sd_alpha     = sd(alpha),
      sd_theta     = sd(theta),
      scale_sigma  = scale_sigma,
      scale_tau    = scale_tau,
      scale_lambda = scale_lambda,
      # Individual-level latent scores
      Alpha_tilde = scale(alpha)[,1],
      Theta_tilde = scale(theta)[,1],
      Sigma_tilde = scale_half_std_normal(sigma),
      Tau_tilde   = scale_half_std_normal(tau),
      Lambda_tilde = scale_half_std_normal(lambda)
    )
  }
  
  # ───────────────────────────── #
  #    Shifted-lognormal model    #
  # ───────────────────────────── #
  
  if(model == "shifted.lognormal") {
    # Minimum RT per subject and plausible delta values
    min_RT <- tapply(df$rt, df$subject, min)
    delta <- truncnorm::rtruncnorm(n = length(unique(df$subject)),
                                   a = 0, b = min_RT,
                                   mean = c(min_RT/1.5), sd = .02)
    
    # Empty values
    alpha <- theta <- sigma <- vector(length = length(unique(df$subject)))
    
    # Ensure subjects goes from 1 to I
    df$subject <- as.integer(as.factor(df$subject))
    
    # Subject-level initial values in log scale
    cnt <- 1
    for(i in unique(df$subject)) {
      alpha[cnt] <- mean(log(df$rt[which(df$subject == i)] - delta[i]))
      theta[cnt] <- mean(log(df$rt[which(df$subject == i & df$condition == incongruent_id)] - delta[i])) - 
        mean(log(df$rt[which(df$subject == i & df$condition == congruent_id)] - delta[i]))
      sigma[cnt] <- sd(log(df$rt[which(df$subject==i)] - delta[i]))
      cnt <- cnt + 1
    }
    
    # Mean and sd of level-1 residual variance
    mu_sigma <- mean(sigma)
    sd_sigma <- sd(sigma)
    
    # Compute offsets and scale values
    scale_sigma <- sd_sigma / sqrt(1 - 2/pi)
    shift_sigma <- mu_sigma - scale_sigma * sqrt(2/pi)
    
    # Gaussian starting values
    initial_values <- list(
      # Population means
      mu_alpha    = mean(alpha),
      mu_theta    = mean(theta),
      mu_delta    = mean(delta),
      shift_sigma = shift_sigma,
      # Population standard deviations
      sd_alpha    = sd(alpha),
      sd_theta    = sd(theta),
      sd_delta    = sd(delta),
      scale_sigma = scale_sigma,
      # Individual-level latent scores
      Alpha_tilde = scale(alpha)[,1],
      Theta_tilde = scale(theta)[,1],
      Delta = delta,
      Sigma_tilde = scale_half_std_normal(sigma)
    )
    
  }
  
  # ─────────────────────────── #
  #    Return initial values    #
  # ─────────────────────────── # 
  
  return(initial_values)
}

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 4: Fit bayesian (generalized) linear mixed models
# ─────────────────────────────────────────────────────────────────────────────

# Estimate generalized linear mixed models (gaussian, ex-gaussian and shifted-lognormal) in Stan
Stan_LMM_estimate <- function(data, model, path = "Stan models/Linear Mixed Models/", method = "MCMC"){
  colnames(data) <- c("subject", "condition", "rt")
  sdata <- Stan.array.data(data)
  
  # Gaussian linear mixed model
  if(model == "gaussian"){
    LMM_model <- cmdstan_model(paste0(path, "BLMM_gaussian.stan"), cpp_options = list(stan_threads = TRUE))
    parnames <- c("mu_alpha", "mu_theta", "mu_sigma", "sd_alpha", "sd_theta", 
                  "sd_sigma", "gamma2", "gamma", "logsd_alpha", "logsd_theta", 
                  "logmu_sigma", "logsd_sigma", "shift_sigma", "scale_sigma", 
                  "logshift_sigma", "logscale_sigma")
    # Ex-gaussian linear mixed model
  } else if(model == "exgaussian_tau"){
    LMM_model <- cmdstan_model(paste0(path, "BLMM_exgaussian_tau.stan"), cpp_options = list(stan_threads = TRUE))
    parnames <- c("mu_alpha", "mu_theta", "mu_tau", "mu_sigma", "sd_alpha", 
                  "sd_theta", "sd_tau", "sd_sigma", "gamma2", "gamma", 
                  "logsd_alpha", "logsd_theta", "logmu_tau", "logmu_sigma",
                  "logsd_tau", "logsd_sigma", "shift_tau", "shift_sigma", 
                  "logshift_tau", "logshift_sigma", "scale_tau", "scale_sigma",
                  "logscale_tau", "logscale_sigma")
    # Shifted-lognormal linear mixed model
  } else if(model == "exgaussian_rate"){
    LMM_model <- cmdstan_model(paste0(path, "BLMM_exgaussian_rate.stan"), cpp_options = list(stan_threads = TRUE))
    parnames <- c("mu_alpha", "mu_theta", "mu_lambda", "mu_sigma", "sd_alpha", 
                  "sd_theta", "sd_lambda", "sd_sigma", "gamma2", "gamma", 
                  "logsd_alpha", "logsd_theta", "logmu_lambda", "logmu_sigma",
                  "logsd_lambda", "logsd_sigma", "shift_lambda", "shift_sigma", 
                  "logshift_lambda", "logshift_sigma", "scale_lambda", "scale_sigma",
                  "logscale_lambda", "logscale_sigma")
    # Shifted-lognormal linear mixed model
  } else if(model == "shifted.lognormal"){
    LMM_model <- cmdstan_model(paste0(path, "BLMM_shlognormal.stan"), cpp_options = list(stan_threads = TRUE))
    parnames <- c("mu_alpha", "mu_theta", "mu_sigma", "mu_delta", "sd_alpha", "sd_theta", 
                  "sd_sigma", "sd_delta", "gamma2", "gamma", "logsd_alpha", "logsd_theta", 
                  "logmu_sigma", "logsd_sigma", "logsd_delta", "shift_sigma", "scale_sigma", 
                  "logshift_sigma", "logscale_sigma")
    # Stop if it's not one of this models
  } else {
    stop("Model must be 'gaussian', 'exgaussian_tau', 'exgaussian_rate', or 'shifted.lognormal'")
  }
  
  # Classic bayesian estimate via MCMC draws
  if(method == "MCMC") {
    # Extract samples from the model
    fit_LMM <- LMM_model$sample(
      sdata, 
      iter_sampling = 4000, 
      iter_warmup = 2000, 
      chains = 3,
      parallel_chains = 3,
      threads_per_chain = 4,
      adapt_delta = .8,
      max_treedepth = 10,
      refresh = 50,
      seed = 2025,
      init = replicate(3, make.inits.1task(df = data, model = model), simplify = F)
    )
    
    # Extract posterior means and sds
    posterior_means <- unlist(fit_LMM$summary(parnames)[,"mean"])
    posterior_sds   <- unlist(fit_LMM$summary(parnames)[,"sd"])
    # Check model convergence to ensure valid inferences
    weak.check   <- any(unlist(fit_LMM$summary(parnames)[,"rhat"])>=1.10)
    strong.check <- any(unlist(fit_LMM$summary(parnames)[,"rhat"])>=1.05)
    # Change names
    names(posterior_means) <- paste0(parnames, "_postMean")
    names(posterior_sds)   <- paste0(parnames, "_postSd")
    # Return posterior means, posterior sds and diagnostic summery
    return(c(posterior_means, posterior_sds, unlist(fit_LMM$diagnostic_summary()), 
             weak.rhat = weak.check, strong.rhat = strong.check, 
             avg_rhat  = mean(unlist(fit_LMM$summary(parnames)[,"rhat"])),
             disp_rhat = sd(unlist(fit_LMM$summary(parnames)[,"rhat"])),
             avg_essb  = mean(unlist(fit_LMM$summary(parnames)[,"ess_bulk"])),
             disp_essb = sd(unlist(fit_LMM$summary(parnames)[,"ess_bulk"])),
             avg_esst  = mean(unlist(fit_LMM$summary(parnames)[,"ess_tail"])),
             disp_esst = sd(unlist(fit_LMM$summary(parnames)[,"ess_tail"])),
             time = fit_LMM$time()$total))
  }
  
  # Estimation via laplace posterior approximation
  if(method == "Laplace") {
    # First, estimate maximum a posteriori estimates
    MAP_ests <- LMM_model$optimize(
      sdata, iter = 1e5,
      jacobian = TRUE,
      threads = 10,
      seed = 2025,
      init = list(make.inits.1task(df = data, model = model))
    )

    # Check for maximum iteration limit without converge
    max_iter <- any(grepl(pattern = "  Maximum number of iterations hit, may not be at an optima", 
                          x = capture.output(MAP_ests$output())))
    
        
    # Then, compute model parameters with laplace approximation
    fit_laplace <- LMM_model$laplace(
      data = sdata, 
      mode = MAP_ests, 
      threads = 6,
      draws = 30000)
    # Model summaries
    posterior_means <- unlist(fit_laplace$summary(parnames)[,"mean"])
    posterior_sds   <- unlist(fit_laplace$summary(parnames)[,"sd"])
    # Change names
    names(posterior_means) <- paste0(parnames, "_postMean")
    names(posterior_sds)   <- paste0(parnames, "_postSd")
    # Return posterior means and sds
    return(c(posterior_means, posterior_sds, max_iter = max_iter))
  } 
  
  # Estimation via bayesian variational inference algorithm
  if(method == "ADVI"){
    # Variational inference estimate
    ADVI_est <- LMM_model$variational(
      sdata, 
      iter = 2e5,
      grad_samples = 10,
      elbo_samples = 1000,
      seed = 2025,
      draws = 10000,
      init = list(make.inits.1task(df = data, model = model))
    )
    # Model summaries
    posterior_means <- unlist(ADVI_est$summary(parnames)[,"mean"])
    posterior_sds   <- unlist(ADVI_est$summary(parnames)[,"sd"])
    # Change names
    names(posterior_means) <- paste0(parnames, "_postMean")
    names(posterior_sds)   <- paste0(parnames, "_postSd")
    # Return posterior means and sds
    return(c(posterior_means, posterior_sds))
  }
  
  # Estimation via pathfinder algorithm
  if(method == "pathfinder"){
    # Pathfinder estimates
    pathf_est <- LMM_model$pathfinder(
      sdata, 
      num_paths = 5,
      max_lbfgs_iters = 1e4, 
      history_size = 50,
      init = replicate(5, make.inits.1task(df = data, model = model), simplify = F)
    )
    # Model summaries
    posterior_means <- unlist(pathf_est$summary(parnames)[,"mean"])
    posterior_sds   <- unlist(pathf_est$summary(parnames)[,"sd"])
    # Change names
    names(posterior_means) <- paste0(parnames, "_postMean")
    names(posterior_sds)   <- paste0(parnames, "_postSd")
    # Return posterior means and sds
    return(c(posterior_means, posterior_sds)) 
  }
}

# Function to fit all models in parallel
fit_LMM_models <- function(data, model, method, cores = NULL){
  
  # Ensure future package
  ensure_packages(c("future", "parallel"))
  
  # Number of individual experiments
  n.exps <- length(unique(all.df$dataset_id))
  
  # Prepare parallelized settings
  plan(sequential)
  
  handlers(list(handler_progress(format = "Estimation progress: [:bar] :percent :elapsedfull",
                                 clear = FALSE)))
  
  # Parallelized estimation with progress bar
  with_progress({
    # Create a progress bar with the desired format
    p <- progressor(steps = n.exps)
    
    # Save results
    fit.blmm <- foreach(
      i = 1:n.exps, 
      .combine = "rbind",
      .options.future = list(
        seed = TRUE, 
        packages = c("cmdstanr", "tidyverse"))) %dofuture% {
          # Progress bar
          p()
          
          # Clean data frame
          temp.df <- clean.data(dat = data[which(data$dataset_id == i),])
          
          # Bayesian LMM estimation
          Stan_LMM_estimate(data = temp.df, model = model, method = method)
        }
  })
  
  # Clean caché memory
  gc()
  
  # Return parameter estimates
  return(fit.blmm)
}


# ─────────────────────────────────────────────────────────────────────────────