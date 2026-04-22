functions {
    // Shifted lognormal with a uniform contamination component
    real shiftlnorm_mix_lpdf(real x,
                            real delta,
                             real mu,
                             real sigma,
                             real lambda,
                             real max_RT) {
      if (x <= delta) {
        return log(lambda) + uniform_lpdf(x | 0, max_RT);
      } else {
        return log_mix(
          lambda,
          uniform_lpdf(x | 0, max_RT),
          lognormal_lpdf(x - delta | mu, sigma)
        );
      }
    }
    // Compute sign adjustments for confirmatory models
    vector compute_sign_adjustment(int J, int M, matrix Lambda_b, matrix Lambda_index) {
        vector[M] sign_adjustment = rep_vector(1.0, M);
        array[J] int task_used = rep_array(0, J);

        for (m in 1:M) {
            int j_m = 0;
            for (j in 1:J) {
                if (Lambda_index[j, m] != 0 && task_used[j] == 0) {
                    j_m = j;
                    task_used[j] = 1;
                    break;
                }
            }
            if (j_m > 0 && Lambda_b[j_m, m] < 0) {
                sign_adjustment[m] = -1.0;
            }
        }
        return sign_adjustment;
    }
    // Scaled beta log-probability density function
    real scaled_beta_lpdf(vector y, real alpha, real beta) {
      int N = num_elements(y);
      real log_prob = 0;
      for (n in 1:N) {
        if (y[n] < -1 || y[n] > 1) {
          reject("y must be between -1 and 1, but found y = ", y[n]);
          }
        }
        if (alpha <= 0 || beta <= 0) {
          reject("alpha and beta must be positive, but found alpha = ", alpha, ", beta = ", beta);
        }
        
        log_prob = sum((alpha - 1) * log1p(y) + (beta - 1) * log1m(y)) 
                   - N * (alpha + beta - 1) * log(2) 
                   - N * lbeta(alpha, beta);
        return log_prob;
    }
}

data {
    int<lower=0> I;                       // Number of subjects
    int<lower=0> J;                       // Number of tasks
    int<lower=0> K;                       // Number of experimental conditions
    int<lower=0> L_max;                   // Max number of trials across subjects
    int<lower=1> M;                       // Number of latent factors
    array[I, K, J] int T_subj;            // Number of trials per subject, condition, and task
    array[I, K, J] vector [L_max] RT;     // Reaction time for each subject, condition, task and trial
    array[I, K, J] real <lower=0> RT_max; // Largest RT for each pearson in each condition and task
    matrix[J, M] Lambda_index;            // Indicator matrix for factor loadings (0 or 1)
    // Prior distribution parameters: means
    matrix[J,2] pr_mu_alpha;
    matrix[J,2] pr_mu_beta;
    matrix[J,2] pr_mu_delta;
    matrix[J,2] pr_shift_sigma;
    matrix[J,2] pr_mu_pi;
    // Prior distribution parameters: Sds
    matrix<lower=0>[J,3] pr_sd_alpha;
    matrix<lower=0>[J,3] pr_sd_beta;
    matrix<lower=0>[J,3] pr_sd_delta;
    matrix<lower=0>[J,3] pr_scale_sigma;
    matrix<lower=0>[J,3] pr_sd_pi;
    // Prior distribution parameters: R_alpha and R_Phi
    real<lower=0> pr_L_R_alpha;
    real<lower=0> pr_L_Phi;
    // Prior distribution parameters: factor loadings
    vector[2] pr_L_raw;
    // Mixture-IS to compute LOO-CV?
    int<lower=0, upper=1> loo_mix_IS;
    // Save log_lik in generated quantities (only in small datasets)
    int<lower=0, upper=1> save_log_lik;        
}

transformed data{
  // Maximu, RT per subject and task
  array[I, J] real UL_d;    
  // Maximum of max_RT_subj
  array[J] real UL_mu;      
  // Compute both values
  for(j in 1:J){
      for(i in 1:I){
          UL_d[i,j] = max(RT_max[i,,j]);
      }
      UL_mu[j] = max(UL_d[,j]);
  }
}

parameters {
    // Population means
    vector[J] mu_alpha;                   // Average reaction time per task
    vector[J] mu_beta;                    // Average experimental effect per task
    vector<lower=0>[J] shift_sigma;       // Shift parameter for level-1 standard deviation
    array[J] real <lower=0, upper=UL_mu> mu_delta; // Average shift parameter
    vector[J] mu_pi;                      // Mean contamination (pre-probit transformation)
    // Population standard deviations
    vector<lower=0>[J] sd_alpha;          // Std dev of average reaction time per task
    vector<lower=0>[J] sd_beta;           // Std dev of experimental effects per task
    vector<lower=0>[J] scale_sigma;       // Scale parameter for level-1 standard deviation
    vector<lower=0>[J] sd_delta;          // Std dev of (unconstrained) shift-parameter distribution
    vector<lower=0>[J] sd_pi;             // Std dev of contamination (pre-probit transformation)
    // Cholesky factor decomposition
    cholesky_factor_corr[J] L_R_alpha;    // Cholesky decomposition for R_alpha
    cholesky_factor_corr[M] L_Phi;        // Cholesky factor of Phi (common factors correlation matrix)
    // Psychometric factor loadings
    matrix<lower=-1,upper=1>[J,M] L_raw;  // Standardized factor loadings (needs extra constraint)
    // Individual-level latent score
    matrix[J, I] alpha_tilde;             // Standardized version of alpha
    matrix[J, I] beta_tilde;              // Standardized version of beta
    matrix<lower=0>[J, I] sigma_tilde;    // Level-1 residual variance 
    array[I, J] real <lower=0, upper=UL_d> Delta; // Shift parameter
    matrix[J, I] pi_tilde;                // Subject-task contamination latent z-score
}

transformed parameters {
// Fix unrelated variables to zero
    matrix[J, M] Lambda_b = L_raw .* Lambda_index;
    // Model-implied psychometric common-factors correlation matrix
    corr_matrix[M] Phi_b = (M == 1) ? rep_matrix(1.0, 1, 1) : multiply_lower_tri_self_transpose(L_Phi);
    // Model-Implied (proportion of) uniqueness variance 
    vector<lower=0, upper=1>[J] Psi_std =  1 - diagonal(quad_form(Phi_b, Lambda_b'));
    // Model-implied correlation matrix
    corr_matrix[J] Rho_beta = add_diag(quad_form(Phi_b, Lambda_b'), Psi_std);
    // Cholesky decomposition of model-implied correlation matrix
    cholesky_factor_corr[J] L_R_beta = cholesky_decompose(Rho_beta);
    
    // Compute non-centered individual-level latent scores 
    matrix[I, J] alpha = rep_matrix(mu_alpha', I) + (diag_pre_multiply(sd_alpha, L_R_alpha) * alpha_tilde)';
    matrix[I, J] beta  = rep_matrix(mu_beta', I) + (diag_pre_multiply(sd_beta, L_R_beta) * beta_tilde)';
    matrix[I, J] sigma = rep_matrix(shift_sigma', I) + diag_pre_multiply(scale_sigma, sigma_tilde)';
    
    // Contamination probabilities
    matrix<lower=0, upper=1>[I, J] pi_cont;
    for (j in 1:J) {
      for (i in 1:I) {
        pi_cont[i, j] = Phi_approx(mu_pi[j] + sd_pi[j] * pi_tilde[j, i]);
      }
    }
}

model {
    // Model priors: Population means
    mu_alpha    ~ normal(pr_mu_alpha[,1],    pr_mu_alpha[,2]);
    mu_beta     ~ normal(pr_mu_beta[,1],     pr_mu_beta[,2]);
    shift_sigma ~ normal(pr_shift_sigma[,1], pr_shift_sigma[,2]);
    mu_delta    ~ normal(pr_mu_delta[,1],   pr_mu_delta[,2]);
    mu_pi       ~ normal(pr_mu_pi[,1],       pr_mu_pi[,2]);
    // Model priors: Population standard deviations
    sd_alpha    ~ student_t(pr_sd_alpha[,1],    pr_sd_alpha[,2],    pr_sd_alpha[,3]);
    sd_beta     ~ student_t(pr_sd_beta[,1],     pr_sd_beta[,2],     pr_sd_beta[,3]);
    scale_sigma ~ student_t(pr_scale_sigma[,1], pr_scale_sigma[,2], pr_scale_sigma[,3]);
    sd_delta    ~ student_t(pr_sd_delta[,1],    pr_sd_delta[,2],    pr_sd_delta[,3]);
    sd_pi       ~ student_t(pr_sd_pi[,1],       pr_sd_pi[,2],       pr_sd_pi[,3]);
    // Model priors: Cholesky factor decomposition
    L_R_alpha ~ lkj_corr_cholesky(pr_L_R_alpha);
    L_Phi     ~ lkj_corr_cholesky(pr_L_Phi);
    // Model priors: Factor loadings
    to_vector(L_raw) ~ scaled_beta(pr_L_raw[1], pr_L_raw[2]);
    // Model priors: Latent variables
    to_vector(alpha_tilde) ~ std_normal();
    to_vector(beta_tilde)  ~ std_normal();
    to_vector(sigma_tilde) ~ std_normal();
    to_vector(pi_tilde)    ~ std_normal();
    // Model priors: shift parameters
    for (j in 1:J) {
      for (i in 1:I) {
        Delta[i,j] ~ normal(mu_delta[j], sd_delta[j]);
      }
    }
    
    // True model estimation
    if (loo_mix_IS == 0) {
      for (i in 1:I) {
        for (j in 1:J) {
          for (k in 1:K) {
            real mu_ijk = alpha[i, j] + beta[i, j] * (k - 1.5);
            for (l in 1:T_subj[i, k, j]) {
              target += shiftlnorm_mix_lpdf(
                RT[i, k, j][l] | Delta[i, j], mu_ijk, sigma[i, j], pi_cont[i, j], UL_d[i, j]
              );
            }
          }
        }
      }
    }
  
    // Mixture model: only useful for LOO
    if (loo_mix_IS == 1) {
      vector[sum(to_array_1d(T_subj))] log_lik;
      int n_id = 1;
      for (i in 1:I) {
        for (j in 1:J) {
          for (k in 1:K) {
            real mu_ijk = alpha[i, j] + beta[i, j] * (k - 1.5);
            for (l in 1:T_subj[i, k, j]) {
              log_lik[n_id] = shiftlnorm_mix_lpdf(
                RT[i, k, j][l] | Delta[i, j], mu_ijk, sigma[i, j], pi_cont[i, j], UL_d[i, j]
              );
              n_id += 1;
            }
          }
        }
      }
  
      target += sum(log_lik);
      target += log_sum_exp(-log_lik);
    }
}

generated quantities {
    // Compute sign adjustment
    vector[M] sign_adjustment = compute_sign_adjustment(J, M, Lambda_b, Lambda_index);
    // Adjusted factor loadings
    matrix[J, M] Lambda_std = Lambda_b;
    for (m in 1:M) { Lambda_std[, m] *= sign_adjustment[m]; }
    // Adjusted latent factor correlation matrix
    corr_matrix[M] Phi_lv = Phi_b;
    for(m1 in 1:M){for(m2 in 1:M){ Phi_lv[m1, m2] *= sign_adjustment[m1] * sign_adjustment[m2]; }}
    // Model implied correlation matrix: random intercepts
    corr_matrix[J] Rho_alpha = multiply_lower_tri_self_transpose(L_R_alpha);
    // Model implied covariance matrix: random intercepts
    cov_matrix[J] Sigma_alpha = quad_form_diag(Rho_alpha, sd_alpha);
    // Unstandardized factor loadings
    matrix[J, M] Lambda = diag_pre_multiply(sd_beta, Lambda_std);
    // Model-Implied covariance matrix: random slopes
    cov_matrix[J] Sigma_beta = quad_form_diag(Rho_beta, sd_beta);
    // Model-Implied uniqueness variance
    vector[J] Psi = Psi_std .* square(sd_beta);
    // Compute level-1 residual variance mean and std.dev
    vector[J] mu_sigma = shift_sigma + scale_sigma .* sqrt(2.0/pi());
    vector[J] sd_sigma = scale_sigma .* sqrt(1 - 2.0/pi());
    // Signal-to-noise ratios
    vector[J] gamma2 = square(sd_beta) ./ (square(mu_sigma) + square(sd_sigma));
    vector[J] gamma  = sqrt(gamma2);
    // Reliability estimate
    vector[J] reliability = gamma2 ./ (gamma2 + 2.0/L_max);
    // Save log-likelihood values if required
    array[save_log_lik * sum(to_array_1d(T_subj))] real log_lik;
    if(save_log_lik == 1) {
        int n_id = 1;
        // Model log-likelihood
        for (i in 1:I) {
          for (j in 1:J) {
            for (k in 1:K) {
              real mu_ijk = alpha[i, j] + beta[i, j] * (k - 1.5);
              for (l in 1:T_subj[i, k, j]) {
                log_lik[n_id] = shiftlnorm_mix_lpdf(
                  RT[i, k, j][l] | Delta[i, j], mu_ijk, sigma[i, j], pi_cont[i, j], UL_d[i, j]
                );
              n_id += 1;
            }
          }
        }
      }
    }
}
