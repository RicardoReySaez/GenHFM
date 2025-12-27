data {
    int<lower=0> I;                       // Number of subjects
    int<lower=0> J;                       // Number of tasks
    int<lower=0> K;                       // Number of experimental conditions
    int<lower=0> L_max;                   // Max number of trials across subjects
    int<lower=1> M;                       // Number of latent factors
    array[I, K, J] int T_subj;            // Number of trials per subject, condition, and task
    array[I, K, J] vector [L_max] RT;     // Reaction time for each subject, condition, task and trial
    array[I, K, J] real <lower=0> RT_min; // Smallest RT for each pearson in each condition and task
    real<lower=0> xi;                     // Unit-vector repulsion parameter value
    // Prior distribution parameters: means
    matrix[J,2] pr_mu_alpha;
    matrix[J,2] pr_mu_beta;
    matrix[J,2] pr_mu_delta;
    matrix[J,2] pr_shift_sigma;
    // Prior distribution parameters: Sds
    matrix<lower=0>[J,3] pr_sd_alpha;
    matrix<lower=0>[J,3] pr_sd_beta;
    matrix<lower=0>[J,3] pr_sd_delta;
    matrix<lower=0>[J,3] pr_scale_sigma;
    // Prior distribution parameters: R_alpha and h2
    real<lower=0> pr_L_R_alpha;
    vector<lower=0>[2] pr_h2;
    // Mixture-IS to compute LOO-CV?
    int<lower=0, upper=1> loo_mix_IS;
    // Save log_lik in generated quantities (only in small datasets)
    int<lower=0, upper=1> save_log_lik;        
}

transformed data{
  // Minimum RT per subject and task
  array[I, J] real UL_d;    
  // Maximum of min_RT_subj
  array[J] real UL_mu;      
  // Compute both values
  for(j in 1:J){
      for(i in 1:I){
          UL_d[i,j] = min(RT_min[i,,j]);
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
    // Population standard deviations
    vector<lower=0>[J] sd_alpha;          // Std dev of average reaction time per task
    vector<lower=0>[J] sd_beta;           // Std dev of experimental effects per task
    vector<lower=0>[J] scale_sigma;       // Scale parameter for level-1 standard deviation
    vector<lower=0>[J] sd_delta;          // Std dev of (unconstrained) shift-parameter distribution
    // Cholesky factor decomposition
    cholesky_factor_corr[J] L_R_alpha;    // Cholesky decomposition for R_alpha
    // Psychometric factor loadings   
    vector<lower=0, upper=1>[J] h2;       // Communality
    matrix[J, M] Z;                       // Unconstrained unit-vector values
    // Individual-level latent score
    matrix[J, I] alpha_tilde;             // Standardized version of alpha
    matrix[J, I] beta_tilde;              // Standardized version of beta
    matrix<lower=0>[J, I] sigma_tilde;    // Level-1 residual variance 
    array[I, J] real <lower=0, upper=UL_d> Delta; // Shift parameter
}

transformed parameters {
    // Unit vector constraint: compute squared euclidean norm
    vector[J] Z_norm = sqrt(diagonal(tcrossprod(Z))); 
    // Apply unit vector constraint in all J rows
    matrix[J, M] Z_UV = diag_pre_multiply(inv(Z_norm), Z);
    // Compute standardized factor loadings using h2 and unit vector
    matrix[J, M] Lambda_std = diag_pre_multiply(sqrt(h2), Z_UV);
    // Model-Implied (proportion of) uniqueness variance
    vector<lower=0, upper=1>[J] Psi_std = 1 - diagonal(tcrossprod(Lambda_std));
    // Model-implied correlation matrix
    corr_matrix[J] Rho_beta = add_diag(tcrossprod(Lambda_std), Psi_std);
    // Cholesky decomposition of model implied correlation matrix
    cholesky_factor_corr[J] L_R_beta = cholesky_decompose(Rho_beta);
    
    // Compute non-centered individual-level latent scores 
    matrix[I, J] alpha = rep_matrix(mu_alpha', I) + (diag_pre_multiply(sd_alpha, L_R_alpha) * alpha_tilde)';
    matrix[I, J] beta  = rep_matrix(mu_beta', I) + (diag_pre_multiply(sd_beta, L_R_beta) * beta_tilde)';
    matrix[I, J] sigma = rep_matrix(shift_sigma', I) + diag_pre_multiply(scale_sigma, sigma_tilde)';
}

model {
    // Model priors: Population means
    mu_alpha    ~ normal(pr_mu_alpha[,1],    pr_mu_alpha[,2]);
    mu_beta     ~ normal(pr_mu_beta[,1],     pr_mu_beta[,2]);
    shift_sigma ~ normal(pr_shift_sigma[,1], pr_shift_sigma[,2]);
    mu_delta    ~ normal(pr_mu_delta[,1],    pr_mu_delta[,2]);
    // Model priors: Population standard deviations
    sd_alpha    ~ student_t(pr_sd_alpha[,1],    pr_sd_alpha[,2],    pr_sd_alpha[,3]);
    sd_beta     ~ student_t(pr_sd_beta[,1],     pr_sd_beta[,2],     pr_sd_beta[,3]);
    scale_sigma ~ student_t(pr_scale_sigma[,1], pr_scale_sigma[,2], pr_scale_sigma[,3]);
    sd_delta    ~ student_t(pr_sd_delta[,1],    pr_sd_delta[,2],    pr_sd_delta[,3]);
    // Model priors: Cholesky factor decomposition
    L_R_alpha ~ lkj_corr_cholesky(pr_L_R_alpha);
    // Model priors: Communalities
    h2 ~ beta(pr_h2[1], pr_h2[2]);
    // Model priors: unit vectors
    target += -0.5 * dot_self(Z_norm) + xi * sum(log(Z_norm));
    // Model priors: Latent variables
    to_vector(alpha_tilde) ~ std_normal();
    to_vector(beta_tilde)  ~ std_normal();
    to_vector(sigma_tilde) ~ std_normal();
    // Model priors: shift parameters
    for(j in 1:J){
      target += normal_lpdf(Delta[,j] | mu_delta[j], sd_delta[j]); 
    }
    
    // True model estimation
    if(loo_mix_IS == 0) {
        // Model Log-Likelihood
        for (i in 1:I) {
            for (j in 1:J) {
                for (k in 1:K) {
                    // Linear prediction
                    real mu_ijk = alpha[i, j] + beta[i, j] * (k - 1.5);
                    // Log-likelihood
                    RT[i, k, j, 1:T_subj[i, k, j]] - Delta[i, j] ~ lognormal(mu_ijk, sigma[i,j]);
                }
            }
        }
    }

    // Mixture model: only useful for LOO
    if(loo_mix_IS == 1) {
        // Empty log-likelihood vector
        vector[sum(to_array_1d(T_subj))] log_lik;
        int n_id = 1;
        // Model log-likelihood
        for (i in 1:I) {
            for (j in 1:J) {
                for (k in 1:K) {
                    // Linear prediction
                    real mu_ijk = alpha[i, j] + beta[i, j] * (k - 1.5);
                    for(l in 1:T_subj[i, k, j]) {
                        // Save log-likelihood values for each trial
                        log_lik[n_id] = lognormal_lpdf(RT[i, k, j, l] - Delta[i, j] | mu_ijk, sigma[i,j]);
                        n_id += 1;
                    }
                }
            }
        }
        // Add log-likelihood and mixture component
        target += sum(log_lik);
        target += log_sum_exp(-log_lik);
    }
}

generated quantities {
    // Model implied correlation matrix: random intercepts
    corr_matrix[J] Rho_alpha = multiply_lower_tri_self_transpose(L_R_alpha);
    // Model implied covariance matrix: random intercepts
    cov_matrix[J] Sigma_alpha = quad_form_diag(Rho_alpha, sd_alpha);
    // Unstandardized factor loadings
    matrix[J, M] Lambda = diag_pre_multiply(sd_beta, Lambda_std);
    // Model-Implied covariance matrix: random slopes
    cov_matrix[J] Sigma_beta = quad_form_diag(Rho_beta, sd_beta);
    // Model-Implied uniqueness variance matrix
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
                    // Linear prediction
                    real mu_ijk = alpha[i, j] + beta[i, j] * (k - 1.5);
                    for(l in 1:T_subj[i, k, j]) {
                        // Save log-likelihood values for each trial
                        log_lik[n_id] = lognormal_lpdf(RT[i, k, j, l] - Delta[i, j] | mu_ijk, sigma[i,j]);
                        n_id += 1;
                    }
                }
            }
        }
    }
}
