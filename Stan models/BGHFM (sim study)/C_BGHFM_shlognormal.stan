functions {
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
    array[I, K, J] real <lower=0> RT_min; // Smallest RT for each pearson in each condition and task
    matrix[J, M] Lambda_index;            // Indicator matrix for factor loadings (0 or 1)
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
    vector[J] mu_theta;                   // Average experimental effect per task
    vector<lower=0>[J] shift_sigma;       // Shift parameter for level-1 standard deviation
    array[J] real <lower=0, upper=UL_mu> mu_delta; // Average shift parameter
    // Population standard deviations
    vector<lower=0>[J] sd_alpha;          // Std dev of average reaction time per task
    vector<lower=0>[J] sd_theta;          // Std dev of experimental effects per task
    vector<lower=0>[J] scale_sigma;       // Scale parameter for level-1 standard deviation
    vector<lower=0>[J] sd_delta;          // Std dev of (unconstrained) shift-parameter distribution
    // Cholesky factor decomposition
    cholesky_factor_corr[J] L_R_Alpha;    // Cholesky decomposition for R_Alpha
    cholesky_factor_corr[M] L_Phi;        // Cholesky factor of Phi (common factors correlation matrix)
    // Psychometric factor loadings
    matrix<lower=-1,upper=1>[J,M] L_raw;  // Standardized factor loadings (needs extra constraint)
    // Individual-level latent score
    matrix[J, I] Alpha_tilde;             // Standardized version of Alpha
    matrix[J, I] Theta_tilde;             // Standardized version of Theta
    matrix<lower=0>[J, I] sigma_tilde;    // Level-1 residual variance 
    array[I, J] real <lower=0, upper=UL_d> Delta; // Shift parameter
}

transformed parameters {
    // Fix unrelated variables to zero
    matrix[J, M] Lambda_b = L_raw .* Lambda_index;
    // Model-implied psychometric common-factors correlation matrix
    corr_matrix[M] Phi_b = (M == 1) ? rep_matrix(1.0, 1, 1) : multiply_lower_tri_self_transpose(L_Phi);
    // Model-Implied (proportion of) uniqueness variance 
    vector<lower=0, upper=1>[J] Psi_std =  1 - diagonal(quad_form(Phi_b, Lambda_b'));
    // Model-implied correlation matrix
    corr_matrix[J] Rho_theta = add_diag(quad_form(Phi_b, Lambda_b'), Psi_std);
    // Cholesky decomposition of model-implied correlation matrix
    cholesky_factor_corr[J] L_R_Theta = cholesky_decompose(Rho_theta);
    
    // Compute non-centered individual-level latent scores 
    matrix[I, J] Alpha = rep_matrix(mu_alpha', I) + (diag_pre_multiply(sd_alpha, L_R_Alpha) * Alpha_tilde)';
    matrix[I, J] Theta = rep_matrix(mu_theta', I) + (diag_pre_multiply(sd_theta, L_R_Theta) * Theta_tilde)';
    matrix[I, J] sigma = rep_matrix(shift_sigma', I) + diag_pre_multiply(scale_sigma, sigma_tilde)';
}

model {
    // Model priors: Population means
    mu_alpha    ~ normal(0, 5);
    mu_theta    ~ std_normal();
    shift_sigma ~ normal(0, 5);
    mu_delta    ~ std_normal();
    // Model priors: Population standard deviations
    sd_alpha    ~ student_t(4, 0, 0.5);
    sd_theta    ~ student_t(4, 0, 0.5);
    scale_sigma ~ student_t(4, 0, 0.5);
    sd_delta    ~ student_t(4, 0, 0.5);
    // Model priors: Cholesky factor decomposition
    L_R_Alpha ~ lkj_corr_cholesky(1);
    L_Phi ~ lkj_corr_cholesky(1);
    // Model priors: Factor loadings
    to_vector(L_raw) ~ scaled_beta(1, 1);
    // Model priors: Latent variables
    to_vector(Alpha_tilde) ~ std_normal();
    to_vector(Theta_tilde) ~ std_normal();
    to_vector(sigma_tilde) ~ std_normal();
    // Model priors: shift parameters
    for(j in 1:J){
      target += normal_lpdf(Delta[,j] | mu_delta[j], sd_delta[j]); 
    }
    
    // Model Log-Likelihood
    for (i in 1:I) {
        for (j in 1:J) {
            for (k in 1:K) {
                // Linear prediction
                real mu_ijk = Alpha[i, j] + Theta[i, j] * (k - 1.5);
                // Log-likelihood
                RT[i, k, j, 1:T_subj[i, k, j]] - Delta[i, j] ~ lognormal(mu_ijk, sigma[i,j]);
            }
        }
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
    corr_matrix[J] Rho_alpha = multiply_lower_tri_self_transpose(L_R_Alpha);
    // Model implied covariance matrix: random intercepts
    cov_matrix[J] Sigma_alpha = quad_form_diag(Rho_alpha, sd_alpha);
    // Unstandardized factor loadings
    matrix[J, M] Lambda = diag_pre_multiply(sd_theta, Lambda_std);
    // Model-Implied covariance matrix: random slopes
    cov_matrix[J] Sigma_theta = quad_form_diag(Rho_theta, sd_theta);
    // Model-Implied uniqueness variance
    vector[J] Psi = Psi_std .* square(sd_theta);
    // Compute level-1 residual variance mean and std.dev
    vector[J] mu_sigma = shift_sigma + scale_sigma .* sqrt(2.0/pi());
    vector[J] sd_sigma = scale_sigma .* sqrt(1 - 2.0/pi());
    // Signal-to-noise ratios
    vector[J] gamma2 = square(sd_theta)./(square(mu_sigma) + square(sd_sigma));
    vector[J] gamma  = sqrt(gamma2);
    // Reliability estimate
    vector[J] reliability = gamma2 ./ (gamma2 + 2.0/L_max);
}
