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
    // Psychometric factor loadings   
    vector<lower=0, upper=1>[J] h2;       // Communality
    matrix[J, M] Z;                       // Unconstrained unit-vector values
    // Individual-level latent score
    matrix[J, I] Alpha_tilde;             // Standardized version of Alpha
    matrix[J, I] Theta_tilde;             // Standardized version of Theta
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
    corr_matrix[J] Rho_theta = add_diag(tcrossprod(Lambda_std), Psi_std);
    // Cholesky decomposition of model implied correlation matrix
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
    // Model priors: Communalities
    h2 ~ beta(1, 1);
    // Model priors: unit vectors
    target += -0.5 * dot_self(Z_norm) + xi * sum(log(Z_norm));
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
    // Model implied correlation matrix: random intercepts
    corr_matrix[J] Rho_alpha = multiply_lower_tri_self_transpose(L_R_Alpha);
    // Model implied covariance matrix: random intercepts
    cov_matrix[J] Sigma_alpha = quad_form_diag(Rho_alpha, sd_alpha);
    // Unstandardized factor loadings
    matrix[J, M] Lambda = diag_pre_multiply(sd_theta, Lambda_std);
    // Model-Implied covariance matrix: random slopes
    cov_matrix[J] Sigma_theta = quad_form_diag(Rho_theta, sd_theta);
    // Model-Implied uniqueness variance matrix
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
