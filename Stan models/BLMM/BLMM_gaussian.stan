functions {
  // Log-likelihood partial sum for within-chain parallelization
  real partial_sum_log_lik(array[] int slice_ID, int start, int end, array [,,] real RT,
                           vector Alpha, vector Theta, array [,] int T_subj,
                           int K, vector sigma) {
        // Start partial sum at zero                    
        real partial_lpdf = 0.0;
        // Partial log-likelihood
        for (i in start:end) {         // Loop over subjects in the slice
            for (k in 1:K) {           // Loop over conditions
            // Linear prediction
            real mu_ikl = Alpha[i] + (k - 1.5) * Theta[i];
            // Log-likelihood
            partial_lpdf += normal_lpdf(RT[i, k, 1:T_subj[i, k]] | mu_ikl, sigma[i]);
            }
        }
        return partial_lpdf;
    }
}

data {
    int<lower=1> I;                     // Number of subjects
    int<lower=1> K;                     // Number of experimental conditions
    int<lower=1> L_max;                 // Max number of trials across subjects
    array[I, K] int T_subj;             // Number of trials per subject and condition
    array[I, K, L_max] real RT;         // Trial response times per subject and condition
}

transformed data {
    array[I] int Slice_ID;              // Slice ID for within-chain parallelization
    for (i in 1:I){ Slice_ID[i] = i; }
}

parameters {
    // Population means
    real <lower=0> mu_alpha;            // Average response time
    real mu_theta;                      // Average experimental effect
    real <lower=0> shift_sigma;         // Shift parameter for level-1 standard deviation
    // Population standard deviations
    real <lower=0> sd_alpha;            // Std dev of average response time
    real <lower=0> sd_theta;            // Std dev of experimental effects
    real <lower=0> scale_sigma;         // Scale parameter for level-1 standard deviation
    // Individual-level latent score (NCP)
    vector [I] Alpha_tilde;
    vector [I] Theta_tilde;
    // Individual-level residual std dev parameter
    vector <lower=0> [I] Sigma_tilde;
}

transformed parameters{
  // Non-centered parameterization for all latent variables
  vector [I] Alpha = mu_alpha + sd_alpha * Alpha_tilde;
  vector [I] Theta = mu_theta + sd_theta * Theta_tilde;
  vector [I] sigma = shift_sigma + scale_sigma * Sigma_tilde;
}

model {
    // Model priors: Population means
    mu_alpha    ~ normal(0, 10);
    mu_theta    ~ normal(0, 10);
    shift_sigma ~ normal(0, 10);
    // Model priors: Population standard deviations
    sd_alpha    ~ inv_gamma(0.5, 0.5);
    sd_theta    ~ inv_gamma(0.5, 0.5);
    scale_sigma ~ inv_gamma(0.5, 0.5);
    // Model priors: latent variables
    Alpha_tilde ~ std_normal();
    Theta_tilde ~ std_normal();
    Sigma_tilde ~ std_normal();
    
    // Model log-likelihood: within-chain parallelization
    target += reduce_sum(partial_sum_log_lik, Slice_ID, 1, RT, 
                         Alpha, Theta, T_subj, K, sigma);
}

generated quantities {
    // Standard deviations in log scale
    real logsd_alpha = log(sd_alpha);        
    real logsd_theta = log(sd_theta);
    // Shift and scale parameters in log scale
    real logshift_sigma = log(shift_sigma);
    real logscale_sigma = log(scale_sigma);
    // Expected value and std dev for level-1 standard deviation
    real mu_sigma = shift_sigma + scale_sigma * sqrt(2.0/pi());
    real sd_sigma = sqrt(square(scale_sigma) * (1 - 2.0/pi()));
    // Expected value and std dev for level-1 standard deviation (log scale)
    real logmu_sigma  = log(mu_sigma);
    real logsd_sigma  = log(sd_sigma);
    // Signal-to-noise ratio
    real gamma2 = square(sd_theta) / (square(mu_sigma) + square(sd_sigma));
    real gamma  = sqrt(gamma2);
}
