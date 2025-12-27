functions {
  // Log-likelihood partial sum for within-chain parallelization
  real partial_sum_log_lik(array[] int slice_ID, int start, int end, array [,] vector RT,
                           vector Alpha, vector Theta, array [,] int T_subj,
                           int K, vector sigma, array[] real Delta) {
        // Start partial sum at zero                    
        real partial_lpdf = 0.0;
        // Partial log-likelihood
        for (i in start:end) {         // Loop over subjects in the slice
            for (k in 1:K) {           // Loop over conditions
            // Linear prediction
            real mu_ikl = Alpha[i] + (k - 1.5) * Theta[i];
            // Log-likelihood
            partial_lpdf += lognormal_lpdf(RT[i, k, 1:T_subj[i, k]] - Delta[i] | mu_ikl, sigma[i]);
            }
        }
        return partial_lpdf;
    }
}

data {
    int<lower=1> I;                             // Number of subjects
    int<lower=1> K;                             // Number of experimental conditions
    int<lower=1> L_max;                         // Max number of trials across subjects
    array[I, K] int T_subj;                     // Number of trials per subject and condition
    array[I, K] vector [L_max] RT;              // Trial response times per subject and condition
    array[I, K] real <lower=0> RT_min;          // Smallest RT for each pearson in each condition
}

transformed data{
  // Minimum RT per subject (shift parameter upper bound per subject)
  array[I] real ub_RT_i;
  // Slice ID for within-chain parallelization
  array[I] int Slice_ID;
  for(i in 1:I){
      ub_RT_i[i] = min(RT_min[i,]);
      Slice_ID[i] = i;
  }
}

parameters {
    // Population means
    real mu_alpha;                               // Average response time
    real mu_theta;                               // Average experimental effect
    real <lower=0, upper=max(ub_RT_i)> mu_delta; // Average shift parameter
    real <lower=0> shift_sigma;                  // Shift parameter for level-1 standard deviation
    // Population standard deviations
    real <lower=0> sd_alpha;                     // Std dev of average response time
    real <lower=0> sd_theta;                     // Std dev of experimental effects
    real <lower=0> sd_delta;                     // Std dev of shift parameter
    real <lower=0> scale_sigma;                  // Scale parameter for level-1 standard deviation
    // Individual-level latent score (NCP)
    vector [I] Alpha_tilde;
    vector [I] Theta_tilde;
    // Individual-level shift and residual std dev parameter
    array[I] real <lower=0, upper=ub_RT_i> Delta;
    vector <lower=0> [I] Sigma_tilde;
}

transformed parameters{
    // Non-centered parameterization for all latent variables
    vector [I] Alpha  = mu_alpha + sd_alpha * Alpha_tilde;
    vector [I] Theta  = mu_theta + sd_theta * Theta_tilde;
    vector [I] sigma  = shift_sigma + scale_sigma * Sigma_tilde;
}

model {
    // Model priors: Population means
    mu_alpha    ~ normal(0, 10);
    mu_theta    ~ normal(0, 10);
    mu_delta    ~ normal(0, 10);
    shift_sigma ~ normal(0, 10);
    // Model priors: Population standard deviations
    sd_alpha    ~ inv_gamma(0.5, 0.5);
    sd_theta    ~ inv_gamma(0.5, 0.5);
    sd_delta    ~ inv_gamma(0.5, 0.5);
    scale_sigma ~ inv_gamma(0.5, 0.5);
    // Model priors: latent variables
    Alpha_tilde ~ std_normal();
    Theta_tilde ~ std_normal();
    Delta ~ normal(mu_delta, sd_delta);
    Sigma_tilde ~ std_normal();
    
    // Model log-likelihood: within-chain parallelization
    target += reduce_sum(partial_sum_log_lik, Slice_ID, 1, RT, 
                         Alpha, Theta, T_subj, K, sigma, Delta);
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
    real logmu_sigma = log(mu_sigma);
    real logsd_sigma = log(sd_sigma);
    // Non-decision time log-standard deviation
    real logsd_delta = log(sd_delta);
    // Signal-to-noise ratio
    real gamma2 = square(sd_theta) / (square(mu_sigma) + square(sd_sigma));
    real gamma  = sqrt(gamma2);
}
