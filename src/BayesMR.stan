/**
 * @author Ioan Gabriel Bucur
 * @description This is a Stan implementation of the Bayesian model described in Section 3 of the paper
 * "Inferring the direction of a causal link and estimating its effect via a Bayesian Mendelian randomization approach"
 * by Ioan Gabriel Bucur, Tom Claassen, and Tom Heskes, published in "Statistical Methods in Medical Research", 2019.
 * @url https://journals.sagepub.com/doi/10.1177/0962280219851817
 */

data {
  int N; // number of observations
  int J; // number of genetic variants
  int n; // number of Bernoulli trials
  real v_spike; // variance of the 'spike' component in the prior 
  real v_slab; // variance of the 'slab' component in the prior
  
  // Sufficient statistics (First and second order moments)
  matrix[J + 3, J + 3] SS; 
  
  // Effect allele frequencies
  vector<lower = 0, upper = 1>[J] eaf;
  
  // Estimates for mu_X and mu_Y (Equation 4 in 3.1)
  real mu_X;
  real mu_Y;
}

transformed data {
  real sd_spike = sqrt(v_spike);
  real sd_slab = sqrt(v_slab);
  
  // mean and standard deviation of genetic variants
  vector[J] mG;
  vector[J] sigma_G;

  for (j in 1:J) {
    mG[j] = SS[1, j + 1];
  }
  
  sigma_G = sqrt(eaf .* (1 - eaf) * n);
}


parameters {
  
  // pooled mixture weights for spike-and-slab priors (see 3.4.3)
  real<lower = 0, upper = 1> wgamma;
  real<lower = 0, upper = 1> walpha;
  
  // standardized structural parameters (see 3.4.1)
  real sbeta;
  vector[J] sgamma;
  vector[J] salpha;
  real<lower = 0> skappa_X;
  real skappa_Y;
  
  // standard deviation of X (exposure) and Y (outcome)
  real<lower = 0> sigma_X;
  real<lower = 0> sigma_Y;
}

transformed parameters {
  
  // original structural parameters, before rescaling (see 3.4.1)
  real beta = sbeta * sigma_Y / sigma_X;
  vector[J] gamma = sgamma * sigma_X ./ sigma_G;
  vector[J] alpha = salpha * sigma_Y ./ sigma_G;
  real kappa_X = skappa_X * sigma_X;
  real kappa_Y = skappa_Y * sigma_Y;
}


model {

  // define derived simplified matrices for simplifying expressions
  matrix[J + 3, 2] f;
  matrix[2, 2] sigma;
  matrix[2, 2] S;
  
  /* Spike-and-slab Prior on Structural Parameters */
  for (j in 1:J) {
    target += log_sum_exp(log(wgamma) + normal_lpdf(sgamma[j] | 0, sd_spike), log1m(wgamma) + normal_lpdf(sgamma[j] | 0, sd_slab));
    target += log_sum_exp(log(walpha) + normal_lpdf(salpha[j] | 0, sd_spike), log1m(walpha) + normal_lpdf(salpha[j] | 0, sd_slab));
  }
  target += log_sum_exp(log(0.5) +  normal_lpdf(sbeta | 0, sd_spike), log(0.5) +  normal_lpdf(sbeta | 0, sd_slab));
  
  target += log_sum_exp(log(0.5) +  normal_lpdf(skappa_X | 0, sd_spike), log(0.5) +  normal_lpdf(skappa_X | 0, sd_slab));
  target += log_sum_exp(log(0.5) +  normal_lpdf(skappa_Y | 0, sd_spike), log(0.5) +  normal_lpdf(skappa_Y | 0, sd_slab));
    
  /* Uninformative Gaussian prior on Scale Parameters */
  sigma_X ~ normal(0, 10);
  sigma_Y ~ normal(0, 10);

  // Derive f
  f[1, 1] = - mu_X;
  f[1, 2] = - (mu_Y + beta * mu_X);
  for (j in 1:J) {
    f[j + 1, 1] = - sigma_X * sgamma[j] / sigma_G[j];
    f[j + 1, 2] = - sigma_Y * (salpha[j] + sbeta * sgamma[j]) / sigma_G[j];
  }
  f[J + 2, 1] = 1;
  f[J + 2, 2] = 0;
  f[J + 3, 1] = 0;
  f[J + 3, 2] = 1;

  // Derive sigma (the covariance of X, Y | G)
  sigma[1, 1] = sigma_X^2 * (1 + skappa_X^2);
  sigma[1, 2] = sigma_X * sigma_Y * (sbeta * (1 + skappa_X^2) + skappa_X * skappa_Y);
  sigma[2, 1] = sigma[1, 2];
  sigma[2, 2] = sigma_Y^2 * (1 + sbeta^2 + (skappa_Y + sbeta * skappa_X)^2);
  
  // S is defined in equation (6) in Section 3.3
  S = f' * SS * f;

  // log-likelihood for conditional gaussian model
  target += N * (- 0.5 * (log(4 * pi() * pi()) + log_determinant(sigma) + trace(inverse(sigma) * S)));
}


