data {
  int<lower=0> N_obs;         // Number of time points (T)
  int<lower=0> D;             // Maximum delay (D)
  array[N_obs, D] int<lower=0> Y;   // Reported cases (T x D matrix)
  //vector[D] q;                // Reporting probabilities (q[1], q[2], ..., q[D])
}

parameters {
  real<lower=0> alpha_lambda;  // Gamma prior shape parameter for lambda
  real<lower=0> beta_lambda;   // Gamma prior rate parameter for lambda
  vector<lower=0>[N_obs] lambda_t;   // Poisson intensities (λ[t]) at each time point
  vector<lower=0>[N_obs] b_t;             // Parameter for delay function q(d)
}

model {
  // Priors
  alpha_lambda ~ uniform(0, 10);
  beta_lambda ~ uniform(0, 10);
  b_t ~ uniform(0, 10);  // Prior for parameter b
  
  // Gamma prior on Poisson intensities (lambda_t)
  lambda_t ~ gamma(alpha_lambda, beta_lambda);
  
  // Likelihood: Marginalized Poisson likelihood for N_t and Binomial for Y
  for (t in 1:N_obs) {
    for (d in 1:D) {
      real q_d = 1 - exp(-b_t[t] * d);  // q(d) as a function of delay d
      Y[t,d] ~ poisson(lambda_t[t] * q_d);
    }
  }
}


