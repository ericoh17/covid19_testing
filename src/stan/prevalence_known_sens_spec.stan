data {
  int<lower = 0> status;
  int<lower = 0> n_samp;
  real<lower = 0, upper = 1> sens;
  real<lower = 0, upper = 1> spec;
  real<lower = 0> logit_prev_prior_scale;
  real mu_logit_prev_mean;
}
parameters {
  real mu_logit_prev;
  real<lower = 0> sigma_logit_prev;
  real<offset = mu_logit_prev, multiplier = sigma_logit_prev> logit_prev;
}
transformed parameters {
  real<lower = 0, upper = 1> prev = inv_logit(logit_prev);
}
model {
  //likelihood
  real sample_positive = prev * sens + (1 - prev) * (1 - spec);
  status ~ binomial(n_samp, sample_positive);
  // priors
  logit_prev ~ normal(mu_logit_prev, sigma_logit_prev);
  sigma_logit_prev ~ normal(0, logit_prev_prior_scale);
  mu_logit_prev ~ normal(mu_logit_prev_mean, 1); 
}