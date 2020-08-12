data {
  int<lower = 0> status;
  int<lower = 0> n_samp;
  real<lower = 0> logit_spec_prior_scale;
  real<lower = 0> logit_sens_prior_scale;
}
parameters {
  real<lower = 0, upper = 1> prev;
  real mu_logit_spec;
  real mu_logit_sens;
  real<lower = 0> sigma_logit_spec;
  real<lower = 0> sigma_logit_sens;
  real<offset = mu_logit_spec, multiplier = sigma_logit_spec> logit_spec;
  real<offset = mu_logit_sens, multiplier = sigma_logit_sens> logit_sens;
}
transformed parameters {
  real spec = inv_logit(logit_spec);
  real sens = inv_logit(logit_sens);
}
model {
  //likelihood
  real sample_positive = prev * sens + (1 - prev) * (1 - spec);
  status ~ binomial(n_samp, sample_positive);
  // priors
  logit_spec ~ normal(mu_logit_spec, sigma_logit_spec);
  logit_sens ~ normal(mu_logit_sens, sigma_logit_sens);
  sigma_logit_spec ~ normal(0, logit_spec_prior_scale);
  sigma_logit_sens ~ normal(0, logit_sens_prior_scale);
  mu_logit_spec ~ normal(4, 2); 
  mu_logit_sens ~ normal(4, 2); 
}