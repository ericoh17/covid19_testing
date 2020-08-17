data {
  int<lower = 0> status;
  int<lower = 0> n_pool;
  int<lower = 1> n_per_pool;
  real<lower = 0> logit_spec_prior_scale;
  real<lower = 0> logit_sens_prior_scale;
  real<lower = 0> logit_prev_prior_scale;
}
parameters {
  // specificity
  real mu_logit_spec;
  real<lower = 0> sigma_logit_spec;
  real<offset = mu_logit_spec, multiplier = sigma_logit_spec> logit_spec;
  // sensitivity
  real mu_logit_sens;
  real<lower = 0> sigma_logit_sens;
  real<offset = mu_logit_sens, multiplier = sigma_logit_sens> logit_sens;
  // prevalence
  real mu_logit_prev;
  real<lower = 0> sigma_logit_prev;
  real<offset = mu_logit_prev, multiplier = sigma_logit_prev> logit_prev;
}
transformed parameters {
  real spec = inv_logit(logit_spec);
  real sens = inv_logit(logit_sens);
  real prev = inv_logit(logit_prev);
}
model {
  //likelihood
  real sample_positive = (1- pow((1 - prev), n_per_pool)) * sens + pow((1 - prev), n_per_pool) * (1 - spec);
  status ~ binomial(n_pool, sample_positive);
  // priors
  // specificity
  logit_spec ~ normal(mu_logit_spec, sigma_logit_spec);
  sigma_logit_spec ~ normal(0, logit_spec_prior_scale);
  mu_logit_spec ~ normal(4, 2);
  // sensitivity
  logit_sens ~ normal(mu_logit_sens, sigma_logit_sens);
  sigma_logit_sens ~ normal(0, logit_sens_prior_scale);
  mu_logit_sens ~ normal(4, 2); 
  // prevalence
  logit_prev ~ normal(mu_logit_prev, sigma_logit_prev);
  sigma_logit_prev ~ normal(0, logit_prev_prior_scale);
  mu_logit_prev ~ normal(-5.5, 2); 
}