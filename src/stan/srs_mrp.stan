data {
  int<lower = 0> n_samp;  
  int<lower = 0, upper = 1> status[n_samp]; 
  int<lower = 1, upper = 4> race[n_samp];  /
  int<lower = 1, upper = 4> age[n_samp];  
  int<lower = 0> num_ps;  
  vector<lower = 0>[num_ps] ps_pop;
  real<lower = 0> coef_prior_scale;
  real beta_mu;
  real<lower = 0, upper = 1> sens;
  real<lower = 0, upper = 1> spec;
}
parameters {
  real beta;  
  real<lower = 0> sigma_race;
  real<lower = 0> sigma_age;
  vector<multiplier = sigma_race>[4] alpha_race;  
  vector<multiplier = sigma_age>[4] alpha_age;
}
model {
  vector[n_samp] prev = inv_logit(beta + 
                                  alpha_race[race] + 
                                  alpha_age[age]
                                  );
  vector[n_samp] sample_positive = prev * sens + (1 - prev) * (1 - spec);
  status ~ bernoulli(sample_positive);
  // random intercept priors
  alpha_race ~ normal(0, sigma_race);
  alpha_age ~ normal(0, sigma_age);
  // beta priors
  beta ~ normal(beta_mu, 1);
  // hyperpriors for random intercepts
  sigma_race ~ normal(0, coef_prior_scale);
  sigma_age ~ normal(0, coef_prior_scale);
}
generated quantities {
  // poststratification
  real prev_total;
  vector[num_ps] prev_ps; 
  int ps_ind;
  ps_ind = 1;
  for (age_ind in 1:4) {
    for (race_ind in 1:4) {
      prev_ps[ps_ind] = inv_logit(beta +
                                  alpha_race[race_ind] + 
                                  alpha_age[age_ind]
                                  );
      ps_ind += 1;
    }
  }
  prev_total = sum(ps_pop .* prev_ps) / sum(ps_pop);
  // posterior predictive samples
  vector[n_samp] y_rep_prev = inv_logit(beta + 
                                        alpha_race[race] + 
                                        alpha_age[age]);
  vector[n_samp] y_rep_pos = y_rep_prev * sens + (1 - y_rep_prev) * (1 - spec);
  int<lower = 0, upper = 1> y_rep[n_samp] = bernoulli_rng(y_rep_pos);
}