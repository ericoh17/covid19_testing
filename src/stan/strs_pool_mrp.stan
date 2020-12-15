data {
  int<lower = 0> n_pool;  
  int<lower = 1> n_per_pool; 
  int<lower = 0> n_samp; 
  int<lower = 0, upper = 1> status[n_pool]; 
  vector<lower = 0, upper = 1>[n_samp] gender;  
  int<lower = 1, upper = 4> race[n_samp]; 
  int<lower = 0> n_geo; 
  int<lower = 1, upper = n_geo> geo[n_samp];  
  vector[n_geo] x_geo_scale;  
  int<lower = 0> num_ps; 
  vector<lower = 0>[num_ps] ps_pop;
  real<lower = 0> coef_prior_scale;
  real beta_1_mu;
  real<lower = 0, upper = 1> sens;
  real<lower = 0, upper = 1> spec;
}
parameters {
  vector[3] beta; 
  real<lower = 0> sigma_race;
  real<lower = 0> sigma_geo;
  vector<multiplier = sigma_race>[4] alpha_race; 
  vector<multiplier = sigma_geo>[n_geo] alpha_geo; 
}
model {
  vector[n_samp] prev = inv_logit(beta[1] + 
                                    beta[2] * gender +
                                    beta[3] * x_geo_scale[geo] +
                                    alpha_race[race] + 
                                    alpha_geo[geo]
                                    );
  int pool_ind;
  vector[n_pool] no_pos_prob;
  pool_ind = 1;
  for (j in 1:n_pool) {
    no_pos_prob[j] = exp(sum(log(1 - prev[pool_ind:(pool_ind + n_per_pool - 1)])));
    pool_ind = pool_ind + n_per_pool;
  }
  vector[n_pool] sample_positive = ((1 - no_pos_prob) * sens) + (no_pos_prob * (1 - spec));
  status ~ bernoulli(sample_positive);
  // random intercept priors
  alpha_race ~ normal(0, sigma_race);
  alpha_geo ~ normal(0, sigma_geo);
  // beta priors
  beta[1] ~ normal(beta_1_mu, 1);
  beta[2] ~ normal(0, 1);
  beta[3] ~ normal(0, 1); 
  // hyperpriors for random intercepts
  sigma_race ~ normal(0, coef_prior_scale);
  sigma_geo ~ normal(0, coef_prior_scale);
}
generated quantities {
  // poststratification
  real prev_total;
  vector[num_ps] prev_ps; 
  int ps_ind;
  ps_ind = 1;
  for (geo_ind in 1:n_geo) {
    for (race_ind in 1:4) {
      for (gender_ind in 0:1) {
        prev_ps[ps_ind] = inv_logit(beta[1] +
                                      beta[2] * gender_ind + 
                                      beta[3] * x_geo_scale[geo_ind] +
                                      alpha_race[race_ind] + 
                                      alpha_geo[geo_ind]
        );
        ps_ind += 1;
      }
    }
  }
  prev_total = sum(ps_pop .* prev_ps) / sum(ps_pop);
}