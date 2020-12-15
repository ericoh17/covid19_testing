data {
  int<lower = 0> n_samp;  
  int<lower = 0, upper = 1> status[n_samp]; 
  vector<lower = 0, upper = 1>[n_samp] gender; 
  int<lower = 1, upper = 4> race[n_samp]; 
  int<lower = 1, upper = 4> age[n_samp]; 
  int<lower = 0> n_clust;  
  int<lower = 1, upper = n_clust> clust[n_samp];  
  vector[n_clust] x_clust_scale; 
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
  real<lower = 0> sigma_age;
  real<lower = 0> sigma_clust;
  vector<multiplier = sigma_race>[4] alpha_race;  
  vector<multiplier = sigma_age>[4] alpha_age;  
  vector<multiplier = sigma_clust>[n_clust] alpha_clust; 
}
model {
  vector[n_samp] prev = inv_logit(beta[1] + 
                                    beta[2] * gender +
                                    beta[3] * x_clust_scale[clust] +
                                    alpha_race[race] + 
                                    alpha_age[age] + 
                                    alpha_clust[clust]
                                    );
  vector[n_samp] sample_positive = prev * sens + (1 - prev) * (1 - spec);
  status ~ bernoulli(sample_positive);
  // random intercept priors
  alpha_race ~ normal(0, sigma_race);
  alpha_age ~ normal(0, sigma_age);
  alpha_clust ~ normal(0, sigma_clust);
  // beta priors
  beta[1] ~ normal(beta_1_mu, 1);
  beta[2] ~ normal(0, 1);
  beta[3] ~ normal(0, 1);  // prior on scaled coef
  // hyperpriors for random intercepts
  sigma_race ~ normal(0, coef_prior_scale);
  sigma_age ~ normal(0, coef_prior_scale);
  sigma_clust ~ normal(0, coef_prior_scale);
}
generated quantities {
  // poststratification
  real prev_total;
  vector[num_ps] prev_ps; 
  int ps_ind;
  ps_ind = 1;
  for (clust_ind in 1:n_clust) {
    for (age_ind in 1:4) {
      for (race_ind in 1:4) {
        for (gender_ind in 0:1) {
          prev_ps[ps_ind] = inv_logit(beta[1] +
                                      beta[2] * gender_ind + 
                                      beta[3] * x_clust_scale[clust_ind] +
                                      alpha_race[race_ind] + 
                                      alpha_age[age_ind] + 
                                      alpha_clust[clust_ind]
                                      );
          ps_ind += 1;
        }
      }
    }
  }
  prev_total = sum(ps_pop .* prev_ps) / sum(ps_pop);
  // posterior predictive samples
  vector[n_samp] y_rep_prev = inv_logit(beta[1] + 
                                        beta[2] * gender +
                                        beta[3] * x_clust_scale[clust] +
                                        alpha_race[race] + 
                                        alpha_age[age] + 
                                        alpha_clust[clust]);
  vector[n_samp] y_rep_pos = y_rep_prev * sens + (1 - y_rep_prev) * (1 - spec);
  int<lower = 0, upper = 1> y_rep[n_samp] = bernoulli_rng(y_rep_pos);
}