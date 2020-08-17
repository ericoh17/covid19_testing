library(cmdstanr)
library(rstan)
library(dplyr)
library(tidyr)
library(posterior)
library(bayesplot)

# function for calculating probability of a pool testing positive
# assume pool is positive if at least one sample in the pool
# is positive
pool_positive <- function(num_pool, prev, sens, spec) {
  
  ((1 - (1 - prev)^num_pool) * sens) + (((1 - prev)^num_pool) * (1 - spec))
  
}


# vector of true prevalences to consider
#true_prevalence <- c(seq(0.001, 0.01, 0.001), seq(0.012, 0.05, 0.002))
true_prevalence <- 0.005

true_sens <- 0.95
true_spec <- 0.995

# number of pools to test
# can be changed as desired
n_test <- 200

# vector of pool sizes to consider
#num_pool_vec <- seq(2, 30, 2)
num_pool <- 5

# give each sampled pool a covid status
meas_status <- rbinom(n = n_test,
                      size = 1,
                      prob = pool_positive(num_pool,
                                           true_prevalence,
                                           true_sens,
                                           true_spec))

# fit model assuming known 
# sensitivity and specificity
pool_prev_model_known <- cmdstan_model("src/stan/pool_prevalence_known_sens_spec.stan")

pool_prev_model_known_data <- list(status = sum(meas_status),
                                   n_pool = n_test,
                                   n_per_pool = num_pool,
                                   sens = true_sens,
                                   spec = true_spec,
                                   logit_prev_prior_scale = 0.5)

pool_prev_fit_known <- pool_prev_model_known$sample(data = pool_prev_model_known_data, 
                                                    seed = 954,
                                                    refresh = 0, 
                                                    parallel_chains = 4, 
                                                    iter_warmup = 10000, 
                                                    iter_sampling = 10000,
                                                    step_size = 0.5,
                                                    adapt_delta = 0.9)

pool_prev_known_results <- pool_prev_fit_known$summary() %>% as.data.frame()

mcmc_hist(pool_prev_fit_known$draws("prev"))

# fit model assuming unknown 
# sensitivity and specificity
pool_prev_model_unknown <- cmdstan_model("src/stan/pool_prevalence_unknown_sens_spec.stan")

pool_prev_model_unknown_data <- list(status = sum(meas_status),
                                     n_pool = n_test,
                                     n_per_pool = num_pool,
                                     logit_spec_prior_scale = 0.25,
                                     logit_sens_prior_scale = 0.25,
                                     logit_prev_prior_scale = 0.25)

pool_prev_fit_unknown <- pool_prev_model_unknown$sample(data = pool_prev_model_unknown_data, 
                                                        seed = 478,
                                                        refresh = 0, 
                                                        parallel_chains = 4, 
                                                        iter_warmup = 10000, 
                                                        iter_sampling = 10000,
                                                        step_size = 0.5,
                                                        adapt_delta = 0.9)

pool_prev_unknown_results <- pool_prev_fit_unknown$summary() %>% as.data.frame()

mcmc_hist(pool_prev_fit_unknown$draws("prev"))

color_scheme_set("gray")
mcmc_scatter(pool_prev_fit_unknown$draws(), 
             pars = c("sens", "prev"), 
             size = 1.5, alpha = 0.5)

mcmc_scatter(pool_prev_fit_unknown$draws(), 
             pars = c("spec", "prev"), 
             size = 1.5, alpha = 0.5)
