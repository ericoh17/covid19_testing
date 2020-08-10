
library(cmdstanr)
library(rstan)
library(dplyr)
library(tidyr)

# function for calculating the probability of a
# sample testing positive
sample_positive <- function(prev, sens, spec) {
  
  (sens * prev) + ((1 - spec) * (1 - prev))
  
}


# vector of true prevalences to consider
#true_prevalence <- c(seq(0.001, 0.01, 0.001), seq(0.012, 0.05, 0.002))
true_prevalence <- 0.1

true_sens <- 0.95
true_spec <- 0.995

# number of samples to test
# can be changed as desired
n_samp <- 200

# give each sampled subject a covid status
meas_status <- rbinom(n = n_samp,
                       size = 1,
                       prob = sample_positive(true_prevalence,
                                              true_sens,
                                              true_spec))

prev_model <- cmdstan_model("src/stan/prevalence.stan")

prev_model_data <- list(status = sum(meas_status), 
                        n_samp = n_samp,
                        logit_spec_prior_scale = 0.5,
                        logit_sens_prior_scale = 0.5)

prev_fit <- prev_model$sample(data = prev_model_data, 
                              refresh = 0, 
                              parallel_chains = 4, 
                              iter_warmup = 1e4, 
                              iter_sampling = 1e4)

rstan::read_stan_csv(prev_fit$output_files())

