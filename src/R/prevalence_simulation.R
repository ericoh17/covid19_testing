
library(cmdstanr)
library(rstan)
library(dplyr)
library(tidyr)
library(posterior)
library(bayesplot)

# function for calculating the probability of a
# sample testing positive
sample_positive <- function(prev, sens, spec) {
  
  (sens * prev) + ((1 - spec) * (1 - prev))
  
}


# vector of true prevalences to consider
#true_prevalence <- c(seq(0.001, 0.01, 0.001), seq(0.012, 0.05, 0.002))
true_prevalence <- 0.005

true_sens <- 0.95
true_spec <- 0.995

# number of samples to test
# can be changed as desired
n_samp <- 2000

# give each sampled subject a covid status
meas_status <- rbinom(n = n_samp,
                       size = 1,
                       prob = sample_positive(true_prevalence,
                                              true_sens,
                                              true_spec))

# fit model assuming known 
# sensitivity and specificity
prev_model_known <- cmdstan_model("src/stan/prevalence_known_sens_spec.stan")

prev_model_known_data <- list(status = sum(meas_status),
                              n_samp = n_samp,
                              sens = true_sens,
                              spec = true_spec,
                              logit_prev_prior_scale = 0.5)

prev_fit_known <- prev_model_known$sample(data = prev_model_known_data, 
                                          seed = 5048,
                                          refresh = 0, 
                                          parallel_chains = 4, 
                                          iter_warmup = 10000, 
                                          iter_sampling = 10000,
                                          step_size = 0.1,
                                          adapt_delta = 0.9)

prev_known_results <- prev_fit_known$summary() %>% as.data.frame()

mcmc_hist(prev_fit_known$draws("prev"))


# fit model assuming unknown
# sensitivity and specificity
prev_model_unknown <- cmdstan_model("src/stan/prevalence_unknown_sens_spec.stan")

prev_model_unknown_data <- list(status = sum(meas_status),
                                n_samp = n_samp,
                                logit_spec_prior_scale = 0.25,
                                logit_sens_prior_scale = 0.25,
                                logit_prev_prior_scale = 0.25)

prev_fit_unknown <- prev_model_unknown$sample(data = prev_model_unknown_data, 
                                              seed = 4829,
                                              refresh = 0, 
                                              parallel_chains = 4, 
                                              iter_warmup = 10000, 
                                              iter_sampling = 10000,
                                              step_size = 0.1,
                                              adapt_delta = 0.9)

prev_unknown_results <- prev_fit_unknown$summary() %>% as.data.frame()

mcmc_hist(prev_fit_unknown$draws("prev"))

color_scheme_set("gray")
mcmc_scatter(prev_fit_unknown$draws(), 
             pars = c("sens", "prev"), 
             size = 1.5, alpha = 0.5)

mcmc_scatter(prev_fit_unknown$draws(), 
             pars = c("spec", "prev"), 
             size = 1.5, alpha = 0.5)
