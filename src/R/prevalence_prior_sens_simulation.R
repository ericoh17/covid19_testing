library(cmdstanr)
library(rstan)
library(dplyr)
library(tidyr)
library(posterior)
library(bayesplot)

options(scipen = 999)

# function for calculating the probability of a
# sample testing positive
sample_positive <- function(prev, sens, spec) {
  
  (sens * prev) + ((1 - spec) * (1 - prev))
  
}

# shortest posterior interval (spin) from Gelman, Carpenter 2020
spin <- function(x, lower = NULL, upper = NULL, conf = 0.95){
  x <- sort(as.vector(x))
  if (!is.null(lower)) {
    if (lower > min(x)) stop("lower bound is not lower than all the data")
    else x <- c(lower, x)
  }
  if (!is.null(upper)) {
    if (upper < max(x)) stop("upper bound is not higher than all the data")
    else x <- c(x, upper)
  }
  n <- length(x)
  gap <- round(conf*n)
  width <- x[(gap+1):n] - x[1:(n-gap)]
  index <- min(which(width==min(width)))
  x[c(index, index + gap)]
}

logit <- function(p) { log(p/(1-p)) }

# vector of true prevalences to consider
true_prevalence <- c(0.001, 0.0025, 0.005, 0.0075, 
                     0.01, 0.02, 0.03, 0.04, 0.05)

# parameter for shifting the prior distribution
# for prevalence. Can change to be larger or smaller
perc_shift <- -0.8

true_sens <- 0.95
true_spec <- 0.995

# number of samples to test
n_samp <- 200

# stan model assuming known sens/spec
prev_model_known <- cmdstan_model("src/stan/prevalence_known_sens_spec.stan")

# stan model assuming unknown sens/spec
prev_model_unknown <- cmdstan_model("src/stan/prevalence_unknown_sens_spec.stan")

prev_mat <- matrix(NA, 
                   nrow = length(true_prevalence), 
                   ncol = 11)
colnames(prev_mat) <- c("true_prevalence",
                        "prev_known_median", "prev_known_mean", 
                        "prev_known_mad", "prev_known_sd",
                        "CP_known",
                        "prev_unknown_median", "prev_unknown_mean", 
                        "prev_unknown_mad", "prev_unknown_sd",
                        "CP_unknown")
prev_mat[,1] <- true_prevalence

for (p in 1:length(true_prevalence)) {
  
  prev_known_median_vec <- rep(NA, 500)
  prev_known_mean_vec <- rep(NA, 500)
  prev_known_mad_vec <- rep(NA, 500)
  prev_known_sd_vec <- rep(NA, 500)
  prev_unknown_median_vec <- rep(NA, 500)
  prev_unknown_mean_vec <- rep(NA, 500)
  prev_unknown_mad_vec <- rep(NA, 500)
  prev_unknown_sd_vec <- rep(NA, 500)
  cp_known_vec <- rep(NA, 500)
  cp_unknown_vec <- rep(NA, 500)
  
  for (i in 1:500) {
    
    # give each sampled subject a covid status
    meas_status <- rbinom(n = n_samp,
                          size = 1,
                          prob = sample_positive(true_prevalence[p],
                                                 true_sens,
                                                 true_spec))
    
    # fit model assuming known 
    # sensitivity and specificity
    prev_model_known_data <- list(status = sum(meas_status),
                                  n_samp = n_samp,
                                  sens = true_sens,
                                  spec = true_spec,
                                  logit_prev_prior_scale = 0.5,
                                  mu_logit_prev_mean = logit(true_prevalence[p]) + perc_shift)
    
    prev_fit_known <- prev_model_known$sample(data = prev_model_known_data, 
                                              seed = 50,
                                              refresh = 0, 
                                              parallel_chains = 4, 
                                              iter_warmup = 10000, 
                                              iter_sampling = 10000,
                                              step_size = 0.1,
                                              adapt_delta = 0.9)
    
    prev_known_results <- prev_fit_known$summary() %>% as.data.frame()
    
    prev_known_ci <- spin(prev_fit_known$draws()[,,"prev"], 
                          lower = 0, 
                          upper = 1, 
                          conf = 0.95)
    
    prev_known_mean_vec[i] <-  prev_known_results[5,2]
    prev_known_median_vec[i] <- prev_known_results[5,3]
    prev_known_sd_vec[i] <-  prev_known_results[5,4]
    prev_known_mad_vec[i] <-  prev_known_results[5,5]
    cp_known_vec[i] <- ifelse(true_prevalence[p] > prev_known_ci[1] & 
                                true_prevalence[p] < prev_known_ci[2], 1, 0)
    
    
    # fit model assuming unknown
    # sensitivity and specificity
    prev_model_unknown_data <- list(status = sum(meas_status),
                                    n_samp = n_samp,
                                    logit_spec_prior_scale = 0.5,
                                    logit_sens_prior_scale = 0.5,
                                    logit_prev_prior_scale = 0.5,
                                    mu_logit_prev_mean = logit(true_prevalence[p]) + perc_shift)
    
    prev_fit_unknown <- prev_model_unknown$sample(data = prev_model_unknown_data, 
                                                  seed = 4829,
                                                  refresh = 0, 
                                                  parallel_chains = 4, 
                                                  iter_warmup = 10000, 
                                                  iter_sampling = 10000,
                                                  step_size = 0.1,
                                                  adapt_delta = 0.9)
    
    prev_unknown_results <- prev_fit_unknown$summary() %>% as.data.frame()
    
    prev_unknown_ci <- spin(prev_fit_unknown$draws()[,,"prev"], 
                            lower = 0, 
                            upper = 1, 
                            conf = 0.95)
    
    prev_unknown_mean_vec[i] <-  prev_unknown_results[13,2]
    prev_unknown_median_vec[i] <- prev_unknown_results[13,3]
    prev_unknown_sd_vec[i] <-  prev_unknown_results[13,4]
    prev_unknown_mad_vec[i] <-  prev_unknown_results[13,5]
    cp_unknown_vec[i] <- ifelse(true_prevalence[p] > prev_unknown_ci[1] & 
                                  true_prevalence[p] < prev_unknown_ci[2], 1, 0)
    
  }
  
  prev_mat[p,] <- c(true_prevalence[p],
                    mean(prev_known_median_vec), 
                    mean(prev_known_mean_vec), 
                    mean(prev_known_mad_vec),
                    mean(prev_known_sd_vec),
                    mean(cp_known_vec),
                    mean(prev_unknown_median_vec), 
                    mean(prev_unknown_mean_vec), 
                    mean(prev_unknown_mad_vec),
                    mean(prev_unknown_sd_vec),
                    mean(cp_unknown_vec))
  
}

write.csv(prev_mat,
          paste0("prevalence_prior_shift_", perc_shift, "_sim_results.csv"),
          row.names = FALSE)


