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

true_sens <- 0.95
true_spec <- 0.995

# number of pools to test
n_test <- 200

# vector of pool sizes to consider
num_pool_vec <- seq(5, 30, 5)

# stan model assuming known sens/spec
pool_prev_model_known <- cmdstan_model("src/stan/pool_prevalence_known_sens_spec.stan")

# stan model assuming unknown sens/spec
pool_prev_model_unknown <- cmdstan_model("src/stan/pool_prevalence_unknown_sens_spec.stan")

prev_mat <- matrix(NA, 
                   nrow = length(true_prevalence) * length(num_pool_vec), 
                   ncol = 12)
colnames(prev_mat) <- c("true_prevalence",
                        "num_pool",
                        "prev_known_median", "prev_known_mean", 
                        "prev_known_mad", "prev_known_sd",
                        "CP_known",
                        "prev_unknown_median", "prev_unknown_mean", 
                        "prev_unknown_mad", "prev_unknown_sd",
                        "CP_unknown")

res_ind <- 1

for (p in 1:length(true_prevalence)) {
  
  for (num_pool in num_pool_vec) {
  
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
    
    # give each sampled pool a covid status
    meas_status <- rbinom(n = n_test,
                          size = 1,
                          prob = pool_positive(num_pool,
                                               true_prevalence[p],
                                               true_sens,
                                               true_spec))
    
    # fit model assuming known 
    # sensitivity and specificity
    pool_prev_model_known_data <- list(status = sum(meas_status),
                                       n_pool = n_test,
                                       n_per_pool = num_pool,
                                       sens = true_sens,
                                       spec = true_spec,
                                       logit_prev_prior_scale = 0.5,
                                       mu_logit_prev_mean = logit(true_prevalence[p]))
    
    pool_prev_fit_known <- pool_prev_model_known$sample(data = pool_prev_model_known_data, 
                                                        seed = 954,
                                                        refresh = 0, 
                                                        parallel_chains = 4, 
                                                        iter_warmup = 10000, 
                                                        iter_sampling = 10000,
                                                        step_size = 0.5,
                                                        adapt_delta = 0.9)
    
    pool_prev_known_results <- pool_prev_fit_known$summary() %>% as.data.frame()
    
    pool_prev_known_ci <- spin(pool_prev_fit_known$draws()[,,"prev"], 
                               lower = 0, 
                               upper = 1, 
                               conf = 0.95)
    
    prev_known_mean_vec[i] <-  pool_prev_known_results[5,2]
    prev_known_median_vec[i] <- pool_prev_known_results[5,3]
    prev_known_sd_vec[i] <-  pool_prev_known_results[5,4]
    prev_known_mad_vec[i] <-  pool_prev_known_results[5,5]
    cp_known_vec[i] <- ifelse(true_prevalence[p] > pool_prev_known_ci[1] & 
                                true_prevalence[p] < pool_prev_known_ci[2], 1, 0)
    
    
    # fit model assuming unknown 
    # sensitivity and specificity
    pool_prev_model_unknown_data <- list(status = sum(meas_status),
                                         n_pool = n_test,
                                         n_per_pool = num_pool,
                                         logit_spec_prior_scale = 0.5,
                                         logit_sens_prior_scale = 0.5,
                                         logit_prev_prior_scale = 0.5,
                                         mu_logit_prev_mean = logit(true_prevalence[p]))
    
    pool_prev_fit_unknown <- pool_prev_model_unknown$sample(data = pool_prev_model_unknown_data, 
                                                            seed = 478,
                                                            refresh = 0, 
                                                            parallel_chains = 4, 
                                                            iter_warmup = 10000, 
                                                            iter_sampling = 10000,
                                                            step_size = 0.5,
                                                            adapt_delta = 0.9)
    
    pool_prev_unknown_results <- pool_prev_fit_unknown$summary() %>% as.data.frame()
    
    pool_prev_unknown_ci <- spin(pool_prev_fit_unknown$draws()[,,"prev"], 
                               lower = 0, 
                               upper = 1, 
                               conf = 0.95)
    
    prev_unknown_mean_vec[i] <-  pool_prev_unknown_results[13,2]
    prev_unknown_median_vec[i] <- pool_prev_unknown_results[13,3]
    prev_unknown_sd_vec[i] <-  pool_prev_unknown_results[13,4]
    prev_unknown_mad_vec[i] <-  pool_prev_unknown_results[13,5]
    cp_unknown_vec[i] <- ifelse(true_prevalence[p] > pool_prev_unknown_ci[1] & 
                                  true_prevalence[p] < pool_prev_unknown_ci[2], 1, 0)
    
    }
    
    prev_mat[res_ind,] <- c(true_prevalence[p],
                            num_pool,
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
    
    res_ind <- res_ind + 1
  
  }
  
}

write.csv(prev_mat,
          "pool_prevalence_sim_results.csv",
          row.names = FALSE)


