args <- commandArgs(TRUE)

# number of samples to test
n_test <- as.numeric(args[1])

# prob of positive test for avg 
# person in sample
beta_1 <- as.numeric(args[2])

# pooling indicator
pool_ind <- args[3]

if (pool_ind == TRUE) {
  
  # pool sizes
  num_pool_vec <- c(2, 5)
  
  # Stan model
  prev_model_known <- cmdstan_model("../stan/pool_prevalence_srs_mrp_known_sens_spec.stan")
  
} else {
  
  # pool sizes
  num_pool_vec <- c(1)
  
  # Stan model
  prev_model_known <- cmdstan_model("../stan/prevalence_srs_mrp_known_sens_spec.stan")
  
}

library(renv)
renv::restore()

library(cmdstanr)
library(posterior)

set.seed(754893)
options(scipen = 999)

source("utils.R")

# true sensitivity/specificity
true_sens <- 0.95
true_spec <- 0.995

# total size of population
n_pop <- 1000

# generate population characteristics #

# race: 1 (white), 2 (black), 3 (hispanic), 4 (other)
n_race <- 4
race <- sample(rep(1:n_race, c(500, 300, 100, 100)))

# age: 1 (21 - 40), 2 (41 - 60), 3 (> 60)
n_age <- 4
age <- sample(rep(1:n_age, c(100, 200, 600, 100)))

pop_cov <- data.frame(race, age)

# get population counts for poststratification table
pop_cov$count <- 1
pop_cov_count <- aggregate(count ~ race + age, data = pop_cov, FUN = length)
pop_cov <- pop_cov[,-ncol(pop_cov)]

# prevalence
pop_cov$true_prev <- invlogit(beta_1)
true_prev_mean <- mean(pop_cov$true_prev)

# create poststratification table  
num_ps <- n_race * n_age

ps_pop <- rep(NA, num_ps)
ps_prev <- rep(NA, num_ps)
ps_ind <- 1
for (age_ind in 1:n_age) {
  for (race_ind in 1:n_race) {
    ps_pop[ps_ind] <- pop_cov_count$count[(pop_cov_count$race == race_ind & 
                                             pop_cov_count$age == age_ind)]
    
    ps_prev[ps_ind] <- mean(pop_cov$true_prev[(pop_cov$race == race_ind & 
                                                 pop_cov$age == age_ind)])
    
    ps_ind <- ps_ind + 1
  }
}


prev_mat <- matrix(NA, nrow = length(num_pool_vec), ncol = 10)
colnames(prev_mat) <- c("true_prevalence", "num_pool",
                        "prev_known_median", "prev_known_mean", 
                        "prev_known_mad", "prev_known_sd",
                        "spin_lower", "spin_upper",
                        "CP_known", "Power_known")
prev_mat_ind <- 1

for (num_pool in num_pool_vec) {
  
  n_samp <- n_test * num_pool
  
  prev_known_median_mat <- matrix(NA, nrow = 500, ncol = 1)
  prev_known_mean_mat <- matrix(NA, nrow = 500, ncol = 1)
  prev_known_mad_mat <- matrix(NA, nrow = 500, ncol = 1)
  prev_known_sd_mat <- matrix(NA, nrow = 500, ncol = 1)
  spin_lower_mat <- matrix(NA, nrow = 500, ncol = 1)
  spin_upper_mat <- matrix(NA, nrow = 500, ncol = 1)
  cp_known_mat <- matrix(NA, nrow = 500, ncol = 1)
  power_known_mat <- matrix(NA, nrow = 500, ncol = 1)
  
  for (i in 1:500) {
    
    # get n_samp samples from full cohort
    samp_pop <- select_srs(pop_cov, n_samp)
    
    if (num_pool == 1) {
      
      # give each sampled subject a test status
      meas_status <- rbinom(n = n_samp,
                            size = 1,
                            prob = sample_positive(samp_pop$true_prev,
                                                   true_sens,
                                                   true_spec))
      
      # model data for Stan
      prev_model_known_data <- list(n_samp = nrow(samp_pop),
                                    status = meas_status,
                                    race = samp_pop$race,
                                    age = samp_pop$age,
                                    num_ps = num_ps,
                                    ps_pop = ps_pop,
                                    coef_prior_scale = 0.5,
                                    beta_mu = beta,
                                    sens = true_sens,
                                    spec = true_spec)
      
    } else {
      
      # randomly assign people to pools
      samp_pop$pool <- sample(rep(1:(nrow(samp_pop)/num_pool), each = num_pool))
      
      # give each sampled subject a test status
      meas_status <- rbinom(n = length(unique(samp_pop$pool)),
                            size = 1,
                            prob = unlist(lapply(split(samp_pop$true_prev, samp_pop$pool),
                                                 pool_positive,
                                                 sens = true_sens,
                                                 spec = true_spec)))
      
      # model data for Stan
      prev_model_known_data <- list(n_pool = length(meas_status),
                                    n_per_pool = num_pool,
                                    n_samp = nrow(samp_pop),
                                    status = meas_status,
                                    race = samp_pop$race,
                                    age = samp_pop$age,
                                    num_ps = num_ps,
                                    ps_pop = ps_pop,
                                    coef_prior_scale = 0.5,
                                    beta_mu = beta,
                                    sens = true_sens,
                                    spec = true_spec)
      
    }
    
    prev_fit_known <- prev_model_known$sample(data = prev_model_known_data, 
                                              seed = 50,
                                              refresh = 0, 
                                              parallel_chains = 4, 
                                              iter_warmup = 10000, 
                                              iter_sampling = 10000,
                                              step_size = 0.1,
                                              adapt_delta = 0.95)
    
    prev_known_results <- as.data.frame(summarise_draws(prev_fit_known$draws("prev_total"), mean, median, sd, mad))
    
    prev_known_mean_mat[i,] <-  prev_known_results[,2]
    prev_known_median_mat[i,] <- prev_known_results[,3]
    prev_known_sd_mat[i,] <-  prev_known_results[,4]
    prev_known_mad_mat[i,] <-  prev_known_results[,5]
    
    spin_cp <- calc_spin_cp(true_prev_mean, prev_fit_known$draws("prev_total"))
    
    spin_lower_mat[i,] <- spin_cp[[1]]
    spin_upper_mat[i,] <- spin_cp[[2]]
    cp_known_mat[i,] <- spin_cp[[3]]
    
    power_known_mat[i,] <- 1 - calc_spin_cp(0, prev_fit_known$draws("prev_total"))[[3]]

  }
  
  prev_mat[prev_mat_ind,] <- cbind(c(true_prev_mean),
                                   num_pool,
                                   apply(prev_known_median_mat, 2, mean_na_rm),
                                   apply(prev_known_mean_mat, 2, mean_na_rm),
                                   apply(prev_known_mad_mat, 2, mean_na_rm),
                                   apply(prev_known_sd_mat, 2, mean_na_rm),
                                   apply(spin_lower_mat, 2, mean_na_rm),
                                   apply(spin_upper_mat, 2, mean_na_rm),
                                   apply(cp_known_mat, 2, mean_na_rm),
                                   apply(power_known_mat, 2, mean_na_rm))
  prev_mat_ind <- prev_mat_ind + 1
  
}

if (pool_ind == TRUE) {
  
  write.csv(prev_mat,
            paste0("../../results/pooling/srs_n_test_", n_test, 
                   "_true_prev_", round(true_prev_mean, 3), "_results.csv"),
            row.names = FALSE)
  
} else {
  
  write.csv(prev_mat,
            paste0("../../results/no_pooling/srs_n_test_", n_test, 
                   "_true_prev_", round(true_prev_mean, 3), "_results.csv"),
            row.names = FALSE)
  
}

