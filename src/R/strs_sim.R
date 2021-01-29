args <- commandArgs(TRUE)

# number of samples to test
n_test <- as.numeric(args[1])

# prob of positive test for avg 
# person in sample
beta_0 <- as.numeric(args[2])

# pooling indicator
pool_ind <- args[3]

library(renv)
renv::restore()

library(cmdstanr)
library(posterior)

if (pool_ind == TRUE) {
  
  # pool sizes
  num_pool_vec <- c(5, 10)
  
  # Stan model
  prev_model_known <- cmdstan_model("../stan/strs_pool_mrp.stan")
  
} else {
  
  # pool sizes
  num_pool_vec <- c(1)
  
  # Stan model
  prev_model_known <- cmdstan_model("../stan/strs_mrp.stan")
  
}

set.seed(754893)
options(scipen = 999)

source("utils.R")

# true sensitivity/specificity
true_sens <- 0.95
true_spec <- 0.995

# total size of population
n_pop <- 20000

# generate population characteristics #

# gender: 0 (female), 1 (male)
gender <- sample(rep(c(0,1), c(10000, 10000)))

# race: 1 (white), 2 (black), 3 (hispanic), 4 (other)
n_race <- 4
race <- sample(rep(1:n_race, c(10000, 6000, 2000, 2000)))

# geographic variable
# (e.g. university buildings of residence)
n_geo <- 10
geo <- sample(1:n_geo, n_pop, replace = TRUE,
              prob = seq(0.3, 0.1, length.out = 10))

pop_cov <- data.frame(gender, race, geo)

# get population counts for poststratification table
pop_cov$count <- 1
pop_cov_count <- aggregate(count ~ gender + race + geo, data = pop_cov, FUN = length)
pop_cov <- pop_cov[,-ncol(pop_cov)]

# add id
pop_cov$id <- 1:nrow(pop_cov)

# geographic level covariate
# (e.g. % of students in building that are athletes)
x_geo <- rnorm(n_geo, 50, 20)

# standardize covariate
x_geo_scale <- (x_geo - mean(x_geo))/sd(x_geo)

pop_cov$x_geo <- x_geo
pop_cov$x_geo_scale <- x_geo_scale

# prevalence parameters
beta_1 <- 0.25
beta_2 <- 0.5

true_prev <- beta_0 + beta_1 * gender

race_re <- rnorm(n_race, 0, 0.5)
geo_re <- rnorm(n_geo, 0, 0.5)

for (geo_ind in 1:n_geo) {
  true_prev[geo == geo_ind] <- true_prev[geo == geo_ind] + beta_2 * x_geo_scale[geo_ind] + geo_re[geo_ind]
}

for (race_ind in 1:n_race) {
  true_prev[race == race_ind] <- true_prev[race == race_ind] + race_re[race_ind]
}

pop_cov$true_prev <- invlogit(true_prev)
true_prev_mean <- mean(pop_cov$true_prev)

# create poststratification table  
num_ps <- length(unique(gender)) * n_race * n_geo

ps_pop <- rep(NA, num_ps)
ps_prev <- rep(NA, num_ps)
ps_ind <- 1
for (geo_ind in 1:n_geo) {
  for (race_ind in 1:n_race) {
    for (gender_ind in 0:1) {
      ps_pop[ps_ind] <- pop_cov_count$count[(pop_cov_count$gender == gender_ind & 
                                               pop_cov_count$race == race_ind & 
                                               pop_cov_count$geo == geo_ind)]
      
      ps_prev[ps_ind] <- mean(pop_cov$true_prev[(pop_cov$gender == gender_ind & 
                                                   pop_cov$race == race_ind & 
                                                   pop_cov$geo == geo_ind)])
      
      ps_ind <- ps_ind + 1
    }
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
  
  # calculate number of subjects to sample from
  # each strata
  # proportional stratified sampling samples
  # proportional to the size of the strata
  num_in_strata <- c(table(pop_cov$geo))
  n_samp_strata <- floor((n_samp/n_pop) * num_in_strata)
  
  if (num_pool == 1) {
    
    # if total number of samples less than n_samp,
    # add remaining from random strata
    num_samp <- sum(n_samp_strata)
    
    if (num_samp < n_samp) {
      
      samp_diff <- n_samp - num_samp
      extra_samp_strata <- sample(1:n_geo, samp_diff, prob = num_in_strata/sum(num_in_strata))
      
      for (i in 1:length(extra_samp_strata)) {
        
        n_samp_strata[extra_samp_strata[i]] <- n_samp_strata[extra_samp_strata[i]] + 1
        
      }
    }
    
  } else {
    
    # round number to sample in each strata so that 
    # they are divisible by the number in each pool
    n_samp_strata <- sapply(n_samp_strata, find_divisible, divisor = num_pool)
    
  }
  
  prev_known_median_mat <- matrix(NA, nrow = 500, ncol = 1)
  prev_known_mean_mat <- matrix(NA, nrow = 500, ncol = 1)
  prev_known_mad_mat <- matrix(NA, nrow = 500, ncol = 1)
  prev_known_sd_mat <- matrix(NA, nrow = 500, ncol = 1)
  spin_lower_mat <- matrix(NA, nrow = 500, ncol = 1)
  spin_upper_mat <- matrix(NA, nrow = 500, ncol = 1)
  cp_known_mat <- matrix(NA, nrow = 500, ncol = 1)
  power_known_mat <- matrix(NA, nrow = 500, ncol = 1)
  
  for (i in 1:500) {
    
    # get n_samp_strat samples from each strata
    samp_pop_lst <- mapply(select_stratified_srs, 
                           split(pop_cov, pop_cov$geo),
                           n_samp_strata)
    
    if (num_pool == 1) {
      
      samp_pop <- pop_cov[pop_cov$id %in% unlist(samp_pop_lst),]
      
      # give each sampled subject a test status
      meas_status <- rbinom(n = n_samp,
                            size = 1,
                            prob = sample_positive(samp_pop$true_prev,
                                                   true_sens,
                                                   true_spec))
      
      # model data for Stan
      prev_model_known_data <- list(n_samp = n_samp,
                                    status = meas_status,
                                    gender = samp_pop$gender,
                                    race = samp_pop$race,
                                    n_geo = n_geo,
                                    geo = samp_pop$geo,
                                    x_geo_scale = x_geo_scale,
                                    num_ps = num_ps,
                                    ps_pop = ps_pop,
                                    coef_prior_scale = 0.5,
                                    beta_0_mu = beta_0,
                                    sens = true_sens,
                                    spec = true_spec)
      
    } else {
      
      samp_pop <- data.frame()
      for (j in 1:length(samp_pop_lst)) {
        
        temp_samp_pop <- merge(data.frame(id = samp_pop_lst[[j]]), 
                               pop_cov, 
                               by = "id", 
                               sort = FALSE)
        samp_pop <- rbind(samp_pop, temp_samp_pop)
        
      }
      samp_pop$pool <- rep(1:(nrow(samp_pop)/num_pool), each = num_pool)
      
      # give each sampled pool a test status
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
                                    gender = samp_pop$gender,
                                    race = samp_pop$race,
                                    n_geo = n_geo,
                                    geo = samp_pop$geo,
                                    x_geo_scale = x_geo_scale,
                                    num_ps = num_ps,
                                    ps_pop = ps_pop,
                                    coef_prior_scale = 0.5,
                                    beta_0_mu = beta_0,
                                    sens = true_sens,
                                    spec = true_spec)
      
    }
    
    # MCMC
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
                                   apply(prev_known_median_mat, 2, mean),
                                   apply(prev_known_mean_mat, 2, mean),
                                   apply(prev_known_mad_mat, 2, mean),
                                   apply(prev_known_sd_mat, 2, mean),
                                   apply(spin_lower_mat, 2, mean),
                                   apply(spin_upper_mat, 2, mean),
                                   apply(cp_known_mat, 2, mean),
                                   apply(power_known_mat, 2, mean))
  prev_mat_ind <- prev_mat_ind + 1
  
}


if (pool_ind == TRUE) {
  
  write.csv(prev_mat,
            paste0("../../strs_pooling_n_test_", n_test, 
                   "_true_prev_", round(true_prev_mean, 3), "_results.csv"),
            row.names = FALSE)
  
} else {
  
  write.csv(prev_mat,
            paste0("../../strs_no_pooling_n_test_", n_test, 
                   "_true_prev_", round(true_prev_mean, 3), "_results.csv"),
            row.names = FALSE)
  
}


