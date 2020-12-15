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
  num_pool_vec <- c(5, 10)
  
  # Stan model
  prev_model_known <- cmdstan_model("../stan/mcrs_pool_mrp.stan")
  
} else {
  
  # pool sizes
  num_pool_vec <- c(1)
  
  # Stan model
  prev_model_known <- cmdstan_model("../stan/mcrs_mrp.stan")
  
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
n_pop <- 200000

# generate population characteristics #

# gender: 0 (female), 1 (male)
gender <- sample(rep(c(0,1), c(100000, 100000)))

# race: 1 (white), 2 (black), 3 (hispanic), 4 (other)
n_race <- 4
race <- sample(rep(1:n_race, c(100000, 60000, 20000, 20000)))

# age: 1 (< 21), 2 (21 - 40), 3 (41 - 60), 4 (> 60)
n_age <- 4
age <- sample(rep(1:n_age, c(20000, 40000, 120000, 20000)))

# clusters (e.g. census blocks groups)
# assume 0 pop census blocks groups not included
n_clust <- 80
clust <- sample(1:n_clust, n_pop, replace = TRUE)

pop_cov <- data.frame(gender, race, age, clust)

# get population counts for poststratification table
pop_cov$count <- 1
pop_cov_count <- aggregate(count ~ gender + race + age + clust, data = pop_cov, FUN = length)
pop_cov <- pop_cov[,-ncol(pop_cov)]

# add id
pop_cov$id <- 1:nrow(pop_cov)

# geographic level covariate
# (e.g. median income in block group)
x_clust <- rnorm(n_clust, 50000, 20000)

# standardize covariate
x_clust_scale <- (x_clust - mean(x_clust))/sd(x_clust)

pop_cov$x_clust <- x_clust
pop_cov$x_clust_scale <- x_clust_scale

# prevalence parameters
beta_2 <- 0.5
beta_3 <- 0.3

true_prev <- beta_1 + beta_2 * gender

age_re <- rnorm(n_age, 0, 0.5)
race_re <- rnorm(n_race, 0, 0.5)
clust_re <- rnorm(n_clust, 0, 0.5)

for (clust_ind in 1:n_clust) {
  true_prev[clust == clust_ind] <- true_prev[clust == clust_ind] + beta_3 * x_clust_scale[clust_ind] + clust_re[clust_ind]
}

for (age_ind in 1:n_age) {
  true_prev[age == age_ind] <- true_prev[age == age_ind] + age_re[age_ind]
}

for (race_ind in 1:n_race) {
  true_prev[race == race_ind] <- true_prev[race == race_ind] + race_re[race_ind]
}

pop_cov$true_prev <- invlogit(true_prev)
true_prev_mean <- mean(pop_cov$true_prev)

# create poststratification table  
num_ps <- length(unique(gender)) * n_race * n_age * n_clust

ps_pop <- rep(NA, num_ps)
ps_prev <- rep(NA, num_ps)
ps_ind <- 1
for (clust_ind in 1:n_clust) {
  for (age_ind in 1:n_age) {
    for (race_ind in 1:n_race) {
      for (gender_ind in 0:1) {
        ps_pop[ps_ind] <- pop_cov_count$count[(pop_cov_count$gender == gender_ind & 
                                                 pop_cov_count$age == age_ind & 
                                                 pop_cov_count$race == race_ind & 
                                                 pop_cov_count$clust == clust_ind)]
        
        ps_prev[ps_ind] <- mean(pop_cov$true_prev[(pop_cov$gender == gender_ind & 
                                                     pop_cov$age == age_ind &
                                                     pop_cov$race == race_ind & 
                                                     pop_cov$clust == clust_ind)])
        
        ps_ind <- ps_ind + 1
      }
    }
  }
}

prev_mat <- matrix(NA, nrow = length(num_pool_vec), ncol = 10)
colnames(prev_mat) <- c("true_prevalence", "num_pool",
                        "prev_median", "prev_mean", 
                        "prev_mad", "prev_sd",
                        "spin_lower", "spin_upper",
                        "CP", "Power")
prev_mat_ind <- 1

for (num_pool in num_pool_vec) {
  
  n_samp <- n_test * num_pool
  
  # randomly sample clusters
  n_samp_clust <- 10
  clust_samp <- sort(sample(1:n_clust, n_samp_clust, replace = FALSE))
  
  # calculate number of people to sample in each cluster
  n_per_clust <- n_samp/n_samp_clust * c(table(clust_samp))
  
  # get population data from sampled clusters
  samp_clust_lst <- vector("list", n_samp_clust)
  c_ind <- 1
  for (c in clust_samp) {
    
    samp_clust_lst[[c_ind]] <- pop_cov[(pop_cov$clust == c),]
    c_ind <- c_ind + 1
    
  }
  
  pop_cov_samp_clust <- do.call(rbind, samp_clust_lst)
  pop_cov_samp_clust <- pop_cov_samp_clust[order(pop_cov_samp_clust$clust),]
  
  prev_known_median_mat <- matrix(NA, nrow = 500, ncol = 1)
  prev_known_mean_mat <- matrix(NA, nrow = 500, ncol = 1)
  prev_known_mad_mat <- matrix(NA, nrow = 500, ncol = 1)
  prev_known_sd_mat <- matrix(NA, nrow = 500, ncol = 1)
  spin_lower_mat <- matrix(NA, nrow = 500, ncol = 1)
  spin_upper_mat <- matrix(NA, nrow = 500, ncol = 1)
  cp_known_mat <- matrix(NA, nrow = 500, ncol = 1)
  power_known_mat <- matrix(NA, nrow = 500, ncol = 1)
  
  for (i in 1:500) {
    
    # get n_per_clust samples from each cluster
    # stratified by age
    samp_pop_lst <- mapply(select_samp_in_cluster_by_age,
                           split(pop_cov_samp_clust, pop_cov_samp_clust$clust),
                           n_per_clust,
                           SIMPLIFY = FALSE)
    
    if (num_pool == 1) {
      
      samp_pop <- pop_cov_samp_clust[pop_cov_samp_clust$id %in% unlist(samp_pop_lst),]
      
      # give each sampled subject a test status
      meas_status <- rbinom(n = nrow(samp_pop),
                            size = 1,
                            prob = sample_positive(samp_pop$true_prev,
                                                   true_sens,
                                                   true_spec))
      
      # model data for Stan
      prev_model_known_data <- list(n_samp = nrow(samp_pop),
                                    status = meas_status,
                                    gender = samp_pop$gender,
                                    race = samp_pop$race,
                                    age = samp_pop$age,
                                    n_clust = n_clust,
                                    clust = samp_pop$clust,
                                    x_clust_scale = x_clust_scale,
                                    num_ps = num_ps,
                                    ps_pop = ps_pop,
                                    coef_prior_scale = 0.5,
                                    beta_1_mu = beta_1,
                                    sens = true_sens,
                                    spec = true_spec)
      
    } else {
      
      samp_pop <- data.frame()
      
      for (j in 1:length(samp_pop_lst)) {
        for (k in 1:length(samp_pop_lst[[j]])) {
          
          temp_samp_pop <- merge(data.frame(id = samp_pop_lst[[j]][[k]]), 
                                 pop_cov_samp_clust, 
                                 by = "id", 
                                 sort = FALSE)
          samp_pop <- rbind(samp_pop, temp_samp_pop)
          
        }
      }
      
      # assign samples to pools within age strata in cluster
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
                                    age = samp_pop$age,
                                    n_clust = n_clust,
                                    clust = samp_pop$clust,
                                    x_clust_scale = x_clust_scale,
                                    num_ps = num_ps,
                                    ps_pop = ps_pop,
                                    coef_prior_scale = 0.5,
                                    beta_1_mu = beta_1,
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
            paste0("../../results/pooling/mcrs_n_test_", n_test, 
                   "_true_prev_", round(true_prev_mean, 3), "_results.csv"),
            row.names = FALSE)
  
} else {
  
  write.csv(prev_mat,
            paste0("../../results/no_pooling/mcrs_n_test_", n_test, 
                   "_true_prev_", round(true_prev_mean, 3), "_results.csv"),
            row.names = FALSE)

}
