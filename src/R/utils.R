# function for calculating the probability of a
# sample testing positive
sample_positive <- function(prev, sens, spec) {
  
  (sens * prev) + ((1 - spec) * (1 - prev))
  
}

# function for calculating probability of a pool testing positive
# assume pool is positive if at least one sample in the pool
# is positive
pool_positive <- function(prev, sens, spec) {
  
  # (1-p1)(1-p2)... on log scale to avoid overflow
  no_pos_prob <- exp(sum(log(1 - prev)))
  
  ((1 - no_pos_prob) * sens) + (no_pos_prob * (1 - spec))
  
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

# function that returns spin interval
# bounds and whether it covers truth
calc_spin_cp <- function(prev_true, draws) {
  
  ci <- spin(draws, 
             lower = 0, 
             upper = 1, 
             conf = 0.95)
  
  cp <- ifelse(prev_true > ci[1] & prev_true < ci[2], 1, 0)
  
  return(list(ci[1], ci[2], cp))
}

logit <- function(p) { log(p/(1-p)) }

invlogit <- function(x) { exp(x)/(1+exp(x))}

# select srs
select_srs <- function(df, samp_n) {
  
  df[sample(nrow(df), samp_n, replace = FALSE),]
  
}

# helper function for strs
select_stratified_srs <- function(df, samp_n) {
  
  strat_df <- df[sample(nrow(df), samp_n, replace = FALSE),]
  return(strat_df$id)
  
}

# get strs by age
# helper for mcrs
select_samp_in_cluster_by_age <- function(df, samp_n) {
  
  clust_by_age_lst <- mapply(select_stratified_srs, 
                             split(df, df$age),
                             rep(samp_n/4, 4),
                             SIMPLIFY = FALSE)
  
  return(clust_by_age_lst)

}

# find closest integer divisible by divisor
find_divisible <- function(dividend, divisor) {
  
  if (dividend %% divisor == 0) {
    
    return(dividend)
    
  } else {
    
    close_less <- dividend - (dividend %% divisor)
    close_more <- (dividend + divisor) - (dividend %% divisor)
    
    if (dividend - close_less > close_more - dividend) {
      
      return(close_more)
      
    } else {
      
      return(close_less)
      
    }
  }
}

mean_na_rm <- function(x) {
  mean(x, na.rm = TRUE)
}



