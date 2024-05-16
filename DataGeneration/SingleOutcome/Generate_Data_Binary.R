
library(tidyverse)
library(lme4)

logistic_inverse = function(log_odds) {
  return( exp(log_odds) / ( 1 + exp(log_odds) ) )
}

logistic_function = function( p ) {
  return( log(p / (1-p)) )
}

# Function to generate data for Parallel CRT Data for Equal Allocation and equal cluster sizes
# Generates data for binary outcome
# Assumes Exhangeable Correlation Structure
# Citation to Rabideau and Wang (2021)
generate_parallel_crt_equal_alloc_binary_data <- function(k,
                                                          m,
                                                          icc,
                                                          sigma2,
                                                          beta_0,
                                                          beta_star) {

  # Calculate variance components
  sigma2_B = sigma2 * icc
  sigma2_W = sigma2 - sigma2_B

  tot_k = 2 * k
  study_n = tot_k * m

  cluster_rand_seq = sample(1:tot_k, tot_k, replace = FALSE)

  trial_grp_membership = factor(rep(cluster_rand_seq, each = m))
  group_trt  = if_else( trial_grp_membership %in% cluster_rand_seq[1:k], 1, 0)

  X = tibble(
    clustername = trial_grp_membership,
    intercept = rep(1, study_n),
    trt = group_trt
  )

  random_effect_seq = rnorm(tot_k, mean = 0, sd = sqrt(sigma2_B))
  random_effect_study_n = rep(random_effect_seq, each = m)

  # Now generate Predictor Matrix
  X_mat = as.matrix(dplyr::select(X, -clustername))

  # Coefficient Vector
  beta_vec = c(beta_0, beta_star)

  # Now generate Binary Outcome
  log_odds_y = X_mat %*% beta_vec + random_effect_study_n
  prob_y = logistic_inverse( log_odds_y )
  outcome = rbinom(study_n, size = 1, prob = prob_y)

  sim_data = bind_cols(X, outcome = outcome)

  # Return data
  return(sim_data)
}



# Example usage
set.seed(123)
k=10
m=200
icc = 0.05

control_prob_event = 0.3
beta_0 = logistic_function(control_prob_event)

beta_star_odds_ratio = 1.2
beta_star = log(beta_star_odds_ratio)

sim_data = generate_parallel_crt_equal_alloc_binary_data(
  k = k,
  m = m,
  icc = icc,
  sigma2 = .001,
  beta_0 = beta_0,
  beta_star = beta_star
)


# Fit GLMM
model <- glmer(outcome ~ trt + (1 | clustername),
               family = binomial, data = sim_data)

# Summarize the model
summary(model)


exp( summary(model)$coefficients )


