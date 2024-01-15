# Install and load necessary packages
install.packages(c("tidyverse", "lme4", "MASS"))
install.packages(c("lme4", "Matrix"), dependencies = TRUE)
library(tidyverse)
library(lme4)


# Function to generate data for Normal Parallel CRT Data for Equal Allocation and equal cluster sizes
# Assumes Exhangeable Correlation Structure
generate_parallel_crt_equal_alloc_data <- function(k, 
                                       m,
                                       icc,
                                       sigma2,
                                       beta_0,
                                       beta_star) {
  # Calculate variance components
  # Find Between Variance
  sigma2_B = sigma2 * icc
  
  #Find Within Variance
  sigma2_W = sigma2 - sigma2_B
  
  # Total Number of Clusters given a Binary Treatment by Cluster
  tot_k = 2*k
  study_n = tot_k*m
  
  # Generate Group/Cluster affiliations
  cluster_rand_seq = sample(1:tot_k, tot_k, replace = FALSE)
  
  trial_grp_membership = factor(rep(cluster_rand_seq, each = m))
  group_trt  = if_else( trial_grp_membership %in% cluster_rand_seq[1:k], 1, 0 )
  
  
  X = tibble( 
    clustername = trial_grp_membership,
    intercept = rep(1, study_n),
    trt = group_trt
  )
  
  random_effect_seq = rnorm( tot_k, 
                             mean = 0,
                             sd = sqrt( sigma2_B )  )
  
  
  random_effect_study_n = rep(random_effect_seq, each = m)
  
  
  # Now generate Predictor Matrix
  X_mat = X %>% select(-clustername)
  X_mat = as.matrix(X_mat)
  
  # Coefficient Vector
  beta_vec = c(beta_0, beta_star)
  
  # Now generate Outcome
  outcome = X_mat %*% beta_vec  + 
    random_effect_study_n + 
    rnorm( study_n, mean = 0,sd = sqrt(sigma2_W) )
    
  sim_data = bind_cols( X, outcome = outcome[,1] )
  
  # Return data
  return(sim_data)
}

# Example usage
set.seed(123)
k=4
m=100
sigma2=5
icc = 0.05
beta_0 = 140
beta_star = 10

sim_data = generate_parallel_crt_equal_alloc_data(
  k = k, 
  m = m,
  icc = icc,
  sigma2 = sigma2,
  beta_0 = beta_0,
  beta_star = beta_star
)

# Fit mixed effect model
mixed_model <- lme4::lmer(outcome ~ trt + (1 | clustername), data = sim_data)

# View model summary
summary(mixed_model)
