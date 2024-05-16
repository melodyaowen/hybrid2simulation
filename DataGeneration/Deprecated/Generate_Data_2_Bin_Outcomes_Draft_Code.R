# Script with data generation functions for two binary endpoints

# logistic_fun() ---------------------------------------------------------------
logistic_fun <- function(p){
  return(log(p / (1-p)))
}

# logistic_inverse_fun() -------------------------------------------------------
logistic_inverse_fun <- function(log_odds){
  return(exp(log_odds) / (1 + exp(log_odds)))
}

# gen_crt_coprimary_data_bin() ------------------------------------------------
# Function to generate data for Normal Parallel CRT Data for Equal Allocation
# and equal cluster sizes
# Assumes Exchangeable Correlation Structure
# Recall notation: j = 1,...,m is index of cluster size
#                  k = 1,...,2K is index of cluster (K is number in trt arm)
#                  q = 1, 2 is index of outcome (Q = 2)
#                  i = 1, 2 is treatment group index
gen_crt_coprimary_data_bin <- function(K, # Number of clusters in treatment arm
                                       m, # Number of individuals in each cluster
                                       beta1, # Int effect for Y1
                                       beta2, # Int effect for Y2
                                       rho01, # ICC for Y1
                                       rho02, # ICC for Y2
                                       rho1, # Corr(Y1, Y2) for two diff individuals
                                       rho2, # Corr(Y1, Y2) for same individual
                                       varY1, # Var(Y1)
                                       varY2, # Var(Y2)
                                       r = 1
                                       ){

  # Variance of treatment assignment
  r_alt <- 1/(r + 1)
  sigmaz.square <- r_alt*(1 - r_alt)

  # Total Number of Clusters given a Binary Treatment by Cluster
  K1 <- ceiling(K)
  K2 <- ceiling(r*K)
  K_total <- K1 + K2
  study_n <- K_total*m

  # Intercepts are assumed to be zero
  beta01 <- 0
  beta02 <- 0

  # Random assignment of clusters to treatment
  groups <- tibble(cluster_id = 1:K_total,
                   trt_group = sample(c(rep(1, K1), rep(0, K2)), replace = FALSE))

  # R is the correlation matrix of Y_k
  R11 <- (1 - rho01)*diag(m) + rho01*matrix(1, nrow = m, ncol = m)
  R22 <- (1 - rho02)*diag(m) + rho02*matrix(1, nrow = m, ncol = m)
  R12 <- (rho2 - rho1)*diag(m) + rho1*matrix(1, nrow = m, ncol = m)
  R <- rbind(cbind(R11, R12), cbind(R12, R22))

  rmvbin(n = 100, margprob, commonprob = diag(margprob),
         bincorr = diag(length(margprob)),
         sigma = diag(length(margprob)),
         colnames = NULL, simulvals = NULL)

  # Make sure matrix is positive definite
  if(is.positive.definite(R) == FALSE){
    stop("R matrix is not positive definite. Check input values for correlations and variances.")
  }

  # Outcome 1 mean
  mu1_k_0 <- theta01 <- exp(beta01)/(exp(beta01) + 1) # Cluster K, treatment 0
  mu1_k_1 <- theta11 <- exp(beta01 + beta1)/
    (exp(beta01 + beta1) + 1) # Cluster K, treatment 1

  # Outcome 2 mean
  mu2_k_0 <- theta02 <- exp(beta02)/(exp(beta02) + 1)
  mu2_k_1 <- theta12 <- exp(beta02 + beta2)/
    (exp(beta02 + beta2) + 1)

  # V matrix for each outcome for each treatment group
  V1_k_0 <- mu1_k_0*(1-mu1_k_0)*diag(1, nrow = m, ncol = m) # outcome 1, trt 0
  V1_k_1 <- mu1_k_1*(1-mu1_k_1)*diag(1, nrow = m, ncol = m) # outcome 1, trt 1
  V2_k_0 <- mu2_k_0*(1-mu2_k_0)*diag(1, nrow = m, ncol = m) # outcome 1, trt 0
  V2_k_1 <- mu2_k_1*(1-mu2_k_1)*diag(1, nrow = m, ncol = m) # outcome 1, trt 1

  # Variance matrices
  Var_Y1_k_0 <- sqrt(V1_k_0) %*% R11 %*% sqrt(V1_k_0)
  Var_Y1_k_1 <- sqrt(V1_k_1) %*% R11 %*% sqrt(V1_k_1)
  Var_Y2_k_0 <- sqrt(V2_k_0) %*% R22 %*% sqrt(V2_k_0)
  Var_Y2_k_1 <- sqrt(V2_k_1) %*% R22 %*% sqrt(V2_k_1)



  # Covariance matrix
  Cov_Y1_Y2_0 <- sqrt(V1_k_0) %*% R12 %*% sqrt(V2_k_0)
  Cov_Y1_Y2_1 <- sqrt(V1_k_1) %*% R12 %*% sqrt(V2_k_1)







}

K = 5 # Number of clusters in treatment arm
m = 100 # Number of individuals in each cluster
beta1 = 0.1 # Int effect for Y1
beta2 = 0.2 # Int effect for Y2
rho01 = 0.03 # ICC for Y1
rho02 = 0.1 # ICC for Y2
rho1 = 0.05 # Corr(Y1, Y2) for two diff individuals
rho2 = 0.2 # Corr(Y1, Y2) for same individual
varY1 = 0.035 # Var(Y1)
varY2 = 0.05 # Var(Y2)
r = 1

k <- 6
m <- 20
icc <- 0.05
sigma2 <- 0.1
beta_0 <- 0
beta_star <- 2
