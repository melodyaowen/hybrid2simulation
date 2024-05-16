# Script with data generation functions for two binary endpoints

# logistic_fun() ---------------------------------------------------------------
logistic_fun <- function(p){
  return(log(p / (1-p)))
}

# logistic_inverse_fun() -------------------------------------------------------
logistic_inverse_fun <- function(log_odds){
  return(exp(log_odds) / (1 + exp(log_odds)))
}

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
                                       r = 1,
                                       prop_Y1, # Proportion of individuals in
                                       # a single cluster with Y1 = 1
                                       prop_Y2 # Proportion of individuals in
                                       # a single cluster with Y2 = 1
){

  # Intercepts are assumed to be zero
  beta01 <- 0
  beta02 <- 0

  # Variance of treatment assignment
  r_alt <- 1/(r + 1)
  sigmaz.square <- r_alt*(1 - r_alt)

  # Total Number of Clusters given a Binary Treatment by Cluster
  K1 <- ceiling(K)
  K2 <- ceiling(r*K)
  K_total <- K1 + K2
  study_n <- K_total*m

  # Random assignment of clusters to treatment
  groups <- tibble(cluster_id = 1:K_total,
                   trt_group = sample(c(rep(1, K1), rep(0, K2)),
                                      replace = FALSE))

  # Calculate the within-cluster component of variance
  sigma2w_1 <- varY1*(1 - rho01) # For outcome 1
  sigma2w_2 <- varY2*(1 - rho02) # For outcome 2

  # Check that the proportion of individuals in a single cluster with Y1 = 1
  # and Y2 = 1 (separately) is valid given the within-cluster component of var
  if(round(prop_Y1*(1 - prop_Y1), 2) != round(sigma2w_1, 2)){
    stop("Check value for the proportion of individuals in a single cluster with Y1 = 1; it does not match up with the calculated within-cluster component of variance calculated from Var(Y1) and rho01.")
  }
  if(round(prop_Y2*(1 - prop_Y2), 2) != round(sigma2w_2, 2)){
    stop("Check value for the proportion of individuals in a single cluster with Y2 = 1; it does not match up with the calculated within-cluster component of variance calculated from Var(Y2) and rho02.")
  }

  # Calculate the between-cluster standard variance
  sigma2b_1 <- (rho01*sigma2w_1)/(1-rho01) # For outcome 1
  sigma2b_2 <- (rho02*sigma2w_2)/(1-rho02) # For outcome 2
  #sigma2b_1 <- varY1*rho01 # gives same thing
  #sigma2b_2 <- varY2*rho02 # gives same thing

  vars = c(varY1, varY2)

  sigmaz.square = sigmaz.square
  m = m
  Q = 2 # Number of outcomes

  # Calculate working correlation matrix
  calCorr <- function(beta1s = c(0, 0),
                      beta2s = c(beta1, beta2),
                      delta2s = c(0,0),
                      rho01 = matrix(c(rho01, rho1,
                                       rho1, rho02),
                                     2, 2),
                      rho2 = matrix(c(1, rho2,
                                      rho2, 1),
                                    2, 2),
                      N = K_total,
                      r = r_alt,
                      m = m,
                      K = 2){

    sigma2ks <- sapply(1:K, function(x) calsigma2ksq(r, m, beta1s[x], beta2s[x], rho01[x,x]))

    meanVector <- sqrt(N)*(beta2s - delta2s)/sqrt(sigma2ks)

    wCor <- diag(K)
    for(k1 in 1:K){
      for(k2 in 1:K){
        if(k1 != k2){
          beta1k1 <- beta1s[k1]
          beta2k1 <- beta2s[k1]
          beta1k2 <- beta1s[k2]
          beta2k2 <- beta2s[k2]
          rho2k1k2 <- rho2[k1,k2]
          rho0k1 <- rho01[k1,k1]
          rho0k2 <- rho01[k2,k2]
          rho1k1k2 <- rho01[k1,k2]
          wCor[k1,k2] <- calCorWk1Wk2(r, m, beta1k1, beta2k1, beta1k2, beta2k2,
                                      rho0k1, rho0k2, rho1k1k2, rho2k1k2)
        }
      }
    }
    return(list(sigma2ks, meanVector, wCor))
  }

  tempVals <- calCorr(beta1s = c(0, 0),
                      beta2s = c(beta1, beta2),
                      delta2s = c(0,0),
                      rho01 = matrix(c(rho01, rho1,
                                       rho1, rho02),
                                     2, 2),
                      rho2 = matrix(c(1, rho2,
                                      rho2, 1),
                                    2, 2),
                      N = K_total,
                      r = r_alt,
                      m = m,
                      K = 2)

  sigma2ks <- tempVals[[1]]
  meanVector <- tempVals[[2]]
  wCorMat <- tempVals[[3]]

  R_k <- constrRi(rho01 = matrix(c(rho01, rho1,
                                 rho1, rho02),
                               2, 2),
                  rho2 = matrix(c(1, rho2,
                                  rho2, 1),
                                2, 2),
                  m = m,
                  K = 2)

  R_11 <- R_k[1:m, 1:m]
  R_12 <- R_k[1:m, (m+1):(2*m)]
  R_22 <- R_k[(m+1):(2*m), (m+1):(2*m)]

  Gamma1 <- calGammaK(r = r_alt, m = m, beta1k = beta01, beta2k = beta1)
  Gamma2 <- calGammaK(r = r_alt, m = m, beta1k = beta02, beta2k = beta2)

  Omega1 <- calOmegak(r = r_alt, m = m,
                      beta1k = beta01, beta2k = beta1, rho0k = rho01)
  Omega2 <- calOmegak(r = r_alt, m = m,
                      beta1k = beta02, beta2k = beta2, rho0k = rho02)

  Sigma1 <- calSigmaK(r = r_alt, m = m,
                      beta1k = beta01, beta2k = beta1, rho0k = rho01)
  Sigma2 <- calSigmaK(r = r_alt, m = m,
                      beta1k = beta02, beta2k = beta2, rho0k = rho02)

  sigma2_1sq <- calsigma2ksq(r = r_alt, m = m,
                             beta1k = beta01, beta2k = beta1, rho0k = rho01)
  sigma2_2sq <- calsigma2ksq(r = r_alt, m = m,
                             beta1k = beta02, beta2k = beta2, rho0k = rho02)

  Omega12 <- calOmegak1k2(r = r_alt, m = m,
                          beta1k1 = beta01, beta2k1 = beta1,
                          beta1k2 = beta02, beta2k2 = beta2,
                          rho1k1k2 = rho1, rho2k1k2 = rho2)

  Cov2Betas <- calCov2Betas(r = r_alt, m = m,
                            beta1k1 = beta01, beta2k1 = beta1,
                            beta1k2 = beta02, beta2k2 = beta2,
                            rho1k1k2 = rho1, rho2k1k2 = rho2)

  sigma2_12sq <- calsigma2k1k2sq(r = r_alt, m = m,
                                 beta1k1 = beta01, beta2k1 = beta1,
                                 beta1k2 = beta02, beta2k2 = beta2,
                                 rho1k1k2 = rho1, rho2k1k2 = rho2)

  CorW1W2 <- calCorWk1Wk2(r = r_alt, m = m,
                          beta1k1 = beta01, beta2k1 = beta1,
                          beta1k2 = beta02, beta2k2 = beta2,
                          rho0k1 = rho01, rho0k2 = rho02,
                          rho1k1k2 = rho1, rho2k1k2 = rho2)
  wCorMat





  R11 <- (1 - rho01)*diag(m) + rho01*matrix(rep(1, m*m), nrow = m, ncol = m)
  R22 <- (1 - rho02)*diag(m) + rho02*matrix(rep(1, m*m), nrow = m, ncol = m)
  R12 <- (rho2 - rho1)*diag(m) + rho1*matrix(rep(1, m*m), nrow = m, ncol = m)
  R <- rbind(cbind(R11, R12), cbind(R12, R22))

  testData <- generate.binary(nObs = 5*2,
                              prop.vec.bin = c(prop_Y1, prop_Y2),
                              corr.mat = R[1:2, 1:2])





  # Generate matrix G, a positive-definite matrix with diagonal entries as
  # the between-cluster variance, and off-diagonal entries as the covariances
  rhob <- 0 # rho_b is described as the between-cluster correlation, not sure
  # what this value would be
  G <- matrix(c(sigma2b_1, # 2x2 matrix because we have Q = 2 outcomes
                rhob*sqrt(sigma2b_1)*sqrt(sigma2b_2),
                rhob*sqrt(sigma2b_1)*sqrt(sigma2b_2),
                sigma2b_2), nrow = 2, ncol = 2)

  # Generate vector of random intercepts for cluster k across Q = 2 endpoints
  b_k <- mvrnorm(n = K_total, mu = rep(0, 2), Sigma = G) %>%
    as.data.frame() %>%
    dplyr::select(b1_k = V1, b2_k = V2) %>%
    mutate(cluster_id = 1:K_total)

  # Outcome 1 mean
  mu1_k_0 <- exp(beta01)/(exp(beta01) + 1) # Cluster K, trt 0
  mu1_k_1 <- exp(beta01 + beta1)/(exp(beta01 + beta1) + 1) # Cluster K, trt 1

  # Outcome 2 mean
  mu2_k_0 <- exp(beta02)/(exp(beta02) + 1)
  mu2_k_1 <- exp(beta02 + beta2)/(exp(beta02 + beta2) + 1)

  # Generate matrix C, a positive definite matrix of correlations for the
  # within-cluster errors e_k
  # off-diagonals can be 0 (this assumes observations from same cluster are
  # conditionally independent given the random effects in univariate and
  # multivariate models and across outcome measures in the multivariate setting)
  R11 <- (1 - rho01)*diag(m) + rho01*matrix(rep(1, m*m), nrow = m, ncol = m)
  R22 <- (1 - rho02)*diag(m) + rho02*matrix(rep(1, m*m), nrow = m, ncol = m)
  R12 <- (rho2 - rho1)*diag(m) + rho1*matrix(rep(1, m*m), nrow = m, ncol = m)
  R <- rbind(cbind(R11, R12), cbind(R12, R22))

  # Generate matrix V, a diagonal matrix of known variance functions
  ## V matrix for each outcome for each treatment group
  V1_k_0 <- mu1_k_0*(1-mu1_k_0)*diag(1, nrow = m, ncol = m) # outcome 1, trt 0
  V1_k_1 <- mu1_k_1*(1-mu1_k_1)*diag(1, nrow = m, ncol = m) # outcome 1, trt 1
  V2_k_0 <- mu2_k_0*(1-mu2_k_0)*diag(1, nrow = m, ncol = m) # outcome 1, trt 0
  V2_k_1 <- mu2_k_1*(1-mu2_k_1)*diag(1, nrow = m, ncol = m) # outcome 1, trt 1

  # Initialize dataset for random errors
  e_kj <- data.frame(e_kj1 = c(),
                     e_kj2 = c(),
                     cluster_id = c())
  for(k in 1:nrow(groups)){ # Loop through each cluster,
    # Specify V_k based on trt group
    if(groups$trt_group[k] == 1){ # Treatment for this cluster is 1
      V_temp <- as.matrix(rbind(cbind(V1_k_1, matrix(0, 1, 1)),
                                cbind(matrix(0, 1, 1), V2_k_1)))
    } else if(groups$trt_group[k] == 0){ # Trt for this cluster is 0
      V_temp <- as.matrix(rbind(cbind(V1_k_0, matrix(0, 1, 1)),
                                cbind(matrix(0, 1, 1), V2_k_0)))
    }

    R_temp <- sqrt(V_temp) %*% R %*% sqrt(V_temp)
    if(is.positive.definite(R_temp) == FALSE){
      stop("R_k matrix is not positive definite. Check input values for correlations and variances.")
    }

    e_temp <- mvrnorm(m, mu = rep(0, 2), Sigma = R_temp) %>%
      as.data.frame() %>%
      dplyr::select(e1_kj = V1, e2_kj = V2) %>%
      mutate(cluster_id = sort(rep(k, m)))

    e_kj <- rbind(e_kj, e_temp)
  }

  # Add subject ID to dataframe
  e_kj <- e_kj %>%
    mutate(subject_id = 1:study_n)

  # Generate dataset based on specifications
  simDat <- tibble(subject_id = 1:study_n,
                   cluster_id = sort(rep(1:K_total, m))) %>%
    left_join(., groups, by = "cluster_id") %>%
    left_join(., b_k, by = "cluster_id") %>%
    left_join(., e_kj, by = c("cluster_id", "subject_id")) %>%
    # Now create outcome columns based on inputs and calculated intercepts
    mutate(log_odds_Y1 = beta1*trt_group + b1_k + e1_kj,
           log_odds_Y2 = beta2*trt_group + b2_k + e2_kj) %>%
    mutate(prob_Y1 = logistic_inverse_fun(log_odds_Y1),
           prob_Y2 = logistic_inverse_fun(log_odds_Y2)) %>%
    mutate(Y1 = rbinom(m*K*2, size = 1, prob = prob_Y1),
           Y2 = rbinom(m*K*2, size = 1, prob = prob_Y2)) %>%
    dplyr::rename(subject_j = subject_id,
                  cluster_k = cluster_id,
                  treatment_z_k = trt_group)

  return(simDat)
}




K = 2 # Number of clusters in treatment arm
m = 5 # Number of individuals in each cluster
beta1 = 0.1 # Trt effect for Y1
beta2 = 0.1 # Trt effect for Y2
rho01 = 0.025 # ICC for Y1
rho02 = 0.025 # ICC for Y2
rho1 = 0.01 # Corr(Y1, Y2) for two diff individuals
rho2 = 0.05 # Corr(Y1, Y2) for same individual
varY1 = 0.23 # Var(Y1)
varY2 = 0.25 # Var(Y2)
r = 1
prop_Y1 = 0.66
prop_Y2 = 0.60



beta1s = c(0, 0)
beta2s = c(beta1, beta2)
delta2s = c(0,0)
rho01 = matrix(c(rho01, rho1,
                 rho1, rho02),
               2, 2)
rho2 = matrix(c(1, rho2,
                rho2, 1),
              2, 2)
N = K_total
r = r_alt
m = m
K = 2
