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
                                       r = 1#,
                                       #prop_Y1, # Proportion of individuals in
                                       # a single cluster with Y1 = 1
                                       #prop_Y2 # Proportion of individuals in
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

  # Calculate the between-cluster standard variance
  sigma2b_1 <- (rho01*sigma2w_1)/(1-rho01) # For outcome 1
  sigma2b_2 <- (rho02*sigma2w_2)/(1-rho02) # For outcome 2
  #sigma2b_1 <- varY1*rho01 # gives same thing
  #sigma2b_2 <- varY2*rho02 # gives same thing

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
  C11 <- (1 - rho01)*diag(1) + rho01*matrix(1, nrow = 1, ncol = 1)
  C22 <- (1 - rho02)*diag(1) + rho02*matrix(1, nrow = 1, ncol = 1)
  C12 <- (rho2 - rho1)*diag(1) + rho1*matrix(1, nrow = 1, ncol = 1)
  C <- rbind(cbind(C11, C12), cbind(C12, C22))

  # Generate matrix V, a diagonal matrix of known variance functions
  ## V matrix for each outcome for each treatment group
  V1_k_0 <- mu1_k_0*(1-mu1_k_0)*diag(1, nrow = 1, ncol = 1) # outcome 1, trt 0
  V1_k_1 <- mu1_k_1*(1-mu1_k_1)*diag(1, nrow = 1, ncol = 1) # outcome 1, trt 1
  V2_k_0 <- mu2_k_0*(1-mu2_k_0)*diag(1, nrow = 1, ncol = 1) # outcome 1, trt 0
  V2_k_1 <- mu2_k_1*(1-mu2_k_1)*diag(1, nrow = 1, ncol = 1) # outcome 1, trt 1

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

    R_temp <- sqrt(V_temp) %*% C %*% sqrt(V_temp)
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
} # End gen_crt_coprimary_data_bin()

# GEE.est() --------------------------------------------------------------------
# Function to estimate dataset statistics, takes in generated binary dataset
# from gen_crt_coprimary_data_bin()
GEE.est <- function(myData){

  ## Main dataset stats --------------------------------------------------------

  # Cluster size of inputted data frame
  m_sim <- nrow(dplyr::filter(myData, cluster_k == 1))

  # Number of clusters in treatment group
  K1 <- dplyr::select(myData, cluster_k, treatment_z_k) %>%
    unique() %>%
    filter(treatment_z_k == 1) %>%
    nrow()
  K_Treatment_sim <- K1

  # Number of clusters in control group
  K2 <- dplyr::select(myData, cluster_k, treatment_z_k) %>%
    unique() %>%
    filter(treatment_z_k == 0) %>%
    nrow()

  # Total K
  K_Total_sim <- K1 + K2

  # Treatment groups
  groups <- myData %>%
    dplyr::select(cluster_k, treatment_z_k) %>%
    unique()

  # Estimate of outcome specific variance
  varY1_sim <- var(myData$Y1)
  varY2_sim <- var(myData$Y2)

  # Probability of a cluster receiving the experimental intervention
  z_bar <- K1/(K1 + K2)

  # Treatment Allocation Ratio, recall r = K2/K1 where K1 = # in treatment group
  r_sim <- nrow(filter(myData, treatment_z_k == 0))/
    nrow(filter(myData, treatment_z_k == 1))

  ## Method 1: Separate GEE models using geepack::geeglm() ---------------------

  # Calculating outcome specific ICC
  m1 <- geepack::geeglm(formula = Y1 ~ treatment_z_k, # creating a geeglm model
                        id = factor(cluster_k),
                        family = binomial(link = "logit"),
                        data = myData, corstr = "exchangeable")
  m1_summary <- summary(m1) # summary output for m1
  rho01_geepack <- m1_summary[12]$corr[[1]]
  beta1_intercept_geepack <- m1_summary$coefficients[1,1]
  beta1_geepack <- m1_summary$coefficients[2,1]

  m2 <- geepack::geeglm(formula = Y2 ~ treatment_z_k, # creating a geeglm model
                        id = factor(cluster_k),
                        family = binomial(link = "logit"),
                        data = myData, corstr = "exchangeable")
  m2_summary <- summary(m2) # summary output for m1
  rho02_geepack <- m2_summary[12]$corr[[1]]
  beta2_intercept_geepack <- m2_summary$coefficients[1,1]
  beta2_geepack <- m2_summary$coefficients[2,1]

  ## Method 2: Using mmm package -----------------------------------------------

  # Fitting generalized multivariate linear mixed model (mmm package)
  fit <- mmm(formula = cbind(Y1, Y2) ~ treatment_z_k,
             id = myData$cluster_k, data = myData,
             family = binomial(link = logit), corStruct = "exchangeable")

  beta1_mmm <- round(fit$coefficients[["Y1.treatment_z_k"]], 4)
  beta2_mmm <- round(fit$coefficients[["Y2.treatment_z_k"]], 4)

  beta1_intercept_mmm <- round(fit$coefficients[["Y1.Intercept"]], 4)
  beta2_intercept_mmm <- round(fit$coefficients[["Y2.Intercept"]], 4)

  wCorr_mmm <- fit$working.correlation # dim (m*2) x (m*2) # R in D.Li Paper
  C11_mmm <- wCorr_mmm[1:m_sim, 1:m_sim]
  C12_mmm <- wCorr_mmm[1:m_sim, (m_sim+1):(2*m_sim)]
  C21_mmm <- wCorr_mmm[(m_sim+1):(2*m_sim), 1:m_sim]
  C22_mmm <- wCorr_mmm[(m_sim+1):(2*m_sim), (m_sim+1):(2*m_sim)]

  # Estimates for rho1 and rho2
  rho1_mmm <- C12_mmm[1,2]
  rho2_mmm <- C12_mmm[1,1]

  # ICC estimates
  rho01_mmm <- C11_mmm[1,2]
  rho02_mmm <- C22_mmm[1,2]

  ## Method 3: Direct GEE calculations using D. Li equations -------------------
  # Summing outcomes per cluster
  sumDat <- myData %>%
    group_by(cluster_k) %>%
    dplyr::summarize(Y1_sum = sum(Y1),
                     Y2_sum = sum(Y2),
                     trt = sum(treatment_z_k)/n()) %>%
    ungroup() %>%
    dplyr::select(cluster_k, trt, Y1_sum, Y2_sum) %>%
    mutate(z1_y1 = trt*Y1_sum,
           z1_y2 = trt*Y2_sum,
           z0_y1 = (1 - trt)*Y1_sum,
           z0_y2 = (1 - trt)*Y2_sum)

  # Estimating betas based on GEE solution
  beta0_cal <- c(log((sum(sumDat$z0_y1))/(m_sim*K2 - sum(sumDat$z0_y1))),
                 log((sum(sumDat$z0_y2))/(m_sim*K2 - sum(sumDat$z0_y2))))
  beta1_intercept_cal <- beta0_cal[1]
  beta2_intercept_cal <- beta0_cal[2]

  beta1_vec_cal <- c(log(sum(sumDat$z1_y1)/(m_sim*K1 - sum(sumDat$z1_y1))) -
                       log(sum(sumDat$z0_y1)/(m_sim*K2 - sum(sumDat$z0_y1))),
                     log(sum(sumDat$z1_y2)/(m_sim*K1 - sum(sumDat$z1_y2))) -
                       log(sum(sumDat$z0_y2)/(m_sim*K2 - sum(sumDat$z0_y2))))
  beta1_cal <- beta1_vec_cal[1]
  beta2_cal <- beta1_vec_cal[2]

  # Thetas and zeta values
  theta0 <- c(exp(beta0_cal[1])/(exp(beta0_cal[1]) + 1),
              exp(beta0_cal[2])/(exp(beta0_cal[2]) + 1))
  theta1 <- c(exp(beta0_cal[1] + beta1_vec_cal[1])/
                (exp(beta0_cal[1] + beta1_vec_cal[1]) + 1),
              exp(beta0_cal[2] + beta1_vec_cal[2])/
                (exp(beta0_cal[2] + beta1_vec_cal[2]) + 1))

  zeta_0 <- c(theta0[1]*(1-theta0[1]), theta0[2]*(1-theta0[2]))
  zeta_1 <- c(theta1[1]*(1-theta1[1]), theta1[2]*(1-theta1[2]))

  mu1_C <- exp(beta1_intercept_cal)/(exp(beta1_intercept_cal) + 1)
  mu2_C <- exp(beta2_intercept_cal)/(exp(beta2_intercept_cal) + 1)

  mu1_T <- exp(beta1_intercept_cal + beta1_cal)/(exp(beta1_intercept_cal +
                                                       beta1_cal) + 1)
  mu2_T <- exp(beta2_intercept_cal + beta2_cal)/(exp(beta2_intercept_cal +
                                                       beta2_cal) + 1)

  # Calculate various matrices, uses outcome specific ICCs from geepack
  # Matches D. Li Code
  Omega_1 <- (m_sim + m_sim*(m_sim - 1)*rho01_geepack)*
    matrix(data = c((1 - z_bar)*zeta_0[1] + z_bar*zeta_1[1],
                                                       z_bar*zeta_1[1],
                                                       z_bar*zeta_1[1],
                                                       z_bar*zeta_1[1]),
                                              nrow = 2, ncol = 2)
  # Variance of S_2(Beta_2) (outcome 2)
  Omega_2 <- (m_sim + m_sim*(m_sim - 1)*rho02_geepack)*
    matrix(data = c((1 - z_bar)*zeta_0[2] + z_bar*zeta_1[2],
                                                       z_bar*zeta_1[2],
                                                       z_bar*zeta_1[2],
                                                       z_bar*zeta_1[2]),
                                              nrow = 2, ncol = 2)

  # Gamma matrices for each outcome (matches D. Li code)
  Gamma_1 <- m_sim*matrix(c((1 - z_bar)*zeta_0[1] + z_bar*zeta_1[1],
                        z_bar*zeta_1[1],
                        z_bar*zeta_1[1],
                        z_bar*zeta_1[1]),
                      nrow = 2, ncol = 2)

  Gamma_2 <- m_sim*matrix(c((1 - z_bar)*zeta_0[2] + z_bar*zeta_1[2],
                        z_bar*zeta_1[2],
                        z_bar*zeta_1[2],
                        z_bar*zeta_1[2]),
                      nrow = 2, ncol = 2)

  # Covariance of S_1(Beta_1) and S_2(Beta_2) Matches D. Li Code
  Sigma_1 <- solve(Gamma_1) %*% Omega_1 %*% solve(Gamma_1)
  Sigma_2 <- solve(Gamma_2) %*% Omega_2 %*% solve(Gamma_2)

  # (2,2)th element of Sigma_K, matches D. Li Code
  sigma2sq_1 <- ((1 + (m_sim - 1)*rho01_geepack)/(m_sim*(1-z_bar)*zeta_0[1])) +
    ((1 + (m_sim - 1)*rho01_geepack)/(m_sim*z_bar*zeta_1[1]))

  sigma2sq_2 <- ((1 + (m_sim - 1)*rho02_geepack)/(m_sim*(1-z_bar)*zeta_0[2])) +
    ((1 + (m_sim - 1)*rho02_geepack)/(m_sim*z_bar*zeta_1[2]))


  ## Creating output specification table ---------------------------------------
  sim_data_stats <- tibble(Parameter = c("K Total", "K Treatment", "m",
                                         "beta1", "beta1 intercept",
                                         "beta2", "beta2 intercept",
                                         "rho01", "rho02",
                                         "rho1", "rho2",
                                         "varY1", "varY2",
                                         "r"),
                           `Estimated Value` = c(K_Total_sim,
                                                 K_Treatment_sim,
                                                 m_sim,
                                                 NA,
                                                 NA,
                                                 NA,
                                                 NA,
                                                 NA,
                                                 NA,
                                                 NA,
                                                 NA,
                                                 round(varY1_sim, 4),
                                                 round(varY2_sim, 4),
                                                 r_sim),
                           `Estimated Value (mmm package)` = c(NA,
                                                               NA,
                                                               NA,
                                                               round(beta1_mmm, 4),
                                                               round(beta1_intercept_mmm, 4),
                                                               round(beta2_mmm, 4),
                                                               round(beta2_intercept_mmm, 4),
                                                               round(rho01_mmm, 4),
                                                               round(rho02_mmm, 4),
                                                               round(rho1_mmm, 4),
                                                               round(rho2_mmm, 4),
                                                               NA,
                                                               NA,
                                                               NA),
                           `Estimated Value (geepack package)` = c(NA,
                                                                   NA,
                                                                   NA,
                                                                   round(beta1_geepack, 4),
                                                                   round(beta1_intercept_geepack, 4),
                                                                   round(beta2_geepack, 4),
                                                                   round(beta2_intercept_geepack, 4),
                                                                   round(rho01_geepack, 4),
                                                                   round(rho02_geepack, 4),
                                                                   NA,
                                                                   NA,
                                                                   NA,
                                                                   NA,
                                                                   NA),
                           `Estimated Value (direct calculation)` = c(NA,
                                                                      NA,
                                                                      NA,
                                                                      round(beta1_cal, 4),
                                                                      round(beta1_intercept_cal, 4),
                                                                      round(beta2_geepack, 4),
                                                                      round(beta2_intercept_cal, 4),
                                                                      NA,
                                                                      NA,
                                                                      NA,
                                                                      NA,
                                                                      NA,
                                                                      NA,
                                                                      NA)
  )

  return(sim_data_stats)
} # End GEE.est()

create_all_bin_sim_dats <- function(n = 100,
                                    K_input,
                                    m_input,
                                    beta1_input,
                                    beta2_input,
                                    rho01_input,
                                    rho02_input,
                                    rho1_input,
                                    rho2_input,
                                    varY1_input,
                                    varY2_input,
                                    r_input = 1){

  simDataList <- vector(mode = 'list', length = n)
  simStatList <- vector(mode = 'list', length = n)

  for(i in 1:n){
    mySimData <- gen_crt_coprimary_data_bin(K = K_input,
                                            m = m_input,
                                            beta1 = beta1_input,
                                            beta2 = beta2_input,
                                            rho01 = rho01_input,
                                            rho02 = rho02_input,
                                            rho1 = rho1_input,
                                            rho2 = rho2_input,
                                            varY1 = varY1_input,
                                            varY2 = varY2_input,
                                            r = r_input)

    simStats <- GEE.est(myData = mySimData) %>%
      mutate(`True Value` = c(K_input + r_input*K_input,
                              K_input, m_input,
                              beta1_input, 0,
                              beta2_input, 0,
                              rho01_input, rho02_input,
                              rho1_input, rho2_input,
                              varY1_input, varY2_input,
                              r_input)) %>%
      relocate(`Parameter`, `True Value`)

    simDataList[[i]] <- mySimData
    simStatList[[i]] <- simStats
  }

  return(list(simDataList, simStatList))
} # End create_all_bin_sim_dats()
