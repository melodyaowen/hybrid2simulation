# calVectorP() -----------------------------------------------------
# Calculates p1 and p2, where p1 is the proportion of individuals in a single cluster with Y1 = 1 and p2 is the proportion of individuals in a single cluster with Y2 = 1
# calTotalVarianceBinary <- function(VarY1, VarY2,
#                                    rho01, rho02){
#
#   # Calculating within-cluster component of variance for Y1 and Y2
#   sigma2W_1 <- VarY1 - VarY1*rho01
#   sigma2W_2 <- VarY2 - VarY2*rho02
#
#   # Calculating p1 and p2
#   p_1 <- c((1 + sqrt(1 - 4*sigma2W_1))/2, (1 - sqrt(1 - 4*sigma2W_1))/2)
#   p_2 <- c((1 + sqrt(1 - 4*sigma2W_2))/2, (1 - sqrt(1 - 4*sigma2W_2))/2)
#
#   return(c(p_1, p_2))
# }
#
# calTotalVarianceBinary(VarY1 = 0.23, VarY2 = 0.25,
#                        rho01 = 0.025, rho02 = 0.025)


# calTotalVarianceBinary() -----------------------------------------------------
# Calculates the total variance from the within-cluster variance and ICC
# rho0 is the outcome specific ICC
# p is the proportion of individuals in a single cluster with Y = 1 (binary)
# calTotalVarianceBinary <- function(rho0, p){
#   sigma2W <- p*(1-p) # within-cluster component of variance
#   sigma2B <- (rho0*sigma2W)/(1-rho0) # between-cluster component of variance
#   totalVar <- sigma2W + sigma2B
#   return(totalVar)
# }


simDataBinLong <- simDataBin %>%
  gather(c(Y1, Y2), key = "Y", value = "Value")

m1 <- MCMCglmm(cbind(Y1, Y2) ~ treatment_z_k, family = rep("categorical", 2),
               data = as.data.frame(simDataBin))

# Using 'MultiOrd' package
ordPmat1 <- matrix( c(0.15, 0.70,
                      0.55, 0.10,
                      0.25, 0.10), 4, 3, byrow = TRUE)
## End(Not run)
cmat1 <- matrix(c(1,0.2,0.2,
                  0.2,1,0.2,
                  0.2,0.2,1), 3, 3, byrow = TRUE)

p <- find.binary.prob(ordPmat1)
finalCorr <- simBinCorr(ordPmat1, cmat1, no.rows = 100000)
y <- generate.binary(1000, p$p, finalCorr$del.next)



# Script with data generation functions for two binary endpoints

# logistic_fun() ---------------------------------------------------------------
logistic_fun <- function(p){
  return(log(p / (1-p)))
}

# logistic_inverse_fun() -------------------------------------------------------
logistic_inverse_fun <- function(log_odds){
  return(exp(log_odds) / (1 + exp(log_odds)))
}

# calTotalVarianceBinary() -----------------------------------------------------
# Calculates the total variance from the within-cluster variance and ICC
# rho0 is the outcome specific ICC
# p is the proportion of individuals in a single cluster with Y = 1 (binary)
# calTotalVarianceBinary <- function(rho0, p){
#   sigma2W <- p*(1-p) # within-cluster component of variance
#   sigma2B <- (rho0*sigma2W)/(1-rho0) # between-cluster component of variance
#   totalVar <- sigma2W + sigma2B
#   return(totalVar)
# }

# calVectorP() -----------------------------------------------------
# Calculates p1 and p2, where p1 is the proportion of individuals in a single cluster with Y1 = 1 and p2 is the proportion of individuals in a single cluster with Y2 = 1
# calTotalVarianceBinary <- function(VarY1, VarY2,
#                                    rho01, rho02){
#
#   # Calculating within-cluster component of variance for Y1 and Y2
#   sigma2W_1 <- VarY1 - VarY1*rho01
#   sigma2W_2 <- VarY2 - VarY2*rho02
#
#   # Calculating p1 and p2
#   p_1 <- c((1 + sqrt(1 - 4*sigma2W_1))/2, (1 - sqrt(1 - 4*sigma2W_1))/2)
#   p_2 <- c((1 + sqrt(1 - 4*sigma2W_2))/2, (1 - sqrt(1 - 4*sigma2W_2))/2)
#
#   return(c(p_1, p_2))
# }
#
# calTotalVarianceBinary(VarY1 = 0.23, VarY2 = 0.25,
#                        rho01 = 0.025, rho02 = 0.025)

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
}

#beta1_input <- 0.8
#beta2_input <- 0.1

simDataBin <- gen_crt_coprimary_data_bin(K = 10, # Number of clusters in treatment arm
                                         m = 100, # Number of individuals in each cluster
                                         beta1 = 0.8, # Int effect for Y1
                                         beta2 = 0.1, # Int effect for Y2
                                         rho01 = 0.07, # ICC for Y1
                                         rho02 = 0.03, # ICC for Y2
                                         rho1 = 0.01, # Corr(Y1, Y2) for two diff individuals
                                         rho2 = 0.05, # Corr(Y1, Y2) for same individual
                                         varY1 = 0.5, # Var(Y1)
                                         varY2 = 0.2, # Var(Y2)
                                         r = 1
)


fit <- mmm(formula = cbind(Y1, Y2) ~ treatment_z_k,
           id = simDataBin$cluster_k, data = simDataBin,
           family = binomial(link = logit), corStruct = "exchangeable")


GEE.est <- function(myData){

  # Cluster size of inputted data frame
  m <- nrow(dplyr::filter(myData, cluster_k == 1))

  # Number of clusters in treatment group
  K1 <- dplyr::select(myData, cluster_k, treatment_z_k) %>%
    unique() %>%
    filter(treatment_z_k == 1) %>%
    nrow()

  # Number of clusters in control group
  K2 <- dplyr::select(myData, cluster_k, treatment_z_k) %>%
    unique() %>%
    filter(treatment_z_k == 0) %>%
    nrow()

  # Total K
  K_total <- K1 + K2

  # Treatment groups
  groups <- myData %>%
    dplyr::select(cluster_k, treatment_z_k) %>%
    unique()

  # Estimate of outcome specific variance
  varY1_est <- var(myData$Y1)
  varY2_est <- var(myData$Y2)

  # Probability of a cluster receiving the experimental intervention
  z_bar <- K1/(K1 + K2)

  # Calculating outcome specific ICC
  m1 <- geepack::geeglm(formula = Y1 ~ treatment_z_k, # creating a geeglm model
                        id = factor(cluster_k),
                        family = binomial(link = "logit"),
                        data = myData, corstr = "exchangeable")
  m1_summary <- summary(m1) # summary output for m1
  rho01_est <- m1_summary[12]$corr[[1]]

  m2 <- geepack::geeglm(formula = Y2 ~ treatment_z_k, # creating a geeglm model
                        id = factor(cluster_k),
                        family = binomial(link = "logit"),
                        data = myData, corstr = "exchangeable")
  m2_summary <- summary(m2) # summary output for m1
  rho02_est <- m2_summary[12]$corr[[1]]

  # Fitting generalized multivariate linear mixed model
  fit <- mmm(formula = cbind(Y1, Y2) ~ treatment_z_k,
             id = simDataBin$cluster_k, data = simDataBin,
             family = binomial(link = logit), corStruct = "exchangeable")

  beta1_est <- round(fit$coefficients[["Y1.treatment_z_k"]], 4)
  beta2_est <- round(fit$coefficients[["Y2.treatment_z_k"]], 4)

  beta1_intercept_est <- round(fit$coefficients[["Y1.Intercept"]], 4)
  beta2_intercept_est <- round(fit$coefficients[["Y2.Intercept"]], 4)

  wCorr_est <- fit$working.correlation # dim (m*2) x (m*2) # R in D.Li Paper
  C11_est <- wCorr_est[1:m, 1:m]
  C12_est <- wCorr_est[1:m, (m+1):(2*m)]
  C21_est <- wCorr_est[(m+1):(2*m), 1:m]
  C22_est <- wCorr_est[(m+1):(2*m), (m+1):(2*m)]

  identical(C12_est, C21_est)
  V

} # End GEE.est



GEE.est <- function(myData){

  # Cluster size of inputted data frame
  m <- nrow(dplyr::filter(myData, cluster_k == 1))

  # Number of clusters in treatment group
  K1 <- dplyr::select(myData, cluster_k, treatment_z_k) %>%
    unique() %>%
    filter(treatment_z_k == 1) %>%
    nrow()

  # Number of clusters in control group
  K2 <- dplyr::select(myData, cluster_k, treatment_z_k) %>%
    unique() %>%
    filter(treatment_z_k == 0) %>%
    nrow()

  # Total K
  K_total <- K1 + K2

  # Treatment groups
  groups <- myData %>%
    dplyr::select(cluster_k, treatment_z_k) %>%
    unique()

  # Estimate of outcome specific variance
  varY1_est <- var(myData$Y1)
  varY2_est <- var(myData$Y2)

  # Probability of a cluster receiving the experimental intervention
  z_bar <- K1/(K1 + K2)

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
  beta0_est <- c(log((sum(sumDat$z0_y1))/(m*K2 - sum(sumDat$z0_y1))),
                 log((sum(sumDat$z0_y2))/(m*K2 - sum(sumDat$z0_y2))))

  beta1_est <- c(log(sum(sumDat$z1_y1)/(m*K1 - sum(sumDat$z1_y1))) -
                   log(sum(sumDat$z0_y1)/(m*K2 - sum(sumDat$z0_y1))),
                 log(sum(sumDat$z1_y2)/(m*K1 - sum(sumDat$z1_y2))) -
                   log(sum(sumDat$z0_y2)/(m*K2 - sum(sumDat$z0_y2))))

  fit <- mmm(formula = cbind(Y1, Y2) ~ treatment_z_k,
             id = simDataBin$cluster_k, data = simDataBin,
             family = binomial(link = logit), corStruct = "exchangeable")

  wCorr_est <- fit$working.correlation # dim (m*2) x (m*2) # R in D.Li Paper
  C11_est <- wCorr_est[1:m, 1:m]
  C12_est <- wCorr_est[1:m, (m+1):(2*m)]
  C21_est <- wCorr_est[(m+1):(2*m), 1:m]
  C22_est <- wCorr_est[(m+1):(2*m), (m+1):(2*m)]



  # Theta's also are notated as mu in the paper
  theta0 <- c(exp(beta0_est[1])/(exp(beta0_est[1]) + 1),
              exp(beta0_est[2])/(exp(beta0_est[2]) + 1))
  theta1 <- c(exp(beta0_est[1] + beta1_est[1])/(exp(beta0_est[1] + beta1_est[1]) + 1),
              exp(beta0_est[2] + beta1_est[2])/(exp(beta0_est[2] + beta1_est[2]) + 1))

  zeta_0 <- c(theta0[1]*(1-theta0[1]), theta0[2]*(1-theta0[2]))
  zeta_1 <- c(theta1[1]*(1-theta1[1]), theta1[2]*(1-theta1[2]))

  m1 <- geepack::geeglm(formula = Y1 ~ treatment_z_k, # creating a geeglm model
                        id = factor(cluster_k),
                        family = binomial(link = "logit"),
                        data = myData, corstr = "exchangeable")
  m1_summary <- summary(m1) # summary output for m1
  rho01_est <- m1_summary[12]$corr[[1]]

  m2 <- geepack::geeglm(formula = Y2 ~ treatment_z_k, # creating a geeglm model
                        id = factor(cluster_k),
                        family = binomial(link = "logit"),
                        data = myData, corstr = "exchangeable")
  m2_summary <- summary(m2) # summary output for m1
  rho02_est <- m2_summary[12]$corr[[1]]

  # Variance of S_1(Beta_1) (outcome 1)
  Omega_1 <- (m + m*(m - 1)*rho01_est)*matrix(data = c((1 - z_bar)*zeta_0[1] + z_bar*zeta_1[1],
                                                       z_bar*zeta_1[1],
                                                       z_bar*zeta_1[1],
                                                       z_bar*zeta_1[1]),
                                              nrow = 2, ncol = 2)

  # Variance of S_2(Beta_2) (outcome 2)
  Omega_2 <- (m + m*(m - 1)*rho02_est)*matrix(data = c((1 - z_bar)*zeta_0[2] + z_bar*zeta_1[2],
                                                       z_bar*zeta_1[2],
                                                       z_bar*zeta_1[2],
                                                       z_bar*zeta_1[2]),
                                              nrow = 2, ncol = 2)

  # Gamma matrices for each outcome
  Gamma_1 <- m*matrix(c((1 - z_bar)*zeta_0[1] + z_bar*zeta_1[1],
                        z_bar*zeta_1[1],
                        z_bar*zeta_1[1],
                        z_bar*zeta_1[1]),
                      nrow = 2, ncol = 2)

  Gamma_2 <- m*matrix(c((1 - z_bar)*zeta_0[2] + z_bar*zeta_1[2],
                        z_bar*zeta_1[2],
                        z_bar*zeta_1[2],
                        z_bar*zeta_1[2]),
                      nrow = 2, ncol = 2)

  # Covariance of S_1(Beta_1) and S_2(Beta_2)
  Sigma_1 <- solve(Gamma_1) %*% Omega_1 %*% solve(Gamma_1)
  Sigma_2 <- solve(Gamma_2) %*% Omega_2 %*% solve(Gamma_2)

  # Calculating Sigma_12 matrix
  mu1_k_0 <- theta0[1] # First outcome, k cluster, control group
  mu1_k_1 <- theta0[2] # First outcome, k cluster, treatment group
  mu2_k_0 <- theta1[1] # Second outcome, k cluster, control group
  mu2_k_1 <- theta1[2] # Second outcome, k cluster, treatment group

  # Diagonal V matrix per trt group, per outcome
  V1_k_0 <- mu1_k_0*(1-mu1_k_0)*diag(1, nrow = 1, ncol = 1) # outcome 1, trt 0
  V1_k_1 <- mu1_k_1*(1-mu1_k_1)*diag(1, nrow = 1, ncol = 1) # outcome 1, trt 1
  V2_k_0 <- mu2_k_0*(1-mu2_k_0)*diag(1, nrow = 1, ncol = 1) # outcome 1, trt 0
  V2_k_1 <- mu2_k_1*(1-mu2_k_1)*diag(1, nrow = 1, ncol = 1) # outcome 1, trt 1

  as.matrix(rbind(cbind(V1_k_1, matrix(0, 1, 1)),
                  cbind(matrix(0, 1, 1), V2_k_1))) %*% c(1, z_bar)

  as.matrix(rbind(cbind(V1_k_1, matrix(0, 1, 1)),
                  cbind(matrix(0, 1, 1), V2_k_1))) %*% matrix(c(1, z_bar),
                                                              nrow = 2,
                                                              ncol = 1)

  # Generate matrix C, a positive definite matrix of correlations for the
  # within-cluster errors e_k
  # off-diagonals can be 0 (this assumes observations from same cluster are
  # conditionally independent given the random effects in univariate and
  # multivariate models and across outcome measures in the multivariate setting)
  C11 <- (1 - rho01_est)*diag(1) + rho01_est*matrix(1, nrow = 1, ncol = 1)
  C22 <- (1 - rho02_est)*diag(1) + rho02_est*matrix(1, nrow = 1, ncol = 1)
  #C12 <- (rho2 - rho1)*diag(1) + rho1*matrix(1, nrow = 1, ncol = 1)
  #C <- rbind(cbind(C11, C12), cbind(C12, C22))

  sum_vec <- c()
  for(k in 1:nrow(groups)){ # Loop through each cluster,
    # Specify V_k based on trt group
    if(groups$treatment_z_k[k] == 1){ # Treatment for this cluster is 1
      V_temp <- as.matrix(rbind(cbind(V1_k_1, matrix(0, 1, 1)),
                                cbind(matrix(0, 1, 1), V2_k_1)))
      D_temp <- as.matrix(rbind(cbind(V1_k_1, matrix(0, 1, 1)),
                                cbind(matrix(0, 1, 1), V2_k_1))) %*% c(1, z_bar)
    } else if(groups$trt_group[k] == 0){ # Trt for this cluster is 0
      V_temp <- as.matrix(rbind(cbind(V1_k_0, matrix(0, 1, 1)),
                                cbind(matrix(0, 1, 1), V2_k_0)))
      D_temp <- as.matrix(rbind(cbind(V1_k_0, matrix(0, 1, 1)),
                                cbind(matrix(0, 1, 1), V2_k_0))) %*% c(1, z_bar)
    }

    #R_temp <- sqrt(V_temp) %*% C %*% sqrt(V_temp)

    D_temp[1] %*% c(beta0_est[1], beta1_est[1]) %*% solve(V_temp) * v12_kk

    sum_vec[k] <- 1

    if(is.positive.definite(R_temp) == FALSE){
      stop("R_k matrix is not positive definite. Check input values for correlations and variances.")
    }

  }



  V1 <- rbind()

  D_1k # outcome 1, per cluster
  D_2k # outcome 2, per cluster


  beta1_est_vector <- c(beta0_est[1], beta1_est[1])
  beta2_est_vector <- c(beta0_est[2], beta1_est[2])
  V_1k # V matrix outcome 1, per cluster
  V_2k # V matrix outcome 2, per cluster

  Sigma_12_Sum <- myData %>%
    group_by(cluster_k)



  cov((K_total^(1/2))*(beta1_est[1] - beta1_input), (K_total^(1/2))*(beta1_est[2] - beta2_input))



  # Same as (2,2) element in Sigma_q matrices above
  sigma2_21 <- (1 + (m-1)*rho01_est)/(m*(1-z_bar)*zeta_0[1]) +
    (1 + (m-1)*rho01_est)/(m*z_bar*zeta_1[1])
  sigma2_22 <- (1 + (m-1)*rho02_est)/(m*(1-z_bar)*zeta_0[2]) +
    (1 + (m-1)*rho02_est)/(m*z_bar*zeta_1[2])

  cor(myData$Y1, myData$Y2)


  var(myData$Y1)
  var(myData$Y2)

}

simDataBinLong <- simDataBin %>%
  gather(c(Y1, Y2), key = "Y", value = "Value")

m1 <- MCMCglmm(cbind(Y1, Y2) ~ treatment_z_k, family = rep("categorical", 2),
               data = as.data.frame(simDataBin))

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
#prop_Y1 <- 0.66
#prop_Y2 <- 0.6

exp(0 + beta1)/(1 + exp(0 + beta1))
exp(0 + beta2)/(1 + exp(0 + beta2))

myData <- simDataBin


# Using 'MultiOrd' package
ordPmat1 <- matrix( c(0.15, 0.70,
                      0.55, 0.10,
                      0.25, 0.10), 4, 3, byrow = TRUE)
## End(Not run)
cmat1 <- matrix(c(1,0.2,0.2,
                  0.2,1,0.2,
                  0.2,0.2,1),3,3,byrow=TRUE)

p <- find.binary.prob(ordPmat1)
finalCorr <- simBinCorr(ordPmat1, cmat1, no.rows = 100000)
y <- generate.binary(1000, p$p, finalCorr$del.next)


