# Script with data generation functions for two continuous endpoints, MLMM

# gen_crt_coprimary_data_cont() ------------------------------------------------
# Function to generate data for Normal Parallel CRT Data for Equal Allocation
# and equal cluster sizes
# Assumes Exchangeable Correlation Structure
# Recall notation: j = 1,...,m is index of cluster size
#                  k = 1,...,2K is index of cluster (K is number in trt arm)
#                  q = 1, 2 is index of outcome (Q = 2)
#                  i = 1, 2 is treatment group index
gen_crt_coprimary_data_cont <- function(K, # Number of clusters in treatment arm
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

  # Beta vectors and other important parameters
  betas <- c(beta1, beta2)
  deltas <- c(0, 0)
  Q <- 2

  # Total Number of Clusters given a Binary Treatment by Cluster
  K1 <- ceiling(K)
  K2 <- ceiling(r*K)
  K_total <- K1 + K2
  study_n <- K_total*m

  # Define sigma matrix
  sigmaks.sq <- diag(calCovbetas_eq(vars = c(varY1, varY2),
                                    rho01 = matrix(c(rho01, rho1,
                                                     rho1, rho02),
                                                   2, 2),
                                    rho2 = matrix(c(1, rho2,
                                                    rho2, 1),
                                                  2, 2),
                                    sigmaz.square = sigmaz.square,
                                    m = m,
                                    Q = 2 # Number of outcomes
  )
  )

  # Define mean vector
  meanVector <- sqrt(K_total)*(betas - deltas)/sqrt(sigmaks.sq)

  # May not need
  # # wCor is correlation between test statistics
  # wCor <- calCorWks(vars = c(varY1, varY2),
  #                   rho01 = matrix(c(rho01, rho1,
  #                                    rho1, rho02),
  #                                  2, 2),
  #                   rho2 = matrix(c(1, rho2,
  #                                   rho2, 1),
  #                                 2, 2),
  #                   sigmaz.square,
  #                   m = m,
  #                   Q = 2 # Number of outcomes
  #                   )

  # Covariance matrix Sigma_phi for Y_i
  Sigma_phi <- constrRiP(rho01 = matrix(c(rho01, rho1,
                                          rho1, rho02),
                                        2, 2),
                         Q = 2, vars = c(varY1, varY2))

  # Covariance matrix Sigma_E for Y_i
  Sigma_E <- constrRiE(rho01 = matrix(c(rho01, rho1,
                                        rho1, rho02),
                                      2, 2),
                       rho2 = matrix(c(1, rho2,
                                       rho2, 1),
                                     2, 2),
                       Q = 2, vars = c(varY1, varY2))

  # Check if matrices are positive definite, if not then throw error
  if(is.positive.definite(Sigma_phi) == FALSE){
    stop("Sigma_phi matrix is not positive definite. Check input values for correlations and variances.")
  }
  if(is.positive.definite(Sigma_E) == FALSE){
    stop("Sigma_E matrix is not positive definite. Check input values for correlations and variances.")
  }

  # Random assignment of clusters to treatment
  groups <- tibble(cluster_id = 1:K_total,
                   trt_group = sample(c(rep(1, K1),
                                        rep(0, K2)),
                                      replace = FALSE))

  # Generate vector of random intercepts for cluster k across Q = 2 endpoints
  # phi_k <- mvrnorm(n = K_total, mu = rep(0, 2), Sigma = Sigma_phi) %>%
  #   as.data.frame() %>%
  #   dplyr::select(phi_k1 = V1, phi_k2 = V2) %>%
  #   mutate(cluster_id = 1:K_total)

  # Using a truncated normal (above 0)
  phi_k <- rtmvnorm(n = K_total, mean = rep(0, 2), sigma = Sigma_phi,
                    lower = rep(0, length = 2), upper = rep(Inf, length = 2)) %>%
    as.data.frame() %>%
    dplyr::select(phi_k1 = V1, phi_k2 = V2) %>%
    mutate(cluster_id = 1:K_total)

  # Generate vector of random intercepts for person j across Q = 2 endpoints
  # e_kj <- mvrnorm(n = m*K_total, mu = rep(0, 2), Sigma = Sigma_E) %>%
  #   as.data.frame() %>%
  #   dplyr::select(e_kj1 = V1, e_kj2 = V2) %>%
  #   mutate(subject_id = 1:study_n,
  #          cluster_id = sort(rep(1:K_total, m)))

  # Using a truncated normal (above 0)
  e_kj <- rtmvnorm(n = m*K_total, mean = rep(0, 2), sigma = Sigma_E,
                   lower = rep(0, length = 2), upper = rep(Inf, length = 2)) %>%
    as.data.frame() %>%
    dplyr::select(e_kj1 = V1, e_kj2 = V2) %>%
    mutate(subject_id = 1:study_n,
           cluster_id = sort(rep(1:K_total, m)))

  # Generate dataset based on specifications
  simDat <- tibble(subject_id = 1:study_n,
                   cluster_id = sort(rep(1:K_total, m))) %>%
    left_join(., groups, by = "cluster_id") %>%
    left_join(., phi_k, by = "cluster_id") %>%
    left_join(., e_kj, by = c("cluster_id", "subject_id")) %>%
    # Now create outcome columns based on inputs and calculated intercepts
    mutate(Y1 = beta1*trt_group + phi_k1 + e_kj1,
           Y2 = beta2*trt_group + phi_k2 + e_kj2,
           Yc = Y1 + Y2) %>%
    dplyr::rename(subject_j = subject_id,
                  cluster_k = cluster_id,
                  treatment_z_k = trt_group)

  # Return data
  return(simDat)
}

# constrRiE() ------------------------------------------------------------------
# Function to construct covariance matrix Sigma_E for Y_i
constrRiE <- function(rho01, rho2, Q, vars){
  rho0Q <- diag(rho01)
  SigmaE_Matrix <- diag((1-rho0Q)*vars)
  for(row in 1:Q){
    for(col in 1:Q){
      if(row != col){
        SigmaE_Matrix[row,col] <- sqrt(vars[row])*sqrt(vars[col])*
          (rho2[row,col]-rho01[row,col])
      }
    }
  }
  return(SigmaE_Matrix)
}

# constrRiP() ------------------------------------------------------------------
# Function to construct covariance matrix Sigma_phi for Y_i
constrRiP <- function(rho01, Q, vars){
  rho0Q <- diag(rho01)
  SigmaP_Matrix <- diag(rho0Q*vars)
  for(row in 1:Q){
    for(col in 1:Q){
      if(row != col){
        SigmaP_Matrix[row,col] <- sqrt(vars[row])*sqrt(vars[col])*rho01[row,col]
      }
    }
  }
  return(SigmaP_Matrix)
}

# calCovbetas_eq() -------------------------------------------------------------
# Function to calculate covariance between betas, equal cluster sizes
calCovbetas_eq <- function(vars, rho01, rho2, sigmaz.square, m, Q){
  rho0Q <- diag(rho01)
  sigmaQ.square <-(1 + (m-1)*rho0Q)*vars/(m*sigmaz.square)
  covMatrix <- diag(sigmaQ.square)
  for(row in 1:Q){
    for(col in 1:Q){
      if(row != col){
        covMatrix[row,col] <- sqrt(vars[row])*sqrt(vars[col])*
          (rho2[row,col]+(m-1)*rho01[row,col])/(m*sigmaz.square)
      }
    }
  }
  return(covMatrix)
}

# calCorWks() ------------------------------------------------------------------
# Calculates correlation between test statistics
calCorWks <- function(vars, rho01, rho2, sigmaz.square, m, Q){
  top <- calCovbetas_eq(vars, rho01, rho2, sigmaz.square, m, Q)
  wCor <- diag(Q)
  for(row in 1:Q){
    for(col in 1:Q){
      if(row != col){
        wCor[row,col] <- top[row,col]/sqrt(top[row,row]*top[col,col])
      }
    }
  }
  return(wCor)
}

# MLMM.estim() -------------------------------------------------------------------
# Fit an MLMM to estimate parameters
MLMM.estim <- function(myData){

  # Combined Outcome Assessment
  fmC <- nlme::lme(Yc ~ treatment_z_k, random = ~ 1|cluster_k,
                   data = myData)
  betaC <- as.numeric(fmC$coefficients$fixed) # 2nd is treatment, 1st is int

  # Get Var(BetaC)
  vcov_fmC <- vcov(fmC)
  varBetaC_sim <- vcov_fmC["treatment_z_k", "treatment_z_k"]

  # Calculate outcome specific ICC for Yc
  sigma2_phi_C <- as.numeric(VarCorr(fmC)[1,1])
  sigma2_e_C <- as.numeric(VarCorr(fmC)[2,1])
  rho0C_sim <- sigma2_phi_C/(sigma2_phi_C + sigma2_e_C)

  # Calculate Var(Yc)
  varYc_model_sim <- sigma2_phi_C + sigma2_e_C

  # Fit MLMM, must transform the data first
  myDataStacked <- myData %>%
    tidyr::pivot_longer(cols = c("Y1", "Y2"),
                        names_to = "Outcome_ID", values_to = "Outcome_Value") %>%
    dplyr::select(subject_j, cluster_k, treatment_z_k,
                  Outcome_ID, Outcome_Value, phi_k1, phi_k2, e_kj1, e_kj2) %>%
    dplyr::arrange(cluster_k, subject_j)

  myMLMM <- nlme::lme(Outcome_Value ~ Outcome_ID * treatment_z_k,
                      random = ~ Outcome_ID | cluster_k/subject_j,
                      data = myDataStacked,
                      method = "REML")

  # Get the summary of the model
  myMLMM_summary <- summary(myMLMM)

  # Extract fixed effects coefficients
  fixed_effects <- myMLMM_summary$tTable

  # Extract Beta1, Beta2, and their intercepts
  intercept_Y1 <- fixed_effects["(Intercept)", "Value"]
  intercept_Y2 <- intercept_Y1 + fixed_effects["Outcome_IDY2", "Value"]
  beta1 <- fixed_effects["treatment_z_k", "Value"]
  beta2 <- beta1 + fixed_effects["Outcome_IDY2:treatment_z_k", "Value"]

  # Extract the variances of the betas for getting the test statistics later
  vcov_matrix <- vcov(myMLMM) # variance-covariance matrix of the fixed effects

  var_beta1 <- vcov_matrix["treatment_z_k", "treatment_z_k"]
  var_beta2 <- var_beta1 + vcov_matrix["Outcome_IDY2:treatment_z_k",
                                       "Outcome_IDY2:treatment_z_k"] +
    2 * vcov_matrix["treatment_z_k", "Outcome_IDY2:treatment_z_k"]

  # Get variance components for cluster and individual-level components
  var_components <- VarCorr(myMLMM) # Extract variance components (random efcts)

  cluster_var_Y1 <- as.numeric(var_components[2,1]) # Cluster Intercept, Var Column
  cluster_var_Y2 <- as.numeric(var_components[3,1]) # Cluster Outcome_IDY2, Var Column

  subject_var_Y1 <- as.numeric(var_components[5,1]) # Subject Intercept, Var Column
  subject_var_Y2 <- as.numeric(var_components[6,1]) # Subject OutcomeIDY2, Var Column
  residual_var <- as.numeric(var_components[7,1]) # Residual variance

  # Calculate total variances for Y1 and Y2
  total_var_Y1 <- cluster_var_Y1 + subject_var_Y1 + residual_var
  total_var_Y2 <- cluster_var_Y2 + subject_var_Y2 + residual_var

  # Calculate ICC for Y1 and Y2
  ICC_Y1 <- cluster_var_Y1 / total_var_Y1
  ICC_Y2 <- cluster_var_Y2 / total_var_Y2



  fm1 <- nlme::lme(Y1 ~ treatment_z_k, random = ~ 1|cluster_k, data = myData)
  fm2 <- nlme::lme(Y2 ~ treatment_z_k, random = ~ 1|cluster_k, data = myData)
  Q <- 2 # Number of outcomes
  zeta <- as.numeric(c(fm1$coefficients$fixed, fm2$coefficients$fixed))
  beta1 <- zeta[1:2]
  beta2 <- zeta[3:4]

  m <- as.numeric(table(myData$cluster_k))

  s2phi1 <- VarCorr(fm1)[1,1]
  s2phi2 <- VarCorr(fm2)[1,1]
  SigmaPhi <- diag(c(s2phi1, s2phi2))
  InvS2Phi <- solve(SigmaPhi)

  s2e1 <- VarCorr(fm1)[2,1]
  s2e2 <- VarCorr(fm2)[2,1]
  SigmaE <- diag(c(s2e1, s2e2))
  InvS2E <- solve(SigmaE)

  # Get Var(Beta1) and Var(Beta2)
  vcov_fm1 <- vcov(fm1)
  var_Beta1 <- vcov_fm1["treatment_z_k", "treatment_z_k"]
  vcov_fm2 <- vcov(fm2)
  var_Beta2 <- vcov_fm2["treatment_z_k", "treatment_z_k"]

  # Output formulas
  out1 <- all.vars(as.formula("Y1 ~ treatment_z_k"))[1]
  out2 <- all.vars(as.formula("Y2 ~ treatment_z_k"))[1]
  arm <- all.vars(as.formula("Y1 ~ treatment_z_k"))[2]

  Y <- as.matrix(myData[,c(out1, out2)])
  ID <- as.numeric(myData$cluster_k)
  n <- length(unique(ID)) # K1 + K2

  facz <- as.factor(as.data.frame(myData[,grep(arm, colnames(myData))])[[1]])
  z <- as.numeric(facz) - 1
  X <- as.matrix(cbind(1, z)) # design matrix

  # Calculate rho1 and rho2

  param <- list(theta = list(zeta = zeta, SigmaE = SigmaE, SigmaPhi = SigmaPhi),
                combined = list(betaC = betaC, varBetaC_sim = varBetaC_sim,
                                rho0C_sim = rho0C_sim,
                                varYc_model_sim = varYc_model_sim),
                varBetas = list(var_Beta1 = var_Beta1, var_Beta2 = var_Beta2))

  return(param)
}

# EM.estim() -------------------------------------------------------------------
# Uses EM estimation to estimate sample size calculation parameters
EM.estim <- function(myData, maxiter = 500, epsilon = 1e-4, verbose = FALSE){

  # Fit mixed model to initialize parameters
  fm1 <- nlme::lme(Y1 ~ treatment_z_k, random = ~ 1|cluster_k, data = myData)
  fm2 <- nlme::lme(Y2 ~ treatment_z_k, random = ~ 1|cluster_k, data = myData)
  K <- 2
  zeta <- as.numeric(c(fm1$coefficients$fixed, fm2$coefficients$fixed))
  beta1 <- zeta[1:2]
  beta2 <- zeta[3:4]

  m <- as.numeric(table(myData$cluster_k))

  s2phi1 <- VarCorr(fm1)[1,1]
  s2phi2 <- VarCorr(fm2)[1,1]
  SigmaPhi <- diag(c(s2phi1, s2phi2))
  InvS2Phi <- solve(SigmaPhi)

  s2e1 <- VarCorr(fm1)[2,1]
  s2e2 <- VarCorr(fm2)[2,1]
  SigmaE <- diag(c(s2e1, s2e2))
  InvS2E <- solve(SigmaE)

  out1 <- all.vars(as.formula("Y1 ~ treatment_z_k"))[1]
  out2 <- all.vars(as.formula("Y2 ~ treatment_z_k"))[1]
  arm <- all.vars(as.formula("Y1 ~ treatment_z_k"))[2]

  Y <- as.matrix(myData[,c(out1, out2)])
  ID <- as.numeric(myData$cluster_k)
  n <- length(unique(ID))

  facz <- as.factor(as.data.frame(myData[,grep(arm, colnames(myData))])[[1]])
  z <- as.numeric(facz) - 1
  X <- as.matrix(cbind(1, z)) # design matrix

  ESSphi1 <- matrix(0, n, K)
  ESSphi2 <- array(0,c(K, K, n))

  delta <- 2*epsilon
  max_modi <- 20
  converge <- 0

  # log likelihood
  loglik <- function(theta){
    beta1 <- theta[1:2]
    beta2 <- theta[3:4]
    sphi11 <- theta[5]
    sphi12 <- theta[6]
    sphi22 <- theta[7]
    se11 <- theta[8]
    se12 <- theta[9]
    se22 <- theta[10]
    SigmaPhi <- matrix(c(sphi11,sphi12,sphi12,sphi22),2,2)
    SigmaE <- matrix(c(se11,se12,se12,se22),2,2)
    InvS2Phi <- solve(SigmaPhi)
    InvS2E <- solve(SigmaE)

    temp <- 0
    for(j in unique(ID)){
      j_index <- match(j, unique(ID))
      Yj <- Y[ID == j,,drop=FALSE]
      Xj <- X[ID == j,,drop=FALSE]
      residj <- Yj - cbind(Xj%*%beta1, Xj%*%beta2)
      obs <- c(t(residj))
      tm1 <- (m[j_index] - 1)*log(det(SigmaE)) +
        log(det(SigmaE + m[j_index]*SigmaPhi))
      InvSS2 <- solve(SigmaE + m[j_index]*SigmaPhi) - InvS2E
      Invj <- kronecker(diag(nrow = m[j_index]), InvS2E) +
        kronecker(matrix(1, m[j_index], m[j_index]), InvSS2)/m[j_index]
      tm2 <- c(t(obs) %*% Invj %*% obs)
      temp <- temp - (tm1 + tm2)/2
    }
    temp
  }
  thetah <- c(zeta, c(SigmaPhi[!lower.tri(SigmaPhi)]),
              c(SigmaE[!lower.tri(SigmaE)]))
  theta <- thetah
  LLold <- loglik(thetah)

  niter <- 1
  while((niter <= maxiter) & (abs(delta) > epsilon)){
    # Expectation step
    for(j in unique(ID)){
      j_index <- match(j, unique(ID))
      Yj <- Y[ID == j,,drop=FALSE]
      Xj <- X[ID == j,,drop=FALSE]
      residj <- Yj - cbind(Xj%*%beta1, Xj%*%beta2)
      Vj <- solve(InvS2Phi + m[j_index]*InvS2E)
      Muj <- as.numeric(Vj %*% InvS2E %*% colSums(residj))
      Nujj <- Vj + tcrossprod(Muj)
      ESSphi1[j_index,] <- Muj
      ESSphi2[,,j_index] <- Nujj
    }

    # Maximization step - phi
    SigmaPhi <- apply(ESSphi2, 1:2, sum)/n
    InvS2Phi <- solve(SigmaPhi)

    # Maximization step - zeta
    # Simplify the expression analytically, and obtain simple expression!
    XXt <- crossprod(X)
    Vzeta <-solve(kronecker(InvS2E, XXt))

    temp_ID <- cumsum(c(0, diff(ID)) != 0) + 1

    rzeta1 <- t(X)%*%(Y[,1] - ESSphi1[temp_ID,1])
    rzeta2 <- t(X)%*%(Y[,2] - ESSphi1[temp_ID,2])
    zeta <- Vzeta %*% rbind(InvS2E[1,1]*rzeta1 + InvS2E[1,2]*rzeta2,
                            InvS2E[2,1]*rzeta1 + InvS2E[2,2]*rzeta2)
    zeta <- c(zeta)
    beta1 = zeta[1:2]
    beta2 = zeta[3:4]

    # Maximization step - epsilon
    re <- Y - cbind(X%*%beta1, X%*%beta2)
    rss <- crossprod(re) + rowSums(sweep(ESSphi2, 3, m, FUN = "*"), dims = 2) -
      crossprod(ESSphi1, rowsum(re, ID)) - crossprod(rowsum(re, ID), ESSphi1)
    SigmaE <- rss/sum(m)
    # SigmaE <- diag(diag(SigmaE))
    InvS2E <- solve(SigmaE)

    # whether the algorithm converges
    # thetah = c(zeta,c(SigmaPhi[!lower.tri(SigmaPhi)]),diag(SigmaE))
    thetah <- c(zeta, c(SigmaPhi[!lower.tri(SigmaPhi)]),
                c(SigmaE[!lower.tri(SigmaE)]))
    LLnew <- loglik(thetah)
    delta <- abs(LLnew - LLold)
    LLold <- LLnew
    converge = (abs(delta)<=epsilon)
    niter <- niter + 1
    if(verbose) cat(paste('iter=', niter), '\t',
                    paste('param.error=',epsilon), '\t',
                    paste('loglik=',LLnew), '\n');

    #print(niter)
    #print(zeta)
    #print(SigmaPhi)
    #print(SigmaE)
    #print(LLnew)
  }
  param <- list(theta = list(zeta = zeta, SigmaE = SigmaE, SigmaPhi = SigmaPhi),
                loglik = LLnew, eps = epsilon, iter = niter)

  return(param)
}

# calculateSimDataStats() ------------------------------------------------------
calculateSimDataStats <- function(em_output, mlmm_output, my_data){

  # Calculate study sizes
  K_Total_sim <- length(unique(my_data$cluster_k))
  K_Treatment_sim <- length(unique(filter(my_data, treatment_z_k == 1)$cluster_k))
  m_sim <- nrow(filter(my_data, cluster_k == 1))

  # Treatment Allocation Ratio, recall r = K2/K1 where K1 = # in treatment group
  r_sim <- nrow(filter(my_data, treatment_z_k == 0))/
    nrow(filter(my_data, treatment_z_k == 1))
  r_alt_sim <- 1/(r_sim + 1)
  sigmaz.square_sim <- r_alt_sim*(1 - r_alt_sim)

  # Calculate simulation outcome variances (from data directly)
  varY1_data_sim <- var(my_data$Y1)
  varY2_data_sim <- var(my_data$Y2)
  varYc_data_sim <- var(my_data$Yc)

  # MLMM WITH EM ESTIMATION ----------------------------------------------------

  # Variance component matrices
  SigmaE_em <- em_output$theta$SigmaE
  SigmaPhi_em <- em_output$theta$SigmaPhi

  sigma2_e_1_em <- SigmaE_em[1,1]
  sigma2_e_2_em <- SigmaE_em[2,2]
  sigma_e_12_em <- SigmaE_em[1,2]

  sigma2_phi_1_em <- SigmaPhi_em[1,1]
  sigma2_phi_2_em <- SigmaPhi_em[2,2]
  sigma_phi_12_em <- SigmaPhi_em[1,2]

  # Calculate simulation outcome variances (model based)
  varY1_model_sim_em <- sigma2_e_1_em + sigma2_phi_1_em
  varY2_model_sim_em <- sigma2_e_2_em + sigma2_phi_2_em

  # Calculate outcome specific ICC's
  rho01_sim_em <- sigma2_phi_1_em/(sigma2_phi_1_em + sigma2_e_1_em)
  rho02_sim_em <- sigma2_phi_2_em/(sigma2_phi_2_em + sigma2_e_2_em)

  # Calculate Intra- and Inter- correlations
  rho1_sim_em <- sigma_phi_12_em/(sqrt(sigma2_phi_1_em + sigma2_e_1_em)*
                                    sqrt(sigma2_phi_2_em + sigma2_e_2_em))
  rho2_sim_em <- (sigma_phi_12_em + sigma_e_12_em)/(sqrt(sigma2_phi_1_em + sigma2_e_1_em)*
                                                      sqrt(sigma2_phi_2_em + sigma2_e_2_em))

  # Calculate betas
  beta1_sim_em <- em_output$theta$zeta[2]
  beta2_sim_em <- em_output$theta$zeta[4]

  beta1_intercept_em <- em_output$theta$zeta[1]
  beta2_intercept_em <- em_output$theta$zeta[3]

  # Calculate variance of betas, i.e. Var(sqrt(n)*(beta_hat - beta))
  sigmaks.sq_em <- diag(calCovbetas_eq(vars = c(varY1_model_sim_em,
                                                varY2_model_sim_em),
                                       rho01 = matrix(c(rho01_sim_em,
                                                        rho1_sim_em,
                                                        rho1_sim_em,
                                                        rho02_sim_em),
                                                      2, 2),
                                       rho2 = matrix(c(1, rho2_sim_em,
                                                       rho2_sim_em, 1),
                                                     2, 2),
                                       sigmaz.square = sigmaz.square_sim,
                                       m = m_sim,
                                       Q = 2))

  # Variance of betas
  varBeta1_calc_sim_em <- sigmaks.sq_em[1]/K_Total_sim
  varBeta2_calc_sim_em <- sigmaks.sq_em[2]/K_Total_sim

  # MLMM ESTIMATION ------------------------------------------------------------

  # Variance component matrices
  SigmaE_mlmm <- mlmm_output$theta$SigmaE
  SigmaPhi_mlmm <- mlmm_output$theta$SigmaPhi

  sigma2_e_1_mlmm <- SigmaE_mlmm[1,1]
  sigma2_e_2_mlmm <- SigmaE_mlmm[2,2]
  sigma_e_12_mlmm <- SigmaE_mlmm[1,2]

  sigma2_phi_1_mlmm <- SigmaPhi_mlmm[1,1]
  sigma2_phi_2_mlmm <- SigmaPhi_mlmm[2,2]
  sigma_phi_12_mlmm <- SigmaPhi_mlmm[1,2]

  # Calculate simulation outcome variances (model based)
  varY1_model_sim_mlmm <- sigma2_e_1_mlmm + sigma2_phi_1_mlmm
  varY2_model_sim_mlmm <- sigma2_e_2_mlmm + sigma2_phi_2_mlmm

  # Calculate outcome specific ICC's
  rho01_sim_mlmm <- sigma2_phi_1_mlmm/(sigma2_phi_1_mlmm + sigma2_e_1_mlmm)
  rho02_sim_mlmm <- sigma2_phi_2_mlmm/(sigma2_phi_2_mlmm + sigma2_e_2_mlmm)

  # Calculate Intra- and Inter- correlations
  rho1_sim_mlmm <- sigma_phi_12_mlmm/(sqrt(sigma2_phi_1_mlmm + sigma2_e_1_mlmm)*
                                        sqrt(sigma2_phi_2_mlmm + sigma2_e_2_mlmm))
  rho2_sim_mlmm <- (sigma_phi_12_mlmm + sigma_e_12_mlmm)/(sqrt(sigma2_phi_1_mlmm + sigma2_e_1_mlmm)*
                                                            sqrt(sigma2_phi_2_mlmm + sigma2_e_2_mlmm))

  # Calculate betas
  beta1_sim_mlmm <- mlmm_output$theta$zeta[2]
  beta2_sim_mlmm <- mlmm_output$theta$zeta[4]

  beta1_intercept_mlmm <- mlmm_output$theta$zeta[1]
  beta2_intercept_mlmm <- mlmm_output$theta$zeta[3]

  # Calculate variance of betas, i.e. Var(sqrt(n)*(beta_hat - beta))
  sigmaks.sq_mlmm <- diag(calCovbetas_eq(vars = c(varY1_model_sim_mlmm,
                                                  varY2_model_sim_mlmm),
                                         rho01 = matrix(c(rho01_sim_mlmm,
                                                          rho1_sim_mlmm,
                                                          rho1_sim_mlmm,
                                                          rho02_sim_mlmm),
                                                        2, 2),
                                         rho2 = matrix(c(1, rho2_sim_mlmm,
                                                         rho2_sim_mlmm, 1),
                                                       2, 2),
                                         sigmaz.square = sigmaz.square_sim,
                                         m = m_sim,
                                         Q = 2))

  # Variance of betas calculated from the output
  varBeta1_calc_sim_mlmm <- sigmaks.sq_mlmm[1]/K_Total_sim
  varBeta2_calc_sim_mlmm <- sigmaks.sq_mlmm[2]/K_Total_sim

  # Combined outcome stats
  betaC_sim_mlmm <- mlmm_output$combined$betaC[2]
  betaC_intercept_mlmm <- mlmm_output$combined$betaC[1]
  rho0C_sim_mlmm <- mlmm_output$combined$rho0C_sim
  varYc_model_sim_mlmm <- mlmm_output$combined$varYc_model_sim
  varBetaC_model_sim_mlmm <- mlmm_output$combined$varBetaC_sim

  # Calculating Var(Betac) from beta1 and beta2
  varBetaC_calc_sim_mlmm <- varBeta1_calc_sim_mlmm + varBeta2_calc_sim_mlmm +
    (1+(1/r_sim))*(1/(K_Treatment_sim*m_sim))*(rho2_sim_em + ((m_sim - 1)*rho1_sim_em)*sqrt(varY1_model_sim_mlmm*varY2_model_sim_mlmm))

  # PRINTING OUTPUT ------------------------------------------------------------
  # Output specification table
  sim_data_stats <- tibble(Parameter = c("K Total", "K Treatment", "m",
                                         "beta1", "beta1 intercept",
                                         "beta2", "beta2 intercept",
                                         "betaC", "betaC intercept",
                                         "rho01", "rho02", "rho0C",
                                         "rho1", "rho2",
                                         "varY1_model",
                                         "varY2_model",
                                         "varYc_model",
                                         "varBeta1 (calc from output)",
                                         "varBeta2 (calc from output)",
                                         "varBetaC (from model)",
                                         "varBetaC (calc from 1 and 2)",
                                         "r"),
                           `Estimated Value (MLMM with EM)` = c(K_Total_sim,
                                                                K_Treatment_sim,
                                                                m_sim,
                                                                round(beta1_sim_em, 4),
                                                                round(beta1_intercept_em, 4),
                                                                round(beta2_sim_em, 4),
                                                                round(beta2_intercept_em, 4),
                                                                NA, NA,
                                                                round(rho01_sim_em, 4),
                                                                round(rho02_sim_em, 4),
                                                                NA,
                                                                round(rho1_sim_em, 4),
                                                                round(rho2_sim_em, 4),
                                                                round(varY1_model_sim_em, 4),
                                                                round(varY2_model_sim_em, 4),
                                                                NA,
                                                                round(varBeta1_calc_sim_em, 4),
                                                                round(varBeta2_calc_sim_em, 4),
                                                                NA,
                                                                NA,
                                                                r_sim
                           ),
                           `Estimated Value (MLMM without EM)` = c(K_Total_sim,
                                                                   K_Treatment_sim,
                                                                   m_sim,
                                                                   round(beta1_sim_mlmm, 4),
                                                                   round(beta1_intercept_mlmm, 4),
                                                                   round(beta2_sim_mlmm, 4),
                                                                   round(beta2_intercept_mlmm, 4),
                                                                   round(betaC_sim_mlmm, 4),
                                                                   round(betaC_intercept_mlmm, 4),
                                                                   round(rho01_sim_mlmm, 4),
                                                                   round(rho02_sim_mlmm, 4),
                                                                   round(rho0C_sim_mlmm, 4),
                                                                   round(rho1_sim_mlmm, 4),
                                                                   round(rho2_sim_mlmm, 4),
                                                                   round(varY1_model_sim_mlmm, 4),
                                                                   round(varY2_model_sim_mlmm, 4),
                                                                   round(varYc_model_sim_mlmm, 4),
                                                                   round(varBeta1_calc_sim_mlmm, 4),
                                                                   round(varBeta2_calc_sim_mlmm, 4),
                                                                   round(varBetaC_model_sim_mlmm, 4),
                                                                   round(varBetaC_calc_sim_mlmm,4),
                                                                   r_sim
                           )
  )

  return(sim_data_stats)
}

# create_all_cont_sim_dats()
# Creates "n" amount of simulated datasets
create_all_cont_sim_dats <- function(n = 100,
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
    mySimData <- gen_crt_coprimary_data_cont(K = K_input,
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

    # Theoretical values for Combined Outcome
    rho0C_theoretical <- (rho01_input*varY1_input + rho02_input*varY2_input +
                            2*rho1_input*sqrt(varY1_input)*sqrt(varY2_input))/
      (varY1_input + varY2_input + 2*rho2_input*sqrt(varY1_input)*sqrt(varY2_input))
    varYc_theoretical <- varY1_input + varY2_input +
      2*rho2_input*sqrt(varY1_input)*sqrt(varY2_input)
    betaC_theoretical <- beta1_input + beta2_input

    myParamsEM <- EM.estim(myData = mySimData,
                           maxiter = 500,
                           epsilon = 1e-4,
                           verbose = FALSE)

    myParamsMLMM <- MLMM.estim(myData = mySimData)

    simStats <- calculateSimDataStats(em_output = myParamsEM,
                                      mlmm_output = myParamsMLMM,
                                      my_data = mySimData) %>%
      filter(Parameter != "varBeta1", Parameter != "varBeta2") %>%
      mutate(`Theoretical Value` = c(K_input + r_input*K_input,
                                     K_input, m_input,
                                     beta1_input, 0,
                                     beta2_input, 0,
                                     betaC_theoretical, 0,
                                     rho01_input, rho02_input, rho0C_theoretical,
                                     rho1_input, rho2_input,
                                     varY1_input, varY2_input, varYc_theoretical,
                                     NA, NA, NA, NA,
                                     r_input)) %>%
      relocate(`Parameter`, `Theoretical Value`)

    simDataList[[i]] <- mySimData
    simStatList[[i]] <- simStats
  }
  return(list(simDataList, simStatList))
}

# run_cont_sim()
# Creates "n" amount of simulated datasets
run_cont_sim <- function(n = 1000,
                         scenarioTable,
                         null = FALSE){

  # Setting the seed
  set.seed(1972)

  # Initialize final dataset
  cont_sim_data <- NULL

  # Loop through each simulation parameter scenario
  for(currScenario in 1:nrow(scenarioTable)){
    for(currData in 1:n){
      # Specify current set of scenario true parameters
      currTrueParams <- scenarioTable[currScenario,]

      # Create a raw dataset based on current true parameters
      rawSimData <- gen_crt_coprimary_data_cont(K = currTrueParams$K,
                                                m = currTrueParams$m,
                                                beta1 = currTrueParams$beta1,
                                                beta2 = currTrueParams$beta2,
                                                rho01 = currTrueParams$rho01,
                                                rho02 = currTrueParams$rho02,
                                                rho1 = currTrueParams$rho1,
                                                rho2 = currTrueParams$rho2,
                                                varY1 = currTrueParams$varY1,
                                                varY2 = currTrueParams$varY2,
                                                r = currTrueParams$r
      )

      # Run EM estimation on current raw dataset
      currOutputEM <- EM.estim(myData = rawSimData, maxiter = 500,
                               epsilon = 1e-4, verbose = FALSE)

      # Extract parameters from EM output and raw dataset
      currEstParams <- calculateSimDataStats(currOutputEM, rawSimData) %>%
        mutate(Parameter = str_replace(Parameter, " ", "_")) %>%
        pivot_wider(names_from = Parameter, values_from = 'Estimated Value') %>%
        mutate(Scenario_ID = currScenario, Dataset_ID = currData) %>%
        relocate(Scenario_ID, Dataset_ID)

      cont_sim_data <- cont_sim_data %>%
        bind_rows(., currEstParams)

      # Printing progress for n datasets
      cat("Completed", currData, "/", n,
          "datasets for scenario", currScenario, "/", nrow(scenarioTable), "\n")

    } # End loop 1:n

  } # End loop 1:nrow(scenarioTable)

  if(null == FALSE){
    # Figuring out how many simulation dataframes are in the directory
    # specifically for the continuous case
    fileNumber <- length(list.files("./SimulationOutput",
                                    pattern = "cont_sim_data",
                                    all.files = FALSE, recursive = TRUE,
                                    full.names = TRUE))

    # Writing the CSV to the file
    write.csv(cont_sim_data,
              paste0("./SimulationOutput/",
                     "cont_sim_data", "_", fileNumber, ".csv"))
  } else{
    # Figuring out how many simulation dataframes are in the directory
    # specifically for the continuous case
    fileNumber <- length(list.files("./SimulationOutput",
                                    pattern = "cont_sim_null_data",
                                    all.files = FALSE, recursive = TRUE,
                                    full.names = TRUE))

    # Writing the CSV to the file
    write.csv(cont_sim_data,
              paste0("./SimulationOutput/",
                     "cont_sim_null_data", "_", fileNumber, ".csv"))
  }
}
