# e_temp <- mvrnorm(m, mu = rep(0, 10), Sigma = C) #%>%
# as.data.frame() %>%
#   dplyr::select(e_kj1 = V1, e_kj2 = V2) %>%
#   mutate(subject_id = 1:study_n,
#          cluster_id = sort(rep(1:K_total, m)))


# VarMat1 <- list() # Cov(Y_kq)
# VarMat2 <- list() # Cov(Y_kq)
# CovMat <- list() # Cov(Y_kq, Y_kq')
#
# for(k in 1:nrow(groups)){ # Loop through each cluster,
#   # Specify V_k based on trt group
#   if(groups$trt_group[k] == 1){ # Treatment for this cluster is 1
#     v_k1 <- (V1_k_1^(1/2)) %*% C11 %*% (V1_k_1^(1/2))
#     v_k2 <- (V2_k_1^(1/2)) %*% C22 %*% (V2_k_1^(1/2))
#     v_k12 <- (V1_k_1^(1/2)) %*% C12 %*% (V2_k_1^(1/2))
#   } else if(groups$trt_group[k] == 0){ # Trt for this cluster is 0
#     v_k1 <- (V1_k_0^(1/2)) %*% C11 %*% (V1_k_0^(1/2))
#     v_k2 <- (V2_k_0^(1/2)) %*% C22 %*% (V2_k_0^(1/2))
#     v_k12 <- (V1_k_0^(1/2)) %*% C12 %*% (V2_k_0^(1/2))
#   }
#   VarMat1[[k]] <- v_k1
#   VarMat2[[k]] <- v_k2
#   CovMat[[k]] <- v_k12
# }

# Initialize lists of length K for V_k and R_k
#V_k <- list()
#R_k <- list()




simDat <- tibble(subject_id = 1:study_n,
                 cluster_id = sort(rep(1:K_total, m))) %>%
  left_join(., groups, by = "cluster_id") %>%
  left_join(., phi_k, by = "cluster_id") %>%
  left_join(., e_kj, by = c("cluster_id", "subject_id")) %>%
  # Now create outcome columns based on inputs and calculated intercepts
  mutate(eta1 = beta1*trt_group + phi_k1,
         eta2 = beta2*trt_group + phi_k2) %>%
  mutate(prob_Y1 = logistic_inverse_fun(eta1) + e_kj1,
         prob_Y2 = logistic_inverse_fun(eta2) + e_kj2) %>%
  mutate(Y1 = rbinom(m*K*2, size = 1, prob = prob_Y1),
         Y2 = rbinom(m*K*2, size = 1, prob = prob_Y2)) %>%
  dplyr::rename(subject_j = subject_id,
                cluster_k = cluster_id,
                treatment_z_k = trt_group)
# Return data
return(simDat)

# Variance matrices
Var_Y1_k_0 <- sqrt(V1_k_0) %*% C11 %*% sqrt(V1_k_0)
Var_Y1_k_1 <- sqrt(V1_k_1) %*% C11 %*% sqrt(V1_k_1)
Var_Y2_k_0 <- sqrt(V2_k_0) %*% C22 %*% sqrt(V2_k_0)
Var_Y2_k_1 <- sqrt(V2_k_1) %*% C22 %*% sqrt(V2_k_1)

Cov_Y1Y2_k_0 <- (V1_k_0^(1/2)) %*% C12 %*% (V2_k_0^(1/2))
Cov_Y1Y2_k_1 <- (V1_k_1^(1/2)) %*% C12 %*% (V2_k_1^(1/2))



# Generate matrix R, the variance-covariance matrix for the residual errors
R <- V^(1/2) %>% C %>% V^(1/2)



# Generate vector of random erros for eta_q,kj, the residual error terms
# that captures any remaining variance in the response not accounted for
# by the fixed and random covariates
e_kj <- mvrnorm(m*K_total, mu = rep(0, 2), Sigma = R) %>%
  as.data.frame() %>%
  dplyr::select(e_kj1 = V1, e_kj2 = V2) %>%
  mutate(subject_id = 1:study_n,
         cluster_id = sort(rep(1:K_total, m)))




rmvbin(1, margprob = c(0.3,0.9))
rmvbin(100, sigma = C, margprob = c(0.66, 0.6))
