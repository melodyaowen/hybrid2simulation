source("./DataGeneration/Generate_Data_2_Cont_Outcomes.R")
source("./RequiredPackages.R")

# Table of all Parameters ------------------------------------------------------
simParameters <- expand.grid(K = c(15),
                             m = c(60),
                             betas = c(paste("0.3 0.7"),
                                       paste("0.3 0.3")),
                             vars = c(paste("1 1")),
                             rho0 = c(paste("0.01 0.1"),
                                      paste("0.1 0.1"),
                                      paste("0.1 0.01")),
                             rho1 = c(0.005),
                             rho2 = c(0.2, 0.5),
                             alpha = c(0.05),
                             r = c(1)
                             ) %>%
  separate(betas, sep = " ", into = c("beta1", "beta2")) %>%
  separate(rho0, sep = " ", into = c("rho01", "rho02")) %>%
  separate(vars, sep = " ", into = c("varY1", "varY2")) %>%
  mutate_if(is.character, as.numeric)

# Run power calculations on all true parameters --------------------------------
truePowerTable <- simParameters %>%
  rowwise() %>%
  mutate('P-Value Bonferroni Power' = calc_pwr_pval_adj(K = K, m = m,
                                                        alpha = alpha,
                                                        beta1 = beta1,
                                                        beta2 = beta2,
                                                        varY1 = varY1,
                                                        varY2 = varY2,
                                                        rho01 = rho01,
                                                        rho02 = rho02,
                                                        rho2  = rho2,
                                                        r = r)$'Final Power'[[1]],
         'P-Value Sidak Power' = calc_pwr_pval_adj(K = K, m = m,
                                                   alpha = alpha,
                                                   beta1 = beta1,
                                                   beta2 = beta2,
                                                   varY1 = varY1,
                                                   varY2 = varY2,
                                                   rho01 = rho01,
                                                   rho02 = rho02,
                                                   rho2  = rho2,
                                                   r = r)$'Final Power'[[2]],
         'P-Value D/AP Power' = calc_pwr_pval_adj(K = K, m = m,
                                                  alpha = alpha,
                                                  beta1 = beta1,
                                                  beta2 = beta2,
                                                  varY1 = varY1,
                                                  varY2 = varY2,
                                                  rho01 = rho01,
                                                  rho02 = rho02,
                                                  rho2  = rho2,
                                                  r = r)$'Final Power'[[3]],
         'Combined Outcome Power' = calc_pwr_comb_outcome(K = K, m = m,
                                                          alpha = alpha,
                                                          beta1 = beta1,
                                                          beta2 = beta2,
                                                          varY1 = varY1,
                                                          varY2 = varY2,
                                                          rho01 = rho01,
                                                          rho02 = rho02,
                                                          rho1 = rho1,
                                                          rho2  = rho2,
                                                          r = r),
         'Weighted 1-DF Test Power' = calc_pwr_single_1dftest(K = K, m = m,
                                                              alpha = alpha,
                                                              beta1 = beta1,
                                                              beta2 = beta2,
                                                              varY1 = varY1,
                                                              varY2 = varY2,
                                                              rho01 = rho01,
                                                              rho02 = rho02,
                                                              rho1 = rho1,
                                                              rho2  = rho2,
                                                              r = r),
         'Disjunctive 2-DF Test Power' = calc_pwr_disj_2dftest(K = K, m = m,
                                                               alpha = alpha,
                                                               beta1 = beta1,
                                                               beta2 = beta2,
                                                               varY1 = varY1,
                                                               varY2 = varY2,
                                                               rho01 = rho01,
                                                               rho02 = rho02,
                                                               rho1 = rho1,
                                                               rho2  = rho2,
                                                               r = r),
         'Conjunctive IU Test Power' = calc_pwr_conj_test(K = K, m = m,
                                                          alpha = alpha,
                                                          beta1 = beta1,
                                                          beta2 = beta2,
                                                          varY1 = varY1,
                                                          varY2 = varY2,
                                                          rho01 = rho01,
                                                          rho02 = rho02,
                                                          rho1 = rho1,
                                                          rho2  = rho2,
                                                          r = r)) %>%
  mutate_at(vars(contains('Power')), funs(.*100)) %>%
  dplyr::select(-alpha, -r, -varY1, -varY2)

View(truePowerTable)

# Generate Continuous Datasets -------------------------------------------------
simSetsCont <- create_all_cont_sim_dats(n = 1,
                                        K_input = 10,
                                        m_input = 500,
                                        beta1_input = 0.5,
                                        beta2_input = 0.1,
                                        rho01_input = 0.07,
                                        rho02_input = 0.05,
                                        rho1_input = 0.02,
                                        rho2_input = 0.25,
                                        varY1_input = 0.3,
                                        varY2_input = 0.2,
                                        r_input = 1)

combinedSimStatsCont <- bind_rows(simSetsCont[[2]], .id = "Index") %>%
  mutate(Index = as.integer(Index)) %>%
  mutate(`Absolute Difference` = abs(`True Value` - `Estimated Value`))

View(simSetsCont[[1]][[1]])
View(simSetsCont[[2]][[1]])





