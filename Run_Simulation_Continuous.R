source("./DataGeneration/Generate_Data_2_Cont_Outcomes.R")
source("./RequiredPackages.R")

# Table of all Parameters ------------------------------------------------------
simParameters <- expand.grid(K = c(8),
                             m = c(60),
                             betas = c(paste("0.2 0.4"),
                                       paste("0.2 0.2")),
                             vars = c(paste("1 1")),
                             rho0 = c(paste("0.05 0.1"),
                                      paste("0.1 0.1"),
                                      paste("0.1 0.05")),
                             rho1 = c(0.005),
                             rho2 = c(0.2, 0.5),
                             alpha = c(0.05),
                             r = c(1)
                             ) %>%
  separate(betas, sep = " ", into = c("beta1", "beta2")) %>%
  separate(rho0, sep = " ", into = c("rho01", "rho02")) %>%
  separate(vars, sep = " ", into = c("varY1", "varY2")) %>%
  mutate_if(is.character, as.numeric) %>%
  rowid_to_column(., "Scenario")

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
  mutate_at(vars(contains('Power')), funs(.*100))

View(truePowerTable)

# Check a few cases to make sure model estimates are ok ------------------------
simSetsCont <- create_all_cont_sim_dats(n = 2,
                                        K_input = 8,
                                        m_input = 60,
                                        beta1_input = 0.2,
                                        beta2_input = 0.4,
                                        rho01_input = 0.05,
                                        rho02_input = 0.1,
                                        rho1_input = 0.005,
                                        rho2_input = 0.2,
                                        varY1_input = 1,
                                        varY2_input = 1,
                                        r_input = 1)

combinedSimStatsCont <- bind_rows(simSetsCont[[2]], .id = "Index") %>%
  mutate(Index = as.integer(Index)) %>%
  mutate(`Absolute Difference` = abs(`True Value` - `Estimated Value`))

View(simSetsCont[[1]][[1]])
View(simSetsCont[[2]][[1]])

# Generate Continuous Datasets -------------------------------------------------
run_cont_sim(n = 5, scenarioTable = simParameters)

contResults <- read.csv("./SimulationOutput/cont_sim_data_0.csv") %>%
  dplyr::select(-X)

contResultsPower <- contResults %>%
  rowwise() %>%
  mutate(alpha = 0.05) %>%
  mutate('P-Value Bonferroni Power' = calc_pwr_pval_adj(K = K_Treatment, m = m,
                                                        alpha = alpha,
                                                        beta1 = beta1,
                                                        beta2 = beta2,
                                                        varY1 = varY1,
                                                        varY2 = varY2,
                                                        rho01 = rho01,
                                                        rho02 = rho02,
                                                        rho2  = rho2,
                                                        r = r)$'Final Power'[[1]],
         'P-Value Sidak Power' = calc_pwr_pval_adj(K = K_Treatment, m = m,
                                                   alpha = alpha,
                                                   beta1 = beta1,
                                                   beta2 = beta2,
                                                   varY1 = varY1,
                                                   varY2 = varY2,
                                                   rho01 = rho01,
                                                   rho02 = rho02,
                                                   rho2  = rho2,
                                                   r = r)$'Final Power'[[2]],
         'P-Value D/AP Power' = calc_pwr_pval_adj(K = K_Treatment, m = m,
                                                  alpha = alpha,
                                                  beta1 = beta1,
                                                  beta2 = beta2,
                                                  varY1 = varY1,
                                                  varY2 = varY2,
                                                  rho01 = rho01,
                                                  rho02 = rho02,
                                                  rho2  = rho2,
                                                  r = r)$'Final Power'[[3]],
         'Combined Outcome Power' = calc_pwr_comb_outcome(K = K_Treatment, m = m,
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
         'Weighted 1-DF Test Power' = calc_pwr_single_1dftest(K = K_Treatment, m = m,
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
         'Disjunctive 2-DF Test Power' = calc_pwr_disj_2dftest(K = K_Treatment, m = m,
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
         'Conjunctive IU Test Power' = calc_pwr_conj_test(K = K_Treatment, m = m,
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
  dplyr::select(Scenario_ID, Dataset_ID, contains("Power"))

View(contResultsPower)

