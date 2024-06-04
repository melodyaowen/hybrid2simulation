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

simParameters_null <- expand.grid(K = c(8),
                                  m = c(60),
                                  betas = c(paste("0 0.4"),
                                            paste("0 0.2")),
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
  mutate('method1_bonf' = calc_pwr_pval_adj(K = K, m = m, alpha = alpha,
                                            beta1 = beta1, beta2 = beta2,
                                            varY1 = varY1, varY2 = varY2,
                                            rho01 = rho01, rho02 = rho02,
                                            rho2  = rho2, r = r)$'Final Power'[[1]],
         'method1_sidak' = calc_pwr_pval_adj(K = K, m = m, alpha = alpha,
                                             beta1 = beta1, beta2 = beta2,
                                             varY1 = varY1, varY2 = varY2,
                                             rho01 = rho01, rho02 = rho02,
                                             rho2  = rho2, r = r)$'Final Power'[[2]],
         'method1_dap' = calc_pwr_pval_adj(K = K, m = m, alpha = alpha,
                                           beta1 = beta1, beta2 = beta2,
                                           varY1 = varY1, varY2 = varY2,
                                           rho01 = rho01, rho02 = rho02,
                                           rho2  = rho2, r = r)$'Final Power'[[3]],
         'method2' = calc_pwr_comb_outcome(K = K, m = m, alpha = alpha,
                                           beta1 = beta1, beta2 = beta2,
                                           varY1 = varY1, varY2 = varY2,
                                           rho01 = rho01, rho02 = rho02,
                                           rho1 = rho1, rho2  = rho2, r = r),
         'method3' = calc_pwr_single_1dftest(K = K, m = m, alpha = alpha,
                                             beta1 = beta1, beta2 = beta2,
                                             varY1 = varY1, varY2 = varY2,
                                             rho01 = rho01, rho02 = rho02,
                                             rho1 = rho1, rho2  = rho2, r = r),
         'method4_chi2' = calc_pwr_disj_2dftest(K = K, m = m, alpha = alpha,
                                                beta1 = beta1, beta2 = beta2,
                                                varY1 = varY1, varY2 = varY2,
                                                rho01 = rho01, rho02 = rho02,
                                                rho1 = rho1, rho2  = rho2,
                                                r = r, dist = "Chi2"),
         'method4_F' = calc_pwr_disj_2dftest(K = K, m = m, alpha = alpha,
                                             beta1 = beta1, beta2 = beta2,
                                             varY1 = varY1, varY2 = varY2,
                                             rho01 = rho01, rho02 = rho02,
                                             rho1 = rho1, rho2  = rho2,
                                             r = r, dist = "F"),
         'method5_T' = calc_pwr_conj_test(K = K, m = m, alpha = alpha,
                                          beta1 = beta1, beta2 = beta2,
                                          varY1 = varY1, varY2 = varY2,
                                          rho01 = rho01, rho02 = rho02,
                                          rho1 = rho1, rho2  = rho2,
                                          r = r, dist = "T"),
         'method5_MVN' = calc_pwr_conj_test(K = K, m = m, alpha = alpha,
                                            beta1 = beta1, beta2 = beta2,
                                            varY1 = varY1, varY2 = varY2,
                                            rho01 = rho01, rho02 = rho02,
                                            rho1 = rho1, rho2  = rho2,
                                            r = r, dist = "MVN")) %>%
  mutate_at(vars(contains('method')), funs(.*100))

truePowerTable_null <- simParameters_null %>%
  rowwise() %>%
  mutate('method1_bonf' = calc_pwr_pval_adj(K = K, m = m, alpha = alpha,
                                            beta1 = beta1, beta2 = beta2,
                                            varY1 = varY1, varY2 = varY2,
                                            rho01 = rho01, rho02 = rho02,
                                            rho2  = rho2, r = r)$'Final Power'[[1]],
         'method1_sidak' = calc_pwr_pval_adj(K = K, m = m, alpha = alpha,
                                             beta1 = beta1, beta2 = beta2,
                                             varY1 = varY1, varY2 = varY2,
                                             rho01 = rho01, rho02 = rho02,
                                             rho2  = rho2, r = r)$'Final Power'[[2]],
         'method1_dap' = calc_pwr_pval_adj(K = K, m = m, alpha = alpha,
                                           beta1 = beta1, beta2 = beta2,
                                           varY1 = varY1, varY2 = varY2,
                                           rho01 = rho01, rho02 = rho02,
                                           rho2  = rho2, r = r)$'Final Power'[[3]],
         'method2' = calc_pwr_comb_outcome(K = K, m = m, alpha = alpha,
                                           beta1 = beta1, beta2 = beta2,
                                           varY1 = varY1, varY2 = varY2,
                                           rho01 = rho01, rho02 = rho02,
                                           rho1 = rho1, rho2  = rho2, r = r),
         'method3' = calc_pwr_single_1dftest(K = K, m = m, alpha = alpha,
                                             beta1 = beta1, beta2 = beta2,
                                             varY1 = varY1, varY2 = varY2,
                                             rho01 = rho01, rho02 = rho02,
                                             rho1 = rho1, rho2  = rho2, r = r),
         'method4_chi2' = calc_pwr_disj_2dftest(K = K, m = m, alpha = alpha,
                                                beta1 = beta1, beta2 = beta2,
                                                varY1 = varY1, varY2 = varY2,
                                                rho01 = rho01, rho02 = rho02,
                                                rho1 = rho1, rho2  = rho2,
                                                r = r, dist = "Chi2"),
         'method4_F' = calc_pwr_disj_2dftest(K = K, m = m, alpha = alpha,
                                             beta1 = beta1, beta2 = beta2,
                                             varY1 = varY1, varY2 = varY2,
                                             rho01 = rho01, rho02 = rho02,
                                             rho1 = rho1, rho2  = rho2,
                                             r = r, dist = "F"),
         'method5_T' = calc_pwr_conj_test(K = K, m = m, alpha = alpha,
                                          beta1 = beta1, beta2 = beta2,
                                          varY1 = varY1, varY2 = varY2,
                                          rho01 = rho01, rho02 = rho02,
                                          rho1 = rho1, rho2  = rho2,
                                          r = r, dist = "T"),
         'method5_MVN' = calc_pwr_conj_test(K = K, m = m, alpha = alpha,
                                            beta1 = beta1, beta2 = beta2,
                                            varY1 = varY1, varY2 = varY2,
                                            rho01 = rho01, rho02 = rho02,
                                            rho1 = rho1, rho2  = rho2,
                                            r = r, dist = "MVN")) %>%
  mutate_at(vars(contains('method')), funs(.*100))

View(truePowerTable)
View(truePowerTable_null)

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
#run_cont_sim(n = 2000, scenarioTable = simParameters)
#run_cont_sim(n = 2000, scenarioTable = simParameters, null = TRUE)

# Read in Continuous Datasets --------------------------------------------------

# Reading in the CSV file with sim results (not null)
contResults <- read.csv("./SimulationOutput/cont_sim_data_0.csv") %>%
  dplyr::select(-X)

# Reading in the CSV file with sim results (null hypothesis case)
contResults_null <- read.csv("./SimulationOutput/cont_sim_null_data_0.csv") %>%
  dplyr::select(-X)

# Evaluate Power for Continuous Data -------------------------------------------

# Dataset of true values
trueVals <- simParameters %>%
  dplyr::select(-K, -m, -r, -alpha) %>%
  rename(Scenario_ID = "Scenario") %>%
  rename_at(vars(-Scenario_ID), function(x) paste0(x,"_true"))

# Comparing estimated parameters with true parameters
paramStats <- contResults %>%
  dplyr::select(-Dataset_ID, -K_Total, -K_Treatment, -m, -r,
                -varBeta1, -varBeta2) %>%
  group_by(Scenario_ID) %>%
  summarise(across(everything(), list(mean = mean#,
                                      #min = min,
                                      #max = max,
                                      #sd = sd
                                      ))) %>%
  left_join(trueVals, by = "Scenario_ID") %>%
  mutate(beta1_intercept_true = 0, beta2_intercept_true = 0) %>%
  dplyr::select(order(colnames(.))) %>%
  relocate(Scenario_ID) %>%
  mutate(across(!contains(c("true", "ID")), ~ round(.x, digits = 4)))

# Calculating power based on all the simulated datasets
contResultsPower <- contResults %>%
  rowwise() %>%
  mutate(alpha = 0.05) %>%
  mutate('method1_bonf' = calc_pwr_pval_adj(K = K_Treatment, m = m, alpha = alpha,
                                            beta1 = beta1, beta2 = beta2,
                                            varY1 = varY1, varY2 = varY2,
                                            rho01 = rho01, rho02 = rho02,
                                            rho2  = rho2, r = r)$'Final Power'[[1]],
         'method1_sidak' = calc_pwr_pval_adj(K = K_Treatment, m = m, alpha = alpha,
                                             beta1 = beta1, beta2 = beta2,
                                             varY1 = varY1, varY2 = varY2,
                                             rho01 = rho01, rho02 = rho02,
                                             rho2  = rho2, r = r)$'Final Power'[[2]],
         'method1_dap' = calc_pwr_pval_adj(K = K_Treatment, m = m, alpha = alpha,
                                           beta1 = beta1, beta2 = beta2,
                                           varY1 = varY1, varY2 = varY2,
                                           rho01 = rho01, rho02 = rho02,
                                           rho2  = rho2, r = r)$'Final Power'[[3]],
         'method2' = calc_pwr_comb_outcome(K = K_Treatment, m = m, alpha = alpha,
                                           beta1 = beta1, beta2 = beta2,
                                           varY1 = varY1, varY2 = varY2,
                                           rho01 = rho01, rho02 = rho02,
                                           rho1 = rho1, rho2  = rho2, r = r),
         'method3' = calc_pwr_single_1dftest(K = K_Treatment, m = m, alpha = alpha,
                                             beta1 = beta1, beta2 = beta2,
                                             varY1 = varY1, varY2 = varY2,
                                             rho01 = rho01, rho02 = rho02,
                                             rho1 = rho1, rho2  = rho2, r = r),
         'method4_chi2' = calc_pwr_disj_2dftest(K = K_Treatment, m = m, alpha = alpha,
                                                beta1 = beta1, beta2 = beta2,
                                                varY1 = varY1, varY2 = varY2,
                                                rho01 = rho01, rho02 = rho02,
                                                rho1 = rho1, rho2  = rho2,
                                                r = r, dist = "Chi2"),
         'method4_F' = calc_pwr_disj_2dftest(K = K_Treatment, m = m, alpha = alpha,
                                             beta1 = beta1, beta2 = beta2,
                                             varY1 = varY1, varY2 = varY2,
                                             rho01 = rho01, rho02 = rho02,
                                             rho1 = rho1, rho2  = rho2,
                                             r = r, dist = "F"),
         'method5_T' = calc_pwr_conj_test(K = K_Treatment, m = m, alpha = alpha,
                                          beta1 = beta1, beta2 = beta2,
                                          varY1 = varY1, varY2 = varY2,
                                          rho01 = rho01, rho02 = rho02,
                                          rho1 = rho1, rho2  = rho2,
                                          r = r, dist = "T"),
         'method5_MVN' = calc_pwr_conj_test(K = K_Treatment, m = m, alpha = alpha,
                                            beta1 = beta1, beta2 = beta2,
                                            varY1 = varY1, varY2 = varY2,
                                            rho01 = rho01, rho02 = rho02,
                                            rho1 = rho1, rho2  = rho2,
                                            r = r, dist = "MVN")) %>%
  mutate_at(vars(contains('method')), funs(.*100))

contResultsPowerOnly <- contResultsPower %>%
  dplyr::select(Scenario_ID, Dataset_ID, contains("method"))

# Summary data favstats of power
contResultsPowerSummary <- contResultsPowerOnly %>%
  dplyr::select(-Dataset_ID) %>%
  group_by(Scenario_ID) %>%
  summarise(across(everything(), list(mean = mean,
                                      min = min,
                                      max = max,
                                      sd = sd)))

# Boxplot data
contPlotData <- contResultsPowerOnly %>%
  pivot_longer(cols = starts_with("method"),
               names_to = "Method", values_to = "Power") %>%
  mutate(Scenario_ID = paste0("Scenario ", Scenario_ID)) %>%
  mutate(Method = recode(Method,
                         "method1_bonf" = "P-Value Bonf.",
                         "method1_sidak" = "P-Value Sidak",
                         "method1_dap" = "P-Value DAP",
                         "method2" = "Combined",
                         "method3" = "Single Weighted",
                         "method4_chi2" = "Disjunctive Chi2",
                         "method4_F" = "Disjunctive F",
                         "method5_T" = "Conjunctive T",
                         "method5_MVN" = "Conjunctive MVN")) %>%
  mutate(Method = factor(Method, levels = c("P-Value Bonf.", "P-Value Sidak",
                                            "P-Value DAP", "Combined",
                                            "Single Weighted",
                                            "Disjunctive Chi2", "Disjunctive F",
                                            "Conjunctive T", "Conjunctive MVN")),
         Scenario_ID = factor(Scenario_ID, c("Scenario 1", "Scenario 2",
                                             "Scenario 3", "Scenario 4",
                                             "Scenario 5", "Scenario 6",
                                             "Scenario 7", "Scenario 8",
                                             "Scenario 9", "Scenario 10",
                                             "Scenario 11", "Scenario 12")))

# True points data for plotting
truePowerTablePlot <- truePowerTable %>%
  dplyr::select(Scenario_ID = Scenario, starts_with("method")) %>%
  mutate(Scenario_ID = paste0("Scenario ", Scenario_ID)) %>%
  pivot_longer(cols = starts_with("method"),
               names_to = "Method", values_to = "Power") %>%
  mutate(Method = recode(Method,
                         "method1_bonf" = "P-Value Bonf.",
                         "method1_sidak" = "P-Value Sidak",
                         "method1_dap" = "P-Value DAP",
                         "method2" = "Combined",
                         "method3" = "Single Weighted",
                         "method4_chi2" = "Disjunctive Chi2",
                         "method4_F" = "Disjunctive F",
                         "method5_T" = "Conjunctive T",
                         "method5_MVN" = "Conjunctive MVN")) %>%
  mutate(Method = factor(Method, levels = c("P-Value Bonf.", "P-Value Sidak",
                                            "P-Value DAP", "Combined",
                                            "Single Weighted",
                                            "Disjunctive Chi2", "Disjunctive F",
                                            "Conjunctive T", "Conjunctive MVN")),
         Scenario_ID = factor(Scenario_ID, c("Scenario 1", "Scenario 2",
                                             "Scenario 3", "Scenario 4",
                                             "Scenario 5", "Scenario 6",
                                             "Scenario 7", "Scenario 8",
                                             "Scenario 9", "Scenario 10",
                                             "Scenario 11", "Scenario 12")))

# Boxplots of power for each scenario and method
ggplot(data = contPlotData, aes(x = Method, y = Power)) + geom_boxplot() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  geom_point(data = truePowerTablePlot, aes(x = Method, y = Power),
             color = 'blue') +
  facet_wrap(~Scenario_ID) +
  ggtitle("Power Results from Co-Primary Continuous CRT Simulation\n(Blue points are true power based on true parameters)")

# Table summarizing power results
resultTable <- contResultsPower %>%
  group_by(Scenario_ID) %>%
  summarise(across(starts_with("method"), list(mean = mean))) %>%
  ungroup() %>%
  left_join(truePowerTable, by = c("Scenario_ID" = "Scenario")) %>%
  mutate(betas = paste0("(", beta1, ", ", beta2, ")"),
         rho0s = paste0("(", rho01, ", ", rho02, ")"),
         varYs = paste0("(", varY1, ", ", varY2, ")")) %>%
  dplyr::select(Scenario_ID,
                K, m, betas,
                rho0s, rho1, rho2, varYs,
                method1_bonf_mean, method1_bonf,
                method1_sidak_mean, method1_sidak,
                method1_dap_mean, method1_dap,
                method1_bonf_mean, method1_bonf_mean,
                method2_mean, method2,
                method3_mean, method3,
                method4_F_mean, method4_F,
                method4_chi2_mean, method4_chi2,
                method5_T_mean, method5_T,
                method5_MVN_mean, method5_MVN) %>%
  mutate(across(starts_with("method"), round, 2)) %>%
  dplyr::rename('P-Value Bonf. (mean)' = method1_bonf_mean,
                'P-Value Bonf. (true)' = method1_bonf,
                'P-Value Sidak (mean)' = method1_sidak_mean,
                'P-Value Sidak (true)' = method1_sidak,
                'P-Value DAP (mean)' = method1_dap_mean,
                'P-Value DAP (true)' = method1_dap,
                'Combined (mean)' = method2_mean,
                'Combined (true)' = method2,
                'Single Weighted (mean)' = method3_mean,
                'Single Weighted (true)' = method3,
                'Disjunctive Chi2 (mean)' = method4_chi2_mean,
                'Disjunctive Chi2 (true)' = method4_chi2,
                'Disjunctive F (mean)' = method4_F_mean,
                'Disjunctive F (true)' = method4_F,
                'Conjunctive T (mean)' = method5_T_mean,
                'Conjunctive T (true)' = method5_T,
                'Conjunctive MVN (mean)' = method5_MVN_mean,
                'Conjunctive MVN (true)' = method5_MVN)

# Evaluate Type I Error Rate for Continuous Data -------------------------------
# Type I error, probability of rejecting null when it is true

typeIerror_data <- contResults_null %>%
  add_row(Scenario_ID = 0, Dataset_ID = 0,
          K_Total = 30, K_Treatment = 15, m = 300,
          beta1 = 0.1, beta2 = 0.1, beta1_intercept = 0, beta2_intercept = 0,
          rho01 = 0.025, rho02 = 0.025, rho1 = 0.01, rho2 = 0.05,
          varY1 = 0.23, varY2 = 0.25, varBeta1 = NA, varBeta2 = NA, r = 1) %>%
  mutate(VIF1 = 1+(m-1)*rho01,
         VIF2 = 1+(m-1)*rho02,
         VIF12 = rho2 + (m-1)*rho1,
         betaC = abs(beta1) + abs(beta2),
         varYc = varY1 + varY2 + 2*rho2*sqrt(varY1)*sqrt(varY2),
         rho0c = (rho01*varY1 + rho02*varY2 + 2*rho1*sqrt(varY1)*sqrt(varY2))/(varY1 + varY2 + 2*rho2*sqrt(varY1)*sqrt(varY2))) %>%
  mutate(method1_stat1 = (beta1^2)/varBeta1,
         method1_stat2 = (beta2^2)/varBeta2,
         #method2_stat = ((abs(beta1) + abs(beta2))^2)/(varBeta1 + varBeta2 + (1 + 1/r)*(1/(K_Treatment*m))*((VIF12)*sqrt(varY1*varY2))),
         method2_stat = (betaC^2)/(2*(varYc/(K_Treatment*m))*(1 + (m-1)*rho0c)),
         method3_stat = ((sqrt((beta1^2)/(2*(varY1/(K_Treatment*m))*VIF1))+sqrt((beta2^2)/(2*(varY2/(K_Treatment*m))*VIF2)))^2)/(2*(1 + VIF12/(sqrt(VIF1*VIF2)))),
         method4_stat = (K_Treatment*m*((beta1^2)*varY2*VIF2 - 2*beta1*beta2*VIF12*sqrt(varY1)*sqrt(varY2) + (beta2^2)*varY1*VIF1))/(2*2*varY1*varY2*(VIF1*VIF2 - VIF12^2)),
         method5_stat1 = (beta1*sqrt(2*K_Treatment))/sqrt((4*varY1*VIF1)/m),
         method5_stat2 = (beta2*sqrt(2*K_Treatment))/sqrt((4*varY2*VIF2)/m)
         )

