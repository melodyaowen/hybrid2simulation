source("./DataGeneration/Generate_Data_2_Cont_Outcomes.R")
source("./RequiredPackages.R")

# Table of all Parameters ------------------------------------------------------
simParameters <- expand.grid(K = c(8),
                             m = c(60),
                             betas = c(paste("0.25 0.4"), # considered 0.25
                                       paste("0.25 0.25")), # and and 0.2
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
                                  betas = c(paste("0 0"),
                                            paste("0 0.25"), # 0.25 or 0.2
                                            paste("0 0.4")),
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

#View(simSetsCont[[1]][[1]])
#View(simSetsCont[[2]][[1]])

# Generate Continuous Datasets -------------------------------------------------
#run_cont_sim(n = 2000, scenarioTable = simParameters)
#run_cont_sim(n = 2000, scenarioTable = simParameters_null, null = TRUE)

# Read in Continuous Datasets --------------------------------------------------

# Reading in the CSV file with sim results (not null)
# contResults <- read.csv("./SimulationOutput/cont_sim_data_0.csv") %>%
#   dplyr::select(-X)

contResults <- read.csv("./SimulationOutput/cont_sim_data_1.csv") %>%
  dplyr::select(-X)

# Reading in the CSV file with sim results (null hypothesis case)
# contResults_null <- read.csv("./SimulationOutput/cont_sim_null_data_0.csv") %>%
#   dplyr::select(-X)

contResults_null <- read.csv("./SimulationOutput/cont_sim_null_data_1.csv") %>%
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
powerResultTable <- contResultsPower %>%
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
  dplyr::rename('1a. P-Value Bonf. (mean)' = method1_bonf_mean,
                '1a. P-Value Bonf. (true)' = method1_bonf,
                '1b. P-Value Sidak (mean)' = method1_sidak_mean,
                '1b. P-Value Sidak (true)' = method1_sidak,
                '1c. P-Value DAP (mean)' = method1_dap_mean,
                '1c. P-Value DAP (true)' = method1_dap,
                '2. Combined (mean)' = method2_mean,
                '2. Combined (true)' = method2,
                '3. Single Weighted (mean)' = method3_mean,
                '3. Single Weighted (true)' = method3,
                '4a. Disjunctive F (mean)' = method4_F_mean,
                '4a. Disjunctive F (true)' = method4_F,
                '4b. Disjunctive Chi2 (mean)' = method4_chi2_mean,
                '4b. Disjunctive Chi2 (true)' = method4_chi2,
                '5a. Conjunctive T (mean)' = method5_T_mean,
                '5a. Conjunctive T (true)' = method5_T,
                '5b. Conjunctive MVN (mean)' = method5_MVN_mean,
                '5b. Conjunctive MVN (true)' = method5_MVN) %>%
  arrange(betas, rho0s)

# Evaluate Type I Error Rate for Continuous Data -------------------------------
# Type I error, probability of rejecting null when it is true

typeIerror_data <- contResults_null %>%
  left_join(dplyr::select(simParameters_null, Scenario_ID = Scenario,
                          beta1_true = beta1, beta2_true = beta2),
            by = "Scenario_ID") %>%
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
         rho0c = (rho01*varY1 + rho02*varY2 +
                    2*rho1*sqrt(varY1)*sqrt(varY2))/(varY1 + varY2 +
                                          2*rho2*sqrt(varY1)*sqrt(varY2)),
         alpha_bonf = 0.05/2,
         alpha_sidak = 1 - sqrt(1 - 0.05),
         alpha_dap = 1 - (1 - 0.05)^(1/(2^(1 - rho2)))) %>%
  mutate(method1_stat1 = (beta1^2)/varBeta1,
         method1_stat2 = (beta2^2)/varBeta2,
         method2_stat_alt = ((abs(beta1) + abs(beta2))^2)/(varBeta1 + varBeta2 +
                      (1 + 1/r)*(1/(K_Treatment*m))*(VIF12*sqrt(varY1*varY2))),
         method2_stat = (betaC^2)/(2*(varYc/(K_Treatment*m))*(1 + (m-1)*rho0c)),
         method3_stat = ((sqrt((beta1^2)/(2*(varY1/(K_Treatment*m))*VIF1)) +
                        sqrt((beta2^2)/(2*(varY2/(K_Treatment*m))*VIF2)))^2)/
                        (2*(1 + VIF12/(sqrt(VIF1*VIF2)))),
         method4_stat = (K_Treatment*m*((beta1^2)*varY2*VIF2 -
                        2*beta1*beta2*VIF12*sqrt(varY1)*sqrt(varY2) +
                        (beta2^2)*varY1*VIF1))/(2*2*varY1*varY2*(VIF1*VIF2 - VIF12^2)),
         method5_stat1 = (beta1*sqrt(2*K_Treatment))/sqrt((4*varY1*VIF1)/m),
         method5_stat2 = (beta2*sqrt(2*K_Treatment))/sqrt((4*varY2*VIF2)/m)
         ) %>%
  mutate(method1_p1 = pchisq(q = method1_stat1, df = 1, lower.tail = FALSE),
         method1_p2 = pchisq(q = method1_stat2, df = 1, lower.tail = FALSE),
         method2_p_alt = pchisq(q = method2_stat_alt, df = 1, lower.tail = FALSE),
         method2_p = pchisq(q = method2_stat, df = 1, lower.tail = FALSE),
         method3_p = pchisq(q = method3_stat, df = 1, lower.tail = FALSE),
         method4_p_F = pf(method4_stat, df1 = 2, df2 = K_Total - 2*2,
                          lower.tail = FALSE, log.p = FALSE),
         method4_p_Chi2 = pchisq(q = method4_stat, df = 2, lower.tail = FALSE),
         method5_p1_Norm = pnorm(q = method5_stat1),
         method5_p2_Norm = pnorm(q = method5_stat2),
         method5_p1_T = pt(q = method5_stat1, df = (K_Total - 2*2)),
         method5_p2_T = pt(q = method5_stat2, df = (K_Total - 2*2)),
         abs_beta1 = abs(beta1),
         abs_beta2 = abs(beta2)
         ) %>%
  mutate_if(is.numeric, round, digits = 4) %>%
  mutate(method1_bonf_reject = ifelse(method1_p1 <= alpha_bonf| method1_p2 <= alpha_bonf,
                                      1, 0),
         method1_sidak_reject = ifelse(method1_p1 <= alpha_sidak | method1_p2 <= alpha_sidak,
                                       1, 0),
         method1_dap_reject = ifelse(method1_p1 <= alpha_dap | method1_p2 <= alpha_dap,
                                     1, 0),
         method2_reject = ifelse(method2_p <= 0.05,
                                 1, 0),
         method3_reject = ifelse(method3_p <= 0.05,
                                 1, 0),
         method4_F_reject = ifelse(method4_p_F <= 0.05,
                                   1, 0),
         method4_Chi2_reject = ifelse(method4_p_Chi2 <= 0.05,
                                      1, 0),
         method5_Norm_reject = ifelse(method5_p1_Norm <= 0.05 & method5_p1_Norm <= 0.05,
                                      1, 0),
         method5_T_reject = ifelse(method5_p1_T <= 0.05 & method5_p1_T <= 0.05,
                                   1, 0))

typeIerror_shortened <- typeIerror_data %>%
  mutate(betaC_true = abs(beta1_true) + abs(beta2_true)) %>%
  dplyr::select(Scenario_ID, Dataset_ID, K_Total, K_Treatment, m,
                beta1, beta1_true, beta2, beta2_true,
                betaC, betaC_true, contains("reject"))

typeIerror_summary <- typeIerror_shortened %>%
  group_by(Scenario_ID) %>%
  summarize(Scenario_n = n(),
            K_Total = mean(K_Total),
            K_Treatment = mean(K_Treatment),
            m = mean(m),
            beta1_mean = mean(beta1),
            beta2_mean = mean(beta2),
            betaC_mean = mean(betaC),
            beta1_true = mean(beta1_true),
            beta2_true = mean(beta2_true),
            betaC_true = mean(betaC_true),
            method1_bonf_reject = sum(method1_bonf_reject),
            method1_sidak_reject = sum(method1_sidak_reject),
            method1_dap_reject = sum(method1_dap_reject),
            method2_reject = sum(method2_reject),
            method3_reject = sum(method3_reject),
            method4_F_reject = sum(method4_F_reject),
            method4_Chi2_reject = sum(method4_Chi2_reject),
            method5_Norm_reject = sum(method5_Norm_reject),
            method5_T_reject = sum(method5_T_reject)
            ) %>%
  mutate(method1_null_case = ifelse(beta1_true == 0 & beta2_true == 0, 1, 0),
         method2_null_case = ifelse(betaC_true == 0, 1, 0),
         method3_null_case = ifelse(beta1_true == 0 & beta2_true == 0, 1, 0),
         method4_null_case = ifelse(beta1_true == 0 & beta2_true == 0, 1, 0),
         method5_null_case = ifelse(beta1_true != 0 & beta2_true != 0, 0, 1)) %>%
  mutate(method1_bonf_typeI_error_rate = ifelse(method1_null_case == 1,
                                                method1_bonf_reject/Scenario_n, NA),
         method1_sidak_typeI_error_rate = ifelse(method1_null_case == 1,
                                                 method1_sidak_reject/Scenario_n, NA),
         method1_dap_typeI_error_rate = ifelse(method1_null_case == 1,
                                               method1_dap_reject/Scenario_n, NA),
         method2_typeI_error_rate = ifelse(method2_null_case == 1,
                                           method2_reject/Scenario_n, NA),
         method3_typeI_error_rate = ifelse(method3_null_case == 1,
                                           method3_reject/Scenario_n, NA),
         method4_F_typeI_error_rate = ifelse(method4_null_case == 1,
                                             method4_F_reject/Scenario_n, NA),
         method4_Chi2_typeI_error_rate = ifelse(method4_null_case == 1,
                                                method4_Chi2_reject/Scenario_n, NA),
         method5_Norm_typeI_error_rate = ifelse(method5_null_case == 1,
                                           method5_Norm_reject/Scenario_n, NA),
         method5_T_typeI_error_rate = ifelse(method5_null_case == 1,
                                           method5_T_reject/Scenario_n, NA)) %>%
  filter(Scenario_ID != 0)

typeIerrorResultTable <- typeIerror_summary %>%
  dplyr::select(-m) %>%
  left_join(simParameters_null, by = c("Scenario_ID" = "Scenario")) %>%
  mutate(betas = paste0("(", beta1, ", ", beta2, ")"),
         rho0s = paste0("(", rho01, ", ", rho02, ")"),
         varYs = paste0("(", varY1, ", ", varY2, ")")) %>%
  dplyr::select(Scenario_ID,
                K, m, betas,
                rho0s, rho1, rho2, varYs,
                method1_bonf_typeI_error_rate,
                method1_sidak_typeI_error_rate,
                method1_dap_typeI_error_rate,
                method2_typeI_error_rate,
                method3_typeI_error_rate,
                method4_F_typeI_error_rate,
                method4_Chi2_typeI_error_rate,
                method5_Norm_typeI_error_rate,
                method5_T_typeI_error_rate) %>%
  mutate(method1_bonf_typeI_error_rate = round(method1_bonf_typeI_error_rate, digits = 3),
         method1_sidak_typeI_error_rate = round(method1_sidak_typeI_error_rate, digits = 3),
         method1_dap_typeI_error_rate = round(method1_dap_typeI_error_rate, digits = 3),
         method2_typeI_error_rate = round(method2_typeI_error_rate, digits = 3),
         method3_typeI_error_rate = round(method3_typeI_error_rate, digits = 3),
         method4_F_typeI_error_rate = round(method4_F_typeI_error_rate, digits = 3),
         method4_Chi2_typeI_error_rate = round(method4_Chi2_typeI_error_rate, digits = 3),
         method5_T_typeI_error_rate = round(method5_T_typeI_error_rate, digits = 3),
         method5_Norm_typeI_error_rate = round(method5_Norm_typeI_error_rate, digits = 3)) %>%
  dplyr::rename('1a. P-Value Bonf.' = method1_bonf_typeI_error_rate,
                '1b. P-Value Sidak' = method1_sidak_typeI_error_rate,
                '1c. P-Value DAP' = method1_dap_typeI_error_rate,
                '2. Combined' = method2_typeI_error_rate,
                '3. Single Weighted' = method3_typeI_error_rate,
                '4a. Disjunctive F' = method4_F_typeI_error_rate,
                '4b. Disjunctive Chi2' = method4_Chi2_typeI_error_rate,
                '5a. Conjunctive T' = method5_T_typeI_error_rate,
                '5b. Conjunctive MVN' = method5_Norm_typeI_error_rate) %>%
  arrange(betas, rho0s)

# Save Result Tables -----------------------------------------------------------

# Power result table
#write.csv(powerResultTable, paste0("./ResultTables/powerResultTable_0.csv"))
#write.csv(powerResultTable, paste0("./ResultTables/powerResultTable_1.csv"))

# Type I error result table
#write.csv(typeIerrorResultTable, paste0("./ResultTables/typeIerrorResultTable_0.csv"))
#write.csv(typeIerrorResultTable, paste0("./ResultTables/typeIerrorResultTable_1.csv"))
