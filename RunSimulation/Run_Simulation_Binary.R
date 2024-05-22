source("./DataGeneration/Generate_Data_2_Bin_Outcomes.R")
source("./RequiredPackages.R")

# Generate Binary Datasets -----------------------------------------------------
simSetsBin <- create_all_bin_sim_dats(n = 1,
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

combinedSimStatsBin <- bind_rows(simSetsBin[[2]], .id = "Index") %>%
  mutate(Index = as.integer(Index)) %>%
  mutate(`Absolute Difference` = abs(`True Value` - `Estimated Value`))
