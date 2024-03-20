source("/Users/melodyowen/Desktop/GitHub/hybrid2simulation/DataGeneration/Generate_Data_2_Cont_Outcomes.R")

# Generate Continuous Datasets
simSetsCont <- create_all_cont_sim_dats(n = 10,
                                        K_input = 10,
                                        m_input = 500,
                                        beta1_input = 0.2,
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
  mutate(`Absolute Difference` = abs(`Input Value` - Value))

# Generate Binary Datasets



