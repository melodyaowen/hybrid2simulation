source("./DataGeneration/Generate_Data_2_Cont_Outcomes.R")
source("./RequiredPackages.R")

# Table of all Parameters ------------------------------------------------------



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





