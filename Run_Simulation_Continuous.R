source("./DataGeneration/Generate_Data_2_Cont_Outcomes.R")
source("./RequiredPackages.R")

# Table of all Parameters ------------------------------------------------------
sim_betas <- tibble(betas = c(paste("0.3 0.7"),
                              paste("0.7 0.7")))
sim_K <- tibble(K = c(15))
sim_m <- tibble(m = c(80))
sim_rho0 <- tibble(rho0 = c(paste("0.01 0.1"),
                            paste("0.1 0.1"),
                            paste("0.1 0.01")))
sim_rho1 <- tibble(rho1 = c(0.025))
sim_rho2 <- tibble(rho2 = c(0.2, 0.5))

simParameters <- expand.grid(K = c(15),
                             m = c(80),
                             betas = c(paste("0.3 0.7"),
                                       paste("0.7 0.7")),
                             rho0 = c(paste("0.01 0.1"),
                                      paste("0.1 0.1"),
                                      paste("0.1 0.01")),
                             rho1 = c(0.025),
                             rho2 = c(0.2, 0.5)
                             ) %>%
  separate(betas, sep = " ", into = c("beta1", "beta2")) %>%
  separate(rho0, sep = " ", into = c("rho01", "rho02"))

# Run power calculations on all true parameters --------------------------------


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





