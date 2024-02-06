myData <- as.data.frame(simSets[[1]][1])

# Returns design speficiations based on simulated dataset
calculate_data_specifications <- function(myData){

  myData <- mySimData # K = 6,
  # m = 100,
  # beta1 = 0.2,
  # beta2 = 0.1,
  # rho01 = 0.07,
  # rho02 = 0.05,
  # rho1 = 0.02,
  # rho2 = 0.25,
  # varY1 = 0.3,
  # varY2 = 0.2,
  # r = 1

  # Calculate outcome specific ICC's from simulated data
  m1 <- geepack::geeglm(formula = Y1 ~ treatment_z_k, # creating a geeglm model
                        id = factor(cluster_k),
                        family = gaussian(link = "identity"),
                        data = myData, corstr = "exchangeable")
  summary_m1 <- summary(m1) # summary output for m1
  rho01_sim <- summary_m1$corr[1,1]
  m2 <- geepack::geeglm(formula = Y2 ~ treatment_z_k, # creating a geeglm model
                        id = factor(cluster_k),
                        family = gaussian(link = "identity"),
                        data = myData, corstr = "exchangeable")
  summary_m2 <- summary(m2) # summary output for m2
  rho02_sim <- summary_m2$corr[1,1]

  # Calculate variances from simulated data
  varY1_sim <- var(myData$Y1)
  varY2_sim <- var(myData$Y2)

  # Calculating rho1 and rho2 values from simulated data
  fm1 <- nlme::lme(Y1 ~ treatment_z_k, random = ~1|cluster_k, data = myData)
  fm2 <- nlme::lme(Y2 ~ treatment_z_k, random = ~1|cluster_k, data = myData)

  zeta_sim <- as.numeric(c(fm1$coefficients$fixed, fm2$coefficients$fixed))
  beta1_sim <- zeta_sim[1:2]
  beta2_sim <- zeta_sim[3:4]
}
