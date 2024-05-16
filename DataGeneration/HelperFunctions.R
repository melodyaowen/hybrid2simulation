# Calculate the total variance from the within-cluster variance and ICC
# rho0 is the outcome specific ICC
# p is the proportion of individiausl in a single cluster with Y = 1 (binary)
calTotalVarianceBinary <- function(rho0, p){
  sigma2W <- p*(1-p) # within-cluster component of variance
  sigma2B <- (rho0*sigma2W)/(1-rho0) # between-cluster component of variance
  totalVar <- sigma2W + sigma2B
  return(totalVar)
}
