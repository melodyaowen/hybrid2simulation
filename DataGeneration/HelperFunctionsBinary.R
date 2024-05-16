# Helper functions for binary data generation

# Function to construct correlation matrix for Y_i -----------------------------
constrRi <- function(rho01, rho2, m, K){
  rhoMatrix <- as.numeric()
  for(row in 1:K){
    rhoMatrix_row <- as.numeric()
    for(col in 1:K){
      corr1 <- rho01[row,col]
      corr2 <- rho2[row,col]
      corrBlock <- (corr2 - corr1)*diag(m) + matrix(corr1, nrow = m, ncol = m)
      rhoMatrix_row  <- cbind(rhoMatrix_row,corrBlock)
    }
    rhoMatrix <- rbind(rhoMatrix,rhoMatrix_row)
  }
  rhoMatrix <- as.matrix(rhoMatrix)
  return(rhoMatrix)
}

# Function to calculate matrix Gamma_k -----------------------------------------
calGammaK <- function(r, m, beta1k, beta2k){
  muC <- exp(beta1k)/(exp(beta1k) + 1)
  muT <- exp(beta1k + beta2k)/(exp(beta1k + beta2k) + 1)
  piece_C <- (1 - r)*muC*(1 - muC)
  piece_T <- r*muT*(1 - muT)
  output <- m* matrix(c((piece_C + piece_T), piece_T, piece_T, piece_T),
                      nrow = 2, ncol = 2)
  return(output)
}

# Function to calculate matrix Omega_k -----------------------------------------
calOmegak <-  function(r, m, beta1k, beta2k, rho0k){
  cstTerm <- (m + m*(m - 1)*rho0k)
  muC <- exp(beta1k)/(exp(beta1k) + 1)
  muT <- exp(beta1k + beta2k)/(exp(beta1k + beta2k) + 1)
  piece_C <- (1 - r)*muC*(1 - muC)
  piece_T <- r*muT*(1 - muT)
  output <- cstTerm* matrix(c((piece_C + piece_T), piece_T, piece_T, piece_T),
                            nrow = 2, ncol = 2)
  return(output)
}

# Function to calculate covariance matrix SigmaK -------------------------------
calSigmaK <- function(r, m, beta1k, beta2k, rho0k){
  Gammak <- calGammaK(r, m, beta1k, beta2k)
  Omegak <- calOmegak(r, m, beta1k, beta2k, rho0k)
  return(solve(Gammak) %*% Omegak %*% solve(Gammak))
}

calsigma2ksq <- function(r, m, beta1k, beta2k, rho0k){
  cstTerm <- 1 + (m - 1)*rho0k
  muC <- exp(beta1k)/(exp(beta1k) + 1)
  muT <- exp(beta1k + beta2k)/(exp(beta1k + beta2k) + 1)
  piece_C <- cstTerm/m/(1 - r)/(muC*(1 - muC))
  piece_T <- cstTerm/m/r/(muT*(1 - muT))
  return(piece_C + piece_T)
}

# Function to calculate Omega_k1k2 ---------------------------------------------
calOmegak1k2 <- function(r, m, beta1k1, beta2k1, beta1k2,
                         beta2k2, rho1k1k2, rho2k1k2){
  muC1 <-  exp(beta1k1)/(exp(beta1k1) + 1)
  muC2 <- exp(beta1k2)/(exp(beta1k2) + 1)
  muT1 <-  exp(beta1k1 + beta2k1)/(exp(beta1k1 + beta2k1) + 1)
  muT2 <- exp(beta1k2 + beta2k2)/(exp(beta1k2 + beta2k2) + 1)
  xiC <- sqrt(muC1*muC2*(1 - muC1)*(1 - muC2))
  xiT <-  sqrt(muT1*muT2*(1 - muT1)*(1 - muT2))
  cstTerm <- m*rho2k1k2 + m*(m - 1)*rho1k1k2
  piece_C <- (1 - r)*xiC
  piece_T <- r*xiT
  output <- cstTerm*matrix(c((piece_C + piece_T), piece_T, piece_T, piece_T),
                           nrow = 2,ncol = 2)
  return(output)
}

# Function to calculate covariance matrix between beta_k1 and beta_k2 ----------
calCov2Betas <- function(r, m, beta1k1, beta2k1, beta1k2,
                         beta2k2, rho1k1k2, rho2k1k2){
  omegak1k2 <- calOmegak1k2(r, m, beta1k1, beta2k1, beta1k2,
                            beta2k2, rho1k1k2, rho2k1k2)
  gammak1 <- calGammaK(r, m, beta1k1, beta2k1)
  gammak2 <- calGammaK(r, m, beta1k2, beta2k2)
  output <- solve(gammak1) %*%  omegak1k2  %*% solve(gammak2)
  return(output)
}


calsigma2k1k2sq <- function(r, m, beta1k1, beta2k1, beta1k2,
                            beta2k2, rho1k1k2, rho2k1k2){
  cstTerm <- rho2k1k2 + (m - 1)*rho1k1k2
  muC1 <-  exp(beta1k1)/(exp(beta1k1) + 1)
  muC2 <- exp(beta1k2)/(exp(beta1k2) + 1)
  muT1 <-  exp(beta1k1 + beta2k1)/(exp(beta1k1 + beta2k1) + 1)
  muT2 <- exp(beta1k2 + beta2k2)/(exp(beta1k2 + beta2k2) + 1)
  xiC <- sqrt(muC1*muC2*(1 - muC1)*(1 - muC2))
  xiT <-  sqrt(muT1*muT2*(1 - muT1)*(1 - muT2))
  pieceC <- cstTerm/m/(1 - r)/xiC
  pieceT <- cstTerm/m/(r)/xiT
  return(pieceC + pieceT)
}

# Function to calculate correlation between test statistics wk1 and wk2 --------
calCorWk1Wk2 <- function(r, m, beta1k1, beta2k1, beta1k2, beta2k2,
                         rho0k1, rho0k2, rho1k1k2, rho2k1k2){
  top <- calsigma2k1k2sq(r, m, beta1k1, beta2k1, beta1k2,
                         beta2k2, rho1k1k2, rho2k1k2)
  sigmak1 <- sqrt(calsigma2ksq(r, m, beta1k1, beta2k1, rho0k1))
  sigmak2 <- sqrt(calsigma2ksq(r, m, beta1k2, beta2k2, rho0k2))
  result <- top/(sigmak1*sigmak2)
  return(result)
}


