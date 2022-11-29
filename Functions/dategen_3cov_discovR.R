# data generation 3 components for discovR #
# 29-sep-2021 #
# soogeun park #

# this is a quick adaptation of the data generation function used for the logistic regression paper
# but now this is for discovR 

data_3cov <- function(I, J, Jy, d1_active, c_active_block1, c_active_block2, 
                    d2_active, comp_sd, signal_level,
                    Py, seed, strict_uncorrelated){
  
  # Inputs #
  # I: number of observations
  # J: number of predictor variables
  # Jy: number of outcome variables
  # d1_active: indices of non-zero loadings for the covariate distinctive to the first block
  # c_active_block1, c_active_block2: indices of non-zero loadings for the common covariate
  # d2_active: indices of non-zero loadings for the covariate distinctive to the second block
  # comp_sd: magnitude of standard deviation defined for the covariates
  # signal_level: the proportion of variance accounted for by the covariates
  # Py: regression coefficients (will be normalized)
  # seed: random seed value
  # strict_uncorrelated: if TRUE, the covariates are orthonormalized. This can only work for low-dimensionality.
  # thus, it is always FALSE in the paper. The covariates are anyhow generated from multivariate normal distribution,
  # with a diagonal covariance matrix.
  
  library(MASS)
  
  I <- I
  J1 <- J/2
  J2 <- J/2
  
  set.seed(seed)
  
  # generating the covariates and noise variables, such that they are all uncorrelated
  H_noise <- MASS::mvrnorm(n = I, mu = rep(0, 3+J+Jy), Sigma = diag(c(comp_sd, rep(comp_sd[1], J), rep(comp_sd[1], Jy))^2))
  
  H <- H_noise[,1:3]
  
  # if TRUE, the covariates are orthonormalized (usually FALSE)
  if (strict_uncorrelated){
    H <- scale(H,T,F)
    
    H <- qr.Q(qr(H))
    
    H <- H * sqrt(comp_sd[1]^2 * I)
  }
  
  X_noise <- H_noise[,4:(J+3)]
  
  Y_noise <- H_noise[,(J+3+1):(3+J+Jy)]
  
  # defining the loadings matrix Px
  truePx <- matrix(0, nrow = J, ncol = 3)
  
  truePx[d1_active,1] <- 1
  truePx[c_active_block1,2] <- 1
  truePx[c_active_block2,2] <- 1
  truePx[d2_active,3] <- 1
  
  # loadings orthonormalized
  truePx <- qr.Q(qr(truePx))
  
  # Xtrue defined
  Xtrue <- H %*% t(truePx)
  
  # amount of SSE (sum of squared errors) needed
  SS_noise <- (1 - signal_level) / signal_level * sum(Xtrue^2)
  
  # normalize the noise vectors to unit length
  X_noise_normalized <- X_noise %*% diag(1/sqrt(diag(t(X_noise) %*% X_noise)))
  
  # correctly scaling the X_noise to achieve the desired signal_level
  X_noise_scaled <- X_noise_normalized * sqrt(SS_noise/J)
  
  sum(X_noise_scaled^2)
  
  # sum(Xtrue^2) / (sum(Xtrue^2) + sum(X_noise_scaled^2))
  
  # observed X defined
  X_train <- Xtrue + X_noise_scaled
  
  # signal level check
  X_signal_level_train <- sum(Xtrue^2) / sum(X_train^2)
  
  # regression coefficients scaled
  dividing <- diag(1/sqrt(diag(t(Py) %*% Py)))
  
  if (sum(colSums(Py == 0) == nrow(Py)) > 0){
    
    zerocolumn <- colSums(Py == 0) == nrow(Py)
    
    dividing <- diag(1/sqrt(diag(t(Py) %*% Py)))
    
    dividing[zerocolumn,zerocolumn] <- 1
    
  }
  
  Py <- Py %*% dividing 
  
  Ytrue <- H %*% t(Py)
  
  # amount of SSE (sum of squared errors) needed for Y (this is the same as the level for X)
  SS_noise <- (1 - signal_level) / signal_level * sum(Ytrue^2)
  
  # normalize the noise vectors to unit length
  Y_noise_normalized <- Y_noise %*% diag(1/sqrt(diag(t(Y_noise) %*% Y_noise)))
  
  # correctly scaling the Y_noise to achieve the desired signal_level
  Y_noise_scaled <- Y_noise_normalized * sqrt(SS_noise/Jy)
  
  sum(Y_noise_scaled^2)
  
  # observed Y defined
  Y_train <- Ytrue + Y_noise_scaled
  
  Y_signal_level_train <- sum(Ytrue^2) / sum(Y_train^2)
  
  
  # test set: exactly the same procedure as the train set #
  
  # generating the covariates and noise variables, such that they are all uncorrelated
  H_noise <- MASS::mvrnorm(n = I, mu = rep(0, 3+J+Jy), Sigma = diag(c(comp_sd, rep(comp_sd[1], J), rep(comp_sd[1], Jy))^2))
  
  H <- H_noise[,1:3]
  
  # if TRUE, the covariates are orthonormalized (usually FALSE)
  if (strict_uncorrelated){
    H <- scale(H,T,F)
    
    H <- qr.Q(qr(H))
    
    H <- H * sqrt(comp_sd[1]^2 * I)
  }
  
  X_noise <- H_noise[,4:(J+3)]
  
  Y_noise <- H_noise[,(J+3+1):(3+J+Jy)]
  
  # defining the loadings matrix Px
  truePx <- matrix(0, nrow = J, ncol = 3)
  
  truePx[d1_active,1] <- 1
  truePx[c_active_block1,2] <- 1
  truePx[c_active_block2,2] <- 1
  truePx[d2_active,3] <- 1
  
  # loadings orthonormalized
  truePx <- qr.Q(qr(truePx))
  
  # Xtrue defined
  Xtrue <- H %*% t(truePx)
  
  # amount of SSE (sum of squared errors) needed
  SS_noise <- (1 - signal_level) / signal_level * sum(Xtrue^2)
  
  # normalize the noise vectors to unit length
  X_noise_normalized <- X_noise %*% diag(1/sqrt(diag(t(X_noise) %*% X_noise)))
  
  # correctly scaling the X_noise to achieve the desired signal_level
  X_noise_scaled <- X_noise_normalized * sqrt(SS_noise/J)
  
  sum(X_noise_scaled^2)
  
  # sum(Xtrue^2) / (sum(Xtrue^2) + sum(X_noise_scaled^2))
  
  # observed X defined
  X_test <- Xtrue + X_noise_scaled
  
  # signal level check
  X_signal_level_test <- sum(Xtrue^2) / sum(X_test^2)
  
  # regression coefficients scaled
  dividing <- diag(1/sqrt(diag(t(Py) %*% Py)))
  
  if (sum(colSums(Py == 0) == nrow(Py)) > 0){
    
    zerocolumn <- colSums(Py == 0) == nrow(Py)
    
    dividing <- diag(1/sqrt(diag(t(Py) %*% Py)))
    
    dividing[zerocolumn,zerocolumn] <- 1
    
  }
  
  Py <- Py %*% dividing 
  
  Ytrue <- H %*% t(Py)
  
  # amount of SSE (sum of squared errors) needed for Y (this is the same as the level for X)
  SS_noise <- (1 - signal_level) / signal_level * sum(Ytrue^2)
  
  # normalize the noise vectors to unit length
  Y_noise_normalized <- Y_noise %*% diag(1/sqrt(diag(t(Y_noise) %*% Y_noise)))
  
  # correctly scaling the Y_noise to achieve the desired signal_level
  Y_noise_scaled <- Y_noise_normalized * sqrt(SS_noise/Jy)
  
  sum(Y_noise_scaled^2)
  
  # observed Y defined
  Y_test <- Ytrue + Y_noise_scaled
  
  Y_signal_level_test <- sum(Ytrue^2) / sum(Y_test^2)
  
  
  
  
  result <- list(X_train = X_train, Y_train = Y_train, X_test = X_test, Y_test = Y_test, truePx = truePx, 
                 Py = Py, X_signal_level_test = X_signal_level_test, X_signal_level_train = X_signal_level_train,
                 Y_signal_level_test = Y_signal_level_test, Y_signal_level_train = Y_signal_level_train)
  
  return (result)
  
}