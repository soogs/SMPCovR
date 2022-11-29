# discovR simulation study #
# 24th November 2022 #
# Soogeun Park #

# rep 1 #

# 1. functions and libraries ####

library(mixOmics)
library(RegularizedSCA)
library(Rcpp)
library(nFactors)

setwd("E:\\Users\\park\\Desktop\\discovR simulation\\to_blade_20_june_2022\\to_blade_20_june_2022\\")

source("./discovR_cpp_ridge.R")
source("./dategen_3cov_discovR.R")
source("./sPLS_cv.R")

sourceCpp("./updateW_discovR_ridge.cpp")

setwd("E:\\Users\\park\\Desktop\\discovR simulation\\simulation 24_nov_2022\\")

source("./conditions_making_24_nov_2022.R")
source("./discovR_cv_rsq_ridge 9-nov-2022.R")


permutations <- function(n){
  if(n==1){
    return(matrix(1))
  } else {
    sp <- permutations(n-1)
    p <- nrow(sp)
    A <- matrix(nrow=n*p,ncol=n)
    for(i in 1:n){
      A[(i-1)*p+1:p,] <- cbind(i,sp+(sp>=i))
    }
    return(A)
  }
}

# congruence for 2 vectors
congruence <- function(a, b){
  
  if(anyNA(cbind(a, b))){
    return (0)
  }
  
  if (sum(abs(a)) < 1e-10 | sum(abs(b)) < 1e-10){
    result <- 0
    
    return(result)
  }
  
  result <- (t(a) %*% b) / sqrt(sum(a^2) * sum(b^2))
  return(c(result))
}


# evaluation criteria functions #
corrects <- function(estimate, defined){
  
  if (is.matrix(estimate)){
    total <- prod(dim(estimate))
    
    estimate_zero <- sum(colSums(estimate == 0) == nrow(estimate))
    defined_zero <- sum(colSums(defined == 0) == nrow(defined))
    
    if (estimate_zero > defined_zero){
      
      tucker <- RegularizedSCA::TuckerCoef(defined, estimate)
      
      estimate <- estimate[,tucker$perm]
      
    } else {
      
      tucker <- RegularizedSCA::TuckerCoef(estimate, defined)
      
      defined <- defined[,tucker$perm]
      
      
    }
    
    
    ratio <- (sum((abs(estimate) > 1e-7) + (abs(defined) > 1e-7) == 2) + 
                sum((abs(estimate) < 1e-7) + (abs(defined) < 1e-7) == 2)) / (total)
    
  } else {
    total <- length(estimate)
    
    ratio <- (sum((abs(estimate) > 1e-7) + (abs(defined) > 1e-7) == 2) + 
                sum((abs(estimate) < 1e-7) + (abs(defined) < 1e-7) == 2)) / (total)
    
  }
  
  return(ratio)
}


# corrects for a vector

corrects_vec <- function(estimate, defined){
  
  ratio <- (sum((abs(estimate) > 1e-7) + (abs(defined) > 1e-7) == 2) + 
              sum((abs(estimate) < 1e-7) + (abs(defined) < 1e-7) == 2)) / (length(estimate))
  
  return(ratio)
  
}

nonzero <- function(estimate, defined){
  
  if (is.matrix(estimate)){
    nonzeros <-  sum(abs(defined) > 1e-7)
    
    tucker <- RegularizedSCA::TuckerCoef(defined, estimate)
    
    estimate <- estimate[,tucker$perm]
    
    ratio <- sum((abs(estimate) > 1e-7) + 
                   (abs(defined) > 1e-7) == 2) / (nonzeros)
    
  } else {
    nonzeros <-  sum(abs(defined) > 1e-7)
    
    ratio <- sum((abs(estimate) > 1e-7) + 
                   (abs(defined) > 1e-7) == 2) / (nonzeros)
    
  }
  
  return (ratio)
}


# common and distinctive identifier #
comdis <- function(estimate){
  J <- nrow(estimate)
  
  d1zero <- colSums(estimate[1:(J/2),]) == 0
  d2zero <- colSums(estimate[(J/2 + 1):J,]) == 0
  
  commz <- (d1zero + d2zero) == 0
  
  result <- list(d1 = sum(d2zero), d2 = sum(d1zero), common = sum(commz))
  
  return(result)
}

# scale_zero function for standardizing a matrix which contains zero columns
scale_zero <- function(x, center, standardize){
  zerocols <- which(colSums(x == 0) == nrow(x))
  nonzerocols <- which(colSums(x == 0) != nrow(x))
  
  result <- x
  
  result[,nonzerocols] <- scale(result[,nonzerocols], center, standardize)
  
  return(result)
  
}


# Model selection process #
# Model selection according to the paper. First alpha and ridge CV with fixed lasso and glasso values
# -> and then one more run of discovR (no CV) to find lasso and glasso values that return the correct weights.

results <- data.frame(matrix(NA, nrow = nrow(condition_df), ncol = 55))

results[,1] <- factor(results[,1], levels = c("low", "high"))

colnames(results) <- c("dimensions", "signal_level", "Jy", "py_pattern", "reps",
                       "model_seed", "noise_seed",
                       
                       "disc_alpha", "disc_ridge_y", "disc_lasso_y",
                       "disc_lasso", "disc_ridge", 
                       "disc_fit", "disc_pred", 
                       "disc_pred0", "disc_pred_nonzero",
                       "disc_tucker",
                       "disc_correct", 
                       "disc_nonzero", 
                       "disc_py_correct",
                       "disc_zeros1", "disc_zeros2", "disc_zeros3",
                       "disc_Y_included", "disc_Y_correct",
                       "disc_CV_Y_included", "disc_CV_Y_correct",
                       
                       "lasso0_alpha", "lasso0_ridge_y", "lasso0_lasso_y",
                       "lasso0_lasso", "lasso0_ridge", 
                       "lasso0_fit", "lasso0_pred", 
                       "lasso0_pred0", "lasso0_pred_nonzero",
                       "lasso0_tucker",
                       "lasso0_correct", 
                       "lasso0_nonzero", 
                       "lasso0_py_correct",
                       "lasso0_zeros1", "lasso0_zeros2", "lasso0_zeros3",
                       "lasso0_zerocoefs", "lasso0_zerocoefs_firstfive",
                       
                       "diacon_fit", "diacon_pred",
                       "diacon_pred0", "diacon_pred_nonzero",
                       "diacon_tucker",
                       "diacon_correct", 
                       "diacon_nonzero", 
                       "diacon_py_correct",
                       "diacon_zerocoefs", "diacon_zerocoefs_firstfive")


# 2. replication starts ####
set.seed(1111) 
model_seed <- sample(x = 1:100000, size = nrow(condition_df))
noise_seed <- sample(x = 1:100000, size = nrow(condition_df))

dim(condition_df)

unique(condition_df$signal_level)

# dividing the task such that each run takes equal number of high dimensionality problem #
repshigh <- rownames(subset(condition_df, dimension == "high"))[1:50]
repslow <- rownames(subset(condition_df, dimension == "low"))[1:50]

reps_to_do <- c(1, as.numeric(repshigh), as.numeric(repslow))

length(reps_to_do)

for (rrr in reps_to_do){ 
  
  modelseeding <- model_seed[rrr]
  noiseseeding <- noise_seed[rrr]
  
  cond <- condition_df[rrr,]
  
  
  if (cond$dimension == "low"){
    I <- 100
    J <- 30
    
    d1_active <- 1:4
    c_active_block1 <- 5:8
    c_active_block2 <- 16:19
    d2_active <- 20:23
  }
  
  if (cond$dimension == "high"){
    I <- 100
    J <- 200
    
    d1_active <- 1:25
    c_active_block1 <- 26:45
    c_active_block2 <- 101:120
    d2_active <- 121:145
    
  }
  
  if (cond$Jy == 5){
    
    if(cond$py_pattern == "all3"){
      
      Py <- matrix(0, nrow = 5, ncol = 3)
      
      Py[1:2,1] <- 1
      Py[3:4,2] <- 1
      Py[5,3] <- 1
      
      Py[2,] <- 1
      
    } 
    
    
    if(cond$py_pattern == "notall3"){
      
      Py <- matrix(0, nrow = 5, ncol = 3)
      
      Py[1,1] <- 1
      Py[2,2] <- 1
      Py[3,3] <- 1
      
      Py[3,] <- 1
    }
    
    if(cond$py_pattern == "all2"){
      
      Py <- matrix(0, nrow = 5, ncol = 3)
      
      Py[1:3,1] <- 1
      Py[4:5,2] <- 1
      
      Py[3,1:2] <- 1
      
    }
    
    
    if(cond$py_pattern == "notall2"){
      
      Py <- matrix(0, nrow = 5, ncol = 3)
      
      Py[1:2,1] <- 1
      Py[2:3,2] <- 1
      
    }
    
    
  } else if (cond$Jy == 20){
    
    
    if(cond$py_pattern == "all3"){
      
      Py <- matrix(0, nrow = 20, ncol = 3)
      
      Py[1:7,1] <- 1
      Py[8:14,2] <- 1
      Py[15:20,3] <- 1
      
      Py[7,] <- 1
      Py[14,] <- 1
      Py[20,] <- 1
      
    } 
    
    
    if(cond$py_pattern == "notall3"){
      
      Py <- matrix(0, nrow = 20, ncol = 3)
      
      Py[1:4,1] <- 1
      Py[5:8,2] <- 1
      Py[9:12,3] <- 1
      
      Py[4,] <- 1
      Py[8,] <- 1
      Py[12,] <- 1
      
    }
    
    if(cond$py_pattern == "all2"){
      
      Py <- matrix(0, nrow = 20, ncol = 3)
      
      Py[1:10,1] <- 1
      Py[11:20,2] <- 1
      
      Py[9:12,1:2] <- 1
      
    }
    
    
    if(cond$py_pattern == "notall2"){
      
      Py <- matrix(0, nrow = 20, ncol = 3)
      
      Py[1:6,1] <- 1
      Py[7:12,2] <- 1
      
      Py[c(6,7),1:2] <- 1
    }
    
    
  }
  
  
  pyzeros <- which(rowSums(Py) == 0)
  
  
  cond$signal_level <- as.numeric(cond$signal_level)
  cond$reps <- as.numeric(cond$reps)
  cond$Jy <- as.numeric(cond$Jy)
  
  # data generating #
  dat <- data_3cov(I = I, J = J, Jy = cond$Jy, d1_active = d1_active, 
                   c_active_block1 = c_active_block1, c_active_block2 = c_active_block2,
                   d2_active = d2_active, comp_sd = c(50,50,50), 
                   signal_level = cond$signal_level,
                   strict_uncorrelated = F,
                   Py = Py, seed = noiseseeding)
  
  # scaling the data #
  X_train <- dat$X_train
  Y_train <- dat$Y_train
  
  X_train <- scale(X_train,T,T)
  Y_train <- scale(Y_train,T,T)
  
  X_test <- dat$X_test
  Y_test <- dat$Y_test
  
  X_test <- scale(X_test,T,T)
  Y_test <- scale(Y_test,T,T)
  
  
  truePx <- matrix(0, nrow = ncol(X_test), ncol = 3)
  
  truePx[d1_active,1] <- 1
  truePx[c_active_block1,2] <- 1
  truePx[c_active_block2,2] <- 1
  truePx[d2_active,3] <- 1
  
  irr <- which(abs(Py) < 1) 
  # index for the covariate that is irrelevant
  
  # splitting the truePx matrix for use for DIABLO #
  truePx_d1 <- truePx[1:(J/2), 1]
  truePx_c <- truePx[, 2]
  truePx_d2 <- truePx[(J/2 +1):J,3]
  
  # discovR - exhasutive CV ####
  
  R <- 3 # number of covariates
  
  # just in case, i can check for the number of covariates with acceleration factor
  svdd <- svd(cbind(Y_train, X_train))
  
  nScree(svdd$d)
  
  nScree(svdd$d^2)
  
  # exhaustive cross validation: alpha, lasso_y and lasso_w #
  time1 <- Sys.time()
  
  lasso_range <- c(exp(seq(log(0.00001), log(0.5), length.out = 7)))
  
  lasso_y_range <- c(exp(seq(log(0.00001), log(0.5), length.out = 7)))
  
  alpha_range <- seq(0.1, 0.9, 0.2)
  
  ranges <- expand.grid(alpha_range, lasso_y_range, lasso_range)
  
  colnames(ranges) <- c("alpha", "lasso_y", "lasso")
  
  ranges <- ranges[order(ranges$lasso, decreasing = T),]
  
  ranges <- ranges[order(ranges$lasso_y, decreasing = T),]
  
  ranges <- ranges[order(ranges$alpha, decreasing = F),]
  
  scd_cv_results <- c()
  scd_cv_conv <- c()
  
  disc_cv_results <- matrix(NA, ncol = 8+3, nrow = nrow(ranges))
  
  colnames(disc_cv_results) <- c("cve", "se", "conv", "Y_zero", 
                                 "Y_included_CV", "Y_correct_CV",
                                 "W0", "Py0",
                                 "alpha", "lasso_y", "lasso")
  
  disc_cv_results <- as.data.frame(disc_cv_results)
  
  
  colnames(ranges) <- c("alpha", "lasso_y", "lasso")
  
  
  for (i in 1:nrow(ranges)){
    
    disc_cv1 <- discovR_cv_rsq(X = X_train, Y = Y_train, 
                               R = R, 
                               alpha = ranges[i,]$alpha, 
                               lasso_w = rep(ranges[i,]$lasso,R), 
                               ridge_w = 1e-7, 
                               lasso_y = rep(ranges[i,]$lasso_y,3),
                               ridge_y = 1e-7, 
                               seed_cv = i+1, seed_method = 111,
                               inits = "rational", include_rational = TRUE, 
                               MAXITER = 500, nrFolds = 5, nonzero_only = TRUE, stop_value = 1e-3, 
                               true_active_Y = rowSums(abs(Py)) > 1e-5)
    
    # apply the model parameters on an entire dataset to study the sparsity level
    fit1 <- discovR_cpp_ridge(X = X_train, Y = Y_train, R = R, 
                              alpha = ranges[i,]$alpha, 
                              lasso_w =  rep(ranges[i,]$lasso,R), 
                              ridge_w = 1e-7,
                              lasso_y = rep(ranges[i,]$lasso_y,3), 
                              ridge_y = 1e-7, inits = "rational", 
                              nrstart = 1, include_rational = TRUE, seed = 11, MAXITER = 1000, 
                              stop_value = 1e-3)
    
    disc_cv_results[i,] <- c(disc_cv1$cve, disc_cv1$se, sum(disc_cv1$convergence), disc_cv1$zero_Y,
                             mean(disc_cv1$Y_included), mean(disc_cv1$Y_correct),
                             sum(fit1$W == 0), sum(fit1$Py == 0),
                             ranges[i,]$alpha, ranges[i,]$lasso_y, ranges[i,]$lasso)
    
    print(i)
    
  }
  
  time2 <- Sys.time()
  
  disc_cv_results <- disc_cv_results[!is.na(disc_cv_results$cve),]
  
  disc_cv_results$zerocoef <- disc_cv_results$W0 + disc_cv_results$Py0
  
  cv_max_index <- which.max(disc_cv_results[,1]) # largest R2
  
  serule <- disc_cv_results[cv_max_index,1] - disc_cv_results[cv_max_index,2] # 1SE rule
  
  disc_possible <- disc_cv_results[disc_cv_results[,1] >= serule,] # larger than the 1SE rule
  
  
  # new model selection (9 nov 2022) #
  # among the models with the lowest number of outcome varibles,
  # choose the most sparse model #
  
  # subsetting only the models with the lowest number of outcome variables
  disc_possible <- disc_possible[disc_possible$Y_zero == max(disc_possible$Y_zero),]
  
  if (is.vector(disc_possible)){
    disc_chosen <- disc_possible
  } else {
    
    disc_possible <- disc_possible[order(disc_possible$lasso, decreasing = T),]
    
    disc_possible <- disc_possible[order(disc_possible$lasso_y, decreasing = T),]
    
    disc_possible <- disc_possible[order(disc_possible$alpha, decreasing = F),]
    
    disc_possible <- disc_possible[order(disc_possible$W0, decreasing = T),]
    
    disc_possible <- disc_possible[order(disc_possible$Py0, decreasing = T),]
    
    disc_possible <- disc_possible[order(disc_possible$zerocoef, decreasing = T),]
    
    disc_chosen <- disc_possible[1,]
    
  }
  
  
  disc_alpha <- disc_chosen$alpha
  disc_lasso <- disc_chosen$lasso
  disc_lasso_y <- disc_chosen$lasso_y
  
  # since we cross-validated already with alpha, py-ridge and w-lasso,
  # we leave the solution as it is!
  
  # final discovR model fitting #
  disc1 <- discovR_cpp_ridge(X = X_train, Y = Y_train, R = R, 
                             alpha = disc_alpha, 
                             lasso_w = rep(disc_lasso,3), 
                             ridge_w = 1e-7,
                             lasso_y = rep(disc_lasso_y,3), 
                             ridge_y = 1e-7, inits = "rational", 
                             nrstart = 1, include_rational = TRUE, seed = 11, stop_value = 1e-3, MAXITER = 1000)
  
  # ber and classifiation rate #
  disc_fit <- scale_zero(X_train %*% disc1$W %*% t(disc1$Py),F,T)
  
  disc_fit <- 1 - (sum((Y_train - disc_fit)^2) / sum(Y_train^2))
  
  disc_pred_scale <- scale_zero(X_test %*% disc1$W %*% t(disc1$Py),F,T)
  
  disc_pred <- 1 - (sum((Y_test - disc_pred_scale)^2) / sum(Y_test^2))
  
  disc_pred0 <- 1 - (sum((Y_test[,pyzeros] - disc_pred_scale[,pyzeros])^2) / sum(Y_test[,pyzeros]^2))
  
  disc_pred_nonzero <- 1 - (sum((Y_test[,-pyzeros] - disc_pred_scale[,-pyzeros])^2) / sum(Y_test[,-pyzeros]^2))
  
  
  # number of outcome variables included by the final model
  disc_Y_included_index  <- rowSums(abs(disc1$Py)) > 1e-5
  
  disc_Y_included  <- sum(rowSums(abs(disc1$Py)) > 1e-5)
  
  true_Y_included_index <- rowSums(abs(Py)) > 1e-5
  
  disc_Y_correct <- (sum((abs(disc_Y_included_index) > 1e-7) + (abs(true_Y_included_index) > 1e-7) == 2) + 
                       sum((abs(disc_Y_included_index) < 1e-7) + (abs(true_Y_included_index) < 1e-7) == 2)) / (length(true_Y_included_index))
  
  
  # number of outcome variables included in each CV fold, and correct in each CV fold
  disc_CV_Y_included <- disc_chosen$Y_included_CV
  
  disc_CV_Y_correct <- disc_chosen$Y_correct_CV
  
  # tucker congruence #
  tucks <- RegularizedSCA::TuckerCoef(truePx, disc1$W)
  
  disc_tucker <- tucks$tucker_value
  
  # nonzero #
  disc_nonzero <- nonzero(estimate = disc1$W, defined = truePx)
  
  
  # correct classificaiton #
  disc_correct <- corrects(estimate = disc1$W, defined = truePx)
  
  # Py correct classification #
  disc1_py_zeros <- colSums(disc1$Py == 0)
  
  disc1_py_zero_index <- which(disc1_py_zeros == ncol(Y_train))
  
  # we need to watch out for situations where Py columns are estimated as zero
  
  if (length(disc1_py_zero_index) == 3){
    
    disc_py_correct <- sum(disc1$Py == Py) / (prod(dim(Py)))
    
  } else if (length(disc1_py_zero_index == 2)){
    
    Py_column <- disc1$Py[,-disc1_py_zero_index]
    
    disc_py1 <- corrects_vec(estimate = Py_column, Py[,1])
    disc_py2 <- corrects_vec(estimate = Py_column, Py[,2])
    disc_py3 <- corrects_vec(estimate = Py_column, Py[,3])
    
    match_index <- which.max(c(disc_py1, disc_py2, disc_py3))
    
    # the remaining columns are zero columns:
    remaining_corrects <- corrects(disc1$Py[,-match_index], Py[,-match_index])
    
    disc_py_correct <- max(c(disc_py1, disc_py2, disc_py3))/3 + remaining_corrects * 2/3
    
  } else if (length(disc1_py_zero_index == 1)){
    
    Py_column <- disc1$Py[,disc1_py_zero_index]
    
    disc_py1 <- corrects_vec(estimate = Py_column, Py[,1])
    disc_py2 <- corrects_vec(estimate = Py_column, Py[,2])
    disc_py3 <- corrects_vec(estimate = Py_column, Py[,3])
    
    match_index <- which.max(c(disc_py1, disc_py2, disc_py3))
    
    # the remaining columns are zero columns:
    remaining_corrects <- corrects(disc1$Py[,-match_index], Py[,-match_index])
    
    disc_py_correct <- max(c(disc_py1, disc_py2, disc_py3))/3 + remaining_corrects * 2/3
    
  } else {
    
    disc_py_correct <- corrects(estimate = disc1$Py, defined = Py)
    
  }
  
  rm(Py_column, disc_py1, disc_py2, disc_py3, match_index, remaining_corrects)
  
  disc_wzeros <- apply(disc1$W == 0,2,sum)
  
  
  # discovR - sequential CV ####
  
  # The sequential approach doesn't make sense with the new model selection at all,
  # because the model selection is all about fitting the model parameters on the entire dataset
  # to evaluate the sparsity level
  
  
  # # initial lasso glasso values per condition #
  # if (cond$dimension == "low"){
  #   lasso_w <-  rep(0.015,3)
  # } else {
  #   
  #   # for high dimesnionality
  #   lasso_w <- rep(0.003,3)
  #   
  # }
  # 
  # lasso_y_range <- c(exp(seq(log(0.00001), log(0.5), length.out = 7)))
  # 
  # alpha_range <- seq(0.1, 0.9, 0.1)
  # 
  # ranges <- expand.grid(alpha_range, lasso_y_range)
  # 
  # colnames(ranges) <- c("alpha", "lasso_y")
  # 
  # ranges <- ranges[order(ranges$lasso_y, decreasing = T),]
  # 
  # ranges <- ranges[order(ranges$alpha, decreasing = F),]
  # 
  # scd_cv_results <- c()
  # scd_cv_conv <- c()
  # 
  # disc_cv_results <- matrix(NA, ncol = 4+2, nrow = nrow(ranges))
  # 
  # colnames(disc_cv_results) <- c("cve", "se", "conv", "Y_zero", "alpha", "lasso_y")
  # 
  # disc_cv_results <- as.data.frame(disc_cv_results)
  # 
  # colnames(ranges) <- c("alpha", "lasso_y")
  # 
  # for (i in 1:nrow(ranges)){
  #   
  #   disc_cv1 <- discovR_cv_rsq(X = X_train, Y = Y_train, 
  #                              R = R, 
  #                              alpha = ranges[i,]$alpha, 
  #                              lasso_w = lasso_w, 
  #                              ridge_w = 1e-7, 
  #                              lasso_y = rep(ranges[i,]$lasso_y,3),
  #                              ridge_y = 1e-7, 
  #                              seed_cv = i+1, seed_method = 111,
  #                              inits = "rational", include_rational = TRUE, 
  #                              MAXITER = 500, nrFolds = 5, nonzero_only = TRUE, stop_value = 1e-3)
  #   
  #   disc_cv_results[i,] <- c(disc_cv1$cve, disc_cv1$se, sum(disc_cv1$convergence), disc_cv1$zero_Y,
  #                            ranges[i,]$alpha, ranges[i,]$lasso_y)
  #   
  #   print(i)
  # }
  # 
  # time2 <- Sys.time()
  # 
  # disc_cv_results <- disc_cv_results[!is.na(disc_cv_results$cve),]
  # 
  # cv_max_index <- which.max(disc_cv_results[,1]) # largest R2
  # 
  # serule <- disc_cv_results[cv_max_index,1] - disc_cv_results[cv_max_index,2] # 1SE rule
  # 
  # disc_possible <- disc_cv_results[disc_cv_results[,1] >= serule,] # larger than the 1SE rule
  # 
  # if (is.vector(disc_possible)){
  #   disc_chosen <- disc_possible
  # } else {
  #   disc_chosen <- disc_possible[1,]
  # }
  # 
  # 
  # disc_seq_alpha <- disc_chosen$alpha
  # disc_seq_lasso_y <- disc_chosen$lasso_y
  # 
  # # now cross-validation for lasso-w #
  # 
  # lasso_w_range <- c(exp(seq(log(0.00001), log(0.5), length.out = 30)))
  # 
  # ranges <- as.data.frame(matrix(sort(lasso_w_range, decreasing = T), ncol = 1))
  # 
  # colnames(ranges) <- c("lasso_w")
  # 
  # scd_cv_results <- c()
  # scd_cv_conv <- c()
  # 
  # disc_cv_results <- matrix(NA, ncol = 4+1, nrow = nrow(ranges))
  # 
  # colnames(disc_cv_results) <- c("cve", "se", "conv", "Y_zero", "lasso_w")
  # 
  # disc_cv_results <- as.data.frame(disc_cv_results)
  # 
  # for (i in 1:nrow(ranges)){
  #   
  #   disc_cv1 <- discovR_cv_rsq(X = X_train, Y = Y_train, 
  #                              R = R, 
  #                              alpha = disc_seq_alpha, 
  #                              lasso_w = rep(ranges[i,], 3), 
  #                              ridge_w = 1e-7, 
  #                              lasso_y = rep(disc_seq_lasso_y,3),
  #                              ridge_y = 1e-7, 
  #                              seed_cv = i+1, seed_method = 111,
  #                              inits = "rational", include_rational = TRUE, 
  #                              MAXITER = 500, nrFolds = 5, nonzero_only = TRUE, stop_value = 1e-3)
  #   
  #   disc_cv_results[i,] <- c(disc_cv1$cve, disc_cv1$se, sum(disc_cv1$convergence), disc_cv1$zero_Y,
  #                            ranges[i,])
  #   
  #   print(i)
  # }
  # 
  # time2 <- Sys.time()
  # 
  # disc_cv_results <- disc_cv_results[!is.na(disc_cv_results$cve),]
  # 
  # cv_max_index <- which.max(disc_cv_results[,1]) # largest R2
  # 
  # serule <- disc_cv_results[cv_max_index,1] - disc_cv_results[cv_max_index,2] # 1SE rule
  # 
  # disc_possible <- disc_cv_results[disc_cv_results[,1] >= serule,] # larger than the 1SE rule
  # 
  # if (is.vector(disc_possible)){
  #   disc_chosen <- disc_possible
  # } else {
  #   disc_chosen <- disc_possible[1,]
  # }
  # 
  # 
  # disc_seq_lasso_w <- disc_chosen$lasso_w
  # 
  # 
  # # final discovR model fitting #
  # disc_seq <- discovR_cpp_ridge(X = X_train, Y = Y_train, R = R, 
  #                               alpha = disc_seq_alpha, 
  #                               lasso_w = rep(disc_seq_lasso_w,3), 
  #                               ridge_w = 1e-7,
  #                               lasso_y = rep(disc_seq_lasso_y,3), 
  #                               ridge_y = 1e-7, inits = "rational", 
  #                               nrstart = 1, include_rational = TRUE, seed = 11, stop_value = 1e-5)
  # 
  # 
  # # ber and classifiation rate #
  # disc_seq_fit <- scale_zero(X_train %*% disc_seq$W %*% t(disc_seq$Py),F,T)
  # 
  # disc_seq_fit <- 1 - (sum((Y_train - disc_seq_fit)^2) / sum(Y_train^2))
  # 
  # disc_seq_pred <- scale_zero(X_test %*% disc_seq$W %*% t(disc_seq$Py),F,T)
  # 
  # disc_seq_pred <- 1 - (sum((Y_test - disc_seq_pred)^2) / sum(Y_test^2))
  # 
  # 
  # # tucker congruence #
  # tucks <- RegularizedSCA::TuckerCoef(truePx, disc_seq$W)
  # 
  # disc_seq_tucker <- tucks$tucker_value
  # 
  # # nonzero #
  # disc_seq_nonzero <- nonzero(estimate = disc_seq$W, defined = truePx)
  # 
  # 
  # # correct classificaiton #
  # disc_seq_correct <- corrects(estimate = disc_seq$W, defined = truePx)
  # 
  # # Py correct classification #
  # 
  # disc_seq_py_zeros <- colSums(disc_seq$Py == 0)
  # 
  # disc_seq_py_zero_index <- which(disc_seq_py_zeros == ncol(Y_train))
  # 
  # # we need to watch out for situations where Py columns are estimated as zero
  # 
  # if (length(disc_seq_py_zero_index) == 3){
  #   
  #   disc_seq_py_correct <- sum(disc_seq$Py == Py) / (prod(dim(Py)))
  #   
  # } else if (length(disc_seq_py_zero_index == 2)){
  #   
  #   Py_column <- disc_seq$Py[,-disc_seq_py_zero_index]
  #   
  #   disc_py1 <- corrects_vec(estimate = Py_column, Py[,1])
  #   disc_py2 <- corrects_vec(estimate = Py_column, Py[,2])
  #   disc_py3 <- corrects_vec(estimate = Py_column, Py[,3])
  #   
  #   match_index <- which.max(c(disc_py1, disc_py2, disc_py3))
  #   
  #   # the remaining columns are zero columns:
  #   remaining_corrects <- corrects(disc_seq$Py[,-match_index], Py[,-match_index])
  #   
  #   disc_seq_py_correct <- max(c(disc_py1, disc_py2, disc_py3))/3 + remaining_corrects * 2/3
  #   
  # } else if (length(disc_seq_py_zero_index == 1)){
  #   
  #   Py_column <- disc_seq$Py[,disc_seq_py_zero_index]
  #   
  #   disc_py1 <- corrects_vec(estimate = Py_column, Py[,1])
  #   disc_py2 <- corrects_vec(estimate = Py_column, Py[,2])
  #   disc_py3 <- corrects_vec(estimate = Py_column, Py[,3])
  #   
  #   match_index <- which.max(c(disc_py1, disc_py2, disc_py3))
  #   
  #   # the remaining columns are zero columns:
  #   remaining_corrects <- corrects(disc_seq$Py[,-match_index], Py[,-match_index])
  #   
  #   disc_seq_py_correct <- max(c(disc_py1, disc_py2, disc_py3))/3 + remaining_corrects * 2/3
  #   
  # } else {
  #   
  #   disc_seq_py_correct <- corrects(estimate = disc_seq$Py, defined = Py)
  #   
  # }
  # 
  # rm(Py_column, disc_py1, disc_py2, disc_py3, match_index, remaining_corrects)
  # 
  # 
  # # common and distinctive component identification #
  # disc_seq_comdis <- comdis(estimate = disc_seq$W)
  # 
  # disc_seq_d1 <- disc_seq_comdis$d1
  # disc_seq_d2 <- disc_seq_comdis$d2
  # disc_seq_common <- disc_seq_comdis$common
  # 
  # disc_seq_wzeros <- apply(disc_seq$W == 0,2,sum)
  # 
  # 
  # # final number of outcome variables included in the model
  # disc_seq_Y_included <- sum(rowSums(abs(disc_seq$Py)) > 1e-7) # this number should either be 20 or 12, ideally. 
  
  
  
  # previous method SSCOVR (lasso for Py = 0) - exhaustive ####
  
  # (for this method, large lasso and small alpha is selected)
  
  # cross validation for alpha and ridge #
  time1 <- Sys.time()
  
  # ranges of alpha and ridge values considered
  alpha_range <- seq(0.1, 0.9, 0.1)
  
  lasso_w_range <- c(exp(seq(log(0.00001), log(0.5), length.out = 15)))
  
  ranges <- expand.grid(alpha_range, lasso_w_range)
  
  colnames(ranges) <- c("alpha", "lasso_w")
  
  ranges <- ranges[order(ranges$lasso_w, decreasing = T),]
  
  ranges <- ranges[order(ranges$alpha, decreasing = F),]
  
  lasso0_cv_results <- c()
  lasso0_cv_conv <- c()
  
  lasso0_cv_results <- matrix(NA, ncol = 3+2, nrow = nrow(ranges))
  
  colnames(lasso0_cv_results) <- c("R2cv", "SE", "conv", "alpha", "lasso")
  
  for (i in 1:nrow(ranges)){
    
    lasso0_cv1 <- discovR_cv_rsq(X = X_train, Y = Y_train, R = R, 
                                 alpha = ranges[i,]$alpha, 
                                 lasso_w = rep(ranges[i,]$lasso_w, R), 
                                 ridge_w = 1e-7,
                                 lasso_y = rep(0,3),
                                 ridge_y = 1e-7, 
                                 seed_cv = i+1, seed_method = 111,
                                 inits = "rational", include_rational = TRUE, 
                                 MAXITER = 100, nrFolds = 5, stop_value = 1e-3, nonzero_only = FALSE,
                                 true_active_Y = rowSums(abs(Py)) > 1e-5)
    
    lasso0_cv_results[i,] <- c(lasso0_cv1$cve, lasso0_cv1$se, sum(lasso0_cv1$convergence),
                               ranges[i,]$alpha, ranges[i,]$lasso_w)
    
    print(i)
  }
  
  time2 <- Sys.time()
  
  lasso0_cv_results <- lasso0_cv_results[complete.cases(lasso0_cv_results),]
  
  cv_max_index <- which.max(lasso0_cv_results[,1]) # largest R2_cv
  
  serule <- lasso0_cv_results[cv_max_index,1] - lasso0_cv_results[cv_max_index,2] # 1SE rule
  
  lasso0_possible <- lasso0_cv_results[lasso0_cv_results[,1] >= serule,]
  
  if (is.vector(lasso0_possible)){
    lasso0_chosen <- lasso0_possible
  } else {
    lasso0_chosen <- lasso0_possible[1,]
  }
  
  lasso0_alpha <- lasso0_chosen[4]
  lasso0_lasso_w <- lasso0_chosen[5]
  
  # final discovR model fitting #
  lasso0_exhaustive <- discovR_cpp_ridge(X = X_train, Y = Y_train, R = R, 
                                         alpha = lasso0_alpha, 
                                         lasso_w = rep(lasso0_lasso_w, 3), 
                                         ridge_w = 1e-7,
                                         lasso_y = rep(0,3), 
                                         ridge_y = 1e-7, inits = "rational", 
                                         nrstart = 1, include_rational = TRUE, seed = 11, stop_value = 1e-3, MAXITER = 1000)
  
  
  
  # ber and classifiation rate #
  lasso0_fit <- scale_zero(X_train %*% lasso0_exhaustive$W %*% t(lasso0_exhaustive$Py),F,T)
  
  lasso0_fit <- 1 - (sum((Y_train - lasso0_fit)^2) / sum(Y_train^2))
  
  lasso0_pred_scale <- scale_zero(X_test %*% lasso0_exhaustive$W %*% t(lasso0_exhaustive$Py),F,T)
  
  lasso0_pred <- 1 - (sum((Y_test - lasso0_pred_scale)^2) / sum(Y_test^2))
  
  lasso0_pred0 <- 1 - (sum((Y_test[,pyzeros] - lasso0_pred_scale[,pyzeros])^2) / sum(Y_test[,pyzeros]^2))
  
  lasso0_pred_nonzero <- 1 - (sum((Y_test[,-pyzeros] - lasso0_pred_scale[,-pyzeros])^2) / sum(Y_test[,-pyzeros]^2))
  
  
  # tucker congruence #
  tucks <- RegularizedSCA::TuckerCoef(truePx, lasso0_exhaustive$W)
  
  lasso0_tucker <- tucks$tucker_value
  
  # nonzero #
  lasso0_nonzero <- nonzero(estimate = lasso0_exhaustive$W, defined = truePx)
  
  
  # correct classificaiton #
  lasso0_correct <- corrects(estimate = lasso0_exhaustive$W, defined = truePx)
  
  # Py correct classification #
  lasso0_py_correct <- corrects(estimate = lasso0_exhaustive$Py, defined = Py)
  
  # I want to know how the regression coefficients that were meant to be estimated by zero 
  # behaves like
  
  # I sort the regression coefficients in a descending order with respect to magnitude
  lasso0_abs <- sort(abs(lasso0_exhaustive$Py), decreasing = T)
  
  # taking the average of the coefficients that were meant to be estimated as zero:
  lasso0_zerocoefs <- mean(lasso0_abs[-c(1:sum(Py))])
  
  lasso0_zerocoefs_firstfive <- mean((lasso0_abs[-c(1:sum(Py))])[1:5])
  
  lasso0_wzeros <- apply(lasso0_exhaustive$W == 0,2,sum)
  
  
  # previous method SSCOVR (lasso for Py = 0) - sequential ####
  
  # The sequential approach does not make sense in this batch of simulation study,
  # where the model selection was done by studying the sparsity level of applying the parameters to 
  # the dataset entirely
  
  # # initial lasso glasso values per condition #
  # if (cond$dimension == "low"){
  #   lasso_w <-  rep(0.015,3)
  # } else {
  #   # for high dimesnionality
  #   lasso_w <- rep(0.003,3)
  # }
  # 
  # # cross validation for alpha and ridge #
  # time1 <- Sys.time()
  # 
  # # ranges of alpha and ridge values considered
  # alpha_range <- seq(0.1, 0.9, 0.1)
  # 
  # ridge_range <-  c(1e-7)
  # 
  # lasso_y_range <- c(0)
  # 
  # ranges <- expand.grid(alpha_range, ridge_range, lasso_y_range)
  # 
  # colnames(ranges) <- c("alpha", "ridge", "lasso_y")
  # 
  # ranges <- ranges[order(ranges$lasso_y, decreasing = T),]
  # 
  # ranges <- ranges[order(ranges$alpha, decreasing = F),]
  # 
  # lasso0_cv_results <- c()
  # lasso0_cv_conv <- c()
  # 
  # lasso0_cv_results <- matrix(NA, ncol = 3+3, nrow = nrow(ranges))
  # 
  # colnames(ranges) <- c("alpha", "ridge", "lasso_y")
  # 
  # for (i in 1:nrow(ranges)){
  #   
  #   lasso0_cv1 <- discovR_cv_rsq(X = X_train, Y = Y_train, R = R, 
  #                                alpha = ranges[i,]$alpha, 
  #                                lasso_w = lasso_w, 
  #                                ridge_w = 1e-7,
  #                                lasso_y = rep(0,3),
  #                                ridge_y = ranges[i,]$ridge, 
  #                                seed_cv = i+1, seed_method = 111,
  #                                inits = "rational", include_rational = TRUE, 
  #                                MAXITER = 100, nrFolds = 5, stop_value = 1e-3, nonzero_only = TRUE)
  #   
  #   lasso0_cv_results[i,] <- c(lasso0_cv1$cve, lasso0_cv1$se, sum(lasso0_cv1$convergence),
  #                              ranges[i,]$alpha, ranges[i,]$ridge, ranges[i,]$lasso_y)
  #   
  #   print(i)
  # }
  # 
  # time2 <- Sys.time()
  # 
  # cv_max_index <- which.max(lasso0_cv_results[,1]) # smallest CVE
  # 
  # serule <- lasso0_cv_results[cv_max_index,1] - lasso0_cv_results[cv_max_index,2] # 1SE rule
  # 
  # lasso0_possible <- lasso0_cv_results[lasso0_cv_results[,1] >= serule,]
  # 
  # if (is.vector(lasso0_possible)){
  #   lasso0_chosen <- lasso0_possible
  # } else {
  #   lasso0_chosen <- lasso0_possible[1,]
  # }
  # 
  # lasso0_seq_alpha <- lasso0_chosen[4]
  # lasso0_seq_ridge <- lasso0_chosen[5]
  # lasso0_seq_lasso_y <- 0
  # 
  # # now cross-validation for lasso-w #
  # lasso_w_range <- c(exp(seq(log(0.00001), log(0.5), length.out = 30)))
  # 
  # ranges <- as.data.frame(matrix(sort(lasso_w_range, decreasing = T), ncol = 1))
  # 
  # colnames(ranges) <- c("lasso_w")
  # 
  # scd_cv_results <- c()
  # scd_cv_conv <- c()
  # 
  # lasso0_cv_results <- matrix(NA, ncol = 4+1, nrow = nrow(ranges))
  # 
  # colnames(lasso0_cv_results) <- c("cve", "se", "conv", "Y_zero", "lasso_w")
  # 
  # lasso0_cv_results <- as.data.frame(lasso0_cv_results)
  # 
  # for (i in 1:nrow(ranges)){
  #   
  #   lasso0_cv1 <- discovR_cv_rsq(X = X_train, Y = Y_train, 
  #                                R = R, 
  #                                alpha = lasso0_seq_alpha, 
  #                                lasso_w = rep(ranges[i,], 3), 
  #                                ridge_w = 1e-7, 
  #                                lasso_y = rep(lasso0_seq_lasso_y,3),
  #                                ridge_y = 1e-7, 
  #                                seed_cv = i+1, seed_method = 111,
  #                                inits = "rational", include_rational = TRUE, 
  #                                MAXITER = 500, nrFolds = 5, nonzero_only = TRUE, stop_value = 1e-3)
  #   
  #   lasso0_cv_results[i,] <- c(lasso0_cv1$cve, lasso0_cv1$se, sum(lasso0_cv1$convergence), lasso0_cv1$zero_Y,
  #                              ranges[i,])
  #   
  #   print(i)
  # }
  # 
  # time2 <- Sys.time()
  # 
  # lasso0_cv_results <- lasso0_cv_results[!is.na(lasso0_cv_results$cve),]
  # 
  # cv_max_index <- which.max(lasso0_cv_results[,1]) # largest R2
  # 
  # serule <- lasso0_cv_results[cv_max_index,1] - lasso0_cv_results[cv_max_index,2] # 1SE rule
  # 
  # lasso0_possible <- lasso0_cv_results[lasso0_cv_results[,1] >= serule,] # larger than the 1SE rule
  # 
  # if (is.vector(lasso0_possible)){
  #   lasso0_chosen <- lasso0_possible
  # } else {
  #   lasso0_chosen <- lasso0_possible[1,]
  # }
  # 
  # lasso0_seq_lasso_w <- lasso0_chosen$lasso_w
  # 
  # 
  # # final discovR model fitting #
  # lasso0_seq <- discovR_cpp_ridge(X = X_train, Y = Y_train, R = R, 
  #                                 alpha = lasso0_seq_alpha, 
  #                                 lasso_w = rep(lasso0_seq_lasso_w,3), 
  #                                 ridge_w = 1e-7,
  #                                 lasso_y = rep(0,3), 
  #                                 ridge_y = 1e-7, inits = "rational", 
  #                                 nrstart = 1, include_rational = TRUE, seed = 11, stop_value = 1e-5)
  # 
  # 
  # 
  # # ber and classifiation rate #
  # lasso0_seq_fit <- scale_zero(X_train %*% lasso0_seq$W %*% t(lasso0_seq$Py),F,T)
  # 
  # lasso0_seq_fit <- 1 - (sum((Y_train - lasso0_seq_fit)^2) / sum(Y_train^2))
  # 
  # lasso0_seq_pred <- scale_zero(X_test %*% lasso0_seq$W %*% t(lasso0_seq$Py),F,T)
  # 
  # lasso0_seq_pred <- 1 - (sum((Y_test - lasso0_seq_pred)^2) / sum(Y_test^2))
  # 
  # 
  # # tucker congruence #
  # tucks <- RegularizedSCA::TuckerCoef(truePx, lasso0_seq$W)
  # 
  # lasso0_seq_tucker <- tucks$tucker_value
  # 
  # # nonzero #
  # lasso0_seq_nonzero <- nonzero(estimate = lasso0_seq$W, defined = truePx)
  # 
  # 
  # # correct classificaiton #
  # lasso0_seq_correct <- corrects(estimate = lasso0_seq$W, defined = truePx)
  # 
  # # Py correct classification #
  # lasso0_seq_py_correct <- corrects(estimate = lasso0_seq$Py, defined = Py)
  # 
  # # common and distinctive component identification #
  # lasso0_seq_comdis <- comdis(estimate = lasso0_seq$W)
  # 
  # lasso0_seq_d1 <- lasso0_seq_comdis$d1
  # lasso0_seq_d2 <- lasso0_seq_comdis$d2
  # lasso0_seq_common <- lasso0_seq_comdis$common
  # 
  # lasso0_seq_wzeros <- apply(lasso0_seq$W == 0,2,sum)
  
  
  # DIABLO-concatenated block.spls ####
  
  # finding the level of sparsity by CV
  
  # first arranging the objects for the diablo function #
  design <- matrix(1, ncol = 1, nrow = 1)
  
  block1 <- data.frame(X_train)
  
  plsX <- list(block1 = block1)
  
  
  # if (cond$dimension == "low"){
  #   list.keepX = list(block1 = c(4,8,4))
  #   
  # } else {
  #   list.keepX = list(block1 = c(25, 40, 25)) 
  # }
  
  if (cond$dimension == "low"){
    
    comp1_range <- c(4,8, 12, 16, seq(20, ncol(X_train), 8))
    comp2_range <- c(4,8, 12, 16, seq(20, ncol(X_train), 8))
    comp3_range <- c(4,8, 12, 16, seq(20, ncol(X_train), 8))
    
  } else {
    
    comp1_range <- unique(sort(c(seq(25, ncol(X_train), 25), seq(40, ncol(X_train), 40)), decreasing = F))
    comp2_range <- unique(sort(c(seq(25, ncol(X_train), 25), seq(40, ncol(X_train), 40)), decreasing = F))
    comp3_range <- unique(sort(c(seq(25, ncol(X_train), 25), seq(40, ncol(X_train), 40)), decreasing = F))
    
  }
  
  ranges <- expand.grid(comp1_range, comp2_range, comp3_range)
  
  colnames(ranges) <- c("comp1", "comp2", "comp3")
  
  pls_cv_results <- matrix(NA, ncol = 2+3, nrow = nrow(ranges))
  
  colnames(pls_cv_results) <- c("cve", "se",  "comp1", "comp2", "comp3")
  
  pls_cv_results <- as.data.frame(pls_cv_results)
  
  for (i in 1:nrow(ranges)){
    
    list.keepX <- list(block1 = c(ranges[i,]$comp1, ranges[i,]$comp2, ranges[i,]$comp3))
    
    pls_cv1 <- sPLS_cv(plsX = plsX, Y = Y_train, ncomp = R, nrFolds = 5, keepX = list.keepX, design = design, seed_cv = i)
    
    pls_cv_results[i,] <- c(pls_cv1$cve, pls_cv1$se, 
                            ranges[i,])
    
    print(i)
  }
  
  time2 <- Sys.time()
  
  cv_max_index <- which.max(pls_cv_results[,1]) # largest R2
  
  serule <- pls_cv_results[cv_max_index,1] - pls_cv_results[cv_max_index,2] # 1SE rule
  
  pls_possible <- pls_cv_results[pls_cv_results[,1] >= serule,] # larger than the 1SE rule
  
  if (is.vector(pls_possible)){
    pls_chosen <- pls_possible
  } else {
    pls_chosen <- pls_possible[1,]
  }
  
  pls_keep1 <- pls_chosen$comp1
  pls_keep2 <- pls_chosen$comp2
  pls_keep3 <- pls_chosen$comp3
  
  
  list.keepX <- list(block1 = c(pls_keep1, pls_keep2, pls_keep3))
  
  # actual model fitting
  diacon1 <- block.spls(X = plsX, Y = Y_train,
                        ncomp = R, keepX = list.keepX, design = design)
  
  # fit #
  diacon_fit <- predict(object = diacon1, newdata = plsX)
  
  diacon_fit_scale <- scale_zero(diacon_fit$predict$block1[,,3],F,T)
  
  diacon_fit <- 1 - (sum((diacon_fit_scale - Y_train)^2) / sum(Y_train^2))
  
  
  block1test <- data.frame(X_test)
  
  plsXtest <- list(block1 = block1test)
  
  diacon_pred <- predict(object = diacon1,
                         newdata = plsXtest)
  
  diacon_pred_scale <- scale_zero(diacon_pred$predict$block1[,,3], F, T)
  
  diacon_pred <- 1 - (sum((diacon_pred_scale - Y_test)^2) / sum(Y_test^2))
  
  diacon_pred0 <- 1 - (sum((Y_test[,pyzeros] - diacon_pred_scale[,pyzeros])^2) / sum(Y_test[,pyzeros]^2))
  
  diacon_pred_nonzero <- 1 - (sum((Y_test[,-pyzeros] - diacon_pred_scale[,-pyzeros])^2) / sum(Y_test[,-pyzeros]^2))
  
  
  
  # tucker congruence #
  diacon_tucks <- RegularizedSCA::TuckerCoef(truePx,
                                             diacon1$loadings$block1)
  
  diacon_tucker <- diacon_tucks$tucker_value
  
  # correct classification #
  diacon_correct <- corrects(estimate = diacon1$loadings$block1,
                             defined = truePx)
  
  # nonzero classification #
  diacon_nonzero <- nonzero(estimate = diacon1$loadings$block1,
                            defined = truePx)
  
  
  # Py correct classification #
  diacon_py_correct <- corrects(estimate = diacon1$loadings$Y, defined = Py)
  
  # I sort the regression coefficients in a descending order with respect to magnitude
  diacon_abs <- sort(abs(diacon1$loadings$Y), decreasing = T)
  
  # taking the average of the coefficients that were meant to be estimated as zero:
  diacon_zerocoefs <- mean(diacon_abs[-c(1:sum(Py))])
  
  diacon_zerocoefs_firstfive <- mean((diacon_abs[-c(1:sum(Py))])[1:5])
  
  
  
  
  
  
  # recording the results #
  
  performances <- t(c(disc_alpha, 1e-7, disc_lasso_y,                      
                      disc_lasso[1], 1e-7,
                      disc_fit, disc_pred,
                      disc_pred0, disc_pred_nonzero,
                      disc_tucker,
                      disc_correct, 
                      disc_nonzero,
                      disc_py_correct,
                      disc_wzeros[1], disc_wzeros[2], disc_wzeros[3],
                      disc_Y_included, disc_Y_correct,
                      disc_CV_Y_included, disc_CV_Y_correct,
                      
                      lasso0_alpha, 1e-7, 0,                      
                      lasso0_lasso_w[1], 1e-7,
                      lasso0_fit, lasso0_pred,
                      lasso0_pred0, lasso0_pred_nonzero,
                      lasso0_tucker,
                      lasso0_correct, 
                      lasso0_nonzero,
                      lasso0_py_correct,
                      lasso0_wzeros[1], lasso0_wzeros[2], lasso0_wzeros[3],
                      lasso0_zerocoefs, lasso0_zerocoefs_firstfive,
                      
                      diacon_fit, diacon_pred, 
                      diacon_pred0, diacon_pred_nonzero,
                      diacon_tucker,
                      diacon_correct, 
                      diacon_nonzero, 
                      diacon_py_correct,
                      diacon_zerocoefs, diacon_zerocoefs_firstfive
  ))
  
  
  resulting <- data.frame(cond$dimension, cond$signal_level, cond$Jy, cond$py_pattern,
                          cond$reps, modelseeding, noiseseeding)
  
  resulting <- cbind(resulting, performances)
  
  # recording the current result in the data frame that I've defined above
  results[rrr,] <- resulting
  
  # if there was a problem, stop replication
  if(anyNA(results[rrr,])){
    stop("NA resulted")
  }
  
  print(rep(rrr, 10))
  
  
  # saving the current "results" object, in case it crashes in the middle 
  save("results", file = "./run_and_results/results1_24_nov_2022.Rdata")
  write.table(x = results, file = "./run_and_results/results1_24_nov_2022.txt", sep=",")
  
  flush.console()
  
}

