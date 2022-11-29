# discovR toy CV Rsq#
# 24th May 2022 #
# Soogeun Park #

# I try to do toy example where I vary the alpha value to see what happens
# I generate the data based on the simulation study 18-oct-2021 

# rep 2 #

# 1. functions and libraries ####
library(foreign)
library(psych)
library(RegularizedSCA)
library(nFactors)
library(mixOmics)
library(Rcpp)

source( "E:\\Users\\park\\Desktop\\Functions\\discovR_cpp_ridge.R")
source("E:\\Users\\park\\Desktop\\Functions\\dategen_3cov_discovR.R")
source("E:\\Users\\park\\Desktop\\Functions\\discovR_cv_rsq.R")

sourceCpp( "E:\\Users\\park\\Desktop\\Functions\\updateW_discovR_ridge.cpp")

setwd("E:\\Users\\park\\Desktop\\CV Rsq toy experiment ridge/")

# just loading the condition_df script
source("E:\\Users\\park\\Desktop\\simulation 17-may-2022/conditions_making_17_may_2022.R")

# 1. data preparation #####
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


# 2. replication starts ####

# I'm not doing replication #

# set.seed(1111) 
# model_seed <- sample(x = 1:100000, size = nrow(condition_df))
# noise_seed <- sample(x = 1:100000, size = nrow(condition_df))
# 
# dim(condition_df)
# 
# unique(condition_df$signal_level)
# 
# # dividing the task such that each run takes equal number of high dimensionality problem #
# repshigh <- rownames(subset(condition_df, dimension == "high"))[1:80]
# repslow <- rownames(subset(condition_df, dimension == "low"))[1:80]
# 
# reps_to_do <- c(as.numeric(repslow), as.numeric(repshigh))
# 
# length(reps_to_do)



# I think it may be nice to compare between:
# high, vaf = 0.9, Jy = 20, all3
# high, vaf = 0.5, Jy = 20, all3
# high, vaf = 0.9, Jy = 20, notall2
# high, vaf = 0.5, Jy = 20, notall2


rrr <-  32

cond <- condition_df[rrr,]

cond$py_pattern <- "notall2"

cond

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


cond$signal_level <- as.numeric(cond$signal_level)
cond$reps <- as.numeric(cond$reps)
cond$Jy <- as.numeric(cond$Jy)

# data generating #
# VAF050 Jy20 notall2 seed = 333 #
dat <- data_3cov(I = I, J = J, Jy = cond$Jy, d1_active = d1_active, 
                 c_active_block1 = c_active_block1, c_active_block2 = c_active_block2,
                 d2_active = d2_active, comp_sd = c(50,50,50), 
                 signal_level = cond$signal_level,
                 strict_uncorrelated = F,
                 Py = Py, seed = 333)



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

# selecting the number of covariates
svdd <- svd(cbind(Y_train, X_train))

nScree(svdd$d)

nScree(svdd$d^2)

R <- 3


lasso_range <- c(0,exp(seq(log(0.00001), log(0.5), length.out = 14)))

lasso_y_range <- c(0,exp(seq(log(0.00001), log(0.5), length.out = 19)))

alpha_range <- seq(0.1, 0.9, 0.1)

ranges <- expand.grid(alpha_range, lasso_y_range, lasso_range)

colnames(ranges) <- c("alpha", "lasso_y", "lasso")


ranges <- ranges[order(ranges$lasso_y, decreasing = T),]

ranges <- ranges[order(ranges$alpha, decreasing = F),]

ranges <- ranges[order(ranges$lasso, decreasing = T),]



results <- matrix(NA, ncol = 19, nrow = nrow(ranges))

colnames(results) <- c("x_insample", "y_insample", "x_outsample", "y_outsample",
                       "CV_Rsq_y", "SE_y", "CV_Rsq_x", "SE_x", "conv",
                       "w_correct", "w_nonzero", "w_nonzero_level", 
                       "py_correct", "py_nonzero",
                       "py_nonzero_level", "Y_nonzero", "alpha", "lasso_y", "lasso")

reps_total <- 1:nrow(ranges)

reps_to_do <- 1401:2700

for (i in reps_to_do){
  
  disc_cv1 <- discovR_cv_rsq(X = X_train, Y = Y_train, 
                             R = R, 
                             alpha = ranges[i,]$alpha, 
                             lasso_w = rep(ranges[i,]$lasso,R), 
                             ridge_w = 1e-10,
                             lasso_y = rep(ranges[i,]$lasso_y,3),
                             ridge_y = 1e-7, 
                             seed_cv = i+1, seed_method = 111,
                             inits = "rational", include_rational = TRUE, 
                             MAXITER = 500, nrFolds = 5, nonzero_only = TRUE, stop_value = 1e-3)
  
  disc1 <- discovR_cpp_ridge(X = X_train, Y = Y_train, R = R, 
                       alpha = ranges[i,]$alpha, 
                       lasso_w = rep(ranges[i,]$lasso,3), 
                       ridge_w = 1e-10,
                       lasso_y = rep(ranges[i,]$lasso_y,3), 
                       ridge_y = 1e-7, inits = "rational", 
                       nrstart = 1, include_rational = TRUE, seed = 22, stop_value = 1e-5)
  
  
  # R squared: Y #
  fitted_values <- X_train %*% disc1$W %*% t(disc1$Py)
  
  final_zeros <- colSums(abs(fitted_values)) == 0
  
  final_nonzeros <- as.numeric(which(final_zeros == F))
  
  y_insample <- 1 - sum((fitted_values[,final_nonzeros] - Y_train[,final_nonzeros])^2) / sum(Y_train[,final_nonzeros]^2)
  
  Y_nonzeros <- length(final_nonzeros)
  
  # R squared: X #
  fitted_values <- X_train %*% disc1$W %*% t(disc1$Px)
  
  final_zeros <- colSums(abs(fitted_values)) == 0
  
  final_nonzeros <- as.numeric(which(final_zeros == F))
  
  x_insample <- 1 - sum((fitted_values[,final_nonzeros] - X_train[,final_nonzeros])^2) / sum(X_train[,final_nonzeros]^2)
  
  # out-sample R squared: Y #
  fitted_values <- X_test %*% disc1$W %*% t(disc1$Py)
  
  final_zeros <- colSums(abs(fitted_values)) == 0
  
  final_nonzeros <- as.numeric(which(final_zeros == F))
  
  y_outsample <- 1 - sum((fitted_values[,final_nonzeros] - Y_test[,final_nonzeros])^2) / sum(Y_test[,final_nonzeros]^2)
  
  
  # R squared: X #
  fitted_values <- X_test %*% disc1$W %*% t(disc1$Px)
  
  final_zeros <- colSums(abs(fitted_values)) == 0
  
  final_nonzeros <- as.numeric(which(final_zeros == F))
  
  x_outsample <- 1 - sum((fitted_values[,final_nonzeros] - X_test[,final_nonzeros])^2) / sum(X_test[,final_nonzeros]^2)
  
  
  
  # correct classificaiton #
  disc_correct <- tryCatch(corrects(estimate = disc1$W, defined = truePx), error = function(e) NA)
  
  disc_nonzero_hit <- tryCatch(nonzero(estimate = disc1$W, defined = truePx), error = function(e) NA)
  
  # Py correct classification #
  disc_py_correct <- tryCatch(corrects(estimate = disc1$Py, defined = Py), error = function(e) NA)
  
  
  
  disc_py_nonzero_hit <- tryCatch(nonzero(estimate = disc1$Py, defined = Py), error = function(e) NA)
  
  # nonzero coefficients
  Py_nonzero <- sum(disc1$Py != 0)
  
  W_nonzero <- sum(disc1$W != 0)
  
  resulting <- c(x_insample, y_insample, x_outsample, y_outsample,
                 disc_cv1$cve, disc_cv1$se, 
                 disc_cv1$cve_x, disc_cv1$se_x, sum(disc_cv1$convergence),
                 disc_correct, disc_nonzero_hit, W_nonzero,
                 disc_py_correct, disc_py_nonzero_hit, 
                 Py_nonzero, Y_nonzeros,
                 ranges[i,]$alpha, ranges[i,]$lasso_y, ranges[i,]$lasso)
  
  results[i,] <- resulting
  save(results, file = "./run2_ridge_toy_high_vaf050_Jy20_notall2_CV_rsq.Rdata")
  
  print(i)
}



















head(results)

results <- as.data.frame(results)

results$y_insample



# results_lasso0005 <- results[results$alpha > 0.69 & results$alpha < 0.71,]
# 
# results_lasso0005 <- results[results$lasso_y == 0.001,]
# 
# results_lasso0005 <- results[results$lasso == 0.005,]



results$lasso_y

results$Y_nonzero


log_lasso <- log(results$lasso)

log_lasso[log_lasso == -Inf] <- -15

results$log_lasso <- log_lasso


log_lasso_y <- log(results$lasso_y)

log_lasso_y[log_lasso_y == -Inf] <- -15

results$log_lasso_y <- log_lasso_y



col_rsq <- c("y_insample", "y_outsample", "x_insample", "x_outsample")

results_long <- tidyr::gather(results, whichrsq, rsq, col_rsq,
                              factor_key = FALSE)

results_long
# 
# 
# 
plot090_all3_lasso0005 <- ggplot(results_long, aes(x=alpha, y=rsq, group=whichrsq)) +
  # geom_line(aes(color=whichrsq)) +
  geom_point(size = 0.9, aes(color = whichrsq))+
  theme(legend.position = "bottom")

plot090_all3_lasso0005 <- ggplot(results_long, aes(x=log_lasso, y=rsq, group=whichrsq)) +
  geom_point(size = 0.9, aes(color = whichrsq))+
  theme(legend.position = "bottom")


plot090_all3_lasso0005 <- ggplot(results_long, aes(x=log_lasso_y, y=rsq, group=whichrsq)) +
  geom_point(size = 0.9, aes(color = whichrsq))+
  # coord_cartesian(xlim = c(0, 0.02)) +
  theme(legend.position = "bottom")

