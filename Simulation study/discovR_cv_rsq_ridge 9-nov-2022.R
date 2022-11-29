# discovR_cv_Rsq  #
# last modified: 9-nov-2022 #
# Soogeun Park # 

# cross-validation for discovR, but using Rsq measure 
# instead of the SSQ


# 9-nov-2022 update
# I check for the (a) number of outcome variables included in each fold and (b) how many of these are correct

# every other detail of the function is the same as discovR_cv_rsq_ridge in "to_blade_20_june_2022"

discovR_cv_rsq <- function(X, Y, R, alpha, lasso_w, 
                       ridge_w, lasso_y, ridge_y, nrFolds, seed_cv, seed_method, 
                       fold_scale = TRUE,
                       inits, include_rational, MAXITER, stop_value, nonzero_only, true_active_Y){
  
  
  # Input #
  # X : predictor matrix
  # Y : response variables
  # R : number of covariates
  # alpha : weighting parameter between X and Y
  # lasso_w : lasso parameter for W
  # ridge_w : ridge parameter for w
  # lasso_y : lasso parameter for Py
  # ridge_y : ridge parameter for Py
  # inits : the initialization strategy used for model-fitting. "rational", "multistart"
  # seed_cv : seed used for the CV. the manner of folds being divided is dependent on this
  # seed_method : seed used for the model-fitting. This matters if you want to do multistart; the random starting values are dependent on this
  # include_rational : if the rational start is included as part of the multistart procedure
  # MAXITER: maximum number of iterations
  
  set.seed(seed_cv)
  
  # randomly shuffle the data
  sampled <- sample(nrow(X))
  
  X <- X[sampled,]
  
  Y <- Y[sampled,]
  
  #Create equally sized folds
  folds <- cut(seq(1,nrow(X)),breaks=nrFolds,labels=FALSE)
  
  cve_k <- data.frame(error_y = NA)
  
  cve_x_k <- data.frame(error_x = NA)
  
  # vector saving the number of outcome variables included
  Y_included <- c()
  
  # vector saving the correct classification rate of outcome variables included
  Y_correct <- c()
  
  convergence <- c()
  
  for(k in 1:nrFolds){
    
    # defining test and train data  
    test_index <- which(folds==k, arr.ind=TRUE)
    
    X_train <- X[-test_index,]
    X_test <- X[test_index, ]
    
    Y_train <- Y[-test_index,]
    Y_test <- Y[test_index,]
    
    # simulation studies prior to 13 April 2022 have all been  conducted with fold_scale = TRUE
    if (fold_scale){
      # center and scale #
      X_train <- scale(X_train,T,T)
      Y_train <- scale(Y_train,T,T)
      
      X_test <- scale(X_test,T,T)
      Y_test <- scale(Y_test,T,T)
      
    } 
    
    # model fitting #
    scd1 <- discovR_cpp_ridge(X = X_train, Y = Y_train, R = R, 
                        alpha = alpha, lasso_w = lasso_w, ridge_w = ridge_w,
                        lasso_y = lasso_y, ridge_y = ridge_y, inits = inits, 
                        nrstart = 1, include_rational = TRUE, seed = seed_method, 
                        stop_value = stop_value, MAXITER = MAXITER)
    
    if (length(scd1$loss_hist) < MAXITER){
      converged <- T
    } else {
      converged <- F
    }
    
    # out of sample prediction #
    pred <- X_test %*% scd1$W %*% t(scd1$Py)
    
    # out of sample prediction for X #
    pred_x <- X_test %*% scd1$W %*% t(scd1$Px)
    
    # If we want to only compare based on the outcome variables that were predicted as nonzero:
    if (nonzero_only){
      
      nonzero_index <- colSums(abs(pred)) > 1e-5
      
      Rsq_i <- 1 - sum((Y_test[,nonzero_index] - pred[,nonzero_index])^2) / sum(Y_test[,nonzero_index]^2)
      
      cve_k[k,1] <- Rsq_i
      
      # Rsq regarding X
      nonzero_index_x <- colSums(abs(pred_x)) > 1e-5
      
      Rsq_x_i <- 1 - sum((X_test[,nonzero_index_x] - pred_x[,nonzero_index_x])^2) / sum(X_test[,nonzero_index_x]^2)
      
      cve_x_k[k,1] <- Rsq_x_i
    } else {
      
      # Rsq for Y
      Rsq_i <- 1 - sum((Y_test - pred)^2) / sum(Y_test^2)
      
      cve_k[k,1] <- Rsq_i
      
      # Rsq for X
      nonzero_index_x <- colSums(abs(pred_x)) > 1e-5
      
      Rsq_x_i <- 1 - sum((X_test - pred_x)^2) / sum(X_test^2)
      
      cve_x_k[k,1] <- Rsq_x_i
      
    }
    
    # checking the number of outcome variables included in each fold
    Y_included[k] <- sum(colSums(abs(pred)) > 1e-5)
    
    Y_included_index <- colSums(abs(pred)) > 1e-5
    
    Y_correct[k] <- (sum((abs(Y_included_index) > 1e-7) + (abs(true_active_Y) > 1e-7) == 2) + 
                sum((abs(Y_included_index) < 1e-7) + (abs(true_active_Y) < 1e-7) == 2)) / (length(true_active_Y))
    
    convergence[k] <- converged
    
    
  }
  
  cve <- colMeans(cve_k)
  # cve = mean cv error
  
  # cve for X
  cve_x <- colMeans(cve_x_k)
  
  se <- apply(cve_k,2,sd) / sqrt(nrFolds)
  # se = standard deviation of the mean cv errors per k / sqrt(nrFold)
  
  # se for X
  se_x <- apply(cve_x_k,2,sd) / sqrt(nrFolds)
  
  # we check how many outcome variables are predicted as 0
  disc_full <- discovR_cpp_ridge(X = X, Y = Y, R = R, 
                           alpha = alpha, lasso_w = lasso_w, ridge_w = ridge_w,
                           lasso_y = lasso_y, ridge_y = ridge_y, inits = inits, 
                           nrstart = 1, include_rational = TRUE, seed = seed_method, stop_value = stop_value, 
                           MAXITER = MAXITER)
  
  zero_Y <- sum(colSums(abs(X %*% disc_full$W %*% t(disc_full$Py))) == 0)
  
  zero_W <- sum(colSums(abs(disc_full$W)) == 0)
  
  zero_Py <- sum(colSums(abs(disc_full$Py)) == 0)
  
  # also, now that we just fitted the model, we can compute the in-sample Rsq
  fit <- X %*% disc_full$W %*% t(disc_full$Py)
  
  fit_Rsq_all <- 1 - sum((Y - fit)^2) / sum(Y^2)
  
  fit_nonzero_index <- colSums(abs(fit)) > 1e-5
  
  if (sum(fit_nonzero_index) == 0){
    fit_Rsq_nonzero <- 0 
  } else{
    
    fit_Rsq_nonzero <- 1 - sum((Y[,fit_nonzero_index] - fit[,fit_nonzero_index])^2) / sum(Y[,fit_nonzero_index]^2)
    
  }
  
  
  # we can also compute the in-sample Rsq in relation to the predictors X
  fit_x <- X %*% disc_full$W %*% t(disc_full$Px)
  
  Rsq_all_x <- 1 - sum((X - fit_x)^2) / sum(X^2)
  
  fit_nonzero_index_x <- colSums(abs(fit_x)) > 1e-5
  
  if (sum(fit_nonzero_index_x) == 0){
    Rsq_nonzero_x <- 0 
  } else{
    
    Rsq_nonzero_x <- 1 - sum((X[,fit_nonzero_index_x] - fit_x[,fit_nonzero_index_x])^2) / sum(X[,fit_nonzero_index_x]^2)
    
  }
  
  
  return(list(cve = cve, se = se, cve_k = cve_k, 
              cve_x = cve_x, se_x = se_x, cve_x_k = cve_x_k,
              convergence = convergence,
              Y_included = Y_included, Y_correct = Y_correct,
              zero_Y = zero_Y, zero_W = zero_W, zero_Py = zero_Py,
              fit_Rsq_all = fit_Rsq_all, fit_Rsq_nonzero = fit_Rsq_nonzero,
              fit_Rsq_all_x = Rsq_all_x, fit_Rsq_nonzero_x = Rsq_nonzero_x))
  

}

