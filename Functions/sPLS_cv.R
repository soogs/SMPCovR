# sPLS-cv-loo #
# initiated: 19 April 2022 #
# Soogeun Park # 

# cross validation for discovR - but for doing Leave-one-out #



sPLS_cv <- function(plsX, Y, ncomp, nrFolds, keepX, design, seed_cv, fold_scale = FALSE){
  
  set.seed(seed_cv)
  
  # randomly shuffle the data
  sampled <- sample(nrow(plsX$block1))
  
  X <- plsX$block1[sampled,]
  
  Y <- Y[sampled,]
  
  #Create equally sized folds
  folds <- cut(seq(1,nrow(X)),breaks=nrFolds,labels=FALSE)
  
  cve_k <- data.frame(error_y = NA)
  
  convergence <- c()
  
  for(k in 1:nrFolds){
    
    # defining test and train data  
    test_index <- which(folds==k, arr.ind=TRUE)
    
    X_train <- X[-test_index,]
    X_test <- X[test_index, ]
    
    plsX_train <- list(block1 = X_train)
    plsX_test <- list(block1 = X_test)
    
    
    Y_train <- Y[-test_index,]
    Y_test <- Y[test_index,]
    
    
    rownames(Y_train) <- rownames(plsX_train$block1)
    # rownames(Y_test) <- rownames(plsX_test$block1)
    
    
    # simulation studies prior to 13 April 2022 have all been  conducted with fold_scale = TRUE
    if (fold_scale){
      # center and scale #
      X_train <- scale(X_train,T,T)
      Y_train <- scale(Y_train,T,T)
      
      X_test <- scale(X_test,T,T)
      Y_test <- scale(Y_test,T,T)
      
    } 
    
    dia1 <- block.spls(X = plsX_train, Y = Y_train,
               ncomp = R, keepX = keepX, design = design)
    
    dia_pred <- predict(object = dia1, newdata = plsX_test)
    
    Rsq_i <- 1 - sum((dia_pred$predict$block1[,,ncomp] - Y_test)^2) / sum(Y_test^2)
    
    cve_k[k,1] <- Rsq_i
    
    # convergence[k] <- converged
    
    
  }
  
  cve <- colMeans(cve_k)
  # cve = mean cv error
  
  se <- apply(cve_k,2,sd) / sqrt(nrFolds)
  # se = standard deviation of the mean cv errors per k / sqrt(nrFold)
  
  return(list(cve = cve, se = se, cve_k = cve_k))
}

