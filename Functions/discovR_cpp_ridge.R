# discovR regression - cpp #
# initiated: 26th May 2022 #
# Soogeun Park #

# 0. backstory ####
# Previously, this discovR function was also with the grouplasso penalty on W
# but that is now replaced with ridge penalty

discovR_cpp_ridge <- function(X, Y, R, alpha, 
                        lasso_w, ridge_w, 
                        lasso_y, ridge_y,
                        inits = c("rational", "oracle", "multistart"), 
                        nrstart = 10, 
                        include_rational = TRUE,
                        seed,
                        MAXITER = 10000, stop_value = 1e-10){
  
  # Input #
  # X : predictor matrix
  # Y : outcome variables (can be a matrix)
  # R : number of covariates
  # alpha : weighting parameter between X and Y
  # lasso_w : lasso penalty for the weights. vector for multiple components
  # ridge_w : ridge penalty for the weights
  # lasso_y : lasso penalty for the regression coefficients
  # ridge_y : ridge penalty for the regression coefficients
  # inits : starting value specification
  # nrstart : number of multiple starts
  # include_rational : whether rational start is included as a set of the multistart
  # MAXITER : Maximum number of iteration
  # stop_value : tolerance
  # seed : set.seed value
  
  # 1. define sub-functions ####
  
  print("Please sourceCpp the updateW function.")
  
  if (length(lasso_w) != R){
    stop("Vector of length R is required as an input for the lasso penalty")
  }
  
  
  # loss calculation function #
  # the inputs for this function are as raw as possible:
  # no pre-weighting them by alpha or the l2 norm of X or Y
  
  losscal <- function(Y, X, W, Px, Py, alpha, lasso_w, ridge_w, lasso_y, ridge_y){
    lasso_w_mat <- matrix(0, nrow = dim(W)[1], ncol = dim(W)[2])
    
    for (r in 1:length(lasso_w)){
      lasso_w_mat[,r] <- lasso_w[r]
    }
    
    lasso_y_mat <- matrix(0, nrow = dim(Py)[1], ncol = dim(Py)[2])
    
    for (r in 1:length(lasso_y)){
      lasso_y_mat[,r] <- lasso_y[r]
    }
    
    pca_loss <- sum((X - X %*% W %*% t(Px))^2) / sum((X^2))
    
    reg_loss <- sum((Y - X %*% W %*% t(Py))^2) / sum(Y^2)
    
    result <- (alpha) * reg_loss + 
      (1 - alpha) * pca_loss + 
      sum(lasso_w_mat * abs(W)) +  
      ridge_w * sum(W^2) +
      sum(lasso_y_mat * abs(Py)) + 
      ridge_y * sum(Py^2)
    
    return (result)
  }
  
  # pcovr function from Katrijn Van Deun's bmc bioinformatics paper #
  pcovr <- function(X, Y, R, alpha){
    
    # if Y is provided as a vector.. #
    if (is.vector(Y)){
      Y <- matrix(data = Y, ncol = 1)
    }
    
    I <- nrow (X)
    Jx <- ncol (X)
    Jy <- ncol (Y)
    J <- Jx + Jy
    eps <- 1e-12
    
    # [Iy,Jy]=size(Y);
    # if Iy~=I, disp(' size Y and X do not match ');return;end;
    
    # weighting X and Y according to the alpha parameter
    w1 <- (I*J*alpha) / (sum(Y^2))
    w2 <- (I*J*(1-alpha)) / (sum(X^2))
    
    wY <- sqrt(w1) * Y
    wX <- sqrt(w2) * X
    
    # Z = [wY wX]
    Xstar <- cbind(wY, wX)
    
    # SVD data #
    if (J > I){
      XX <- X %*% t(X)
      eigs <- eigen(XX)
      Ux <- eigs$vectors
      Ssq <- eigs$values
      Sxsq <-Ssq
      Sxsq[Sxsq < 1e-6] <- 1e-6
      Sx <- sqrt(Sxsq)
      
      invSx <- Sx^(-1)
      
      Vx <- t(X) %*% Ux %*% diag(invSx)
      
      Uxstar <- t(Ux) %*% Xstar
      
      svd1 <- svd(Uxstar)
      U <- svd1$u[,1:R]
      S <- diag(svd1$d[1:R])
      V <- svd1$v[,1:R]
      
      
      if (R == 1){
        U <- matrix(svd1$u[,1:R], ncol = 1)
        S <- matrix(svd1$d[1:R], nrow = 1, ncol = 1)
        V <- matrix(svd1$v[,1:R], ncol = 1)
      }
      
    } else if (I >= J){
      XX <- t(X) %*% X
      eigs <- eigen(XX)
      
      Vx <- eigs$vectors
      Ssq <- eigs$values
      Sxsq <- Ssq
      Sxsq[Sxsq < 1e-6] <- 1e-6
      
      Sx <- sqrt(Sxsq)
      
      invSx <- Sx^(-1)
      
      Ux <- X %*% Vx %*% diag(invSx)
      Uxstar <- t(Ux) %*% Xstar
      
      svd1 <- svd(Uxstar)
      U <- svd1$u[,1:R]
      S <- diag(svd1$d[1:R])
      V <- svd1$v[,1:R]
      
      if (R == 1){
        U <- matrix(svd1$u[,1:R], ncol = 1)
        S <- matrix(svd1$d[1:R], nrow = 1, ncol = 1)
        V <- matrix(svd1$v[,1:R], ncol = 1)
      }
    }
    
    P <- V
    
    W <- Vx %*% diag(invSx) %*% U %*% S
    
    # reweighting step #
    Py <- P[1:Jy,] / (sqrt(w1))
    Px <- P[-c(1:Jy),] / (sqrt(w2))
    
    # if Py turns out to be a vector,
    # we need to make it into a matrix of one row 
    if (is.vector(Py)){
      Py <- matrix(data = Py, ncol = 1)
    }
    
    # fit measure #
    RsqX <- 1-sum(sum((X - X %*% W %*% t(Px))^2))/(sum(sum(X^2)))
    Rsqy <- 1-sum(sum((Y - X %*% W %*% t(Py))^2))/(sum(sum(Y^2)))
    
    return_list <- list(W = W, Px = Px, Py = Py, RsqX = RsqX, Rsqy = Rsqy)
    
    return (return_list)
  }
  
  # updatePy function #
  updatePy <- function(X, Y, R, W, Py, Px, alpha, 
                       lasso_w, ridge_w, lasso_y, ridge_y){
    
    loss_py <- function(X, Y, W, Py, lasso_y, ridge_y, alpha){
      
      lasso_y_mat <- matrix(0, nrow = dim(Py)[1], ncol = dim(Py)[2])
      
      for (r in 1:length(lasso_y)){
        lasso_y_mat[,r] <- lasso_y[r]
      }
      
      result <- alpha / sum(Y^2) * sum((Y - X %*% W %*% t(Py))^2) + 
        sum(lasso_y_mat * abs(Py)) + ridge_y * sum(Py^2)
      
      return(result)
      
    }
    
    Tmat <- X %*% W
    
    Py_check <- Py
    
    Py_old <- Py
    
    loss_old <- 100000
    
    loss_Py_old <- 100000
    
    conv <- FALSE
    
    conv_counter <- 0
    
    while (!conv){
      
      for (r in 1:R){
        
        for (h in 1:ncol(Y)){
          t_hr <- Y[,h] - Tmat %*% Py[h,] +  Tmat[,r] * Py[h,r] 
          
          left <- t(Tmat[,r]) %*% t_hr
          
          right <- sum(Y^2) * lasso_y[r] / (2 * alpha)
          
          denom <- sum(Tmat[,r]^2) + (sum(Y^2) / alpha) * ridge_y
          
          # soft thresholding #
          if (left > right){
            
            py_hr <- (left - right) / denom
            
          } else if (left < (-right)){
            
            py_hr <- (left + right) / denom
            
          } else if (abs(left) <= right){
            
            py_hr <- 0
          }
          
          
          # update the regression coefficient #
          Py[h,r] <- py_hr
          
          # loss check (the entire loss) #
          loss_new <- losscal(Y = Y, X = X, W = W, alpha = alpha, 
                              Px = Px, Py = Py, lasso_w = lasso_w,
                              ridge_w = ridge_w,
                              lasso_y = lasso_y, ridge_y = ridge_y)
          
          # loss check (only wrt Py)
          loss_Py_new <- loss_py(X = X, Y = Y, W = W, Py = Py, 
                                 lasso_y = lasso_y, 
                                 ridge_y = ridge_y, alpha = alpha)
          
          
          if ((loss_new - loss_old) > 1e-13){
            print("entire loss up")
            
            stop()
          } else{
            loss_old <- loss_new
          }
          
          
          if ((loss_Py_new - loss_Py_old) > 1e-13){
            print("loss wrt Py up")
            
            stop()
          } else{
            loss_Py_old <- loss_Py_new
          }
          
        }
        
      }
      
      
      # convergence check #
      if (sum(abs(Py_old - Py)) < 1e-12){
        conv <- TRUE
      } else {
        
        Py_old <- Py
        
        conv_counter <- conv_counter + 1
      }
      
    }
    return(list(Py = Py, conv_counter = conv_counter))
  }
  
  # updatePx function #
  updatePx <- function(alpha, X, W){
    
    constant <- sqrt((1 - alpha) / sum(X^2))
    
    svdd <- svd(t(constant * X %*% W) %*% (constant * X))
    
    Px_new <- svdd$v %*% t(svdd$u)
    
    return(Px_new)
  }
  
  
  # 2. define a few objects ####
  
  # Y could be provided as a column vector (in a univariate case)
  # turn it into a matrix object in that case
  if (is.vector(Y)){
    Y <- matrix(data = Y, ncol = 1)
  }
  
  # 3. initial values generated ####
  
  if (inits == "rational"){
    # W and P from pcovr
    pcovr_results <- pcovr(X = X, Y = Y, R = R, alpha = alpha)
    
    W <- pcovr_results$W
    
    nrstart <- 1
  } 
  
  if (inits == "oracle"){
    # user-specified matrices for W and P
    nrstart <- 1
  } 
  
  # vector and list to save results from multiple starting values
  multi_results <- list()
  multi_loss <- c()
  
  for (nr in 1:nrstart){
    
    if (inits == "multistart"){
      
      if (nr == nrstart & include_rational){
        # the last set of the multistart approach is always rational start
        
        pcovr_results <- pcovr(X = X, Y = Y, R = R, alpha = alpha)
        
        W <- pcovr_results$W
        
      } else {
        # for the other sets of starting values, initial values are random
        
        set.seed(seed+nr-1)
        
        W <- matrix(stats::runif(n = ncol(X)*R, min = -1, max = 1), nrow = ncol(X), ncol = R)
        
      }
      
    }
    
    # initial Py #
    Tmat <- X %*% W
    
    Py <- t(MASS::ginv(t(Tmat) %*% Tmat) %*% t(Tmat) %*% Y)
    
    colsumX2 <- colSums(X^2)
    
    # initial loss #
    loss0 <- 10000
    
    # convergence starting 
    conv <- 0
    iter <- 1
    
    loss_hist <- c(loss0)
    
    # 5. estimation ####
    while (conv == 0){
      
      # Px update #
      Px <- updatePx(alpha = alpha, X = X, W = W)
      
      # no loss calculation after Px, because I'm very confident that the loss does not increase
      # # loss check #
      # loss1 <- losscal(Y = Y, X = X, W = W, Px = Px, Py = Py, alpha = alpha,
      #                  lasso_w = lasso_w, grouplasso_w = grouplasso_w,
      #                  lasso_y = lasso_y, ridge_y = ridge_y, blockindex = blockindex)
      # 
      # # current loss is always smaller than the previous loss
      # if ((loss1 - loss0) > 1e-10){
      #   print ("ERROR: current loss > previous loss, after Px estimation")
      #   stop()
      # }
      # 
      # # loss update
      # loss0 <- loss1
      
      # Py update #
      Py <- updatePy(X = X, Y = Y, R = R, W = W, 
                     Py = Py, Px = Px, alpha = alpha, 
                     lasso_w = lasso_w, ridge_w = ridge_w, 
                     lasso_y = lasso_y, ridge_y = ridge_y)$Py
      
      loss1 <- losscal(Y = Y, X = X, W = W, Px = Px, Py = Py, 
                       alpha = alpha,
                       lasso_w = lasso_w, ridge_w = ridge_w, 
                       lasso_y = lasso_y, ridge_y = ridge_y)
      
      # current loss is always smaller than the previous loss
      if ((loss0 - loss1) < stop_value){
        conv <- 1
      }
      
      if ((loss1 - loss0) > 1e-10){
        print ("ERROR: current loss > previous loss, after Py estimation")
        stop()
      }
      
      iter <- iter + 1
      
      loss_hist[iter] <- loss1
      
      loss0 <- loss1
      
      
      # W update #
      W <- updateW_cpp_ridge(X = X, Y = Y, W = W, Px = Px, Py = Py, R = R, alpha = alpha, 
                       lasso_w = lasso_w, ridge_w = ridge_w,
                       lasso_y = lasso_y, ridge_y = ridge_y)$W_new
      
      loss1 <- losscal(Y = Y, X = X, W = W, Px = Px, Py = Py, 
                       alpha = alpha,
                       lasso_w = lasso_w, ridge_w = ridge_w,
                       lasso_y = lasso_y, ridge_y = ridge_y)
      
      # current loss is always smaller than the previous loss
      if ((loss0 - loss1) < stop_value){
        conv <- 1
      }
      
      if ((loss1 - loss0) > 1e-10){
        print ("ERROR: current loss > previous loss, after weights estimation")
        stop()
      }
      iter <- iter + 1
      
      if (iter > MAXITER){
        conv <- 1
      }
      
      loss_hist[iter] <- loss1
      
      loss0 <- loss1
    }
    
    
    # providing names for the objects
    rownames(W) <- colnames(X)
    rownames(Px) <- colnames(X)
    rownames(Py) <- colnames(Y)
    
    result_list <- list(W = W, Py = Py, Px = Px, 
                        loss = loss1, loss_hist = loss_hist, 
                        iter = iter)
    
    multi_loss[nr] <- loss1
    multi_results[[nr]] <- result_list
    
  }
  
  lossmin <- which.min(multi_loss)
  multi_results_min <- multi_results[[lossmin]]
  
  return_results <- multi_results_min
  return (return_results)
}
