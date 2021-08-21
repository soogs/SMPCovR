# discovR regression #
# initiated: 15th June 2021 #
# last modified: 21st August 2021 #
# Soogeun Park #

# 0. backstory ####
# So far, the pcovr extensions for multiblock data that I have been developing have 
# been focusing on inducing the sparsity on the weights (W)
# the method would be further extended to also induce sparsity on the regression coefficients

discovR <- function(X, Y, blockcols, R, alpha, 
                    lasso_w, grouplasso_w, lasso_y, 
                    ridge_y,
                    inits = c("rational", "oracle", "multistart"), 
                    nrstart = 10, 
                    include_rational = TRUE,
                    seed,
                    MAXITER = 10000, stop_value = 1e-10){
  
  # Input #
  # X : predictor matrix
  # Y : outcome variables (can be a matrix)
  # blockcols : vector specifying the number of variables that each data block has
  # R : number of covariates
  # alpha : weighting parameter between X and Y
  # lasso_w : lasso penalty for the weights
  # grouplasso_w : group lasso penalty for the weights
  # lasso_y : lasso penalty for the regression coefficients
  # ridge_y : ridge penalty for the regression coefficients
  # inits : starting value specification
  # nrstart : number of multiple starts
  # include_rational : whether rational start is included as a set of the multistart
  # MAXITER : Maximum number of iteration
  # stop_value : tolerance
  # seed : set.seed value
  
  # 1. define sub-functions ####
  
  if (length(lasso_w) != R){
    stop("Vector of length R is required as an input for the lasso penalty")
  }
  
  # updatePx function #
  # Adopted from the Github repository for the Sparse PCovR paper (Van Deun et al., 2018)
  # https://doi.org/10.1186/s12859-018-2114-5
  updatePx <- function(wX, W, X){
    
    if (sum(colSums(W != 0) == 0) > 0){
      print ("ERROR: W matrix has zero-columns. This disallows the P computation. Try a lower l1 penalty.")
      stop()
    }
    
    Tmat <- X %*% W
    K1 <- t(wX) %*% Tmat
    K <- t(Tmat) %*% wX %*% t(wX) %*% Tmat
    eigs <- eigen(K)
    V <- eigs$vectors
    Ssq <- eigs$values
    
    S <- diag(Ssq^(-0.5))
    
    if(ncol(W) == 1){
      S <- matrix(Ssq^(-0.5), ncol = 1, nrow = 1)
    }
    
    Px <- K1 %*% V %*% S %*% t(V)
    
    return (Px)
  }
  
  # loss calculation function #
  # the inputs for this function are as raw as possible:
  # no pre-weighting them by alpha or the l2 norm of X or Y
  
  losscal <- function(Y, X, W, Px, Py, alpha, lasso_w, grouplasso_w, lasso_y, ridge_y, blockindex){
    lasso_w_mat <- matrix(0, nrow = dim(W)[1], ncol = dim(W)[2])
    
    for (r in 1:length(lasso_w)){
      lasso_w_mat[,r] <- lasso_w[r]
    }
  
    
    glasso_norm <- function(x, grouplasso_w, blockindex, R){
      l2norm <- 0
      
      for (r in 1:R){
        for (i in 1:length(blockindex)){
          l2norm <- l2norm + grouplasso_w[r] * sqrt(sum(x[blockindex[[i]],r]^2)) * sqrt(length(blockindex[[i]]))
        }
      }
      return (l2norm)
    }
    
    
    pca_loss <- sum((X - X %*% W %*% t(Px))^2) / sum((X^2))
    
    reg_loss <- sum((Y - X %*% W %*% Py)^2) / sum(Y^2)
    
    l2norm <- glasso_norm(x = W, grouplasso_w = grouplasso_w, blockindex = blockindex, R = ncol(W))
    
    result <- (alpha) * reg_loss + 
      (1 - alpha) * pca_loss + 
      sum(lasso_w_mat * abs(W)) +  
      l2norm + 
      lasso_y * sum(abs(Py)) + 
      ridge_y * sum(Py^2)
    
    return (result)
  }
  
  # pcovr function from bmc bioinformatics #
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
      Py <- matrix(data = Py, nrow = 1)
    }
    
    # fit measure #
    RsqX <- 1-sum(sum((X - X %*% W %*% t(Px))^2))/(sum(sum(X^2)))
    Rsqy <- 1-sum(sum((Y - X %*% W %*% t(Py))^2))/(sum(sum(Y^2)))
    
    return_list <- list(W = W, Px = Px, Py = Py, RsqX = RsqX, Rsqy = Rsqy)
    
    return (return_list)
  }
  
  # updatePy function #
  updatePy <- function(X, W, Y, R, Py, Px, alpha, lasso_w, ridge_w, lasso_y, ridge_y){
    
    loss_py <- function(X, Y, W, Py, lasso_y, ridge_y, alpha){
      
      result <- alpha / sum(Y^2) * sum((Y - X %*% W %*% Py)^2) + lasso_y * sum(abs(Py)) + ridge_y * sum(Py^2)
      
      return(result)
      
    }
    
    I <- nrow(X)
    Jx <- ncol(X)
    Jy <- ncol(Y)
    J <- Jx + Jy
    
    Tmat <- X %*% W
    
    ixt <- diag(ncol(Y)) %x% Tmat
    
    y <- c(Y)
    
    Py_check <- Py
    
    py <- c(Py)
    
    py_old <- py
    
    loss_old <- 100000
    
    conv <- FALSE
    
    conv_counter <- 0
    
    while (!conv){
      
      for (k in 1:(ncol(Y) * R)){
        
        rk <- y - ixt %*% py + ixt[,k] * py[k]
        
        tk <- ixt[,k]
        
        A <- (I*J*alpha * sum(tk^2)) / (sum(Y^2))
        
        left <- A * t(tk) %*% rk / (sum(tk^2))
        
        right <- lasso_y / (2)
        
        # soft thresholding #
        if (left > right){
          
          pyk <- (left - right) / (A + ridge_y)
          
        } else if (left < (-right)){
          
          pyk <- (left + right) / (A + ridge_y)
          
        } else if (abs(left) <= right){
          
          pyk <- 0
        }
        
        # update the regression coefficient #
        py[k] <- pyk 
        
        Py_check <- matrix(py, nrow = R)
        
        # loss check (the entire loss) #
        loss_new <- losscal(Y = Y, X = X, W = W, alpha = alpha, Px = Px, Py = Py_check, lasso_w = lasso_w,
                            ridge_w = ridge_w, lasso_y = lasso_y, ridge_y = ridge_y)
        
        # loss_new <- loss_py(X = X, Y = Y, W = W, Py = Py_check, lasso_y = lasso_y, ridge_y = ridge_y, alpha = alpha)
        
        
        if ((loss_new - loss_old) > 1e-10){
          print("loss up")
        } else{
          loss_old <- loss_new
        }
        
      }
      
      # convergence check #
      if (sum(abs(py_old - py)) < 1e-10){
        conv <- TRUE
      } else {
        
        py_old <- py
        
        conv_counter <- conv_counter + 1
      }
    }
    return(list(py = py, conv_counter = conv_counter))
  }
  
  
  # 2. define a few objects ####
  
  # Y could be provided as a column vector. transfer into matrix
  if (is.vector(Y)){
    Y <- matrix(data = Y, ncol = 1)
  }
  
  I <- nrow(X)
  Jx <- ncol(X)
  Jy <- ncol(Y)
  J <- Jx + Jy
  
  # weighting X and Y according to the alpha parameter and the I*J
  w1 <- (I*J*alpha) / (sum(Y^2))
  w2 <- (I*J*(1-alpha)) / (sum(X^2))
  
  wY <- sqrt(w1) * Y
  wX <- sqrt(w2) * X
  
  # Z = [wY wX]
  Z <- cbind(wY, wX)
  
  # blockindex object needed to allow for esimation (no multiblock)
  blockindex <- list(1:ncol(X))
  
  # common and distinctive structure (only common component)
  cd <- rep(1, R)
  
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
    
  # vector to save results later
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
      
      Py <- MASS::ginv(t(Tmat) %*% Tmat) %*% t(Tmat) %*% Y
      
      colsumX2 <- colSums(X^2)
      
      ssZ <- sum(Z^2)
      
      # initial loss #
      loss0 <- 10000
      
      # convergence starting 
      conv <- 0
      iter <- 1
      
      loss_hist <- c(loss0)
      
      # 5. estimation ####
      while (conv == 0){
        
        # P given W #
        Px <- updatePx(wX = wX, X = X, W = W)
        
        py <- updatePy(X = X, W = W, Y = Y, R = R, Py = Py, Px = Px, alpha = alpha, 
                       lasso_w = lasso_w, ridge_w = ridge_w, 
                       lasso_y = lasso_y, ridge_y = ridge_y)$py
        
        Py <- matrix(py, nrow = R)
        
        # reweighitng the Px for loss calculationg
        Px_weighted <- Px / (sqrt(w2))
        
        loss1 <- losscal(Y = Y, X = X, W = W, Px = Px_weighted, Py = Py, alpha = alpha,
                         lasso_w = lasso_w, ridge_w = ridge_w,
                         lasso_y = lasso_y, ridge_y = ridge_y)
        
        # current loss is always smaller than the previous loss
        if ((loss0 - loss1) < stop_value){
          conv <- 1
        }
        
        if ((loss1 - loss0) > 1e-10){
          print ("ERROR: current loss > previous loss, after loadings estimation")
          # stop()
        }
        
        
        
        iter <- iter + 1
        
        loss_hist[iter] <- loss1
        
        loss0 <- loss1
        
        
        # I think Py should be reweighted before providing as input for the updateW?
        Py_weighted <- Py * sqrt(w1)
        
        P_weighted <- rbind(t(Py_weighted), Px)
        
        colsumP2 <- colSums(P_weighted^2)
        
        # W given P #
        W <- updateW_cpp(Z = Z, W = W, X = X, P = P_weighted, R = R, lambda1 = as.matrix(lasso_w), lambda2 = ridge_w, colsumX2 = as.matrix(colsumX2), colsumP2 = as.matrix(colsumP2), blockindex = blockindex, cd = cd)
        
        # reweighting the Px for loss calculation
        Px_weighted <- Px / (sqrt(w2))
        
        loss1 <- losscal(Y = Y, X = X, W = W, Px = Px_weighted, Py = Py, alpha = alpha,
                         lasso_w = lasso_w, ridge_w = ridge_w,
                         lasso_y = lasso_y, ridge_y = ridge_y)
        
        # current loss is always smaller than the previous loss
        if ((loss0 - loss1) < stop_value){
          conv <- 1
        }
        
        if ((loss1 - loss0) > 1e-10){
          print ("ERROR: current loss > previous loss, after weights estimation")
          # stop()
        }
        iter <- iter + 1
        
        loss_hist[iter] <- loss1
        
        loss0 <- loss1
      }
      
      
      P <- rbind(t(Py), Px)
      
      # reweighting step #
      # Py <- matrix(P[1:Jy,] / (sqrt(w1)), ncol = R)
      Px <- P[-c(1:Jy),] / (sqrt(w2))
      
      if (R == 1){
        Px <- matrix(P[-c(1:Jy),] / (sqrt(w2)), ncol = 1)
      }
      
      result_list <- list(W = W, P = P, loss = loss1, loss_hist = loss_hist, iter = iter, Px = Px, Py = Py)
      multi_loss[nr] <- loss1
      multi_results[[nr]] <- result_list
      
    }
    
    lossmin <- which.min(multi_loss)
    multi_results_min <- multi_results[[lossmin]]
    
  return_results <- result_list
  return (return_results)
}
