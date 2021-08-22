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
  # lasso_w : lasso penalty for the weights. vector for multiple components
  # grouplasso_w : group lasso penalty for the weights. vector for multiple components
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
  
  
  # loss calculation function #
  # the inputs for this function are as raw as possible:
  # no pre-weighting them by alpha or the l2 norm of X or Y
  
  losscal <- function(Y, X, W, Px, Py, alpha, lasso_w, grouplasso_w, lasso_y, ridge_y, blockindex){
    lasso_w_mat <- matrix(0, nrow = dim(W)[1], ncol = dim(W)[2])
    
    for (r in 1:length(lasso_w)){
      lasso_w_mat[,r] <- lasso_w[r]
    }
    
    
    lasso_y_mat <- matrix(0, nrow = dim(Py)[1], ncol = dim(Py)[2])
    
    for (r in 1:length(lasso_y)){
      lasso_y_mat[,r] <- lasso_y[r]
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
    
    reg_loss <- sum((Y - X %*% W %*% t(Py))^2) / sum(Y^2)
    
    l2norm <- glasso_norm(x = W, grouplasso_w = grouplasso_w, blockindex = blockindex, R = ncol(W))
    
    result <- (alpha) * reg_loss + 
      (1 - alpha) * pca_loss + 
      sum(lasso_w_mat * abs(W)) +  
      l2norm + 
      sum(lasso_y_mat * abs(Py)) + 
      ridge_y * sum(Py^2)
    
    return (result)
  }
  
  # pcovr function from Katrijn's bmc bioinformatics paper #
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
  updatePy <- function(X, Y, R, W, Py, Px, alpha, 
                       lasso_w, grouplasso_w, lasso_y, ridge_y,
                       blockindex){
    
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
                              grouplasso_w = grouplasso_w, 
                              lasso_y = lasso_y, ridge_y = ridge_y,
                              blockindex = blockindex)
          
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
  
  
  # updateW function #
  updateW <- function(Y, X, R, W, Py, Px, alpha, 
                      lasso_w, grouplasso_w, lasso_y, 
                      ridge_y, blockindex){
    
    
    # soft-thresholding operator
    soft <- function(x, lambda){
      x2 <- abs(x) - lambda
      
      if(x2 < 0){x2 <- 0}
      
      result <- sign(x) * x2
      
      return(result)
    }
    
    # general loss calculation #
    
    losscal <- function(Y, X, W, Px, Py, alpha, lasso_w, grouplasso_w, lasso_y, ridge_y, blockindex){
      lasso_w_mat <- matrix(0, nrow = dim(W)[1], ncol = dim(W)[2])
      
      for (r in 1:length(lasso_w)){
        lasso_w_mat[,r] <- lasso_w[r]
      }
      
      
      lasso_y_mat <- matrix(0, nrow = dim(Py)[1], ncol = dim(Py)[2])
      
      for (r in 1:length(lasso_y)){
        lasso_y_mat[,r] <- lasso_y[r]
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
      
      reg_loss <- sum((Y - X %*% W %*% t(Py))^2) / sum(Y^2)
      
      l2norm <- glasso_norm(x = W, grouplasso_w = grouplasso_w, blockindex = blockindex, R = ncol(W))
      
      result <- (alpha) * reg_loss + 
        (1 - alpha) * pca_loss + 
        sum(lasso_w_mat * abs(W)) +  
        l2norm + 
        sum(lasso_y_mat * abs(Py)) + 
        ridge_y * sum(Py^2)
      
      return (result)
    }
    
    SSY <- sum(Y^2)
    
    SSX <- sum(X^2)
    
    # W_new will be the new matrix
    W_new <- W
    
    # initial loss calculation #
    loss0 <- losscal(Y = Y, X = X, W = W, Py = Py, 
                     Px = Px, alpha = alpha, lasso_w = lasso_w,
                     grouplasso_w = grouplasso_w, 
                     lasso_y = lasso_y,
                     ridge_y = ridge_y, blockindex = blockindex)
    
    # saving loss history #
    loss_hist <- loss0
    
    # set.seed(seed)
    
    # randomizing the order of coordinate descent
    r_order <- sample(1:R)
    # r_order <- 1:R
    
    conv <- FALSE
    
    loss0_out <- loss0
    
    while (!conv){
      
      for (r in r_order){
        
        k_order <- sample(1:length(blockindex))
        # k_order <- 1:length(blockindex)
        
        for (k in k_order){
          
          Wk <- c(W_new[blockindex[[k]],r])
          
          Xk <- X[,blockindex[[k]]]
          
          # r_k: contribution for Y, except group k & component r
          r_k <- Y - X %*% W_new %*% t(Py) + Xk %*% Wk %*% t(Py[,r])
          
          # s_k: contribution for X, except for group k & component r
          s_k <- X %*% Px[,r] - (X %*% W_new[,r] - Xk %*% Wk)
          
          soft_y <- (2*alpha/SSY) * (r_k %*% Py[,r])
          
          soft_x <- (2 *(1-alpha)/SSX) * s_k
          
          k_soft <- diag(c(soft_y + soft_x)) %*% Xk
          
          k_soft <- colSums(k_soft)
          
          k_soft <- unlist(lapply(k_soft, function(x){soft(x, lasso_w[r])}))
          
          # check if l2-norm of k_soft is smaller than grouplasso penalty #
          J_k <- length(blockindex[[k]])
          
          if(sqrt(sum(k_soft^2)) <= (grouplasso_w[r] * sqrt(J_k))){
            Wk_new <- Wk
            
            Wk_new[] <- 0
            
            # update the entire matrix W_new
            W_new[blockindex[[k]],r] <- Wk_new
            
            # loss W checking #
            loss1 <- losscal(Y = Y, X = X, W = W_new, Py = Py, 
                             Px = Px, alpha = alpha, lasso_w = lasso_w,
                             grouplasso_w = grouplasso_w, lasso_y = lasso_y,
                             ridge_y = ridge_y, blockindex = blockindex)
            
            if (loss1 > loss0){
              print("loss increase after block sparsified")
              stop()
            }
            
            # saving loss history #
            loss_hist <- append(loss_hist, loss1)
            
            loss0 <- loss1
            
          } else {
            # if the group coefficient is not zero vector,
            # iteration of coordinate descent for elementwise sparsity
            
            h_order <- sample(1:length(blockindex[[k]]))
            
            # saving objects prior to loop
            PX <- Px %x% X 
            
            Xk <- X[,blockindex[[k]]]
            
            Xl <- X[,-blockindex[[k]]]
            
            conv_ingroup <- FALSE
            
            Wk <- c(W_new[blockindex[[k]],r])
            
            for (h in h_order){
              
              r_kh <- Y - X %*% W_new %*% t(Py) + (matrix(Xk[,h] * Wk[h], ncol = 1) %*% t(Py[,r]))
              
              s_kh <- X %*% Px[,r] - X %*% W_new[,r] + matrix(Xk[,h] * Wk[h], ncol = 1)
              
              soft_y_h <- (2*alpha/SSY) * (r_kh %*% Py[,r])
              
              soft_x_h <- (2 *(1-alpha)/SSX) * s_kh
              
              kh_soft <- t(soft_y_h + soft_x_h) %*% Xk[,h]
              
              soft_Wkh <- soft(x = kh_soft, lambda = lasso_w[r])
              
              denom <- (((2 * alpha / SSY) * (t(Py[,r]) %*% Py[,r]) + 
                           (2 * (1-alpha) / SSX)) * (t(Xk[,h]) %*% Xk[,h])) +
                (grouplasso_w[r] * sqrt(J_k) / sqrt(sum(Wk^2)))
              
              Wkh_new <- soft_Wkh / denom
              
              # Wk_new is the candidate Wk vectors
              Wk_new <- Wk
              Wk_new[h] <- Wkh_new
              
              # updating W_new matrix
              W_new[blockindex[[k]],r] <- Wk_new
              
              # loss calculation #
              loss1 <- losscal(Y = Y, X = X, W = W_new, Py = Py, 
                               Px = Px, alpha = alpha, lasso_w = lasso_w,
                               lasso_y = lasso_y,
                               grouplasso_w = grouplasso_w, ridge_y = ridge_y, blockindex = blockindex)
              
              if (loss1 > loss0){
                print("loss increase via coordinate descent")
                
                stop()
              }
              
              
              if (Wkh_new != 0){
                if (loss0 - loss1 < 1e-7){
                  conv_ingroup <- TRUE
                }
              }
              
              # saving loss history #
              loss_hist <- append(loss_hist, loss1)
              
              loss0 <- loss1
              
              Wk <- c(W_new[blockindex[[k]], r])
              
            } 
          }
        }
      }
      
      # one loop complete. check loss for convergence 
      loss1_out <- losscal(Y = Y, X = X, W = W_new, Py = Py, 
                           Px = Px, alpha = alpha, lasso_w = lasso_w,
                           lasso_y = lasso_y,
                           grouplasso_w = grouplasso_w, 
                           ridge_y = ridge_y, blockindex = blockindex)
      
      if (loss1_out - loss0_out > 0){
        print("loss increase after one whole loop")
        stop()
      }
      
      if (abs(loss0_out - loss1_out) < 1e-7){
        conv <- TRUE
      } else{
        
        loss0_out <- loss1_out
        
      }
      
      
    }
    
    result <- list(W = W_new, loss_hist = loss_hist)
    
    return(result)
    
  }
  
  
  # updatePx function #
  updatePx <- function(alpha, X, W){
    
    if (sum(colSums(W != 0) == 0) > 0){
      print ("ERROR: W matrix has zero-columns. This disallows the P computation. Try a lower l1 penalty.")
      stop()
    }
    
    constant <- sqrt((1 - alpha) / sum(X^2))
    
    svdd <- svd(t(constant * X %*% W) %*% (constant * X))
    
    Px_new <- svdd$v %*% t(svdd$u)
    
    return(Px_new)
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
