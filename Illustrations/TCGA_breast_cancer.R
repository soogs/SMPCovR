# discovR illustration #
# initiated: 22nd August 2021 #
# last modified: 22nd August 2021 #
# Soogeun Park #

# I will demonstrate the discovR function (regression problem)
# with the breast cancer data, from the mixOmics package

# 1. Data load and pre-process ####

library(mixOmics)
library(RegularizedSCA)

data('breast.TCGA')
# extract training data
data = list(mRNA = breast.TCGA$data.train$mrna, 
            miRNA = breast.TCGA$data.train$mirna, 
            proteomics = breast.TCGA$data.train$protein)

# check dimension
lapply(data, dim) # these are the three data blocks

# let's take the first five variables from the proteomics block as our outcome variables
Y <- data$proteomics[,1:5]

data$proteomics <- data$proteomics[,-(1:5)]

lapply(data, dim)

# concatenating the data, and center-scale each block, such that 
# the sum of squares for each block is the same for all the blocks
# this is done with the pre_process function from the RegularizedSCA package

dat <- cbind(pre_process(data$mRNA, weight = T), 
      pre_process(data$miRNA, weight = T), 
      pre_process(data$proteomics, weight = T))

# also scaling the outcome variables:
Y <- pre_process(Y, weight = T)

# 2. Model selection ####


# 3. discovR estimation ####

time1 <- Sys.time()

hi <- discovR(X = dat, Y = Y, blockcols = c(200, 184, 137), R = 2, 
              alpha = 0.4, lasso_w = c(0.002,0.002), 
              grouplasso_w = c(0.0002, 0.0002),
              lasso_y = c(0.0001, 0.0001), 
              ridge_y = 0.0001, inits = "rational", 
              nrstart = 1, seed = 1111, MAXITER = 10000, w_one_iteration = TRUE, stop_value = 1e-5)

time2 <- Sys.time()

time3 <- Sys.time()

hi2 <- discovR_cpp(X = dat, Y = Y, blockcols = c(200, 184, 137), R = 2, 
              alpha = 0.4, lasso_w = c(0.002,0.002), 
              grouplasso_w = c(0.0002, 0.0002),
              lasso_y = c(0.0001, 0.0001), 
              ridge_y = 0.0001, inits = "rational", 
              nrstart = 1, seed = 1111, MAXITER = 10000, stop_value = 1e-5)

time4 <- Sys.time()


hi <- discovR(X = dat, Y = Y, blockcols = c(200, 184, 137), R = 1, 
              alpha = 0.4, lasso_w = c(0.0002), 
              grouplasso_w = c(0.00002),
              lasso_y = c(0.0001), 
              ridge_y = 0.0001, inits = "rational", 
              nrstart = 1, seed = 1111, MAXITER = 10000, w_one_iteration = TRUE, stop_value = 1e-5)


