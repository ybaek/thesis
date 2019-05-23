#
# Constants for the model
# (Calculation remote)
#

library(abind)
library(dplyr)
library(geosphere)
library(Matrix)
library(parallel)

source("000-utilities.R")
load("data/ReducedData.RData")
# Global variables
data_grid <- real_data
time_seq <- seq(as.POSIXct("2015-01-01 00:00:00", tz="UTC"),
	        as.POSIXct("2017-12-31 23:00:00", tz="UTC"), by=3600)
N <- length(unique(data_grid$id)) 
U <- length(time_seq) ## 26304
P <- 15
R <- 24

# Parallelizing: split and stack matrix rows
nc <- detectCores()
cl <- makeCluster(nc, type="FORK")
indsList <- parallel::splitIndices(U, nc)

# Computing parameters (matrix algebra on global variables)
computeParam <- function(inds) {
  len <- max(inds) - min(inds) + 1 # dimension of the arrays
  y <- matrix(NA, N, len)
  X <- array(NA, c(N, P, len))
  C <- array(NA, c(N, N, len))
  S <- array(NA, c(N, R, len))
  Phi <- array(NA, c(R, R, len))
  etaSigma <- array(NA, c(R, R, len))
  xiSigma <- array(NA, c(R, R, len))
  
  for (t in inds) {
    # 1. Components from the dataset (y, X, C, Q)
    dataT <- data_grid[data_grid$dt==time_seq[t],]
    i <- t - min(inds) + 1
    y[,i] <- dataT$log_dust
    Xmat <- as.matrix(dataT[,c("temp","ws","humid","pres","elevation")])
    Xmat <- scale(Xmat, scale=FALSE)
    Xintmat <- apply(combn(colnames(Xmat)[-5], 2), 2,
                     function(names) Xmat[,names[1]] * Xmat[,names[2]])
    X[,,i] <- cbind(Xmat, Xintmat, Xmat[,-5]^2)
    
    # 2. Similarity matrix and its Laplacian
    distT <- normDist_fn(distMat_fn(as.matrix(dataT[,10:11])))
    C[,,i] <- 1-distT
    C[,,i] <- ifelse(C[,,i] > .5, C[,,i], 0) # temp. cutoff
    diag(C[,,i]) <- 0
  
    # 3. S (basis functions)
    E1 <- moran_fn(C[,,i])
    S[,,i] <- E1$vectors[,!dplyr::near(E1$values, 0)]
    
    # 4. etaSigma (covariance matrix)
    etaSigma[,,i] <- as.matrix(nearPD(t(S[,,i]) %*% (diag(rowSums(C[,,i])) - C[,,i]) %*%
    		     				S[,,i])$mat)
    
    # 5. Phi (transition matrix) 
    if (t>1) {  
      P <- diag(1,R) - t(S[,,i]) %*% hat_fn(X[,,i]) %*% S[,,i]
      E2 <- eigen(P %*% diag(1,R) %*% P)
      Phi[,,i] <-  E2$vectors
    }
  }
  return(list(y=y, X=X, C=C, S=S, Phi=Phi, etaSigma=etaSigma))
}
paramList <- clusterApply(cl, indsList, computeParam)
stopCluster(cl) # Make sure to stop (manually)

param_bind <- function(str, axis) {
  my_bind <- function(a,b) {
    result <- list()
    result[[str]] <- abind(a[[str]], b[[str]], along=axis)
    return(result)
  }
  Reduce(my_bind, paramList)
}


y <- param_bind("y",2)$y
X <- param_bind("X",3)$X
C <- param_bind("C",3)$C
S <- param_bind("S",3)$S
etaSigma <- param_bind("etaSigma", 3)$etaSigma
Phi <- param_bind("Phi", 3)$Phi[,,-1] # Except for 1st t

xiSigma <- array(NA, dim=c(R,R,(U-1)))
for (t in 2:U) {
  xiSigma[,,(t-1)] <- as.matrix(nearPD(etaSigma[,,t] - Phi[,,(t-1)] %*% etaSigma[,,(t-1)] %*% t(Phi[,,(t-1)]))$mat) 
}

save(list=c("y", "X", "C", "S", "Phi", "etaSigma", "xiSigma", "N", "P", "R", "U"), file="modelConstants.RData")
