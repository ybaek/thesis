
# Constants and samplers are sourced from a separate file
load("data/modelConstants.RData")
source("code/sampler_code.R")
source("code/000-utilities.R")
library(Matrix)

start <- dt_search_fn("2015-05-16 00:00:00")
end <- dt_search_fn("2015-05-20 23:00:00")
U <- end - start + 1
N <- 24
y <- y[-25,start:end]
X <- X[-25,,start:end] # intercept
C <- C[-25,-25,start:end]

# Calculate bases ANEW

S <- moran_fn(C[,,1])
S <- S$vectors[,(S$values > -1/24 & !dplyr::near(S$values, 0))]
R <- ncol(S)
etaSigma <- as.matrix(
  nearPD(t(S) %*% (diag(rowSums(C[,,1])) - C[,,1]) %*% S)$mat
)
Phi <- array(NA, dim=c(R,R,U-1))
xiSigma <- array(NA, dim=c(R,R,U-1))
for (i in 2:U) {
  Proj <- diag(1,R) - t(S) %*% hat_fn(X[,,i]) %*% S
  E2 <- eigen(as.matrix(nearPD(Proj %*% diag(1,R) %*% Proj)$mat))
  Phi[,,(i-1)] <-  E2$vectors
  xiSigma[,,(i-1)] <- as.matrix(
    nearPD(
      etaSigma - Phi[,,(i-1)] %*% etaSigma %*% t(Phi[,,(i-1)])
    )$mat
  )
}
rm(list=c("C", "E2", "Proj"))
gc()

# BUGS code as written in NIMBLE 0.7.0.
library(nimble)
code <- nimbleCode({
  # Indexing: Space runs first, Time runs last
  # Space: 1~N; Variable/Spectral: 1~P/R; Time: 1~U.
  # 1. Likelihood
  for (i in 1:N) {
    for (t in 1:U) {
      mu[i,t] <- beta0[t] + sum(X[i,1:P,t] * beta[1:P,t]) + 
        sum(S[i,1:R] %*% eta[1:R,t])
      y[i,t] ~ dnorm(mu[i,t], var=var_eps[i,i])
      ## X: covariates (design matrix)
      ## S: basis function matrix
      ## var_eps: obs. error (time-dependent)
    }
  }
  
  # 2. Latent state/random effects
  zerosR[1:R] <- rep(0,R)
  var_eta_mat[1:R,1:R] <- diag(rep(var_eta, R))
  etaSigma1[1:R,1:R] <- var_eta_mat[1:R,1:R] %*% 
    etaSigma[1:R,1:R]
  eta[1:R,1] ~ dmnorm(zerosR[1:R], 
                      cov=etaSigma1[1:R,1:R])
  for (t in 1:(U-1)) {
    etaMuT[1:R,t] <- Phi[1:R,1:R,t] %*% eta[1:R,t]
    xiSigmaT[1:R,1:R,t] <- var_eta_mat[1:R,1:R] %*% 
      xiSigma[1:R,1:R,t]
    eta[1:R,(t+1)] ~ dmnorm(etaMuT[1:R,t], 
                            cov=xiSigmaT[1:R,1:R,t])
    ## Covariance matrix == var_eta * target matrix
    ## Phi: transition matrix
  }
  
  # 3. Priors
  # 3.1. Regression slope
  zerosP[1:P] <- rep(0, P)
  betaCov[1:P,1:P] <- diag(rep(10^12, P))
  for (t in 1:U) {
    beta0[t] ~ dnorm(0, var=10^12)
    beta[1:P,t] ~ dmnorm(zerosP[1:P], cov=betaCov[1:P,1:P])
  }
    
  # 3.2. Variance parameters
  var_eta ~ dinvgamma(2,1)
  for (i in 1:N) {
    var_eps[i,i] ~ dunif(0,10^12)
    for (j in (i+1):N) {
      var_eps[i,j] <- 0
      var_eps[j,i] <- 0
    }
  }
  
  # Deterministic quantity (for sampling var_eta)
  SSER[1] <- (t(eta[1:R,1]) %*% inverse(etaSigma[1:R,1:R]) %*% (eta[1:R,1]))[1,1]
  for (t in 2:U) {
    SSER[t] <- SSER[t-1] + (
      t(eta[1:R,t] - etaMuT[1:R,(t-1)]) %*%
        inverse(xiSigmaT[1:R,1:R,(t-1)]) %*% (eta[1:R,t] - etaMuT[1:R,(t-1)])
    )[1,1]
  }
})
constants <- list(N=N, U=U, R=R, P=P)
data <- list(y=y, X=X, S=S, Phi=Phi, etaSigma=etaSigma, xiSigma=xiSigma)
Rmodel <- nimbleModel(code, constants, data, check=FALSE, calculate=FALSE)

# Build top-level parameter samplers & filter
conf <- configureMCMC(Rmodel, nodes=c("beta0", "y"))
for (t in 1:U) {
  conf$addSampler(paste0("beta[1:15, ", t, "]"), "betaVec")
}
conf$addSampler("var_eta", "varEta")
diagNodes <- paste0("var_eps[", 1:N, ", ", 1:N, "]")
conf$addSampler(diagNodes, "AF_slice")
conf$addMonitors("y")
Rmcmc <- buildMCMC(conf, useConjugacy=FALSE)
Rfilter <- buildEnsembleKF(Rmodel, "eta")

# Compile
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project=Rmodel)
Cfilter <- compileNimble(Rfilter, project=Rmodel, resetFunctions=TRUE)

# Initializing: where to start to explore?
# slope <- OLS estimate as Mode
library(abind)
y_ols <- c(y)
X_ols <- matrix(NA, N, P)
for (i in 1:U) {
  if (i==1) X_ols[1:N,] <- X[,,i]
  else X_ols <- rbind(X_ols, X[,,i])
}
beta_ols <- as.numeric(coef(lm(y ~ ., data=cbind.data.frame(y=y_ols, X_ols))))
for (t in 1:U) { 
  Cmodel[[paste0("beta0[", t, "]")]] <- beta_ols[1]
  Cmodel[[paste0("beta[1:15, ", t, "]")]] <- beta_ols[2:16]
}
# Run chains -> store
# 3 chains, each 10,000 iters

# t1 <- Sys.time()
# samplesList <- list(chain1 = list(), chain2 = list(), chain3 = list())
samplesList <- list()
set.seed(0)
#for (i in 1:3) {
  Cmodel$var_eps <- diag(runif(N, 0.1, 1))
  Cmodel$var_eta <- rinvgamma(1, 2, 1)
  Cmodel$eta[,1] <- chol(Cmodel$var_eta * Cmodel$etaSigma) %*% rnorm(R) 
  for (t in 2:U) {  
    Cmodel$eta[,t] <- Cmodel$Phi[,,(t-1)] %*% Cmodel$eta[,(t-1)] +
      chol(Cmodel$var_eta * Cmodel$xiSigma[,,(t-1)]) %*% rnorm(R)
  }
  for (j in 1:11) {
    if (j==1) { # burn-in
      Cmcmc$run(2000)
      Cfilter$run(2000)
      Cmcmc$run(1000)
    }
    else {
      Cfilter$run(1000)
      Cmcmc$run(1000)
      #samplesList[[i]][[(j-1)]] <- as.matrix(Cmcmc$mvSamples)
      samplesList[[(j-1)]] <- as.matrix(Cmcmc$mvSamples)
    }
  }
#}
#t2 <- Sys.time()
#(t2 - t1) / 3 ## How much time has been spent??

# chain1 <- Reduce(rbind, samplesList[[1]])
# chain2 <- Reduce(rbind, samplesList[[2]])
# chain3 <- Reduce(rbind, samplesList[[3]])
# samples <- list(chain1=chain1, chain2=chain2, chain3=chain3)
# samples_total <- Reduce(rbind, samples)
# save(samples, file="data/1.RData")

# Node names
samples <- Reduce(rbind, samplesList)
colnames(samples) <- colnames(samplesList[[1]])
betaIndex <- c(paste0("beta0[", 1:U, "]"), 
               paste0(rep(paste0("beta[", 1:P, ", "), each=U), 1:U, "]"))
params <- c(betaIndex, "var_eta", diagNodes)
Cmodel$resetData()
simMat <- array(NA, dim=c(N, U, 100))
set.seed(0)
for (i in 1:100) {
  for (p in params) Cmodel[[p]] <- samples[sample(10000, 1), p]
  while(dplyr::near(0, Cmodel[["var_eta"]])) Cmodel[["var_eta"]] <- samples[sample(10000, 1), "var_eta"]
  Cmodel$calculate()
  Cmodel$simulate("eta")
  Cmodel$calculate()
  Cmodel$simulate("y")
  simMat[,,i] <- Cmodel$y
}
save(simMat, file="data/sim1.RData")

#######

rm(list=ls())
gc()

# Constants and samplers are sourced from a separate file
load("data/modelConstants.RData")
source("code/sampler_code.R")
source("code/000-utilities.R")
library(Matrix)

start <- dt_search_fn("2017-08-16 00:00:00")
end <- dt_search_fn("2017-08-20 23:00:00")
U <- end - start + 1
N <- 24
y <- y[-20,start:end]
X <- X[-20,,start:end] # intercept
C <- C[-20,-20,start:end]

# Calculate bases ANEW

S <- moran_fn(C[,,1])
S <- S$vectors[,(S$values > -1/24 & !dplyr::near(S$values, 0))]
R <- ncol(S)
etaSigma <- as.matrix(
  nearPD(t(S) %*% (diag(rowSums(C[,,1])) - C[,,1]) %*% S)$mat
)
Phi <- array(NA, dim=c(R,R,U-1))
xiSigma <- array(NA, dim=c(R,R,U-1))
for (i in 2:U) {
  Proj <- diag(1,R) - t(S) %*% hat_fn(X[,,i]) %*% S
  E2 <- eigen(as.matrix(nearPD(Proj %*% diag(1,R) %*% Proj)$mat))
  Phi[,,(i-1)] <-  E2$vectors
  xiSigma[,,(i-1)] <- as.matrix(
    nearPD(
      etaSigma - Phi[,,(i-1)] %*% etaSigma %*% t(Phi[,,(i-1)])
    )$mat
  )
}
rm(list=c("C", "E2", "Proj"))
gc()

# BUGS code as written in NIMBLE 0.7.0.
library(nimble)
code <- nimbleCode({
  # Indexing: Space runs first, Time runs last
  # Space: 1~N; Variable/Spectral: 1~P/R; Time: 1~U.
  # 1. Likelihood
  for (i in 1:N) {
    for (t in 1:U) {
      mu[i,t] <- beta0[t] + sum(X[i,1:P,t] * beta[1:P,t]) + 
        sum(S[i,1:R] %*% eta[1:R,t])
      y[i,t] ~ dnorm(mu[i,t], var=var_eps[i,i])
      ## X: covariates (design matrix)
      ## S: basis function matrix
      ## var_eps: obs. error (time-dependent)
    }
  }
  
  # 2. Latent state/random effects
  zerosR[1:R] <- rep(0,R)
  var_eta_mat[1:R,1:R] <- diag(rep(var_eta, R))
  etaSigma1[1:R,1:R] <- var_eta_mat[1:R,1:R] %*% 
    etaSigma[1:R,1:R]
  eta[1:R,1] ~ dmnorm(zerosR[1:R], 
                      cov=etaSigma1[1:R,1:R])
  for (t in 1:(U-1)) {
    etaMuT[1:R,t] <- Phi[1:R,1:R,t] %*% eta[1:R,t]
    xiSigmaT[1:R,1:R,t] <- var_eta_mat[1:R,1:R] %*% 
      xiSigma[1:R,1:R,t]
    eta[1:R,(t+1)] ~ dmnorm(etaMuT[1:R,t], 
                            cov=xiSigmaT[1:R,1:R,t])
    ## Covariance matrix == var_eta * target matrix
    ## Phi: transition matrix
  }
  
  # 3. Priors
  # 3.1. Regression slope
  zerosP[1:P] <- rep(0, P)
  betaCov[1:P,1:P] <- diag(rep(10^12, P))
  for (t in 1:U) {
    beta0[t] ~ dnorm(0, var=10^12)
    beta[1:P,t] ~ dmnorm(zerosP[1:P], cov=betaCov[1:P,1:P])
  }
  
  # 3.2. Variance parameters
  var_eta ~ dinvgamma(2,1)
  for (i in 1:N) {
    var_eps[i,i] ~ dunif(0,10^12)
    for (j in (i+1):N) {
      var_eps[i,j] <- 0
      var_eps[j,i] <- 0
    }
  }
  
  # Deterministic quantity (for sampling var_eta)
  SSER[1] <- (t(eta[1:R,1]) %*% inverse(etaSigma[1:R,1:R]) %*% (eta[1:R,1]))[1,1]
  for (t in 2:U) {
    SSER[t] <- SSER[t-1] + (
      t(eta[1:R,t] - etaMuT[1:R,(t-1)]) %*%
        inverse(xiSigmaT[1:R,1:R,(t-1)]) %*% (eta[1:R,t] - etaMuT[1:R,(t-1)])
    )[1,1]
  }
})
constants <- list(N=N, U=U, R=R, P=P)
data <- list(y=y, X=X, S=S, Phi=Phi, etaSigma=etaSigma, xiSigma=xiSigma)
Rmodel <- nimbleModel(code, constants, data, check=FALSE, calculate=FALSE)

# Build top-level parameter samplers & filter
conf <- configureMCMC(Rmodel, nodes=c("beta0", "y"))
for (t in 1:U) {
  conf$addSampler(paste0("beta[1:15, ", t, "]"), "betaVec")
}
conf$addSampler("var_eta", "varEta")
diagNodes <- paste0("var_eps[", 1:N, ", ", 1:N, "]")
conf$addSampler(diagNodes, "AF_slice")
conf$addMonitors("y")
Rmcmc <- buildMCMC(conf, useConjugacy=FALSE)
Rfilter <- buildEnsembleKF(Rmodel, "eta")

# Compile
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project=Rmodel)
Cfilter <- compileNimble(Rfilter, project=Rmodel, resetFunctions=TRUE)

# Initializing: where to start to explore?
# slope <- OLS estimate as Mode
library(abind)
y_ols <- c(y)
X_ols <- matrix(NA, N, P)
for (i in 1:U) {
  if (i==1) X_ols[1:N,] <- X[,,i]
  else X_ols <- rbind(X_ols, X[,,i])
}
beta_ols <- as.numeric(coef(lm(y ~ ., data=cbind.data.frame(y=y_ols, X_ols))))
for (t in 1:U) { 
  Cmodel[[paste0("beta0[", t, "]")]] <- beta_ols[1]
  Cmodel[[paste0("beta[1:15, ", t, "]")]] <- beta_ols[2:16]
}
# Run chains -> store
# 3 chains, each 10,000 iters

# t1 <- Sys.time()
# samplesList <- list(chain1 = list(), chain2 = list(), chain3 = list())
samplesList <- list()
set.seed(0)
#for (i in 1:3) {
Cmodel$var_eps <- diag(runif(N, 0.1, 1))
Cmodel$var_eta <- rinvgamma(1, 2, 1)
Cmodel$eta[,1] <- chol(Cmodel$var_eta * Cmodel$etaSigma) %*% rnorm(R) 
for (t in 2:U) {  
  Cmodel$eta[,t] <- Cmodel$Phi[,,(t-1)] %*% Cmodel$eta[,(t-1)] +
    chol(Cmodel$var_eta * Cmodel$xiSigma[,,(t-1)]) %*% rnorm(R)
}
for (j in 1:11) {
  if (j==1) { # burn-in
    Cmcmc$run(2000)
    Cfilter$run(2000)
    Cmcmc$run(1000)
  }
  else {
    Cfilter$run(1000)
    Cmcmc$run(1000)
    #samplesList[[i]][[(j-1)]] <- as.matrix(Cmcmc$mvSamples)
    samplesList[[(j-1)]] <- as.matrix(Cmcmc$mvSamples)
  }
}
#}
#t2 <- Sys.time()
#(t2 - t1) / 3 ## How much time has been spent??

# chain1 <- Reduce(rbind, samplesList[[1]])
# chain2 <- Reduce(rbind, samplesList[[2]])
# chain3 <- Reduce(rbind, samplesList[[3]])
# samples <- list(chain1=chain1, chain2=chain2, chain3=chain3)
# samples_total <- Reduce(rbind, samples)
# save(samples, file="data/1.RData")

# Node names
samples <- Reduce(rbind, samplesList)
colnames(samples) <- colnames(samplesList[[1]])
betaIndex <- c(paste0("beta0[", 1:U, "]"), 
               paste0(rep(paste0("beta[", 1:P, ", "), each=U), 1:U, "]"))
params <- c(betaIndex, "var_eta", diagNodes)
Cmodel$resetData()
simMat <- array(NA, dim=c(N, U, 100))
set.seed(0)
for (i in 1:100) {
  for (p in params) Cmodel[[p]] <- samples[sample(10000, 1), p]
  while(dplyr::near(0, Cmodel[["var_eta"]])) Cmodel[["var_eta"]] <- samples[sample(10000, 1), "var_eta"]
  Cmodel$calculate()
  Cmodel$simulate("eta")
  Cmodel$calculate()
  Cmodel$simulate("y")
  simMat[,,i] <- Cmodel$y
}
save(simMat, file="data/sim2.RData")

###########

rm(list=ls())
gc()
# Constants and samplers are sourced from a separate file
load("data/modelConstants.RData")
source("code/sampler_code.R")
source("code/000-utilities.R")
library(Matrix)

start <- dt_search_fn("2016-03-30 00:00:00")
end <- dt_search_fn("2016-04-03 23:00:00")
U <- end - start + 1
N <- 24
y <- y[-20,start:end]
X <- X[-20,,start:end] # intercept
C <- C[-20,-20,start:end]

# Calculate bases ANEW
R <- 4
S <- array(NA, dim=c(N,R,U))
Phi <- array(NA, dim=c(R,R,U-1))
etaSigma <- array(NA, dim=c(R,R,U))
xiSigma <- array(NA, dim=c(R,R,U-1))
for (i in 1:U) {
  E1 <- moran_fn(C[,,i])
  S[,,i] <- E1$vectors[,(E1$values > -1/24 & !dplyr::near(E1$values, 0))]
  etaSigma[,,i] <- as.matrix(
    nearPD(t(S[,,i]) %*% (diag(rowSums(C[,,i])) - C[,,i]) %*% S[,,i])$mat
  )
  if (i > 1) {
    Proj <- diag(1,R) - t(S[,,i]) %*% hat_fn(X[,,i]) %*% S[,,i]
    E2 <- eigen(as.matrix(nearPD(Proj %*% diag(1,R) %*% Proj)$mat))
    Phi[,,(i-1)] <-  E2$vectors
    xiSigma[,,(i-1)] <- as.matrix(
      nearPD(
        etaSigma[,,i] - Phi[,,(i-1)] %*% etaSigma[,,(i-1)] %*% t(Phi[,,(i-1)])
      )$mat
    )
  }
}
rm(list=c("C", "E1", "E2", "Proj"))
gc()

# BUGS code as written in NIMBLE 0.7.0.
library(nimble)
code <- nimbleCode({
  # Indexing: Space runs first, Time runs last
  # Space: 1~N; Variable/Spectral: 1~P/R; Time: 1~U.
  # 1. Likelihood
  for (i in 1:N) {
    for (t in 1:U) {
      mu[i,t] <- beta0[t] + sum(X[i,1:P,t] * beta[1:P,t]) + 
        sum(S[i,1:R,t] %*% eta[1:R,t])
      y[i,t] ~ dnorm(mu[i,t], var=var_eps[i,i])
      ## X: covariates (design matrix)
      ## S: basis function matrix
      ## var_eps: obs. error (time-dependent)
    }
  }
  
  # 2. Latent state/random effects
  zerosR[1:R] <- rep(0,R)
  var_eta_mat[1:R,1:R] <- diag(rep(var_eta, R))
  etaSigma1[1:R,1:R] <- var_eta_mat[1:R,1:R] %*% 
    etaSigma[1:R,1:R,1]
  eta[1:R,1] ~ dmnorm(zerosR[1:R], 
                      cov=etaSigma1[1:R,1:R])
  for (t in 1:(U-1)) {
    etaMuT[1:R,t] <- Phi[1:R,1:R,t] %*% eta[1:R,t]
    xiSigmaT[1:R,1:R,t] <- var_eta_mat[1:R,1:R] %*% 
      xiSigma[1:R,1:R,t]
    eta[1:R,(t+1)] ~ dmnorm(etaMuT[1:R,t], 
                            cov=xiSigmaT[1:R,1:R,t])
    ## Covariance matrix == var_eta * target matrix
    ## Phi: transition matrix
  }
  
  # 3. Priors
  # 3.1. Regression slope
  zerosP[1:P] <- rep(0, P)
  betaCov[1:P,1:P] <- diag(rep(10^12, P))
  for (t in 1:U) {
    beta0[t] ~ dnorm(0, var=10^12)
    beta[1:P,t] ~ dmnorm(zerosP[1:P], cov=betaCov[1:P,1:P])
  }
  
  # 3.2. Variance parameters
  var_eta ~ dinvgamma(2,1)
  for (i in 1:N) {
    var_eps[i,i] ~ dunif(0,10^12)
    for (j in (i+1):N) {
      var_eps[i,j] <- 0
      var_eps[j,i] <- 0
    }
  }
  
  # Deterministic quantity (for sampling var_eta)
  SSER[1] <- (t(eta[1:R,1]) %*% inverse(etaSigma[1:R,1:R,1]) %*% (eta[1:R,1]))[1,1]
  for (t in 2:U) {
    SSER[t] <- SSER[t-1] + (
      t(eta[1:R,t] - etaMuT[1:R,(t-1)]) %*%
        inverse(xiSigmaT[1:R,1:R,(t-1)]) %*% (eta[1:R,t] - etaMuT[1:R,(t-1)])
    )[1,1]
  }
})
constants <- list(N=N, U=U, R=R, P=P)
data <- list(y=y, X=X, S=S, Phi=Phi, etaSigma=etaSigma, xiSigma=xiSigma)
Rmodel <- nimbleModel(code, constants, data, check=FALSE, calculate=FALSE)

# Build top-level parameter samplers & filter
conf <- configureMCMC(Rmodel, nodes=c("beta0", "y"))
for (t in 1:U) {
  conf$addSampler(paste0("beta[1:15, ", t, "]"), "betaVec2")
}
conf$addSampler("var_eta", "varEta")
diagNodes <- paste0("var_eps[", 1:N, ", ", 1:N, "]")
conf$addSampler(diagNodes, "AF_slice")
conf$addMonitors("y")
Rmcmc <- buildMCMC(conf, useConjugacy=FALSE)
Rfilter <- buildEnsembleKF(Rmodel, "eta")

# Compile
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project=Rmodel)
Cfilter <- compileNimble(Rfilter, project=Rmodel, resetFunctions=TRUE)

# Initializing: where to start to explore?
# slope <- OLS estimate as Mode
library(abind)
y_ols <- c(y)
X_ols <- matrix(NA, N, P)
for (i in 1:U) {
  if (i==1) X_ols[1:N,] <- X[,,i]
  else X_ols <- rbind(X_ols, X[,,i])
}
beta_ols <- as.numeric(coef(lm(y ~ ., data=cbind.data.frame(y=y_ols, X_ols))))
for (t in 1:U) { 
  Cmodel[[paste0("beta0[", t, "]")]] <- beta_ols[1]
  Cmodel[[paste0("beta[1:15, ", t, "]")]] <- beta_ols[2:16]
}
# Run chains -> store
# 3 chains, each 10,000 iters

# t1 <- Sys.time()
# samplesList <- list(chain1 = list(), chain2 = list(), chain3 = list())
samplesList <- list()
set.seed(0)
#for (i in 1:3) {
  Cmodel$var_eps <- diag(runif(N, 0.1, 1))
  Cmodel$var_eta <- rinvgamma(1, 2, 1)
  Cmodel$eta[,1] <- chol(Cmodel$var_eta * Cmodel$etaSigma[,,1]) %*% rnorm(R) 
  for (t in 2:U) {  
    Cmodel$eta[,t] <- Cmodel$Phi[,,(t-1)] %*% Cmodel$eta[,(t-1)] +
      chol(Cmodel$var_eta * Cmodel$xiSigma[,,(t-1)]) %*% rnorm(R)
  }
  for (j in 1:11) {
    if (j==1) { # burn-in
      Cmcmc$run(2000, reset=TRUE)
      Cfilter$run(2000)
      Cmcmc$run(1000)
    }
    else {
      Cfilter$run(1000)
      Cmcmc$run(1000)
      # samplesList[[i]][[(j-1)]] <- as.matrix(Cmcmc$mvSamples)
      samplesList[[(j-1)]] <- as.matrix(Cmcmc$mvSamples)
    }
  }
#}
#t2 <- Sys.time()
#(t2 - t1) / 3 ## How much time has been spent??

# chain1 <- Reduce(rbind, samplesList[[1]])
# chain2 <- Reduce(rbind, samplesList[[2]])
# chain3 <- Reduce(rbind, samplesList[[3]])
# samples <- list(chain1=chain1, chain2=chain2, chain3=chain3)
# samples_total <- Reduce(rbind, samples)
# save(samples, file="data/3.RData")

# Node names
samples <- Reduce(rbind, samplesList)
colnames(samples) <- colnames(samplesList[[1]])
betaIndex <- c(paste0("beta0[", 1:U, "]"), 
               paste0(rep(paste0("beta[", 1:P, ", "), each=U), 1:U, "]"))
params <- c(betaIndex, "var_eta", diagNodes)
Cmodel$resetData()
simMat <- array(NA, dim=c(N, U, 100))
for (i in 1:100) {
  for (p in params) Cmodel[[p]] <- samples[sample(10000, 1), p]
  Cmodel$calculate()
  Cmodel$simulate("eta")
  Cmodel$calculate()
  Cmodel$simulate("y")
  simMat[,,i] <- Cmodel$y
}
save(simMat, file="data/sim3.RData")
