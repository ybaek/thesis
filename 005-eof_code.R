library(dplyr)
library(ggplot2)
library(imputeTS)
library(reshape2)
source("code/000-utilities.R")
load("data/modelConstants.RData")
y_impute <- apply(y, 1, function(ts) na.interpolation(ts, option="stine"))
E <- eof_fn(y_impute)$v[,1:9] # bases
R <- ncol(E)


# sim_list <- list()
# set.seed(0)
# for (i in 1:100) {
#   sim_data <- apply(y_impute, 1, function(x) x[sample(26, 26, replace=F)])
#   D <- eof_fn(t(sim_data))$d
#   sim_list[[i]] <-  D / sum(D)
# }
# simEOFs_long <- melt(Reduce(rbind, sim_list))[,2:3]
# png("figures/model/fig2-eof.png", width=7, height=5, units="in", res=250)
# par(mar=c(4.5,4,2,2)+.1)
# plot(E$d/sum(E$d), pch=19, col="red", type="b", 
#      xlab="Eigenvalue", ylab="Relative variance")
# boxplot(simEOFs_long$value ~ simEOFs_long$Var2, axes=F, add=T)
# dev.off()

######

# Let's fit a model (different, simpler)
## Which time window to go to?
dt_search_fn <- function(dt) {
  ## convert POSIX
  dt_n <- as.POSIXct(dt, tz="UTC")
  ## Search here
  start <- as.POSIXct("2015-01-01 00:00:00", tz="UTC")
  end <- as.POSIXct("2017-12-31 23:00:00", tz="UTC")
  time_seq <- seq(start, end, by=3600)
  index <- which(time_seq == dt_n)
  return(index)
}
start <- dt_search_fn("2016-03-30 00:00:00")
end <- dt_search_fn("2016-04-03 23:00:00")
U <- end - start + 1
N <- 24
y <- y[-20,start:end]
X <- X[-20,,start:end]

######

library(nimble)
code <- nimbleCode({
  for (i in 1:N) {
    for (t in 1:U) {
      mu[i,t] <- beta0 + sum(X[i,1:P,t] * beta[1:P]) + 
        sum(S[i,1:R] %*% eta[1:R,t])
      y[i,t] ~ dnorm(mu[i,t], var=var_eps[i,i])
    }
  }
  
  zerosR[1:R] <- rep(0,R)
  var_error_mat[1:R,1:R] <- diag(rep(var_error, R))
  eta[1:R,1] ~ dmnorm(zerosR[1:R], cov=var_error_mat[1:R,1:R])
  for (t in 1:(U-1)) {
    ## "Off-diagonal" model
    etaMuT[1:R,t] <- Phi[1:R,1:R] %*% eta[1:R,t]
    eta[1:R,(t+1)] ~ dmnorm(etaMuT[1:R,t], cov=var_error_mat[1:R,1:R])
  }
  
  beta0 ~ dnorm(0, var=10^12)
  for (p in 1:P) {
    beta[p] ~ dnorm(0, var=10^12)
  }
  
  var_error ~ dinvgamma(2, 1)
  for (i in 1:N) {
    var_eps[i,i] ~ dunif(0,10^12)
    for (j in (i+1):N) {
      var_eps[i,j] <- 0
      var_eps[j,i] <- 0
    }
  }
  
  # Another parameters in this case: Phi's
  ## Off-diagonals
  for (r in 1:(R-1)) {
    Phi[r,(r+1)] ~ dunif(-0.999, 0.999)
    Phi[r,r] ~ dunif(-0.999, 0.999)
    Phi[(r+1),r] ~ dunif(-0.999, 0.999)
  }
  Phi[R,R] ~ dunif(-0.999, 0.999)
})
constants = list(P=P, R=R, U=U, N=N)
data = list(y=y, X=X, S=E)
Rmodel <- nimbleModel(code, constants, data, check=F, calculate=F)
conf <- configureMCMC(Rmodel, nodes=c("var_error", "y"))
conf$addSampler(c("beta0", "beta[1:15]"), "AF_slice")
diagNodes <- paste0("var_eps[", 1:N, ", ", 1:N, "]")
phiDiagNodes <- paste0("Phi[", 1:R, ", ", 1:R, "]")
offDiagNodes <- paste0("Phi[", c(1,2,2,3,3,4,4,5,5,6,6,7,7,8,8,9), ", ",
                       c(2,1,3,2,4,3,5,4,6,5,7,6,8,7,9,8), "]")
conf$addSampler(diagNodes, "AF_slice")
conf$addSampler(c(phiDiagNodes, offDiagNodes), "AF_slice")
Rmcmc <- buildMCMC(conf, useConjugacy=FALSE)
Rfilter <- buildEnsembleKF(Rmodel, "eta")

# Compile
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project=Rmodel)
Cfilter <- compileNimble(Rfilter, project=Rmodel, resetFunctions=TRUE)

library(abind)
y_ols <- c(y)
X_ols <- matrix(NA, N, P)
for (i in 1:U) {
  if (i==1) X_ols[1:N,] <- X[,,i]
  else X_ols <- rbind(X_ols, X[,,i])
}
beta_ols <- as.numeric(coef(lm(y ~ ., data=cbind.data.frame(y=y_ols, X_ols))))
Cmodel$beta0 <- beta_ols[1]
Cmodel$beta <- beta_ols[2:16]
set.seed(0)
Cmodel$var_error <- rinvgamma(1, 2, 1)
Cmodel$var_eps <- diag(runif(N, 0.1, 1))
Cmodel$Phi <- diag(0, R)
for (r in 1:(R-1)) {
  Cmodel[[paste0("Phi[", r+1, ", ", r, "]")]] <- runif(1, -0.999, 0.999)
  Cmodel[[paste0("Phi[", r, ", ", r+1, "]")]] <- runif(1, -0.999, 0.999)
  Cmodel[[paste0("Phi[", r, ", ", r, "]")]] <- runif(1, -0.999, 0.999)
}
Cmodel[["Phi[9, 9]"]] <- runif(1, -0.999, 0.999)
Cmodel$eta[,1] <- rnorm(R) 
for (t in 2:U) {  
  Cmodel$eta[,t] <- Cmodel$Phi %*% Cmodel$eta[,(t-1)] +
    rmnorm_chol(1, rep(0, R), chol(diag(Cmodel$var_error, R)))
}
samples <- list()

t1 <- Sys.time() # Start!
set.seed(0)
for (j in 1:11) {
  if (j==1) { # burn-in
    Cmcmc$run(2000, reset=TRUE)
    Cfilter$run(1000)
    Cmcmc$run(1000)
  }
  else {
    Cfilter$run(1000)
    Cmcmc$run(1000)
    samples[[(j-1)]] <- as.matrix(Cmcmc$mvSamples)
  }
}
t2 <- Sys.time()
t2 - t1 ## How much time has been spent??
samples <- Reduce(rbind, samples)

sim <- array(NA, dim=c(N, U, 100))
Cmodel$resetData()
params <- c("beta0", paste0("beta[", 1:P, "]"), diagNodes, phiDiagNodes, offDiagNodes, "var_error")
set.seed(0)
draw <- sample(10000,1)
for (v in params) Cmodel[[v]] <- samples[draw,v]
for (i in 1:100) {
  Cmodel$calculate()
  Cmodel$simulate("eta")
  Cmodel$calculate()
  Cmodel$simulate("y")
  sim[,,i] <- Cmodel$y
}

rmses3 <- apply(sim, 3, function(X) mean((X-y)^2, na.rm=T))
save(rmses3, file="data/rmses3.RData")
