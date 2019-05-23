
# Constants and samplers are sourced from a separate file
load("data/modelConstants.RData")
source("sampler_code.R")
library(methods)
library(utils)
.Machine$double.eps <- .Machine$double.eps^2

# 2015 Oct. 1st ~ 2015 Oct. 15th (NA's approx. 10%)
start <- 6553
end <- 7296
U <- end - start + 1
y <- y[,start:end]
X <- X[,,start:end]
S <- S[,,start:end]
etaSigma <- etaSigma[,,start:end]
Phi <- Phi[,,start:(end-1)]
xiSigma <- xiSigma[,,start:(end-1)]

rm(C)
gc()

var(c(y), na.rm=T) ## 0.35
set.seed(0)
sim_error <- matrix(rnorm(N*U, sd=sqrt(0.7)), N, U) # high StoN ratio
y_sim <- y + sim_error

# Configuring the model
library(nimble)
code <- nimbleCode({
  # Indexing: Space runs first, Time runs last
  # Space: 1~N; Variable/Spectral: 1~P/R; Time: 1~U.
  # 1. Likelihood
  for (i in 1:N) {
    for (t in 1:U) {
      mu[i,t] <- sum(X[i,1:P,t] * beta[1:P,t]) + sum(S[i,1:R,t] * eta[1:R,t])
      # X: covariates (design matrix)
      # S: basis function matrix
      # var_eps: obs. error (time-dependent)
      y[i,t] ~ dnorm(mu[i,t], var=var_eps[t] + sim_var)
    }
  }
  
  # 2. Latent state/random effects
  zerosR[1:R] <- rep(0,R)
  etaSigma1[1:R,1:R] <- var_eta * etaSigma[1:R,1:R,1]
  # cov. matrix == var_eta * target cov. matrix (constant)
  eta[1:R,1] ~ dmnorm(mean=zerosR[1:R], cov=etaSigma1[1:R,1:R]) # zero mean
  for (t in 1:(U-1)) {
    etaMuT[1:R,t] <- Phi[1:R,1:R,t] %*% eta[1:R,t]
    # Phi: transition matrix
    etaSigmaT[1:R,1:R,t] <- var_eta * xiSigma[1:R,1:R,t]
    eta[1:R,(t+1)] ~ dmnorm(mean=etaMuT[1:R,t], cov=etaSigmaT[1:R,1:R,t])
  }
  
  # 3. Priors
  # 3.1. Regression slope
    zerosP[1:P] <- rep(0, P)
    betaCov[1:P,1:P] <- diag(rep(10^12,P))
    for (t in 1:U) {
      beta[1:P,t] ~ dmnorm(zerosP[1:P], cov=betaCov[1:P,1:P]) # Normal iid
    }
  # 3.2. Variance parameters
  var_eta ~ dinvgamma(2,1) # Inverse gamma
  for (t in 1:U) {
    var_eps[t] ~ dinvgamma(2,1) # Inverse gamma iid
  }
  sim_var ~ dinvgamma(2,1)
  # Deterministic quantity (for sampling var_eta)  
  loopSum[1] <- (t(eta[1:R,1]) %*% inverse(etaSigma1[1:R,1:R]) %*% (eta[1:R,1]))[1,1]
  for (t in 2:U) {
    loopSum[t] <- loopSum[t-1] + (
      t(eta[1:R,t] - Phi[1:R,1:R,(t-1)] %*% eta[1:R,(t-1)]) %*% 
      inverse(etaSigmaT[1:R,1:R,(t-1)]) %*% (eta[1:R,t] - Phi[1:R,1:R,(t-1)] %*% eta[1:R,(t-1)])
    )[1,1]
  }
})

constants <- list(N=N, U=U, R=R, P=P, U=U)
data <- list(y=y_sim, X=X, S=S, Phi=Phi, etaSigma=etaSigma, xiSigma=xiSigma)
Rmodel <- nimbleModel(code, constants, data, check=FALSE, calculate=FALSE)

# Build top-level parameter samplers
conf <- configureMCMC(Rmodel, nodes=c("var_eps", "y", "sim_var"))
conf$addMonitors("y") 
for (t in 1:U) {
  betaNode <- paste0("beta[1:", P, ", ", t, "]")
  conf$addSampler(betaNode, "betaVec")
}
conf$addSampler("var_eta", "varEta")
Rmcmc <- buildMCMC(conf)

# Build latent state filter
Rfilter <- buildEnsembleKF(Rmodel, nodes="eta")

# Compile
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project=Rmodel)
Cfilter <- compileNimble(Rfilter, project=Rmodel, resetFunctions=TRUE)

# Run chains -> store
# 3 chains, each 10,000 iters

samplesList <- list(chain1 = list(), chain2 = list(), chain3 = list())

for (i in 1:3) {
  set.seed(0)
  Cmodel$beta <- matrix(runif(P*U, -.5, .5), nrow=P)
  Cmodel$var_eps <- runif(U, 0.01, 2)
  Cmodel$var_eta <- runif(1, 0.01, 2)
  Cmodel$sim_var <- runif(1, 0.01, 2)
  Cmodel$eta <- matrix(rep(1, R*U), nrow=R)
  for (j in 1:10) {
    if (j==1) {
      Cmcmc$run(1000, reset=TRUE)
      Cfilter$run(1000)
    }
    else {
      Cmcmc$run(1000)
      Cfilter$run(1000)
      samplesList[[i]][[(j-1)]] <- as.matrix(Cmcmc$mvSamples)
    }
  }
}

###############

save(list=c("chain1", "chain2", "chain3"), file="data/chains_sim.RData")

chain1 <- Reduce(rbind, samplesList$chain1)
chain2 <- Reduce(rbind, samplesList$chain2)
chain3 <- Reduce(rbind, samplesList$chain3)
rm(samplesList)
gc()

samples <- rbind(chain1[2000:9000,], chain2[2000:9000,], chain3[2000:9000,])
samples <- samples[,-5401]


#########

# Plots: to test differents similarities
library(dplyr)
library(ggplot2)
load("data/modelConstants.RData")
# source("003.01-distance.R")
real_data_for_plot <- real_data %>% arrange(dt) %>% head(25) %>% select(id, lon, lat)
coords_mat <- real_data_for_plot %>% select(lon,lat) %>% as.matrix()
dist_mat <- distMat_fn(coords_mat)
rownames(dist_mat) <- real_data_for_plot$id
colnames(dist_mat) <- real_data_for_plot$id
dist_mat1 <- normDist_fn(dist_mat)
sim_mat1 <- 1-dist_mat1
diag(sim_mat1) <- 0
sim_mat1 <- ifelse(sim_mat1 > .5, sim_mat1, 0)
rownames(sim_mat1) <- rownames(corMat)
colnames(sim_mat1) <- colnames(corMat)
sim_mat1_for_plot <- sim_mat1
sim_mat1_for_plot[upper.tri(sim_mat1)] <- NA
sim_mat1_for_plot <- melt(sim_mat1_for_plot) %>% mutate(Var1=as.factor(Var1),
                                               Var2=as.factor(Var2), value=value)
g3 <- ggplot(sim_mat1_for_plot, aes(Var1, Var2)) + geom_tile(aes(fill=value)) +
  scale_fill_distiller(palette="Spectral", direction=-1, lim=c(0,1), na.value=NA) +
  labs(x="ID", y="ID", fill="Similarity") +
  theme(axis.text.x = element_text(angle=45, vjust=.7, colour="black"), 
        axis.text.y = element_text(colour="black"),
        panel.background = element_rect(fill="white"),
        legend.margin=margin(0,0,0,0,"pt"))
ggsave(g3, filename="fig3-similarity_heat1.png", device="png", 
       path="figures/model/", width=5, height=5, dpi=320)

