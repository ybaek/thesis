library(Matrix)
load("data/modelConstants.RData")
load("data/1.RData")
load("data/ReducedData.RData")
source("code/000-utilities.R")
start <- dt_search_fn("2015-05-16 00:00:00")
end <- dt_search_fn("2015-05-20 23:00:00")

# To test predictions:
future <- y[-25,(end+1):(end+12)]
future_X <- X[-25,,(end+1):(end+12)]

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

betaIndex <- c(paste0("beta0[", 1:U, "]"),
  paste0("beta[", rep(1:P, each=U), ", ", 1:U, "]"))
diagNodes <- paste0("var_eps[", 1:N, ", ", 1:N, "]")

params <- list()
varNames <- c(betaIndex, diagNodes, "var_eta")
samples_total <- Reduce(rbind, samples)
colnames(samples_total) <- colnames(samples$chain1)
set.seed(0)
for (v in varNames) {
  draw <- sample(10000*3,1)
  if (v=="var_eta") {
    while (dplyr::near(0, samples_total[draw,"var_eta"])) draw <- sample(10000*3, 1)
  }
  params[[v]] <- samples_total[draw,v]
}
library(mvtnorm)
Phi_k <- array(NA, dim=c(R, R, 13))
Phi_k[,,1] <- as.matrix(Phi[,,(U-1)])
xiSigma_k <- array(NA, dim=c(R, R, 12))
for (t in 1:12) {
  Proj <- diag(1,R) - t(S) %*% hat_fn(future_X[,,t]) %*% S
  E2 <- eigen(as.matrix(nearPD(Proj %*% diag(1,R) %*% Proj)$mat))
  Phi_k[,,(t+1)] <-  E2$vectors
  xiSigma_k[,,t] <- as.matrix(
    nearPD(etaSigma - Phi_k[,,t] %*% etaSigma %*% t(Phi_k[,,t]))$mat
  )
}  
Phi_k <- Phi_k[,,2:13]

# Kalman Forecasting
fore <- array(NA, dim=c(N, 12, 100))
set.seed(0)
for (i in 1:100) {
  latents <- matrix(NA, R, U)
  for (t in 1:U) {
    if (t==1) latents[,1] <- rmvnorm(1, rep(0, R), params[["var_eta"]] * etaSigma[,,U])
    else latents[,t] <- rmvnorm(1, 
                                Phi[,,(t-1)] %*% latents[,(t-1)],
                                params[["var_eta"]] * xiSigma[,,(t-1)]
    )
  }
  foreLatents <- matrix(NA, R, 13)
  foreLatents[,1] <- latents[,U]
  for (t in 1:12) {
    foreLatents[,(t+1)] <- rmvnorm(1, Phi[,,t] %*% foreLatents[,t],
                            params[["var_eta"]] * xiSigma[,,t]
    )
    for (n in 1:N) {
      mu <- params[["beta0"]] + t(future_X[n,1:P,t]) %*% unlist(params[betaIndex]) +
        t(S[n,1:R]) %*% foreLatents[,(t+1)]
      fore[n,t,i] <- rnorm(1, mu, sqrt(params[[paste0("var_eps[", n, ", ", n, "]")]]))
    }
  }
}
foreMat <- rowMeans(fore, dims=2)
rmses3 <- apply(fore, 3, function(X) mean((X-future)^2, na.rm=T))
rownames(future) <- unique(real_data$id)[-25]
colnames(future) <- seq(as.POSIXct("2015-05-21 00:00:00", tz="UTC"), 
                   as.POSIXct("2015-05-21 11:00:00", tz="UTC"), by=3600)
rownames(foreMat) <- unique(real_data$id)[-25]
colnames(foreMat) <- seq(as.POSIXct("2016-04-04 00:00:00", tz="UTC"), 
                         as.POSIXct("2016-04-04 11:00:00", tz="UTC"), by=3600)

g1 <- ggplot(melt(future), aes(Var2, as.factor(Var1))) + geom_tile(aes(fill=value)) +
  labs(x="Hour", y="ID", fill=NULL) + scale_x_continuous(breaks=NULL) + 
  scale_fill_distiller(palette="BuPu", direction=1)
## Interpolations
future_interp <- future
future_interp[which(is.na(future))] <- foreMat[which(is.na(future))]
g2 <- ggplot(melt(future_interp), aes(Var2, as.factor(Var1))) + geom_tile(aes(fill=value)) +
  labs(x="Hour", y="ID", fill=NULL) + scale_x_continuous(breaks=NULL) + 
  scale_fill_distiller(palette="BuPu", direction=1, limits=c(1.5, 4.5), oob=scales::squish)
## Residuals
g3 <- ggplot(melt(abs(future-foreMat)), aes(Var2, as.factor(Var1))) + geom_tile(aes(fill=value)) +
  labs(x="Hour", y="ID", fill=NULL) + scale_x_continuous(breaks=NULL) + 
  scale_fill_distiller(palette="BuPu", direction=1)
save(list=c("g1", "g2", "g3"), file="data/plots6.RData")
