#
# Convergence plots
#

load("data/1.RData")
samples_t <- Reduce(rbind, samples)
colnames(samples_t) <- colnames(samples$chain1)

library(basicMCMCplots)
diagNodes <- paste0("var_eps[", 1:24, ", ", 1:24, "]")
## Trace plot
png("images/trace.png", 7, 6, "in", res=300)
par(mfrow=c(2,2), mar=c(3,3,2,1), mgp=c(1.8,1,0))
plot(samples$chain1[-c(1:1000),diagNodes[1]], type="l", col="red",
     xlab="Iteration", ylab=expression(sigma[paste(epsilon,1)]^2))
plot(samples$chain2[-c(1:1000),diagNodes[1]], type="l", col="red",
     xlab="Iteration", ylab=expression(sigma[paste(epsilon,1)]^2))
plot(samples$chain3[-c(1:1000),diagNodes[1]], type="l", col="red",
     xlab="Iteration", ylab=expression(sigma[paste(epsilon,1)]^2))
dev.off()

## Rank plot
ranks1 <- rank(Reduce(c, lapply(samples, function(x) x[-c(1:1000),diagNodes[1]])))

png("images/ranks.png", 7, 6, "in", res=300)
par(mfrow=c(2,2), mar=c(3,3,1,1), mgp=c(1.8,1,0))
hist(ranks1[1:9000], col="lightblue", main=NULL, 
     xlab=expression(paste("Ranks of ", sigma[paste(epsilon,1)]^2)))
hist(ranks1[9001:18000], col="lightblue", main=NULL, 
     xlab=expression(paste("Ranks of ", sigma[paste(epsilon,1)]^2)))
hist(ranks1[18001:27000], col="lightblue", main=NULL, 
     xlab=expression(paste("Ranks of ", sigma[paste(epsilon,1)]^2)))
dev.off()

# Using samplesPlot pkg
samplesPlot(samples$chain1, burnin=1000,
            diagNodes, legend=F, traceplot=F, width=7, height=6,
            file="images/densities")

chainsPlot(lapply(samples, function(x) x[-c(1:1000),]), 
           diagNodes, scale=T, cex=0.8, 
           width=7, height=6, file="images/chains")

library(coda)
samples_t <- as.mcmc.list(lapply(samples, function(x) as.mcmc(x[,c(diagNodes, "var_eta")])))
gelman.diag(samples_t)

######
library(ggplot2)
library(reshape2)
load("data/sim3.RData")
sim <- simMat[,,sample(100,1)]
rownames(sim) <- rownames(y)
colnames(sim) <- colnames(y)
g1 <- ggplot(melt(y), aes(Var2, value)) + geom_line(aes(col=rank(Var1)), show.legend=F) +
  labs(x="Hour", y="log") + 
  scale_x_continuous(breaks=NULL) +
  scale_color_distiller(palette="BuPu", direction=1)
g2 <- ggplot(melt(y-sim), aes(Var2, value)) + geom_line(aes(col=Var1), show.legend=F) +
  labs(x="Hour", y="log") + 
  scale_x_continuous(breaks=NULL) +
  scale_y_continuous(limits=c(-2,2)) +
  scale_color_distiller(palette="BuPu", direction=1)
ggsave("time3-1.png", g1, "png", "images/", width=6, height=3, units="in", dpi=320)
ggsave("time2-2.png", g2, "png", "images/", width=6, height=3, units="in", dpi=320)

mspe <- list()
for (i in 1:100) mspe[[i]] <- (y - simMat[,,i])^2
mspe <- sqrt(Reduce(`+`, mspe)/100)
rownames(mspe) <- rownames(y)
colnames(mspe) <- colnames(y)
g3 <- ggplot(melt(mspe), aes(Var2, value)) + geom_line(aes(col=rank(Var1)), show.legend=F) +
  labs(x="Hour", y="log") + 
  scale_x_continuous(breaks=NULL) +
  scale_color_distiller(palette="BuPu", direction=1)
ggsave("time3-3.png", g3, "png", "images/", width=6, height=3, units="in", dpi=320)
#####

load("data/sim3.RData")
png("images/check3.png", width=7, height=5, units="in", res=320)
par(mfrow=c(2,2), mar=c(4,1.96,1,1), mgp=c(1.9,1,0))
hist(apply(simMat, 3, function(X) max(X, na.rm=T)),
     xlim=c(4.9, 7), col="lightblue", main=NULL, xlab="Max")
abline(v=max(y, na.rm=T), col="red", lwd=2)
hist(apply(simMat, 3, function(X) min(X, na.rm=T)),
      col="lightblue", main=NULL, xlab="Min")
abline(v=min(y, na.rm=T), col="red", lwd=2)
hist(apply(simMat, 3, function(X) median(X, na.rm=T)),
     xlim=c(3.85, 3.92), col="lightblue", main=NULL, xlab="Median")
abline(v=median(y, na.rm=T), col="red", lwd=2)
hist(apply(simMat, 3, function(X) sd(X, na.rm=T)),
     col="lightblue", main=NULL, xlab="SD")
abline(v=sd(y, na.rm=T), col="red", lwd=2)
dev.off()