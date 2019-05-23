#
# Customized Gibbs samplers
# from full conditionals
#

# For compilation, I CAN'T use arrays
# -> use control argument to pass on 3rd dimension

library(nimble)

# 1. Regression slope full conditional
# [ beta[1:P,t] | eta[1:R,t], var_eps[t] ]
betaVec <- nimbleFunction(
  name = "betaVec",
  contains = sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    t <- as.integer(unlist(strsplit(strsplit(target, split=",")[[1]][2], "]")))
    # Access node names
    dataNode <- paste0("y[,", t, "]")
    XNode <- paste0("X[,,", t, "]")
    basisNode <- "S"
    latentNode <- paste0("eta[,", t, "]")
    oeNode <- paste0("var_eps")
    # Access dimensions (unchanged over sampling)
    N <- length(model[[dataNode]])
    P <- dim(model[[XNode]])[2]
    R <- dim(model[[basisNode]])[2]
    
    calcNodes <- model$getDependencies(target)
  },
  run = function() {
    # Access node values
    dataVec <- model[[dataNode]]
    XMat <- model[[XNode]]
    basisMat <- model[[basisNode]]
    latentVec <- model[[latentNode]]
    oeMat <- model[[oeNode]]
    # Acess node parameters: slope covariance matrix.
    priorBetaCov <- model$getParam(target, "cov")
    # Exact full conditional
    postVar <- inverse(t(XMat) %*% inverse(oeMat) %*% XMat + priorBetaCov)
    postMean <- (
      postVar %*% t(XMat) %*% inverse(oeMat) %*% (dataVec - (basisMat %*% latentVec)[,1])
    )[,1] # coerced into vector
    beta_draw <- rmnorm_chol(1, postMean, chol(postVar), prec_param=FALSE)
    # Update
    model[[target]] <<- beta_draw
    nimCopy(from=model, to=mvSaved, row=1, nodes=calcNodes, logProb=TRUE)
  },
  methods = list( reset = function(){} )
)

betaVec2 <- nimbleFunction(
  name = "betaVec2",
  contains = sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    t <- as.integer(unlist(strsplit(strsplit(target, split=",")[[1]][2], "]")))
    # Access node names
    dataNode <- paste0("y[,", t, "]")
    XNode <- paste0("X[,,", t, "]")
    basisNode <- paste0("S[,,", t, "]")
    latentNode <- paste0("eta[,", t, "]")
    oeNode <- paste0("var_eps")
    # Access dimensions (unchanged over sampling)
    N <- length(model[[dataNode]])
    P <- dim(model[[XNode]])[2]
    R <- dim(model[[basisNode]])[2]
    
    calcNodes <- model$getDependencies(target)
  },
  run = function() {
    # Access node values
    dataVec <- model[[dataNode]]
    XMat <- model[[XNode]]
    basisMat <- model[[basisNode]]
    latentVec <- model[[latentNode]]
    oeMat <- model[[oeNode]]
    # Acess node parameters: slope covariance matrix.
    priorBetaCov <- model$getParam(target, "cov")
    # Exact full conditional
    postVar <- inverse(t(XMat) %*% inverse(oeMat) %*% XMat + priorBetaCov)
    postMean <- (
      postVar %*% t(XMat) %*% inverse(oeMat) %*% (dataVec - (basisMat %*% latentVec)[,1])
    )[,1] # coerced into vector
    beta_draw <- rmnorm_chol(1, postMean, chol(postVar), prec_param=FALSE)
    # Update
    model[[target]] <<- beta_draw
    nimCopy(from=model, to=mvSaved, row=1, nodes=calcNodes, logProb=TRUE)
  },
  methods = list( reset = function(){} )
)

# 2. random effects scale factor (time-ind. variance)
# [ var_eta | eta[1:R,1:U] ]

varEta <- nimbleFunction(
  name = "varEta",
  contains = sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    latentNode <- "eta"
    R <- dim(model[[latentNode]])[1]
    U <- dim(model[[latentNode]])[2]
    sumNode <- paste0("SSER[", U, "]") # sum involving ALL latent effect nodes

    calcNodes <- model$getDependencies(target)
  },
  run = function() {
    # Access node parameters: inverse-gamma parameters for var_eta 
    priorAlpha <- model$getParam(target, "shape")
    priorBeta <- model$getParam(target, "rate")
    
    # Access node values (loopSum) / calculate conditional
    sumValue <- model[[sumNode]]
    postAlpha <- priorAlpha + (R*U)/2
    postBeta <- priorBeta + sumValue / 2
    
    # Update
    model[[target]] <<- rinvgamma(1, shape=postAlpha, rate=postBeta)
    nimCopy(from=model, to=mvSaved, row=1, nodes=calcNodes, logProb=TRUE)
  },
  methods = list( reset = function(){} )
)
