#
# Some numerical functions
#

# Reading CSV and merging them for this dataset

## Handling unequal # of columns/seps (big data)
unequal_fn <- function(file) read.csv(file, header=T, fill=T)

mergecsv_fn <- function(path, key, order, sep) {
  # path refers to a path of csv files
  # key == string of year to filter from
  # order == order of the key in the filename
  # sep delimiter to split file name strings
  flist <- list.files(path, full.names=T)
  flist_sep <- strsplit(flist, split=sep)
  l <- list()
  for (i in 1:length(flist)) {
    f <- flist_sep[[i]]
    if (f[order]==key) l <- append(l, flist[i])
  }
  return(do.call(rbind.data.frame, lapply(l, unequal_fn)))
}

# Normalizing Euclidean pairwise distance
normDist_fn <- function(x, ...) (x-min(x, ...)) / (max(x, ...) - min(x, ...))

# Projection
hat_fn <- function(X, tol=.Machine$double.eps^2) return(X%*%solve(t(X)%*%X, tol=tol)%*%t(X))

# Orthogonal complement to C(X)
orthog_proj_fn <- function(X) {
  n <- dim(X)[1]
  I <- diag(1, n)
  return(I-hat_fn(X))
}

# Calculating Moran's I (no restriction on spatial weights)
moran_fn <- function(W, X=NULL) {
  n <- nrow(W)
  ones <- rep(1, n) 
  if (is.null(X)) {
    P <- diag(1, n) - (ones %o% ones / n)
  }
  else P <- orthog_proj_fn(X)
  return(eigen(P%*%W%*%P))
}

# Calculating empirical orthogonal functions (EOFs)
eof_fn <- function(Z) {
  ## The data matrix is indexed as: T (time) x n (space)
  temp <- dim(Z)[1]
  spat <- dim(Z)[2]
  Zscaled <- 1 / sqrt(temp-1) * (Z - outer(rep(1, temp), apply(Z, 2, mean)))
  E <- svd(Zscaled)
  return(E)
}

# Geodesic distance for coordinates data frame
distMat_fn <- function(m) {
  ## m is an array-like object ordered by rows into Long/lat
  result <- matrix(NA, nrow(m), nrow(m))
  for (i in 1:(nrow(m)-1)) {
    result[i,i] <- 0
    for (j in (i+1):nrow(m)) {
      d <- geosphere::distGeo(m[i,1:2], m[j,1:2])
      result[i,j] <- result[j,i] <- d
    }
    result[j,j] <- 0
  }
  return(result)
}

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

