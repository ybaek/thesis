#
# Interpolating missing values in covariates
# (Remote)
#

library(dplyr)
library(imputeTS)
library(purrr)
library(reshape2)

# Run the following for interpolation
load("data/FinalData.RData")

# Handling pressure and precipitation
pres_fill <- list()
for (n in unique(data_grid$id)) {
  sub_pres <- filter(data_grid, id==n)$pres
  sub_flag <- filter(data_grid, id==n)$flag
  unrecorded_inds <- which(is.na(sub_flag))
  for (i in 2:length(sub_pres)) {
    if (!(i %in% unrecorded_inds) & is.na(sub_pres[i])) {
      sub_pres[i] <- sub_pres[i-1]
    }
  }
  pres_fill[[n]] <- sub_pres
}
pres_fill <- unlist(pres_fill)
data_grid$pres <- pres_fill
precip_fill <- data_grid$precip
precip_fill[!is.na(data_grid$flag) & is.na(precip_fill)] <- 0.05
data_grid$precip <- precip_fill

# Cast into arrays
my_acast <- function(var) acast(data_grid, id ~ dt, value.var=var)
Xs <- lapply(list("temp", "humid", "ws", "pres", "precip"), my_acast)

# Cross-validation
# 1. Sampling indices
# set.seed(0)
# inds_list <- lapply(Xs, function(a) apply(a, 1, function(r) which(!is.na(r))))
# inds_list <- lapply(inds_list, function(l) lapply(l, function(x) sample(x, 1050)))
# inds_matlist <- lapply(inds_list, function(l) matrix(unlist(l), nrow=length(l), byrow=T))
# 
# # 2. Training set
# Xs_train <- Xs
# for (i in 1:length(Xs)) {
#   for (j in 1:nrow(Xs[[i]])) {
#     Xs_train[[i]][ j, inds_matlist[[i]][j,] ] <- NA
#   }
# }
# 
# # 3. Interpolation: using seasonally split splines
# my_interp_spl <- function(x) imputeTS::na.interpolation(x, option="spline")
# my_interp_stine <- function(x) imputeTS::na.interpolation(x, option="stine")
# my_interp_ma <- function(x) imputeTS::na.ma(x, k=6, weighting="linear")
# Xs_interp1 <- lapply(Xs_train, function(a) t(apply(a, 1, my_interp_spl)))
# Xs_interp2 <- lapply(Xs_train, function(a) t(apply(a, 1, my_interp_stine)))
# Xs_interp3 <- lapply(Xs_train, function(a) t(apply(a, 1, my_interp_ma)))
# Xs_random <- lapply(Xs_train, function(a) t(apply(a, 1, imputeTS::na.random)))
# 
# # Imputation results
# # RMSE
# RMSE <- matrix(0, 5, 4)
# rownames(RMSE) <- c("temp", "humid", "ws", "pres", "precip")
# colnames(RMSE) <- c("Spline", "Stineman", "Ma", "Random")
# for (i in 1:length(Xs)) {
#   for (j in 1:nrow(Xs[[i]])) {
#     RMSE[i, 1] <- RMSE[i,1] + sum((Xs_interp1[[i]][j, inds_matlist[[i]][j,]] - Xs[[i]][j, inds_matlist[[i]][j,]])^2)
#     RMSE[i, 2] <- RMSE[i,2] + sum((Xs_interp2[[i]][j, inds_matlist[[i]][j,]] - Xs[[i]][j, inds_matlist[[i]][j,]])^2)
#     RMSE[i, 3] <- RMSE[i,3] + sum((Xs_interp3[[i]][j, inds_matlist[[i]][j,]] - Xs[[i]][j, inds_matlist[[i]][j,]])^2)
#     RMSE[i, 4] <- RMSE[i,4] + sum((Xs_random[[i]][j, inds_matlist[[i]][j,]] - Xs[[i]][j, inds_matlist[[i]][j,]])^2)
#   }
# }
# RMSE <- RMSE / 26
# RMSE
# 
my_interp_stine <- function(x) na.seadec(ts(x, frequency=24), algorithm="interpolation",
                                         option="stine")

# Impute using whatever best algorithm
Xs <- lapply(Xs, function(a) t(apply(a, 1, my_interp_stine)))

# Then, back to original data structure:
for (i in 1:length(Xs)) {
  colnames(Xs[[i]]) <- seq(as.POSIXct("2006-01-01 00:00:00", tz="UTC"), 
    as.POSIXct("2017-12-31 23:00:00", tz="UTC"), by=3600)
}
Xs <- pmap(list(Xs, list("temp", "humid", "ws", "pres", "precip")),
     function(a, name) melt(a, varnames=c("dt", "id"), value.name=name))
Xs <- cbind(Xs[[1]], Xs[[2]][,3], Xs[[3]][,3], Xs[[4]][,3], Xs[[5]][,3])
colnames(Xs)[4:7] <- c("humid", "ws", "pres", "precip")
data_grid[, c("temp", "humid", "ws", "pres", "precip")] <- Xs %>% arrange(id, dt) %>% select(-dt, -id)

save(data_grid, file="data/modelData.RData")
