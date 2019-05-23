
library(dplyr)
library(geosphere)
library(ggplot2)
library(lubridate)
library(reshape2)

load("data/ReducedData.RData")
source("code/000-utilities.R")

# Examining distance and distance matrix
## WGS84 ellipsoid distance matrix
stations <- real_data %>% arrange(dt) %>% 
  select(id, lon, lat) %>% unique() %>% head(25)
coords_mat <- stations %>% select(lon, lat) %>% 
  as.matrix()
dist_mat <- matrix(NA, 25,25)

for (i in 1:(nrow(coords_mat)-1)) {
  dist_mat[i,i] <- 0 ## Same point has distance zero
  for (j in (i+1):nrow(coords_mat)) { 
    d <- distGeo(coords_mat[i, c("lon", "lat")],
                 coords_mat[j, c("lon", "lat")])
    dist_mat[i,j] <- d
    dist_mat[j,i] <- d
  }
  dist_mat[j,j] <- 0
}
rownames(dist_mat) <- stations$id
colnames(dist_mat) <- stations$id
dist_mat_for_plot <- dist_mat/100000
dist_mat_for_plot[upper.tri(dist_mat)] <- NA
dist_mat_for_plot <- melt(dist_mat_for_plot) %>% 
  mutate(Var1=as.factor(Var1), Var2=as.factor(Var2), value=value)
  
g1 <- ggplot(dist_mat_for_plot, aes(Var1, Var2)) + geom_tile(aes(fill=value)) +
  scale_fill_distiller(palette="Spectral", direction=-1, na.value=NA) +
  labs(x="ID", y="ID", fill="Distance\n(1e5 m)") +
  theme(axis.text.x = element_text(angle=45, vjust=.7, colour="black"), 
        axis.text.y = element_text(colour="black"),
        panel.background = element_rect(fill="white"),
        legend.margin=margin(0,0,0,0,"pt"))
ggsave(g1, filename="fig5-distance_heat.png", device="png", 
       path="images", width=5, height=5, dpi=320)

# Correlation differs from distance (non-parametric-ness of spatial)
real_mat <- acast(real_data[,c(1,2,13)], id ~ dt)
corMat <- cor(t(real_mat), use = "pairwise.complete.obs", method="spearman")
corMat_for_plot <- corMat
corMat_for_plot[upper.tri(corMat)] <- NA
corMat_for_plot <- melt(corMat_for_plot) %>%
  transmute(Var1=as.factor(Var1), Var2=as.factor(Var2), value=value)
g2 <- ggplot(corMat_for_plot, aes(Var1, Var2)) + geom_tile(aes(fill=value)) +
  scale_fill_distiller(palette="Spectral", direction=1, na.value=NA) +
  labs(x="ID", y="ID", fill="Spearman\nCorrelation") +
  theme(axis.text.x = element_text(angle=45, vjust=.7, colour="black"), 
        axis.text.y = element_text(colour="black"),
        panel.background = element_rect(fill="white"),
        legend.margin=margin(0,0,0,0,"pt"))
ggsave(g2, filename="fig6-corrleation_heat.png", device="png", 
       path="images", width=5, height=5, dpi=320)

#####
# An eigensystem associated to the Moran's I statistic
# Ex. A 5x5 grid.
library(igraph)
grid <- expand.grid(x=seq(.5, 4.5, by=1), y=seq(.5, 4.5, by=1))
adj_grid <- as_adj(make_lattice(c(5,5)), sparse = F) # denotes only neighbordhood structure
dist_grid <- as.matrix(stats::dist(grid)) # denotes the strength of that connection
cutoff <- median(unlist(dist_grid[dist_grid!=0]))
## ~0.3290306: median norm. distance (excluding 0's)
# ... (not using)
sim_grid3 <- exp(-dist_grid)
sim_grid3 <- ifelse(sim_grid3 < exp(-cutoff), 0, sim_grid3)
diag(sim_grid3) <- 0

E1 <- moran_fn(adj_grid) # ~10 (excluding singularities)
E4 <- moran_fn(sim_grid3) # ~7 (excluding singularities)
S1 <- E1$vectors 
S4 <- E4$vectors 
colnames(S1) <- paste0("basis1_", seq(1:25))
colnames(S4) <- paste0("basis4_", seq(1:25))
grid <- cbind(grid, S1, S4)

# Areal plot (mid-points)
rect_x1 = rep(c(0, 1, 2, 3, 4), 5)
rect_x2 = rect_x1 + 1
rect_y1 = rep(c(0, 1, 2, 3, 4), each = 5)
rect_y2 = rect_y1 + 1
g1 <- ggplot() + 
  geom_rect(aes(xmin = rect_x1, xmax = rect_x2, ymin = rect_y1, ymax = rect_y2), fill="white", color="black") +
  geom_point(data=grid, aes(x, y)) +
  theme_bw()
#ggsave("figures/simulation/areal_grid.png", g1, device="png", height=5.18, width=6.19, units="in")

simulPlot_fn <- function(var) {
  g <- ggplot(grid, aes(x, y)) + geom_tile(aes_string(fill=var), show.legend=F) + scale_fill_distiller(direction=1)
  #ggsave(paste0("figures/simulation/", var, ".png"), g, device="png")
}
Map(simulPlot_fn, as.list(c(colnames(S1), colnames(S4))))

## A "pseudo"-Laplacian matrix is guaranteed to be positive semidefinite.
Q <- diag(rowSums(sim_grid3), 25) - sim_grid3

# 
# ## 1.
# set.seed(0)
# for (i in 1:10^5) {
#   y <- rnorm(25)
#   I[i,1] <- moran_fn(y, dist_grid)
# }
# 
# ## 2. Set exp. covariogram parameters
# library(MASS)
# rho <- runif(1,0,10)
# sigma <- runif(1,0,10)
# mu <- rep(0, 25)
# Sigma <- sigma^2 * exp(-dist_grid/rho)
# set.seed(0)
# for (i in 1:10^5) {
#   y <- mvrnorm(1, mu, Sigma)
#   I[i,2] <- moran_fn(y, dist_grid)
# }
# 
# ## 3. Spatial correlation, not necessarily based on distance
# set.seed(0)
# Sigma <- rWishart(1, 25, diag(1, 25))[,,1]
# for (i in 1:10^5) {
#   y <- mvrnorm(1, mu, Sigma)
#   I[i,3] <- moran_fn(y, dist_grid)
# }
# # mean and variance are similar as if no independence
# 
# apply(I, 2, mean) ## theoretical mean== -0.04167
# # -0.04164629 -0.14324302 -0.04793842
# apply(I, 2, var) ## all lesser than theoretical variance of 0.0067
# # 0.0004660707 0.0037051618 0.0006813005
# 
# ## Graphics
# data_grid$xmin = rep(c(0, 1, 2, 3, 4), 5)
# data_grid$xmax = data_grid$xmin + 1
# data_grid$ymin = rep(c(0, 1, 2, 3, 4), each = 5)
# data_grid$ymax = data_grid$ymin + 1
# 
# set.seed(0)
# data_grid$random = rnorm(25)
# data_grid$exp = mvrnorm(1, mu, sigma^2 * exp(-dist_grid/rho))
# data_grid$nondist = mvrnorm(1, mu, Sigma)
# 
# ggplot(data_grid) +
#   geom_rect(aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="white", colour="black") +
#   geom_point(aes(x, y)) + theme_minimal()
# 
# g2 <- ggplot(data_grid) +
#   geom_tile(aes(x, y, fill=scale(random))) + labs(x="", y="", fill="z")
# g3 <- ggplot(data_grid) +
#   geom_tile(aes(x, y, fill=scale(exp))) + labs(x="", y="", fill="z")
# g4 <- ggplot(data_grid) +
#   geom_tile(aes(x, y, fill=scale(nondist))) + labs(x="", y="", fill="z")
# ggsave("figures/fig12-normal_sim.png", g2, "png")
# ggsave("figures/fig13-exp_covgm_sim.png", g3, "png")
# ggsave("figures/fig14-mvn_sim.png", g4, "png")
# 
# 
# g2 <- ggplot() + geom_histogram(aes(I[,1], ..density..), fill="lightblue", colour="black") +
#   xlab("Moran's I")
# g3 <- ggplot() + geom_histogram(aes(I[,2], ..density..), fill="lightblue", colour="black") +
#   xlab("Moran's I")
# g4 <- ggplot() + geom_histogram(aes(I[,3], ..density..), fill="lightblue", colour="black") +
#   xlab("Moran's I")
# ggsave("figures/fig15-normal_hist.png", g2, "png")
# ggsave("figures/fig16-exp_covgm_hist.png", g3, "png")
# ggsave("figures/fig17-mvn_hist.png", g4, "png")



