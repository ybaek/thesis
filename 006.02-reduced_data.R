#
# Reduced dataset
#

library(dplyr)
library(ggplot2)
library(gstat)
library(igraph)
library(imputeTS)
library(lubridate)
library(reshape2)
library(sp)
library(spacetime)
library(SpatioTemporal)

load("data/ReducedData.RData")
dim(real_data[!is.na(real_data$log_dust),]) #596373

stn_data <- real_data %>% group_by(id, lon, lat, elevation) %>% summarise()

korea <- map_data("world", "South Korea")
g1 <- ggplot() +
  geom_polygon(data=korea, aes(long, lat, group=group), fill="gray", col="black") +
  geom_point(data=stn_data, aes(lon, lat, fill=elevation), alpha=.8, shape=21, colour="black", size=2.75) +
  geom_text(data=stn_data, aes(lon, lat, label=id), nudge_y=0.125, size=2.75, fontface="bold") +
  labs(fill="Elevation (m)") + 
  scale_fill_distiller(palette="Reds", direction=1) + theme_dark()
ggsave("fig3-map.png", g1, device="png", path = "figures/eda/")


# Season re-edit
real_data$season <- rep("", dim(real_data)[1])
real_data$season <- ifelse((month(real_data$dt) >= 3 & month(real_data$dt) <=5), "spring", real_data$season)
real_data$season <- ifelse((month(real_data$dt) >= 6 & month(real_data$dt) <=8), "summer", real_data$season)
real_data$season <- ifelse((month(real_data$dt) >= 9 & month(real_data$dt) <=11), "fall", real_data$season)
real_data$season <- ifelse((month(real_data$dt) >= 1 & month(real_data$dt) <=2 | month(real_data$dt)==12), "winter", real_data$season)
real_data$season <- factor(real_data$season, levels=c("spring", "summer", "fall", "winter"))

# Overall Monthly data plot
forPlot <- real_data %>% 
  group_by(id, year=year(dt), month=month(dt)) %>% 
  summarise(y = mean(log_dust, na.rm=T)) %>% 
  ungroup()
forPlot$t <- rep(1:36, 25) # 12 months * 3 years * 25 stns

g2 <- ggplot(forPlot, aes(month, as.factor(id))) +
  geom_tile(aes(fill=y), colour="grey50") +
  facet_wrap(~ year, nrow=3, ncol=1) +
  labs(x="Month", y="", 
       fill=expression(paste("log ", mu, "g/m"^{3}))) +
  scale_x_continuous(breaks=1:12) +
  scale_fill_distiller(palette="YlGnBu", direction=1) +
  theme(axis.text.y = element_text(colour="black", size=6),
        panel.background = element_blank(),
        strip.text.x = element_text(margin=margin(0.05,0,0.05,0, "cm")))
ggsave("fig1-monthly_station.png", path="figures/eda/", g2, device="png",
       width=7, height=6, units="in")  

# Variability across seasons
forPlot <- real_data %>% mutate(year=year(dt)) %>% 
  filter(!is.na(season), log_dust > log(.001))
g4 <- ggplot(forPlot, aes(season, log_dust)) + 
  geom_boxplot(aes(fill=season), outlier.alpha=.1) +
  facet_wrap(~year, nrow=1, ncol=3) + labs(x="", y="log levels") +
  scale_fill_hue(h=c(0, 250)) +
  scale_x_discrete(breaks=NULL) +
  theme(legend.margin=margin(0,0,0,0,"pt"))

ggsave(g4, file="fig9-season_boxplot.png", path="figures/eda/", device="png",
       width=7, height=5, units="in")

# Variogram
wgs84_str <- "+proj=longlat +datum=WGS84"
agg_grid <- real_data %>% group_by(id, year=year(dt), month=month(dt)) %>% 
  summarise(y = mean(log_dust, na.rm=T)) %>% 
  arrange(year, month)
coords <- unique(real_data %>% group_by(id, lon, lat) %>% count())[-c(4,6,9,20,29),1:2]
time <- as.POSIXct(paste(rep(2015:2017, each=3), 
                         rep(c(paste0("0", 1:9), 10:12), 3), 
                         "01", sep="-"), tz="UTC")
time <- time[order(time)]
data <- agg_grid %>% ungroup() %>% select(y)
STFDF <- STFDF(
  sp = SpatialPoints(coords, proj4string=CRS(wgs84_str)),
  time = time, data = data
)
vv <- variogramST(y ~ 1, data=STFDF, assumeRegular=TRUE)
png("figures/eda/fig15-variogram.png", width=5, height=5, units="in", res=250)
plot(vv, main="Monthly Mean Variogram")
dev.off()

# ACF (short-lag)
# How seasonality can be accounted for by predictors
my_acast <- function(var) acast(real_data, id ~ dt, value.var=var)
response <- my_acast("log_dust")
y_acfs <- apply(response, 1, function(r) acf(r, lag.max=150, plot=F, na.action=na.pass))

png("figures/eda/fig10-response_acfs.png", width=8, height=8, units="in", res=300)
par(mfrow=c(5,5), mar=c(3,2,2,0.8), oma=c(1,1.2,1,1), mgp=c(1.5,0.6,0))
lapply(y_acfs, function(x) plot(x, ylab="", main=""))
dev.off()

# Predictors: Air temp, ws, humidity, air pres.
Xs <- lapply(list("temp","ws","humid","pres"), my_acast)
X_acfs <- lapply(Xs, 
                 function(l) apply(l, 1, 
                                   function(r) acf(r, lag.max=150, plot=F, na.action=na.pass)))
png("figures/eda/fig11-air_temp_acfs.png", width=8, height=8, units="in", res=300)
par(mfrow=c(5,5), mar=c(3,2,2,0.8), oma=c(1,1.2,1,1), mgp=c(1.5,0.6,0))
lapply(X_acfs[[1]], function(x) plot(x, ylab="", main=""))
dev.off()

png("figures/eda/fig12-wind_speed_acfs.png", width=8, height=8, units="in", res=300)
par(mfrow=c(5,5), mar=c(3,2,2,0.8), oma=c(1,1.2,1,1), mgp=c(1.5,0.6,0))
lapply(X_acfs[[2]], function(x) plot(x, ylab="", main=""))
dev.off()

png("figures/eda/fig13-humidity_acfs.png", width=8, height=8, units="in", res=300)
par(mfrow=c(5,5), mar=c(3,2,2,0.8), oma=c(1,1.2,1,1), mgp=c(1.5,0.6,0))
lapply(X_acfs[[3]], function(x) plot(x, ylab="", main=""))
dev.off()

png("figures/eda/fig14-air_pres_acfs.png", width=8, height=8, units="in", res=300)
par(mfrow=c(5,5), mar=c(3,2,2,0.8), oma=c(1,1.2,1,1), mgp=c(1.5,0.6,0))
lapply(X_acfs[[4]], function(x) plot(x, ylab="", main=""))
dev.off()

# Association between predictors and response
with(real_data, cor(temp,log_dust,use="na.or.complete",method="spearman")) #-0.116
# Temperature is peculiar.
airtemp_mat <- Xs[[1]]
airtemp_mat <- t(apply(airtemp_mat, 1, function(row) na.interpolation(row, option="stine")))
## otherwise doesn't work w/ NA
temp_seas <- t(apply(airtemp_mat, 1, 
                    function(row) stats::decompose(ts(row, frequency=24))$seasonal))
cor(c(t(temp_seas)), real_data$log_dust, use="na.or.complete", method="spearman") ##-0.137
## Not much difference. Why negative direction?

with(real_data, cor(ws,log_dust,use="na.or.complete",method="spearman")) # ~null
with(real_data, cor(humid,log_dust,use="na.or.complete",method="spearman")) #-0.2
with(real_data, cor(pres,log_dust,use="na.or.complete",method="spearman")) #0.15
## All other directions match, but only very mild strength.

# Graph structure
source("000-utilities.R")
forPlot <- unique(real_data[,c("id","lon","lat")])[-c(3,5,8,19,28),2:3]
distMat <- distMat_fn(as.matrix(forPlot))
rownames(distMat) <- unique(real_data$id)
colnames(distMat) <- unique(real_data$id)

## Similarity Matrix
rho <- 100000
cutoff <- 165000
distMat <- exp(-distMat/rho)
distMat <- ifelse(distMat > exp(-cutoff/rho), distMat, 0) # temp. cutoff
diag(distMat) <- 0
## Forming a graph (data frame==edge list)
el <- as.data.frame(t(combn(unique(real_data$id), 2)))
el$weight <- 0
for (i in 1:nrow(el)) {
  n1 <- as.character(el[i,1]); n2 <- as.character(el[i,2])
  if (distMat[n1,n2]!=0) el[i,3] <- distMat[n1, n2]
}
el <- el[el$weight!=0,]
graph <- graph_from_data_frame(el, directed=FALSE)

png("figures/model/fig1-graph.png", width=8, height=8, units="in", res=250)
plot(graph, layout=layout_(graph, nicely()))
dev.off()

