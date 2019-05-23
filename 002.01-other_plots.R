#
# Some final graphic work
#

load("data/FinalData.RData")
load("data/Stations_ASOS.RData")

library(ggplot2)
library(dplyr)
library(lubridate)
library(reshape2)

korea <- map_data("world", "South Korea")

final_stns <- data_grid %>% group_by(id, lon, lat, elevation) %>% count()
asos_stns_model <- asos_stn %>% filter(id %in% unique(final_stns$id) | id==129 | id==159) %>% 
  select(id, lon, lat, elevation)

g1 <- ggplot() +
  geom_polygon(data=korea, aes(long, lat, group=group), fill="gray") +
  geom_point(data=asos_stns_model, aes(lon, lat), fill="blue",
             colour="black", pch=21, size=2.75) +
  geom_point(data=final_stns, aes(lon, lat, fill=elevation),
             colour="black", size=2.25, pch=21) +
  geom_text(data=final_stns, aes(lon, lat, label=id), 
            nudge_y=0.15, size=2.75, fontface="bold") +
  scale_fill_distiller(palette="Reds", direction=1) +
  labs(fill="Elevation (m)", 
       title="Data Locations Used for Modeling") +
  theme_dark()

ggsave("figures/eda/fig26-final_locations.png", g1, device="png", width=8, height=6,
       units="in", dpi=300)

asos_stns_model$id[12] <- 132 # 서산 => 안면도
asos_stns_model$id[27] <- 160 # 부산 => 구덕산

#####

agg_grid <- data_grid %>% 
  group_by(id, year=year(dt), month=month(dt)) %>% 
  summarise(y = mean(log_dust, na.rm=T)) %>% 
  ungroup()

g2 <- ggplot(agg_grid, aes(month, as.factor(id))) +
  geom_tile(aes(fill=y), colour="grey50") +
  facet_wrap(~ year, nrow=3) +
  labs(x="Month", y="ID", 
       fill=expression(paste("log ", mu, "g/m"^3)),
       title="Monthly Average Levels per Station") +
  scale_x_continuous(breaks=1:12) +
  scale_fill_distiller(palette="Spectral") +
  theme(axis.text.x = element_text(angle=45, vjust=.7, colour="black"),
        axis.text.y = element_text(colour="black", size=7.5),
        panel.background = element_blank(),
        strip.text.x = element_text(margin=margin(0.05,0,0.05,0, "cm")))
ggsave("fig1-monthly_station.png", path="images/", g2, device="png",
       width=8, height=6, units="in")  


#####

# For elevation differences between ASOS and PM10 stations...

#####

# ST-variogram

library(gstat)
library(lubridate)
library(sp)
library(spacetime)

wgs84_str <- "+proj=longlat +datum=WGS84"

agg_grid <- data_grid %>% group_by(id, year=year(dt), month=month(dt)) %>% 
  summarise(y = mean(log_dust, na.rm=T)) %>% 
  arrange(year, month)
coords <- unique(data_grid %>% group_by(id, lon, lat) %>% count())[-c(4,6,9,20,29),1:2]
time <- as.POSIXct(paste(rep(2006:2017, each=12), 
              rep(c(paste0("0", 1:9), 10:12), 12), 
              "01", sep="-"), tz="UTC")
data <- agg_grid %>% select(y)
STFDF <- STFDF(
  sp = SpatialPoints(coords, proj4string=CRS(wgs84_str)),
  time = time, data = data
)
vv <- variogramST(y ~ 1, data=STFDF, assumeRegular=TRUE)
png("figures/model/fig2-variogram.png", width=6, height=6, units="in", res=250)
plot(vv, main="Monthly Mean Variogram")
dev.off()

#####
# Correlation btw. covariates
corMat <- data_grid %>% select(log_dust, temp, precip, ws, humid, pres, elevation) %>% 
  as.matrix() %>% cor(use="na.or.complete", method="spearman")
g3 <- ggplot(data=melt(corMat), aes(Var1, Var2)) +
  geom_tile(aes(fill=value)) + 
  labs(x="", y="", fill="Spearman\nCorrelation", title="Before Interpolation") +
  scale_fill_distiller(palette="RdBu", limits=c(-1,1)) + 
  theme(axis.text.x = element_text(colour="black"), axis.text.y = element_text(colour="black"),
        panel.background = element_blank())
ggsave("fig3-correlation_plot.png", path="figures/model/", g3, device="png",
       width=8, height=6, units="in")  

# After interpolation
load("data/modelData.RData")
corMat <- data_grid %>% select(log_dust, temp, precip, ws, humid, pres, elevation) %>% 
  as.matrix() %>% cor(use="na.or.complete", method="spearman")
g4 <- ggplot(data=melt(corMat), aes(Var1, Var2)) +
  geom_tile(aes(fill=value)) + 
  labs(x="", y="", fill="Spearman\nCorrelation", title="After Interpolation") +
  scale_fill_distiller(palette="RdBu", limits=c(-1,1)) + 
  theme(axis.text.x = element_text(colour="black"), axis.text.y = element_text(colour="black"),
        panel.background = element_blank())
ggsave("fig4-correlation_plot_post.png", path="figures/model/", g4, device="png",
       width=8, height=6, units="in")  

# Space-only variograms
## NOT REPLICABLE 
## vv (Apr. 10th)
library(gstat)
variogram_data1 <- real_data %>% filter(year(dt)==2017, month(dt)==4, day(dt)==10) %>%
  group_by(id, lon, lat) %>%
  summarise(mean = mean(log_dust, na.rm=T)) %>%
  filter(!is.nan(mean))
variogram_data2 <- real_data %>% filter(year(dt)==2017, month(dt)==9, day(dt)==1) %>%
  group_by(id, lon, lat) %>%
  summarise(mean = mean(log_dust, na.rm=T)) %>%
  filter(!is.nan(mean))
variogram_data1 <- variogram_data1[-4,]
variogram_data2 <- variogram_data2[-4,]
spdf1 <- SpatialPointsDataFrame(coords=variogram_data1[,c(2,3)],
                       data=variogram_data1[,4],
                       proj4string = CRS(wgs84_str))
spdf2 <- SpatialPointsDataFrame(coords=variogram_data2[,c(2,3)],
                       data=variogram_data2[,4],
                       proj4string = CRS(wgs84_str))
vv1 <- variogram(mean ~ 1, data=spdf1)
vv2 <- variogram(mean ~ 1, data=spdf2)
vgm1 <- fit.variogram(vv1, vgm(1, "Exp", 500, 1))
vgm2 <- fit.variogram(vv1, vgm(1, "Sph", 500, 1))

png("presentation/images/sp_variograms1.png", width=7, height=6, units="in", res=250)
p1 <- plot(vv1, model=vgm1, pch=19, lwd=2, main="Exponential fit")
p2 <- plot(vv1, model=vgm2, pch=19, lwd=2, main="Spherical fit")
print(p1, split = c(1, 1, 2, 1), more = TRUE)
print(p2, split = c(2, 1, 2, 1))
dev.off()

png("presentation/images/sp_variograms2.png", width=7, height=6, units="in", res=250)
p3 <- plot(vv1, model=vgm1, pch=19, lwd=2, main="Apr. 1")
p4 <- plot(vv2, model=vgm2, pch=19, lwd=2, col="red", main="Sep. 1")
print(p3, split = c(1, 1, 2, 1), more = TRUE)
print(p4, split = c(2, 1, 2, 1))
dev.off()