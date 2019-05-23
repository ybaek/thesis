#
# Preprocessing (before Modeling)
# Interpolating X -- the covariates
# CAUTION MUST BE TAKEN IN INFERENCE
#

library(dplyr)
library(ggplot2)
library(GGally)
library(lubridate)

load("data/Data.RData")
load("data/Stations_dust.RData")
# There must be no NA's for: Latitude + Longitude + Elevation 
for (i in unique(data_grid$id)) {
  if (i %in% c(108, 115, 116, 152, 232)) {
    sub <- filter(data_grid, id==i) %>% select(dt, lon, lat, elevation)
    # dustStn_tbl must be ordered by end dates
    t <- dustStn_tbl[dustStn_tbl$id==i,"end"][[1]][1]
    sub1 <- filter(sub, as.Date(dt) < t)
    sub2 <- filter(sub, as.Date(dt) >= t)
    sub1$lon <- dustStn_tbl[dustStn_tbl$id==i,]$lon[1]
    sub1$lat <- dustStn_tbl[dustStn_tbl$id==i,]$lat[1]
    sub1$elevation <- dustStn_tbl[dustStn_tbl$id==i,]$elevation[1]
    sub2$lon <- dustStn_tbl[dustStn_tbl$id==i,]$lon[2]
    sub2$lat <- dustStn_tbl[dustStn_tbl$id==i,]$lat[2]
    sub2$elevation <- dustStn_tbl[dustStn_tbl$id==i,]$elevation[2]
    data_grid[(data_grid$id==i & as.Date(data_grid$dt) < t), c("lon", "lat", "elevation")] <- sub1[,c("lon", "lat", "elevation")]
    data_grid[(data_grid$id==i & as.Date(data_grid$dt) >= t), c("lon", "lat", "elevation")] <- sub2[,c("lon", "lat", "elevation")]
  }
  else {
    sub <- filter(data_grid, id==i) %>% select(lon, lat, elevation)
    sub$lon <- dustStn_tbl[dustStn_tbl$id==i,]$lon
    sub$lat <- dustStn_tbl[dustStn_tbl$id==i,]$lat
    sub$elevation <- dustStn_tbl[dustStn_tbl$id==i,]$elevation
    data_grid[data_grid$id==i, c("lon", "lat", "elevation")] <- sub
  }
}
# save(data_grid, file="data/Data.RData")

# Handling pressure and precip.
for (i in unique(data_grid$id)) {
  sub <- filter(data_grid, id==i) %>% select(pres, precip, flag_pres, flag_precip)
  for (r in 2:nrow(sub)) {
    pres <- sub[r,"pres"][[1]]
    precip <- sub[r,"precip"][[1]]
    f1 <- sub[r,"flag_pres"][[1]]
    f2 <- sub[r,"flag_precip"][[1]]
    if (is.na(f1) | is.na(f2)) next
    if (is.na(pres) & f1==1) sub[r,]$pres <- sub[(r-1),]$pres
    if (is.na(precip) & f1==1) sub[r,]$pres <- 0
  }
  data_grid[data_grid$id==i, c("pres", "precip")] <- sub[,c("pres", "precip")]
}


# Time series plot may draw out the effect direction between covariate and response
daily_agg <- data_grid %>%
  group_by(id, year=year(dt), month=month(dt), day=day(dt)) %>% 
  summarise(ld = mean(log_dust, na.rm=T), temp = mean(temp, na.rm=T), humid = mean(humid, na.rm=T),
            ws = mean(ws, na.rm=T), pres = mean(pres, na.rm=T), precip = mean(precip, na.rm=T))

# Graphics
base_g <- ggplot((daily_agg %>% filter(year==2007, id!=501) %>% arrange(month, day)), aes(x=paste(month, day), group=id, col=id)) +
  labs(x="") + scale_x_discrete(breaks=NULL, labels=NULL) + scale_color_distiller(palette="YlOrRd", direction=1)

g1 <- base_g + geom_line(aes(y = ld), show.legend=F) 
g2 <- base_g + geom_line(aes(y = temp), show.legend=F)
g3 <- base_g + geom_line(aes(y = humid), show.legend=F)
g4 <- base_g + geom_line(aes(y = ws), show.legend=F)
g5 <- base_g + geom_line(aes(y = log(pres)), show.legend=F)
g6 <- base_g + geom_line(aes(y = precip), show.legend=F)

gList1 <- ggmatrix(list(g1, g5, g3, g4, g6), 5, 1, showAxisPlotLabels=T, 
         yAxisLabels=c("log PM10", "log hPa", "humid (%)", "wind (m/s)", "precip (mm)")) +
  ggtitle("Year 2007")
# Think about detrending temperature

#####
# data_grid$wd <- ifelse(data_grid$wd %in% c("SW", "WSW", "W", "WNW", "NW"),
#                        "W", data_grid$wd) # regime: Westerly vs. not.

## Things to think about
## 1. How to handle station with covariates observed only after 2006
## -> Only one case: 북춘천 (93). How to handle this case?
## -> Drew on an +1 hr data plus some jitters.

## 2. Certain stations have disproportionately large no. of NA's for certain 
## variables only. This will likely cause problems to interpolation.
## (NOT yet fully resolved)

## 3. Data quality issues (e.g., why so many NA's for air pressure?)
## -> Waiting for KMA correspondence. Likely guess: air pres and precip
## Should not be a big issue.

load("data/ASOS_AWS.RData")
sub <- covars %>% filter(id==93 | id==101) %>% 
  filter(dt > as.POSIXct("2016-10-01 00:00:00", tz="UTC"))
base_g <- ggplot(sub, aes(x=dt, group=id, col=as.factor(id)))
g1 <- base_g + geom_line(aes(y=temp))
g2 <- base_g + geom_line(aes(y=humid))
g3 <- base_g + geom_line(aes(y=ws))
g4 <- base_g + geom_line(aes(y=pres))
g5 <- base_g + geom_line(aes(y=precip))
gList2 <- ggmatrix(list(g1, g2, g3, g4, g5), 5, 1, showAxisPlotLabels=T, 
                   yAxisLabels=c("Temp", "Humid", "Speed", "Pres", "Precip")) +
  ggtitle("Comparing weather for IDs 93 and 101")

sub93 <- sub %>% filter(id==93)
sub101 <- sub %>% filter(id==101)

# Lag maximizing CCF was ~ -1 hr (except precipitation)
my_ccf_fn <- function (var) ccf(sub93[[var]], sub101[[var]], lag.max=4500, na.action=na.pass, plot=F) 
l <- list("temp", "humid", "ws", "pres", "precip")
ccfList <- lapply(l, my_ccf_fn)

hist(sub93$temp - sub101$temp) # very std. normal-like
hist(sub93$humid - sub101$humid)
hist(sub93$pres - sub101$pres) # not centered around zero
hist(sub93$ws - sub101$ws) # very std. normal-like
hist(sub93$precip - sub101$precip)

png("figures/eda/fig25-CCF.png", width=7, height=4, units="in", res=300)
par(mfrow=c(2,3), mar=c(3,4,3,2)+.1, mgp=c(2, 1, 0))
plot(ccfList[[1]], main="Temp")
plot(ccfList[[2]], main="Humid")
plot(ccfList[[3]], main="Wind Speed")
plot(ccfList[[4]], main="Air Pres")
plot(ccfList[[5]], main="Precipitation")
dev.off()
# rm(sub)
# rm(sub93)
# rm(sub101)

# When we choose to interpolate, we can use data from y. 2015, Dec.
# ID=101, for starting points
start_93 <- as.POSIXct("2006-01-01 00:00:00", tz="UTC") + 3600
end_93 <- as.POSIXct("2016-10-01 00:00:00", tz="UTC") + 3600
interp101 <- data_grid %>% filter(id==101) %>% select(-lon, -lat, -elevation, -log_dust, -season) %>% 
  filter(dt >= start_93, dt <= end_93) %>% arrange(dt)
interp93 <- data_grid %>% filter(id==93) %>% select(-lon, -lat, -elevation, -log_dust, -season) %>% 
  filter(dt <= as.POSIXct("2016-10-01 00:00:00", tz="UTC"))
interp93[,3:15] <- interp101[,3:15] # 춘천 data => 북춘천 data with lag=+1 hr.
# Need to add jitters here... (see histograms above)
save(interp93, file="data/ID93_interp.RData")
# rm(interp101)

### UPDATE (1/15): NOT taking this appraoch. Just merge them as one station.
