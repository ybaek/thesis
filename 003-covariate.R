#
# Processing dataset of other weather covariates
# (Final data structure is made)
#

library(dplyr)
library(ggplot2)
library(GGally)
library(lubridate)

load("data/dust.RData") ## dim. 2943420
load("data/Stations_dust.RData") ## 32 unique ID's
load("data/Stations_ASOS.RData")
load("data/Stations_AWS.RData")

# ASOS data + AWS data (see notes for location id extrapolation details)
# Slight deviations in latitude/longitude/elevation exist

########
# 1. Reading the data files
source("000-utilities.R")

lasos <- list()
laws <- list()
for (y in 2006:2017) {
  lasos[[(y-2005)]] <- mergecsv_fn("data/ASOS_data/", as.character(y), 7, "_")
  laws[[(y-2005)]] <- mergecsv_fn("data/AWS_data/", as.character(y), 7, "_")
}

# 2. Preprocessing
asos <- do.call(rbind.data.frame, lasos)
asos <- asos[,c(1:7, 9:10)]
asos <- asos %>% 
  transmute(id=지점, dt=as.POSIXct(일시, tz="UTC"), temp=`기온..C.`,
            precip=`강수량.mm.`, ws=`풍속.m.s.`, wd=`풍향.16방위.`,
            humid=`습도...`, dew=`이슬점온도..C.`, pres=현지기압.hPa.) %>% 
  as_tibble()

aws <- do.call(rbind.data.frame, laws)
aws <- aws[,c(1:7, 9)] %>% 
  transmute(id=지점, dt=as.POSIXct(일시, tz="UTC"), temp=`기온..C.`,
            precip=`강수량.mm.`, ws=`풍속.m.s.`, wd=`풍향.deg.`,
            humid=`습도...`, pres=현지기압.hPa.) %>% 
  as_tibble()
aws$dew <- NA # No dew point info. for AWS

covars <- rbind.data.frame(asos, aws) %>% arrange(id, dt)
## Wind direction: 16(+1) directions
covars$wd <- cut(covars$wd, breaks=seq(0, 360+22.5, by=22.5), right=F,
                 labels=c("CV", "NNE", "NE", "ENE", "E", "ESE", "SE", "SSE", "S",
                          "SSW", "SW", "WSW", "W", "WNW", "NW", "NNW", "N"))
save(covars, file="data/ASOS_AWS.RData") ## dim. 3162057
#####

# Collating response and covariance (The whole observation grid)
dust_tbl_sub <- dust_tbl %>% 
  filter(id!=144) %>% # exclude Oseong-san
  filter(year>=2006, year<=2017) ## dim. 2500826

covars$id <- ifelse(covars$id==695, 94, covars$id) # Gwangdeok-san
covars$id <- ifelse(covars$id==129, 132, covars$id) # Seosan => Ahnmyeon-do
covars$id <- ifelse(covars$id==159, 160, covars$id) # Busan => Gudeok-san
covars <- covars %>% arrange(id, dt)

t_start <- as.POSIXct("2006-01-01 00:00:00", tz="UTC")
t_end <- as.POSIXct("2017-12-31 23:00:00", tz="UTC")
ts <- seq(t_start, t_end, 3600)
## 105192==24 * (365 * 9 + 366 * 3)
locs <- unique(dust_tbl_sub$id) ## 31 unique locations (w/ differing elevations)
st_grid <- expand.grid(id=locs, dt=ts) %>% arrange(id, dt) %>% as_tibble() # dim. 3260952


data_grid <- left_join(st_grid, covars, by=c("id", "dt"))
data_grid <- left_join(data_grid, (dust_tbl %>% select(id, dt, lon, lat, elevation, log_dust, season)),
                  by=c("id", "dt"))
save(data_grid, file="data/Data.RData") # Final data structure, with all missing values

all_missing_y <- data_grid %>% group_by(dt) %>% 
  summarise(no_y = all(is.na(log_dust))) %>% filter(no_y) # 25 hrs
# ~23.3% missing values in the response (760127)

# For pressure/precipitation, need to distinguish btw. actual and non-actual
# missing values (e.g., little pres. change/trace values)
covars$flag_temp <- ifelse(is.na(covars$temp), 1, 0)
covars$flag_humid <- ifelse(is.na(covars$humid), 1, 0)
covars$flag_ws <- ifelse(is.na(covars$ws), 1, 0)
covars$flag_wd <- ifelse(is.na(covars$wd), 1, 0)
covars$flag_pres <- ifelse(is.na(covars$pres), 1, 0)
covars$flag_precip <- ifelse(is.na(covars$precip), 1, 0)
tb1 <- left_join((data_grid %>% select(id, dt)), covars, by=c("id", "dt"))
data_grid$flag_temp <- tb1$flag_temp
data_grid$flag_humid <- tb1$flag_humid
data_grid$flag_ws <- tb1$flag_ws
data_grid$flag_wd <- tb1$flag_wd
data_grid$flag_pres <- tb1$flag_pres
data_grid$flag_precip <- tb1$flag_precip
rm(tb1)

# All missing values are different for covariates (most much less severe)
# There are 103976 (~3%) of total hours omitted (will be reduced if we replace
# Bukchuncheon with Chuncheon data).

# Below is OUTDATED--wait till further repsonse

# sum(is.na(data_grid$temp)) ## 105693 (~3%)
# sum(is.na(data_grid$humid)) ## 222865 (~6%)
# sum(is.na(data_grid$ws)) ## 109888 (~3%)
# sum(is.na(data_grid$wd)) ## 109933 (~3%)
# sum(is.na(data_grid$flag_pres)) 
# sum(is.na(data_grid$flag_precip)) 
# 
# # When covariate is missing but response is not
# sum(is.na(data_grid$temp) & !is.na(data_grid$log_dust)) # 8965
# sum(is.na(data_grid$humid) & !is.na(data_grid$log_dust)) # 82066
# sum(is.na(data_grid$ws) & !is.na(data_grid$log_dust)) # 11966
# sum(is.na(data_grid$wd) & !is.na(data_grid$log_dust)) # 12010
# sum(is.na(data_grid$pres) & !is.na(data_grid$log_dust)) # 315745
# sum(is.na(data_grid$flag) & !is.na(data_grid$log_dust)) # 292435
# 
# # When covariate and response are both missing
# sum(is.na(data_grid$temp) & is.na(data_grid$log_dust)) # 96728
# sum(is.na(data_grid$humid) & is.na(data_grid$log_dust)) # 140799
# sum(is.na(data_grid$ws) & is.na(data_grid$log_dust)) # 97922
# sum(is.na(data_grid$wd) & is.na(data_grid$log_dust)) # 97923
# sum(is.na(data_grid$pres) & is.na(data_grid$log_dust)) # 294121
# sum(is.na(data_grid$flag) & is.na(data_grid$log_dust)) # 125120

# Is there a station where all or near-all covariates are missing?
tbl1 <- data_grid %>% group_by(id) %>%
  summarise(temp = sum(is.na(temp)),
            humid = sum(is.na(humid)),
            pres = sum(is.na(pres)),
            ws = sum(is.na(ws)), wd=sum(is.na(wd)),
            precip = sum(is.na(flag)))
## 북춘천(94)은 시작일이 제일 늦으므로 결측값이 제일 많다. (94225)
## 관측 시작 이후부터는 결측값이 없음.
tbl1 %>% arrange(desc(temp)) # 연평도(501)->5013
tbl1 %>% arrange(desc(humid)) # 연평도(501)->61224
tbl1 %>% arrange(desc(pres)) # 관악산(116)->102456
tbl1 %>% arrange(desc(ws)) # 연평도(501)->6041
tbl1 %>% arrange(desc(wd)) # 연평도(501)->6040
tbl1 %>% arrange(desc(precip)) # 광덕산(94)/안면도(132)/구덕산(160)->105192

#####

# Graphics
# (Fixed)


# Effect of westerly wind on high dust levels
## Direction: CV==(wind speed < 0.4m/s)
total_2017_sub <- data_grid %>% filter(year(dt)==2017, log_dust!=log(.001))
high_by_wd <- total_2017_sub %>% filter(log_dust > log(270)) %>% ## 3 sd above mean
  group_by(wd) %>% count() %>% arrange(wd) # All spring

g5 <- ggplot(high_by_wd, aes(wd, n)) +
  geom_col(aes(fill=wd), show.legend=F) + 
  #geom_text(aes(label=scales::percent(n/85))) +
  coord_polar("x", start=-0.235) + labs(x="", y="", title=expression(paste("Wind direction when >270", mu, "g/m"^3))) + 
  scale_y_continuous(breaks=NULL) + scale_color_hue(h=c(0,360))
ggsave("figures/eda/fig16-wind_direction_class.png", g5, device="png")
## West and Southwest winds account for ~42% (Comprise half with CV)

g6 <- ggplot(total_2017_sub, aes(log_dust)) +
  geom_density(aes(group=(wd %in% c("NW", "WNW", "WSW", "SW")), fill=(wd %in% c("NW", "WNW", "WSW", "SW"))), colour="gray", alpha=.5) +
  labs(x=expression(paste("PM"[10], " (log)")), fill="Wind direction") + 
  scale_fill_discrete(labels=c("Non-W", "Westerly"))
ggsave("figures/eda/fig17-density_direction.png", g6, device="png")

#####
load("data/FinalData.RData")
total_2017_sub <- data_grid %>% filter(year(dt)==2017, log_dust!=log(.001))
l <- list()
for (i in 1:4) l[[i]] <-total_2017_sub[total_2017_sub$season==levels(total_2017_sub$season)[i],]

# Pressures can be "cut"


plot_covar <- function(df, x) {
  ggplot(df, aes_string(x, "log_dust")) +
  geom_point() +
  geom_smooth(aes(group=as.factor(wd), col=as.factor(wd)), se=F, show.legend=F)
}

plotList1 <- lapply(l, function(df) plot_covar(df, "temp"))
plotList2 <- lapply(l, function(df) plot_covar(df, "humid"))
plotList3 <- lapply(l, function(df) plot_covar(df, "ws"))
plotList4 <- lapply(l, function(df) plot_covar(df, "pres"))
plotList5 <- lapply(l, function(df) plot_covar(df, "precip"))

g7 <- ggmatrix(plotList1, 2, 2, title=expression(paste("Air temperature(", degree, "C)")))
g8 <- ggmatrix(plotList2, 2, 2, title="Humidity (%)")
g9 <- ggmatrix(plotList3, 2, 2, title="Wind speed (m/s)")
g10 <- ggmatrix(plotList4, 2, 2, title="Air pressure (hPa)")
g11 <- ggmatrix(plotList5, 2, 2, title="Precipitation (mm)")
ggsave("figures/eda/fig18-temp_17.png", g7, device="png", dpi=150)
ggsave("figures/eda/fig19-humid_17.png", g8, device="png", dpi=150)
ggsave("figures/eda/fig20-ws_17.png", g9, device="png", dpi=150)
ggsave("figures/eda/fig21-pres_17.png", g10, device="png", dpi=150)
ggsave("figures/eda/fig22-precip_17.png", g11, device="png", dpi=150)


### Southwesterly wind seems to have a peculiar effect.
### Interestingly, West + Northwesterly winds had more prominent effect in 2007
### (For > 800, accounted for ~64% of all episodes)
### It seems consistent that westerly wind has significant effect on increased PM10 levels



