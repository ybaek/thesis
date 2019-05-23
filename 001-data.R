#
# Raw code for exploring dust level data in S. Korea
#


library(dplyr)
library(ggplot2)
library(readr)
library(lubridate)

dust_tbl <- read_csv("data/Dust_data.csv")
dustStn_tbl <- read_csv("data/Dust_stn_info.csv")

#####

#1. Processing data for PM10 lvls (response)

# Station IDs with changed positions
## It seems only names have changed for 132 and 229 (NEED TO CONFIRM)
dustStn_tbl <- dustStn_tbl[-c(17,34),]
dustStn_tbl[c(17,33),3] <- NA
changed_ids <- c(102, 108, 115, 116, 152, 232)

dustStn_tbl <- dustStn_tbl %>% arrange(지점, 종료일)

for (i in 1:nrow(dustStn_tbl)) {
  if (dustStn_tbl[i,1] %in% changed_ids & is.na(dustStn_tbl[i,3])) {
    dustStn_tbl[i,1] <- dustStn_tbl[i,1] - .5
  }
}

inds <- which(dust_tbl$지점 %in% changed_ids)
temp_df <- dust_tbl[inds,c(1,2)]
for (i in 1:nrow(temp_df)) {
  if (as.Date(temp_df[[i,2]]) > dustStn_tbl[dustStn_tbl$지점==temp_df[[i,1]],]$종료일[1]) {
    temp_df[i,1] <- temp_df[[i,1]] - .5
  }
}
dust_tbl[inds,1:2] <- temp_df

dust_tbl <- left_join(dust_tbl, dustStn_tbl, by="지점")
dust_tbl <- dust_tbl %>% 
  rename(dust = `1시간평균 미세먼지농도(㎍/㎥)`, elevation = `노장해발고도(m)`, dt = 일시,
         lon = 경도, lat = 위도, id = 지점, start = 시작일, end = 종료일, name = 지점명) %>% 
  mutate(dust = ifelse(dust <= -999, NA, dust))
dust_tbl <- dust_tbl %>% 
  mutate(log_dust = log(dust), year=as.integer(year(dt)), month=as.integer(month(dt)), 
         day=day(dt), hour = hour(dt))
dust_tbl$flag <- ifelse(dust_tbl$dust==0, "2", "0") # 0=="2"
dust_tbl$dust <- ifelse(dust_tbl$dust==0, 0.001, dust_tbl$dust) # add 0.001
dust_tbl$flag <- ifelse(is.na(dust_tbl$dust), "1", dust_tbl$flag) # NA=="1"
dust_tbl$log_dust <- log(dust_tbl$dust)
dust_tbl <- dust_tbl[,-c(7, 11:14)]
dust_tbl$id <- as.integer(ceiling(dust_tbl$id)) # back to orig. IDs

# seasonality
dust_tbl$season <- rep("", dim(dust_tbl)[1])
dust_tbl$season <- ifelse((dust_tbl$month >= 3 & dust_tbl$month <=5), "spring", dust_tbl$season)
dust_tbl$season <- ifelse((dust_tbl$month >= 6 & dust_tbl$month <=8), "summer", dust_tbl$season)
dust_tbl$season <- ifelse((dust_tbl$month >= 9 & dust_tbl$month <=11), "fall", dust_tbl$season)
dust_tbl$season <- ifelse(((dust_tbl$month >= 1 & dust_tbl$month <=2) | (dust_tbl$month == 12)), "winter", dust_tbl$season)
dust_tbl$season <- factor(dust_tbl$season, levels=c("spring", "summer", "fall", "winter"))

save(dust_tbl, file="data/dust.RData")

dustStn_tbl <- dustStn_tbl %>% 
  transmute(id=지점, start=시작일, end=종료일, name=지점명, lon=경도,
         lat=위도, elevation=`노장해발고도(m)`)
dustStn_tbl$id <- as.integer(ceiling(dustStn_tbl$id)) # back to orig. IDs
dustStn_tbl <- dustStn_tbl %>% arrange(id, end) # Order by dates
save(dustStn_tbl, file="data/Stations_dust.RData")

#####

# 2. Processing covariate dataset (ASOS + AWS)

# First, the stations
asos_stn <- read_csv("data/ASOS_stn_info.csv")
asos_stn <- asos_stn %>% 
  transmute(id=지점, start=시작일, end=종료일, name=지점명,
            lon=경도, lat=위도, elevation=`노장해발고도(m)`) %>% 
  arrange(id, end)
save(asos_stn, file="data/Stations_ASOS.RData")
aws_stn <- read_csv("data/AWS_stn_info.csv")
aws_stn <- aws_stn %>% 
  transmute(id=지점, start=시작일, end=종료일, name=지점명,
            lon=경도, lat=위도, elevation=`노장해발고도(m)`) %>% 
  arrange(id, end)
save(aws_stn, file="data/Stations_AWS.RData")

# Location changes (=> see file: 003-covariate.R)

#####

# Missing values
## THERE ARE A LOT (roughly ~36% from the total POSSIBLE data)
## We start from 21:00, March 21st (since that is the earliest we can go)
st_grid <- expand.grid(
  dt = seq(as.POSIXct("2003-03-21 21:00:00", tz="UTC"),
           as.POSIXct("2018-12-31 23:00:00", tz="UTC"), by=3600),
  id = unique(dust_tbl$id)
)
obs_tbl <- dust_tbl %>% select(id, dt, dust)
total_tbl <- left_join(st_grid, obs_tbl, by=c("id", "dt"))
sum(is.na(total_tbl$dust)) ## ~35.5%
all_miss_dt <- total_tbl %>% group_by(dt) %>% summarise(all_miss=all(is.na(dust))) %>% filter(all_miss)
all_miss_year <- all_miss_dt %>% mutate(year = year(dt))
g7 <- ggplot(all_miss_year, aes(year)) +
  geom_bar(fill="orange") + 
  labs(x="year", title="Years with all missing hrs.") +
  scale_x_continuous(breaks=seq(2003, 2018, 1))
ggsave("figures/eda/all_missing.png", g7, device="png")


#####

# Graphics

# 1. Seeing the distribution (by year)
g1 <- ggplot(data=dust_tbl, aes(x=log_dust)) +
  geom_histogram(aes(y = ..density..), fill="orange") +
  facet_wrap(~year) +
  labs(x="Hourly avg. dust (log)", title="Distribution by year") + 
  theme_bw()
ggsave("fig1-yearly_density.png", g1, device="png")


# 2. Geographic map of stations (more detailed map of Korea needed)
korea <- map_data("world", "South Korea")
g4 <- ggplot() +
  geom_polygon(data=korea, aes(long, lat, group=group), fill="gray", col="black") +
  geom_point(data=dustStn_tbl, aes(lon, lat, fill=elevation), alpha=.8, shape=21, colour="black", size=2.75) +
  geom_text(data=dustStn_tbl, aes(lon, lat, label=id), nudge_y=0.125, size=2.75, fontface="bold") +
  labs(title="PM10 station locations", fill="Elevation (m)") + 
  scale_fill_distiller(palette="Reds", direction=1) + theme_dark()
ggsave("figures/fig5-map.png", g4, device="png")

# 3. Summary of location/elevation changes
dupl_ids <- c(102,108, 115, 116, 152, 232)
dust_dupl_stn <- dustStn_tbl %>% arrange(id, end) %>%  # MUST be ordered (arrow direction must match)
  filter(id %in% dupl_ids) %>% select(id, lon, lat, elevation)
# "Normalizing"
for (i in seq(2,12,by=2)) dust_dupl_stn[i,2:4] <- dust_dupl_stn[i,2:4] - dust_dupl_stn[(i-1),2:4]
for (i in seq(1,11,by=2)) dust_dupl_stn[i,2:4] <- 0
g5 <- ggplot(dust_dupl_stn, aes(lon, lat)) +
  geom_path(aes(group=id), arrow=arrow(length=unit(.05, "inches"), type="closed")) +
  geom_text(aes(label=id), nudge_y=.00125, size=2.75) +
  geom_vline(xintercept=0, lwd=.2) + 
  labs(x="lon diff. (deg)", y="lat diff. (deg)") + theme_bw()

g6 <- ggplot(dust_dupl_stn, aes(as.factor(id), elevation)) +
  geom_path(aes(group=id), arrow=arrow(length=unit(.05, "inches"), type="closed")) +
  geom_hline(yintercept=0, lwd=.2) + 
  labs(x="ID", y="elevation diff. (m)") + theme_bw()
ggsave("figures/eda/fig12-dust_stn-loc_change.png", g5, device="png")
ggsave("figures/eda/fig13-dust_stn-elev_change.png", g6, device="png")

# Same for ASOS/AWS stations
g1 <- ggplot() +
  geom_polygon(data=korea, aes(long, lat, group=group), fill="gray", col="black") +
  geom_point(data=asos_stn, aes(lon, lat, fill=elevation), size=2, pch=21, alpha=.8, col="black") +
  scale_fill_distiller(palette="Blues", direction=1) +
  labs(fill="Elevation (m)") + theme_dark()
g2 <- ggplot() +
  geom_polygon(data=korea, aes(long, lat, group=group), fill="gray", col="black") +
  geom_point(data=(dustStn_tbl %>% filter(!(id %in% c(116, 94, 229, 501, 132, 144, 160)))), 
             aes(lon, lat, fill="ASOS present"), pch=21, colour="black", size=2.25) +
  geom_point(data=(dustStn_tbl %>% filter(id %in% c(116, 94, 229, 501))), 
             aes(lon, lat, fill="AWS present"), pch=21, colour="black", size=2.25) +
  geom_point(data=(dustStn_tbl %>% filter(id %in% c(132, 144, 160))), 
             aes(lon, lat, fill="None present"), pch=21, colour="black", size=2.25) +
  geom_text(data=dustStn_tbl, aes(lon, lat, label=id), 
            nudge_y=0.15, size=2.75, fontface="bold") +
  scale_fill_manual(name="Status", values=c("ASOS present"="blue", "AWS present"="green", "None present"="red")) +
  labs(title="Classifying PM10 stations") + theme_dark()
ggsave("fig4-asos_stn.png", g1, device="png", path="figures/eda/")
ggsave("figures/fig7-pm10_stn_class.png", g2, device="png")

## ASOS station location/elevation changes
# asos_dust_ids <- asos_2017$id %>% unique() # IDs applicable to PM10 locations
# asos_dust_stn <- asos_stn %>% filter(id %in% asos_dust_ids) %>% arrange(id, end)
# dupl_ids <- c(100, 102, 108, 136, 140, 143, 146, 152, 192, 232) 
# ## IDS with location changes
# asos_dupl_stn <- asos_dust_stn %>% filter(id %in% dupl_ids) %>% select(id, lon, lat, elevation)
# # "Normalizing"
# for (i in c(seq(2, 14, by=2), c(19,21))) asos_dupl_stn[i,2:4] <- asos_dupl_stn[i,2:4] - asos_dupl_stn[(i-1),2:4]
# for (i in c(seq(1, 13, by=2), c(18,20))) asos_dupl_stn[i,2:4] <- 0
# asos_dupl_stn[16:17,2:4] <- asos_dupl_stn[16:17,2:4] - asos_dupl_stn[c(15,15),2:4]
# asos_dupl_stn[15,2:4] <- 0
# 
# g3 <- ggplot(asos_dupl_stn, aes(lon, lat)) +
#   geom_path(aes(group=id), arrow=arrow(length=unit(.05, "inches"), type="closed")) +
#   geom_text(aes(label=id), nudge_y=.0025, size=2.75) +
#   geom_vline(xintercept=0, lwd=.2) + labs(x="lon diff. (deg)", y="lat diff. (deg)") + theme_bw()
# g4 <- ggplot(asos_dupl_stn, aes(as.factor(id), elevation)) +
#   geom_path(aes(group=id), arrow=arrow(length=unit(.05, "inches"), type="closed")) +
#   geom_hline(yintercept=0, lwd=.2) + labs(x="ID", y="elevation diff. (m)") + theme_bw()
# ggsave("figures/eda/fig14-asos_stn-loc_change.png", g3, device="png")
# ggsave("figures/eda/fig15-asos_stn-elev_change.png", g4, device="png")
