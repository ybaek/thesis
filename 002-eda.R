
library(dplyr)
library(ggplot2)
library(GGally)


# normalizing euclidean distance

load("dust.RData")

########

# episodes of severe PM10 levels 
dust_severe <- dust_tbl %>% filter(dust >= 800)
g2 <- ggplot(dust_severe, aes(year, log_dust)) +
  geom_col() + ylab("count") + scale_x_continuous(breaks=seq(2004, 2016, 1)) + theme_bw() +
  ggtitle(expression(paste("Severe levels (>800", mu, "g/", m^3, ")")))
g3 <- ggplot(dust_severe, aes(year, log_dust)) +
  geom_col(aes(fill=as.factor(season)), position="fill") +
  scale_x_continuous(breaks=seq(2004, 2016, 1)) +
  labs(y="", fill="season", title="Contribution per each season") +
  scale_fill_hue(h=c(0, 250)) + theme_bw()

gm <- ggmatrix(list(g2, g3), 2, 1, "Year", 
         title=expression(paste("Counts of severe episodes (>800", mu, "g/m"^3, ")")),
         legend=2)

ggsave("figures/eda/fig9-severe_per_seasons.png", gm, device="png")

# Variability across seasons
g4 <- ggplot((dust_tbl %>% filter(year>2012, year<2018, dust!=0.001)), aes(season, log_dust)) + 
  geom_boxplot(aes(fill=season), outlier.alpha=.1) +
  facet_wrap(~year, nrow=1, ncol=5) + labs(x="") +
  scale_fill_hue(h=c(0, 250)) +
  scale_x_discrete(breaks=NULL) +
  theme(legend.margin=margin(0,0,0,0,"pt"))

ggsave(g4, file="figures/eda/fig11-season_boxplot.png", device="png", dpi=120)

# Which stations did not have at least one occurrence of high level dust?
no_severe_ids <- dust_tbl %>% group_by(id) %>% summarise() %>% 
  anti_join((dust_severe %>% group_by(id) %>% count()), by="id")
no_severe_tbl <- dust_tbl %>% filter(id %in% no_severe_ids$id) %>% 
  group_by(id, lon, lat) %>% summarise()

korea <- map_data("world", "South Korea")
no_severe_tbl$id <- round(no_severe_tbl$id)
g5 <- ggplot() +
  geom_polygon(data=korea, aes(long, lat, group=group), fill="gray", col="black") +
  geom_point(data=no_severe_tbl, aes(lon, lat), fill="red", colour="black", pch=21, size=3) +
  geom_label(data=no_severe_tbl, aes(lon, lat, label=id), nudge_y=0.25, size=3) +
  labs(title="Stations with low PM levels") + theme_dark()
ggsave(g5, file="figures/fig8-no_severe_pm.png", device="png")



#########

# Further examination of seasonality (using spacetime)
library(spacetime)
library(SpatioTemporal)

dust_sub <- dust_tbl %>% 
  select(id, dt, log_dust, lat, lon, elevation)
obs <- dust_sub %>% 
  transmute(obs = log_dust, 
            date = dt,
            ID = as.character(id))
# For package reasons, time periods must be changed to integers
time_lookup <- select(obs, date) %>% unique() %>% arrange(date)
time_lookup$index <- seq(1, dim(time_lookup)[1])
obs <- inner_join(obs, time_lookup, by="date")
obs <- obs[,-2]
obs <- obs %>% rename(date = index)
###
covars <- dust_sub %>% 
  select(id, lat, lon) %>% unique() %>% rename(ID = id)
STdata <- createSTdata(obs, covars)


# ACF plots reveal long memory + slight yearly seasonality (cusps at 24 and 48-lags)
## Does this have to do with different locations/elevations of different stations?
plot_acf_fn <- function(id) {
  png(paste0("figures/acf-for-station", id, ".png"))
  plot(STdata, "acf", id)
  dev.off()
}
plot_pacf_fn <- function(id) {
  png(paste0("figures/pacf-for-station", id, ".png"))
  plot(STdata, "pacf", id)
  dev.off()
}
# Filter(plot_acf_fn, STdata$covars$ID)
# Filter(plot_pacf_fn, STdata$covars$ID)

# REMARKS: see note.


