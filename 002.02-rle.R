
#
# Run-length beyond 2 hrs of high levels
# (criterion: raw dust of 385 ("3 sigma rule"))
#


library(data.table)
library(dplyr)
library(lubridate)
library(ggplot2)

load("data/modelData.RData")
# 1. Severe counts
stns_years_severe <- data_grid %>% mutate(year = lubridate::year(dt)) %>% 
  group_by(id, year) %>% summarise(severe_count = sum(log_dust>log(800), na.rm=T))
g1 <- ggplot(stns_years_severe) + 
  geom_tile(aes(as.factor(id), as.factor(year), fill=severe_count), colour="grey50") + 
  labs(x="Station", y="Year", fill="count",
       title=expression(paste("Severe Episodes (>800", mu, "g m"^{-3}, ")"))) +
  scale_x_discrete(breaks=NULL) +
  scale_fill_distiller(palette="YlGnBu", direction=1) +
  theme(axis.text.x = element_text(angle=45, vjust=.7, colour="black", size=7.5),
        axis.text.y = element_text(colour="black", size=10),
        panel.background = element_blank(),
        legend.margin=margin(0,0,0,0,"pt")) + coord_flip()
ggsave("fig27-severe.png", path="figures/eda/", g1, device="png",
       width=8, height=6, units="in")

# 2. Prolonged period (>= 2hrs)
IDs <- unique(data_grid$id)
l <- list()
for (I in IDs) {
  dataI <- filter(data_grid, id==I) %>% mutate(severe = log_dust > log(385)) %>% data.table()
  dataI[,rle_severe := rleid(severe)]
  dataI[severe==TRUE, length_severe := nrow(.SD), by = "rle_severe"]
  dataI[length_severe>=2, check := TRUE]
  dataI[,year := lubridate::year(dt)]
  l[[as.character(I)]] <- dataI[,c("id", "year","check")]
}

years_rle <- Reduce(rbind, l) %>% group_by(year) %>% summarise(sum = sum(check, na.rm=T))
stns_rle <- lapply(l, function(df) sum(df$check, na.rm=T)) %>% as.data.frame() %>% 
  tidyr::gather("station", "count")
stns_years_rle <- lapply(l, function(df) (df %>% group_by(id, year) %>% summarise(sum = sum(check, na.rm=T)))) %>% 
  Reduce(rbind, .)

g2 <- ggplot(stns_years_rle) + geom_tile(aes(as.factor(id), as.factor(year), fill=sum), colour="grey50") + 
  labs(x="Station", y="Year", fill="count",
       title=expression(paste("Prolonged Episodes (>3", sigma, ", >=2 hrs.)"))) +
  scale_fill_distiller(palette="YlGnBu", direction=1) +
  scale_x_discrete(breaks=NULL) +
  theme(axis.text.x = element_text(angle=45, vjust=.7, colour="black", size=7.5),
        axis.text.y = element_text(colour="black", size=10),
        panel.background = element_blank(),
        legend.margin=margin(0,0,0,0,"pt")) + coord_flip()
ggsave("fig28-prolonged.png", path="figures/eda/", g2, device="png",
       width=8, height=6, units="in")
