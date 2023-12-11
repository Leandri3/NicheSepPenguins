#----------------------------------------------------------------------------------------
# This code: 
# Reads in divestats and cluster data for dives 5m and deeper
# Filters for only foraging dives
# Extracts the dive characteristics that we want to use for each dive
# Summarizes dive characteristics to a .rds file
#
# Leandri de Kock
# April 2022
#-------------------------------------------------------------------------------------- 

#load packages:
library(tidyverse)

#--------------------------------------------------------------------------------------
# 1. Kopaitic chinstraps
#--------------------------------------------------------------------------------------

# Load divestats file with clusters and dive residuals included:
clus_divestats <- readRDS('./outputs/clusters/divestats_5m_xy_clusters_chinstrap_kop.rds')

names(clus_divestats)

# filter for foraging dives
clus_divestats = clus_divestats %>%
  dplyr::filter(dive_5m_cluster == 1)

# Check that there is the same amount of dives compared to the original divestats file
#clus_divestats = clus_divestats %>%
#  dplyr::filter(ID == "Kopaitic_Bro1_KI01")
# same number of dives as original divestats 
# Advantage: Includes all individuals from species x island for all breeding phases. 

# Add new variableS: "group" for breeding phase; "island" to specify island
clus_divestats$group = substr(clus_divestats$ID, 10, nchar(clus_divestats$ID)-6)
unique(clus_divestats$group)

clus_divestats$island = substr(clus_divestats$ID, 1, nchar(clus_divestats$ID)-10)
unique(clus_divestats$island)

# add in breeding phase: Incubation or chick-rearing: 
# Change Bro and Cre to chick-rearing phase:
#--------------------------------------------------------------
# rename 
clus_divestats$phase = if_else(clus_divestats$group  == "Bro", "Chick",clus_divestats$group)
clus_divestats$phase = if_else(clus_divestats$phase == "Cre", "Chick", clus_divestats$phase)

unique(clus_divestats$phase)
unique(clus_divestats$group)
#--------------------------------------------------------------------

# Any other dive metrics that we can think about? 

# Make a summary for dive metrics:
dive_metrics <- clus_divestats %>%
  dplyr::select(c(begdesc, end.dive, botttim, divetim, bottdist, maxdep, postdive.dur, dives.in.trip, prop.deep.dives,
           prop.botttime, total.bouts, ID, species, x, y, dive.res, dive_5m_cluster, group, island, phase))
names(dive_metrics)

# Other dive metrics to include: 
#------------------------------------------------------------------------------------------------
# 1. number of foraging dives in a track (per individual)  -- we have this by only using foraging dives
# How many foraging dives did this penguin do in the track? 
#foraging.dives = dive_metrics %>%
#                group_by(ID) %>%
#                filter(divecluster == 1) %>%  # cluster 1 = foraging
#                tally() # how many records are there?
#              
#foraging.dives   
#
## Add number of foraging dives per individual track to dive_metrics
#dive_metrics_for = left_join(dive_metrics, foraging.dives, by = c("ID" = "ID"))
#dive_metrics_for
#
## rename n ( number of foraging dives)
#dive_metrics_for = dive_metrics_for %>%
#  rename(for.dives = n)
#
## proportion of foraging dives:
#dive_metrics_for$prop.for.dives = dive_metrics_for$for.dives / dive_metrics_for$dives.in.trip
#dive_metrics_for$prop.for.dives
#
#summary(dive_metrics_for$prop.for.dives)
## mean  = 0.3620 = 36.2 % of dives are foraging dives

#-----------------------------------------------------------------------------------------------
# 2. Attempts of catch per unit effort (ACPUEt): 
# number of possible wiggles (bottdist)/ total time in bottom duration
# can be by dive level (every dive) or by track level (every track), i.e. every individual on its own
# Total number of possible prey capture attempts (with or without success) relative to the total time spent in the
# bottom phase per trip

# By dive level: Proportion of possible prey capture attempts performed for the total time spent 
# at the bottom phase of the dive
dive_metrics$acpue = dive_metrics$bottdist/dive_metrics$botttim
dive_metrics$acpue  # *100 for percentage
summary(dive_metrics$acpue)

#-----------------------------------------------------------------------------------------------
# 3. Calculate solar position (elevation) and determine time of day for dives using maptools::solarpos
# Load packages
library(raster)
library(maptools)

dat_nel = dive_metrics %>%
  dplyr::select(begdesc, ID, x, y)

str(dat_nel)

# make a copy for the correct time zone
# remove rows with NAs in column
dat.gmt = dat_nel[complete.cases(dat_nel$y),]

# make sure date.time is set to GMT (or UTC)
attr(dat.gmt$begdesc, "tzone") # check time zone

# # projections
utm.prj = " +proj=utm +zone=21 +south +datum=WGS84 +units=m +no_defs "   # Chris UTM King George
wgs84 <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"
#
# # create a spatial points object:
# # remove rows with NAs in column
# dat.sp = dat_nel[complete.cases(dat_nel$y),]
#
# #To assign a known CRS to spatial data:
 utm.coord = SpatialPoints(cbind(dat.gmt$x, dat.gmt$y), proj4string=CRS(utm.prj))
 utm.coord
#
# # To transform from one CRS to another:
 wgs.coord <- spTransform(utm.coord, CRS(wgs84))
 wgs.coord
#
 dat.gmt$lon2 <- wgs.coord$coords.x1
 dat.gmt$lat2 <- wgs.coord$coords.x2
 head(dat.gmt)
#
# # First change time zone
# # make sure date.time is set to GMT (or UTC)
# attr(dat.sp$begdesc, "tzone") # check time zone # gmt
#
# # Make a new variable with antarctic peninsula time zone
# # OlsonNames() #Returns  a list of valid time zone names.
# dat.sp$begdesc.ap.tz = lubridate::with_tz(dat.sp$begdesc, "America/Punta_Arenas")   #antarctic peninsula: UTC -3
# attr(dat.sp$begdesc.ap.tz, "tzone") # check time zone
#
# dat.sp$begdesc.ap.tz
#
# # make a second spatial object  to change timezone of points
# dat.sp2 = dat.sp %>%
#   dplyr::select(begdesc.ap.tz, lon2, lat2)
#
# str(dat.sp2)
#
# dat.sp2 = SpatialPoints(cbind(dat.sp2$lon2, dat.sp2$lat2), proj4string=CRS(wgs84))
#
# # Change from one time zone to another? # can use same time zone for both islands
# d = as.POSIXct(dat.sp$begdesc.ap.tz, tz = "America/Punta_Arenas")  # antarctic peninsula: UTC-3
#
# # solarpos returns a matrix with the solar azimuth (in degrees from North), and elevation (in degrees measured vertically
# # from that point on the horizon up to the object) .

dat.gmt$solar.ele.begdesc <- maptools::solarpos(as.matrix(cbind(dat.gmt$lon2, dat.gmt$lat2)),
                                                proj4string=CRS(wgs84),
                                                direction = 'sunset', dat.gmt$begdesc, POSIXct.out = T)[,2]


# copy antarcitc peninsula time zone to original gmt dataset
# dat.gmt$begdesc.ap_tz <- dat.sp$begdesc.ap.tz

# add in the solar elevation to original gmt dataset
# dat.gmt$solar.ele.begdesc <- s[,2]
# dat.gmt

plot(dat.gmt$begdesc, dat.gmt$solar.ele.begdesc)

# classify solar elevation as day, night or twilight
# for us: how low does the sun dip below the horizon?
summary(dat.gmt$solar.ele.begdesc)

# create a day, night, twilight variable
dat.gmt$solar.ele.begdesc.day = dat.gmt$solar.ele.begdesc
dat.gmt$solar.ele.begdesc.day <- replace(dat.gmt$solar.ele.begdesc.day, dat.gmt$solar.ele.begdesc.day > 0, 'day')
dat.gmt$solar.ele.begdesc.day <- replace(dat.gmt$solar.ele.begdesc.day, dat.gmt$solar.ele.begdesc.day < '0', 'twilight')
dat.gmt$solar.ele.begdesc.day

names(dat.gmt)
table(dat.gmt$solar.ele.begdesc.day)

--------------------------------------------------------

dive_metrics_day <- left_join(dive_metrics, dat.gmt, by=c("begdesc" = "begdesc", "ID" = "ID", 
                                                                  "x" = "x", "y" = "y"))
dive_metrics_day

# save the dive metric file
# saveRDS(dive_metrics_day, paste0("./outputs/summary/dive_5m_metrics_foraging_", dive_metrics_day$species[1], "_", dive_metrics_day$island[1], '.rds'))
