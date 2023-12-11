
#------------------------------------------------------------------
# This code  
# 1) imports Kopaitic Island GPS data
# 2) speed filter GPS data (remove unreliable locations) with trip::sda
# 3) crawlWrap locations of multiple animals with a function (momentuHMM)
# to give locs every 5 minutes

# Leandri de Kock, Chris Oosthuizen
# April 2022
# September 2022 

#------------------------------------------------------------------

# load packages
library(tidyverse)
library(momentuHMM)
library(sf)

Sys.setenv(TZ = "GMT")

# -------------------------------------------------------------------------
#### read in the GPS data for a group of penguins
# -------------------------------------------------------------------------
# raw GPS data from tag has already been processed to trips (locations at the nest removed)
gpsdat = readRDS("./outputs/data/combined data/CS_Kop_all_files_combined.rds")

# make sure date.time is set to GMT (or UTC)
attr(gpsdat$date.time, "tzone") # check time zone
gpsdat$date.time = lubridate::with_tz(gpsdat$date.time, "GMT")  
attr(gpsdat$date.time, "tzone") # check time zone

#-----aside---------------------------------------------------------
#----------------------------------------------------------------
## Set deployment latitudes (home.lat) (per trip, not per id)
#----------------------------------------------------------------
# add GPS locations of breeding sites, although I do not use this in the end.

# unique(gpsdat$site)
# 
# PaddyRocks.lat = -62.3109	
# FarNorth.lat = -62.2903
# WindyRidge.lat = -62.3058	
# Sunset.lat = -62.3037	
# 
# PaddyRocks.lon = -59.2115
# FarNorth.lon = 	-59.2273
# WindyRidge.lon = -59.2325
# Sunset.lon = -59.2321
#   
# gpsdat = gpsdat %>%
#   mutate(
#     home.lat = case_when(
#        site == "PaddyRocks"  ~ PaddyRocks.lat,
#         site == "FarNorth" ~ FarNorth.lat,
#         site == "WindyRidge" ~ WindyRidge.lat,
#        site == "Sunset"  ~ Sunset.lat))
# 
# gpsdat = gpsdat %>%
#   mutate(
#     home.lon = case_when(
#        site == "PaddyRocks"  ~ PaddyRocks.lon,
#         site == "FarNorth" ~ FarNorth.lon,
#         site == "WindyRidge" ~ WindyRidge.lon,
#        site == "Sunset"  ~ Sunset.lon))
#-----aside---------------------------------------------------------


# make a summary of track / trip start and end points 

trip.summ = gpsdat %>%
   group_by(id, trip_id) %>%
   summarise(track.start.date = min(track.start.date),   # first trip starts
             track.end.date = max(track.end.date),       # last trip ends
             trip.start.date = min(trip.start.date),  
             trip.end.date = max(trip.end.date),
             track.duration.hr = min(track.duration.hr),   # delta track start & end
             trip.duration.hr = max(trip.duration.hr),
             start = min(start),   # when GPS was switched on
             end = max(end)    # when GPS was switched off
             )

trip.summ 
#write.csv(trip.summ, "summ.csv")
# trip.sum still include trip_id == 0, which is when the penguin is at the nest. Remove those

trip.sum = trip.summ %>% 
  dplyr::filter(trip_id != 0) %>%   # remove locations not on trips
  dplyr::select(id,	trip_id, trip.start.date,	trip.end.date)
trip.sum$crawl_id = paste0(trip.sum$id, "_", trip.sum$trip_id)

trip.sum

#----aside---------------------------------------------------------
# check = gpsdat %>%
#    group_by(id, trip_id) %>%
#    summarise(track.start.date = n_distinct(track.start.date),
#              track.end.date = n_distinct(track.end.date),
#              trip.start.date = n_distinct(trip.start.date),
#              trip.end.date = n_distinct(trip.end.date),
#              track.duration.hr = n_distinct(track.duration.hr),
#              trip.duration.hr = n_distinct(trip.duration.hr),
#              start = n_distinct(start),
#              end = n_distinct(end)
#              )
#----aside---------------------------------------------------------

# Process data
gpsdat = gpsdat %>%
  as_tibble(gpsdat) %>%
  dplyr::select(date.time, lat, lon, tracks, id, trip_id)  %>%   # remove unnecessary columns
  dplyr::filter(trip_id != 0)  %>%   # remove points at nests
  dplyr::group_by(id) %>%
  arrange(date.time) %>% # put the data in date order
  distinct(date.time, .keep_all=TRUE) %>%   # remove duplicate date.time records for individuals
  dplyr::ungroup()

gpsdat

unique(gpsdat$tracks)

gpsdat = gpsdat %>% mutate(stage =
                     case_when(tracks == 'Inc1'~ "Incubation", 
                               tracks == 'Inc2' ~ "Incubation",
                               tracks == 'Bro1' ~ "Brood",
                               tracks == 'Bro2' ~ "Brood",
                               tracks == 'Bro3' ~ "Brood",
                               tracks == 'Bro4' ~ "Brood",
                               tracks == 'Cre1' ~ "Creche",
                               tracks == 'Cre2' ~ "Creche")
)


gpsdat
unique(gpsdat$id)
length(unique(gpsdat$id))

gpsdat$crawl_id = paste0(gpsdat$id, "_", gpsdat$trip_id)
n_distinct(gpsdat$crawl_id)
unique(gpsdat$crawl_id)

# clean up
gpsdat = gpsdat %>% dplyr::select(-tracks)


#------------------------------------------------------------------------
# Add nest coordinates to data
#------------------------------------------------------------------------
# I want to add the lat lon coordinates of the island at the start and end points
# of each trip, i.e., when the penguin left the island and returned to it
# as determined from TDR data. This makes sure every penguin departs from and arrives
# back at the island

HP.lon = -57.9162
HP.lat = -63.3160

trip.sum2 = trip.sum %>%
              dplyr::select(id, crawl_id, trip.start.date, trip.end.date, trip_id)
trip.sum2$lat = HP.lat    # use as a proxy for all nests
trip.sum2$lon = HP.lon    # use as a proxy for all nests
trip.sum2$stage = "home"
trip.sum2

# Add the coordinates to the trip start points
trip.starts = trip.sum2 %>%
              dplyr::select(trip.start.date, lat, lon, id, trip_id, stage, crawl_id)%>%
              dplyr::rename(date.time = trip.start.date)
trip.starts

# Add the coordinates to the trip end points
trip.ends = trip.sum2 %>%
              dplyr::select(trip.end.date, lat, lon, id, trip_id, stage, crawl_id)%>%
              dplyr::rename(date.time = trip.end.date)
trip.ends

# add the start and end GPS points and times to the main data
dim(gpsdat)
gpsdat = rbind(gpsdat, trip.starts)
#gpsdat = rbind(gpsdat, trip.ends)    # don't add end point as this results in unnatural tracks for those tags that ran out of battery
dim(gpsdat)

#-------------------------------------------------------------------------------
####  Speed filter
#-------------------------------------------------------------------------------
# quick plot to see whether there are outliers

# first define the number of colors you want
library(RColorBrewer)
nb.cols <- length(unique(gpsdat$id))
mycolors <- colorRampPalette(brewer.pal(8, "Set2"))(nb.cols)

ggplot() + 
  geom_point(data = gpsdat, size = 0.25, alpha = 1, aes(x = lon, y = lat, color = id)) +
  scale_x_continuous(expand = expansion(mult = c(.1, .1))) +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.1))) +
  scale_fill_manual(values = mycolors) +
  theme_bw() +
  theme(legend.position="none")


# use trip::sda to filter data based on speed, distance, angle
# define coordinates and CRS 
sp::coordinates(gpsdat) <- ~lon+lat
sp::proj4string(gpsdat) <- sp::CRS("+init=epsg:4326", doCheckCRSArgs = FALSE)

library(trip)
# create object "trip"
gpsdat.tr <- trip::trip(gpsdat, c("date.time", "id"))
summary(gpsdat.tr)

# speed filter
gpsdat.tr$speedfilter <- trip::sda(gpsdat.tr, 
                                   smax = 15,       # maximum speed, in km/h
                                   ang = c(15, 25),    # Freitas et al 2008's minimum turning angles in degrees
                                   distlim = c(2.5, 5)) 	# Freitas et al 2008's maximum step lengths in km

table(gpsdat.tr$speedfilter) # FALSE points are being removed

# remove any false locations 
gpsdat = subset(gpsdat.tr, gpsdat.tr$speedfilter == "TRUE")  
gpsdat = as_tibble(gpsdat)  # change back from trip object to tibble

# quick plot to see new data with no outliers
ggplot() + 
  geom_point(data = gpsdat, size = 0.25, alpha = 1, aes(x = lon, y = lat, color = id)) +
 # facet_grid(~ stage) + 
  scale_x_continuous(expand = expansion(mult = c(.1, .1))) +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.1))) +
  scale_fill_manual(values = mycolors) +
  theme_bw()+
  theme(legend.position="none")

unique(gpsdat$trip_id)
n_distinct(gpsdat$id)  # should be 62 for Kopaitic
n_distinct(gpsdat$crawl_id) # should be 226 for Kopaitic

#Look at sample size (GPS uplinks) and tracking time per indiv
samplesize = gpsdat %>% 
             count(crawl_id)
samplesize

min(samplesize$n)     
hist(samplesize$n)

gpsdat = merge(gpsdat, samplesize, 
                 by.x=c("crawl_id", "crawl_id"))

# now exclude birds with few data points
# If it has 2 data points it is only the 2 points at home - don't have any GPS data from the trip (only TDR data)
# Now it includes at least x GPS points at sea
gpsdat = gpsdat %>%
         dplyr::filter(n > 4)             # keeping 2 home locations plus 3 at sea reduces the data to 224 trips
       
# Some tracks can be too short in time too, not only in GPS points:
tracktime = gpsdat %>% 
  group_by(crawl_id) %>% 
  mutate(timediff =  as.numeric(difftime(last(date.time), first(date.time)), units="mins")) %>%
  slice(1) %>% 
  dplyr::select(crawl_id, timediff)
min(tracktime$timediff)

# Make the minimum track time 15 minutes:
gpsdat = merge(gpsdat, tracktime, 
               by.x=c("crawl_id", "crawl_id"))

gpsdat = gpsdat %>%
  dplyr::filter(timediff > 15)             # reduces the data to 222 trips

min(gpsdat$timediff)

# --------------------------------------------------------
# run crawl model to predict locations at 5 min intervals
# --------------------------------------------------------
crawldat = gpsdat %>%
  dplyr::select(crawl_id, date.time, lat, lon)  %>%
  dplyr::rename(ID = crawl_id) %>%  # crawl requires variable called ID
  sf::st_as_sf(coords=c("lon", "lat"), crs = 'EPSG:4326') %>% # gpsdata (lat lon) is in wgs84
  sf::st_transform('EPSG:32721') # transform lat lon data to UTM zone 21S https://epsg.io/32721 (crawl requires UTM)

print(crawldat)

# specify error around GPS points: assume a 50 m isotropic error ellipse for the measurement error model 
# source: September 2, 2021 momentuHMM vignette
lnError <- crawl::argosDiag2Cov(50,50,0) 
crawldat$ln.sd.x = lnError$ln.sd.x
crawldat$ln.sd.y = lnError$ln.sd.y
crawldat$error.corr = lnError$error.corr

# Run Crawl models to predict REGULAR track locations (at-sea to include foraging, resting and swimming)
crw2min = crawlWrap(crawldat, 
                    theta =c(5.5,-0.2), 
                    Time.name = "date.time",
                    timeStep = "5 min",
                    fixPar = c(1,1,NA,NA),
                    err.model = list(x = ~ln.sd.x-1,
                                     y = ~ln.sd.y-1,
                                    rho = ~error.corr),
                    attempts = 100)

# save model output
crw.fit = crw2min$crwFits

# get predictions from crawl and transform from UTM back to lat lon
crw.pred = as_tibble(prepData(crw2min)) %>%
  dplyr::select(ID, date.time, x, y) %>%
  dplyr::rename(id = ID, lon = x, lat = y) %>%
  sf::st_as_sf(coords = c("lon", "lat"), crs = 'EPSG:32721') %>%
  sf::st_transform(crs = 'EPSG:4326')

# plot the predictions
crw.preddat = as_tibble(prepData(crw2min))

ggplot(data = crw.preddat, aes(x = x, y = y)) +
  geom_point(size = 1, alpha = 1) +
  theme_bw() 


# save model outputs
crw.preddat$round = substr(as.character(crw.preddat$ID) , 1, 13)
unique(crw.preddat$round)

crw.preddat$stage = substr(as.character(crw.preddat$ID) , 1, 12)
unique(crw.preddat$stage)
crw.preddat = crw.preddat %>% 
                   mutate(stage =
                             case_when(stage == 'Kopaitic_Bro'~ 'Kopaitic_Bro',
                                       stage == 'Kopaitic_Cre'~ 'Kopaitic_Cre',
                                       stage == 'Kopaitic_Inc'~ 'Kopaitic_Inc'))

unique(crw.preddat$stage)

saveRDS(crw.fit, paste0("./outputs/crawl/5 min crawl_model_fit_","CS_Kopaitic", '.rds'))
saveRDS(crw.preddat, paste0("./outputs/crawl/5 min crawl_model_preddat_", "CS_Kopaitic", '.rds'))
