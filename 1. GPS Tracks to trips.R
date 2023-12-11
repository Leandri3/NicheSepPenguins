
#----------------------------------------------------------
# Leandri de Kock, Chris Oosthuizen
# Sept 2021

# Divide GPS track into separate trips to sea based on diving data
# Calculate maximum and cumulative trip distances and trip durations based on GPS

# Works for a single penguin at a time

#----------------------------------------------------------
setwd("C:/Users/cwoosthuizen/Desktop/krilltokt penguins")

Sys.setenv(TZ = "GMT")

library(tidyverse)
library(data.table)
library(lubridate)
library(chron)

# Import chinstrap penguin accellerometer data (rds file)
df = readRDS("./data/GPS/Kopaitic chinstrap/rds data/CS_Kop_R1_Inc1_KI01.rds")
df

# combine id with round to give unique ID's
df$id = paste(df$island, df$tracks, df$tag, sep= "_")

# save a copy
df.input = df

# Fill in missing values for lat and lon with previous value
df = df %>% fill(lat) %>%
            fill(lon) 

# Fill in missing values for lat and lon with later value (fill first set of NA's)
df = df %>% fill(lat, .direction = "up")  %>%
            fill(lon, .direction = "up")  

# only accelerometer values exist for every row - e.g. for depth, "NA" fills the gaps
# remove all those rows that have "NA" in them, leaving just depth values. 
dat  = df %>% tidyr::drop_na(depth)
dat 

# make depth negative
dat$depth = dat$depth * -1 

#-------------------------------------------------------------------------------------------------
# Select periods that the bird was at the nest (no diving) and remove those points
# This is done manually by setting the following:
# - date.time where you want to start remove data
# - istart (lag after date.time to start removing data)
# - iend - how many HOURS of data to remove?

# https://stackoverflow.com/questions/31856643/removing-data-between-specific-dates-in-r
#-------------------------------------------------------------------------------------------------

# set the date based on the input data. There can be any number of entries here:
# there will be 0 if the bird did only 1 trip
# there will be 3 if the bird did 4 separate trips

# the point is to create gaps in the data: not to remove all points at the nest
# that happens later

# plot
plot(dat$date.time, dat$depth, pch = ".")

# select the rough start date of 'land' periods here, and the number of hours the penguin 
# was at the nest. Start dates need not be very precise, but it should not include 
# diving data

# ============================================
# ============================================
# Multiple trips 
# ============================================
# ============================================

# You can use the locator() to help you identify the starting date to specify in the code:
#mouseclicks = locator()
 mouseclicks = as.data.frame(mouseclicks$x)
 colnames(mouseclicks) = "date"
 mydates = mouseclicks$date
 class(mydates) = c('POSIXt','POSIXct')
 as.data.frame((mydates))

 dd <- data.frame(date.time = c('2018-12-11 07:52:07'),
                   istart = c(0),  
                  iend = c(1))  #hours 
 
 dd$date.time = as.POSIXct(format(dd$date.time))
 dd
 dates <- paste(dd[,1])
 dates
 istart <- dd$date.time + dd[,2]
 iend <- dd$date.time + dd[,3]*60*60
 
# # draw polygon - the data inside the polygon will be removed (should not contain dives)
 for (i in 1:length(iend)){
   polygon(c(istart[i],iend[i],iend[i],istart[i]),c(0,0,-110,-110),
           col=rgb(1, 0, 0, 0.4), border= "red")
 }
 
# #-------------------------------------------
# # Now remove the data inside the polygons
# #-------------------------------------------
 
 ret <- rep(FALSE, NROW(dat))
 for (i in seq_along(istart)) {
     ret <- ret | ((dat$date.time >= istart[i]) & (dat$date.time <= iend[i]))
 }
 
 dat_split <- dat[!ret, ]

# # there should be gaps where the polygons were
plot(dat_split$date.time, dat_split$depth, pch = ".")
# 
 for (i in 1:length(iend)){
   polygon(c(istart[i],iend[i],iend[i],istart[i]),c(0,0,-110,-110),
           col=rgb(0, 1, 0,0.5), border= NA)
 }

rstudioapi::savePlotAsImage(paste0("./outputs/plots/chinstrap/Kopaitic GPS processing/", "CS_dive_", dat$id[1],".png"), width=900,height=600)

#----------------------------------------------------------------
# Further refine splitting track into trips:
#----------------------------------------------------------------
# Now work with the data where positions at the nest have been removed.
dat = dat_split   # It is OK to throw an error here if there was a single trip

# first calculate time difference between points:
# order
dat = dat[order(dat$date.time),]

# add Date and Time
dat$Date <- as.Date(dat$date.time) 
dat$Time <- format(as.POSIXct(dat$date.time), format = "%H:%M:%S") 

# convert to chron date 
dat$chron.date_time <- chron(as.character(dat$Date),dat$Time, format=c(dates="Y-m-d",times="h:m:s")) 

# calculate time diff in minutes
timediff <- diff(dat$chron.date_time)*24*60

# remove first entry without any difference
dat <- dat[-1,]

# assign timediff column
dat$timediff <- as.numeric(timediff)          

# absolute time differences
dat$timelag <-as.numeric(abs(timediff))

range(dat$timediff)
range(dat$timelag)

#What are the 10 largest gaps in timediff?
tail(sort(dat$timelag), 10 )
format(tail(sort(dat$timelag), 10 ), scientific = FALSE) 

# Assume there are X large gaps (i.e. X + 1 foraging trips)
# that you made above. 
# Then you want to edit the timelag value below, so that you 
# also select X number of the highest values. 

# add a column for trip start date
# when there is a timelag of more than X minutes, assign a 1, otherwise assign a 0
dat = dat %>% 
     mutate(tripstart = case_when(timelag > 50 ~ 1, 
                                          TRUE ~ 0))

sum(dat$tripstart)  # should be the number of trips, less 1 (the very first start date)

# now also make first point in data a 1 (first start date)
dat$tripstart[1] <- 1

# subset to only starting points of each trip
starts  = subset(dat, dat$tripstart == 1)

# Plot start points
# ggplot(data = dat, aes(date.time, depth, colour =  Temp)) +
#     geom_point(shape= ".") +
#     scale_colour_gradientn(colours = terrain.colors(10))+
#     geom_point(data = starts, aes(x = date.time, y = depth), 
#                colour = 'red')

# add a column for the end date of trips
dat = dat %>%  
      mutate(tripend = lead(tripstart))   # 'lead' takes the previous value (see lead and lag from dplyr)

dat$tripend[length(dat$tripend)] <- 1   # make the very last data point also a '1' for the end of that trip

ends = subset(dat, dat$tripend == 1)

# plot start and end points
# ggplot(data = dat, aes(date.time, depth, colour =  Temp)) +
#   geom_point(shape = ".") +
#   scale_colour_gradientn(colours = terrain.colors(10))+
#   geom_point(data = starts, aes(x = date.time, y = depth), colour = 'red')+
#   geom_point(data = ends, aes(x = date.time, y = depth), colour = 'blue')


# Assign trip numbers:
dat$trip_id <- cumsum(ifelse(dat$tripstart == 1, 1, 0))   
# 1 = value to look for
# 1 =  number that each PCE should 'advance' with if yes
# 0 =  number that each PCE should 'advance' with if no 

unique(dat$trip_id)

ggplot(data = dat, aes(date.time, depth, colour = as.factor(trip_id)))+
 geom_point(shape = ".") +
  #scale_colour_gradientn(colours = terrain.colors(10))+
  geom_point(data = starts, aes(x = date.time, y = depth), colour = 'red')+
  geom_point(data = ends, aes(x = date.time, y = depth), colour = 'blue') 


#------------------------------------------------------------
# cut out leading zeros of depth for each trip
#------------------------------------------------------------
# https://stackoverflow.com/questions/67950197/dplyr-find-first-non-zero-element-and-last-non-zero-element-and-trim-vector-by

# Each trip still has a period at the start where depth is 0 (penguin on land)
# also, each trip ends with many records of depth = 0 (penguin back on land)
# This code removes the leading and trailing zeros, by trip

# Because of measurement error, etc. there may be small non-zero movements 
# in the depth data even when the penguin is close to land.
# make a second depth variable (depth_zero) where all depths shallower than 
# 1 m is set to 0. This 'smooths' near-surface depth data
dat = dat %>% mutate(depth_zero = ifelse(depth > -1 , 0, depth))

# now use this to cut trips
dat.trips = dat %>%
        group_by(trip_id) %>% 
        filter(cumsum(depth_zero != 0) != 0) %>%                
        filter(rev(cumsum(rev(depth_zero != 0))) != 0 ) %>%      
     ungroup()

dim(dat)
dim(dat.trips)  # some leading/trailing zero's were removed

ggplot(data = dat.trips, aes(date.time, depth, colour = as.factor(trip_id)))+
    geom_point(shape = ".") +
    #scale_colour_gradientn(colours = terrain.colors(10))+
    geom_point(data = starts, aes(x = date.time, y = depth), colour = 'red')+
    geom_point(data = ends, aes(x = date.time, y = depth), colour = 'orange') 


# ---------------Aside------------------------------------
# look at each trip (zoom in to see if the code worked)
# trip 1
trip1 = subset(dat, dat$trip_id == 1)
ggplot(data = trip1, aes(date.time, depth, colour = as.factor(trip_id) ))+
geom_line() 

trip1.trips = subset(dat.trips, dat.trips$trip_id == 1)
ggplot(data = trip1.trips, aes(date.time, depth, colour = as.factor(trip_id) ))+
geom_line() 

# ---
ggplot(data = trip1[1:20000,], aes(date.time, depth, colour = as.factor(trip_id) ))+
geom_line() 

ggplot(data = trip1.trips[1:20000,], aes(date.time, depth, colour = as.factor(trip_id) ))+
geom_line() 

# works!
# ---------------Aside------------------------------------

# trip start and trip end was based on the old splits.
# you've now cleaned that based on depth.
# Need to re-define trip start and end dates with the new data:

dat.trips =  dat.trips %>%
         group_by(trip_id) %>%
         mutate(tripstart = if_else(duplicated(trip_id) == FALSE, 1, 0)) %>%
         ungroup()

dat.trips =  dat.trips %>%
         group_by(trip_id) %>%
         mutate(tripend = if_else(rev(duplicated(trip_id)) == FALSE, 1, 0)) %>%
         ungroup()


# redefine and plot the new trip and end dates:
starts_cleaned  = subset(dat.trips, dat.trips$tripstart == 1)
ends_cleaned = subset(dat.trips, dat.trips$tripend == 1)

ggplot(data = dat.trips, aes(date.time, depth, colour = as.factor(trip_id)))+
    geom_point(shape = ".") +
    #scale_colour_gradientn(colours = terrain.colors(10))+
    geom_point(data = starts_cleaned, aes(x = date.time, y = depth), colour = 'red')+
    geom_point(data = ends_cleaned, aes(x = date.time, y = depth), colour = 'purple') 

rstudioapi::savePlotAsImage(paste0("./outputs/plots/chinstrap/Kopaitic GPS processing/", "CS_trip_", dat$id[1],".png"), width=900,height=600)

# Should now look good.
unique(dat.trips$trip_id)
unique(dat.trips$id)

# renumber trip_id if there are gaps (should not be the case)
# dat.trips <- dat.trips %>%
# mutate(trip_id = as.numeric(factor(trip_id, levels = unique(trip_id)))) 

head(as.data.frame(dat.trips))
dim(dat.trips)


# ----------------------------------------------------------------------------------------
# Summary stats: start and end date of every animal's tracking period (over multiple trips)
# ----------------------------------------------------------------------------------------

# Mutate a column with track start dates 
dat.trips =  dat.trips %>%
     mutate(track.start.date = date.time[1]) %>%
     ungroup() 

# Select the last record for every tag id (1 row per animal)
temp =  dat.trips %>%
         slice(n())  %>%
         mutate(track.end.date = date.time) %>%
         dplyr::select(id, track.end.date) %>%
         ungroup() 

# Merge end dates with data 
dat.trips = left_join(dat.trips, temp, by = "id")   # merge: left_join keeps all rows of tripdata, and matches 'id' 

# ------------------------------------------------------------------------------
# Summary stats: start and end of each foraging trip (within tracks)
# ------------------------------------------------------------------------------

# Mutate a column with trip start dates 
dat.trips = dat.trips %>%
              group_by(trip_id) %>%
              mutate(trip.start.date =  date.time[1]) %>%
              ungroup() 
 
# Select the last record for every tag id (1 row per animal)
temp2 = dat.trips %>%
              group_by(trip_id) %>%
              slice(n())  %>%
              mutate(trip.end.date = date.time) %>%
              dplyr::select(trip_id, trip.end.date) %>%
              ungroup() 
 
# Merge end dates with data 
dat.trips = left_join(dat.trips, temp2, by = "trip_id")   # merge: left_join keeps all rows of tripdata, and matches 'id' 

# ------------------------------------------------------------------------------
# Summary stats: duration 
# ------------------------------------------------------------------------------
dat.trips = mutate(dat.trips, track.duration.hr = 
                     abs(difftime(track.end.date, track.start.date, units="hours")))

dat.trips = mutate(dat.trips, trip.duration.hr = 
                     abs(difftime(trip.end.date, trip.start.date, units="hours")))

head(as.data.frame(dat.trips))
dim(dat.trips)

tripdata.sum = dat.trips %>% 
               dplyr::distinct(trip_id, .keep_all = TRUE)

as.data.frame(tripdata.sum)
dim(tripdata.sum)

write.csv(tripdata.sum, "./outputs/summary/CS_Kop_tripdata_summary.csv")

# Append other summary data to this file using the following code:
# https://stackoverflow.com/questions/45930711/appending-a-new-line-into-an-existing-csv-file
# write.table(tripdata.sum,  
#              file = "./outputs/summary/cs_Kop_tripdata_summary.csv", 
#              append = T, 
#              sep = ',', 
#              row.names = T, 
#              col.names = F )

# ----------------------------------------------------------------------------
# Calculate distance moved per trip from GPS data
# This is 'probably' not for analysis, but it is useful to plot
# ----------------------------------------------------------------------------

library(ggspatial)
library(sp)  

# ---------------------------------------------------------------
# Kopaitic had a 1-sec sampling rate for the GPS.
# Subsample - make a column that summarizes time as every 30 seconds
# --------------------------------------------------------------- 

dat.gps = df.input
 
#Subsample - make a column that summarizes time as every 30 seconds
dat.gps$date.time.breaks = cut(dat.gps$date.time, breaks = "30 sec")

#order
dat.gps <- dat.gps[with(dat.gps, order(date.time)),]

#remove duplicates from 30 sec time interval (based on date.time.breaks)
dat.gps = dat.gps[!duplicated(dat.gps$date.time.breaks),]

# reduce to GPS points
dat.gps = dat.gps %>% tidyr::drop_na(lat)

# Merge trip metrics with gps data 
temp <- dat.trips %>%
          dplyr::select(date.time, trip_id,
                track.start.date,  track.end.date,
                trip.start.date, trip.end.date,
                track.duration.hr, trip.duration.hr)
 
dat.gps = left_join(dat.gps, temp, by = "date.time")   # merge: left_join keeps all rows of dat.gps, and matches 'date.time' 
dim(dat.gps) 

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Simple exploratory plotting of tracks
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# set a projection
prj = " +proj=longlat + ellps = WGS84 + datum=WGS84"

# Change NA to 0 so that you can plot it 
dat.gps$trip_id[is.na(dat.gps$trip_id)] <- 0

# Convert df object to an sf object
sf_locs <- sf::st_as_sf(dat.gps, coords = c("lon","lat")) %>%
  sf::st_set_crs(prj)

ggplot() +
  layer_spatial(sf_locs, size = 0.75, alpha = 1, aes(color = as.factor(trip_id))) +  # specified stage as ID above
  scale_x_continuous(expand = expansion(mult = c(.6, .6))) +
  scale_color_brewer(palette = "Dark2") +
  theme(legend.position = "bottom") +
  ggtitle("Breeding Gentoo penguins",
          subtitle = "tracks_colours represent trips")


#----------------------------------------------------------------------------------
# add deployment duration and distance: Ryan Reisinger (IYNA paper / Github)
#----------------------------------------------------------------------------------
# Calculate distance using GPS data (only the unique points) & draw plot

# blue lines on plot are from first and last GPS points 
# red lines on plot are from splitting tracks into trips based on diving data
# grey points in plot (NA) have been removed from dat.trips 

library(geosphere)

dt <- list()
ids2 <- unique(dat.gps$TagID)

for(i in 1:length(ids2)){
  
  d <- dat.gps[dat.gps$TagID == ids2[i],]
  d.start <- min(d$date.time)
  d.end <- max(d$date.time)
#  dur <- d.end - d.start
#  dur <- as.character(dur)
  d$start <- rep(d.start, nrow(d))
  d$end <- rep(d.end, nrow(d))
#  d$dur <- rep(dur, nrow(d))
#  slon <- rep(d$lon[1], nrow(d))   # use if you do not have nest locations
#  slat <- rep(d$lat[1], nrow(d))   # use if you do not have nest locations
  slon <- rep(d$nest.lon[1], nrow(d))   # use if you have nest locations
  slat <- rep(d$nest.lat[1], nrow(d))   # use if you have nest locations
  s <- as.matrix(cbind(slon, slat))
  e <- as.matrix(cbind(d$lon, d$lat))
  d$col.dist <- distGeo(s, e)/1000     # distance, in km
  d$dist.max <- max(d$col.dist)
 # d = subset(d, d$col.dist < 200)
  t = ggplot(d, aes(x = date.time, y = col.dist, col =  as.factor(trip_id)))+
    geom_point() + 
    geom_hline(yintercept = 0)+ 
    geom_vline(xintercept = d$start, col = 'blue')+
    geom_vline(xintercept = d$end, col = 'blue')+
    geom_vline(xintercept = d$trip.start.date, col = 'red', linetype="dotted")+     # individual trips start
    geom_vline(xintercept = d$trip.end.date, col = 'red',  linetype="dotted")+       # individual trips end
    geom_vline(xintercept = d$track.start.date, col = 'red', linetype="dotted")+        # track start
    geom_vline(xintercept = d$track.end.date, col = 'red', linetype="dotted")+        #track end
    ggtitle(d$id[1])
   print(t)  
  dt[i] <- list(d)
}

gps.tripdata <- do.call(rbind, dt)
as.data.frame(head(gps.tripdata))
dim(gps.tripdata)

saveRDS(gps.tripdata, (paste0("./outputs/data/GPS points/", "CS_clean_gps_", dat$id[1], '.rds')))

rstudioapi::savePlotAsImage(paste0("./outputs/plots/chinstrap/Kopaitic GPS processing/", "CS_distance_", dat$id[1],".png"), width=900,height=600)


#------------------------------------------------------------------------------
# Calculate duration and distance using dat.trips, drawing no plot to save time  
#------------------------------------------------------------------------------

dt <- list()
ids2 <- unique(dat.trips$TagID)

for(i in 1:length(ids2)){
  
  d <- dat.trips[dat.trips$TagID == ids2[i],]
  d.start <- min(d$date.time)
  d.end <- max(d$date.time)
#  dur <- d.end - d.start
#  dur <- as.character(dur)
  d$start <- rep(d.start, nrow(d))
  d$end <- rep(d.end, nrow(d))
#  d$dur <- rep(dur, nrow(d))
#  slon <- rep(d$lon[1], nrow(d))   # use if you do not have nest locations
#  slat <- rep(d$lat[1], nrow(d))   # use if you do not have nest locations
  slon <- rep(d$nest.lon[1], nrow(d))   # use if you have nest locations
  slat <- rep(d$nest.lat[1], nrow(d))   # use if you have nest locations
  s <- as.matrix(cbind(slon, slat))
  e <- as.matrix(cbind(d$lon, d$lat))
  d$col.dist <- distGeo(s, e)/1000
  d$dist.max <- max(d$col.dist)
  # t = ggplot(d, aes(x = date.time, y = col.dist, col =  as.factor(trip_id)))+
  #   geom_point() + 
  #   geom_hline(yintercept = 0)+ 
  #   geom_vline(xintercept = d$start, col = 'blue')+
  #   geom_vline(xintercept = d$end, col = 'blue')+
  #   geom_vline(xintercept = d$trip.start.date, col = 'red', linetype="dotted")+     # individual trips start
  #   geom_vline(xintercept = d$trip.end.date, col = 'red',  linetype="dotted")+       # individual trips end
  #   geom_vline(xintercept = d$track.start.date, col = 'red', linetype="dotted")+        # track start
  #   geom_vline(xintercept = d$track.end.date, col = 'red', linetype="dotted")+        #track end
  #   ggtitle(d$id[1])
  #  print(t) 
  dt[i] <- list(d)
}

all.tripdata <- do.call(rbind, dt)
as.data.frame(head(all.tripdata))

saveRDS(all.tripdata, (paste0("./outputs/data/", "CS_clean_", dat$id[1], '.rds')))

# ----------------------------------------------
# Plot with mapview
# ----------------------------------------------
# convert data to an sf object:

library(sf)
library(mapview)
mapviewOptions(fgb = FALSE)

# set map projection
wgs84 = "+proj=longlat +ellps=WGS84 +datum=WGS84"

# make simple feature geometry for trip locations
sf <- 
   gps.tripdata  %>%
   st_as_sf(coords = c('lon', 'lat'), crs= wgs84)

# make simple feature geometry for nest location
sf.nest <-  
        gps.tripdata  %>%
        st_as_sf(coords = c('nest.lon', 'nest.lat'), crs= wgs84)

movementmap = 
   sf %>% 
     mapview(zcol= "trip_id",               
             cex = 5, lwd = 0.5, alpha = 1,
            # col.regions=list("red", "green", "yellow", "blue"), # select trip_id colours here if you want
             na.color = "transparent",
             map.types = 'Esri.WorldImagery',
             crs = wgs84) +
   sf.nest %>% 
        mapview(cex = 5, lwd = 0.5, alpha = 1, pch = 15,
                color = "white",     # line colour for nest
                col.regions = "red")  # fill colour for nest

movementmap

#Save it - but takes about 5 min

mapshot(movementmap, file = 
      paste0("./outputs/plots/chinstrap/Kopaitic GPS processing/", "CS_mapview_", dat$id[1],".png"), width=90,height=60)


#works!
#------------------------------------------------------------
# Plot distance and depth together
#------------------------------------------------------------
gps.tripdataNA = gps.tripdata
library(naniar)
gps.tripdataNA = gps.tripdataNA  %>% replace_with_na(replace = list(trip_id = 0))
unique(gps.tripdataNA$trip_id)

plot1 = ggplot()+
    geom_point(data = gps.tripdataNA, aes(x = date.time, y = col.dist, col =  as.factor(trip_id))) + 
    geom_hline(yintercept = 0)+ 
 #   geom_vline(xintercept = unique(gps.tripdata$track.start.date), col = 'blue')+
  #  geom_vline(xintercept = unique(gps.tripdata$track.end.date), col = 'blue')+
    geom_vline(xintercept = unique(gps.tripdata$trip.start.date), col = 'red', linetype="dotted")+     # individual trips start
    geom_vline(xintercept = unique(gps.tripdata$trip.end.date), col = 'red',  linetype="dotted")+       # individual trips end
    ggtitle(d$id[1]) +
    theme_minimal() +
    theme(axis.title.x = element_blank(),
        axis.text.x = element_blank())

plot2 = ggplot(all.tripdata, aes(x = date.time, y = depth, col =  as.factor(trip_id)))+
    geom_point() + 
    geom_hline(yintercept = 0)+ 
    theme_minimal() +
    geom_vline(xintercept = unique(all.tripdata$trip.start.date), col = 'red', linetype="dotted")+     # individual trips start
    geom_vline(xintercept = unique(all.tripdata$trip.end.date), col = 'red',  linetype="dotted")       # individual trips end

# first method
mx <- ymd_hms(max(gps.tripdata$date.time))
mn <- ymd_hms(min(gps.tripdata$date.time))

library(gridExtra)

grid.arrange(plot1 + scale_x_datetime(limits=c(mn,mx)),
             plot2 + scale_x_datetime(limits=c(mn,mx)))

rstudioapi::savePlotAsImage(paste0("./outputs/plots/chinstrap/Kopaitic GPS processing/", "CS_distance_dive_", dat$id[1],".png"), width=900,height=600)



# -----------------sanity check-----------------------------
# Compare analysis on all data vs gps data only:

mean(all.tripdata$dist.max)
mean(gps.tripdata$dist.max)

mean(all.tripdata$col.dist)
mean(gps.tripdata$col.dist)

dim(all.tripdata)
dim(gps.tripdata)

# Are there any NA in the all.tripdata df?
# colnames(all.tripdata)[!complete.cases(t(all.tripdata))]

# -----------------------------------
# trip package analysis
# getting maximum distance per trip 
# -----------------------------------


library(trip)

d <- gps.tripdata[ , c("id", "date.time", "lon", "lat", "trip_id")]
coordinates(d) <- ~lon+lat
proj4string(d) <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")
tr <- trip(d, c("date.time", "trip_id")) #create trip object
tr

grd <- makeGridTopology(tr, cellsize = c(.25, .25)) #create grid topology
trg <- tripGrid(tr, grid = grd) #create time spent sp data.frame
library(raster)
tsr <- raster(trg) #create a raster
tsr <- tsr/(60*60) #convert to hours
#writeRaster(tsr, "tsr.grd", format = "raster", overwrite = TRUE) #write to file

test <- rasterToPoints(tsr)

tsdf <- as.data.frame(trg) #create a data frame
names(tsdf) <- c("t", "lon", "lat")
tsdf$t <- tsdf$t/(60*60)

#distances
dsts <- dplyr::select(gps.tripdata, date.time, trip_id)
dsts$d <- trackDistance(tr, prev = F, longlat = T)

dis.id <- unique(dsts$trip_id)
dis.l <- list()

for(i in 1:length(dis.id)){
  a <- dsts[dsts$trip_id == dis.id[i], ]
  b <- sum(a$d)
  dis.l[i] <- list(b)
}

dis <- cbind(dis.id, unlist(dis.l))
dis <- as.data.frame(dis)
names(dis) <- c("id", "d")
dis$d = as.numeric(dis$d)
dis$distance <- round(dis$d , digits = 0)
dis

# maximum distance from home, per trip:
homedist(tr, home = NULL)

# -----------------sanity check-----------------------------


