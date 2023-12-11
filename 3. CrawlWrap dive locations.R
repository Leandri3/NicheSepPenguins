#------------------------------------------------------------------
# This code  
# speed filter GPS data (remove unreliable locations) with trip::sda
# crawlWrap dive locations of multiple animals with a function (momentuHMM)
# extract temperature in bins from tdr dive data

# Leandri de Kock, Chris Oosthuizen, Andy Lowther
# December 2021

#------------------------------------------------------------------

# load packages
library(momentuHMM)
library(sf)
library(tidyverse)

Sys.setenv(TZ = "GMT")

# -------------------------------------------------------------------------
#### read in the GPS data for a group of penguins
# -------------------------------------------------------------------------
# raw GPS data was processed to trips (locations at the nest removed)
gpsdat = readRDS("./outputs/data/combined data/CS_Kop_all_Bro_files_combined.rds")

# make sure date.time is set to GMT (or UTC)
attr(gpsdat$date.time, "tzone") # check time zone
gpsdat$date.time = lubridate::with_tz(gpsdat$date.time, "GMT")  
attr(gpsdat$date.time, "tzone") # check time zone

# Process data
gpsdat = gpsdat %>%
  as_tibble(gpsdat) %>%
  dplyr::select(-X, -Y, -Z, -depth, -Temp, -TagID, -tag)  %>%   # remove unnecessary columns
  mutate(species = "chinstrap" ) %>% # update species labels if needed ("gentoo" or "chinstrap")
  dplyr::group_by(id) %>%
  arrange(date.time) %>% # put the data in date order
  distinct(date.time, .keep_all=TRUE) %>%   # remove duplicate date.time records for individuals
  dplyr::ungroup()

gpsdat
unique(gpsdat$id)

#-------------------------------------------------------------------------------
####  Speedfilter
#-------------------------------------------------------------------------------
# quick plot to see whether there are outliers

# first automatically define the number of colors you need
library(RColorBrewer)
nb.cols <- length(unique(gpsdat$id))
mycolors <- colorRampPalette(brewer.pal(8, "Set2"))(nb.cols)

ggplot() + 
  geom_point(data = gpsdat, size = 0.25, alpha = 1, aes(x = lon, y = lat, color = id)) +
  scale_x_continuous(expand = expansion(mult = c(.1, .1))) +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.1))) +
  scale_fill_manual(values = mycolors) +
  theme_bw()

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
  scale_x_continuous(expand = expansion(mult = c(.1, .1))) +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.1))) +
  scale_fill_manual(values = mycolors) +
  theme_bw()


#-----------------------------------------------------------------------------------
# Import divestats
#-----------------------------------------------------------------------------------
# List divestats files in the directory that you want
divefiles <- list.files(path = "C:/Users/cwoosthuizen/Desktop/krilltokt penguins/outputs/divemove/Kopaitic chinstrap",
          pattern = '.*divestats_chinstrap_Kopaitic_Bro.*\\.rds$',
          full.names = TRUE)

print(divefiles)

# make empty list
divestatslist = list()

for (i in divefiles) {             # to run on all individuals selected in divefiles
# Import divestats data and save in a list
divestatslist[[i]] = readRDS(i)
}

divestats = as_tibble(bind_rows(divestatslist))
divestats

# make sure date.time is set to GMT (or UTC)
attr(divestats$begdesc, "tzone") # check time zone
divestats$begdesc = lubridate::with_tz(divestats$begdesc, "GMT")  
attr(divestats$begdesc, "tzone") # check time zone

# check that it is the same animal IDs:
unique(gpsdat$id)
unique(divestats$ID)

# --------------------------------------------------------
# run crawl model as a function, to predict latlon locations at dives
# --------------------------------------------------------

crawldat = gpsdat %>%
  dplyr::select(id, date.time, lat, lon)  %>%
  dplyr::rename(ID = id) %>%  # crawl requires variable called ID
  sf::st_as_sf(coords=c("lon", "lat"), crs=4326) %>% # gpsdata (lat lon) is in wgs84 (4326)
  sf::st_transform(32721) # transform lat lon data to UTM zone 21S https://epsg.io/32721 (crawl requires UTM)

print(crawldat)

# specify error around GPS points: assume a 50 m isotropic error ellipse for the measurement error model 
# source: September 2, 2021 momentuHMM vignette
lnError <- crawl::argosDiag2Cov(50,50,0) 
crawldat$ln.sd.x = lnError$ln.sd.x
crawldat$ln.sd.y = lnError$ln.sd.y
crawldat$error.corr = lnError$error.corr


# set up crawlwrap function to predict locations for dives
crawl_dive = function(gps = x, dives = y) {

 crawlpredLoc_dives = crawlWrap(gps,
                    theta =c(5.5,-0.2), 
                    Time.name = "date.time",
                    predTime = dives$begdesc,
                    fixPar = c(1,1,NA,NA),
                    err.model = list(x = ~ln.sd.x-1,
                                     y = ~ln.sd.y-1,
                                    rho = ~error.corr),
                    attempts = 10)
                          
 predLoc_dives = prepData(crawlpredLoc_dives)

 return(predLoc_dives)
 
  }


# Split data frames into lists of individuals (required for pmap)
gpslist  = crawldat %>%
  group_split(ID)
gpslist

divelist = divestats %>%
  group_split(ID)
divelist

# combine lists into 1 object
gps_divelist <- list(gpslist, divelist)

# run crawl model on each individual
cr.pred = purrr::pmap(gps_divelist, crawl_dive)         
 
# convert lists to data tibble
crwdat = dplyr::bind_rows(cr.pred) # Output from prepdata (crawlWrap), (includes PREDICTED locations, x, y, steps and turning angles) for each individual separately with separate data frames

# save model outputs
crawldat$group = substr(crawldat$ID , 1, nchar(crawldat$ID)-5)
unique(crawldat$group)
crawldat$group = substr(crawldat$ID , 1, nchar(crawldat$ID)-6)
unique(crawldat$group)

saveRDS(crwdat, paste0("./outputs/crawl/crawl_divelocs_", gpsdat$species[1],"_", crawldat$group[1], '.rds'))

# add lat lon locations from crawl model to divestats table
divestats_xy = left_join(divestats, crwdat, by = c("ID" = "ID", "begdesc" = "date.time"))

# save divestats table with lat lon added
saveRDS(divestats_xy, paste0("./outputs/crawl/divestats_xy_", gpsdat$species[1],"_", crawldat$group[1], '.rds'))

dim(crwdat)    # should be the same nr of rows
dim(divestats)
dim(divestats_xy)

# plot dive locations
p =  ggplot(data = divestats_xy, aes(x = x, y = y, color = ID)) +
  geom_point(size = 1, alpha = 1) +
  scale_x_continuous(expand = expansion(mult = c(.1, .1))) +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.1))) +
  scale_fill_manual(values = mycolors) +
  theme_bw() 

print(p)

