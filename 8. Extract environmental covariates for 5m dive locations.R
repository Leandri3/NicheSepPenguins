
# Nelson and Kopaitic Island 2018/2019
# Extracting remote sensed environmental covariates 
# Leandri de Kock, Chris Oosthuizen
# August 2022

#library(rgdal)
library(raster)
library(terra)
library(marmap)
library(RColorBrewer)
library(rasterVis)
library(fields)
library(sp)
library(tidyverse)
library(viridis)
library(tidync)
library(ncdf4)
library(tidyr)
library(sf)

Sys.setenv(TZ = "GMT")

#---------------------------------------------------------
# Set up spatial limits of area of interest and projections
#---------------------------------------------------------
# Set projections
utm.prj = " +proj=utm +zone=21 +south +datum=WGS84 +units=m +no_defs "   # Chris UTM King George
wgs84 <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"


# set limits of study area in lat lon
minx <- -62.0       # extends map to west  
maxx <- -56.0       # extends map to east
miny <- -67.0       # extends map to south
maxy <- -58.0       # extends map to north

# make sf object of limits (in lat/lon)
lim = data.frame(x = c(minx, maxx),  y = c(miny, maxy)) %>% 
  st_as_sf(coords = 1:2, crs = 4326)

# transform limits to polar stereo, to crop the bathemetry data with these limits
lim_polar = st_transform(lim, 9354)
lim_polar


#---------------------------------------------------------
# 1. Bathymetry
#---------------------------------------------------------
#~~~
#UPGRADED to a new bathymetry version from GEBCO to IBSCO
#~~~

# IBSCO 500 m Southern Ocean bathymetry data
# https://doi.pangaea.de/10.1594/PANGAEA.937574?format=html#download

# This imports IBSCO bathymetery values for the entire Southern Ocean: 
# The digital bathymetric model (DBM) of IBCSO Version 1.0 has a 500m x 500m resolution 
# based on a Polar Stereographic projection for the area south of 60Â° S.
# https://www.scar.org/science/ibcso/resources/

#---------------------------------------------------------
# 1. Load tif of bathymetry data and crop to study area 
#---------------------------------------------------------
# The key thing is to download the correct tif (not the RGB tif)

bath = raster("./data/Environmental inputs/IBCSO_v2_bed.tif")   
bath
#plot(bath)

# crop to study area     
e <- as(extent(lim_polar), 'SpatialPolygons')
class(e)
proj4string(e) <- CRS("+init=epsg:9354")

# crop     
ibsco <- crop(bath, e)
plot(ibsco)
ibsco[ibsco > 0] <- 0
ibsco

#---------------------------------------------------------
# 1a. create a landmask from bathemetry
#---------------------------------------------------------
landmask <- ibsco
landmask[landmask >= 0] <- NA
class(landmask)
plot(landmask)

# reproject raster to UTM
# make a new landmask spatraster
landmask2 = terra::rast(landmask)
ibsco_utm = terra::project(landmask2, "EPSG:32721", method = "near")
plot(ibsco_utm)

# reproject raster to WGS84
ibsco_wgs = terra::project(landmask2, "EPSG:4326", method = "near")
plot(ibsco_wgs)

#---------------------------------------------------------
# 2. Slope
#---------------------------------------------------------

slope <- terrain(landmask, opt ="slope")
#plot(slope)
slope = projectRaster(slope, crs = wgs84)
plot(slope)

#---------------------------------------------------------
# 3a. distance to shelf
#---------------------------------------------------------
# takes a few minutes to run

shelf <- rasterToContour(landmask, levels = -200) # bathymetry contour to select
shelf <- rasterize(shelf, y = raster(ext = extent(landmask), crs = crs(landmask)))
shelf.dist <- distance(shelf)
shelf.dist <- shelf.dist/1000   # convert to km
landmask1 <- landmask
landmask1 <- resample(x = landmask1, y = shelf.dist)
shelf.dist <- mask(shelf.dist, mask=landmask1)
plot(shelf.dist, colNA= "black")

saveRDS(shelf.dist, "./outputs/environmental space covariates/shelf.dist 200 m.rds")

## Read in shelf once saved (so that you don't have to run it all the time)
#shelf.dist = read_rds("./outputs/environmental space covariates/shelf.dist 200 m.rds")
#shelf.dist = projectRaster(shelf.dist, crs = wgs84)
#plot(shelf.dist)

#----------------------------------------------------------------------------------------------
# 3b. Distance to shelf with negative values over the shelf and positive values in the deep
#----------------------------------------------------------------------------------------------

# Must multiply depth raster by distance to shelf raster
# Distance to shelf raster is in a different resolution, so have to make the resolutions the same
# (different resolutions because of resample)
#
 landmask   # check resolution of bathymetry raster
 shelf.dist # check resolution
 
 # Make resoltions the same:
 shelf.dist.l <- resample(shelf.dist, landmask)
 shelf.dist.l  # same res as landmask
 
 # set function
 rc <- function(landmask, shelf.dist.l)
    {    ifelse(landmask >= -200,  -1,   # if landmask is shallower or = to -200, * -1 (On-shelf locations)
         ifelse(landmask <  -200,  1) )   # if landmask is deeper than -200 * 1 (off-shelf locations)
    }
 
 #Apply function to raster stack
 r.class <- overlay(landmask, fun=rc)
 plot(r.class, colNA = "black")       # plot with on shelf / off-shelf dictomy (values are either 1 or -1)
 
 # Now we have a raster 'mask' that says where the bathymetry is deeper than 100 m (= 1) and
 # where it is not (= -1)
 # Now multiply this raster with the shelf.dist raster (with appropriate resolution):
 
 shelf.dist.posneg  = r.class * shelf.dist.l
 plot(shelf.dist.posneg, colNA = "black")
 shelf.dist.posneg 
 
 saveRDS(shelf.dist.posneg, "./outputs/environmental space covariates/shelf.dist.posneg 200 m.rds")
 
## Read in shelf once saved (so that you don't have to run it all the time)
#shelf.dist.posneg = read_rds("./outputs/environmental space covariates/shelf.dist.posneg 200 m.rds")
#shelf.dist.posneg = projectRaster(shelf.dist.posneg, crs = wgs84)
#plot(shelf.dist.posneg)

#---------------------------------------------------------
# 4. distance to coast for each track location 
#---------------------------------------------------------
 coast <- rasterToContour(ibsco, levels = 0)  # bathymetry contour to select (sea-level)
 coast <- rasterize(coast, y = raster(ext = extent(ibsco), crs = crs(ibsco)))
 coast.dist <- distance(coast)
 coast.dist <- coast.dist/1000   # convert to km
 landmask1 <- landmask
 landmask1 <- resample(x = landmask1, y = coast.dist)
 coast.dist <- mask(coast.dist, mask=landmask1)
 plot(coast.dist, colNA= "black")
 
 saveRDS(coast.dist, "./outputs/environmental space covariates/distancetocoast.rds")

#coast.dist = read_rds("./outputs/environmental space covariates/distancetocoast.rds")
#coast.dist = projectRaster(coast.dist, crs = wgs84)
#plot(coast.dist)

#---------------------------------------------------------
# 5. SST
#---------------------------------------------------------
# Import SST data from 2018-11-25 to 2019-02-14

# SST data was downloaded using the script Downloading podaac nasa 1 km SST data.R
# https://cran.r-project.org/web/packages/heatwaveR/vignettes/OISST_preparation.html

#OISST_data <- readRDS("./data/Environmental inputs/krilltokt_OISSTdata20181101_20190228.Rds")

#~~~
# UPGRADED to a newer version of OISST data: 1 km ultra high resolution SST data from Erddap and podaac
#~~~~
# with a grid resolution of 0.01 (1 km x 1 km). 
# Data downloaded with: Downloading podaac nasa 1 km SST data script
# Data source: GHRSST Level 4 MUR Global Foundation Sea Surface Temperature Analysis (v4.1)
# https://podaac.jpl.nasa.gov/dataset/MUR-JPL-L4-GLOB-v4.1
#~~~

OISST_data <- readRDS("./data/Environmental inputs/krilltokt_OISSTdata_1km_25nov_14feb.Rds")
class(OISST_data)
OISST_data

# check if the dates are correct
summary(OISST_data)   # should be from 25 November 2018 to 14 February 2019

# Look at data 
# ggplot() +
#    geom_tile(data = OISST_data  , aes(x = lon, y = lat, fill = temp))  +
#    scale_fill_viridis(direction = 1, option = "plasma",
#    name = "SST")

#ggplot(OISST_data, aes(lon, lat, fill = temp)) + geom_raster() + facet_wrap(~t)
# it clearly becomes warmer in the Bransfield Strait from January

# convert sst to raster
# First transform tibble from long to wide, so that each 'date' is a column. 
# This is needed to make a raster stack with a raster per day
# Check: all 'days' have 2270 data points. 
Nn = OISST_data  %>%
   group_by(t) %>%
   summarise(n = length(t))
Nn 
unique(Nn$n)
unique(Nn$t) # 82 unique days

OISST.wide = OISST_data  %>% 
  spread(t, temp)
OISST.wide

# Now create a rasterBrick from regularly gridded points: rasterFromXYZ
OISST.r <- rasterFromXYZ(OISST.wide, crs = wgs84, digits = 2)
OISST.r

indx <- as.data.frame(unique(OISST_data$t))
colnames(indx) = "t"
indx$date.time <- as.POSIXct(indx$t, format = "%Y-%m-%d")
#indx$Date <- as.Date(indx$date.time, tz = "GMT")
indx$date.names <- format(indx$date.time, "%Y%m%d")

# Make a layer for each day
OISST.rt <- setZ(OISST.r, indx$date.time)
OISST.rt 
names(OISST.rt) <- indx$date.names 
OISST.rt 
plot(OISST.rt)

plot(OISST.rt[[1]])

# # Get a mean SST value for the period you are working with: 
# # Statistics across cells: 
# #mean_sst <-cellStats((OISST.rt[[40:48]]), stat = mean) # for all days
# #mean_sst
# 
# # Use calc() function to calculate values for a new raster* object from another raster* object, using a formula
# # CS nel bro: 26 December - 14 February
# mean_sst_chick <- calc(OISST.rt[[32:82]], mean)
# mean_sst_chick
# 
# # now we have a mean_sst value for each cell during the brood period
# plot(mean_sst_chick)

#---------------------------------------------------------
# 6. CURRENT
#---------------------------------------------------------
# Import Current data from 2018-11-01 to 2019-02-28

# Current data was downloaded using the script Krilltokt download current (CURRENT).R
# NOAA Near Real Time Geostrophic Currents Data

CURRENT_data <- readRDS("./data/Environmental inputs/krilltokt_CURRENTdata.Rds")
class(CURRENT_data)
CURRENT_data

#plot(CURRENT_data$lon, CURRENT_data$lat)

# check if the dates are correct
summary(CURRENT_data)   # should be from November 2018 to February 2019

#CURRENT_data$u_current = round(CURRENT_data$u_current, 2)
#CURRENT_data$v_current = round(CURRENT_data$v_current, 2)
CURRENT_data$lat = round(CURRENT_data$lat, 2)
CURRENT_data$lon = round(CURRENT_data$lon, 2)
CURRENT_data

# Look at data 
 ggplot() +
    geom_tile(data = CURRENT_data  , aes(x = lon, y = lat, fill = v_current))  +
    scale_fill_viridis(direction = 1, option = "plasma",
    name = "V-curent")
#
#ggplot(CURRENT_data, aes(lon, lat, fill = temp)) + geom_raster() + facet_wrap(~t)
# it clearly becomes warmer in the Bransfield Strait from January

# convert CURRENT to raster
# First transform tibble from long to wide, so that each 'date' is a column. 
# This is needed to make a raster stack with a raster per day
# Check: all 'days' have 2611 data points. 
Nn = CURRENT_data  %>%
   group_by(t) %>%
   summarise(n = length(t))
Nn 
unique(Nn$n)
unique(Nn$t) # 120 unique days

#---------------------------------------------------------
# 6a. U_current: Eastward 
#---------------------------------------------------------
CURRENT_data_u = CURRENT_data %>%
  dplyr::select(lon, lat, t, u_current)
CURRENT_data_u

CURRENT.wide_u = CURRENT_data_u  %>% 
  spread(t, u_current)
CURRENT.wide_u

# Now create a rasterBrick from regularly gridded points: rasterFromXYZ
CURRENT.r_u <- rasterFromXYZ(CURRENT.wide_u, res = c(0.2, 0.2), crs = wgs84)
CURRENT.r_u
# Error: x cell sizes are not regular == fixed with rounding all variables to 2 decimals


indx <- as.data.frame(unique(CURRENT_data_u$t))
colnames(indx) = "t"
indx$date.time <- as.POSIXct(indx$t, format = "%Y-%m-%d")
#indx$Date <- as.Date(indx$date.time, tz = "GMT")
indx$date.names <- format(indx$date.time, "%Y%m%d")

# Make a layer for each day
CURRENT.rt_u <- setZ(CURRENT.r_u, indx$date.time)
CURRENT.rt_u 
names(CURRENT.rt_u) <- indx$date.names 
CURRENT.rt_u 
plot(CURRENT.rt_u)

plot(CURRENT.rt_u[[1]])

# # Get a mean current_U value for the period you are working with: 
# # Statistics across cells: 
# #mean_current_u <-cellStats((CURRENT.rt_u[[49:97]]), stat = mean) # for all days
# #mean_current_u
# 
# # Use calc() function to calculate values for a new raster* object from another raster* object, using a formula
# # Brood period for cs_Nelson = 26 December - 14 February
# mean_current_u_chick <- calc(CURRENT.rt_u[[56:106]], mean)
# mean_current_u_chick
# 
# # now we have a mean_sst value for each cell during the brood period
# plot(mean_current_u_chick)


#--------------------------------------------------------------------------
# 6b. V_current: Northward 
#------------------------------------------------------------------------
CURRENT_data_v = CURRENT_data %>%
  dplyr::select(lon, lat, t, v_current)
CURRENT_data_v

CURRENT.wide_v = CURRENT_data_v  %>% 
  spread(t, v_current)
CURRENT.wide_v

# Now create a rasterBrick from regularly gridded points: rasterFromXYZ
CURRENT.r_v <- rasterFromXYZ(CURRENT.wide_v, res = c(0.2, 0.2), crs = wgs84)
CURRENT.r_v
# Error: x cell sizes are not regular == fixed with rounding all variables to 2 decimals

indx <- as.data.frame(unique(CURRENT_data_v$t))
colnames(indx) = "t"
indx$date.time <- as.POSIXct(indx$t, format = "%Y-%m-%d")
#indx$Date <- as.Date(indx$date.time, tz = "GMT")
indx$date.names <- format(indx$date.time, "%Y%m%d")

# Make a layer for each day
CURRENT.rt_v <- setZ(CURRENT.r_v, indx$date.time)
CURRENT.rt_v 
names(CURRENT.rt_v) <- indx$date.names 
CURRENT.rt_v 
plot(CURRENT.rt_v)

plot(CURRENT.rt_v[[1]])

#  Get a mean current_U value for the period you are working with: 
#  Statistics across cells: 
# mean_current_v <-cellStats((CURRENT.rt_v[[49:97]]), stat = mean) # for all days
# mean_current_v

# # Use calc() function to calculate values for a new raster* object from another raster* object, using a formula
# # Brood period for cs_Nelson = 26 December - 14 February
# mean_current_v_chick <- calc(CURRENT.rt_v[[56:106]], mean)
# mean_current_v_chick
# 
# # now we have a mean_sst value for each cell during the brood period
# plot(mean_current_v_chick)

#---------------------------------------------------------
# 7. fsle 
#---------------------------------------------------------

# https://www.aviso.altimetry.fr/en/data/products/value-added-products/fsle-finite-size-lyapunov-exponents.html

# Log into Myaviso
# wcoosthuizen@zoology.up.ac.za
# SO1sit
# Go to my products 
# https://www.aviso.altimetry.fr/en/my-aviso-plus/my-products.html

# Navigate to value added products, and FSLE, and se Extraction 
# Should open a website that will load (a bit slow) 
# -----------------------------------------------------------------------------
# Download 
# DATASET-DUACS-DT-GLOBAL-ALLSAT-MADT-FSLE--AGGREGATED-DATASET
# (delayed time) 
# for a specific block. 
# DESCRIPTION
# FSLEs based on the maximum eigenvalue of the Cauchy-Green strain tensor
# UNIT	
# days-1
# -----------------------------------------------------------------------------
# This is better than using ftp (fileZilla) where you cannot select a geographic area.

#----------------------------------
## Import via raster package 
#----------------------------------
# use 'raster' to get a rasterlayer with multiple 'bands'
# But best us to use 'brick' and not 'stack' to get rasterbrick.

fsle = brick("./data/Environmental inputs/dataset-duacs-dt-global-allsat-madt-fsle_1647859809653.nc")
fsle
# plot first 4 days
plot(fsle[[1:4]])

# # Notice that longitude is not given as negative (but as 0 to 360)
# # https://gis.stackexchange.com/questions/266535/change-a-raster-from-longitude-display-180-180-to-0-360
# # Use rotate to change a raster from longitude display (-180, 180) to (0, 360)

fsle = raster::rotate(fsle)
plot(fsle[[1:4]])
fsle 

# Get a mean current_U value for the period you are working with: 
# Statistics across cells: 
#mean_fsle <-cellStats((fsle[[49:97]]), stat = mean) # for all days
#mean_fsle

## Use calc() function to calculate values for a new raster* object from another raster* object, using a formula
## Brood period for cs_Nelson = 26 December - 14 February
#mean_fsle_chick <- calc(fsle[[56:106]], mean)
#mean_fsle_chick
#
## now we have a mean_sst value for each cell during the brood period
#plot(mean_fsle_chick)

#---------------------------------------------------------
# 9. Distance to nest
#---------------------------------------------------------
#xy <- c(-58.42, -62.21)  # lat lon of Nelson or Nelson
#nest.dist <- distanceFromPoints(landmask, xy) 
#nest.dist <- nest.dist/1000   # convert to km
#nest.dist <- mask(nest.dist, mask=landmask)
#plot(nest.dist)


#----------------------------------------------------------
# 10. Salinity and potential temperature 
#-----------------------------------------------------------

# Data was downloaded from CMEMS (Copernicus Marine service): 
# https://resources.marine.copernicus.eu/product-detail/GLOBAL_REANALYSIS_PHY_001_031/INFORMATION 

# Products: Global Ocean Ensemble Physics Reanalysis
# GLOBAL_REANALYSIS_PHY_001_031


# At the surface? To characterise the water masses
#--------------------------------------------------------------------
sal <-brick("./data/Environmental inputs/Nelson_global-reanalysis-phy-001-031-grepv2-daily_1659881110354.nc", varname = "so_foam", level =1)
sal

plot(sal[[10:18]])

class(sal)

# Get a mean sal value for the period you are working with: 
# Statistics across cells: 
#mean_sal_chick <-cellStats((sal[[19:67]]), stat = mean) # for all days
#mean_sal_chick

# # Use calc() function to calculate values for a new raster* object from another raster* object, using a formula
# # Brood period for cs_Nelson = 26 December - 14 February
# # ! different dates to subset to for salinity, because we have data for a shorter period of time
# mean_sal_chick <- calc(sal[[26:76]], mean)
# mean_sal_chick
# 
# # now we have a mean_sst value for each cell during the brood period
# plot(mean_sal_chick)

#---------------------------------------------------------
# 11. Density ocean mixed layer thickness (units = m)
#---------------------------------------------------------

# Data was downloaded from CMEMS (Copernicus Marine service): 
# https://resources.marine.copernicus.eu/product-detail/GLOBAL_REANALYSIS_PHY_001_031/INFORMATION 

# Products: Global Ocean Ensemble Physics Reanalysis
# GLOBAL_REANALYSIS_PHY_001_031

# - #ocean_mixed_layer_thickness_defined_by_sigma_theta (MLD) - also from CMEMS
#--------------------------------------------------------
mld <-brick("./data/Environmental inputs/Nelson_global-reanalysis-phy-001-031-grepv2-daily_1659881110354.nc", varname = "mlotst_foam", level =1)
mld

plot(mld[[10:18]])
class(mld)

# Get a mean mld value for the period you are working with: 
# Statistics across cells: 
#mean_mld_chick <-cellStats((mld[[19:67]]), stat = mean) # for all days
#mean_mld_chick

## Use calc() function to calculate values for a new raster* object from another raster* object, using a formula
## Brood period for 31 December to 07 February
## ! different dates to subset to for salinity, because we have data for a shorter period of time
#mean_mld_chick <- calc(mld[[26:76]], mean)
#mean_mld_chick

## now we have a mean_MLD value for each cell during the brood period
#plot(mean_mld_chick)

#---------------------------------------------------------
# 12. Potential temperature
#---------------------------------------------------------

# Data was downloaded from CMEMS (Copernicus Marine service): 
# https://resources.marine.copernicus.eu/product-detail/GLOBAL_REANALYSIS_PHY_001_031/INFORMATION 

# Products: Global Ocean Ensemble Physics Reanalysis
# GLOBAL_REANALYSIS_PHY_001_031

# - sea_water_potential_temperature (T)- also from CMEMS
#--------------------------------------------------------
pot<-brick("./data/Environmental inputs/Nelson_global-reanalysis-phy-001-031-grepv2-daily_1659881110354.nc", varname = "thetao_foam", level =1)
pot

plot(pot[[10:18]])
class(pot)

## Get a mean pot value for the period you are working with: 
## Statistics across cells: 
##mean_pot_chick <-cellStats((pot[[19:67]]), stat = mean) # for all days
##mean_pot_chick
#
## Use calc() function to calculate values for a new raster* object from another raster* object, using a formula
## Brood period for cs_Nelson = 31 December to 07 February
## ! different dates to subset to for salinity, because we have data for a shorter period of time
#mean_pot_chick <- calc(pot[[26:76]], mean)
#mean_pot_chick
#
## now we have a mean_pot value for each cell during the brood period
#plot(mean_pot_chick)


#---------------------------------------------------------
# Extractions - Gentoos
#---------------------------------------------------------
# Get data of penguin lat lon positions for which you want an environmental value

# 1. Nelson Gentoos
gent_nel <- readRDS("./outputs/summary/dive_5m_metrics_foraging_gentoo_Nelson.rds")

# gent_nel <-list.files(path = "C:/Users/cwoosthuizen/Desktop/krilltokt penguins/outputs/crawl",
#           pattern = '.*divestats_temp_gentoo_Nelson.*\\.rds$',
#           full.names = TRUE)
# 
# print(gent_nel)
# 
# # make empty list
# gentslist = list()
# 
# for (i in gent_nel) {             # to run on all individuals selected in divefiles
# # Import divetemps data and save in a list
# gentslist[[i]] = readRDS(i)
# }
# 
# gent_nel = as_tibble(bind_rows(gentslist))
str(gent_nel)

# make sure date.time is set to GMT (or UTC)
attr(gent_nel$begdesc, "tzone") # check time zone

# gent_nel = gent_nel %>%
#   mutate(species = "gentoo")%>%
#   mutate(island = "nelson")
# 
# # create a new column for each breeding stage
# gent_nel$stage = substr(gent_nel$ID , 8, nchar(gent_nel$ID)-6)
# unique(gent_nel$stage)

# Change Gua to Bro
#--------------------------------------------------------------
# rename Gua = Bro
# gent_nel$stage = if_else(gent_nel$stage == "Gua", "Bro", gent_nel$stage )
# unique(gent_nel$stage)

# Look at sample size (uplinks) per stage # bro should be 96906
gent_nel %>%
  group_by(group) %>%
  dplyr::summarise(n= n()) # Bro = 90433
# works!
#----------------------------------------------------------------

plot(gent_nel$x, gent_nel$y, pch = ".", cex = 2)

# 2. Kopaitic gentoos
gent_kop <- readRDS("./outputs/summary/dive_5m_metrics_foraging_gentoo_Kopaitic.rds")

# gent_kop <-list.files(path = "C:/Users/cwoosthuizen/Desktop/krilltokt penguins/outputs/crawl",
#           pattern = '.*divestats_temp_gentoo_Kopaitic.*\\.rds$',
#           full.names = TRUE)
# 
# print(gent_kop)
# 
# # make empty list
# gentslist = list()
# 
# for (i in gent_kop) {             # to run on all individuals selected in divefiles
# # Import divetemps data and save in a list
# gentslist[[i]] = readRDS(i)
# }
# 
# gent_kop = as_tibble(bind_rows(gentslist))
str(gent_kop)

# make sure date.time is set to GMT (or UTC)
attr(gent_kop$begdesc, "tzone") # check time zone

# gent_kop = gent_kop %>%
#   mutate(species = "gentoo")%>%
#   mutate(island = "kopaitic")

# create a new column for each breeding stage
# gent_kop$stage = substr(gent_kop$ID , 10, nchar(gent_kop$ID)-6)
# unique(gent_kop$stage)

# Look at sample size (uplinks) per stage # bro should be 96906
gent_kop %>%
  group_by(group) %>%
  dplyr::summarise(n= n()) # Bro = 90433
# works!
#----------------------------------------------------------------
plot(gent_kop$x, gent_kop$y, pch = ".", cex = 2)

locdat = as_tibble(bind_rows(gent_nel, gent_kop))
head(locdat)

locdat = locdat %>%
  rename(lon.x = x) %>%
  rename(lat.y = y) %>%
  filter(maxdep > 5)

# remove rows with NA values in lat column 
locdat = locdat[complete.cases(locdat$lat.y),] # there were two rows with NaN values

#To assign a known CRS to spatial data:
utm.coord = SpatialPoints(cbind(locdat$lon.x, locdat$lat.y), proj4string=CRS(utm.prj))
utm.coord

# To transform from one CRS to another:
wgs.coord <- spTransform(utm.coord, CRS(wgs84))
wgs.coord

locdat$lon2 <- wgs.coord$coords.x1
locdat$lat2 <- wgs.coord$coords.x2
head(locdat)
plot(locdat$lon2,locdat$lat2, pch = ".", cex = 2)

dat = as_tibble(locdat) #%>%  # add in important dive parameters? 
    # dplyr::select(ID, begdesc, botttim, divetim, bottdist, maxdep, 
    #                lon.x, lat.y, dive.res, divecluster, group, lon2, lat2,
    #                island, species,  
    #                for.dives, prop.for.dives)

dat$Date = as.Date(dat$begdesc)
dat$Date = as.POSIXct(dat$Date, format="%Y-%m-%d",tz= "GMT")

dat
str(dat)

## Extractions (for environmental variables in utm)
coord.proj <- cbind(dat$lon.x, dat$lat.y)
coord.proj <- as.data.frame(coord.proj)
names(coord.proj) <- c("lon.x", "lat.y")
dim(coord.proj)
head(coord.proj)

## Extractions (for environmental variables in lat lon)
coord.wgs84 <- cbind(dat$lon2, dat$lat2)
coord.wgs84 <- as.data.frame(coord.wgs84)
names(coord.wgs84) <- c("lon.wgs84", "lat.wgs84")
dim(coord.wgs84)
head(coord.wgs84)

#------------------------------------------------------
# Extract static environmental variables
#------------------------------------------------------
dat$depth <- raster::extract(x = ibsco_wgs,
                             y = cbind(coord.wgs84$lon.wgs84, coord.wgs84$lat.wgs84))

dat$slope <- raster::extract(x = slope,
                             y = cbind(coord.wgs84$lon.wgs84, coord.wgs84$lat.wgs84))

dat$shelf <- raster::extract(x = shelf.dist,
                             y = cbind(coord.wgs84$lon.wgs84, coord.wgs84$lat.wgs84))

dat$shelf_posneg <- raster::extract(x = shelf.dist.posneg,
                             y = cbind(coord.wgs84$lon.wgs84, coord.wgs84$lat.wgs84))

dat$coast <- raster::extract(x = coast.dist,
                             y = cbind(coord.wgs84$lon.wgs84, coord.wgs84$lat.wgs84))

#------------------------------------------------------
# Extract dynamic environmental variables
#------------------------------------------------------

xtract<-function(stk, coords, dates_track){
  rasterdates <- getZ(stk)
  df<-raster::extract(x=stk, y=coords) #extract values for all rasters  
  vals<-vector("numeric") #empty vector
  for(i in 1:nrow(df)){
    x = which(rasterdates == dates_track[i]) #gives mathcing in stk
    vals[i] <- df[i, x]
  }
  return(vals)
}

dat$sst <- xtract(stk = OISST.rt,
                  coords = cbind(coord.wgs84$lon.wgs84, coord.wgs84$lat.wgs84),
                  dates_track = dat$Date)

dat$u_current <- xtract(stk = CURRENT.rt_u,
                  coords = cbind(coord.wgs84$lon.wgs84, coord.wgs84$lat.wgs84),
                  dates_track = dat$Date)

dat$v_current <- xtract(stk = CURRENT.rt_v,
                  coords = cbind(coord.wgs84$lon.wgs84, coord.wgs84$lat.wgs84),
                  dates_track = dat$Date)

dat$fsle <- xtract(stk = fsle,
                  coords = cbind(coord.wgs84$lon.wgs84, coord.wgs84$lat.wgs84),
                  dates_track = dat$Date)

dat$sal <- xtract(stk = sal,
                  coords = cbind(coord.wgs84$lon.wgs84, coord.wgs84$lat.wgs84),
                  dates_track = dat$Date)

as.data.frame(dat)
dat$depth <- dat$depth$IBCSO_v2_bed

# ------------------------------------------
# Sanity check
# ------------------------------------------

#library(viridis)

ggplot() +
   geom_point(data = dat , aes(x = lon.x, y = lat.y, col = depth), size = 0.5)  +
   scale_color_viridis(option = "D")

ggplot() +
   geom_point(data = dat , aes(x = lon.x, y = lat.y, col = slope), size = 0.5)  +
    scale_color_viridis(option = "D")

ggplot() +
   geom_point(data = dat , aes(x = lon.x, y = lat.y, col = shelf), size = 0.5)  +  # SHELF AT 200m
   scale_color_viridis(option = "D")

ggplot() +
   geom_point(data = dat , aes(x = lon.x, y = lat.y, col = shelf_posneg), size = 0.5)  + #shelf at 200 m
   scale_color_viridis(option = "D")


ggplot() +
   geom_point(data = dat , aes(x = lon.x, y = lat.y, col = coast), size = 0.5)  +
   scale_color_viridis(option = "D")

ggplot() +
   geom_point(data = dat , aes(x = lon.x, y = lat.y,  col = sst), size = 0.5)  +
    scale_color_viridis(option = "D")

ggplot() +
   geom_point(data = dat , aes(x = lon.x, y = lat.y,  col = u_current), size = 0.5)  +
    scale_color_viridis(option = "D")

ggplot() +
   geom_point(data = dat , aes(x = lon.x, y = lat.y,  col = v_current), size = 0.5)  +
    scale_color_viridis(option = "D")

ggplot() +
   geom_point(data = dat , aes(x = lon.x, y = lat.y,  col = fsle), size = 0.5)  +
    scale_color_viridis(option = "D")

ggplot() +
   geom_point(data = dat , aes(x = lon.x, y = lat.y,  col =sal), size = 0.5)  +
    scale_color_viridis(option = "D")

# Are there NA's in the ENVIRONMENTAL data?  
colnames(dat)[!complete.cases(t(dat))]  

# dat[is.na(dat$ice),]  missing values on multiple days
# How many NA's are there?  
sum(is.na(dat$slope))
#sum(is.na(dat$shelf))
sum(is.na(dat$depth))
sum(is.na(dat$shelf_posneg))
#sum(is.na(dat$coast))
sum(is.na(dat$sst))
sum(is.na(dat$u_current))
sum(is.na(dat$v_current))
sum(is.na(dat$sal))

# Hardly any for sst. Simply use the closest available data:

# Funtion to select the nearest data point:
# https://stackoverflow.com/questions/10077415/replacing-nas-in-r-with-nearest-value
#---start ----------------------------------------------------- 
f1 <- function(dat) {
  N <- length(dat)
  na.pos <- which(is.na(dat))
  if (length(na.pos) %in% c(0, N)) {
    return(dat)
  }
  non.na.pos <- which(!is.na(dat))
  intervals  <- findInterval(na.pos, non.na.pos,
                             all.inside = TRUE)
  left.pos   <- non.na.pos[pmax(1, intervals)]
  right.pos  <- non.na.pos[pmin(N, intervals+1)]
  left.dist  <- na.pos - left.pos
  right.dist <- right.pos - na.pos

  dat[na.pos] <- ifelse(left.dist <= right.dist,
                        dat[left.pos], dat[right.pos])
  return(dat)
}
#--end ------------------------------------------------------

dat$sst <- f1(dat$sst)
dat$slope <- f1(dat$slope)
dat$depth <- f1(dat$depth)
#dat$shelf <- f1(dat$shelf)
dat$shelf_posneg <- f1(dat$shelf_posneg)
#dat$coast <- f1(dat$coast)
dat$u_current <- f1(dat$u_current)
dat$v_current <- f1(dat$v_current)
dat$sal <- f1(dat$sal)

# Check: Are there NA's in the data?  # Only dive metrics
colnames(dat)[!complete.cases(t(dat))]  

# Exactly the same plot; sst is just filled in for the 11 NA values near land. 
ggplot() +
   geom_point(data = dat , aes(x = lon.x, y = lat.y,  col = slope))  +
    scale_color_viridis(option = "D")

ggplot() +
   geom_point(data = dat , aes(x = lon.x, y = lat.y,  col = depth))  +
    scale_color_viridis(option = "D")

ggplot() +
   geom_point(data = dat , aes(x = lon.x, y = lat.y,  col = shelf_posneg))  +
    scale_color_viridis(option = "D")

ggplot() +
   geom_point(data = dat , aes(x = lon.x, y = lat.y,  col = shelf))  +
    scale_color_viridis(option = "D")

ggplot() +
   geom_point(data = dat , aes(x = lon.x, y = lat.y,  col = sst))  +
    scale_color_viridis(option = "D")

ggplot() +
   geom_point(data = dat , aes(x = lon.x, y = lat.y,  col = coast))  +
    scale_color_viridis(option = "D")

ggplot() +
   geom_point(data = dat , aes(x = lon.x, y = lat.y,  col = u_current))  +
    scale_color_viridis(option = "D")

ggplot() +
   geom_point(data = dat , aes(x = lon.x, y = lat.y,  col = v_current))  +
    scale_color_viridis(option = "D")
dat

#------------------------------------------------------------------
# Benthic/ pelagic diving: 
# Is a dive pelagic or benthic?
#-----------------------------------------------------------------
#dat$ibcso <- dat$depth$IBCSO_v2_bed
#dat$ibcso
#str(dat$ibcso)
str(dat$depth)

names(dat)

# how to work out percentage use of water column: 
pwc = dat %>%
  dplyr::select(lon.x, lat.y, begdesc, ID, maxdep, depth)
head(pwc)

#pwc$bathy <- raster::extract(x = gebco,
#                             y = cbind(coord.wgs84$lon.wgs84, coord.wgs84$lat.wgs84))

# make depth positive: 
pwc$bathy = pwc$depth*-1

pwc$bathy = replace(pwc$bathy, pwc$bathy == 0, 0.1)

pwc$bathy.prop = pwc$maxdep/pwc$bathy
pwc$bathy.prop

plot(pwc$bathy.prop)
summary(pwc$bathy.prop)

#pwc.5 = pwc %>%
#  filter(bathy.prop > 5)

ggplot() +
   geom_point(data = pwc , aes(x = lon.x, y = lat.y,  col = bathy.prop), size = 0.5)  +
    scale_color_viridis(option = "D")

# how many observations above 100%?
#pwc.over100 = pwc %>%
#  filter(bathy.prop >1)  # 17813 observations

pwc$bathy = round(pwc$bathy.prop, 2)
pwc$bathy <- replace(pwc$bathy, pwc$bathy > 0.80, 'benthic')
pwc$bathy <- replace(pwc$bathy, pwc$bathy <= '0.80', 'pelagic')
pwc$bathy

ggplot() +
   geom_point(data = pwc , aes(x = lon.x, y = lat.y,  col = bathy), size = 0.2)  #+
    #scale_color_viridis(option = "D")

feed = pwc %>%
  group_by(bathy) %>%
  summarise(n = n())
feed

head(pwc)

full_dive_env<- left_join(dat, pwc, by=c("ID" = "ID", "begdesc" = "begdesc", 
                                     "lon.x" = "lon.x", "lat.y" = "lat.y", 
                                     "maxdep" = "maxdep", "depth" = "depth"))
full_dive_env

# Save lat lon extracted environmental data
saveRDS(full_dive_env, './outputs/environmental covariates/gentoos_foraging_5m_dive metrics_and_extracted_environmental_covariates.rds')

