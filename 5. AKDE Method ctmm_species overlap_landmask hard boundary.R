#----------------------------------------------------------------
# Fit ctmm utilization distribution to krilltokt GPS data.
# Exclude land areas from UDs using a land mask
# Use crawl GPS data

#!!!
# REDo ctmm on the same grid so that we can compare between species
# Kopaitic - Brood


# Leandri 
# July 2023
#------------------------------------------------------------------

#------------------------------------
# 1. Make landmask for ctmm UD
#------------------------------------

# ibcso 500 m Southern Ocean bathemetry data
# https://doi.pangaea.de/10.1594/PANGAEA.937574?format=html#download

library(raster)
# 
# # map data of entire Southen Ocean
# ibcso = raster("./data/Environmental inputs/IBCSO_v2_bed.tif")   
# ibcso
# 
# #---------------------------------------------------------
# # Set up spatial limits of area of interest and projections
# #---------------------------------------------------------
# 
# # set limits of study area in lat lon
# minx <- -69.0       # extends map to west  
# maxx <- -55.0       # extends map to east
# miny <- -68.0       # extends map to south
# maxy <- -56.0       # extends map to north
# 
# library(sf)
# 
# # read in bathymetry raster of WAP (from IBCSO)
ibcso = raster('./outputs/data/Environmental inputs/IBCSO_v2_WAP.tif')
raster::plot(ibcso)


library(RColorBrewer)
blue.col <- colorRampPalette(brewer.pal(9, "Blues"))
plot(ibcso, col = rev(blue.col(255)))

# set cells above sea level to NA and others to 1
ibcso[ibcso >= 0] <- NA
# ibcso[ibcso < 0] <- 1
raster::plot(ibcso, col = rev(blue.col(255)))

# 
# # ---------------------------------------------
# # make a SpatialPolygonsDataFrame of raster
# # ---------------------------------------------
# 
# # see here to remove the outer boundary:
# # https://gis.stackexchange.com/questions/441888/raster-to-spatialpolygonsdataframe-close-crop-polygons-at-the-raster-boundary/441891#441891
# # but this does not work to model the UDs, so keep the landmask that I had in the question, not in the answer. 
# 
library(stars)
library(sf)
# 
# # make a SpatialPolygonsDataFrame 
# spdf = as_Spatial(st_as_sf(stars::st_as_stars(ibcso),
#                         as_points = F, merge = T))
# 
# raster::plot(spdf, col="red", border="green")
# # transform projection
# spdfUTM = spTransform(spdf, CRS("+init=epsg:32721")) # to match that of ctmm::projection(DATA)
# raster::plot(spdfUTM, col="red", border="green")
# 
# saveRDS(spdfUTM, "./data/Environmental inputs/ctmmlandmask.rds")

spdfUTM = readRDS("./outputs/data/Environmental inputs/ctmmlandmask.rds")
raster::plot(spdfUTM, col="red", border="green")

# make the grid the same for both species 
grid = raster(spdfUTM, res = 500)
grid

##---------------
# Gentoos
##---------------
#------------------------------------
#  2. Fit ctmm UDs with landmask 
#------------------------------------

# ctmm model based on code from Chris Fleming: https://groups.google.com/g/ctmm-user/c/JUJLZI-1GL8 
# parallelization based on code from Ingo  https://groups.google.com/g/ctmm-user/c/22EkOdfNLdM/m/JGqds6qxAAAJ 

library(ctmm)
library(doParallel)
library(foreach) 
library(tidyverse)

# read in penguin GPS data, with locations predicted every two minutes
# This is Kopaitic Island
dat = readRDS("./outputs/crawl/5 min crawl_model_preddat_GT_Kopaitic.rds")

plot(dat$x, dat$y, pch = '.')
unique(dat$stage)

# clean up and select a single breeding stage
dat = dat %>% 
  dplyr::select(ID, date.time, x, y, stage) %>% 
  dplyr::mutate(island = "kop") %>% 
  dplyr::rename(trip_id = ID, lon.x = x, lat.y = y) %>%
  #  dplyr::filter(stage == "Kopaitic_Inc")
  dplyr::filter(stage == "Kopaitic_Bro")
#  dplyr::filter(stage == "Kopaitic_Cre")
# create ID variable

# Make an animal ID from the trip_id data. ctmm must run over animals (tracks), not individual foraging trips
dat$id = substr(dat$trip_id, 1, 18) 
unique(dat$id)

# Input data should be in the Movebank format. This requires LAT LON (wgs84) coordinates!
names(dat)

# Back transform to wgs84
library(sp)
#To assign a known CRS to spatial data:
utm.coord = SpatialPoints(cbind(dat$lon.x, dat$lat.y), proj4string=CRS('EPSG:32721'))
# To transform from one CRS to another:
wgs.coord <- spTransform(utm.coord, CRS('EPSG:4326'))
dat$longitude <- wgs.coord$coords.x1
dat$latitude <- wgs.coord$coords.x2
str(dat)

unique(dat$id)

dat = dat %>% 
  dplyr::select(id, date.time, longitude, latitude, lon.x, lat.y) %>% 
  dplyr::rename(x = lon.x, y = lat.y) %>%
  dplyr::group_by(id) %>% 
  dplyr::arrange(date.time) %>%
  dplyr::distinct(date.time, .keep_all = TRUE) %>%
  ungroup()

# define projection
DATA <- as.telemetry(dat, projection = 'EPSG:32721')
ctmm::projection(DATA)

## back-end for foreach function (this works for Windows)
cl <- parallel::makeCluster(detectCores(), outfile = "")
doParallel::registerDoParallel(cl)

#---------------------------------------------------------------
# Import SpatialPolygonsDataframes to use as a hard boundary
#--------------------------------------------------------------

# rename the landmask made from bathymetry
COAST = spdfUTM

GRID = raster(COAST)

FITS = readRDS("./outputs/ctmm/FITS_GT_Kop_Bro.rds")


# plot with tracks
length(DATA)
library(RColorBrewer)
# Define the number of colors you want
nb.cols <- length(DATA)
COL <- colorRampPalette(brewer.pal(8, "Dark2"))(nb.cols)
plot(DATA, col = COL, SP = COAST, col.SP = 'grey') 

#----------------
# Run akde 
#----------------

# Use weights = F for evenly spaced data; weights = T for uneven times

# create a 'akde' function
ud_function <- function(j){
  akde(DATA[[j]], FITS[[j]][[1]], SP = COAST, SP.in = T,  weights = F, grid = grid)
}

UDS <- foreach(j=1:length(DATA), .packages='ctmm') %dopar% 
  {ud_function(j)}

summary(UDS[[1]])

# Mean home range for gentoos based on individual HR
# this is a straight mean of the individual densities that doesn't model population variance
MEAN = mean(UDS)

# # plot coastal boundary in the background
# plot(DATA, MEAN, 
#      col = COL, 
#      col.DF = 'orange', 
#      lwd.level = 2, 
#      SP = COAST, col.SP = 'grey') 
# 
# # plot bathymetry in the background - not entirely correct
# plot(DATA, MEAN, 
#      col = COL, 
#      col.DF = 'darkorange', 
#      lwd.level = 2,
#      R = ibcso, col.R = rev(blue.col(255)))

extent(UDS)
# OVER = overlap(UDS) # pairwise comparisons between individuals

# OVER$CI

# #individual level
# cs_uds = readRDS("./outputs/ctmm/UDS_GT_Kop_Bro.rds")
# summary(cs_uds[[1]])
# ind_ovr <- overlap(cs_uds)
# mean(cs_uds)

## population kernel density estimate (paper coming)
# Whether SP.in = T or SP,in = F does not seem to matter
PKDE <- pkde(DATA, UDS, SP = COAST, SP.in = F, smooth = T, grid = grid) # distribution of the population, grid = grid
# ?pkde
summary(PKDE)
plot(PKDE, SP=COAST, col.SP='grey') 
# # plot coastal boundary in the background
# plot(DATA, PKDE, 
#      col = COL, 
#      col.DF = 'darkorange', 
#      lwd.level = 2, 
#      SP = COAST, col.SP = 'grey') 
# 
# # plot bathymetry in the background - not entirely correct
# plot(DATA, PKDE, 
#      col = COL, 
#      col.DF = 'darkorange', 
#      lwd.level = 2,
#      R = ibcso, col.R = rev(blue.col(255)))


extent(PKDE)

# > extent(PKDE)
# x       y
# min -439537.1 2409034
# max 1079962.9 3792034
# > 

##-------------------------------
# Chinstraps
##---------------
# read in penguin GPS data, with locations predicted every two minutes
# This is Kopaitic Island
dat2 = readRDS("./outputs/crawl/5 min crawl_model_preddat_CS_Kopaitic.rds")

plot(dat2$x, dat2$y, pch = '.')
unique(dat2$stage)

# clean up and select a single breeding stage
dat2 = dat2 %>% 
  dplyr::select(ID, date.time, x, y, stage) %>% 
  dplyr::mutate(island = "kop") %>% 
  dplyr::rename(trip_id = ID, lon.x = x, lat.y = y) %>%
  #  dplyr::filter(stage == "Kopaitic_Inc")
  dplyr::filter(stage == "Kopaitic_Bro")
#  dplyr::filter(stage == "Kopaitic_Cre")
# create ID variable

# Make an animal ID from the trip_id data. ctmm must run over animals (tracks), not individual foraging trips
dat2$id = substr(dat2$trip_id, 1, 18) 
unique(dat2$id)

# Input data should be in the Movebank format. This requires LAT LON (wgs84) coordinates!
names(dat2)

# Back transform to wgs84
library(sp)
#To assign a known CRS to spatial data:
utm.coord = SpatialPoints(cbind(dat2$lon.x, dat2$lat.y), proj4string=CRS('EPSG:32721'))
# To transform from one CRS to another:
wgs.coord <- spTransform(utm.coord, CRS('EPSG:4326'))
dat2$longitude <- wgs.coord$coords.x1
dat2$latitude <- wgs.coord$coords.x2
str(dat2)

unique(dat2$id)

dat2 = dat2 %>% 
  dplyr::select(id, date.time, longitude, latitude, lon.x, lat.y) %>% 
  dplyr::rename(x = lon.x, y = lat.y) %>%
  dplyr::group_by(id) %>% 
  dplyr::arrange(date.time) %>%
  dplyr::distinct(date.time, .keep_all = TRUE) %>%
  ungroup()

# define projection
DATA2 <- as.telemetry(dat2, projection = 'EPSG:32721')
ctmm::projection(DATA2)

# ## back-end for foreach function (this works for Windows)
# cl <- parallel::makeCluster(detectCores(), outfile = "")
# doParallel::registerDoParallel(cl)

#---------------------------------------------------------------
# Import SpatialPolygonsDataframes to use as a hard boundary
#--------------------------------------------------------------

# rename the landmask made from bathymetry
#COAST = spdfUTM

#GRID = raster(COAST)

FITS2 = readRDS("./outputs/ctmm/FITS_CS_Kop_Bro.rds")


# plot with tracks
length(DATA2)
library(RColorBrewer)
# Define the number of colors you want
nb.cols <- length(DATA2)
COL2 <- colorRampPalette(brewer.pal(8, "Dark2"))(nb.cols)
plot(DATA2, col = COL2, SP = COAST, col.SP = 'grey') 

#----------------
# Run akde 
#----------------

# Use weights = F for evenly spaced data; weights = T for uneven times

# create a 'akde' function
ud_function <- function(j){
  akde(DATA2[[j]], FITS2[[j]][[1]], SP = COAST, SP.in = T,  weights = F, grid = grid)
}

UDS2 <- foreach(j=1:length(DATA2), .packages='ctmm') %dopar% 
  {ud_function(j)}

summary(UDS2[[1]])

# Mean home range for gentoos based on individual HR
# this is a straight mean of the individual densities that doesn't model population variance
MEAN2 = mean(UDS2)
MEAN2
# # plot coastal boundary in the background
# plot(DATA2, MEAN2, 
#      col = COL2, 
#      col.DF = 'orange', 
#      lwd.level = 2, 
#      SP = COAST, col.SP = 'grey') 
# 
# # plot bathymetry in the background - not entirely correct
# plot(DATA2, MEAN2, 
#      col = COL2, 
#      col.DF = 'darkorange', 
#      lwd.level = 2,
#      R = ibcso, col.R = rev(blue.col(255)))

extent(UDS2)
OVER2 = overlap(UDS2) # pairwise comparisons between individuals

OVER2$CI

# #individual level
# cs_uds = readRDS("./outputs/ctmm/UDS_GT_Kop_Bro.rds")
# summary(cs_uds[[1]])
# ind_ovr <- overlap(cs_uds)
# mean(cs_uds)

## population kernel density estimate (paper coming)
# Whether SP.in = T or SP,in = F does not seem to matter
PKDE2 <- pkde(DATA2, UDS2, SP = COAST, SP.in = F, smooth = T, grid = grid) # distribution of the population

summary(PKDE2)
plot(PKDE2, SP=COAST, col.SP='grey') 
# # plot coastal boundary in the background
# plot(DATA2, PKDE2, 
#      col = COL2, 
#      col.DF = 'darkorange', 
#      lwd.level = 2, 
#      SP = COAST, col.SP = 'grey') 
# 
# # plot bathymetry in the background - not entirely correct
# plot(DATA2, PKDE2, 
#      col = COL2, 
#      col.DF = 'darkorange', 
#      lwd.level = 2,
#      R = ibcso, col.R = rev(blue.col(255)))


extent(PKDE2)
#extent(UDS2)
extent(PKDE)
# extent(UDS)

# #set extent so that it is the same despite hr-sizes
# EXT = extent(PKDE)
# 
# par(mfrow = c(1,2))
# plot(PKDE, col.DF = 'darkorange', SP=COAST, col.SP='lightblue', level = 0.95, level.UD = 0.95, ext = EXT)
# title ('gentoo')
# plot(PKDE2, col.DF = 'green', SP=COAST, col.SP='lightblue', level = 0.95, level.UD = 0.95, ext = EXT)
# title('chinstraps')
# 
# par(mfrow = c(1,1))

# how to combine pkdes of both species to assess overlap? 
#://groups.google.com/g/ctmm-user/c/NQn95t_M8tg/m/J_sOaDlABQAJ
# https://groups.google.com/g/ctmm-user/c/L7911xHoXCE

# combine lists of gentoo and chinstrap from the same island and breeding stage into one list
kop_bro <- c(list(gt=PKDE), list(cs=PKDE2))
#projection(kop_bro) <-median(kop_bro)
# then you can use the overlap() function on that compiled list 
# 95% overlap - BA
kop_bro_overlap_95 <- overlap(kop_bro, level = 0.95)  #, ext = EXT
kop_bro_overlap_95

ba_95 = kop_bro_overlap_95$CI

write.csv(ba_95, './outputs/ctmm/BA species 95_overlap results_Kopaitic_Bro.csv')

# 50% overlap - BA
kop_bro_overlap_50 <- overlap(kop_bro, level = 0.50)  #, ext = EXT
kop_bro_overlap_50

ba_50 = kop_bro_overlap_50$CI

write.csv(ba_50, './outputs/ctmm/BA species 50_overlap results_Kopaitic_Bro.csv')





