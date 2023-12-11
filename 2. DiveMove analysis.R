
#----------------------------------------------------------
# Leandri de Kock, Chris Oosthuizen 
# October 2021

# diveMove analysis 
# Zero-offset correction
# Calculate dive statistics for all individuals in a deployment round
#----------------------------------------------------------

Sys.setenv(TZ = "GMT")

library(tidyverse)
library(data.table)
library(lubridate)
library(trip)
library(sf)
library(diveMove)

#-----------------------------------------------------------------------------

# List files in the directory
filenames <- list.files(path = "C:/Users/cwoosthuizen/Desktop/krilltokt penguins/outputs/data",
           pattern = "*rds",
           full.names = F)
filenames

# This is gentoos and chinstraps, and in our cleaned data, we have not actually put the 
# species in the data!
# Subset here to only list 1 species for 1 island, to break up the analysis a bit:

# List KOPAITIC CHINSTRAP files in the directory
all_files <- list.files(path = "C:/Users/cwoosthuizen/Desktop/krilltokt penguins/outputs/data",
          pattern = '.*CS_clean_Kopaitic.*\\.rds$',
          full.names = TRUE)

# all_files <- list.files(path = "C:/Users/cwoosthuizen/Desktop/krilltokt penguins/outputs/data",
#           pattern = '.*Kopaitic_Bro1_KI08.*\\.rds$',
#           full.names = TRUE)


# How many files are there in the directory? 
dplyr::n_distinct(all_files)

# print the name of each file
for (file in all_files) {
  print(file)
}


#------------------------------------------
# Import the data from all files:
#------------------------------------------

for (file in all_files) {             # to run on all individuals in the folder
  
# Import chinstrap penguin accelerometer data 
divedat = readRDS(file)

# make depth positive for diveMove
divedat$depth = divedat$depth*-1

# calculate sampling rate in seconds and add it to the plot below as a heading
sampling.rate = unique(round(unique(divedat$timediff)*60))

# plot and save, for quick reference
plot(divedat$date.time, divedat$depth, pch = ".", 
           main = c(paste0("CS_dive_", divedat$id[1]),
           paste0("Sampling rate = ", sampling.rate, " sec")))

rstudioapi::savePlotAsImage(paste0("./outputs/plots/chinstrap/diveMove_Kopaitic/", "CS_dive_", divedat$id[1],".png"), width=900,height=600)


#====================================================================
#====================================================================
# DONT SET depth TO NEGATIVE BEFORE TDR ANALYSIS! 
# The calibrateDepth function works with positive depths! 
#====================================================================
#====================================================================

# diveMove:

#now start creating the TDR dive object
filename = "Dive data"

tdr <- createTDR(time = divedat$date.time, 
                 depth = divedat$depth, 
                 concurrentData = divedat[, !names(divedat) %in% c("date.time","depth")],#OR [, -c(1,2)], #the - means it excludes variables 1 & 2  #or dives[, c(3:7, 9:15)],
                 speed = FALSE, 
                 dtime = 1,   #  sampling interval used in seconds
                 file = filename)

#plotTDR(tdr)
# ****   Make sure 'depth' is positive and starts at zero !!! **** 
show(tdr)

# start detecting periods of activity (i.e individual dives)
tdr.calib <- calibrateDepth(tdr,
                  dive.thr = 5,       # only select dives deeper than 5 m.
                  zoc.method = 'filter',
                  k=c(3, 5760),
                  probs=c(0.5, 0.02),
                  knot.factor=20, 
                  descent.crit.q=0.01, ascent.crit.q=0.0,
                  na.rm=T)

#plotTDR(tdr.calib, diveNo=7:8, what="phases", depth.lim=c(0, 80))
#plotDiveModel(tdr.calib, diveNo=16)

#extractDive(tdr.calib, diveNo=2:8)
# tdr.calib@call
# tdr.calib@tdr
# tdr.calib@gross.activity
# tdr.calib@dive.activity  # getDAct(tdr.calib)

# extract zero-offset corrected data
tdr.dat = as.data.frame(tdr.calib@tdr)
tdr.dat = rename(tdr.dat, date.time = time)  # rename to date.time
head(tdr.dat)

# make depth negative
tdr.dat$depth = tdr.dat$depth * -1

# save zero-offset corrected data
saveRDS(tdr.dat, (paste0("./outputs/divemove/Kopaitic chinstrap/", "ZOC_tdr_chinstrap_", divedat$id[1], '.rds')))

#------------------------------------------------
# create dive summary metrics for each dive 
#------------------------------------------------

dive.stats <- diveStats(tdr.calib)        
head(dive.stats)

# add time stamps columns:
# this is for the track start and end date and probably not very useful for a dive-level
# summary table, so excluded
# stamps <- stampDive(tdr.calib, ignoreZ=F)
# dive.stats <- data.frame(stamps, dive.stat)
# dive.stats

# drop some columns that you will probably not use in analysis:
dive.stats = dive.stats %>% dplyr::select (-c(#phase.no, activity,  # these columns are from 'stampDive
                          descD.1stqu, descD.3rdqu, ascD.1stqu,
                          bottD.1stqu, bottD.3rdqu, ascD.3rdqu))

# select only dives with divetim less than 16 minutes (to exclude any extreme errors)
dive.stats = dive.stats %>% dplyr::filter(divetim < 1000)

# make max dive depth negative
dive.stats$maxdep = dive.stats$maxdep * -1

# how to id dives? Add row numbers to dive data - each row is a new dive
dive.stats$dive.nr <- seq.int(nrow(dive.stats))

# Make a new vector, dive.midpoint, which is the time half-way between the end of descent 
# and start of ascent (i.e., the middle of the dive)
dive.stats$dive.midpoint = dive.stats$begasc - dive.stats$enddesc  # = bottom time 
dive.stats$dive.midpoint = dive.stats$enddesc + 0.5*dive.stats$dive.midpoint
  
# when do dives end? (time at end of dive)
dive.stats$end.dive <- as.POSIXct(dive.stats$asctim, origin = dive.stats$begasc, tz = "GMT") 

# Did this work? Calculate dive duration in secs:
# dive.stats$dive.dur = (dive.stats$end.dive - dive.stats$begdesc) + 0.5
# dive.stats$divetim - dive.stats$dive.dur   # differ with 0.5 secs from that reported from divestats. Add 0.5, then all the same. 
# So yes, it works!
# Don't have to add dive.dur to the df as divemove reports divetim (which is the same)

# How many dives did this penguin do? Add this number to dive.stats
dive.stats$dives.in.trip = as.numeric(length(dive.stats$begdesc))

# How many dives did this penguin do deeper than a certain max depth??
deep.dives = dive.stats %>% 
                filter(maxdep < -5) %>% 
                tally()   # how many records are there?
deep.dives 

#add number of deep dives to dive.stats
dive.stats$deep.dives = as.numeric(rep(deep.dives[1], length(dive.stats$begdesc)))

# proportion of deep dives:
dive.stats$prop.deep.dives = dive.stats$deep.dives / dive.stats$dives.in.trip

# proportion bottom time
dive.stats$prop.botttime = (dive.stats$botttim/as.numeric(dive.stats$divetim))

#------------------------
# dive bouts
#------------------------

# Calculate time difference between end of dive x and start of dive x+1
# Verified that the initial 3 calculations are correct time differences of inter-dive time.
# This is the gap BEFORE the dive
# diveMove reports postdiveDuration which is the gap AFTER the dive.
dive.stats$dive.gap <- difftime(dive.stats$begdesc, lag(dive.stats$end.dive, 
                         default = dive.stats$begdesc[1]), units = "min")

# Make a dive bout, by classifying all dives with a time period of less than 5 minutes into a bout.
# Add 1 at the end, to make the bouts start at bout 1, and not at bout 0.
dive.stats$bout <- cumsum(ifelse(dive.stats$dive.gap > 5, 1, 0)) + 1
# 5 =   how many minutes for a new bout to start?
# 1 =  number that each bout should 'advance' with if yes
# 0 =  number that each bout should 'advance' with if no 

dive.stats$total.bouts <- max(dive.stats$bout)

dive.stats$ID = divedat$id[1]

# add a column: trip_id
 unique(tdr.dat$trip_id)
 
 trip_id = tdr.dat %>%
   select(id, trip_id, date.time)%>%
   group_by(id) %>%
   rename(ID = id)
 trip_id
 
dive.stats = left_join(dive.stats, trip_id, by = c("ID" = "ID", "begdesc" = "date.time"))

unique(dive.stats$trip_id)

saveRDS(dive.stats, (paste0("./outputs/divemove/Kopaitic chinstrap/", "divestats_chinstrap_", divedat$id[1], '.rds')))

# ------------------------------------------------------------------------
# Calculate the mean values per bout, for selected columns in dive.stats
# ------------------------------------------------------------------------

bout.stats <- dive.stats %>%
                 replace(is.na(.), 0)  %>%    # replace NA with zero, otherwise cannot calculate means
                 group_by(bout) %>%
                 summarise(bout.starttime = min(begdesc),
                           bout.endtime = max(end.dive),
                           bout.duration = bout.endtime - bout.starttime,
                           mean.botttim = mean(botttim),
                           mean.divetim = mean(divetim),
                           mean.bottdep = mean(bottdep.mean),
                           mean.maxdep = mean(maxdep),
                           bout.n.dives = sum(unique(dive.nr)))
                           
                             
head(as.data.frame(bout.stats))

saveRDS(bout.stats, (paste0("./outputs/divemove/Kopaitic chinstrap/", "boutstats_chinstrap_", divedat$id[1], '.rds')))

# plot and save, for quick reference
plot(dive.stats$maxdep, dive.stats$divetim, pch = ".",
           main = c(paste0("CS_dive_", divedat$id[1]),
           paste0("Dive time against max depth")))

rstudioapi::savePlotAsImage(paste0("./outputs/plots/chinstrap/diveMove_Kopaitic/", "CS_divetim_maxdep", divedat$id[1],".png"), width=900,height=600)


}    # loop ends 


# test the position of 'dive.midpoint' - works
tst = tdr.dat[110000:112050,]
plot(tst$date.time, tst$depth, type = "l")

# now add max depth from diveMove 
points(dive.stats$enddesc , dive.stats$maxdep, col = "red", pch = 20)   
points(dive.stats$dive.midpoint , dive.stats$maxdep, col = "blue", pch = 20)    # nearer the time at depth from tdr
points(dive.stats$begdesc , rep(0, length(dive.stats$begdesc)) , col = "darkgreen", pch = 20)    # nearer the time at depth from tdr
points(dive.stats$end.dive , rep(0, length(dive.stats$end.dive)) , col = "orange", pch = 20)    # nearer the time at depth from tdr
points(dive.stats$begasc, dive.stats$maxdep, col = "purple", pch = 20)    # nearer the time at depth from tdr

#tst = tdr.dat[9000:11250,]
#plot(tst$date.time, tst$depth, type = "l")

