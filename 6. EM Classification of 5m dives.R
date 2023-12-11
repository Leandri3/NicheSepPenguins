#----------------------------------------------------------------------------------------
# This code:
# Estimates dive residuals (using linear mixed-effects models)
# Applies an Expectation-Maximization algorithm to 
# Classify diving types into 2 clusters (foraging and non-foraging) per species, per island
#
# Using the following predictors:Bottom time
#                                Dive residuals 
#                                Maximum depth
#
# Separates into three different Clusters: 1 = Foraging
#
# Nelson Island

# Leandri de Kock, Chris Oosthuizen, Andy Lowther
#  August 2022
#-----------------------------------------------------------------------------------------
#setwd("C:/Users/cwoosthuizen/Desktop/krilltokt penguins")

# load packages
library(tidyverse)
library(lme4)
library(ggplot2)
library(Rmixmod)

Sys.setenv(TZ = "GMT")

##-------------------------------------------------------------------------
# 1. Nelson, Chinstraps 
##-------------------------------------------------------------------------

# Import divestats
chins <- list.files(path = "C:/Users/cwoosthuizen/Desktop/krilltokt penguins/outputs/crawl",
          pattern = '.*divestats_xy_chinstrap_Nelson.*\\.rds$',
          full.names = TRUE)
print(chins)

# make empty list
chinslist = list()

for (i in chins) {             # to run on all individuals selected in divefiles
# Import divestats data and save in a list
chinslist[[i]] = readRDS(i)
}

divestats = as_tibble(bind_rows(chinslist))
divestats

divestats = divestats %>%
  mutate(species = "chinstrap")

head(divestats)
names(divestats)

# Make maxdep positive:
divestats$maxdep = divestats$maxdep*-1
divestats$maxdep 

# filter for dives deeper than 5m
divestats_5m = divestats %>%
  filter(maxdep > 5)
# filters out about 60 000 dives 

#-------------------------------------------------------------------------------------------------
# Get dive residuals for classification using linear mixed effects models
#------------------------------------------------------------------------------------------------
# Now use a linear mixed effects model to account for random effects of different individuals: 
#library(lme4)

# Rescaling variables: around a mean = 0 and sd = 1

#divetim.stand (stand = standardised)
divestats_5m$divetim.stand = (divestats_5m$divetim - (mean(divestats_5m$divetim)))/sd(divestats_5m$divetim)
mean(divestats_5m$divetim.stand)
sd(divestats_5m$divetim.stand)

#maxdep.stand
divestats_5m$maxdep.stand = (divestats_5m$maxdep - (mean(divestats_5m$maxdep)))/sd(divestats_5m$maxdep)
mean(divestats_5m$maxdep.stand)
sd(divestats_5m$maxdep.stand)

# random effects: ID and depth
model1 <- lmer(divetim.stand ~ maxdep.stand + (maxdep.stand|ID), data = divestats_5m, REML = T) 
# when rescaling variables beforehand, there are no warnings about model failed to converge, model unidentifiable
model1
summary(model1)
AIC(model1)
plot(model1)

# Get the dive residuals from the model
residuals(model1,  type = 'pearson', scaled = T)
divestats_5m$dive.res = residuals(model1, type = 'pearson')
divestats_5m

# Plot the dive residuals against depth
ggplot(divestats_5m, aes(x = maxdep.stand, y = dive.res)) + geom_point() + geom_smooth(method = "lm")
ggplot(divestats_5m, aes(x = maxdep, y = dive.res)) + geom_point() + geom_smooth(method = "lm")

# Negative residual values: Relatively shorter dives at a given depth that may indicate relatively lower effort
# Positive residual values: Relatively longer dives at a given depth that may indicate relatively lower effort

# To see correlation matrix code: see EM Penguins_Kop_CS.R script

#===========================   APPLY EM ALGORITHM  ============================ 

#library(Rmixmod) # see instruction of the package 
# Select variables to be included in the EM 
# (variables that will be used to classify dives into different types) 
# Here you need to decide which of the divestats variables are important to classify dives.
# see ??diveStats for meaning of the dive parameters

# # Including 3 variables to classify dives:  

data_cluster = divestats_5m %>%
                   dplyr::select(botttim, # bottom time
                                 #divetim, # dive duration
                                #bottdist, # indication of wiggles
                                 dive.res,  #,)
                                 maxdep)
data_cluster
# check for NA
colnames(data_cluster)[!complete.cases(t(data_cluster))]  

# change NA values to zero, otherwise algorithm fails (NA for bottom time = no bottom time, = 0 seconds)
data_cluster[is.na(data_cluster)] = 0
colnames(data_cluster)[!complete.cases(t(data_cluster))]  

# Run the EM model with 3 cluster types:  -1 = foraging; 2 & 3 = non-foraging

system.time(        # system.time is useful when running big data. Use <- in association with system.time as an = gives an error
    emdive <- mixmodCluster(data = data_cluster, 
                            nbCluster = 3,    # how many dive types are there?
                            models = mixmodGaussianModel (), 
                            criterion= c("BIC","ICL"),
                            strategy = mixmodStrategy(algo = "EM", 
                                                      nbTry = 1,
                                                      initMethod = "random", 
                                                      nbTryInInit = 100,
                                                      nbIterationInInit = 5, 
                                                      nbIterationInAlgo = 200,
                                                      epsilonInInit = 0.001, 
                                                      epsilonInAlgo = 0.001)))

# show a summary of the best model containing the estimated parameters , the likelihood
summary(emdive)

#dev.new()  # open a window for plot. Does not always plot properly in the plotting box
plot(emdive) 

# Can identify 3 clusters: 
# Foraging (Cluster 1, Longest bottom time, Shallow to deep dives, high effort)
# Non-foraging (Cluster 2, Short bottom time, medium to deep dives, expected effort)
# Non-foraging (Cluster 3, short bottom time, shallow dives,less than expected effort)

# Save plot 
rstudioapi::savePlotAsImage(paste0("./outputs/plots/chinstrap/Nelson clusters/", "CS_EM_5m_3 Clusters.png"), width=900,height=600)

# add clusters to data 
divestats_5m$dive_5m_cluster<-(emdive@bestResult@partition)
head(divestats_5m)

# save divestats_5m table with clusters added
saveRDS(divestats_5m, paste0("./outputs/clusters/divestats_5m_xy_clusters_", divestats_5m$species[1],"_nel",'.rds'))

# Exploration only:
# make own plots: max depth by cluster (dive type)
ggplot(data = divestats_5m, aes(x = seq(1:length(maxdep)), y = maxdep, color = as.factor(dive_5m_cluster))) +
  geom_point(size = 1, alpha = 0.9, shape=16) +
  theme(plot.margin = margin()) +
  theme_bw()

# make own plots: bottom time by cluster (dive type)
ggplot(data = divestats_5m, aes(x = seq(1:length(maxdep)), y = botttim, color = as.factor(dive_5m_cluster))) +
  geom_point(size = 1, alpha = 0.9, shape=16) +
  theme(plot.margin = margin()) +
  theme_bw()

# make own plots: bottom time by cluster (dive type)
ggplot(data = divestats_5m, aes(x = seq(1:length(maxdep)), y = bottdist, color = as.factor(dive_5m_cluster))) +
  geom_point(size = 1, alpha = 0.9, shape=16) +
  theme(plot.margin = margin()) +
  theme_bw()

ggplot(data = divestats_5m, aes(x = x, y = y, color = as.factor(dive_5m_cluster))) +
  geom_point(size = 1, alpha = 0.9, shape='.') +
  theme(plot.margin = margin()) +
  theme_bw()

# Save plot 
rstudioapi::savePlotAsImage(paste0("./outputs/plots/chinstrap/Nelson clusters/", "CS_5m_Tracks by 3 Clusters.png"), width=900,height=600)


################################################################################################################

#-------------------------------------------------------------------
# 2. Nelson, Gentoos
#----------------------------------------------------------------
# Import divestats
gents <- list.files(path = "C:/Users/cwoosthuizen/Desktop/krilltokt penguins/outputs/crawl",
          pattern = '.*divestats_xy_gentoo_Nelson.*\\.rds$',
          full.names = TRUE)
print(gents)

# make empty list
gentslist = list()

for (i in gents) {             # to run on all individuals selected in divefiles
# Import divestats data and save in a list
gentslist[[i]] = readRDS(i)
}

divestats = as_tibble(bind_rows(gentslist))
divestats

divestats = divestats %>%
  mutate(species = "gentoo")

head(divestats)
names(divestats)

# Make maxdep positive:
divestats$maxdep = divestats$maxdep*-1
divestats$maxdep 

# filter for dives deeper than 5m
divestats_5m = divestats %>%
  filter(maxdep > 5)
# filters out about 14 000 dives 

#-------------------------------------------------------------------------------------------------
# Get dive residuals for classification using linear mixed effects models
#------------------------------------------------------------------------------------------------
# Now use a linear mixed effects model to account for random effects of different individuals: 
#library(lme4)

# Rescaling variables: around a mean = 0 and sd = 1

#divetim.stand (stand = standardised)
divestats_5m$divetim.stand = (divestats_5m$divetim - (mean(divestats_5m$divetim)))/sd(divestats_5m$divetim)
mean(divestats_5m$divetim.stand)
sd(divestats_5m$divetim.stand)

#maxdep.stand
divestats_5m$maxdep.stand = (divestats_5m$maxdep - (mean(divestats_5m$maxdep)))/sd(divestats_5m$maxdep)
mean(divestats_5m$maxdep.stand)
sd(divestats_5m$maxdep.stand)

# random effects: ID and depth
model1 <- lmer(divetim.stand ~ maxdep.stand + (maxdep.stand|ID), data = divestats_5m, REML = T) 
# when rescaling variables beforehand, there are no warnings about model failed to converge, model unidentifiable
model1
summary(model1)
AIC(model1)
plot(model1)

# Get the dive residuals from the model
residuals(model1,  type = 'pearson', scaled = T)
divestats_5m$dive.res = residuals(model1, type = 'pearson')
divestats_5m

# Plot the dive residuals against depth
ggplot(divestats_5m, aes(x = maxdep.stand, y = dive.res)) + geom_point() + geom_smooth(method = "lm")
ggplot(divestats_5m, aes(x = maxdep, y = dive.res)) + geom_point() + geom_smooth(method = "lm")

# Negative residual values: Relatively shorter dives at a given depth that may indicate relatively lower effort
# Positive residual values: Relatively longer dives at a given depth that may indicate relatively lower effort

# To see correlation matrix code: see EM Penguins_Kop_CS.R script

#===========================   APPLY EM ALGORITHM  ============================ 

#library(Rmixmod) # see instruction of the package 
# Select variables to be included in the EM 
# (variables that will be used to classify dives into different types) 
# Here you need to decide which of the divestats variables are important to classify dives.
# see ??diveStats for meaning of the dive parameters

# # Including 3 variables to classify dives:  

data_cluster = divestats_5m %>%
                   dplyr::select(botttim, # bottom time
                                 #divetim, # dive duration
                                #bottdist, # indication of wiggles
                                 dive.res,  #,)
                                 maxdep)
data_cluster
# check for NA
colnames(data_cluster)[!complete.cases(t(data_cluster))]  

# change NA values to zero, otherwise algorithm fails (NA for bottom time = no bottom time, = 0 seconds)
data_cluster[is.na(data_cluster)] = 0
colnames(data_cluster)[!complete.cases(t(data_cluster))]  

# Run the EM model with 3 cluster types:  -1 = foraging; 2 & 3 = non-foraging

system.time(        # system.time is useful when running big data. Use <- in association with system.time as an = gives an error
    emdive <- mixmodCluster(data = data_cluster, 
                            nbCluster = 3,    # how many dive types are there?
                            models = mixmodGaussianModel (), 
                            criterion= c("BIC","ICL"),
                            strategy = mixmodStrategy(algo = "EM", 
                                                      nbTry = 1,
                                                      initMethod = "random", 
                                                      nbTryInInit = 100,
                                                      nbIterationInInit = 5, 
                                                      nbIterationInAlgo = 200,
                                                      epsilonInInit = 0.001, 
                                                      epsilonInAlgo = 0.001)))

# show a summary of the best model containing the estimated parameters , the likelihood
summary(emdive)

#dev.new()  # open a window for plot. Does not always plot properly in the plotting box
plot(emdive) 

# Can identify 3 clusters: 
# Foraging (Cluster 1, Long bottom time, Shallow to deep dives, high effort)
# Non-foraging(Cluster 2, Short bottom time, shallow to deep dives, expected effort)

# Save plot 
rstudioapi::savePlotAsImage(paste0("./outputs/plots/gentoo/Nelson clusters/", "GT_EM_5m_3 Clusters.png"), width=900,height=600)

# add clusters to data 
divestats_5m$dive_5m_cluster<-(emdive@bestResult@partition)
head(divestats_5m)
#divestats_5m$dive_5m_cluster_new = as.factor(divestats_5m$dive_5m_cluster)

## change number of clusters so that 1 = foraging 
#new_cluster <- divestats_5m %>%
#  select(ID, begdesc, dive_5m_cluster)%>%
#  mutate(dive_5m_cluster_new = as.factor(dive_5m_cluster)) 
#
#new_cluster <- new_cluster %>%
#  mutate(dive_5m_cluster_new = replace(dive_5m_cluster_new, dive_5m_cluster_new == 1, 3)) %>% 
#  mutate(dive_5m_cluster_new = replace(dive_5m_cluster_new, dive_5m_cluster_new == 2, 1))
#
## everything = 1 foraging, others  = non-foraging
#new_cluster
#
#divestats_5m_new <- left_join (divestats_5m, new_cluster, by = c('ID' = "ID", "begdesc" = "begdesc", "dive_5m_cluster" = "dive_5m_cluster" ))
#names(divestats_5m_new)


# save divestats_5m table with clusters added
saveRDS(divestats_5m, paste0("./outputs/clusters/divestats_5m_xy_clusters_", divestats_5m$species[1],"_nel",'.rds'))

# Exploration only:
# make own plots: max depth by cluster (dive type)
ggplot(data = divestats_5m, aes(x = seq(1:length(maxdep)), y = maxdep, color = as.factor(dive_5m_cluster))) +
  geom_point(size = 1, alpha = 0.9, shape=16) +
  theme(plot.margin = margin()) +
  theme_bw()

# make own plots: bottom time by cluster (dive type)
ggplot(data = divestats_5m, aes(x = seq(1:length(maxdep)), y = botttim, color = as.factor(dive_5m_cluster))) +
  geom_point(size = 1, alpha = 0.9, shape=16) +
  theme(plot.margin = margin()) +
  theme_bw()

# make own plots: bottom time by cluster (dive type)
ggplot(data = divestats_5m, aes(x = seq(1:length(maxdep)), y = bottdist, color = as.factor(dive_5m_cluster))) +
  geom_point(size = 1, alpha = 0.9, shape=16) +
  theme(plot.margin = margin()) +
  theme_bw()

ggplot(data = divestats_5m, aes(x = x, y = y, color = as.factor(dive_5m_cluster))) +
  geom_point(size = 1, alpha = 0.9, shape=16) +
  theme(plot.margin = margin()) +
  theme_bw()

# Save plot 
rstudioapi::savePlotAsImage(paste0("./outputs/plots/gentoo/Nelson clusters/", "GT_5m_Tracks by 3 Clusters.png"), width=900,height=600)
####################################################################################################################

