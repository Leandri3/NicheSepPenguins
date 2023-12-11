# load packages
library(tidyverse)
library(lmerTest)
library(car)
library(ggplot2)
library(effects)
library(sjPlot)
library(lattice)
library(ggeffects)
library(arm)
library(mgcv)
library(glmmTMB)
library(DHARMa)
library(MetBrewer)
library(RColorBrewer)
library(splines)
library(performance)
library(MuMIn)



#-----------------------------------------------------------------
# read in dive metrics
gent <- readRDS('./outputs/environmental covariates/gentoos_foraging_5m_dive metrics_and_env and bathy_covariates_23_05_2023.rds')
chin <- readRDS('./outputs/environmental covariates/chinstraps_foraging_5m_dive metrics_and_env and bathy_covariates_23_05_2023.rds')

# include all dive data from both populations
dive <- as_tibble(bind_rows(gent, chin))

head(dive)
names(dive)


### 1. Sanity check ---------------------------------------------

# Breeding stage check
unique(dive$phase)
unique(dive$group)


names(dive)

# rename some columns
dive = dive %>%
  rename(feed= feeddive)
dive = dive %>%
  rename(bdepth = new_bdepth)   

names(dive)
str(dive)


# make sure date.time is set to GMT (or UTC)
attr(dive$begdesc, "tzone") # check time zone

# how many dive observations per species?
dive.sp = dive %>%
  group_by(species)%>%
  summarise(n = n())
dive.sp

# create a new ID variable with 'species_name'

dive$ids = paste0(dive$species, '_', dive$ID)
dive$ids

unique(dive$ids)  # 221 is correct. 


# make a day variable since 1 Dec
library(lubridate)

dive$date = date(dive$begdesc)
dive$date

summary(dive$date)

# how many days since 1 December? 
dive$num.day=julian(dive$date, origin = as.Date('2018-12-01'))
dive$num.day

summary(dive$num.day)

# # change bathymetry to all values >300 to 300m 
# dive$bathy300 = dive$bdepth
# dive$bathy300 = with(dive, replace(bathy300, bdepth > 300, 300))
# plot(dive$bathy300)

# standardize predictor variables eg. sst 

dive$z.sst = scale(dive$sst)
plot(dive$z.sst)
plot(dive$sst)

dive$z.bdepth = scale(dive$bdepth)

dive$z.solar.ele.begdesc = scale(dive$solar.ele.begdesc)

dive$z.num.day = scale(dive$num.day)


# Calculate 'time' for correlation structure:
g = dive %>% 
  group_by(ids) %>%   
  dplyr::mutate(t = difftime(begdesc, begdesc[1], units='hours'))  

g$t2 = plyr::round_any(as.numeric(g$t), 0.5, f = ceiling) 
length(unique(g$t2))

g$times = glmmTMB::numFactor(as.numeric(g$t2))

# change group levels so that reference level can be Incubation. Then we compare diving depths between Inc 
# and Chick-rearing phases

g$group = relevel(as.factor(g$group), ref = 'Inc') 

## quick correlation test  
# library(corrplot)
# dive.cor = dive %>%
#   dplyr::select(maxdep, z.num.day, z.solar.ele.begdesc, z.bdepth, z.sst)
# str(dive.cor)
# n = cor(dive.cor)
# corrplot(n, method = 'number')
# Num.day & sst = positive correlation 0.60 Thus it is smaller than 0.7 and can still be included in models. 
# all the other variables are not correlated.

# Improve time by using parallelization
parallel::detectCores() #8 

nt <- min(parallel::detectCores(),8)

#------------------------------------------------------------------------
# Model selection
# - make sure that REML = F, so that model selection can be done with ML. 
#------------------------------------------------------------------------

# Modelset 1: With 'group'

# GROUP
mod1 = glmmTMB(maxdep ~  
                 species +  # 'known'  
                 island +   # CH1 
                 group +    # CH1
                 species:group +    # depth diff between CS and GT vary by INC/BRO/CRE
                 island:group +     # depth diff between islands vary by INC/BRO/CRE
                 (1|ids) +
                 ou(times + 0|ids),
               data= g,
               family = Gamma(link = "log"),
               control = glmmTMBControl(parallel = nt),
               REML = F)


mod2 = glmmTMB(maxdep ~  
                 species +  # 'known'  
                 island +   # CH1
                 species:island +   # test: do the species differ by the same amount between the two sites  
                 group +    # CH1
                 species:group +    # depth diff between CS and GT vary by INC/BRO/CRE
                 island:group +
                 (1|ids) +
                 ou(times + 0|ids),
               data= g,
               family = Gamma(link = "log"),
               control = glmmTMBControl(parallel = nt),
               REML = F)

mod3 = glmmTMB(maxdep ~  
                 species +  # 'known'  
                 island +   # CH1s  
                 group +    # CH1
                 island:group +
                 (1|ids) +
                 ou(times + 0|ids),
               data= g,
               family = Gamma(link = "log"),
               control = glmmTMBControl(parallel = nt),
               REML = F)

mod4 = glmmTMB(maxdep ~  
                 species +  # 'known'  
                 island +   # CH1
                 species:island +   # test: do the species differ by the same amount between the two sites  
                 group +    # CH1
                 island:group +
                 (1|ids) +
                 ou(times + 0|ids),
               data= g,
               family = Gamma(link = "log"),
               control = glmmTMBControl(parallel = nt),
               REML = F)

mod5 = glmmTMB(maxdep ~  
                 species +  # 'known'  
                 island +   # CH1 
                 group +    # CH1
                 species:group +
                 (1|ids) +
                 ou(times + 0|ids),
               data= g,
               family = Gamma(link = "log"),
               control = glmmTMBControl(parallel = nt),
               REML = F)

mod6 = glmmTMB(maxdep ~  
                 species +  # 'known'  
                 island +   # CH1
                 species:island +   # test: do the species differ by the same amount between the two sites  
                 group +    # CH1
                 species:group + 
                 (1|ids) +
                 ou(times + 0|ids),
               data= g,
               family = Gamma(link = "log"),
               control = glmmTMBControl(parallel = nt),
               REML = F)

mod7 = glmmTMB(maxdep ~  
                 species +  # 'known'  
                 island +   # CH1  
                 group +    # CH1
                 (1|ids) +
                 ou(times + 0|ids),
               data= g,
               family = Gamma(link = "log"),
               control = glmmTMBControl(parallel = nt),
               REML = F)

mod8 = glmmTMB(maxdep ~  
                 species +  # 'known'  
                 island +   # CH1
                 species:island +   # test: do the species differ by the same amount between the two sites  
                 group +    # CH1
                 (1|ids) +
                 ou(times + 0|ids),
               data= g,
               family = Gamma(link = "log"),
               control = glmmTMBControl(parallel = nt),
               REML = F)

###------------------------------------------------------------------------------------------------------------------

# Group is also proxy for: 1) space = where they are (bathymetry can be an indication, Inc - deeper waters, Bro & Cre - on the shelf)
#                       2) Solar availability (more solar light during Dec + almost no twilight vs Feb- more twilight) - 
#                        we would expect CS to dive shallower more often during Cre
#                       3) SST (environmental conditions over the season)

# Can these environmental factors better explain the pattern in vertical overlap that exists compared to the group effect? 

# Model set 2: With environmental factors
# replace group with other environmental factors everywhere; 
# (group and island:group, species:group),  but no island*bdepth interaction 

# start with a full a priori model: 
mod9 = glmmTMB(maxdep ~  
                 species +  # 'known'  
                 island +   # CH1
                 z.solar.ele.begdesc +
                 species:z.solar.ele.begdesc + # test: do the species differ in how sun elev affect maxdep
                 island:z.solar.ele.begdesc+
                 z.sst +  
                 species:z.sst +  # test: do the species differ in how sst affect maxdep
                 island:z.sst+
                 splines::ns(z.bdepth,4)+  
                 species:splines::ns(z.bdepth,4) +  
                 (1|ids) +
                 ou(times + 0|ids),
               data= g,
               control = glmmTMBControl(parallel = nt),
               REML = F)


mod10 = glmmTMB(maxdep ~  
                  species +  # 'known'  
                  island +   # CH1
                  z.solar.ele.begdesc +
                  species:z.solar.ele.begdesc + # test: do the species differ in how sun elev affect maxdep
                  island:z.solar.ele.begdesc+
                  z.sst +  
                  island:z.sst+
                  splines::ns(z.bdepth,4)+  
                  species:splines::ns(z.bdepth,4) +  
                  (1|ids) +
                  ou(times + 0|ids),
                data= g,
                control = glmmTMBControl(parallel = nt),
                REML = F)

mod11 = glmmTMB(maxdep ~  
                  species +  # 'known'  
                  island +   # CH1
                  z.solar.ele.begdesc +
                  species:z.solar.ele.begdesc + # test: do the species differ in how sun elev affect maxdep
                  island:z.solar.ele.begdesc+
                  z.sst +  
                  species:z.sst + 
                  splines::ns(z.bdepth,4)+  
                  species:splines::ns(z.bdepth,4) +  
                  (1|ids) +
                  ou(times + 0|ids),
                data= g,
                control = glmmTMBControl(parallel = nt),
                REML = F)

mod12 = glmmTMB(maxdep ~  
                  species +  # 'known'  
                  island +   # CH1
                  z.solar.ele.begdesc +
                  species:z.solar.ele.begdesc + # test: do the species differ in how sun elev affect maxdep
                  island:z.solar.ele.begdesc+
                  z.sst + 
                  splines::ns(z.bdepth,4)+  
                  species:splines::ns(z.bdepth,4) +  
                  (1|ids) +
                  ou(times + 0|ids),
                data= g,
                control = glmmTMBControl(parallel = nt),
                REML = F)

mod13 = glmmTMB(maxdep ~  
                  species +  # 'known'  
                  island +   # CH1
                  z.solar.ele.begdesc +
                  species:z.solar.ele.begdesc + # test: do the species differ in how sun elev affect maxdep
                  island:z.solar.ele.begdesc+
                  splines::ns(z.bdepth,4)+  
                  species:splines::ns(z.bdepth,4) +  
                  (1|ids) +
                  ou(times + 0|ids),
                data= g,
                control = glmmTMBControl(parallel = nt),
                REML = F)

mod14 = glmmTMB(maxdep ~  
                  species +  # 'known'  
                  island +   # CH1
                  z.solar.ele.begdesc +
                  species:z.solar.ele.begdesc + 
                  z.sst +  
                  species:z.sst +  # test: do the species differ in how sst affect maxdep
                  island:z.sst+
                  splines::ns(z.bdepth,4)+  
                  species:splines::ns(z.bdepth,4) +  
                  (1|ids) +
                  ou(times + 0|ids),
                data= g,
                control = glmmTMBControl(parallel = nt),
                REML = F)

mod15 = glmmTMB(maxdep ~  
                  species +  # 'known'  
                  island +   # CH1
                  z.solar.ele.begdesc +
                  species:z.solar.ele.begdesc + 
                  z.sst + 
                  island:z.sst+
                  splines::ns(z.bdepth,4)+  
                  species:splines::ns(z.bdepth,4) +  
                  (1|ids) +
                  ou(times + 0|ids),
                data= g,
                control = glmmTMBControl(parallel = nt),
                REML = F)


# compare models with different predictors
AIC(mod1, mod2, mod3, mod4, mod5, mod6, mod7, mod8, mod9, mod10, mod11, mod12, mod13, mod14, mod15) 

models <- c('mod1', 'mod2', 'mod3', 'mod4', 'mod5', 'mod6', 'mod7', 'mod8', 'mod9', 'mod10', 'mod11', 'mod12', 'mod13', 'mod14', 'mod15')
aics <- AIC(mod1, mod2, mod3, mod4, mod5, mod6, mod7, mod8, mod9, mod10, mod11, mod12, mod13, mod14, mod15)
delta.aics <- aics$AIC - min(aics$AIC) # To work out the change in delta AIC #see definitions in book (from Anderson (2008))
exp.delta <- exp(-0.5*delta.aics)
wi <- exp.delta/sum(exp.delta)# these are the Akaike weights for each model #The probability that model is the actual (fitted) k-l best model in the set (Anderson 2008) (See Burnham et al (2011) paper) 
(modtable <- data.frame(models, numpar=aics$df, aics$AIC, delta.aics, wi))

## write to a CSV file

write.csv(modtable, "./outputs/interspecific lmms/model selection output_breed and environmental for chapter.csv")

######################################################################
# Breeding stage effects
#####################################################################


# Plot group effects: 

group.mod = glmmTMB(maxdep ~  
                      species +  # 'known'  
                      island +   # CH1
                      # drop1 # species:island +   # test: do the species differ by the same amount between the two sites  
                      group +    # CH1
                      species:group +    # depth diff between CS and GT vary by INC/BRO/CRE
                      island:group +
                      (1|ids) +
                      ou(times + 0|ids),
                    data= g,
                    family = Gamma(link = "log"),
                    control = glmmTMBControl(parallel = nt),
                    REML = T)

summary(group.mod)

# Conditional model:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)             3.59672    0.05515   65.22  < 2e-16 ***
#   speciesgentoo           0.41731    0.05727    7.29 3.18e-13 ***
#   islandNelson           -0.36670    0.06030   -6.08 1.19e-09 ***
#   groupBro                0.23255    0.06617    3.51 0.000441 ***
#   groupCre                0.11467    0.08366    1.37 0.170457    
# speciesgentoo:groupBro -0.23339    0.07451   -3.13 0.001735 ** 
#   speciesgentoo:groupCre -0.17064    0.08089   -2.11 0.034884 *  
#   islandNelson:groupBro   0.08778    0.07334    1.20 0.231308    
# islandNelson:groupCre   0.37713    0.09075    4.16 3.24e-05 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# > 

plot(allEffects(group.mod))
# save image 
rstudioapi::savePlotAsImage(paste0("./outputs/interspecific lmms/plots/vertical overlap_group alleffects.png"), width=900,height=600)


performance::r2(group.mod) #R2c = 0.39, R2m = 0.29

tab_model(group.mod, CSS = css_theme("cells"), show.re.var = F)
acf(resid(group.mod))

#------------------------------------------
# check model performance: 
library(DHARMa)

res = simulateResiduals(group.mod)
plot(res, asFactor = T)


# 8.  Summary of model output in an HTML table

#output for one model
tab_model(group.mod, CSS = css_theme("cells"), 
          file = paste0('./outputs/interspecific lmms/glmmTMB_vertical overlap_group effects model summary.doc'))

# don't show all the random effects
tab_model(group.mod, CSS = css_theme("cells"), show.re.var = F)

# variance = residual variance
# tau00 ID = Variance from ID as random effect intercept
# ICC (Inter-correlation coefficient) = repeatability of individuals; captures the within-class similarity of the covariate adjusted data values
# N ID = Number of individuals
# Observations: Number of max dive observations
# Marginal R-squared: Variance of fixed effects
# Conditional R-squared: Variance of fixed and random effects


# Customizing HTML tables - Daniel L?decke
# see https://cran.r-project.org/web/packages/sjPlot/vignettes/table_css.html


####--------------------------------------------------------------------------------------------------------
# 2. Plotting group effects: 

summary(group.mod)

p = ggpredict(group.mod, terms = c('species','group', 'island'))
p
plot(p)


#observed data
ggplot() +
  geom_point(g, mapping = aes(group, maxdep, col = species), 
             position = position_dodge(.5), shape =1, alpha = 0.5)+ 
  facet_wrap(island~.)+
  scale_color_manual(values=met.brewer("Java", 2, direction = 1), name = 'Species')

pred <- as.data.frame(p)
head(pred)

pred = pred %>%
  dplyr::rename(group = x, 
                maxdep = predicted, 
                species= group, 
                island = facet)

# specify facet plot order - predicted
pred$species_f = factor(pred$species, levels=c('gentoo','chinstrap'))
levels(pred$species_f) <- list(gentoo = 'gentoo', chinstrap = 'chinstrap') 
pred$species_f

# specify facet plot order - observed
g$species_f = factor(g$species, levels=c('gentoo','chinstrap'))
levels(g$species_f) <- list(gentoo = 'gentoo', chinstrap = 'chinstrap') 
g$species_f

# specify facet plot order - predicted
pred$group_f = factor(pred$group, levels=c('Inc','Bro', 'Cre'))
levels(pred$group_f) <- list(Incubation = 'Inc', Brood = 'Bro', Creche = 'Cre') 
pred$group_f

# specify facet plot order - observed
g$group_f = factor(g$group, levels=c('Inc','Bro', 'Cre'))
levels(g$group_f) <- list(Incubation = 'Inc', Brood = 'Bro', Creche = 'Cre') 
g$group_f

head(pred)

# violin plots? 
s = ggplot(data =g, mapping = aes(group_f, maxdep, col =species_f))+
  geom_point(data = g, position = position_jitterdodge(dodge.width = .9, jitter.width = .3), 
             shape = 20, alpha = 0.1)+
  facet_wrap(island~.)+
  geom_violin(trim = T, alpha = 0.2, linewidth = 0.8) +  
  # coord_cartesian(ylim= c(0.0, 0.5))+ 
  geom_point(data= pred, aes(group_f, maxdep, group = species_f), #, col=species_f
             size = 2, alpha = 1, position = position_dodge(.9),col = 'grey20')+#
  geom_errorbar(data= pred, aes(group_f, maxdep, ymin = conf.low, ymax = conf.high, group = species_f), #, col = species_f
                position = position_dodge(.9),  alpha = 0.2, linewidth = 1.2, col = 'grey20')+#
  scale_color_manual(values=met.brewer("Java", 2, direction = -1), name = 'Species')+
  scale_y_reverse()
#scale_colour_manual(values = c('day'= 'orange', 'twilight' = 'purple'))
#geom_point(size = 1.5, alpha = 3, position = position_dodge(.9))+
#geom_boxplot(width=0.1, position = position_dodge(.9), size = 0.1)
#stat_summary(fun.data=mean_sdl, mult=1, 
#              geom="pointrange", color="red")
s

s <- s + theme_bw() + 
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black"),
        strip.text.x = element_text(size = 14), 
        strip.text.y = element_text(size = 14), 
        #panel.grid.major.x = element_text(size = 14), 
        # panel.grid.major.y = element_text(size = 14),
        legend.text = element_text(size = 14), 
        legend.title = element_text(size = 16), 
        legend.position = 'bottom', 
        axis.text.x = element_text(size = 12, colour = "black"), 
        axis.text.y = element_text(size = 12, colour = "black"), 
        axis.title = element_text(size = 16))+
  #ggtitle("Gentoo penguins")+
  xlab("Breeding stage") +
  ylab("Maximum depth (m)")
s

# save image 
rstudioapi::savePlotAsImage(paste0("./outputs/interspecific lmms/plots/vertical overlap_group effects.png"), width=900,height=600)

######################################################################
# Environmental covariate effects
#####################################################################
###------------------------------------------------------------------------------------------------------------------

# Group is also proxy for: 1) space = where they are (bathymetry can be an indication, Inc - deeper waters, Bro & Cre - on the shelf)
#                       2) Solar availability (more solar light during Dec + almost no twilight vs Feb- more twilight) - 
#                        we would expect CS to dive shallower more often during Cre
#                       3) SST (environmental conditions over the season)

# Can these environmental factors better explain the pattern in vertical overlap that exists compared to the group effect? 

# Model set 2: With environmental factors
# replace group with other environmental factors everywhere; 
# (group and island:group, species:group),  but no island*bdepth interaction 

full.mod = glmmTMB(maxdep ~  
                     species +  # 'known'  
                     island +   # CH1
                     ### species:island +   # test: do the species differ by the same amount between the two sites  
                     #group +    # CH1
                     #species:group +    # depth diff between CS and GT vary by INC/BRO/CRE
                     #island:group +
                     z.solar.ele.begdesc +
                     species:z.solar.ele.begdesc + # test: do the species differ in how sun elev affect maxdep
                     island:z.solar.ele.begdesc+
                     z.sst +  
                     species:z.sst +  # test: do the species differ in how sst affect maxdep
                     island:z.sst+
                     splines::ns(z.bdepth,4)+  
                     species:splines::ns(z.bdepth,4) +  # test: do the species differ in how bath affect maxdep
                     ###  island:splines::ns(z.bdepth,4) +
                     #znum.day +
                     (1|ids) +
                     ou(times + 0|ids),
                   data= g,
                   control = glmmTMBControl(parallel = nt),
                   REML = T)

summary(full.mod)
plot(allEffects(full.mod))

# save image 
rstudioapi::savePlotAsImage(paste0("./outputs/interspecific lmms/plots/vertical overlap_environmental alleffects.png"), width=900,height=600)

performance::r2(full.mod)  #R2c = 0.69, # R2m = 0.56 

tab_model(group.mod, CSS = css_theme("cells"), show.re.var = F)
acf(resid(group.mod))

#------------------------------------------
# check model performance: 
library(DHARMa)

res = simulateResiduals(full.mod)
plot(res, asFactor = T)


# 8.  Summary of model output in an HTML table

#output for one model
tab_model(full.mod, CSS = css_theme("cells"), 
          file = paste0('./outputs/interspecific lmms/glmmTMB_vertical overlap_environmental effects model summary.doc'))

# don't show all the random effects
tab_model(full.mod, CSS = css_theme("cells"), show.re.var = F)

# variance = residual variance
# tau00 ID = Variance from ID as random effect intercept
# ICC (Inter-correlation coefficient) = repeatability of individuals; captures the within-class similarity of the covariate adjusted data values
# N ID = Number of individuals
# Observations: Number of max dive observations
# Marginal R-squared: Variance of fixed effects
# Conditional R-squared: Variance of fixed and random effects


# Customizing HTML tables - Daniel L?decke
# see https://cran.r-project.org/web/packages/sjPlot/vignettes/table_css.html



### ----------------------------------
# 3. Plotting environmental effects with full model: 

#---
# 3.1 - Bathymetric depth
#---
min(g$z.bdepth)

p = ggpredict(full.mod, terms = c('z.bdepth[all]','species', 'island')) # Takes very long to run
p
plot(p)

pred <- as.data.frame(p)
head(pred)

pred = pred %>%
  dplyr::rename(bdepth = x, 
                maxdep = predicted, 
                species = group, 
                island = facet)

head(pred)

ggplot() + 
  geom_point(data = g, aes(x = z.bdepth, y = maxdep), shape = 20, alpha = 0.01, col = 'grey20')+
  facet_wrap(island~.)+
  coord_cartesian(xlim = c(NA, 1.5))+ scale_y_reverse()+
  geom_line(data = pred, aes(x=bdepth, y=maxdep, col = species), linewidth = 1)+
  geom_ribbon(pred, mapping = aes(x=bdepth, y=maxdep, ymin = conf.low, ymax = conf.high, col = species),
              alpha = 0.1)+
  scale_color_manual(values=met.brewer("Java", 2, direction =-1), name = 'Species')+
  # scale_x_continuous(breaks = seq(0,2000, 100))+
  theme_bw() + 
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black"),
        strip.text.x = element_text(size = 14), 
        strip.text.y = element_text(size = 14), 
        #panel.grid.major.x = element_text(size = 14), 
        # panel.grid.major.y = element_text(size = 14),
        legend.text = element_text(size = 14), 
        legend.title = element_text(size = 16),
        legend.position = 'bottom',
        axis.text.x = element_text(size = 12, colour = "black", angle = 90), 
        axis.text.y = element_text(size = 12, colour = "black"), 
        axis.title = element_text(size = 16))+ 
  xlab("Scaled bathymetric depth (m)") +
  ylab("Maximum depth (m)") 

# save image 
rstudioapi::savePlotAsImage(paste0("./outputs/interspecific lmms/plots/vertical overlap_bathy effects.png"), width=900,height=600)
# model predictions look better, but confidence intervals are still large for gentoos, because data beyond 1000m is unavailable. 

# what if you restrict plotting to +_ 1.5
q = ggpredict(full.mod, terms = c('z.bdepth[-0.5:1.5]','species', 'island')) # Takes very long to run
q
plot(q)


# 3.2 - Solar elevation
min(g$z.solar.ele.begdesc)

s = ggpredict(full.mod, terms = c('z.solar.ele.begdesc[all]','species', 'island')) # Takes a long time to run 
s
plot(s)


pred <- as.data.frame(s)
head(pred)

pred = pred %>%
  dplyr::rename(solar.ele.begdesc= x, 
                maxdep = predicted, 
                species = group, 
                island = facet)

head(pred)

ggplot() + 
  geom_point(data = g, aes(x = z.solar.ele.begdesc, y = maxdep), shape = '.', alpha = 0.1, col = 'grey20')+
  facet_wrap(island~.)+
  geom_line(data = pred, aes(x=solar.ele.begdesc, y=maxdep, col = species), linewidth = 1)+
  geom_ribbon(pred, mapping = aes(x=solar.ele.begdesc, y=maxdep, ymin = conf.low, ymax = conf.high, col = species),
              alpha = 0.1)+
  scale_color_manual(values=met.brewer("Java", 2, direction =-1), name = 'Species')+
  scale_y_reverse()+
  xlab("Solar elevation (degrees above the horizon)") +
  ylab("Maximum depth (m)")+ 
  # scale_x_continuous(breaks = seq(0,2000, 100))+
  theme_classic()+
  theme(legend.position = "bottom", 
        axis.text.x = element_text(angle = 90))

# specify facet plot order - predicted
pred$species_f = factor(pred$species, levels=c('gentoo','chinstrap'))
levels(pred$species_f) <- list(gentoo = 'gentoo', chinstrap = 'chinstrap') 
pred$species_f

# specify facet plot order - observed
g$species_f = factor(g$species, levels=c('gentoo','chinstrap'))
levels(g$species_f) <- list(gentoo = 'gentoo', chinstrap = 'chinstrap') 
g$species_f

head(pred)

# instead of violin plots, do side panels: 
# https://www.business-science.io/code-tools/2021/05/18/marginal_distributions.html?utm_content=buffer9350f&utm_medium=social&utm_source=twitter.com&utm_campaign=buffer
library(ggside)

s = ggplot(data=g, mapping = aes(z.solar.ele.begdesc, maxdep, col =species_f))+
  geom_point(data = g, #position = position_jitterdodge(dodge.width = .9), 
             shape = 20, alpha = 0.02)+ facet_wrap(island~.)+
  # coord_cartesian(ylim = c(120, NA), xlim = c(-12, 55))+
  geom_line(data= pred, aes(solar.ele.begdesc, maxdep, col=species_f), 
            size = 1, alpha = 1)+#,col = 'grey20'
  geom_ribbon(data= pred, aes(solar.ele.begdesc, maxdep, ymin = conf.low, ymax = conf.high, col = species_f), #
              alpha = 0.2, size = 1)+ 
  scale_y_reverse()+
  # add in the side panel densities;
  geom_xsidedensity(data = g,  aes(y = after_stat(density), # adds a side density plot (top panel)
                                   fill = species_f), 
                    alpha = 0.4) +              # position = 'stack'
  geom_ysidedensity(data = g,  aes(x = after_stat(density),   # adds a side density plot (right panel)
                                   fill = species_f), 
                    alpha = 0.4)+
  scale_color_manual(values=met.brewer("Java", 2, direction = -1), name = 'Species')+
  scale_fill_manual(values=met.brewer("Java", 2, direction = -1), name = 'Species') +
  scale_ysidex_continuous(n.breaks = 2)
s

s <- s + theme_bw() + 
  theme(ggside.panel.scale.x = 0.2,             # increase size of side panels
        ggside.panel.scale.y = 0.2,
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black"),
        strip.text.x = element_text(size = 14), 
        strip.text.y = element_text(size = 14), 
        #panel.grid.major.x = element_text(size = 14), 
        # panel.grid.major.y = element_text(size = 14),
        legend.text = element_text(size = 14), 
        legend.title = element_text(size = 16),
        legend.position = 'bottom',
        axis.text.x = element_text(size = 12, colour = "black"), 
        axis.text.y = element_text(size = 12, colour = "black"), 
        axis.title = element_text(size = 16))+
  labs(title = "Solar elevation effects on maximum depths of species by island", 
       x = 'Scaled solar elevation (degrees above the horizon)', 
       y = 'Maximum depth (m)')
s

# save image 
rstudioapi::savePlotAsImage(paste0("./outputs/interspecific lmms/plots/vertical overlap_solar elevation effects.png"), width=900,height=600)

#------------
# 3.3 - SST 
#--------------
t = ggpredict(full.mod, terms = c('z.sst[all]','species', 'island'))
t
plot(t)


pred <- as.data.frame(t)
head(pred)

pred = pred %>%
  dplyr::rename(sst= x, 
                maxdep = predicted, 
                species = group, 
                island = facet)

head(pred)

s = ggplot(data =pred, aes(sst, maxdep, col = species)) +
  geom_line(size = 1, alpha = 1) +  
  #stat_smooth(method = "lm", se = FALSE)
  geom_ribbon(
    aes(ymin = conf.low, ymax = conf.high), 
    alpha = 0.3, linewidth = 1.0
  ) +
  scale_y_reverse()+ facet_wrap(island~.)+
  scale_color_manual(values=met.brewer("Java", 2, direction = -1), name = 'Species')
#geom_point(trip, mapping = aes(island, trip.duration, col = species), position = position_dodge(.5), 
#shape = 1) #+ facet_wrap(~group, ncol =3)
s

# try adding observed points to the predicted estimates
s <- s + geom_point(g, mapping = aes(x = sst, y = maxdep), 
                    shape = 20, alpha = 0.01, col = 'grey20')+
  coord_cartesian(ylim= c(150, NA))
s  

s <- s + theme_bw() + 
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black"),
        strip.text.x = element_text(size = 14), 
        strip.text.y = element_text(size = 14), 
        #panel.grid.major.x = element_text(size = 14), 
        # panel.grid.major.y = element_text(size = 14),
        legend.text = element_text(size = 14), 
        legend.title = element_text(size = 16), 
        axis.text.x = element_text(size = 12, colour = "black"), 
        axis.text.y = element_text(size = 12, colour = "black"), 
        axis.title = element_text(size = 16))+
  #ggtitle("Kopaitic Island")+
  xlab("Sea surface temperature (degrees Celsius)") +
  ylab("Maximum depth (m)")
s


# specify facet plot order - predicted
pred$species_f = factor(pred$species, levels=c('gentoo','chinstrap'))
levels(pred$species_f) <- list(gentoo = 'gentoo', chinstrap = 'chinstrap') 
pred$species_f

# specify facet plot order - observed
g$species_f = factor(g$species, levels=c('gentoo','chinstrap'))
levels(g$species_f) <- list(gentoo = 'gentoo', chinstrap = 'chinstrap') 
g$species_f

# instead of violin plots, do side panels: 
# https://www.business-science.io/code-tools/2021/05/18/marginal_distributions.html?utm_content=buffer9350f&utm_medium=social&utm_source=twitter.com&utm_campaign=buffer
library(ggside)

s = ggplot(data=g, mapping = aes(sst, maxdep, col =species_f))+
  geom_point(data = g, #position = position_jitterdodge(dodge.width = .9), 
             shape = 20, alpha = 0.01)+
  # coord_cartesian(ylim = c(120, NA), xlim = c(-12, 55))+
  geom_line(data= pred, aes(sst, maxdep, col=species_f), 
            size = 1, alpha = 1)+#,col = 'grey20'
  geom_ribbon(data= pred, aes(sst, maxdep, ymin = conf.low, ymax = conf.high, col = species_f), #
              alpha = 0.2, size = 1)+ 
  scale_y_reverse()+facet_wrap(island~.)+
  # add in the side panel densities;
  geom_xsidedensity(data = g,  aes(y = after_stat(density), # adds a side density plot (top panel)
                                   fill = species_f), 
                    alpha = 0.4) +              # position = 'stack'
  geom_ysidedensity(data = g,  aes(x = after_stat(density),   # adds a side density plot (right panel)
                                   fill = species_f), 
                    alpha = 0.4)+
  scale_color_manual(values=met.brewer("Java", 2, direction = -1), name = 'Species')+
  scale_fill_manual(values=met.brewer("Java", 2, direction = -1), name = 'Species') +
  scale_ysidex_continuous(n.breaks = 2)
s

s <- s + theme_bw() + 
  theme(ggside.panel.scale.x = 0.2,             # increase size of side panels
        ggside.panel.scale.y = 0.3,
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black"),
        strip.text.x = element_text(size = 14), 
        strip.text.y = element_text(size = 14), 
        #panel.grid.major.x = element_text(size = 14), 
        # panel.grid.major.y = element_text(size = 14),
        legend.text = element_text(size = 14), 
        legend.title = element_text(size = 16),
        legend.position = 'bottom',
        axis.text.x = element_text(size = 12, colour = "black"), 
        axis.text.y = element_text(size = 12, colour = "black"), 
        axis.title = element_text(size = 16))+
  labs(x = 'Scaled sea surface temperature (degrees Celsius)', 
       y = 'Maximum depth (m)')
s

# save image 
rstudioapi::savePlotAsImage(paste0("./outputs/interspecific lmms/plots/vertical overlap_sst effects.png"), width=900,height=600)


