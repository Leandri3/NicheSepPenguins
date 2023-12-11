# NicheSepPenguins
This repository contains examples of the important R scripts that I used in my data analysis for my MSc dissertation, Chapter 3: Central-place foraging constraints and environmental conditions modify niche partitioning in sympatric Pygoscelis penguin populations. I did not include scripts where I plotted the figures.

1) GPS tracks to trips
   -Divide GPS track into separate trips to sea based on diving data
   -Calculate maximum and cumulative trip distances and trip durations based on GPS

2) DiveMove Analysis
   - Zero-offset correction
   - Calculate dive statistics for all individuals in a deployment round

3) CrawlWrap for dive locations
   -speed filter GPS data (remove unreliable locations) with trip::sda
   -crawlWrap dive locations of multiple animals with a function (momentuHMM)
   
4) CrawlWrap - 5 min
   -crawlWrap 5 minute locations of multiple animals with a function (momentuHMM)

5) Autocorrelation kernel density estimation method
  - Fit ctmm utilization distribution to GPS data.
  - Exclude land areas from UDs using a land mask
  - Use 5 - min crawl GPS data
  - Estimate overlap using Bhattacharyya (BA) index

6) EM classification of 5m dives
  - Estimates dive residuals (using linear mixed-effects models)
  - Applies an Expectation-Maximization algorithm to Classify diving types into 2 clusters (foraging and non-foraging) per species, per island
  - Using the following predictors: Bottom time, Dive residuals, Maximum dive depth
  - Separates into three different Clusters: 1 = Foraging

7) Extract solar elevation for dive locations
  - Calculate solar elevation for every foraging dive location

8) Extract environmental covariates for 95% UDs
  - Imports remote-sensed environmental variables from various sources
  - Extracts remote-sensed environmental information for each dive location

9) Dive behaviours with generalised linear-mixed effects models
For maximum dive depth behaviours
  - Do model selection using - glmmTMB
  - Includes autocorrelation structures
  - family = Gamma(link = "log"),
  - model set 1 = breeding stage (group)
  - model set 2 = environmental effects
  - Predicts effects from best model sets
  - Plots figures of model output
