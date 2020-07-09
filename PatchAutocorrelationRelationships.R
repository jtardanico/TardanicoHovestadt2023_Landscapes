#================================================
# Relationships with autocorrelation
#================================================
# Statistics and organism traits broken down by
# scenario and autocorrelation degree.

# Autocorrelation values: Uniform, H=0, H=0.5, H=1
# Values apply to patch temperature and habitat type

#=============================================
# Scenario key
#=============================================
# Four climate scenarios
# 0: static climate
# 1: Static mean with random variation, st.dev. = 1
# 2: Rising mean w/ no random variation
# 3: Rising mean w/ random variation, st.dev. = 1
#=============================================
# Relevant parameters and details
#=============================================

#=============================================
# Main script
#=============================================

# Clear workspace
rm(list=ls())

#-------------------------------------------

library(ggplot2)
library(readr)
library(lattice)
library(tidyverse)

patches <- read.csv("~/Documents/PhD/LandklifModel/Landscapes/Model1/Experiment1/G1Patches/patchstats_h1_test12.txt")
patches$H_t <- as.factor(patches$H_t)
patches$H_h <- as.factor(patches$H_h)

levels(patches$H_t)
levels(patches$H_h)

setwd("Documents/PhD/LandklifModel/Landscapes/Model1/Experiment1/test10/plots1")


#============================================================
# Whole landscape density plots by trait and climate scenario
#============================================================

# Subsetting

i <- subset(patches,patches$Timestep==-1) # i for 'initial'
b <- subset(patches,patches$Timestep==0) # b for 'burnt in'
e <- subset(patches,patches$Timestep==50) # e for 'end'


  scen0i <- subset(i, i$clim_scen==0)
  scen1i <- subset(i, i$clim_scen==1)
  scen2i <- subset(i, i$clim_scen==2)
  scen3i <- subset(i, i$clim_scen==3)
  
  scen0b <- subset(b, b$clim_scen==0)
  scen1b <- subset(b, b$clim_scen==1)
  scen2b <- subset(b, b$clim_scen==2)
  scen3b <- subset(b, b$clim_scen==3)
  
  scen0e <- subset(e, e$clim_scen==0)
  scen1e <- subset(e, e$clim_scen==1)
  scen2e <- subset(e, e$clim_scen==2)
  scen3e <- subset(e, e$clim_scen==3)



# Population size

ggplot(scen0i,aes(pop)) + geom_histogram(binwidth = 1) + facet_grid(H_t~H_h,labeller = label_both) #+ xlim(0,370) + ylim(0,100)
ggplot(scen1i,aes(pop)) + geom_histogram(binwidth = 1) + facet_grid(H_t~H_h,labeller = label_both) #+ xlim(0,370) + ylim(0,100)
ggplot(scen2i,aes(pop)) + geom_histogram(binwidth = 1) + facet_grid(H_t~H_h,labeller = label_both) #+ xlim(0,370) + ylim(0,100)
ggplot(scen3i,aes(pop)) + geom_histogram(binwidth = 1) + facet_grid(H_t~H_h,labeller = label_both) #+ xlim(0,370) + ylim(0,100)

ggplot(scen0b,aes(pop)) + geom_histogram(binwidth = 1) + facet_grid(H_t~H_h,labeller = label_both) + xlim(0,380) #+ ylim(0,100)
ggplot(scen1b,aes(pop)) + geom_histogram(binwidth = 1) + facet_grid(H_t~H_h,labeller = label_both) + xlim(0,380) #+ ylim(0,100)
ggplot(scen2b,aes(pop)) + geom_histogram(binwidth = 1) + facet_grid(H_t~H_h,labeller = label_both) + xlim(0,380) #+ ylim(0,100)
ggplot(scen3b,aes(pop)) + geom_histogram(binwidth = 1) + facet_grid(H_t~H_h,labeller = label_both) + xlim(0,380) #+ ylim(0,100)

ggplot(scen0e,aes(pop)) + geom_histogram(binwidth = 1) + facet_grid(H_t~H_h,labeller = label_both) + xlim(0,380) #+ ylim(0,100)
ggplot(scen1e,aes(pop)) + geom_histogram(binwidth = 1) + facet_grid(H_t~H_h,labeller = label_both) + xlim(0,380) #+ ylim(0,100)
ggplot(scen2e,aes(pop)) + geom_histogram(binwidth = 1) + facet_grid(H_t~H_h,labeller = label_both) + xlim(0,380) #+ ylim(0,100)
ggplot(scen3e,aes(pop)) + geom_histogram(binwidth = 1) + facet_grid(H_t~H_h,labeller = label_both) + xlim(0,380) #+ ylim(0,100)

# No. of lineages
ggplot(scen0i,aes(richness)) + geom_histogram() + facet_grid(H_t~H_h,labeller = label_both)
ggplot(scen1i,aes(richness)) + geom_histogram() + facet_grid(H_t~H_h,labeller = label_both)
ggplot(scen2i,aes(richness)) + geom_histogram() + facet_grid(H_t~H_h,labeller = label_both)
ggplot(scen3i,aes(richness)) + geom_histogram() + facet_grid(H_t~H_h,labeller = label_both)

ggplot(scen0b,aes(richness)) + geom_histogram() + facet_grid(H_t~H_h,labeller = label_both) + xlim(0,15)
ggplot(scen1b,aes(richness)) + geom_histogram() + facet_grid(H_t~H_h,labeller = label_both) + xlim(0,15)
ggplot(scen2b,aes(richness)) + geom_histogram() + facet_grid(H_t~H_h,labeller = label_both) + xlim(0,15)
ggplot(scen3b,aes(richness)) + geom_histogram() + facet_grid(H_t~H_h,labeller = label_both) + xlim(0,15)

ggplot(scen0e,aes(richness)) + geom_histogram() + facet_grid(H_t~H_h,labeller = label_both) + xlim(0,15)
ggplot(scen1e,aes(richness)) + geom_histogram() + facet_grid(H_t~H_h,labeller = label_both) + xlim(0,15)
ggplot(scen2e,aes(richness)) + geom_histogram() + facet_grid(H_t~H_h,labeller = label_both) + xlim(0,15)
ggplot(scen3e,aes(richness)) + geom_histogram() + facet_grid(H_t~H_h,labeller = label_both) + xlim(0,15)

# Shannon index
ggplot(scen0i,aes(shannon)) + geom_histogram() + facet_grid(H_t~H_h,labeller = label_both)
ggplot(scen1i,aes(shannon)) + geom_histogram() + facet_grid(H_t~H_h,labeller = label_both)
ggplot(scen2i,aes(shannon)) + geom_histogram() + facet_grid(H_t~H_h,labeller = label_both)
ggplot(scen3i,aes(shannon)) + geom_histogram() + facet_grid(H_t~H_h,labeller = label_both)

ggplot(scen0b,aes(shannon)) + geom_histogram() + facet_grid(H_t~H_h,labeller = label_both)
ggplot(scen1b,aes(shannon)) + geom_histogram() + facet_grid(H_t~H_h,labeller = label_both)
ggplot(scen2b,aes(shannon)) + geom_histogram() + facet_grid(H_t~H_h,labeller = label_both)
ggplot(scen3b,aes(shannon)) + geom_histogram() + facet_grid(H_t~H_h,labeller = label_both)

ggplot(scen0e,aes(shannon)) + geom_histogram() + facet_grid(H_t~H_h,labeller = label_both)
ggplot(scen1e,aes(shannon)) + geom_histogram() + facet_grid(H_t~H_h,labeller = label_both)
ggplot(scen2e,aes(shannon)) + geom_histogram() + facet_grid(H_t~H_h,labeller = label_both)
ggplot(scen3e,aes(shannon)) + geom_histogram() + facet_grid(H_t~H_h,labeller = label_both)

# Simpson index
ggplot(scen0i,aes(simpson)) + geom_histogram() + facet_grid(H_t~H_h,labeller = label_both)
ggplot(scen1i,aes(simpson)) + geom_histogram() + facet_grid(H_t~H_h,labeller = label_both)
ggplot(scen2i,aes(simpson)) + geom_histogram() + facet_grid(H_t~H_h,labeller = label_both)
ggplot(scen3i,aes(simpson)) + geom_histogram() + facet_grid(H_t~H_h,labeller = label_both)

ggplot(scen0b,aes(simpson)) + geom_histogram() + facet_grid(H_t~H_h,labeller = label_both)
ggplot(scen1b,aes(simpson)) + geom_histogram() + facet_grid(H_t~H_h,labeller = label_both)
ggplot(scen2b,aes(simpson)) + geom_histogram() + facet_grid(H_t~H_h,labeller = label_both)
ggplot(scen3b,aes(simpson)) + geom_histogram() + facet_grid(H_t~H_h,labeller = label_both)

ggplot(scen0e,aes(simpson)) + geom_histogram() + facet_grid(H_t~H_h,labeller = label_both)
ggplot(scen1e,aes(simpson)) + geom_histogram() + facet_grid(H_t~H_h,labeller = label_both)
ggplot(scen2e,aes(simpson)) + geom_histogram() + facet_grid(H_t~H_h,labeller = label_both)
ggplot(scen3e,aes(simpson)) + geom_histogram() + facet_grid(H_t~H_h,labeller = label_both)

# Habitat fitness (niche match)
ggplot(scen0i,aes(fh)) + geom_histogram(binwidth = 0.01) + facet_grid(H_t~H_h,labeller = label_both) + xlim(0,1)
ggplot(scen1i,aes(fh)) + geom_histogram(binwidth = 0.01) + facet_grid(H_t~H_h,labeller = label_both) + xlim(0,1)
ggplot(scen2i,aes(fh)) + geom_histogram(binwidth = 0.01) + facet_grid(H_t~H_h,labeller = label_both) + xlim(0,1)
ggplot(scen3i,aes(fh)) + geom_histogram(binwidth = 0.01) + facet_grid(H_t~H_h,labeller = label_both) + xlim(0,1)

ggplot(scen0b,aes(fh)) + geom_histogram(binwidth = 0.01) + facet_grid(H_t~H_h,labeller = label_both) + xlim(0,1)
ggplot(scen1b,aes(fh)) + geom_histogram(binwidth = 0.01) + facet_grid(H_t~H_h,labeller = label_both) + xlim(0,1)
ggplot(scen2b,aes(fh)) + geom_histogram(binwidth = 0.01) + facet_grid(H_t~H_h,labeller = label_both) + xlim(0,1)
ggplot(scen3b,aes(fh)) + geom_histogram(binwidth = 0.01) + facet_grid(H_t~H_h,labeller = label_both) + xlim(0,1)

ggplot(scen0e,aes(fh)) + geom_histogram(binwidth = 0.01) + facet_grid(H_t~H_h,labeller = label_both) + xlim(0,1)
ggplot(scen1e,aes(fh)) + geom_histogram(binwidth = 0.01) + facet_grid(H_t~H_h,labeller = label_both) + xlim(0,1)
ggplot(scen2e,aes(fh)) + geom_histogram(binwidth = 0.01) + facet_grid(H_t~H_h,labeller = label_both) + xlim(0,1)
ggplot(scen3e,aes(fh)) + geom_histogram(binwidth = 0.01) + facet_grid(H_t~H_h,labeller = label_both) + xlim(0,1)

# Temperature fitness (niche match)
ggplot(scen0i,aes(ft)) + geom_histogram(binwidth = 0.01) + facet_grid(H_t~H_h,labeller = label_both) + xlim(0,1)
ggplot(scen1i,aes(ft)) + geom_histogram(binwidth = 0.01) + facet_grid(H_t~H_h,labeller = label_both) + xlim(0,1)
ggplot(scen2i,aes(ft)) + geom_histogram(binwidth = 0.01) + facet_grid(H_t~H_h,labeller = label_both) + xlim(0,1)
ggplot(scen3i,aes(ft)) + geom_histogram(binwidth = 0.01) + facet_grid(H_t~H_h,labeller = label_both) + xlim(0,1)

ggplot(scen0b,aes(ft)) + geom_histogram(binwidth = 0.01) + facet_grid(H_t~H_h,labeller = label_both) + xlim(0,1)
ggplot(scen1b,aes(ft)) + geom_histogram(binwidth = 0.01) + facet_grid(H_t~H_h,labeller = label_both) + xlim(0,1)
ggplot(scen2b,aes(ft)) + geom_histogram(binwidth = 0.01) + facet_grid(H_t~H_h,labeller = label_both) + xlim(0,1)
ggplot(scen3b,aes(ft)) + geom_histogram(binwidth = 0.01) + facet_grid(H_t~H_h,labeller = label_both) + xlim(0,1)

ggplot(scen0e,aes(ft)) + geom_histogram(binwidth = 0.01) + facet_grid(H_t~H_h,labeller = label_both) + xlim(0,1)
ggplot(scen1e,aes(ft)) + geom_histogram(binwidth = 0.01) + facet_grid(H_t~H_h,labeller = label_both) + xlim(0,1)
ggplot(scen2e,aes(ft)) + geom_histogram(binwidth = 0.01) + facet_grid(H_t~H_h,labeller = label_both) + xlim(0,1)
ggplot(scen3e,aes(ft)) + geom_histogram(binwidth = 0.01) + facet_grid(H_t~H_h,labeller = label_both) + xlim(0,1)

# Overall fitness 
ggplot(scen0i,aes(f)) + geom_histogram(binwidth = 0.01) + facet_grid(H_t~H_h,labeller = label_both) + xlim(0,1)
ggplot(scen1i,aes(f)) + geom_histogram(binwidth = 0.01) + facet_grid(H_t~H_h,labeller = label_both) + xlim(0,1)
ggplot(scen2i,aes(f)) + geom_histogram(binwidth = 0.01) + facet_grid(H_t~H_h,labeller = label_both) + xlim(0,1)
ggplot(scen3i,aes(f)) + geom_histogram(binwidth = 0.01) + facet_grid(H_t~H_h,labeller = label_both) + xlim(0,1)

ggplot(scen0b,aes(f)) + geom_histogram(binwidth = 0.01) + facet_grid(H_t~H_h,labeller = label_both) + xlim(0,1)
ggplot(scen1b,aes(f)) + geom_histogram(binwidth = 0.01) + facet_grid(H_t~H_h,labeller = label_both) + xlim(0,1)
ggplot(scen2b,aes(f)) + geom_histogram(binwidth = 0.01) + facet_grid(H_t~H_h,labeller = label_both) + xlim(0,1)
ggplot(scen3b,aes(f)) + geom_histogram(binwidth = 0.01) + facet_grid(H_t~H_h,labeller = label_both) + xlim(0,1)

ggplot(scen0e,aes(f)) + geom_histogram(binwidth = 0.01) + facet_grid(H_t~H_h,labeller = label_both) + xlim(0,1)
ggplot(scen1e,aes(f)) + geom_histogram(binwidth = 0.01) + facet_grid(H_t~H_h,labeller = label_both) + xlim(0,1)
ggplot(scen2e,aes(f)) + geom_histogram(binwidth = 0.01) + facet_grid(H_t~H_h,labeller = label_both) + xlim(0,1)
ggplot(scen3e,aes(f)) + geom_histogram(binwidth = 0.01) + facet_grid(H_t~H_h,labeller = label_both) + xlim(0,1)

# Mean temperature optimum
ggplot(scen0i,aes(T_opt)) + geom_histogram(binwidth = 0.05) + facet_grid(H_t~H_h,labeller = label_both) #+ xlim(10,15)
ggplot(scen1i,aes(T_opt)) + geom_histogram(binwidth = 0.05) + facet_grid(H_t~H_h,labeller = label_both) #+ xlim(10,15)
ggplot(scen2i,aes(T_opt)) + geom_histogram(binwidth = 0.05) + facet_grid(H_t~H_h,labeller = label_both) #+ xlim(10,15)
ggplot(scen3i,aes(T_opt)) + geom_histogram(binwidth = 0.05) + facet_grid(H_t~H_h,labeller = label_both) #+ xlim(10,15)

ggplot(scen0b,aes(T_opt)) + geom_histogram(binwidth = 0.05) + facet_grid(H_t~H_h,labeller = label_both) #+ xlim(10,15)
ggplot(scen1b,aes(T_opt)) + geom_histogram(binwidth = 0.05) + facet_grid(H_t~H_h,labeller = label_both) #+ xlim(10,15)
ggplot(scen2b,aes(T_opt)) + geom_histogram(binwidth = 0.05) + facet_grid(H_t~H_h,labeller = label_both) #+ xlim(10,15)
ggplot(scen3b,aes(T_opt)) + geom_histogram(binwidth = 0.05) + facet_grid(H_t~H_h,labeller = label_both) #+ xlim(10,15)

ggplot(scen0e,aes(T_opt)) + geom_histogram(binwidth = 0.05) + facet_grid(H_t~H_h,labeller = label_both) #+ xlim(10,15)
ggplot(scen1e,aes(T_opt)) + geom_histogram(binwidth = 0.05) + facet_grid(H_t~H_h,labeller = label_both) #+ xlim(10,15)
ggplot(scen2e,aes(T_opt)) + geom_histogram(binwidth = 0.05) + facet_grid(H_t~H_h,labeller = label_both) #+ xlim(10,15)
ggplot(scen3e,aes(T_opt)) + geom_histogram(binwidth = 0.05) + facet_grid(H_t~H_h,labeller = label_both) #+ xlim(10,15)

# Mean temperature tolerance
ggplot(scen0i,aes(T_sd)) + geom_histogram(binwidth = 0.05) + facet_grid(H_t~H_h,labeller = label_both) #+ xlim(0,3)
ggplot(scen1i,aes(T_sd)) + geom_histogram(binwidth = 0.05) + facet_grid(H_t~H_h,labeller = label_both) #+ xlim(0,3)
ggplot(scen2i,aes(T_sd)) + geom_histogram(binwidth = 0.05) + facet_grid(H_t~H_h,labeller = label_both) #+ xlim(0,3)
ggplot(scen3i,aes(T_sd)) + geom_histogram(binwidth = 0.05) + facet_grid(H_t~H_h,labeller = label_both) #+ xlim(0,3)

ggplot(scen0b,aes(T_sd)) + geom_histogram(binwidth = 0.05) + facet_grid(H_t~H_h,labeller = label_both) + xlim(0,5)
ggplot(scen1b,aes(T_sd)) + geom_histogram(binwidth = 0.05) + facet_grid(H_t~H_h,labeller = label_both) + xlim(0,5)
ggplot(scen2b,aes(T_sd)) + geom_histogram(binwidth = 0.05) + facet_grid(H_t~H_h,labeller = label_both) + xlim(0,5)
ggplot(scen3b,aes(T_sd)) + geom_histogram(binwidth = 0.05) + facet_grid(H_t~H_h,labeller = label_both) + xlim(0,5)

ggplot(scen0e,aes(T_sd)) + geom_histogram(binwidth = 0.05) + facet_grid(H_t~H_h,labeller = label_both) + xlim(0,5)
ggplot(scen1e,aes(T_sd)) + geom_histogram(binwidth = 0.05) + facet_grid(H_t~H_h,labeller = label_both) + xlim(0,5)
ggplot(scen2e,aes(T_sd)) + geom_histogram(binwidth = 0.05) + facet_grid(H_t~H_h,labeller = label_both) + xlim(0,5)
ggplot(scen3e,aes(T_sd)) + geom_histogram(binwidth = 0.05) + facet_grid(H_t~H_h,labeller = label_both) + xlim(0,5)


# Mean habitat optimum
ggplot(scen0i,aes(H_opt)) + geom_histogram(binwidth = 0.01) + facet_grid(H_t~H_h,labeller = label_both) + xlim(0,1)
ggplot(scen1i,aes(H_opt)) + geom_histogram(binwidth = 0.01) + facet_grid(H_t~H_h,labeller = label_both) + xlim(0,1)
ggplot(scen2i,aes(H_opt)) + geom_histogram(binwidth = 0.01) + facet_grid(H_t~H_h,labeller = label_both) + xlim(0,1)
ggplot(scen3i,aes(H_opt)) + geom_histogram(binwidth = 0.01) + facet_grid(H_t~H_h,labeller = label_both) + xlim(0,1)

ggplot(scen0b,aes(H_opt)) + geom_histogram(binwidth = 0.01) + facet_grid(H_t~H_h,labeller = label_both) + xlim(-1,2)
ggplot(scen1b,aes(H_opt)) + geom_histogram(binwidth = 0.01) + facet_grid(H_t~H_h,labeller = label_both) + xlim(-1,2)
ggplot(scen2b,aes(H_opt)) + geom_histogram(binwidth = 0.01) + facet_grid(H_t~H_h,labeller = label_both) + xlim(-1,2)
ggplot(scen3b,aes(H_opt)) + geom_histogram(binwidth = 0.01) + facet_grid(H_t~H_h,labeller = label_both) + xlim(-1,2)

ggplot(scen0e,aes(H_opt)) + geom_histogram(binwidth = 0.01) + facet_grid(H_t~H_h,labeller = label_both) + xlim(-1,2)
ggplot(scen1e,aes(H_opt)) + geom_histogram(binwidth = 0.01) + facet_grid(H_t~H_h,labeller = label_both) + xlim(-1,2)
ggplot(scen2e,aes(H_opt)) + geom_histogram(binwidth = 0.01) + facet_grid(H_t~H_h,labeller = label_both) + xlim(-1,2)
ggplot(scen3e,aes(H_opt)) + geom_histogram(binwidth = 0.01) + facet_grid(H_t~H_h,labeller = label_both) + xlim(-1,2)

# Mean habitat tolerance
ggplot(scen0i,aes(H_sd)) + geom_histogram(binwidth = 0.05) + facet_grid(H_t~H_h,labeller = label_both) + xlim(0,5)
ggplot(scen1i,aes(H_sd)) + geom_histogram(binwidth = 0.05) + facet_grid(H_t~H_h,labeller = label_both) + xlim(0,5)
ggplot(scen2i,aes(H_sd)) + geom_histogram(binwidth = 0.05) + facet_grid(H_t~H_h,labeller = label_both) + xlim(0,5)
ggplot(scen3i,aes(H_sd)) + geom_histogram(binwidth = 0.05) + facet_grid(H_t~H_h,labeller = label_both) + xlim(0,5)

ggplot(scen0b,aes(H_sd)) + geom_histogram(binwidth = 0.05) + facet_grid(H_t~H_h,labeller = label_both) + xlim(0,5)
ggplot(scen1b,aes(H_sd)) + geom_histogram(binwidth = 0.05) + facet_grid(H_t~H_h,labeller = label_both) + xlim(0,5)
ggplot(scen2b,aes(H_sd)) + geom_histogram(binwidth = 0.05) + facet_grid(H_t~H_h,labeller = label_both) + xlim(0,5)
ggplot(scen3b,aes(H_sd)) + geom_histogram(binwidth = 0.05) + facet_grid(H_t~H_h,labeller = label_both) + xlim(0,5)

ggplot(scen0e,aes(H_sd)) + geom_histogram(binwidth = 0.05) + facet_grid(H_t~H_h,labeller = label_both) + xlim(0,5)
ggplot(scen1e,aes(H_sd)) + geom_histogram(binwidth = 0.05) + facet_grid(H_t~H_h,labeller = label_both) + xlim(0,5)
ggplot(scen2e,aes(H_sd)) + geom_histogram(binwidth = 0.05) + facet_grid(H_t~H_h,labeller = label_both) + xlim(0,5)
ggplot(scen3e,aes(H_sd)) + geom_histogram(binwidth = 0.05) + facet_grid(H_t~H_h,labeller = label_both) + xlim(0,5)

# Mean dispersal probability
ggplot(scen0i,aes(disp_l)) + geom_histogram(binwidth = 0.01) + facet_grid(H_t~H_h,labeller = label_both) + xlim(0,1)
ggplot(scen1i,aes(disp_l)) + geom_histogram(binwidth = 0.01) + facet_grid(H_t~H_h,labeller = label_both) + xlim(0,1)
ggplot(scen2i,aes(disp_l)) + geom_histogram(binwidth = 0.01) + facet_grid(H_t~H_h,labeller = label_both) + xlim(0,1)
ggplot(scen3i,aes(disp_l)) + geom_histogram(binwidth = 0.01) + facet_grid(H_t~H_h,labeller = label_both) + xlim(0,1)

ggplot(scen0b,aes(disp_l)) + geom_histogram(binwidth = 0.01) + facet_grid(H_t~H_h,labeller = label_both) + xlim(0,1)
ggplot(scen1b,aes(disp_l)) + geom_histogram(binwidth = 0.01) + facet_grid(H_t~H_h,labeller = label_both) + xlim(0,1)
ggplot(scen2b,aes(disp_l)) + geom_histogram(binwidth = 0.01) + facet_grid(H_t~H_h,labeller = label_both) + xlim(0,1)
ggplot(scen3b,aes(disp_l)) + geom_histogram(binwidth = 0.01) + facet_grid(H_t~H_h,labeller = label_both) + xlim(0,1)

ggplot(scen0e,aes(disp_l)) + geom_histogram(binwidth = 0.01) + facet_grid(H_t~H_h,labeller = label_both) + xlim(0,1)
ggplot(scen1e,aes(disp_l)) + geom_histogram(binwidth = 0.01) + facet_grid(H_t~H_h,labeller = label_both) + xlim(0,1)
ggplot(scen2e,aes(disp_l)) + geom_histogram(binwidth = 0.01) + facet_grid(H_t~H_h,labeller = label_both) + xlim(0,1)
ggplot(scen3e,aes(disp_l)) + geom_histogram(binwidth = 0.01) + facet_grid(H_t~H_h,labeller = label_both) + xlim(0,1)

# Mean global dispersal probability
ggplot(scen0i,aes(disp_g)) + geom_histogram(binwidth = 0.01) + facet_grid(H_t~H_h,labeller = label_both) + xlim(0,1)
ggplot(scen1i,aes(disp_g)) + geom_histogram(binwidth = 0.01) + facet_grid(H_t~H_h,labeller = label_both) + xlim(0,1)
ggplot(scen2i,aes(disp_g)) + geom_histogram(binwidth = 0.01) + facet_grid(H_t~H_h,labeller = label_both) + xlim(0,1)
ggplot(scen3i,aes(disp_g)) + geom_histogram(binwidth = 0.01) + facet_grid(H_t~H_h,labeller = label_both) + xlim(0,1)

ggplot(scen0b,aes(disp_g)) + geom_histogram(binwidth = 0.01) + facet_grid(H_t~H_h,labeller = label_both) + xlim(0,1)
ggplot(scen1b,aes(disp_g)) + geom_histogram(binwidth = 0.01) + facet_grid(H_t~H_h,labeller = label_both) + xlim(0,1)
ggplot(scen2b,aes(disp_g)) + geom_histogram(binwidth = 0.01) + facet_grid(H_t~H_h,labeller = label_both) + xlim(0,1)
ggplot(scen3b,aes(disp_g)) + geom_histogram(binwidth = 0.01) + facet_grid(H_t~H_h,labeller = label_both) + xlim(0,1)

ggplot(scen0e,aes(disp_g)) + geom_histogram(binwidth = 0.01) + facet_grid(H_t~H_h,labeller = label_both) + xlim(0,1)
ggplot(scen1e,aes(disp_g)) + geom_histogram(binwidth = 0.01) + facet_grid(H_t~H_h,labeller = label_both) + xlim(0,1)
ggplot(scen2e,aes(disp_g)) + geom_histogram(binwidth = 0.01) + facet_grid(H_t~H_h,labeller = label_both) + xlim(0,1)
ggplot(scen3e,aes(disp_g)) + geom_histogram(binwidth = 0.01) + facet_grid(H_t~H_h,labeller = label_both) + xlim(0,1)


ggplot(patches,aes(H_sd)) + geom_histogram()
ggplot(patches,aes(T_sd)) + geom_histogram()

#=====================================================
# Landscape raster plots by trait and climate scenario
#=====================================================

a <- scale_fill_gradientn(colours = hcl.colors(100),limits = c(1,15))
b <- scale_fill_gradientn(colours = hcl.colors(100),limits = c(0,2))
c <- scale_fill_gradientn(colours = hcl.colors(100),limits = c(1,370))
d <- scale_fill_gradientn(colours = hcl.colors(100),limits = c(0,1))
e <- scale_fill_gradientn(colours = hcl.colors(100),limits = c(-1,1.8))


# Population size

ggplot(scen0i,aes(x,y)) + geom_raster(aes(fill=pop)) + facet_grid(H_t~H_h,labeller = label_both) + c
ggplot(scen1i,aes(x,y)) + geom_raster(aes(fill=pop)) + facet_grid(H_t~H_h,labeller = label_both) + c
ggplot(scen2i,aes(x,y)) + geom_raster(aes(fill=pop)) + facet_grid(H_t~H_h,labeller = label_both) + c
ggplot(scen3i,aes(x,y)) + geom_raster(aes(fill=pop)) + facet_grid(H_t~H_h,labeller = label_both) + c

ggplot(scen0b,aes(x,y)) + geom_raster(aes(fill=pop)) + facet_grid(H_t~H_h,labeller = label_both) + c
ggplot(scen1b,aes(x,y)) + geom_raster(aes(fill=pop)) + facet_grid(H_t~H_h,labeller = label_both) + c
ggplot(scen2b,aes(x,y)) + geom_raster(aes(fill=pop)) + facet_grid(H_t~H_h,labeller = label_both) + c
ggplot(scen3b,aes(x,y)) + geom_raster(aes(fill=pop)) + facet_grid(H_t~H_h,labeller = label_both) + c

ggplot(scen0e,aes(x,y)) + geom_raster(aes(fill=pop)) + facet_grid(H_t~H_h,labeller = label_both) + c
ggplot(scen1e,aes(x,y)) + geom_raster(aes(fill=pop)) + facet_grid(H_t~H_h,labeller = label_both) + c
ggplot(scen2e,aes(x,y)) + geom_raster(aes(fill=pop)) + facet_grid(H_t~H_h,labeller = label_both) + c
ggplot(scen3e,aes(x,y)) + geom_raster(aes(fill=pop)) + facet_grid(H_t~H_h,labeller = label_both) + c

# No. of Lineages

ggplot(scen0i,aes(x,y)) + geom_raster(aes(fill=richness)) + facet_grid(H_t~H_h,labeller = label_both) 
ggplot(scen1i,aes(x,y)) + geom_raster(aes(fill=richness)) + facet_grid(H_t~H_h,labeller = label_both)
ggplot(scen2i,aes(x,y)) + geom_raster(aes(fill=richness)) + facet_grid(H_t~H_h,labeller = label_both)
ggplot(scen3i,aes(x,y)) + geom_raster(aes(fill=richness)) + facet_grid(H_t~H_h,labeller = label_both)

ggplot(scen0b,aes(x,y)) + geom_raster(aes(fill=richness)) + facet_grid(H_t~H_h,labeller = label_both) + a
ggplot(scen1b,aes(x,y)) + geom_raster(aes(fill=richness)) + facet_grid(H_t~H_h,labeller = label_both) + a
ggplot(scen2b,aes(x,y)) + geom_raster(aes(fill=richness)) + facet_grid(H_t~H_h,labeller = label_both) + a
ggplot(scen3b,aes(x,y)) + geom_raster(aes(fill=richness)) + facet_grid(H_t~H_h,labeller = label_both) + a

ggplot(scen0e,aes(x,y)) + geom_raster(aes(fill=richness)) + facet_grid(H_t~H_h,labeller = label_both) + a
ggplot(scen1e,aes(x,y)) + geom_raster(aes(fill=richness)) + facet_grid(H_t~H_h,labeller = label_both) + a
ggplot(scen2e,aes(x,y)) + geom_raster(aes(fill=richness)) + facet_grid(H_t~H_h,labeller = label_both) + a
ggplot(scen3e,aes(x,y)) + geom_raster(aes(fill=richness)) + facet_grid(H_t~H_h,labeller = label_both) + a

# Shannon index

ggplot(scen0i,aes(x,y)) + geom_raster(aes(fill=shannon)) + facet_grid(H_t~H_h,labeller = label_both) + b
ggplot(scen1i,aes(x,y)) + geom_raster(aes(fill=shannon)) + facet_grid(H_t~H_h,labeller = label_both) + b
ggplot(scen2i,aes(x,y)) + geom_raster(aes(fill=shannon)) + facet_grid(H_t~H_h,labeller = label_both) + b
ggplot(scen3i,aes(x,y)) + geom_raster(aes(fill=shannon)) + facet_grid(H_t~H_h,labeller = label_both) + b

ggplot(scen0b,aes(x,y)) + geom_raster(aes(fill=shannon)) + facet_grid(H_t~H_h,labeller = label_both) + b
ggplot(scen1b,aes(x,y)) + geom_raster(aes(fill=shannon)) + facet_grid(H_t~H_h,labeller = label_both) + b
ggplot(scen2b,aes(x,y)) + geom_raster(aes(fill=shannon)) + facet_grid(H_t~H_h,labeller = label_both) + b
ggplot(scen3b,aes(x,y)) + geom_raster(aes(fill=shannon)) + facet_grid(H_t~H_h,labeller = label_both) + b

ggplot(scen0e,aes(x,y)) + geom_raster(aes(fill=shannon)) + facet_grid(H_t~H_h,labeller = label_both) + b
ggplot(scen1e,aes(x,y)) + geom_raster(aes(fill=shannon)) + facet_grid(H_t~H_h,labeller = label_both) + b
ggplot(scen2e,aes(x,y)) + geom_raster(aes(fill=shannon)) + facet_grid(H_t~H_h,labeller = label_both) + b
ggplot(scen3e,aes(x,y)) + geom_raster(aes(fill=shannon)) + facet_grid(H_t~H_h,labeller = label_both) + b



# Simpson index

ggplot(scen0i,aes(x,y)) + geom_raster(aes(fill=simpson)) + facet_grid(H_t~H_h,labeller = label_both) + c
ggplot(scen1i,aes(x,y)) + geom_raster(aes(fill=simpson)) + facet_grid(H_t~H_h,labeller = label_both) + c
ggplot(scen2i,aes(x,y)) + geom_raster(aes(fill=simpson)) + facet_grid(H_t~H_h,labeller = label_both) + c
ggplot(scen3i,aes(x,y)) + geom_raster(aes(fill=simpson)) + facet_grid(H_t~H_h,labeller = label_both) + c

ggplot(scen0b,aes(x,y)) + geom_raster(aes(fill=simpson)) + facet_grid(H_t~H_h,labeller = label_both) + c
ggplot(scen1b,aes(x,y)) + geom_raster(aes(fill=simpson)) + facet_grid(H_t~H_h,labeller = label_both) + c
ggplot(scen2b,aes(x,y)) + geom_raster(aes(fill=simpson)) + facet_grid(H_t~H_h,labeller = label_both) + c
ggplot(scen3b,aes(x,y)) + geom_raster(aes(fill=simpson)) + facet_grid(H_t~H_h,labeller = label_both) + c

ggplot(scen0e,aes(x,y)) + geom_raster(aes(fill=simpson)) + facet_grid(H_t~H_h,labeller = label_both) + c
ggplot(scen1e,aes(x,y)) + geom_raster(aes(fill=simpson)) + facet_grid(H_t~H_h,labeller = label_both) + c
ggplot(scen2e,aes(x,y)) + geom_raster(aes(fill=simpson)) + facet_grid(H_t~H_h,labeller = label_both) + c
ggplot(scen3e,aes(x,y)) + geom_raster(aes(fill=simpson)) + facet_grid(H_t~H_h,labeller = label_both) + c

# Mean habitat fitness (niche match)

ggplot(scen0i,aes(x,y)) + geom_raster(aes(fill=fh)) + facet_grid(H_t~H_h,labeller = label_both) + d
ggplot(scen1i,aes(x,y)) + geom_raster(aes(fill=fh)) + facet_grid(H_t~H_h,labeller = label_both) + d
ggplot(scen2i,aes(x,y)) + geom_raster(aes(fill=fh)) + facet_grid(H_t~H_h,labeller = label_both) + d
ggplot(scen3i,aes(x,y)) + geom_raster(aes(fill=fh)) + facet_grid(H_t~H_h,labeller = label_both) + d

ggplot(scen0b,aes(x,y)) + geom_raster(aes(fill=fh)) + facet_grid(H_t~H_h,labeller = label_both) + d
ggplot(scen1b,aes(x,y)) + geom_raster(aes(fill=fh)) + facet_grid(H_t~H_h,labeller = label_both) + d
ggplot(scen2b,aes(x,y)) + geom_raster(aes(fill=fh)) + facet_grid(H_t~H_h,labeller = label_both) + d
ggplot(scen3b,aes(x,y)) + geom_raster(aes(fill=fh)) + facet_grid(H_t~H_h,labeller = label_both) + d

ggplot(scen0e,aes(x,y)) + geom_raster(aes(fill=fh)) + facet_grid(H_t~H_h,labeller = label_both) + d
ggplot(scen1e,aes(x,y)) + geom_raster(aes(fill=fh)) + facet_grid(H_t~H_h,labeller = label_both) + d
ggplot(scen2e,aes(x,y)) + geom_raster(aes(fill=fh)) + facet_grid(H_t~H_h,labeller = label_both) + d
ggplot(scen3e,aes(x,y)) + geom_raster(aes(fill=fh)) + facet_grid(H_t~H_h,labeller = label_both) + d

# Mean temperature fitness (niche match )

ggplot(scen0i,aes(x,y)) + geom_raster(aes(fill=ft)) + facet_grid(H_t~H_h,labeller = label_both) + d
ggplot(scen1i,aes(x,y)) + geom_raster(aes(fill=ft)) + facet_grid(H_t~H_h,labeller = label_both) + d
ggplot(scen2i,aes(x,y)) + geom_raster(aes(fill=ft)) + facet_grid(H_t~H_h,labeller = label_both) + d
ggplot(scen3i,aes(x,y)) + geom_raster(aes(fill=ft)) + facet_grid(H_t~H_h,labeller = label_both) + d

ggplot(scen0b,aes(x,y)) + geom_raster(aes(fill=ft)) + facet_grid(H_t~H_h,labeller = label_both) + d
ggplot(scen1b,aes(x,y)) + geom_raster(aes(fill=ft)) + facet_grid(H_t~H_h,labeller = label_both) + d
ggplot(scen2b,aes(x,y)) + geom_raster(aes(fill=ft)) + facet_grid(H_t~H_h,labeller = label_both) + d
ggplot(scen3b,aes(x,y)) + geom_raster(aes(fill=ft)) + facet_grid(H_t~H_h,labeller = label_both) + d

ggplot(scen0e,aes(x,y)) + geom_raster(aes(fill=ft)) + facet_grid(H_t~H_h,labeller = label_both) + d
ggplot(scen1e,aes(x,y)) + geom_raster(aes(fill=ft)) + facet_grid(H_t~H_h,labeller = label_both) + d
ggplot(scen2e,aes(x,y)) + geom_raster(aes(fill=ft)) + facet_grid(H_t~H_h,labeller = label_both) + d
ggplot(scen3e,aes(x,y)) + geom_raster(aes(fill=ft)) + facet_grid(H_t~H_h,labeller = label_both) + d

# Mean overall fitness

ggplot(scen0i,aes(x,y)) + geom_raster(aes(fill=f)) + facet_grid(H_t~H_h,labeller = label_both) + d
ggplot(scen1i,aes(x,y)) + geom_raster(aes(fill=f)) + facet_grid(H_t~H_h,labeller = label_both) + d
ggplot(scen2i,aes(x,y)) + geom_raster(aes(fill=f)) + facet_grid(H_t~H_h,labeller = label_both) + d
ggplot(scen3i,aes(x,y)) + geom_raster(aes(fill=f)) + facet_grid(H_t~H_h,labeller = label_both) + d

ggplot(scen0b,aes(x,y)) + geom_raster(aes(fill=f)) + facet_grid(H_t~H_h,labeller = label_both) + d
ggplot(scen1b,aes(x,y)) + geom_raster(aes(fill=f)) + facet_grid(H_t~H_h,labeller = label_both) + d
ggplot(scen2b,aes(x,y)) + geom_raster(aes(fill=f)) + facet_grid(H_t~H_h,labeller = label_both) + d
ggplot(scen3b,aes(x,y)) + geom_raster(aes(fill=f)) + facet_grid(H_t~H_h,labeller = label_both) + d

ggplot(scen0e,aes(x,y)) + geom_raster(aes(fill=f)) + facet_grid(H_t~H_h,labeller = label_both) + d
ggplot(scen1e,aes(x,y)) + geom_raster(aes(fill=f)) + facet_grid(H_t~H_h,labeller = label_both) + d
ggplot(scen2e,aes(x,y)) + geom_raster(aes(fill=f)) + facet_grid(H_t~H_h,labeller = label_both) + d
ggplot(scen3e,aes(x,y)) + geom_raster(aes(fill=f)) + facet_grid(H_t~H_h,labeller = label_both) + d

# Mean temperature optimum

ggplot(scen0i,aes(x,y)) + geom_raster(aes(fill=T_opt)) + facet_grid(H_t~H_h,labeller = label_both) + e
ggplot(scen1i,aes(x,y)) + geom_raster(aes(fill=T_opt)) + facet_grid(H_t~H_h,labeller = label_both) + e
ggplot(scen2i,aes(x,y)) + geom_raster(aes(fill=T_opt)) + facet_grid(H_t~H_h,labeller = label_both) + e
ggplot(scen3i,aes(x,y)) + geom_raster(aes(fill=T_opt)) + facet_grid(H_t~H_h,labeller = label_both) + e

ggplot(scen0b,aes(x,y)) + geom_raster(aes(fill=T_opt)) + facet_grid(H_t~H_h,labeller = label_both) + e 
ggplot(scen1b,aes(x,y)) + geom_raster(aes(fill=T_opt)) + facet_grid(H_t~H_h,labeller = label_both) + e
ggplot(scen2b,aes(x,y)) + geom_raster(aes(fill=T_opt)) + facet_grid(H_t~H_h,labeller = label_both) + e
ggplot(scen3b,aes(x,y)) + geom_raster(aes(fill=T_opt)) + facet_grid(H_t~H_h,labeller = label_both) + e

ggplot(scen0e,aes(x,y)) + geom_raster(aes(fill=T_opt)) + facet_grid(H_t~H_h,labeller = label_both) + e
ggplot(scen1e,aes(x,y)) + geom_raster(aes(fill=T_opt)) + facet_grid(H_t~H_h,labeller = label_both) + e
ggplot(scen2e,aes(x,y)) + geom_raster(aes(fill=T_opt)) + facet_grid(H_t~H_h,labeller = label_both) + e
ggplot(scen3e,aes(x,y)) + geom_raster(aes(fill=T_opt)) + facet_grid(H_t~H_h,labeller = label_both) + e

# Mean temperature tolerance

ggplot(scen0i,aes(x,y)) + geom_raster(aes(fill=T_sd)) + facet_grid(H_t~H_h,labeller = label_both)# + b
ggplot(scen1i,aes(x,y)) + geom_raster(aes(fill=T_sd)) + facet_grid(H_t~H_h,labeller = label_both)# + b
ggplot(scen2i,aes(x,y)) + geom_raster(aes(fill=T_sd)) + facet_grid(H_t~H_h,labeller = label_both)# + b
ggplot(scen3i,aes(x,y)) + geom_raster(aes(fill=T_sd)) + facet_grid(H_t~H_h,labeller = label_both)# + b

ggplot(scen0b,aes(x,y)) + geom_raster(aes(fill=T_sd)) + facet_grid(H_t~H_h,labeller = label_both) + b
ggplot(scen1b,aes(x,y)) + geom_raster(aes(fill=T_sd)) + facet_grid(H_t~H_h,labeller = label_both) + b
ggplot(scen2b,aes(x,y)) + geom_raster(aes(fill=T_sd)) + facet_grid(H_t~H_h,labeller = label_both) + b
ggplot(scen3b,aes(x,y)) + geom_raster(aes(fill=T_sd)) + facet_grid(H_t~H_h,labeller = label_both) + b

ggplot(scen0e,aes(x,y)) + geom_raster(aes(fill=T_sd)) + facet_grid(H_t~H_h,labeller = label_both) + b
ggplot(scen1e,aes(x,y)) + geom_raster(aes(fill=T_sd)) + facet_grid(H_t~H_h,labeller = label_both) + b
ggplot(scen2e,aes(x,y)) + geom_raster(aes(fill=T_sd)) + facet_grid(H_t~H_h,labeller = label_both) + b
ggplot(scen3e,aes(x,y)) + geom_raster(aes(fill=T_sd)) + facet_grid(H_t~H_h,labeller = label_both) + b

# Mean habitat optimum

ggplot(scen0i,aes(x,y)) + geom_raster(aes(fill=H_opt)) + facet_grid(H_t~H_h,labeller = label_both) + d
ggplot(scen1i,aes(x,y)) + geom_raster(aes(fill=H_opt)) + facet_grid(H_t~H_h,labeller = label_both) + d
ggplot(scen2i,aes(x,y)) + geom_raster(aes(fill=H_opt)) + facet_grid(H_t~H_h,labeller = label_both) + d
ggplot(scen3i,aes(x,y)) + geom_raster(aes(fill=H_opt)) + facet_grid(H_t~H_h,labeller = label_both) + d

ggplot(scen0b,aes(x,y)) + geom_raster(aes(fill=H_opt)) + facet_grid(H_t~H_h,labeller = label_both) + d
ggplot(scen1b,aes(x,y)) + geom_raster(aes(fill=H_opt)) + facet_grid(H_t~H_h,labeller = label_both) + d
ggplot(scen2b,aes(x,y)) + geom_raster(aes(fill=H_opt)) + facet_grid(H_t~H_h,labeller = label_both) + d
ggplot(scen3b,aes(x,y)) + geom_raster(aes(fill=H_opt)) + facet_grid(H_t~H_h,labeller = label_both) + d

ggplot(scen0e,aes(x,y)) + geom_raster(aes(fill=H_opt)) + facet_grid(H_t~H_h,labeller = label_both) + d
ggplot(scen1e,aes(x,y)) + geom_raster(aes(fill=H_opt)) + facet_grid(H_t~H_h,labeller = label_both) + d
ggplot(scen2e,aes(x,y)) + geom_raster(aes(fill=H_opt)) + facet_grid(H_t~H_h,labeller = label_both) + d
ggplot(scen3e,aes(x,y)) + geom_raster(aes(fill=H_opt)) + facet_grid(H_t~H_h,labeller = label_both) + d

# Mean habitat tolerance

ggplot(scen0i,aes(x,y)) + geom_raster(aes(fill=H_sd)) + facet_grid(H_t~H_h,labeller = label_both) + c
ggplot(scen1i,aes(x,y)) + geom_raster(aes(fill=H_sd)) + facet_grid(H_t~H_h,labeller = label_both) + c
ggplot(scen2i,aes(x,y)) + geom_raster(aes(fill=H_sd)) + facet_grid(H_t~H_h,labeller = label_both) + c
ggplot(scen3i,aes(x,y)) + geom_raster(aes(fill=H_sd)) + facet_grid(H_t~H_h,labeller = label_both) + c

ggplot(scen0b,aes(x,y)) + geom_raster(aes(fill=H_sd)) + facet_grid(H_t~H_h,labeller = label_both) + c
ggplot(scen1b,aes(x,y)) + geom_raster(aes(fill=H_sd)) + facet_grid(H_t~H_h,labeller = label_both) + c
ggplot(scen2b,aes(x,y)) + geom_raster(aes(fill=H_sd)) + facet_grid(H_t~H_h,labeller = label_both) + c
ggplot(scen3b,aes(x,y)) + geom_raster(aes(fill=H_sd)) + facet_grid(H_t~H_h,labeller = label_both) + c

ggplot(scen0e,aes(x,y)) + geom_raster(aes(fill=H_sd)) + facet_grid(H_t~H_h,labeller = label_both) + c
ggplot(scen1e,aes(x,y)) + geom_raster(aes(fill=H_sd)) + facet_grid(H_t~H_h,labeller = label_both) + c
ggplot(scen2e,aes(x,y)) + geom_raster(aes(fill=H_sd)) + facet_grid(H_t~H_h,labeller = label_both) + c
ggplot(scen3e,aes(x,y)) + geom_raster(aes(fill=H_sd)) + facet_grid(H_t~H_h,labeller = label_both) + c

# Mean general dispersal probability

ggplot(scen0i,aes(x,y)) + geom_raster(aes(fill=disp_l)) + facet_grid(H_t~H_h,labeller = label_both) + d
ggplot(scen1i,aes(x,y)) + geom_raster(aes(fill=disp_l)) + facet_grid(H_t~H_h,labeller = label_both) + d
ggplot(scen2i,aes(x,y)) + geom_raster(aes(fill=disp_l)) + facet_grid(H_t~H_h,labeller = label_both) + d
ggplot(scen3i,aes(x,y)) + geom_raster(aes(fill=disp_l)) + facet_grid(H_t~H_h,labeller = label_both) + d

ggplot(scen0b,aes(x,y)) + geom_raster(aes(fill=disp_l)) + facet_grid(H_t~H_h,labeller = label_both) + d
ggplot(scen1b,aes(x,y)) + geom_raster(aes(fill=disp_l)) + facet_grid(H_t~H_h,labeller = label_both) + d
ggplot(scen2b,aes(x,y)) + geom_raster(aes(fill=disp_l)) + facet_grid(H_t~H_h,labeller = label_both) + d
ggplot(scen3b,aes(x,y)) + geom_raster(aes(fill=disp_l)) + facet_grid(H_t~H_h,labeller = label_both) + d

ggplot(scen0e,aes(x,y)) + geom_raster(aes(fill=disp_l)) + facet_grid(H_t~H_h,labeller = label_both) + d
ggplot(scen1e,aes(x,y)) + geom_raster(aes(fill=disp_l)) + facet_grid(H_t~H_h,labeller = label_both) + d
ggplot(scen2e,aes(x,y)) + geom_raster(aes(fill=disp_l)) + facet_grid(H_t~H_h,labeller = label_both) + d
ggplot(scen3e,aes(x,y)) + geom_raster(aes(fill=disp_l)) + facet_grid(H_t~H_h,labeller = label_both) + d

# Mean global dispersal probability

ggplot(scen0i,aes(x,y)) + geom_raster(aes(fill=disp_g)) + facet_grid(H_t~H_h,labeller = label_both) + d
ggplot(scen1i,aes(x,y)) + geom_raster(aes(fill=disp_g)) + facet_grid(H_t~H_h,labeller = label_both) + d
ggplot(scen2i,aes(x,y)) + geom_raster(aes(fill=disp_g)) + facet_grid(H_t~H_h,labeller = label_both) + d
ggplot(scen3i,aes(x,y)) + geom_raster(aes(fill=disp_g)) + facet_grid(H_t~H_h,labeller = label_both) + d

ggplot(scen0b,aes(x,y)) + geom_raster(aes(fill=disp_g)) + facet_grid(H_t~H_h,labeller = label_both) + d
ggplot(scen1b,aes(x,y)) + geom_raster(aes(fill=disp_g)) + facet_grid(H_t~H_h,labeller = label_both) + d
ggplot(scen2b,aes(x,y)) + geom_raster(aes(fill=disp_g)) + facet_grid(H_t~H_h,labeller = label_both) + d
ggplot(scen3b,aes(x,y)) + geom_raster(aes(fill=disp_g)) + facet_grid(H_t~H_h,labeller = label_both) + d

ggplot(scen0e,aes(x,y)) + geom_raster(aes(fill=disp_g)) + facet_grid(H_t~H_h,labeller = label_both) + d
ggplot(scen1e,aes(x,y)) + geom_raster(aes(fill=disp_g)) + facet_grid(H_t~H_h,labeller = label_both) + d
ggplot(scen2e,aes(x,y)) + geom_raster(aes(fill=disp_g)) + facet_grid(H_t~H_h,labeller = label_both) + d
ggplot(scen3e,aes(x,y)) + geom_raster(aes(fill=disp_g)) + facet_grid(H_t~H_h,labeller = label_both) + d

#------------------------------------------------------
# Scatter plots

ggplot(scen0i,aes(T_opt,T_sd)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)
ggplot(scen1i,aes(T_opt,T_sd)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)
ggplot(scen2i,aes(T_opt,T_sd)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)
ggplot(scen3i,aes(T_opt,T_sd)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)

ggplot(scen0b,aes(T_opt,T_sd)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)
ggplot(scen1b,aes(T_opt,T_sd)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)
ggplot(scen2b,aes(T_opt,T_sd)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)
ggplot(scen3b,aes(T_opt,T_sd)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)

ggplot(scen0e,aes(T_opt,T_sd)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)
ggplot(scen1e,aes(T_opt,T_sd)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)
ggplot(scen2e,aes(T_opt,T_sd)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)
ggplot(scen3e,aes(T_opt,T_sd)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)

#-------------------------------------------------------

ggplot(scen0i,aes(T_opt,H_opt)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)
ggplot(scen1i,aes(T_opt,H_opt)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)
ggplot(scen2i,aes(T_opt,H_opt)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)
ggplot(scen3i,aes(T_opt,H_opt)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)

ggplot(scen0b,aes(T_opt,H_opt)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)
ggplot(scen1b,aes(T_opt,H_opt)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)
ggplot(scen2b,aes(T_opt,H_opt)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)
ggplot(scen3b,aes(T_opt,H_opt)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)

ggplot(scen0e,aes(T_opt,H_opt)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)
ggplot(scen1e,aes(T_opt,H_opt)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)
ggplot(scen2e,aes(T_opt,H_opt)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)
ggplot(scen3e,aes(T_opt,H_opt)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)

#-------------------------------------------------------

ggplot(scen0i,aes(T_opt,disp_l)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)
ggplot(scen1i,aes(T_opt,disp_l)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)
ggplot(scen2i,aes(T_opt,disp_l)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)
ggplot(scen3i,aes(T_opt,disp_l)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)

ggplot(scen0b,aes(T_opt,disp_l)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)
ggplot(scen1b,aes(T_opt,disp_l)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)
ggplot(scen2b,aes(T_opt,disp_l)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)
ggplot(scen3b,aes(T_opt,disp_l)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)

ggplot(scen0e,aes(T_opt,disp_l)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)
ggplot(scen1e,aes(T_opt,disp_l)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)
ggplot(scen2e,aes(T_opt,disp_l)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)
ggplot(scen3e,aes(T_opt,disp_l)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)

#------------------------------------------------------

ggplot(scen0i,aes(T_sd,disp_l)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)
ggplot(scen1i,aes(T_sd,disp_l)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)
ggplot(scen2i,aes(T_sd,disp_l)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)
ggplot(scen3i,aes(T_sd,disp_l)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)

ggplot(scen0b,aes(T_sd,disp_l)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)
ggplot(scen1b,aes(T_sd,disp_l)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)
ggplot(scen2b,aes(T_sd,disp_l)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)
ggplot(scen3b,aes(T_sd,disp_l)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)

ggplot(scen0e,aes(T_sd,disp_l)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)
ggplot(scen1e,aes(T_sd,disp_l)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)
ggplot(scen2e,aes(T_sd,disp_l)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)
ggplot(scen3e,aes(T_sd,disp_l)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)

#-------------------------------------------------------

ggplot(scen0i,aes(H_opt,disp_l)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)
ggplot(scen1i,aes(H_opt,disp_l)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)
ggplot(scen2i,aes(H_opt,disp_l)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)
ggplot(scen3i,aes(H_opt,disp_l)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)

ggplot(scen0b,aes(H_opt,disp_l)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)
ggplot(scen1b,aes(H_opt,disp_l)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)
ggplot(scen2b,aes(H_opt,disp_l)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)
ggplot(scen3b,aes(H_opt,disp_l)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)

ggplot(scen0e,aes(H_opt,disp_l)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)
ggplot(scen1e,aes(H_opt,disp_l)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)
ggplot(scen2e,aes(H_opt,disp_l)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)
ggplot(scen3e,aes(H_opt,disp_l)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)

#-------------------------------------------------------

ggplot(scen0i,aes(H_sd,disp_l)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)
ggplot(scen1i,aes(H_sd,disp_l)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)
ggplot(scen2i,aes(H_sd,disp_l)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)
ggplot(scen3i,aes(H_sd,disp_l)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)

ggplot(scen0b,aes(H_sd,disp_l)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)
ggplot(scen1b,aes(H_sd,disp_l)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)
ggplot(scen2b,aes(H_sd,disp_l)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)
ggplot(scen3b,aes(H_sd,disp_l)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)

ggplot(scen0e,aes(H_sd,disp_l)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)
ggplot(scen1e,aes(H_sd,disp_l)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)
ggplot(scen2e,aes(H_sd,disp_l)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)
ggplot(scen3e,aes(H_sd,disp_l)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)

#--------------------------------------------------------

ggplot(scen0i,aes(T_opt,disp_g)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)
ggplot(scen1i,aes(T_opt,disp_g)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)
ggplot(scen2i,aes(T_opt,disp_g)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)
ggplot(scen3i,aes(T_opt,disp_g)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)

ggplot(scen0b,aes(T_opt,disp_g)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)
ggplot(scen1b,aes(T_opt,disp_g)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)
ggplot(scen2b,aes(T_opt,disp_g)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)
ggplot(scen3b,aes(T_opt,disp_g)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)

ggplot(scen0e,aes(T_opt,disp_g)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)
ggplot(scen1e,aes(T_opt,disp_g)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)
ggplot(scen2e,aes(T_opt,disp_g)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)
ggplot(scen3e,aes(T_opt,disp_g)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)

#--------------------------------------------------------

ggplot(scen0i,aes(T_sd,disp_g)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)
ggplot(scen1i,aes(T_sd,disp_g)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)
ggplot(scen2i,aes(T_sd,disp_g)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)
ggplot(scen3i,aes(T_sd,disp_g)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)

ggplot(scen0b,aes(T_sd,disp_g)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)
ggplot(scen1b,aes(T_sd,disp_g)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)
ggplot(scen2b,aes(T_sd,disp_g)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)
ggplot(scen3b,aes(T_sd,disp_g)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)

ggplot(scen0e,aes(T_sd,disp_g)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)
ggplot(scen1e,aes(T_sd,disp_g)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)
ggplot(scen2e,aes(T_sd,disp_g)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)
ggplot(scen3e,aes(T_sd,disp_g)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)

#--------------------------------------------------------

ggplot(scen0i,aes(H_opt,disp_g)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)
ggplot(scen1i,aes(H_opt,disp_g)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)
ggplot(scen2i,aes(H_opt,disp_g)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)
ggplot(scen3i,aes(H_opt,disp_g)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)

ggplot(scen0b,aes(H_opt,disp_g)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)
ggplot(scen1b,aes(H_opt,disp_g)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)
ggplot(scen2b,aes(H_opt,disp_g)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)
ggplot(scen3b,aes(H_opt,disp_g)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)

ggplot(scen0e,aes(H_opt,disp_g)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)
ggplot(scen1e,aes(H_opt,disp_g)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)
ggplot(scen2e,aes(H_opt,disp_g)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)
ggplot(scen3e,aes(H_opt,disp_g)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)

#--------------------------------------------------------

ggplot(scen0i,aes(H_sd,disp_g)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)
ggplot(scen1i,aes(H_sd,disp_g)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)
ggplot(scen2i,aes(H_sd,disp_g)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)
ggplot(scen3i,aes(H_sd,disp_g)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)

ggplot(scen0b,aes(H_sd,disp_g)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)
ggplot(scen1b,aes(H_sd,disp_g)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)
ggplot(scen2b,aes(H_sd,disp_g)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)
ggplot(scen3b,aes(H_sd,disp_g)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)

ggplot(scen0e,aes(H_sd,disp_g)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)
ggplot(scen1e,aes(H_sd,disp_g)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)
ggplot(scen2e,aes(H_sd,disp_g)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)
ggplot(scen3e,aes(H_sd,disp_g)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)

#-------------------------------------------------------

ggplot(scen0i,aes(temp_t,disp_l)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)
ggplot(scen1i,aes(temp_t,disp_l)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)
ggplot(scen2i,aes(temp_t,disp_l)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)
ggplot(scen3i,aes(temp_t,disp_l)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)

ggplot(scen0b,aes(temp_t,disp_l)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)
ggplot(scen1b,aes(temp_t,disp_l)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)
ggplot(scen2b,aes(temp_t,disp_l)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)
ggplot(scen3b,aes(temp_t,disp_l)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)

ggplot(scen0e,aes(temp_t,disp_l)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)
ggplot(scen1e,aes(temp_t,disp_l)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)
ggplot(scen2e,aes(temp_t,disp_l)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)
ggplot(scen3e,aes(temp_t,disp_l)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)

#------------------------------------------------------

ggplot(scen0,aes(temp_t,disp_g)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)
ggplot(scen1,aes(temp_t,disp_g)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)
ggplot(scen2,aes(temp_t,disp_g)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)
ggplot(scen3,aes(temp_t,disp_g)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)

#-----------------------------------------------------

ggplot(scen0i,aes(habitat,disp_l)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)
ggplot(scen1i,aes(habitat,disp_l)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)
ggplot(scen2i,aes(habitat,disp_l)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)
ggplot(scen3i,aes(habitat,disp_l)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)

ggplot(scen0b,aes(habitat,disp_l)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)
ggplot(scen1b,aes(habitat,disp_l)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)
ggplot(scen2b,aes(habitat,disp_l)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)
ggplot(scen3b,aes(habitat,disp_l)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)

ggplot(scen0e,aes(habitat,disp_l)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)
ggplot(scen1e,aes(habitat,disp_l)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)
ggplot(scen2e,aes(habitat,disp_l)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)
ggplot(scen3e,aes(habitat,disp_l)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)

#-----------------------------------------------------

ggplot(scen0i,aes(habitat,disp_g)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)
ggplot(scen1i,aes(habitat,disp_g)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)
ggplot(scen2i,aes(habitat,disp_g)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)
ggplot(scen3i,aes(habitat,disp_g)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)

ggplot(scen0b,aes(habitat,disp_g)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)
ggplot(scen1b,aes(habitat,disp_g)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)
ggplot(scen2b,aes(habitat,disp_g)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)
ggplot(scen3b,aes(habitat,disp_g)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)

ggplot(scen0e,aes(habitat,disp_g)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)
ggplot(scen1e,aes(habitat,disp_g)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)
ggplot(scen2e,aes(habitat,disp_g)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)
ggplot(scen3e,aes(habitat,disp_g)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)

#-----------------------------------------------------
  
ggplot(scen0i,aes(T_sd,H_sd)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both) + xlim(0,3) + ylim(0,3)
ggplot(scen1i,aes(T_sd,H_sd)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both) + xlim(0,3) + ylim(0,3)
ggplot(scen2i,aes(T_sd,H_sd)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both) + xlim(0,3) + ylim(0,3)
ggplot(scen3i,aes(T_sd,H_sd)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both) + xlim(0,3) + ylim(0,3)

ggplot(scen0b,aes(T_sd,H_sd)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both) + xlim(0,3) + ylim(0,3)
ggplot(scen1b,aes(T_sd,H_sd)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both) + xlim(0,3) + ylim(0,3)
ggplot(scen2b,aes(T_sd,H_sd)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both) + xlim(0,3) + ylim(0,3)
ggplot(scen3b,aes(T_sd,H_sd)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both) + xlim(0,3) + ylim(0,3)

ggplot(scen0e,aes(T_sd,H_sd)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both) + xlim(0,3) + ylim(0,3)
ggplot(scen1e,aes(T_sd,H_sd)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both) + xlim(0,3) + ylim(0,3)
ggplot(scen2e,aes(T_sd,H_sd)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both) + xlim(0,3) + ylim(0,3)
ggplot(scen3e,aes(T_sd,H_sd)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both) + xlim(0,3) + ylim(0,3)

ggplot(scen0i,aes(disp_l,disp_g)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)
ggplot(scen1i,aes(disp_l,disp_g)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)
ggplot(scen2i,aes(disp_l,disp_g)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)
ggplot(scen3i,aes(disp_l,disp_g)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)

ggplot(scen0b,aes(disp_l,disp_g)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)
ggplot(scen1b,aes(disp_l,disp_g)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)
ggplot(scen2b,aes(disp_l,disp_g)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)
ggplot(scen3b,aes(disp_l,disp_g)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)

ggplot(scen0e,aes(disp_l,disp_g)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)
ggplot(scen1e,aes(disp_l,disp_g)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)
ggplot(scen2e,aes(disp_l,disp_g)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)
ggplot(scen3e,aes(disp_l,disp_g)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)


ggplot(scen0i,aes(f,pop)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)
ggplot(scen1i,aes(f,pop)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)
ggplot(scen2i,aes(f,pop)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)
ggplot(scen3i,aes(f,pop)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)

ggplot(scen0b,aes(f,pop)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)
ggplot(scen1b,aes(f,pop)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)
ggplot(scen2b,aes(f,pop)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)
ggplot(scen3b,aes(f,pop)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)

ggplot(scen0e,aes(f,pop)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)
ggplot(scen1e,aes(f,pop)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)
ggplot(scen2e,aes(f,pop)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)
ggplot(scen3e,aes(f,pop)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)
  