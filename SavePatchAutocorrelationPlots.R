
#=============================================
# Main script
#=============================================

# Clear workspace
rm(list=ls())

#-------------------------------------------

# Packages
library(readr)
library(dplyr)
library(vegan)
library(data.table)
library(ggplot2)


p1 <- "Documents/PhD/LandklifModel/Landscapes/Model1/Experiment1/G1Patches"
files <- list.files(path = p1, pattern = "*.txt",full.names = TRUE,recursive = FALSE)


lapply(files, function(x){
  print(x)
  y <- sub(".*_","",sub(".txt","",x))
  print(y)
  pdir <- paste(p1,"/",y,"/",sep="")
  print(pdir)
  if(dir.exists(pdir)==TRUE){
    savedir <- pdir
    print(savedir)
  } else {
    dir.create(pdir)
    savedir <- pdir
    print(savedir)
  }
  
  
  patches <- read_csv(x, col_types = cols(H_h = col_factor(levels = c("0", "0.5", "1", "Uniform")),
                      H_t = col_factor(levels = c("0","0.5", "1", "Uniform"))))
  patches$H_t <- as.factor(patches$H_t)
  patches$H_h <- as.factor(patches$H_h)
  
  levels(patches$H_t)
  levels(patches$H_h)
  
  # Subsetting
  scen0 <- subset(patches, patches$clim_scen==0)
  scen1 <- subset(patches, patches$clim_scen==1)
  scen2 <- subset(patches, patches$clim_scen==2)
  scen3 <- subset(patches, patches$clim_scen==3)
  
  
  # Population size
  
  
  
  ggplot(scen0,aes(pop)) + geom_histogram(binwidth = 1) + facet_grid(H_t~H_h,labeller = label_both) + xlim(0,370) + ylim(0,100)
  ggsave(paste(savedir,"test6_pop_hist_0.png",sep=""))
  ggplot(scen1,aes(pop)) + geom_histogram(binwidth = 1) + facet_grid(H_t~H_h,labeller = label_both) + xlim(0,370) + ylim(0,100)
  ggsave(paste(savedir,"test6_pop_hist_1.png",sep=""))
  ggplot(scen2,aes(pop)) + geom_histogram(binwidth = 1) + facet_grid(H_t~H_h,labeller = label_both) + xlim(0,370) + ylim(0,100)
  ggsave(paste(savedir,"test6_pop_hist_2.png",sep=""))
  ggplot(scen3,aes(pop)) + geom_histogram(binwidth = 1) + facet_grid(H_t~H_h,labeller = label_both) + xlim(0,370) + ylim(0,100)
  ggsave(paste(savedir,"test6_pop_hist_4.png",sep=""))
  
  # No. of lineages
  ggplot(scen0,aes(lineages)) + geom_histogram() + facet_grid(H_t~H_h,labeller = label_both)
  ggsave(paste(savedir,"test6_lineages_hist_0.png",sep=""))
  ggplot(scen1,aes(lineages)) + geom_histogram() + facet_grid(H_t~H_h,labeller = label_both)
  ggsave(paste(savedir,"test6_lineages_hist_1.png",sep=""))
  ggplot(scen2,aes(lineages)) + geom_histogram() + facet_grid(H_t~H_h,labeller = label_both)
  ggsave(paste(savedir,"test6_lineages_hist_2.png",sep=""))
  ggplot(scen3,aes(lineages)) + geom_histogram() + facet_grid(H_t~H_h,labeller = label_both)
  ggsave(paste(savedir,"test6_lineages_hist_3.png",sep=""))
  
  # Shannon index
  ggplot(scen0,aes(shannon)) + geom_density() + facet_grid(H_t~H_h,labeller = label_both)
  ggplot(scen1,aes(shannon)) + geom_density() + facet_grid(H_t~H_h,labeller = label_both)
  ggplot(scen2,aes(shannon)) + geom_density() + facet_grid(H_t~H_h,labeller = label_both)
  ggplot(scen3,aes(shannon)) + geom_density() + facet_grid(H_t~H_h,labeller = label_both)
  
  # Simpson index
  ggplot(scen0,aes(simpson)) + geom_density() + facet_grid(H_t~H_h,labeller = label_both)
  ggplot(scen1,aes(simpson)) + geom_density() + facet_grid(H_t~H_h,labeller = label_both)
  ggplot(scen2,aes(simpson)) + geom_density() + facet_grid(H_t~H_h,labeller = label_both)
  ggplot(scen3,aes(simpson)) + geom_density() + facet_grid(H_t~H_h,labeller = label_both)
  
  # Habitat fitness (niche match)
  ggplot(scen0,aes(fit.h)) + geom_histogram(binwidth = 0.01) + facet_grid(H_t~H_h,labeller = label_both) + xlim(0,1)
  ggsave(paste(savedir,"test6_fith_hist_0.png",sep=""))
  ggplot(scen1,aes(fit.h)) + geom_histogram(binwidth = 0.01) + facet_grid(H_t~H_h,labeller = label_both) + xlim(0,1)
  ggsave(paste(savedir,"test6_fith_hist_1.png",sep=""))
  ggplot(scen2,aes(fit.h)) + geom_histogram(binwidth = 0.01) + facet_grid(H_t~H_h,labeller = label_both) + xlim(0,1)
  ggsave(paste(savedir,"test6_fith_hist_2.png",sep=""))
  ggplot(scen3,aes(fit.h)) + geom_histogram(binwidth = 0.01) + facet_grid(H_t~H_h,labeller = label_both) + xlim(0,1)
  ggsave(paste(savedir,"test6_fith_hist_3.png",sep=""))
  
  # Temperature fitness (niche match)
  ggplot(scen0,aes(fit.t)) + geom_histogram(binwidth = 0.01) + facet_grid(H_t~H_h,labeller = label_both) + xlim(0,1)
  ggsave(paste(savedir,"test6_fitt_hist_0.png",sep=""))
  ggplot(scen1,aes(fit.t)) + geom_histogram(binwidth = 0.01) + facet_grid(H_t~H_h,labeller = label_both) + xlim(0,1)
  ggsave(paste(savedir,"test6_fitt_hist_1.png",sep=""))
  ggplot(scen2,aes(fit.t)) + geom_histogram(binwidth = 0.01) + facet_grid(H_t~H_h,labeller = label_both) + xlim(0,1)
  ggsave(paste(savedir,"test6_fitt_hist_2.png",sep=""))
  ggplot(scen3,aes(fit.t)) + geom_histogram(binwidth = 0.01) + facet_grid(H_t~H_h,labeller = label_both) + xlim(0,1)
  ggsave(paste(savedir,"test6_fitt_hist_3.png",sep=""))
  
  # Overall fitness 
  ggplot(scen0,aes(fitness)) + geom_histogram(binwidth = 0.01) + facet_grid(H_t~H_h,labeller = label_both) + xlim(0,1)
  ggsave(paste(savedir,"test6_fitness_hist_0.png",sep=""))
  ggplot(scen1,aes(fitness)) + geom_histogram(binwidth = 0.01) + facet_grid(H_t~H_h,labeller = label_both) + xlim(0,1)
  ggsave(paste(savedir,"test6_fitness_hist_1.png",sep=""))
  ggplot(scen2,aes(fitness)) + geom_histogram(binwidth = 0.01) + facet_grid(H_t~H_h,labeller = label_both) + xlim(0,1)
  ggsave(paste(savedir,"test6_fitness_hist_2.png",sep=""))
  ggplot(scen3,aes(fitness)) + geom_histogram(binwidth = 0.01) + facet_grid(H_t~H_h,labeller = label_both) + xlim(0,1)
  ggsave(paste(savedir,"test6_fitness_hist_3.png",sep=""))
  
  # Mean temperature optimum
  ggplot(scen0,aes(T_opt)) + geom_histogram(binwidth = 0.05) + facet_grid(H_t~H_h,labeller = label_both) + xlim(10,15)
  ggsave(paste(savedir,"test6_topt_hist_0.png",sep=""))
  ggplot(scen1,aes(T_opt)) + geom_histogram(binwidth = 0.05) + facet_grid(H_t~H_h,labeller = label_both) + xlim(10,15)
  ggsave(paste(savedir,"test6_topt_hist_1.png",sep=""))
  ggplot(scen2,aes(T_opt)) + geom_histogram(binwidth = 0.05) + facet_grid(H_t~H_h,labeller = label_both) + xlim(10,15)
  ggsave(paste(savedir,"test6_topt_hist_2.png",sep=""))
  ggplot(scen3,aes(T_opt)) + geom_histogram(binwidth = 0.05) + facet_grid(H_t~H_h,labeller = label_both) + xlim(10,15)
  ggsave(paste(savedir,"test6_topt_hist_3.png",sep=""))
  
  # Mean temperature tolerance
  ggplot(scen0,aes(T_sd)) + geom_histogram(binwidth = 0.05) + facet_grid(H_t~H_h,labeller = label_both) + xlim(0,3)
  ggsave(paste(savedir,"test6_tsd_hist_0.png",sep=""))
  ggplot(scen1,aes(T_sd)) + geom_histogram(binwidth = 0.05) + facet_grid(H_t~H_h,labeller = label_both) + xlim(0,3)
  ggsave(paste(savedir,"test6_tsd_hist_1.png",sep=""))
  ggplot(scen2,aes(T_sd)) + geom_histogram(binwidth = 0.05) + facet_grid(H_t~H_h,labeller = label_both) + xlim(0,3)
  ggsave(paste(savedir,"test6_tsd_hist_2.png",sep=""))
  ggplot(scen3,aes(T_sd)) + geom_histogram(binwidth = 0.05) + facet_grid(H_t~H_h,labeller = label_both) + xlim(0,3)
  ggsave(paste(savedir,"test6_tsd_hist_3.png",sep=""))
  
  # Mean habitat optimum
  ggplot(scen0,aes(H_opt)) + geom_histogram(binwidth = 0.01) + facet_grid(H_t~H_h,labeller = label_both) + xlim(0,1)
  ggsave(paste(savedir,"test6_hopt_hist_0.png",sep=""))
  ggplot(scen1,aes(H_opt)) + geom_histogram(binwidth = 0.01) + facet_grid(H_t~H_h,labeller = label_both) + xlim(0,1)
  ggsave(paste(savedir,"test6_hopt_hist_1.png",sep=""))
  ggplot(scen2,aes(H_opt)) + geom_histogram(binwidth = 0.01) + facet_grid(H_t~H_h,labeller = label_both) + xlim(0,1)
  ggsave(paste(savedir,"test6_hopt_hist_2.png",sep=""))
  ggplot(scen3,aes(H_opt)) + geom_histogram(binwidth = 0.01) + facet_grid(H_t~H_h,labeller = label_both) + xlim(0,1)
  ggsave(paste(savedir,"test6_hopt_hist_3.png",sep=""))
  
  # Mean habitat tolerance
  ggplot(scen0,aes(H_sd)) + geom_histogram(binwidth = 0.05) + facet_grid(H_t~H_h,labeller = label_both) + xlim(0,3)
  ggsave(paste(savedir,"test6_hsd_hist_0.png",sep=""))
  ggplot(scen1,aes(H_sd)) + geom_histogram(binwidth = 0.05) + facet_grid(H_t~H_h,labeller = label_both) + xlim(0,3)
  ggsave(paste(savedir,"test6_hsd_hist_1.png",sep=""))
  ggplot(scen2,aes(H_sd)) + geom_histogram(binwidth = 0.05) + facet_grid(H_t~H_h,labeller = label_both) + xlim(0,3)
  ggsave(paste(savedir,"test6_hsd_hist_2.png",sep=""))
  ggplot(scen3,aes(H_sd)) + geom_histogram(binwidth = 0.05) + facet_grid(H_t~H_h,labeller = label_both) + xlim(0,3)
  ggsave(paste(savedir,"test6_hsd_hist_3.png",sep=""))
  
  # Mean dispersal probability
  ggplot(scen0,aes(disp_l)) + geom_histogram(binwidth = 0.01) + facet_grid(H_t~H_h,labeller = label_both) + xlim(0,1)
  ggsave(paste(savedir,"test6_displ_hist_0.png",sep=""))
  ggplot(scen1,aes(disp_l)) + geom_histogram(binwidth = 0.01) + facet_grid(H_t~H_h,labeller = label_both) + xlim(0,1)
  ggsave(paste(savedir,"test6_displ_hist_1.png",sep=""))
  ggplot(scen2,aes(disp_l)) + geom_histogram(binwidth = 0.01) + facet_grid(H_t~H_h,labeller = label_both) + xlim(0,1)
  ggsave(paste(savedir,"test6_displ_hist_2.png",sep=""))
  ggplot(scen3,aes(disp_l)) + geom_histogram(binwidth = 0.01) + facet_grid(H_t~H_h,labeller = label_both) + xlim(0,1)
  ggsave(paste(savedir,"test6_displ_hist_3.png",sep=""))
  
  # Mean global dispersal probability
  ggplot(scen0,aes(disp_g)) + geom_histogram(binwidth = 0.01) + facet_grid(H_t~H_h,labeller = label_both) + xlim(0,1)
  ggsave(paste(savedir,"test6_dispg_hist_0.png",sep=""))
  ggplot(scen1,aes(disp_g)) + geom_histogram(binwidth = 0.01) + facet_grid(H_t~H_h,labeller = label_both) + xlim(0,1)
  ggsave(paste(savedir,"test6_dispg_hist_1.png",sep=""))
  ggplot(scen2,aes(disp_g)) + geom_histogram(binwidth = 0.01) + facet_grid(H_t~H_h,labeller = label_both) + xlim(0,1)
  ggsave(paste(savedir,"test6_dispg_hist_2.png",sep=""))
  ggplot(scen3,aes(disp_g)) + geom_histogram(binwidth = 0.01) + facet_grid(H_t~H_h,labeller = label_both) + xlim(0,1)
  ggsave(paste(savedir,"test6_dispg_hist_3.png",sep=""))

  
  #=====================================================
  # Landscape raster plots by trait and climate scenario
  #=====================================================
  
  a <- scale_fill_gradientn(colours = hcl.colors(100),limits = c(1,6))
  b <- scale_fill_gradientn(colours = hcl.colors(100),limits = c(0,1.6))
  c <- scale_fill_gradientn(colours = hcl.colors(100),limits = c(1,370))
  d <- scale_fill_gradientn(colours = hcl.colors(100),limits = c(0,1))
  e <- scale_fill_gradientn(colours = hcl.colors(100),limits = c(11,14))
  
  
  # Population size
  
  ggplot(scen0,aes(x,y)) + geom_raster(aes(fill=pop)) + facet_grid(H_t~H_h,labeller = label_both) + c
  ggsave(paste(savedir,"test6_pop_map_0.png",sep=""))
  ggplot(scen1,aes(x,y)) + geom_raster(aes(fill=pop)) + facet_grid(H_t~H_h,labeller = label_both) + c
  ggsave(paste(savedir,"test6_pop_map_1.png",sep=""))
  ggplot(scen2,aes(x,y)) + geom_raster(aes(fill=pop)) + facet_grid(H_t~H_h,labeller = label_both) + c
  ggsave(paste(savedir,"test6_pop_map_2.png",sep=""))
  ggplot(scen3,aes(x,y)) + geom_raster(aes(fill=pop)) + facet_grid(H_t~H_h,labeller = label_both) + c
  ggsave(paste(savedir,"test6_pop_map_3.png",sep=""))
  
  # No. of Lineages
  
  ggplot(scen0,aes(x,y)) + geom_raster(aes(fill=lineages)) + facet_grid(H_t~H_h,labeller = label_both) + a
  ggsave(paste(savedir,"test6_lineages_map_0.png",sep=""))
  ggplot(scen1,aes(x,y)) + geom_raster(aes(fill=lineages)) + facet_grid(H_t~H_h,labeller = label_both) + a
  ggsave(paste(savedir,"test6_lineages_map_1.png",sep=""))
  ggplot(scen2,aes(x,y)) + geom_raster(aes(fill=lineages)) + facet_grid(H_t~H_h,labeller = label_both) + a
  ggsave(paste(savedir,"test6_lineages_map_2.png",sep=""))
  ggplot(scen3,aes(x,y)) + geom_raster(aes(fill=lineages)) + facet_grid(H_t~H_h,labeller = label_both) + a
  ggsave(paste(savedir,"test6_lineages_map_3.png",sep=""))
  
  # Mean habitat fitness (niche match)
  
  ggplot(scen0,aes(x,y)) + geom_raster(aes(fill=fit.h)) + facet_grid(H_t~H_h,labeller = label_both) + d
  ggsave(paste(savedir,"test6_fith_map_0.png",sep=""))
  ggplot(scen1,aes(x,y)) + geom_raster(aes(fill=fit.h)) + facet_grid(H_t~H_h,labeller = label_both) + d
  ggsave(paste(savedir,"test6_fith_map_1.png",sep=""))
  ggplot(scen2,aes(x,y)) + geom_raster(aes(fill=fit.h)) + facet_grid(H_t~H_h,labeller = label_both) + d
  ggsave(paste(savedir,"test6_fith_map_2.png",sep=""))
  ggplot(scen3,aes(x,y)) + geom_raster(aes(fill=fit.h)) + facet_grid(H_t~H_h,labeller = label_both) + d
  ggsave(paste(savedir,"test6_fith_map_3.png",sep=""))
  
  # Mean temperature fitness (niche match )
  
  ggplot(scen0,aes(x,y)) + geom_raster(aes(fill=fit.t)) + facet_grid(H_t~H_h,labeller = label_both) + d
  ggsave(paste(savedir,"test6_fitt_map_0.png",sep=""))
  ggplot(scen1,aes(x,y)) + geom_raster(aes(fill=fit.t)) + facet_grid(H_t~H_h,labeller = label_both) + d
  ggsave(paste(savedir,"test6_fitt_map_1.png",sep=""))
  ggplot(scen2,aes(x,y)) + geom_raster(aes(fill=fit.t)) + facet_grid(H_t~H_h,labeller = label_both) + d
  ggsave(paste(savedir,"test6_fitt_map_2.png",sep=""))
  ggplot(scen3,aes(x,y)) + geom_raster(aes(fill=fit.t)) + facet_grid(H_t~H_h,labeller = label_both) + d
  ggsave(paste(savedir,"test6_fitt_map_3.png",sep=""))
  
  # Mean overall fitness
  
  ggplot(scen0,aes(x,y)) + geom_raster(aes(fill=fitness)) + facet_grid(H_t~H_h,labeller = label_both) + d
  ggsave(paste(savedir,"test6_fitness_map_0.png",sep=""))
  ggplot(scen1,aes(x,y)) + geom_raster(aes(fill=fitness)) + facet_grid(H_t~H_h,labeller = label_both) + d
  ggsave(paste(savedir,"test6_fitness_map_1.png",sep=""))
  ggplot(scen2,aes(x,y)) + geom_raster(aes(fill=fitness)) + facet_grid(H_t~H_h,labeller = label_both) + d
  ggsave(paste(savedir,"test6_fitness_map_2.png",sep=""))
  ggplot(scen3,aes(x,y)) + geom_raster(aes(fill=fitness)) + facet_grid(H_t~H_h,labeller = label_both) + d
  ggsave(paste(savedir,"test6_fitness_map_3.png",sep=""))
  
  # Mean temperature optimum
  
  ggplot(scen0,aes(x,y)) + geom_raster(aes(fill=T_opt)) + facet_grid(H_t~H_h,labeller = label_both) + e
  ggsave(paste(savedir,"test6_topt_map_0.png",sep=""))
  ggplot(scen1,aes(x,y)) + geom_raster(aes(fill=T_opt)) + facet_grid(H_t~H_h,labeller = label_both) + e
  ggsave(paste(savedir,"test6_topt_map_1.png",sep=""))
  ggplot(scen2,aes(x,y)) + geom_raster(aes(fill=T_opt)) + facet_grid(H_t~H_h,labeller = label_both) + e
  ggsave(paste(savedir,"test6_topt_map_2.png",sep=""))
  ggplot(scen3,aes(x,y)) + geom_raster(aes(fill=T_opt)) + facet_grid(H_t~H_h,labeller = label_both) + e
  ggsave(paste(savedir,"test6_topt_map_3.png",sep=""))
  
  # Mean temperature tolerance
  
  ggplot(scen0,aes(x,y)) + geom_raster(aes(fill=T_sd)) + facet_grid(H_t~H_h,labeller = label_both) + b
  ggsave(paste(savedir,"test6_tsd_map_0.png",sep=""))
  ggplot(scen1,aes(x,y)) + geom_raster(aes(fill=T_sd)) + facet_grid(H_t~H_h,labeller = label_both) + b
  ggsave(paste(savedir,"test6_tsd_map_1.png",sep=""))
  ggplot(scen2,aes(x,y)) + geom_raster(aes(fill=T_sd)) + facet_grid(H_t~H_h,labeller = label_both) + b
  ggsave(paste(savedir,"test6_tsd_map_2.png",sep=""))
  ggplot(scen3,aes(x,y)) + geom_raster(aes(fill=T_sd)) + facet_grid(H_t~H_h,labeller = label_both) + b
  ggsave(paste(savedir,"test6_tsd_map_3.png",sep=""))
  
  # Mean habitat optimum
  
  ggplot(scen0,aes(x,y)) + geom_raster(aes(fill=H_opt)) + facet_grid(H_t~H_h,labeller = label_both) + d
  ggsave(paste(savedir,"test6_hopt_map_0.png",sep=""))
  ggplot(scen1,aes(x,y)) + geom_raster(aes(fill=H_opt)) + facet_grid(H_t~H_h,labeller = label_both) + d
  ggsave(paste(savedir,"test6_hopt_map_1.png",sep=""))
  ggplot(scen2,aes(x,y)) + geom_raster(aes(fill=H_opt)) + facet_grid(H_t~H_h,labeller = label_both) + d
  ggsave(paste(savedir,"test6_hopt_map_2.png",sep=""))
  ggplot(scen3,aes(x,y)) + geom_raster(aes(fill=H_opt)) + facet_grid(H_t~H_h,labeller = label_both) + d
  ggsave(paste(savedir,"test6_hopt_map_3.png",sep=""))
  
  # Mean habitat tolerance
  
  ggplot(scen0,aes(x,y)) + geom_raster(aes(fill=H_sd)) + facet_grid(H_t~H_h,labeller = label_both) + c
  ggsave(paste(savedir,"test6_hsd_map_0.png",sep=""))
  ggplot(scen1,aes(x,y)) + geom_raster(aes(fill=H_sd)) + facet_grid(H_t~H_h,labeller = label_both) + c
  ggsave(paste(savedir,"test6_hsd_map_1.png",sep=""))
  ggplot(scen2,aes(x,y)) + geom_raster(aes(fill=H_sd)) + facet_grid(H_t~H_h,labeller = label_both) + c
  ggsave(paste(savedir,"test6_hsd_map_2.png",sep=""))
  ggplot(scen3,aes(x,y)) + geom_raster(aes(fill=H_sd)) + facet_grid(H_t~H_h,labeller = label_both) + c
  ggsave(paste(savedir,"test6_hsd_map_3.png",sep=""))
  
  # Mean general dispersal probability
  
  ggplot(scen0,aes(x,y)) + geom_raster(aes(fill=disp_l)) + facet_grid(H_t~H_h,labeller = label_both) + d
  ggsave(paste(savedir,"test6_displ_map_0.png",sep=""))
  ggplot(scen1,aes(x,y)) + geom_raster(aes(fill=disp_l)) + facet_grid(H_t~H_h,labeller = label_both) + d
  ggsave(paste(savedir,"test6_displ_map_1.png",sep=""))
  ggplot(scen2,aes(x,y)) + geom_raster(aes(fill=disp_l)) + facet_grid(H_t~H_h,labeller = label_both) + d
  ggsave(paste(savedir,"test6_displ_map_2.png",sep=""))
  ggplot(scen3,aes(x,y)) + geom_raster(aes(fill=disp_l)) + facet_grid(H_t~H_h,labeller = label_both) + d
  ggsave(paste(savedir,"test6_displ_map_3.png",sep=""))
  
  # Mean global dispersal probability
  
  ggplot(scen0,aes(x,y)) + geom_raster(aes(fill=disp_g)) + facet_grid(H_t~H_h,labeller = label_both) + d
  ggsave(paste(savedir,"test6_dispg_map_0.png",sep=""))
  ggplot(scen1,aes(x,y)) + geom_raster(aes(fill=disp_g)) + facet_grid(H_t~H_h,labeller = label_both) + d
  ggsave(paste(savedir,"test6_dispg_map_1.png",sep=""))
  ggplot(scen2,aes(x,y)) + geom_raster(aes(fill=disp_g)) + facet_grid(H_t~H_h,labeller = label_both) + d
  ggsave(paste(savedir,"test6_dispg_map_2.png",sep=""))
  ggplot(scen3,aes(x,y)) + geom_raster(aes(fill=disp_g)) + facet_grid(H_t~H_h,labeller = label_both) + d
  ggsave(paste(savedir,"test6_dispg_map_3.png",sep=""))
  
  #========================================================
  # Scatterplots
  #========================================================
  
  ggplot(scen0,aes(T_opt,T_sd)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)
  ggsave(paste(savedir,"test6_topt_tsd_scat_0.png",sep=""))
  ggplot(scen1,aes(T_opt,T_sd)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)
  ggsave(paste(savedir,"test6_topt_tsd_scat_1.png",sep=""))
  ggplot(scen2,aes(T_opt,T_sd)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)
  ggsave(paste(savedir,"test6_topt_tsd_scat_2.png",sep=""))
  ggplot(scen3,aes(T_opt,T_sd)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)
  ggsave(paste(savedir,"test6_topt_tsd_scat_3.png",sep=""))
  
  ggplot(scen0,aes(T_opt,H_opt)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)
  ggsave(paste(savedir,"test6_topt_hopt_scat_0.png",sep=""))
  ggplot(scen1,aes(T_opt,H_opt)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)
  ggsave(paste(savedir,"test6_topt_hopt_scat_1.png",sep=""))
  ggplot(scen2,aes(T_opt,H_opt)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)
  ggsave(paste(savedir,"test6_topt_hopt_scat_2.png",sep=""))
  ggplot(scen3,aes(T_opt,H_opt)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)
  ggsave(paste(savedir,"test6_topt_hopt_scat_3.png",sep=""))
  
  ggplot(scen0,aes(T_opt,disp_l)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)
  ggsave(paste(savedir,"test6_topt_displ_scat_0.png",sep=""))
  ggplot(scen1,aes(T_opt,disp_l)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)
  ggsave(paste(savedir,"test6_topt_displ_scat_1.png",sep=""))
  ggplot(scen2,aes(T_opt,disp_l)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)
  ggsave(paste(savedir,"test6_topt_displ_scat_2.png",sep=""))
  ggplot(scen3,aes(T_opt,disp_l)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)
  ggsave(paste(savedir,"test6_topt_displ_scat_3.png",sep=""))
  
  ggplot(scen0,aes(T_sd,disp_l)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)
  ggsave(paste(savedir,"test6_tsd_displ_scat_0.png",sep=""))
  ggplot(scen1,aes(T_sd,disp_l)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)
  ggsave(paste(savedir,"test6_tsd_displ_scat_1.png",sep=""))
  ggplot(scen2,aes(T_sd,disp_l)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)
  ggsave(paste(savedir,"test6_tsd_displ_scat_2.png",sep=""))
  ggplot(scen3,aes(T_sd,disp_l)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)
  ggsave(paste(savedir,"test6_tsd_displ_scat_3.png",sep=""))
  
  ggplot(scen0,aes(H_opt,disp_l)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)
  ggsave(paste(savedir,"test6_hopt_displ_scat_0.png",sep=""))
  ggplot(scen1,aes(H_opt,disp_l)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)
  ggsave(paste(savedir,"test6_hopt_displ_scat_1.png",sep=""))
  ggplot(scen2,aes(H_opt,disp_l)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)
  ggsave(paste(savedir,"test6_hopt_displ_scat_2.png",sep=""))
  ggplot(scen3,aes(H_opt,disp_l)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)
  ggsave(paste(savedir,"test6_hopt_displ_scat_3.png",sep=""))
  
  ggplot(scen0,aes(H_sd,disp_l)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)
  ggsave(paste(savedir,"test6_hsd_displ_scat_0.png",sep=""))
  ggplot(scen1,aes(H_sd,disp_l)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)
  ggsave(paste(savedir,"test6_hsd_displ_scat_1.png",sep=""))
  ggplot(scen2,aes(H_sd,disp_l)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)
  ggsave(paste(savedir,"test6_hsd_displ_scat_2.png",sep=""))
  ggplot(scen3,aes(H_sd,disp_l)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)
  ggsave(paste(savedir,"test6_hsd_displ_scat_3.png",sep=""))
  
  
  ggplot(scen0,aes(T_opt,disp_g)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)
  ggsave(paste(savedir,"test6_topt_dispg_scat_0.png",sep=""))
  ggplot(scen1,aes(T_opt,disp_g)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)
  ggsave(paste(savedir,"test6_topt_dispg_scat_1.png",sep=""))
  ggplot(scen2,aes(T_opt,disp_g)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)
  ggsave(paste(savedir,"test6_topt_dispg_scat_2.png",sep=""))
  ggplot(scen3,aes(T_opt,disp_g)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)
  ggsave(paste(savedir,"test6_topt_dispg_scat_3.png",sep=""))
  
  ggplot(scen0,aes(T_sd,disp_g)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)
  ggsave(paste(savedir,"test6_tsd_dispg_scat_0.png",sep=""))
  ggplot(scen1,aes(T_sd,disp_g)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)
  ggsave(paste(savedir,"test6_tsd_dispg_scat_1.png",sep=""))
  ggplot(scen2,aes(T_sd,disp_g)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)
  ggsave(paste(savedir,"test6_tsd_dispg_scat_2.png",sep=""))
  ggplot(scen3,aes(T_sd,disp_g)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)
  ggsave(paste(savedir,"test6_tsd_dispg_scat_3.png",sep=""))
  
  ggplot(scen0,aes(H_opt,disp_g)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)
  ggsave(paste(savedir,"test6_hopt_dispg_scat_0.png",sep=""))
  ggplot(scen1,aes(H_opt,disp_g)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)
  ggsave(paste(savedir,"test6_hopt_dispg_scat_1.png",sep=""))
  ggplot(scen2,aes(H_opt,disp_g)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)
  ggsave(paste(savedir,"test6_hopt_dispg_scat_2.png",sep=""))
  ggplot(scen3,aes(H_opt,disp_g)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)
  ggsave(paste(savedir,"test6_hopt_dispg_scat_3.png",sep=""))
  
  ggplot(scen0,aes(H_sd,disp_g)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)
  ggsave(paste(savedir,"test6_hsd_dispg_scat_0.png",sep=""))
  ggplot(scen1,aes(H_sd,disp_g)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)
  ggsave(paste(savedir,"test6_hsd_dispg_scat_1.png",sep=""))
  ggplot(scen2,aes(H_sd,disp_g)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)
  ggsave(paste(savedir,"test6_hsd_dispg_scat_2.png",sep=""))
  ggplot(scen3,aes(H_sd,disp_g)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)
  ggsave(paste(savedir,"test6_hsd_dispg_scat_3.png",sep=""))
  
  
  ggplot(scen0,aes(temp_t,disp_l)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)
  ggsave(paste(savedir,"test6_tempt_displ_scat_0.png",sep=""))
  ggplot(scen1,aes(temp_t,disp_l)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)
  ggsave(paste(savedir,"test6_tempt_displ_scat_1.png",sep=""))
  ggplot(scen2,aes(temp_t,disp_l)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)
  ggsave(paste(savedir,"test6_tempt_displ_scat_2.png",sep=""))
  ggplot(scen3,aes(temp_t,disp_l)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)
  ggsave(paste(savedir,"test6_tempt_displ_scat_3.png",sep=""))
  
  ggplot(scen0,aes(temp_t,disp_g)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)
  ggsave(paste(savedir,"test6_tempt_dispg_scat_0.png",sep=""))
  ggplot(scen1,aes(temp_t,disp_g)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)
  ggsave(paste(savedir,"test6_tempt_dispg_scat_1.png",sep=""))
  ggplot(scen2,aes(temp_t,disp_g)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)
  ggsave(paste(savedir,"test6_tempt_dispg_scat_2.png",sep=""))
  ggplot(scen3,aes(temp_t,disp_g)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)
  ggsave(paste(savedir,"test6_tempt_dispg_scat_3.png",sep=""))
  
  ggplot(scen0,aes(habitat,disp_l)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)
  ggsave(paste(savedir,"test6_habitat_displ_scat_0.png",sep=""))
  ggplot(scen1,aes(habitat,disp_l)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)
  ggsave(paste(savedir,"test6_habitat_displ_scat_1.png",sep=""))
  ggplot(scen2,aes(habitat,disp_l)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)
  ggsave(paste(savedir,"test6_habitat_displ_scat_2.png",sep=""))
  ggplot(scen3,aes(habitat,disp_l)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)
  ggsave(paste(savedir,"test6_habitat_displ_scat_3.png",sep=""))
  
  ggplot(scen0,aes(habitat,disp_g)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)
  ggsave(paste(savedir,"test6_habitat_dispg_scat_0.png",sep=""))
  ggplot(scen1,aes(habitat,disp_g)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)
  ggsave(paste(savedir,"test6_habitat_dispg_scat_0.png",sep=""))
  ggplot(scen2,aes(habitat,disp_g)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)
  ggsave(paste(savedir,"test6_habitat_dispg_scat_0.png",sep=""))
  ggplot(scen3,aes(habitat,disp_g)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both)
  ggsave(paste(savedir,"test6_habitat_dispg_scat_0.png",sep=""))
  
  
  ggplot(scen0,aes(T_sd,H_sd)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both) + xlim(0,3) + ylim(0,3)
  ggsave(paste(savedir,"test6_tsd_hsd_scat_0.png",sep=""))
  ggplot(scen1,aes(T_sd,H_sd)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both) + xlim(0,3) + ylim(0,3)
  ggsave(paste(savedir,"test6_tsd_hsd_scat_0.png",sep=""))
  ggplot(scen2,aes(T_sd,H_sd)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both) + xlim(0,3) + ylim(0,3)
  ggsave(paste(savedir,"test6_tsd_hsd_scat_0.png",sep=""))
  ggplot(scen3,aes(T_sd,H_sd)) + geom_point() + facet_grid(H_t~H_h,labeller = label_both) + xlim(0,3) + ylim(0,3)
  ggsave(paste(savedir,"test6_tsd_hsd_scat_0.png",sep=""))
  # Make and save some graphs
  
  })
  