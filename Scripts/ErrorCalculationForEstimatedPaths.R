############################################################################################
# Title: Error calculation and plotting
# Author: Denisse Fierro Arcos
# Version: 1
# Date last updated: 2020-05-01
# 
############################################################################################

#Clear global environment
rm(list = ls())

#Upload useful libraries
library(HMMoce)
library(tidyverse)
library(ggforce)
library(raster)
library(geosphere)

#Get project (current) directory
MainDir <- getwd()
#Identify folder where shapefile and csv outputs will be saved
OutFolder <- "../Spatial/MiniPATs_NewProcess/"
#Getting a list of file directories of datasets containing HMMoce results
HMM_results <- list.files(pattern = "HMMoce_results.csv", recursive = T)

#Create a loop using list of directories to save best results as identified by the model
for(i in HMM_results){
  #Load HMMoce result for each tag
  TagResults <- read.csv(i) %>% arrange(nll)
  #Identify best model (first row as dataset has been arranged by nll) and load results
  #into environment
  #Manually identifying the second best nll as the best track for two tags in particular
  if(grepl("157566", i) | grepl("178974", i)){
    load(list.files(pattern = glob2rx(paste0(TagResults$name[2], "-HMMoce.rda$")), 
                  recursive = T))
  }else{load(list.files(pattern = glob2rx(paste0(TagResults$name[1], "-HMMoce.rda$")), 
                        recursive = T))}
  
  # Plotting errors for calculated paths
  ## for a 90% confidence interval use threshold = 10
  ## see ?getCtr for description of each input but these are smoother, calc.track, and 
  #setup.grid outputs, respectively
  bnds <- getCtr(res$s, res$tr, res$g, threshold = 10, makePlot = F)
  
  ##Converting CIs into a lon and lat error metrics
  df <- cbind(res$tr, t(data.frame(lapply(bnds, 
                                          FUN = function(x) c(x$yDist, x$xDist))))) %>% 
    plyr::rename(c("1" = "ydist", "2" = "xdist")) %>% 
    ##some contours will likely not calculate so you can fill this in with the mean error 
    #(or something else) if you wish
    mutate(xdist = case_when(is.na(xdist) ~ mean(xdist, na.rm = T), 
                             TRUE ~ xdist),
           ydist = case_when(is.na(ydist) ~ mean(ydist, na.rm = T), 
                             TRUE ~ ydist))
  
  #Saving data as csv file - Names depend on best estimated track for the tag
  if(grepl("157566", i) | grepl("178974", i)){
    write.csv(df, paste0(OutFolder, "csvFiles/", TagResults$name[2], ".csv"), 
              row.names = F)
  }else{write.csv(df, paste0(OutFolder, "csvFiles/", TagResults$name[1], ".csv"), 
                    row.names = F)}
  
  ### Creating and exporting shapefile with CRS: WGS84
  shp <- df
  coordinates(shp) = ~lon+lat
  proj4string(shp) = CRS("++proj=longlat +datum=WGS84")
  #Saving shapefile with name of best estimated track for the tag
  if(grepl("157566", i) | grepl("178974", i)){
    shapefile(shp, paste0(OutFolder, TagResults$name[2], ".shp"), overwrite = T)
  }else{shapefile(shp, paste0(OutFolder, TagResults$name[1], ".shp"), overwrite = T)}
  
  ### Calculating distance between points in km
  #Creating dataset which includes dates, coordinates and distances for each estimated
  #location for that tag 
  distPts <- tibble(date = shp@data$date, 
                    lon = shp@coords[,"lon"], 
                    lat = shp@coords[,"lat"],
                    #distHaversine() provides distance in m, but it is divided by 1000 to 
                    #present distance in Km
                    dist_km =  c(NA, round(distHaversine(shp)/1000, 3)),
                    #Calculating speed as m per second
                    vel_ms = dist_km*1000/86400)
  
  #Saving distance and speed calculations with the name of best estimated track for the tag
  if(grepl("157566", i) | grepl("178974", i)){
    write.csv(distPts, paste0("../Outputs/BehaviourCalcs/", TagResults$name[2], ".csv"), 
              row.names = F)
  }else{write.csv(distPts, paste0( "../Outputs/BehaviourCalcs/", TagResults$name[1], 
                                   ".csv"), row.names = F)}
 
  ###Plotting directly in R
  ## make a plot
  world <- map_data('world')
  xl <- c(min(df$lon) - 2, max(df$lon) + 2)
  yl <- c(min(df$lat) - 2, max(df$lat) + 2)
  
  ## simple map of move data
  m1 <- ggplot() + geom_polygon(data = world, aes(x = long, y = lat, group = group)) + 
    coord_fixed(xlim = xl, ylim = yl, ratio = 1.3) + xlab('') + ylab('')
  
  ## add confidence intervals (use geom_ellipse from ggforce package)
  m1 <- m1 + geom_ellipse(data = df, aes(x0 = lon, y0 = lat, a = xdist, b = ydist, 
                                         angle = 0), alpha = 0.1, fill = 'grey', 
                          colour = 'grey')+
    geom_point(data = df, aes(x = lon, y = lat, color = date))+theme_bw()
  m1
  }

