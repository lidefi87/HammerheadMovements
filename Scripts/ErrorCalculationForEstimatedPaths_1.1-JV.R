############################################################################################
# Title: Error calculation and plotting
# Author: Denisse Fierro Arcos and Jeremy Vaudo
# Version: 1.1
# Date created: 2020-05-02
# Date last updated: 2020-05-08 JV
# 
############################################################################################

#Clear global environment
rm(list = ls())

################################ CALLING RELEVANT LIBRARIES ################################
{ library(HMMoce)
  library(tidyverse)
  library(raster)
  library(geosphere)
  library(rgdal)}

############################ MANUALLY SET VARIABLES FOR ANALYSIS ###########################
#Get project (current) directory
MainDir <- getwd()

#Identify folder where shapefile and csv outputs will be saved
OutFolder <- "../Spatial/MiniPATs_NewProcess/"
if(!dir.exists(OutFolder)){
  dir.create(OutFolder, recursive = T)}

# Set confidence interval to be plotted
CI <- 90

########################### GET CONFIDENCE INTERVALS FOR TRACKS #####################################
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
  
  # Getting error for each estimated location
  ## see ?getCtr for description of each input but these are smoother, calc.track, and 
  #setup.grid outputs, respectively
  bnds <- getCtr(res$s, res$tr, res$g, threshold = (100 - CI), makePlot = F)
  
  ## gets polygons for error distributions
  counter <- 0
  df <- do.call(rbind, lapply(bnds, FUN = function(x){
    counter <<- counter + 1
    data.frame(x[["ctr"]]@lines[[1]]@Lines[[1]]@coords) %>%
      mutate(day = counter)}))
  
  ###Plotting directly in R
  ## make a plot
  world <- map_data('world')
  xl <- c(min(res$tr$lon) - 2, max(res$tr$lon) + 2)
  yl <- c(min(res$tr$lat) - 2, max(res$tr$lat) + 2)
  
  ## simple map of movement data
  worldmap <- ggplot() + geom_polygon(data = world, aes(x = long, y = lat, group = group)) + 
    coord_fixed(xlim = xl, ylim = yl, ratio = 1.3) + xlab('') + ylab('')
  
  trackplot <- worldmap +
    geom_polygon(data = df, aes(x = x, y = y, group = day), alpha = 0.1, fill = 'grey', colour = 'grey') + 
    geom_point(data = res$tr, aes(x = lon, y = lat)) +
    theme_bw()
  
  #Saving resulting plot
  if(grepl("157566", i) | grepl("178974", i)){
    ggsave(paste0(OutFolder, TagResults$name[2], ".tiff"), trackplot, "tiff", dpi = 400)
  }else{ggsave(paste0(OutFolder, TagResults$name[1], ".tiff"), trackplot, "tiff", dpi = 400)}
  
  #Saving data as csv file - Names depend on best estimated track for the tag
  #Create a csv folder if one does not yet exist
  if(!dir.exists(paste0(OutFolder,"csvFiles/"))){
    dir.create(paste0(OutFolder,"csvFiles/"), recursive = T)}
  #Save files
  if(grepl("157566", i) | grepl("178974", i)){
    write.csv(res$tr, paste0(OutFolder, "csvFiles/", TagResults$name[2], ".csv"), 
              row.names = F)
  }else{write.csv(res$tr, paste0(OutFolder, "csvFiles/", TagResults$name[1], ".csv"), 
                  row.names = F)}
  
  ### Creating and exporting shapefile with CRS: WGS84
  shp <- res$tr
  coordinates(shp) = ~lon+lat
  proj4string(shp) = CRS("++proj=longlat +datum=WGS84")
  
  #Saving shapefile with name of best estimated track for the tag
  if(grepl("157566", i) | grepl("178974", i)){
    shapefile(shp, paste0(OutFolder, TagResults$name[2], ".shp"), overwrite = T)
  }else{shapefile(shp, paste0(OutFolder, TagResults$name[1], ".shp"), overwrite = T)}
  
  #Creating and saving shapefile of confidence intervals for best fit track
  ConInt <- SpatialPolygons(mapply(function(x, id){
  xy <- dplyr::select(x, -day) %>% Polygon() 
  Polygons(list(Polygon(xy)), id)}, 
  group_split(df, day), unique(df$day)), 
  proj4string = CRS("++proj=longlat +datum=WGS84"))

  ConInt.df <- SpatialPolygonsDataFrame(ConInt, data.frame(id = unique(df$day), 
                                                           row.names = unique(df$day)))
  
  if(grepl("157566", i) | grepl("178974", i)){
    shapefile(ConInt.df, paste0(OutFolder, TagResults$name[2], "_CI.shp"), 
              overwrite = T)
  }else{shapefile(ConInt.df, paste0(OutFolder, TagResults$name[1], "_CI.shp"), overwrite = T)}
  
  ### Calculating distance between points in km
  #Creating dataset which includes dates, coordinates and distances for each estimated
  #location for that tag 
  distPts <- tibble(date = shp@data$date, 
                    lon = shp@coords[,"lon"], 
                    lat = shp@coords[,"lat"],
                    p_resident = shp@data$p,
                    #distGeo() provides distance in m, but it is divided by 1000 to 
                    #present distance in Km
                    dist_km =  head(c(NA, round(distGeo(shp)/1000, 3)), -1),
                    #Calculating speed as m per second
                    vel_ms = dist_km*1000/86400)
  
  #Saving distance and speed calculations with the name of best estimated track for the tag
  if(grepl("157566", i) | grepl("178974", i)){
    write.csv(distPts, paste0(OutFolder, TagResults$name[2], ".csv"), 
              row.names = F)
  }else{write.csv(distPts, paste0(OutFolder, TagResults$name[1], 
                                   ".csv"), row.names = F)}
  }

