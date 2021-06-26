############################################################################################
# Title: Error calculation and plotting
# Author: Denisse Fierro Arcos and Jeremy Vaudo
# Version: 1.2
# Date created: 2020-05-02
# Date last updated: 2020-05-12 DFA
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
#Removing getting current directory as project script is saved within project - DFA

#Identify folder where shapefile and csv outputs will be saved - Back folder added because
#this script is linked to an R project with a set path - DFA
OutFolder <- "../Spatial/MiniPATs_NewProcess/"
if(!dir.exists(OutFolder)){
  dir.create(OutFolder, recursive = T)}

# Set confidence interval to be plotted
CI <- 90

########################### GET CONFIDENCE INTERVALS FOR TRACKS #####################################
#Creating variable to save total distance covered by each tagged shark
# x <- data.frame()

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
  #Counter has been removed, instead dates have been extracted and included in final dataset
  #because we will use them to check when the animal left the GMR - DFA. Section initially
  #added by JV.
  df <- do.call(rbind, lapply(bnds, FUN = function(x){
    data.frame(x[["ctr"]]@lines[[1]]@Lines[[1]]@coords,
               #Including a new date column - DFA
               date =  x[["loc"]]$date)}))
  
  ###Plotting directly in R
  ## make a plot
  world <- map_data('world')
  xl <- c(min(res$tr$lon) - 2, max(res$tr$lon) + 2)
  yl <- c(min(res$tr$lat) - 2, max(res$tr$lat) + 2)
  
  ## simple map of movement data
  worldmap <- ggplot() + geom_polygon(data = world, aes(x = long, y = lat, group = group)) + 
    coord_fixed(xlim = xl, ylim = yl, ratio = 1.3) + xlab('') + ylab('')
  
  trackplot <- worldmap +
    geom_polygon(data = df, aes(x = x, y = y, group = date), alpha = 0.1, fill = 'grey', 
                 colour = 'grey') + 
    geom_point(data = res$tr, aes(x = lon, y = lat)) +
    theme_bw()
  #Visualising newly produced map - DFA
  trackplot
  
  #Create a ggplotMaps folder if one does not yet exist - DFA
  if(!dir.exists(paste0(OutFolder, "ggplotMaps/"))){
    dir.create(paste0(OutFolder, "ggplotMaps/"), recursive = T)}
  #Saving resulting plot
  if(grepl("157566", i) | grepl("178974", i)){
    ggsave(paste0(OutFolder, "ggplotMaps/", TagResults$name[2], ".tiff"), trackplot, 
           "tiff", dpi = 400)
  }else{ggsave(paste0(OutFolder, "ggplotMaps/", TagResults$name[1], ".tiff"), trackplot, 
               "tiff", dpi = 400)}
  
  ### Creating and exporting shapefile with CRS: WGS84
  shp <- res$tr
  coordinates(shp) = ~lon+lat
  proj4string(shp) = CRS("+proj=longlat +datum=WGS84")
  
  #Saving shapefile with name of best estimated track for the tag
  if(grepl("157566", i) | grepl("178974", i)){
    shapefile(shp, paste0(OutFolder, TagResults$name[2], ".shp"), overwrite = T)
  }else{shapefile(shp, paste0(OutFolder, TagResults$name[1], ".shp"), overwrite = T)}
  
  #Creating and saving shapefile of confidence intervals for best fit track - JV
  ConInt <- SpatialPolygons(mapply(function(x, id){
    xy <- dplyr::select(x, -date) %>% Polygon()
    Polygons(list(Polygon(xy)), id)}, 
    group_split(df, date), unique(df$date)), 
    proj4string = CRS("++proj=longlat +datum=WGS84"))

  ConInt.df <- SpatialPolygonsDataFrame(ConInt, data.frame(id = unique(df$date), 
                                                           row.names = unique(df$date)))
  #Saving shapefile - JV
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
  
  #Adding up the distance covered by each tagged shark
  # x <- distPts %>% summarise(TotDist = sum(dist_km, na.rm = T)) %>% 
  #   mutate(Tag = str_extract(i, "[0-9]{6}")) %>% 
  #   bind_rows(x, .)
  
  #Create a csv folder if one does not yet exist - JV
  if(!dir.exists(paste0(OutFolder,"csvFiles/"))){
    dir.create(paste0(OutFolder,"csvFiles/"), recursive = T)} 
  #Saving distance and speed calculations with the name of best estimated track for the tag
  if(grepl("157566", i) | grepl("178974", i)){
    write.csv(distPts, paste0(OutFolder, "csvFiles/", TagResults$name[2], ".csv"), 
              row.names = F)
  }else{write.csv(distPts, paste0(OutFolder, "csvFiles/", TagResults$name[1], 
                                   ".csv"), row.names = F)}
  }