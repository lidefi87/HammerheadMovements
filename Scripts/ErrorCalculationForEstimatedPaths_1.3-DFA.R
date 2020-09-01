############################################################################################
# Title: Error calculation and plotting
# Author: Denisse Fierro Arcos and Jeremy Vaudo
# Version: 1.3
# Date created: 2020-05-02
# Date last updated: 2020-08-04 DFA
# 
############################################################################################

#Clear global environment
rm(list = ls())

################################ CALLING RELEVANT LIBRARIES ################################
{ library(HMMoce)
  library(tidyverse)
  library(raster)
  library(geosphere)
  library(rgdal)
  library(sf)} #sf added to manipulate data points easier - DFA

############################ MANUALLY SET VARIABLES FOR ANALYSIS ###########################
#Removing getting current directory as project script is saved within project - DFA

#Importing shapefile with GMR boundaries for data visualisation and to categorise points
#as being inside or outside the reserve - DFA
Gal <- read_sf("Data/Layers/ReservaMarinaWGS84.shp")

#Identify folder where shapefile and csv outputs will be saved - Back folder added because
#this script is linked to an R project with a set path - DFA
OutFolder <- "../Spatial/MiniPATs_NewProcess/"
if(!dir.exists(OutFolder)){
  dir.create(OutFolder, recursive = T)}

#Uploading file with information about actual track end
EndTrack <- openxlsx::read.xlsx("Data/MiniPAT_Tags/HH_MiniPAT_ActualTrackEnd_2016-19.xlsx",
                                detectDates = T) %>% janitor::clean_names()

# Set confidence interval to be plotted
CI <- 90

########################### GET CONFIDENCE INTERVALS FOR TRACKS #####################################
#Creating variable to save total distance covered by each tagged shark
x <- data.frame()
#Creating variable to save mean change in long when sharks were in and out of the GMR 
DiffLong <- data.frame()
DiffLongAll <- data.frame()

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
  #Extracting actual track end for tag - DFA
  end <- EndTrack %>% filter(tag_id == str_extract(i, "[0-9]{6}"))
  #Removing data after actual track end - DFA
  df <- df %>% filter(date <= end$track_end_date)
  
  #Extracting actual track from res data prior to converting to shapefile - DFA
  shp <- res$tr %>% filter(date <= end$track_end_date) %>% 
    #Creating new column for months and ordering it based on month name - DFA
    mutate(month = factor(months(date), levels = month.name, ordered = T),
           #Second column containing the lunar phase for location estimation date - DFA
           moonPhase = lunar::lunar.phase(date, name = 8)) %>% 
    #Reordering months based on the year of tracking - DFA
    mutate(month = forcats::fct_reorder(month, lubridate::year(date))) %>% 
    #Calculate absolute change in longitude between points - DFA
    mutate(DifLong = c(NA, abs(diff(lon))))
  
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
                 colour = 'grey') + theme_bw() +
    #Adding colour based on moon phase at time of location estimation - DFA
    geom_point(data = shp, aes(x = lon, y = lat, colour = DifLong))
    
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
  #Changing factor column into character
  shp <- shp %>% mutate(month = as.character(month)) %>% 
    st_as_sf(coords = c("lon", "lat"), crs = 4326)
  
  #Visualising differences in longitude change
  ggplot()+geom_sf(data = Gal)+geom_sf(data = shp, aes(color = DifLong))
  
  #Creating new GMR column which specifies if estimate was in or out of the GMR to 
  #calculate mean values of change - DFA
  xx <- shp %>% 
    # mutate(GMR = case_when(st_intersects(shp, Gal, sparse = F) == T ~ "In", 
    #                              st_intersects(shp, Gal, sparse = F) == F ~ "Out")) %>% 
    #Extracting coordinates from geometry to include them dataset
    mutate(lon = st_coordinates(.)[,1], lat = st_coordinates(.)[,2]) %>% 
    #Classifying points based on their location relative to W91.75
    mutate(GMR = case_when(lon <= -91.75 ~ "In",
                           lon > -91.75 ~ "Out")) %>% 
    #Essentially turning shapefile into data frame
    st_drop_geometry() %>% 
    #Add an extra column with infromation about the ptt (tag)
    mutate(ptt = str_extract(i, "[0-9]{6}"))
  DiffLongAll <- rbind(DiffLongAll, xx)
  
  yy <- xx %>% 
    #Grouping variables based on location of point in relation to GMR
    group_by(GMR) %>% 
    summarise(meanDifLong = mean(DifLong, na.rm = T))
  
  #Adding this information to empty data frame so a summary of all tags can be saved at
  #the end - DFA
  DiffLong <- rbind(DiffLong, yy)
  #Removing unused variable
  rm(yy)
  
  #Saving shapefile with name of best estimated track for the tag
  if(grepl("157566", i) | grepl("178974", i)){
    write_sf(shp, paste0(OutFolder, TagResults$name[2], ".shp"), overwrite = T)
  }else{write_sf(shp, paste0(OutFolder, TagResults$name[1], ".shp"), overwrite = T)}
  
  #Creating and saving shapefile of confidence intervals for best fit track - JV
  ConInt <- SpatialPolygons(mapply(function(x, id){
    xy <- dplyr::select(x, -date) %>% Polygon()
    Polygons(list(Polygon(xy)), id)}, 
    group_split(df, date), unique(df$date)), 
    proj4string = CRS("+proj=longlat +datum=WGS84"))

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
  distPts <- res$tr %>% filter(date <= end$track_end_date) %>% 
    as_tibble() %>%
    #distGeo() provides distance in m, but it is divided by 1000 to 
    #present distance in Km
    mutate(dist_km =  head(c(NA, round(distGeo(.[, c("lon", "lat")])/1000, 3)), -1),
    #Calculating speed as m per second
    vel_ms = dist_km*1000/86400) %>% 
    #Calculating absolute change in longitude between points
    mutate(DifLong = c(NA, abs(diff(lon))))
  
  #Adding up the distance covered by each tagged shark - DFA
  x <- distPts %>% summarise(TotDist = sum(dist_km, na.rm = T)) %>%
    mutate(Tag = str_extract(i, "[0-9]{6}")) %>%
    bind_rows(x, .)
  
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

#Save the file containing distances covered per shark and changes in longitude - DFA
write.csv(x, "../Outputs/Distance_covered_per_shark.csv")
write.csv(DiffLong, "../Outputs/DailyLongitudeDifference.csv")



# Statistical analysis ----------------------------------------------------
#Uploading relevant libraries
library(ggpubr)
library(rstatix)

#Tags with return trips
TagsRet <- c(157566:157567, 174049:174051, 178974)

#Extracting these tags from the combined dataste
TagsInt <- DiffLongAll %>% 
  filter(ptt %in% TagsRet) %>% 
  #Remove rows with NA values
  drop_na(DifLong)

TagsInt %>% 
  group_by(GMR, ptt) %>% 
  #Calculating summary stats (mean and SD)
  get_summary_stats(DifLong, type = "mean_sd") %>% 
  #Keep tags with 3 or more observations
  filter(n >= 3) %>%
  #Keeping only tags that appear in and out of the GMR for comparison
  group_by(ptt) %>%
  #Save the tag numbers to filter dataset
  filter(n() > 1) %>% 
  distinct(ptt) -> Tags

#Keeping tags with enough points for comparison
TagsComp <- TagsInt %>% filter(ptt %in% Tags$ptt)

#Visualising results
TagsComp %>% 
  ggboxplot(x = "GMR", y = "DifLong", color = "ptt", add = "point")
TagsComp %>% 
  ggboxplot(x = "GMR", y = "DifLong", add = "point")

#Identifying outliers
ExtOutliers <- TagsComp %>% 
  #Selecting relevant columns
  select(ptt, DifLong, GMR) %>% 
  #Grouping by location and tag
  group_by(GMR) %>% 
  #Identifying outliers per location and tag
  identify_outliers(., DifLong) %>% 
  #Only keeping extreme outliers
  filter(is.extreme == TRUE)

#Checking normality
#QQ plots
TagsComp %>% 
  ggqqplot("DifLong", ggtheme = theme_bw())+facet_grid(~GMR)
#Shapiro Wilk's test
TagsComp %>% 
  group_by(GMR) %>% 
  shapiro_test(DifLong)

#PERMANOVA
library(vegan)
#Using data where a tag has at least three observations
TagsComp %>% 
  #Grouping by location and tag
  group_by(GMR, ptt) %>% 
  #Testing for differences in location, while controlling for tag (PTT)
  adonis(DifLong ~ GMR, ., permutations = 9999, strata = .$ptt)

TagsComp %>% group_by(GMR) %>% 
  summarise(MeanLonDif = mean(DifLong, na.rm = T))

#End