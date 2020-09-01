# Title: Analysing SPOT data from Hammerhead sharks (S. lewini) in the TEP
# Author: Denisse Fierro Arcos
# Version: 4
# Date last updated: 2020-08-10

#Clear global environment
rm(list = ls())

################################### Uploading relevant libraries ##########################################
# library(rgdal)
library(raster)
library(sf)
library(rgeos)
library(tidyverse)
# library(bsam)
library(openxlsx)

############################### Uploading relevant layers for analysis ####################################
#Uploading shapefiles needed for analysis - Galapagos Islands and GMR
# GalIslWGS84 <-  shapefile("Data/Layers/galapagos_island_wgs84.shp")
GMRWGS84 <- read_sf("Data/Layers/ReservaMarinaWGS84.shp")

################################### SATELLITE DATA PREPARATION ##########################################
#Uploading supporting info for tag deployments
SuppInfo <- read.xlsx("Data/SPOT_Tags/SPOT6TagList_Deployments.xlsx", sheet = "ReleaseDates", 
                      detectDates = T)

#Getting a list of the directories within the folder - Keep files with tag numbers only
TagNumb <- dir(path = "Data/SPOT_Tags/", pattern = "[0-9]{6}$", full.names = T)
TagData <- dir(path = "Data/SPOT_Tags/", pattern = "-RawArgos.csv$", full.names = T, 
               recursive = T)

{#Variable containing location classes related to data quality (Classes: G > 3 > 2 > 1 > 0 > A > B > Z)
#Change this variable if you want to include more or less location classes
classInt = c(0:3)

#Create empty dataset to store data for all tags
combinedData <- data.frame()

#Cleaning data within each folder
for(i in seq_along(TagData)){
  #Upload file containing raw Argos data
  rawData <- read.csv(TagData[i], stringsAsFactors = F)
  #Keep only location classes specified above
  tagData <- rawData %>% as_tibble() %>% filter(Class %in% classInt) %>% 
    #Keep only variables of interest: Tag ID, date, ARGOS location class (lc), longitude, latitude
    select(PTT,  Class, PassDate, PassTime, Longitude, Latitude) %>% 
    #Eliminate duplicate entries for the same date and time
    distinct(Class, PassDate, PassTime, .keep_all = T) %>%
    #Change variable type of date column to date
    mutate(PassDate = lubridate::parse_date_time(PassDate, orders = "dmy"),
           PTT = factor(PTT)) %>% 
    #Filter out detections after the release date specified in the SuppInfo dataset
    filter(PassDate <= SuppInfo$ReleaseDate[i]) %>% 
    #Unite date and time columns and convert them to date format
    unite(Date, PassDate, PassTime, sep = " ", remove = F) %>% 
    mutate(Date = as.POSIXct(Date)) %>% 
    arrange(Date) %>% 
    #distGeo() provides distance in m, but it is divided by 1000 to calculate distance in Km
    mutate(dist_km =  head(c(NA, round(geosphere::distGeo(.[c("Longitude", "Latitude")])/1000, 
                                 3)), -1),
     #Calculating time difference between rows in seconds
     TimeDiff = as.numeric(Date-lag(Date), units = "secs"),
     #Calculating speed as m per second
     vel_ms = round((dist_km*1000)/TimeDiff, 3),
     #Calculating differences in longitude between points
     DifLon = c(NA, abs(diff(Longitude))))
  
  #Checking if there are any observations with calculated speeds > 7 m/s and removing them
  if(nrow(tagData)!= 0){
    ind <- vector()
    for(j in 2:nrow(tagData)){
      #If the class before observation is better, then removing current observation
      if(tagData$vel_ms[j] > 7 & (tagData$Class[j-1] > tagData$Class[j])){
        ind <- append(ind, j)
        #If the class before observation is worse, then removing previous observation
      }else if(tagData$vel_ms[j] > 7 & (tagData$Class[j-1] < tagData$Class[j])){
        ind <- append(ind, j-1)
        #If the class before observation is the same, then removing current observation
      }else if(tagData$vel_ms[j] > 7 & (tagData$Class[j-1] == tagData$Class[j])){
        ind <- append(ind, j)
      }}}
  #Remove vectors identified above
  if(length(ind) != 0){
    tagData <- tagData[-ind,]}
  #Recalculate velocity and longitude differences
  tagData <- tagData %>% mutate(dist_km =  head(c(NA, round(geosphere::distGeo(.[c("Longitude", "Latitude")])/1000, 
                                     3)), -1),
         #Calculating time difference between rows in seconds
         TimeDiff = as.numeric(Date-lag(Date), units = "secs"),
         #Calculating speed as m per second
         vel_ms = round((dist_km*1000)/TimeDiff, 3),
         #Calculating differences in longitude between points
         DifLon = c(NA, abs(diff(Longitude))))
  #Save clean data as csv file for each tag
  write.csv(tagData, paste0(TagNumb[i], "_clean.csv"), row.names = F)
  #Adding tag data to combined dataframe
  combinedData <- rbind(combinedData, tagData)
  #Delete variables that are no longer needed
  rm(tagData, rawData)}
  
#Remove variables no longer in use
rm(i, ind, j)

#Add a column with lunar phases to the combined dataset with all SPOT information
combinedData <- combinedData %>% mutate(moonPhase = lunar::lunar.phase(Date, name = 8),
                        #Add month column to order data
                        month = factor(months(Date), levels = month.name, ordered = T)) %>% 
  #Reordering months based on the year of tracking
  mutate(month = forcats::fct_reorder(month, lubridate::year(Date)))
#Save data
write_csv(combinedData, "../Outputs/SPOTdata.csv")

#Create new variable calculating the daily change in longitud per tag
DailySum <- combinedData %>% 
  #Select only relevant variables
  select(PTT, Date, PassDate, Longitude, Latitude, vel_ms) %>% 
  #Calculate mean daily velocity per tag
  group_by(PTT, PassDate) %>% 
  mutate(vel_ms = mean(vel_ms, na.rm = T)) %>% 
  #Select the first and last daily location estimate per tag to calculate the total
  #daily change in longitude
  filter(Date == min(Date) | Date == max(Date)) %>% 
  group_by(PTT) %>% 
  mutate(DifLon = c(NA, abs(diff(Longitude)))) %>% 
  #Remove all groupings
  ungroup()

#Creating a shapefile using the above daily summaries
Dailyshp <- DailySum %>% 
  #Assigning WGS84 coordinate reference system
  st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326) 

#Draw map
ggplot()+geom_sf(data = GMRWGS84)+
  geom_sf(data = Dailyshp, aes(colour = PTT))

#Creating a new column GMR which will indicate whether each observation is inside or
#outside the GMR
InOut <- Dailyshp %>% 
  mutate(GMR = case_when(st_intersects(Dailyshp, GMRWGS84, sparse = F) == T ~ "In",
                         st_intersects(Dailyshp, GMRWGS84, sparse = F) == F ~ "Out")) %>% 
  #Change shapefile into data frame
  st_drop_geometry() %>% 
  #Calculate mean daily values per tag
  group_by(PTT, PassDate, GMR) %>% 
  summarise(across(c(vel_ms, DifLon), mean, na.rm = T)) %>% 
  drop_na(DifLon)

write_csv(InOut, "../Outputs/SPOTs_MeanDistLon.csv")
}





# Applying time SSM to SPOT data with foieGras package --------------------
library(foieGras)

#Create empty dataset to store data for all tags
combinedFG <- data.frame()

#Compiling all raw data for all tags
for(i in seq_along(TagData)){
  #Upload file containing raw Argos data
  rawData <- read.csv(TagData[i], stringsAsFactors = F)
  #Extract data needed for analysis
  tagData <- rawData %>% as_tibble() %>% 
    #Remove NA values
    drop_na(Longitude, Latitude) %>% 
    #Keep only variables of interest: Tag ID, date, ARGOS location class (lc), 
    #longitude, latitude
    select(id = PTT,  lc = Class, PassDate, PassTime, lon = Longitude, 
           lat = Latitude) %>% 
    #Change variable type of date column to date
    mutate(PassDate = lubridate::parse_date_time(PassDate, orders = "dmy"),
           id = factor(id)) %>% 
    #Filter out detections after the release date specified in the SuppInfo dataset
    filter(PassDate <= SuppInfo$ReleaseDate[i]) %>% 
    #Unite date and time columns and convert them to date format
    unite(date, PassDate, PassTime, sep = " ", remove = F) %>% 
    mutate(date = as.POSIXct(date)) %>% 
    arrange(date) 
  #Adding tag data to combined dataframe
  combinedFG <- rbind(combinedFG, tagData)
  #Delete variables that are no longer needed
  rm(tagData, rawData)}

#Divide tracks into smaller sections if there are gaps longer than 7 days
combinedFG <- combinedFG %>% group_by(id) %>% 
  mutate(TimeDiff = as.numeric(date-lag(date), units = "days"),
         id = as.character(id)) 
which(combinedFG$TimeDiff > 7)
which(is.na(combinedFG$TimeDiff))

#Manually changing the name of id to divide track into smaller sections
combinedFG$id[145:196] = paste(combinedFG$id[145], 2, sep = "_")
combinedFG$id[708:1136] = paste(combinedFG$id[708], 2, sep = "_")
combinedFG$id[1137:1171] = paste(combinedFG$id[1137], 3, sep = "_")
combinedFG$id[1172:1226] = paste(combinedFG$id[1172], 4, sep = "_")

#Fitting a continuous-time SSM to filter Argos satellite geolocation data
#Maximum travel rate set at 7 m/s
#Correlated random walk model used
#Predictions to be calculated every 12 hours
fit <- combinedFG %>% ungroup() %>% 
  select(id, date, lc, lon, lat) %>%
  fit_ssm(., vmax = 7, model = "crw", time.step = 24)
#Plotting all predictions (lat and lon to appear in the same graph per tag)
plot(fit, what = "predicted", type = 2)

#Get shapefile containing predicted points for all tags
shp <- grab(fit, "predicted") %>% 
  st_transform(crs(GMRWGS84))

#Plotting map
ggplot()+geom_sf(data = GMRWGS84)+
  geom_sf(data = shp, aes(colour = id))

#Classify data points relative to its position to the GMR
PredPts <- shp %>% 
  # mutate(GMR = case_when(st_intersects(., GMRWGS84, sparse = F) == T ~ "In",
  #                        st_intersects(., GMRWGS84, sparse = F) == F ~ "Out")) %>% 
  #Extracting coordinates from geometry to include them dataset
  mutate(lon = st_coordinates(.)[,1], lat = st_coordinates(.)[,2]) %>% 
  #Classifying points based on their location relative to W91.75
  mutate(GMR = case_when(lon <= -91.75 ~ "In",
                         lon > -91.75 ~ "Out")) %>% 
  #Change shapefile into data frame
  st_drop_geometry() %>% 
  select(id, date, GMR, lon, lat) %>% 
  group_by(id) %>% 
  #Calculating differences in longitude between points
  mutate(DifLon = c(NA, abs(diff(lon)))) %>% 
  #Including final point for 198362 - model did not converge
  bind_rows(combinedFG %>% filter(id == "198362_2") %>% 
              filter(date == max(date)) %>%
              distinct(lc, date, .keep_all = T) %>% 
              select(-c("lc", "PassDate", "PassTime", "TimeDiff")) %>% 
              #Removing the track section from the ID
              mutate(GMR = "Out", id = str_extract(id, "[0-9]{6}"))) %>% 
  #Removing the track section from the ID
  mutate(id = str_extract(id, "[0-9]{6}"))


# Testing for differences in means using Analysis of Variance -------------
library(ggpubr)
library(rstatix)

#Getting summary statistics only for groups that have observations inside and outside
#the GMR
PredPts %>% 
  #Grouping by location and tag
  group_by(GMR, id) %>% 
  summarise(N = n()) %>% 
  #Keep tags with 3 or more observations
  filter(N >= 3) %>% 
  print() %>% 
  #Keeping only tags that appear in and out of the GMR for comparison
  group_by(id) %>%
  filter(n() >1) %>%
  #Save the tag numbers to filter dataset
  distinct(id) -> Tags

#Keeping only tags with locations inside and outside the GMR
CompTags <- PredPts %>% 
  filter(id %in% Tags$id)

#Visualise distribution of data points inside tags being compared
CompTags %>%
  ggboxplot(x = "GMR", y = "DifLon", color = "id", add = "point")
CompTags %>%
  ggboxplot(x = "GMR", y = "DifLon", add = "point")

#Identifying outliers
ExtOutliers <- CompTags %>% 
  #Selecting relevant columns
  select(id, DifLon, GMR) %>% 
  #Grouping by location and tag
  group_by(GMR, id) %>% 
  #Identifying outliers per location and tag
  identify_outliers(., DifLon) %>% 
  #Only keeping extreme outliers
  filter(is.extreme == TRUE)

#Checking normality
#QQ plots
CompTags %>% 
  ggqqplot("DifLon", ggtheme = theme_bw())+facet_grid(~GMR)
#Shapiro Wilk's test
CompTags %>% 
  group_by(GMR) %>% 
  shapiro_test(DifLon)

#Removing extreme outliers
NoOut <- CompTags %>% 
  anti_join(ExtOutliers, by = c("DifLon", "GMR")) %>% 
  mutate(id = factor(id),
         GMR = factor(GMR))
NoOut %>% 
  group_by(GMR) %>% 
  shapiro_test(DifLon)
#Visualising data without outliers
NoOut %>% 
  ggqqplot("DifLon", ggtheme = theme_bw())+facet_grid(~GMR)
#Cannot run ANOVA - Too few observations

#PERMANOVA
library(vegan)
#Using data where a tag has at least three observations
NoOut %>% 
  #Grouping by location and tag
  group_by(GMR, id) %>% 
  drop_na() %>% 
  #Testing for differences in location, while controlling for tag (PTT)
  adonis(DifLon ~ GMR, ., permutations = 9999, strata = .$id)

NoOut %>% group_by(GMR) %>% 
  summarise(MeanLonDif = mean(DifLon, na.rm = T))

NoOut %>% distinct(id)

#Save shapefiles individually
for(i in seq_along(unique(PredPts$id))){
  x <- PredPts %>% 
    mutate(moonPhase = lunar::lunar.phase(date, name = 8)) %>% 
    st_as_sf(coords = c("lon", "lat"), crs = 4326) 
  x %>% filter(id == unique(x$id)[i]) %>% 
    st_write(paste0("../Spatial/SPOT_Tags/SSM/T", unique(x$id)[i], ".shp"), 
             delete_layer = T)
}
rm(x)

#Calculating speed for tag 198368
Speed198368 <- PredPts %>% filter(id == 198368) %>% 
  mutate(dist_km =  head(c(NA, round(geosphere::distGeo(.[c("lon", "lat")])/1000, 
                                     3)), -1),
         #Calculating time difference between rows in seconds
         TimeDiff = as.numeric(date-lag(date), units = "secs"),
         #Calculating speed as m per second
         vel_ms = round((dist_km*1000)/TimeDiff, 3))
Speed198368 %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326) %>% 
  st_write("../Spatial/SPOT_Tags/SSM/T198368_Speed.shp", delete_layer = T)

# AMT movement analysis ---------------------------------------------------
{#First we duplicate the GMR layer
GMR <- GMRWGS84

#We then create a variable with the new bounding box values
new_bb <-  c(-100, -6, -76, 14)
names(new_bb) <- c("xmin", "ymin", "xmax", "ymax")
#Assign this variable as a bounding box attribute
attr(new_bb, "class") = "bbox"

#We assign the newly created bounding box to the GMR layer
attr(st_geometry(GMR), "bbox") = new_bb
#Checking bounding box has been correctly assigned
st_bbox(GMR)

#Create a raster using the GMR layer with the new boundaries
rGMR <- raster(GMR, res = 0.1)
rGMR <- fasterize::fasterize(GMR, rGMR)
#Check the raster
plot(rGMR)
rm(GMR)

#Reclassifying into 1 representing the GMR and 0 representing outside areas
#First create reclassification matrix
rclmat <- matrix(c(0.9, 1, 1, NA, NA, 0), ncol = 3, byrow = T)
rcGMR <- reclassify(rGMR, rclmat)
#Check resulting raster
plot(rcGMR)
#Changing name of variable inside the raster to GMR
names(rcGMR) <- "GMR"
rm(rGMR)

#Prepare location data
combSPOT <- CompPts %>% 
  select(x = lon, y = lat, t = date, id = id) %>% 
  nest(-id) %>% mutate(trk = map(data, function(d) {
    amt::make_track(d, x, y, t, crs = sp::CRS("+init=epsg:4326"))
  }))

#Obtaining a summarise of the time distribution intervals between successive locations 
combSPOT %>% mutate(sr = lapply(trk, summarize_sampling_rate)) %>%
  select(id, sr) %>% unnest


m1 <- combSPOT %>%
  mutate(ssf = lapply(trk, function(x) {
    #We have chosen to resample once per day
    x %>% track_resample(rate = hours(24)) %>% 
      filter_min_n_burst(min_n = 3) %>%
      steps_by_burst() %>% 
      random_steps() %>% 
      extract_covariates(rcGMR, where = "start")
    }))

}

