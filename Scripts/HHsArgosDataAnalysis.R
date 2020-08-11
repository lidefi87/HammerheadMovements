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

#Variable containing location classes related to data quality (Classes: G > 3 > 2 > 1 > 0 > A > B > Z)
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
  tagData <- tagData[-ind,]
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
rm(classInt, i, TagData, TagNumb, SuppInfo, ind, j)

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


# Testing for differences in means using Analysis of Variance -------------
library(ggpubr)
library(rstatix)

#Getting summary statistics only for groups that have observations inside and outside
#the GMR
InOut %>% 
  #Grouping by location and tag
  group_by(GMR, PTT) %>% 
  #Calculating summary stats (mean and SD)
  get_summary_stats(DifLon, type = "mean_sd") %>% 
  #Keep tags with 3 or more observations
  filter(n >= 3) %>% 
  print() %>% 
  #Keeping only tags that appear in and out of the GMR for comparison
  group_by(PTT) %>%
  filter(n() >1) %>%
  #Save the tag numbers to filter dataset
  distinct(PTT) -> Tags

#Keep tags previously identified
CompTags <- InOut %>% 
  #Kepp only tags previously identified as suitable for comparison
  filter(PTT %in% Tags$PTT)
  
#Visualise distribution of data points inside tags being compared
CompTags %>%
  ggboxplot(x = "GMR", y = "DifLon", color = "PTT", add = "point")
CompTags %>%
  ggboxplot(x = "GMR", y = "DifLon", add = "point")

#Identifying outliers
ExtOutliers <- CompTags %>% 
  #Selecting relevant columns
  select(PTT, DifLon, GMR) %>% 
  #Grouping by location and tag
  group_by(GMR) %>% 
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
  mutate(PTT = factor(PTT),
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
InOut %>% 
  #Grouping by location and tag
  group_by(GMR, PTT) %>% 
  #Keep tags with 3 or more observations
  filter(n() >= 3) %>%
  #Testing for differences in location, while controlling for tag (PTT)
  adonis(DifLon ~ GMR, ., permutations = 9999, strata = .$PTT)





#All locations are in water - Section commented out
######################## CREATING POINT SHAPEFILE WITH CLEAN DATA FOR SSM ##############################
#Creating shapefile that only includes water points using the Galapagos layer
#Create a vector of unique datapoints with CRS: WGS84
# combinedDataLyr <-  combinedData
# coordinates(combinedDataLyr) <-  ~Longitude+Latitude
# proj4string(combinedDataLyr) <-  CRS("++proj=longlat +datum=WGS84")
# crs(combinedDataLyr) <-  crs(GalIslWGS84)
#
# #Find points that fall on land using the island layer and delete them from the main database
# x = over(combinedDataLyr, GalIslWGS84)
# #Use previous data frame to find water only points (i.e., rows with NA values) in combined dataframe
# combinedData = combinedData[which(is.na(x$id)),]
# #Remove variables no longer needed
# rm(x)
# #Update point shapefile
# combinedDataLyr = combinedData
# coordinates(combinedDataLyr) = ~Longitude+Latitude
# proj4string(combinedDataLyr) = CRS("++proj=longlat +datum=WGS84")
# crs(combinedDataLyr) = crs(GalIslWGS84)
# 
# #Find points that fall within GMR borders and create a data frame
# filteredPts = combinedDataLyr[GMRWGS84,]
# rm(combinedDataLyr)
# 
# #Plotting maps to check if any points are located on land
# plot(filteredPts)
# plot(GalIslWGS84, add = T)
# plot(GMRWGS84, add = T)
#Zoom option can also be used as follows:
# zoom(GalIslWGS84, ext = extent(combinedData))
# plot(combinedData, add = T)

# #Saving final output as shapefile
# writeOGR(filteredPts, dsn = "C:/Users/Denisse/Documents/TigerSharksSatelliteDataKernels/Layers",
#          layer = "cleanCombinedData", driver="ESRI Shapefile", overwrite_layer = T)
#
#
################################# BAYESIAN STATE-SPACE MODEL (SSM) #####################################
#Section not applicable - commented out
# #Extracting filtered points to dataframe
# combinedData = cbind(filteredPts@data, filteredPts@coords)

#Rename dataset columns as needed for SSM application
colnames(combinedData) <- c("id", "lc", "date", "lon", "lat")
# rm(filteredPts)

#Plot tracks with coordinates and varying colour based on date
combinedData %>% ggplot(aes(lat, lon, col = date))+geom_path()+geom_point()+facet_wrap(~id)

#Calculate number of days between observations initiating from the same tag
combinedData <- combinedData %>% group_by(id) %>% mutate(days = (date-lag(date))/86400) %>% 
  #Find gaps in observations of seven days or longer
  mutate(split = case_when(days >=7, ~ "split"))

#Find gaps in observations of seven days or longer
which(combinedData$days >= 7)

#Keep tags with over 60 observations
#Get number of observations per tag number
combinedData60 <- combinedData %>% right_join(combinedData %>% group_by(id) %>% 
                                              summarise(Obs = n()) %>% filter(Obs > 60) %>% 
                                              select(id), by = "id")

#Fit state-space model
fit <- fit_ssm(combinedData, model = "hDCRWS", tstep = 0.5, adapt = 10000)

#Checking fit of model
diag_ssm(fit)
map_ssm(fit)
plot_fit(fit)
dev.off()

#Extracting modelled values and creating a point shapefile
result = get_summary(fit)
coordinates(result) = ~lon+lat
proj4string(result) = CRS("++proj=longlat +datum=WGS84")
crs(result) = crs(GalIslWGS84)
rm(combinedData)
#Find points that fall on land using the island layer and delete them from the main database
x = over(result, GalIslWGS84)
#Use previous data frame to find water only points (i.e., rows with NA values) in combined dataframe
cleanData = fit$summary[which(is.na(x$id)),]
#Remove variables no longer needed
rm(x)
#Update point shapefile
coordinates(cleanData) = ~lon+lat
proj4string(cleanData) = CRS("++proj=longlat +datum=WGS84")
crs(cleanData) = crs(GalIslWGS84)
rm(result)
# #Save result as a shapefile
# shapefile(cleanData, "C:/Users/Denisse/Documents/TigerSharksSatelliteDataKernels/Layers/SSMresults.shp",
#           overwrite = T)


############################# KERNEL DENSITY ESTIMATION USING SSM DATA #################################
library(aspace)
library(spatialEco)

#Calculate standard distance (SDD) to estimate bandwidth (h)
coords = data.frame(lat = cleanData$lat, lon = cleanData$lon)
calc_sdd(points = coords) #Calculate SDD using coordinates
#Plot SDD with Galapagos map as background to visually check results
plot_sdd(plotnew = T)
plot(GalIslWGS84, add = T)
#Save SDD calculation as variable to be used in bandwidth estimation
SDDpts = r.SDD$SDD
#Using Silverman rule of thumb: h = 1.06*SD*n^-0.2 to calculate bandwith
h = 1.06*SDDpts*(nrow(coords)^(-0.2))
#Remove variables no longer in use
rm(coords, r.SDD, SDDpts, sddatt, sddloc)

#Unweighted KDE calculation
kde = sp.kde(x = cleanData, bw = h, n = 10000)
#Remove any land areas from KDE raster
kdeClip = raster::mask(kde, GalIslWGS84, inverse = T)
rm(kde)

############################ CALCULATING PERCENTAGE VOLUME CONTOURS ###################################
#Create percentage volume contours (50%, 75%, 95%)
p95 = raster.vol(kdeClip, p = 0.95)
p75 = raster.vol(kdeClip, p = 0.75)
p50 = raster.vol(kdeClip, p = 0.50)
#Merge all contours together in one raster
PVC = p95+p75+p50
rm(p95, p75, p50)
#Change zero values into NA - Reclassify all other values
PVC[PVC == 0] = NA
PVC[PVC == 3] = 0.5
PVC[PVC == 2] = 0.75
PVC[PVC == 1] = 0.95
#Plot final raster with heat colour ramp
plot(PVC, col = heat.colors(3))
plot(GalIslWGS84, add = T)

#Calculate area in km^2 covered by each contour
tapply(area(PVC), PVC[], sum)

# #Save final raster as a TIF file
# writeRaster(PVC, filename = "C:/Users/Denisse/Documents/TigerSharksSatelliteDataKernels/Layers/KDEresult.tif",
#             format = "GTiff", overwrite = T)

#######################################################################################################

#Calculate proportion of SSM locations within 5 and 10 km from turtle nesting beaches
#Import locations of turtle nesting beaches
Tort = shapefile("C:/Users/Denisse/Documents/TigerSharksSatelliteDataKernels/Layers/PlayasTortugas.shp")
#Create buffers (5 and 10 km) using turtle nesting beaches (unit is in meters)
x = buffer(Tort, width = 5000)
y = buffer(Tort, width = 10000)
#Clip out land areas using Galapagos layer
Tort_5km = erase(x, GalIslWGS84)
Tort_10km = erase(y, GalIslWGS84)
#Remove unused variables
rm(Tort, x, y)

#Extract points within buffers
Pts_5km = intersect(cleanData, Tort_5km)
Pts_10km = intersect(cleanData, Tort_10km)

#Calculate proportion of points within each buffer area
prop5km = nrow(Pts_5km)/nrow(cleanData)
prop10km = nrow(Pts_10km)/nrow(cleanData)

# #Save all shapefiles created
# shapefile(Pts_5km, filename = "C:/Users/Denisse/Documents/TigerSharksSatelliteDataKernels/Layers/5kmPts.shp",
#           overwrite = T)
# shapefile(Pts_10km, filename = "C:/Users/Denisse/Documents/TigerSharksSatelliteDataKernels/Layers/10kmPts.shp",
#           overwrite = T)
# shapefile(Tort_5km, filename = "C:/Users/Denisse/Documents/TigerSharksSatelliteDataKernels/Layers/5kmBuff.shp",
#           overwrite = T)
# shapefile(Tort_10km, filename = "C:/Users/Denisse/Documents/TigerSharksSatelliteDataKernels/Layers/10kmBuff.shp",
#           overwrite = T)
# shapefile(result, filename = "C:/Users/Denisse/Documents/TigerSharksSatelliteDataKernels/Layers/SSMresults.shp", 
          # overwrite = T)
