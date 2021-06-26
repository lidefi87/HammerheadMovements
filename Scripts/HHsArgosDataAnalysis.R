# Title: Analysing SPOT data from Hammerhead sharks (S. lewini) in the TEP
# Author: Denisse Fierro Arcos
# Version: 5
# Date last updated: 2021-06-26

# Clearing global environment
rm(list = ls())


###########################################################################
# Loading relevant libraries ----------------------------------------------
# library(raster)
library(sf)
# library(rgeos)
library(tidyverse)
library(lubridate)
library(foieGras)


###########################################################################
# Loading relevant shapefiles ---------------------------------------------
# Uploading shapefiles needed for analysis - Galapagos Islands and GMR
# GalIslWGS84 <-  shapefile("Data/Layers/galapagos_island_wgs84.shp")
GMRWGS84 <- read_sf("Data/Layers/ReservaMarinaWGS84.shp")


###########################################################################
# Satellite data preparation ----------------------------------------------
# Uploading supporting info for tag deployments
SuppInfo <- openxlsx::read.xlsx("Data/SPOT_Tags/SPOT6TagList_Deployments.xlsx", sheet = "ReleaseDates", 
                                detectDates = T) %>% 
  #Keeping only animals whose remained tagged for at least one week
  filter(DaysTagged >= 7)

#Getting a list of the directories within the folder - Keep files with tag numbers only
TagNumb <- dir(path = "Data/SPOT_Tags/", pattern = "[0-9]{6}[A-Z]{0,1}$", full.names = T)
TagData <- dir(path = "Data/SPOT_Tags/", pattern = "-RawArgos.csv$", full.names = T, 
               recursive = T)

#Keeping tags only the tags with over a week of data
overOneWeek <- str_match(TagNumb, "[0-9]{6}[A-Z]{0,1}") %in% SuppInfo$ID
TagNumb <- TagNumb[overOneWeek]
TagData <- TagData[overOneWeek]
#Removing unused variables
rm(overOneWeek)

#Create empty list to store data for all tags
tagData <- list()

#Cleaning data within each folder
for(i in seq_along(TagData)){
  #Upload file containing raw Argos data
  tagData[[i]] <- read_csv(TagData[i]) %>%
    #Removing rows with no longitude and/or latitude information
    drop_na(c(Longitude, Latitude)) %>% 
    #Keep only variables of interest: Tag ID, date, ARGOS location class (lc), longitude, latitude
    select(PTT,  Class, PassDate, PassTime, Longitude, Latitude) %>% 
    #Eliminate duplicate entries for the same date and time
    distinct() %>%
    #Unite data and time prior to transforming to time date format
    unite(date, PassDate, PassTime, sep = " ") %>% 
    #Changing data to correct format
    mutate(date = dmy_hms(date),
           PTT = factor(PTT)) %>% 
    #Remove any observations after the release date
    filter(date <= SuppInfo$ReleaseDate[i]) %>% 
    #Ensure data is arranged chronologically
    arrange(date) %>% 
    #distGeo() provides distance in m, but it is divided by 1000 to calculate distance in Km
    mutate(dist_km =  head(c(NA, round(geosphere::distGeo(.[c("Longitude", "Latitude")])/1000, 3)), -1),
         #Calculating time difference between rows in seconds
         TimeDiff = as.numeric(date-lag(date), units = "secs"),
         #Calculating speed as m per second
         vel_ms = round((dist_km*1000)/TimeDiff, 3),
         #Calculating differences in longitude between points
         DifLon = c(NA, abs(diff(Longitude))),
         DifLat = c(NA, abs(diff(Latitude))))
  }

#Transform list to dataframe and prepare prior to application of SSM
tagData <- tagData %>% 
  bind_rows() %>% 
  relocate(Class, .after = date) %>% 
  rename(id = PTT, lc = Class, lon = Longitude, lat = Latitude) %>% 
  #Removing outliers too far west and on land
  filter(lon > -100 & lon < -75)

#Saving clean tagging data
write_csv(tagData, "../Outputs/SPOTdata.csv")


###########################################################################
# Fitting a continuous-time state-space model to filter Argos data --------
#Reconstruct track with one point every 12 hours
fit <- tagData %>% 
  select(id:lat) %>% 
  fit_ssm(vmax = 7, model = "crw", time.step = 12, control = ssm_control(verbose = 0))

#Plotting reconstructed tracks
plot(fit, type = 2, what = "predicted", ask = F)


###########################################################################
# Fitting a move persitence model to estimate behaviour -------------------
# Analysis based on reconstructed tracks
fmp <- fit_mpm(fit, what = "predicted", model = "mpm")

# Plotting persistence of movements
plot(fmp, ncol = 4, pal = "Zissou1", rev = T, ask = F)

# Plotting reconstructed tracks with persistence given colour of points
fmap(fit, fmp, what = "predicted")


###########################################################################
# Saving reconstructed tracks with movement persistence data --------------
#Extracting data and including this in a list
routes <- list(predicted = grab(fit, "predicted", as_sf = F),
               behaviours = grab(fmp, what = "fitted", as_sf = F))

#Match tracks with behaviour based on the id and date columns
routes <- routes %>% 
  reduce(full_join, by = c("id", "date"))

#Saving results
write_csv(routes, "../Outputs/Reconstructed_SPOT_tracks.csv")


###########################################################################
# Linearity calculations --------------------------------------------------
# Calculating distance from start to end point
linearity <- list (dist_st_end = routes %>% 
                     #Grouping by tag ID
                     group_by(id) %>% 
                     #Selecting relevant columns only
                     select(id:lat) %>% 
                     #Selecting the first and last location
                     filter(date == min(date) | date == max(date)) %>% 
                     #Splitting data per tag ID prior to distance calculation
                     split(.$id) %>% 
                     #Calculating distance bertween start and end point in Km 
                     map_dbl(~head(geosphere::distGeo(.[c("lon", "lat")])/1000, 1)) %>%
                     #Saving result as data frame
                     data.frame(start_end_dist = .) %>% 
                     rownames_to_column("id"),
                   dist_tot =  routes %>% 
                     group_by(id) %>% 
                     #Selecting relevant columns only
                     select(id:lat) %>% 
                     #Splitting data per tag ID prior to distance calculation
                     split(.$id) %>% 
                     #Calculating distance bertween start and end point in Km 
                     map(~head(c(NA, geosphere::distGeo(.[c("lon", "lat")])/1000), -1)) %>% 
                     map(sum, na.rm = T) %>%
                     unlist() %>% 
                     data.frame(tot_dist = .) %>% 
                     rownames_to_column("id")) %>% 
  #Join distance from start to end point and total distance
  reduce(full_join, by = "id") %>% 
  #Calculate linearity
  mutate(linearity = start_end_dist/tot_dist)

#Saving results
write_csv(linearity, "../Outputs/Linearity_Reconstructed_SPOT_tracks.csv")



# Plotting reconstructed routes -------------------------------------------
# Creating a shapefile with reconstructed tracks
routes_shp <- routes %>% 
  st_as_sf(coords = c("lon", "lat"), crs = 4326)

# Accesssing world map
library(rnaturalearth)
world <- ne_countries(scale = "medium", returnclass = "sf")

# Plotting using ggplot
ggplot(data = routes_shp)+
  geom_sf(aes(col = g))+
  geom_sf(data = world)+
  scale_color_distiller(type = "div")+
  coord_sf(xlim = c(-100, -75), ylim = c(-2, 16), expand = FALSE)+
  facet_wrap(~id)
ggsave("../Outputs/Maps/Reconstructed_Facets.pdf")
#   
# # Testing for differences in means using Analysis of Variance -------------
# library(ggpubr)
# library(rstatix)
# 
# #Getting summary statistics only for groups that have observations inside and outside
# #the GMR
# PredPts %>% 
#   #Grouping by location and tag
#   group_by(GMR, id) %>% 
#   summarise(N = n()) %>% 
#   #Keep tags with 3 or more observations
#   filter(N >= 3) %>% 
#   print() %>% 
#   #Keeping only tags that appear in and out of the GMR for comparison
#   group_by(id) %>%
#   filter(n() >1) %>%
#   #Save the tag numbers to filter dataset
#   distinct(id) -> Tags
# 
# #Keeping only tags with locations inside and outside the GMR
# CompTags <- PredPts %>% 
#   filter(id %in% Tags$id)
# 
# #Visualise distribution of data points inside tags being compared
# CompTags %>%
#   ggboxplot(x = "GMR", y = "DifLon", color = "id", add = "point")
# CompTags %>%
#   ggboxplot(x = "GMR", y = "DifLon", add = "point")
# 
# #Identifying outliers
# ExtOutliers <- CompTags %>% 
#   #Selecting relevant columns
#   select(id, DifLon, GMR) %>% 
#   #Grouping by location and tag
#   group_by(GMR, id) %>% 
#   #Identifying outliers per location and tag
#   identify_outliers(., DifLon) %>% 
#   #Only keeping extreme outliers
#   filter(is.extreme == TRUE)
# 
# #Checking normality
# #QQ plots
# CompTags %>% 
#   ggqqplot("DifLon", ggtheme = theme_bw())+facet_grid(~GMR)
# #Shapiro Wilk's test
# CompTags %>% 
#   group_by(GMR) %>% 
#   shapiro_test(DifLon)
# 
# #Removing extreme outliers
# NoOut <- CompTags %>% 
#   anti_join(ExtOutliers, by = c("DifLon", "GMR")) %>% 
#   mutate(id = factor(id),
#          GMR = factor(GMR))
# NoOut %>% 
#   group_by(GMR) %>% 
#   shapiro_test(DifLon)
# #Visualising data without outliers
# NoOut %>% 
#   ggqqplot("DifLon", ggtheme = theme_bw())+facet_grid(~GMR)
# #Cannot run ANOVA - Too few observations
# 
# #PERMANOVA
# library(vegan)
# #Using data where a tag has at least three observations
# NoOut %>% 
#   #Grouping by location and tag
#   group_by(GMR, id) %>% 
#   drop_na() %>% 
#   #Testing for differences in location, while controlling for tag (PTT)
#   adonis(DifLon ~ GMR, ., permutations = 9999, strata = .$id)
# 
# NoOut %>% group_by(GMR) %>% 
#   summarise(MeanLonDif = mean(DifLon, na.rm = T))
# 
# NoOut %>% distinct(id)
# 
# #Save shapefiles individually
# for(i in seq_along(unique(PredPts$id))){
#   x <- PredPts %>% 
#     mutate(moonPhase = lunar::lunar.phase(date, name = 8)) %>% 
#     st_as_sf(coords = c("lon", "lat"), crs = 4326) 
#   x %>% filter(id == unique(x$id)[i]) %>% 
#     st_write(paste0("../Spatial/SPOT_Tags/SSM/T", unique(x$id)[i], ".shp"), 
#              delete_layer = T)
# }
# rm(x)
# 
# #Calculating speed for tag 198368
# Speed198368 <- PredPts %>% filter(id == 198368) %>% 
#   mutate(dist_km =  head(c(NA, round(geosphere::distGeo(.[c("lon", "lat")])/1000, 
#                                      3)), -1),
#          #Calculating time difference between rows in seconds
#          TimeDiff = as.numeric(date-lag(date), units = "secs"),
#          #Calculating speed as m per second
#          vel_ms = round((dist_km*1000)/TimeDiff, 3))
# Speed198368 %>%
#   st_as_sf(coords = c("lon", "lat"), crs = 4326) %>% 
#   st_write("../Spatial/SPOT_Tags/SSM/T198368_Speed.shp", delete_layer = T)
# 
# # AMT movement analysis ---------------------------------------------------
# {#First we duplicate the GMR layer
# GMR <- GMRWGS84
# 
# #We then create a variable with the new bounding box values
# new_bb <-  c(-100, -6, -76, 14)
# names(new_bb) <- c("xmin", "ymin", "xmax", "ymax")
# #Assign this variable as a bounding box attribute
# attr(new_bb, "class") = "bbox"
# 
# #We assign the newly created bounding box to the GMR layer
# attr(st_geometry(GMR), "bbox") = new_bb
# #Checking bounding box has been correctly assigned
# st_bbox(GMR)
# 
# #Create a raster using the GMR layer with the new boundaries
# rGMR <- raster(GMR, res = 0.1)
# rGMR <- fasterize::fasterize(GMR, rGMR)
# #Check the raster
# plot(rGMR)
# rm(GMR)
# 
# #Reclassifying into 1 representing the GMR and 0 representing outside areas
# #First create reclassification matrix
# rclmat <- matrix(c(0.9, 1, 1, NA, NA, 0), ncol = 3, byrow = T)
# rcGMR <- reclassify(rGMR, rclmat)
# #Check resulting raster
# plot(rcGMR)
# #Changing name of variable inside the raster to GMR
# names(rcGMR) <- "GMR"
# rm(rGMR)
# 
# #Prepare location data
# combSPOT <- CompPts %>% 
#   select(x = lon, y = lat, t = date, id = id) %>% 
#   nest(-id) %>% mutate(trk = map(data, function(d) {
#     amt::make_track(d, x, y, t, crs = sp::CRS("+init=epsg:4326"))
#   }))
# 
# #Obtaining a summarise of the time distribution intervals between successive locations 
# combSPOT %>% mutate(sr = lapply(trk, summarize_sampling_rate)) %>%
#   select(id, sr) %>% unnest
# 
# 
# m1 <- combSPOT %>%
#   mutate(ssf = lapply(trk, function(x) {
#     #We have chosen to resample once per day
#     x %>% track_resample(rate = hours(24)) %>% 
#       filter_min_n_burst(min_n = 3) %>%
#       steps_by_burst() %>% 
#       random_steps() %>% 
#       extract_covariates(rcGMR, where = "start")
#     }))
# 
# }

