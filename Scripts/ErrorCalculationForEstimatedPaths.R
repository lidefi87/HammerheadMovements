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

### Set directory where all data will be saved
setwd("C:/Users/Denisse/Documents/Manuscripts/miniPAT_SatTags/Spatial/MiniPATs_NewProcess")

# Plotting errors for calculated paths
## for a 90% confidence interval use threshold = 10
## see ?getCtr for description of each input but these are smoother, calc.track, and 
#setup.grid outputs, respectively
bnds <- getCtr(res$s, res$tr, res$g, threshold = 10, makePlot = F)

##this converts those CIs into a lon and lat error metrics
df <- cbind(res$tr, t(data.frame(lapply(bnds, FUN = function(x) c(x$yDist, x$xDist))))) %>% 
  plyr::rename(c("1" = "ydist", "2" = "xdist")) %>% 
  ##some contours will likely not calculate so you can fill this in with the mean error 
  #(or something else) if you wish
  mutate(xdist = case_when(is.na(xdist) ~ mean(xdist, na.rm = T), 
                           TRUE ~ xdist),
         ydist = case_when(is.na(ydist) ~ mean(ydist, na.rm = T), 
                           TRUE ~ ydist))

#save data in csv format for plotting in GIS software
write.csv(df, "csvFiles/157567_idx2_bndNA_par2.25.csv", row.names = F)

### Creating and exporting shapefile
#Creating shapefile with CRS: WGS84
coordinates(df) = ~lon+lat
proj4string(df) = CRS("++proj=longlat +datum=WGS84")
shapefile(df, "157567_idx2_bndNA_par2.25.shp", overwrite = T)

### Calculating distance between points in km
#Creating dataset which includes dates, coordinates and distances for each estimated
#location for that tag 
#distHaversine() provides distance in m, but they are divided by 1000 to present 
#distance in Km
distPts <- cbind(df@data$date, df@coords, dist_km =  c(NA, distHaversine(df)/1000))

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
  geom_point(data = df, aes(x = lon, y = lat))

