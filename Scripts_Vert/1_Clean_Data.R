# Analysis of Hammerhead sharks (Sphyrna lewini) vertical movements
# Step 1: Accessing and cleaning the data
# Date of latest update: 2021-01-08
# Version: 1
# Function: Acessing and cleaning vertical movement data of S. lewini. It will also
# save some density distribution plots
# Prepared for the Sharks Ecology Project of the Charles Darwin Foundation (CDF)
# Script related to publication entitled "xxxx" in xxxx magazine.

# Uploading relevant libraries --------------------------------------------
library(tidyverse)
library(stringr)
library(suncalc)
library(ggpubr)
library(factoextra)
library(ggsci)
library(diveMove)

# Accessing data ----------------------------------------------------------
#Getting full paths for diving data - A total of 11 tags have useful information
VertMov <- list.files(path = "../../HH_Movements/Data/MiniPAT_Tags/", full.names = T,
                      pattern = "-Series.csv$", recursive = T)
#Getting list of tag numbers from list of depth data 
TagID <- str_extract(VertMov, pattern =  "[0-9]{6}")
#Getting full paths for horizontal movements
HorMov <- list.files(path = "../../Spatial/MiniPATs_NewProcess/csvFiles/", full.names = T)

# Merging horizontal and vertical data ------------------------------------
#Initialising empty dataframe to store all movement data for all tags
Movs <- data.frame()
plot_list <- list()
TDRs <- list()
TDRs_Calib <- list()

#Creating loop to include all 11 tags used in the study
for(i in seq_along(TagID)){
  #Reading vertical data for each tag
  Vert <- read_csv(VertMov[which(str_detect(VertMov, TagID[i]))], ) %>% 
    #Selecting only the columns relevant to this process 
    select(Ptt, Day, Time, Depth, DRange, Temperature, TRange) %>% 
    #Changing day column to Date format to match horizontal database
    mutate(Day = parse_date(Day, format = "%d-%b-%Y")) %>% 
    #Changing the name of the column to match horizontal data
    rename("date" = "Day") %>% 
    #Removing rows with no depth info as it cannot be used in further analyses,
    drop_na(Depth)
  
  #Histograms of frequency distribution of depth values per tag
  p <- Vert %>% ggplot(aes(x = Depth))+
    #Histogram using density values
    geom_histogram(binwidth = 100, aes(y = ..density.., fill = ..density..), size = 0.1, 
                   colour = "grey20")+
    #Labelling each plot with the tag name it comes from
    labs(title = TagID[i])+
    #Creating a scale color for the mean and median values to be shown in the plot
    scale_color_manual(name = "Statistics", values = c(Median = "purple", Mean = "green4"),
                       guide = guide_legend(title.position = "top", direction = "horizontal"))+
    #Changing the color palette for the density values
    scale_fill_distiller(palette = "YlOrRd", 
                         #Change breaks, so only two appear to avoid overlap in the labels
                         breaks = c(0.002, 0.006))+
    #Customising the colorbar
    guides(fill = guide_colorbar(title = "Density", title.position = "top", 
                                 direction = "horizontal",
                                 ticks.colour = "black", frame.colour = "black"))+
    #Plotting the density distribution
    stat_function(fun = dnorm, args = list(mean = mean(Vert$Depth), sd = sd(Vert$Depth)), 
                  lwd = 0.2,  color = "black")+
    #Include a vertical line showing the mean and median depth for that particular tag
    geom_vline(aes(color = "Mean", xintercept = mean(Depth)), lty = 4)+
    geom_vline(aes(color = "Median", xintercept = median(Depth)), lty = 2)+
    theme_bw()+
    #Changing the font and size for all text in the plot
    theme(plot.title = element_text(hjust = 0.5, family = "sans", size = 12),
          axis.text = element_text(family = "sans", size = 12),
          legend.text = element_text(family = "sans", size = 12),
          legend.title = element_text(family = "sans", size = 12),
          #Removing the axes titles
          axis.title = element_blank())
  #Saving plots under one list ready to create a composite plot
  plot_list[[i]] <- p
  
  #Reading horizontal data for each tag
  Hor <- read_csv(HorMov[which(str_detect(HorMov, TagID[i]))]) %>% 
    #Selecting only relevant columns
    select(date, lon, lat)
  
  #Sunset and sunrise data
  Hor <- Hor %>% 
    #Obtaining time zone based on lat and lon information
    mutate(tz = lutz::tz_lookup_coords(lat, lon, warn = F),
           #When GMT is given as timezone it has the wrong symbol, fixing this manually
           tz = stringr::str_replace(tz, "\\+", "-")) 
  ##Obtaining sunrise and sunset times for daily locations
  y <- data.frame()
  for(j in 1:nrow(Hor)){
    #One day added to date to correct reported error in suncalc
    z <- getSunlightTimes(Hor$date[j]+1, Hor$lat[j], Hor$lon[j], tz = Hor$tz[j], 
                          keep = c("sunrise", "sunset"))
    y <- rbind(y, z)}
  #Add sunrise and sunset times to Hor data frame
  Hor <- Hor %>% 
    right_join(y %>% select(-date), by = c("lat", "lon"))
  
  #Removing unused variables
  rm(y, z)
  
  #Merging vertical and horizontal data 
  x <- Vert %>% left_join(Hor, by = c("date")) %>% 
    janitor::clean_names() %>% 
    #Removing rows with no longitude, as this means the tag was no longer on the shark
    drop_na(lon) %>% 
    #Merging date and time into one column
    unite(dateTime, date, time, sep = " ", remove = F) %>% 
    #Changing dateTime and date columns into Date objects
    mutate(dateTime = parse_datetime(dateTime),
           date = parsedate::parse_date(date),
           #Calculating time difference (in seconds) between each observation
           timeDiff = as.numeric(dateTime-lag(dateTime), units = "secs"),
           ptt = factor(ptt))
  
  #Merge all data from all tags into one data frame
  Movs <- rbind(Movs, x)
  
  #DiveMove analysis
  TDRs <- append(TDRs, createTDR(time = x$dateTime, depth = x$depth, 
                                 concurrentData = x[,c(7,9,10)], file = VertMov[i]))
  TDRs_Calib <- append(TDRs_Calib, calibrateDepth(TDRs[[i]], dive.thr = 5, zoc.method = "offset", 
                                                  offset = 0.5, dive.model = "smooth.spline", 
                                                  smooth.par = NULL, knot.factor = 20))
}

#Removing variables no longer needed
rm(Vert, Hor, x, p, i, j)

# Classifying data --------------------------------------------------------
#Classiying depths into predefined bins
breaksDepth <- c(0, 10, 25, 50, 200, 400, 1000, max(Movs$depth))
tagsDepth <- c("<10", "10-25", "25-50", "50-200", "200-400", "400-1000", ">1000")
#Classifying temperatures into predefined bins
breaksTemp <- c(0, seq(from = 10, to = 26, by = 2), max(Movs$temperature, na.rm = T))
tagsTemp <- c("<10", "10-12", "12-14", "14-16", "16-18", "18-20", "20-22", "22-24", 
              "24-26", ">26")

#Reclassifying depth and temperature
Movs <- Movs %>% 
  #Correcting time zone of dateTime column
  mutate(dateTime = timechange::time_force_tz(dateTime, tz = tz),
         DepthBin = cut(depth, breaks = breaksDepth, 
                        include.lowest = T, right = F, labels = tagsDepth),
         TempBin = cut(temperature, breaks = breaksTemp, 
                       include.lowest = T, right = F, labels = tagsTemp)) %>% 
  #Changing newly created columns into ordered factors
  mutate(DepthBin = factor(DepthBin, levels = tagsDepth, ordered = T),
         TempBin = factor(TempBin, levels = tagsTemp, ordered = T)) %>% 
  #Classifying whether observations occurred at night or during the morning
  mutate(diel = case_when(dateTime < sunrise ~ "Night", 
                          dateTime >= sunrise & dateTime < sunset ~ "Morning",
                          dateTime >= sunset ~ "Night"))

#Removing variables no longer in use
rm(breaksDepth, breaksTemp, tagsDepth, tagsTemp)
