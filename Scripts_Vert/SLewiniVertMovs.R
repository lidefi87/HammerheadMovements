# Analysis of Hammerhead sharks (Sphyrna lewini) vertical movements -------
# Date of latest update: 2021-01-04
# Version: 1
# Prepared for the Sharks Ecology Project of the Charles Darwin Foundation (CDF)
# Script related to publication entitled "xxxx" in xxxx magazine.


# Uploading relevant libraries --------------------------------------------
{library(tidyverse)
library(stringr)
library(suncalc)
library(ggpubr)
library(factoextra)
library(ggsci)
library(diveMove)}

# Accessing data ----------------------------------------------------------
#Getting full paths for diving data - A total of 11 tags have useful information
VertMov <- list.files(path = "../../HH_Movements/Data/MiniPAT_Tags/", full.names = T,
                   pattern = "-Series.csv$", recursive = T)
#Getting list of tag numbers from list of depth data 
TagID <- str_extract(VertMov, pattern =  "[0-9]{6}")
#Getting full paths for horizontal movements
HorMov <- list.files(path = "../../Spatial/MiniPATs_NewProcess/csvFiles/", full.names = T)


# Getting metadata about tagged sharks ------------------------------------
miniPat <- readxl::read_excel("../../SupportingData/HH_MiniPats_2016_2019.xlsx", 
                              sheet = "MiniPATS") %>% 
  janitor::clean_names()
#Calculating summary statistics
miniPat %>% 
  #drop rows with no data
  drop_na(days_at_liberty) %>% 
  #calculating summary statistics
  rstatix::get_summary_stats(estimated_total_length_m, 
                             days_at_liberty, 
                             distance_covered_km, 
                             type = "common")
#Calculating distance per day
miniPat %>% 
  #drop rows with no data
  drop_na(days_at_liberty) %>% 
  select(tag_id, distance_covered_km, days_at_liberty) %>% 
  mutate(dist_day = distance_covered_km/days_at_liberty) %>% 
  rstatix::get_summary_stats(dist_day, type = "common")

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

#Making composite figure with depth histograms for each tag
#Extracting legend from first figure in the list
leg <- ggpubr::get_legend(plot_list[[1]])
#Including the extracted legend at the end of the plot list as a plot prior to creating 
#composite plot
plot_list[[i+1]] <- as_ggplot(leg)
#Creating composite image. Legend turned off because it is being included at the end
g <- ggarrange(plotlist = plot_list, ncol = 4, nrow = 3, legend = F)
#Annotate figure with the axis labels
g <- annotate_figure(g, bottom = "Depth (m)", left = "Density")
#Saving composite image to disk
ggsave("Scripts_Vert/Figures/DepthHistograms.tiff", g, device = "tiff", dpi = 400, width = 25, 
       height = 23.44, units = "cm")

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


# Summarising movement data -----------------------------------------------
Movs %>% group_by(ptt, diel) %>% 
  rstatix::get_summary_stats(depth, temperature, type = "mean_se")

#Obtaining summary for depth and temperature per tag
Movs %>% select(ptt, diel, date, depth, temperature) %>% 
  group_by(ptt, diel, date) %>% 
  mutate(minDepth = min(depth), 
         maxDepth = max(depth)) %>% 
  group_by(ptt, diel) %>% 
  rstatix::get_summary_stats(minDepth, maxDepth, temperature, type = "mean_se") %>% 
  select(-n) %>% 
  pivot_wider(names_from = variable, values_from = c(mean, se)) %>% 
  unite("MaxDepth", mean_maxDepth, se_maxDepth, sep = " ± ") %>% 
  unite("MinDepth", mean_minDepth, se_minDepth, sep = " ± ") %>% 
  unite("Temp", mean_temperature, se_temperature, sep = " ± ") %>% 
  pivot_wider(names_from = diel, values_from = c(MaxDepth, MinDepth, Temp)) %>% 
  write.csv("Scripts_Vert/Outputs/minMaxDepth.csv", row.names = F)


#Testing for differences in mean temp, min and max depths between morning and night
#Extracting data
SumDepTemp <- Movs %>% select(ptt, diel, date, depth, temperature) %>% 
  group_by(ptt, diel, date) %>% 
  summarise(minDepth = min(depth), 
         maxDepth = max(depth),
         meanTemp = mean(temperature)) %>% 
  mutate(year = lubridate::year(date),
         month = lubridate::month(date))

#Univariate PERMANOVA
#Minimum depth
#Differences between day and night - Non significant (p = 0.342)
vegan::adonis(minDepth ~ diel, data = SumDepTemp, permutations = 999, 
              method = "euclidean")
#Differences among animals - Non significant (p = 0.342)
vegan::adonis(minDepth ~ ptt, data = SumDepTemp, permutations = 999, 
              method = "euclidean")
#Differences across years - Non significant (p = 0.209)
vegan::adonis(minDepth ~ year, data = SumDepTemp, permutations = 999, 
              method = "euclidean")
#Differences across months - Non significant (p = 0.974)
vegan::adonis(minDepth ~ month, data = SumDepTemp, permutations = 999, 
              method = "euclidean")

#Maximum depth
#Differences between day and night - Significant (p = 0.001)
vegan::adonis(maxDepth ~ diel, data = SumDepTemp, permutations = 999, 
              method = "euclidean")
#Differences among animals - Significant (p = 0.001)
vegan::adonis(maxDepth ~ ptt, data = SumDepTemp, permutations = 999, 
              method = "euclidean")
#Differences across years - Non significant (p = 0.641)
vegan::adonis(maxDepth ~ year, data = SumDepTemp, permutations = 999, 
              method = "euclidean")
#Differences across months - Significant (p = 0.001)
vegan::adonis(maxDepth ~ month, data = SumDepTemp, permutations = 999, 
              method = "euclidean")
#Differences between day and night across individuals - Interaction non-sig (p = 0.721)
vegan::adonis(maxDepth ~ diel*ptt, data = SumDepTemp, permutations = 999, 
              method = "euclidean")
#Differences between day and night across months - Interaction non-sig (p = 0.798)
vegan::adonis(maxDepth ~ diel*month, data = SumDepTemp, permutations = 999, 
              method = "euclidean")
#Differences among individuals across months - Interaction sig (p = 0.001)
vegan::adonis(maxDepth ~ ptt*month, data = SumDepTemp, permutations = 999, 
              method = "euclidean")
#Differences between day and night across individuals and months
vegan::adonis(maxDepth ~ diel*ptt*month, data = SumDepTemp, permutations = 999, 
              method = "euclidean")

#Mean temperatures
#Differences between day and night - Non significant (p = 0.421)
vegan::adonis(meanTemp ~ diel, data = SumDepTemp %>% drop_na(meanTemp), 
              permutations = 999, method = "euclidean")
#Differences among animals - Significant (p = 0.001)
vegan::adonis(meanTemp ~ ptt, data = SumDepTemp %>% drop_na(meanTemp),  
              method = "euclidean")
#Differences across years - Non significant (p = 0.092)
vegan::adonis(meanTemp ~ year, data = SumDepTemp %>% drop_na(meanTemp), 
              method = "euclidean")
#Differences across months - Non significant (p = 0.738)
vegan::adonis(meanTemp ~ month, data = SumDepTemp %>% drop_na(meanTemp), 
              method = "euclidean")

## Diel diving patterns ----------------------------------------------------
#Maximum daily depth reached per animal during the morning and the evening
DielDives <- {Movs %>%
    group_by(ptt, date, diel) %>% summarise(maxD = max(depth, na.rm = T))%>% 
    ggplot(aes(x = date, y = maxD, colour = diel))+
    geom_line()+
    scale_color_manual(values = c("#7899BC", "black"))+
    scale_y_reverse()+
    facet_wrap(ptt~., scales = "free_x")+
    theme_bw()+
    labs(x = "", y = "Maximum depth\n (m)")+
    theme(panel.grid = element_blank(),
          panel.grid.major.y = element_line(colour = "#edeeeb"),
          legend.title = element_blank(),
          legend.text = element_text(family = "sans", size = 12),
          legend.position = "top",
          axis.title = element_text(family = "sans", size = 12),
          axis.text = element_text(family = "sans", size = 12))}
#Saving graph
ggsave("Scripts_Vert/Figures/MaxDepthDiel.tiff", DielDives, device = "tiff", dpi = 400, width = 35, 
       height = 32.85, units = "cm")

# Depth-time distribution -------------------------------------------------
x <- Movs %>% 
  #Create a year category
  mutate(year = lubridate::year(date),
         #Pooling 2016 and 2017 together because data belongs to the same sharks
         year = case_when(year == 2016 | year == 2017 ~ "2016-2017",
                          T ~ as.character(year))) %>% 
  #Calculate the total time each animal was tracked during the day and night per year
  group_by(ptt, diel, year) %>% mutate(timeT = sum(timeDiff, na.rm = T)) %>% 
  #Calculate the proportion of time each animal was tracked at each depth bin (day and night)
  group_by(ptt, diel, year, DepthBin) %>% 
  summarise(TotTime = mean(timeT, na.rm = T),
            BinTime = sum(timeDiff, na.rm = T),
            PropTime = (BinTime/TotTime)*100) %>%
  #Calculate the mean and SE of the proportion of time spent at each depth bin
  group_by(diel, year, DepthBin) %>% 
  summarise(MeanProp = mean(PropTime, na.rm = T),
            SE_Prop = plotrix::std.error(PropTime, na.rm = T))

#Create graph with the depth bins reordered in reverse so the shallower depths appear on the top
DepthTime <- x %>% ggplot(aes(reorder(DepthBin, desc(DepthBin))))+ 
  #Facet per year (as columns)
  facet_grid(.~year)+
  #Subsetting the data based on diel patterns
  geom_bar(data = subset(x, diel == "Night"), aes(y = MeanProp, fill = diel), stat = "identity",
           position = "dodge", color = "black")+
  geom_errorbar(data = subset(x, diel == "Night"), aes(ymin = MeanProp, ymax = MeanProp+SE_Prop), 
                width = 0.2, position = position_dodge(.9))+
  #Second subset with negative y axis (MeanProp) so it appears on the left hand side of the plot
  geom_bar(data = subset(x, diel == "Morning"), aes(y = -MeanProp, fill = diel), stat = "identity",
           position = "dodge", color = "black")+
  geom_errorbar(data = subset(x, diel == "Morning"), aes(ymin = -MeanProp, ymax = -MeanProp-SE_Prop), 
                width = 0.2, position = position_dodge(.9))+
  #Flipping coordinates so depth appears on the y axis
  coord_flip()+
  #Changing the labels along the y axis (MeanProp), so they all appear as positive numbers
  scale_y_continuous(limits = c(-75, 75),
                     breaks = seq(-75, 75, 25), 
                     labels = c(seq(75, 25, -25), seq(0, 75, 25)))+
  #Change the fill of boxes - colourblind safe
  scale_fill_manual(values = c("#f5f5f5", "#5ab4ac"))+
  labs(x = "Depth classes (m)", y = "Time at depth (%)")+
  theme_bw()+
  theme(axis.title = element_text(family = "sans", size = 12),
        axis.text = element_text(family = "sans", size = 12),
        axis.ticks.length.y = unit(.15, "cm"),
        legend.position = "none", 
        strip.text = element_text(family = "sans", size = 12), 
        panel.grid = element_blank())

#Saving figure
ggsave("Scripts_Vert/Figures/DepthTime.tiff", DepthTime, device = "tiff", dpi = 400, width = 12,
       height = 5.5)
rm(x)

#Proportion of time spent at depth per month
x <- Movs %>% 
  mutate(year = lubridate::year(date),
         month = lubridate::month(date),
         monthName = factor(month.name[month], levels = month.name, ordered = T)) %>% 
  #Calculate the total time each animal was tracked during the day and night (per month)
  group_by(ptt, diel, monthName) %>% mutate(timeT = sum(timeDiff, na.rm = T)) %>% 
  #Calculate the proportion of time each animal was tracked at each depth bin (day and night)
  group_by(ptt, diel, monthName, DepthBin) %>% 
  summarise(TotTime = mean(timeT, na.rm = T),
            BinTime = sum(timeDiff, na.rm = T),
            PropTime = (BinTime/TotTime)*100) %>%
  #Calculate the mean and SE of the proportion of time spent at each depth bin
  group_by(diel, monthName, DepthBin) %>% 
  summarise(MeanProp = mean(PropTime, na.rm = T),
            SE_Prop = plotrix::std.error(PropTime, na.rm = T))

#Create graph with the depth bins reordered in reverse so the shallower depths appear on the top
DepthTimeMonth <- x %>% ggplot(aes(reorder(DepthBin, desc(DepthBin))))+ 
  facet_wrap(~monthName)+
  #Subsetting the data based on diel patterns
  geom_bar(data = subset(x, diel == "Night"), aes(y = MeanProp, fill = diel), stat = "identity",
           position = "dodge", color = "black")+
  geom_errorbar(data = subset(x, diel == "Night"), aes(ymin = MeanProp, ymax = MeanProp+SE_Prop), 
                width = 0.2, position = position_dodge(.9))+
  #Second subset with negative y axis (MeanProp) so it appears on the left hand side of the plot
  geom_bar(data = subset(x, diel == "Morning"), aes(y = -MeanProp, fill = diel), stat = "identity",
           position = "dodge", color = "black")+
  geom_errorbar(data = subset(x, diel == "Morning"), aes(ymin = -MeanProp, ymax = -MeanProp-SE_Prop), 
                width = 0.2, position = position_dodge(.9))+
  #Flipping coordinates so depth appears on the y axis
  coord_flip()+
  #Changing the labels along the y axis (PropTime), so they all appear as positive numbers
  scale_y_continuous(limits = c(-75, 75),
                     breaks = seq(-75, 75, 25), 
                     labels = c(seq(75, 25, -25), seq(0, 75, 25)))+
  #Change box colour fill - colourblind safe
  scale_fill_manual(values = c("#f5f5f5", "#5ab4ac"))+
  labs(x = "Depth classes (m)", y = "Time at depth (%)")+
  theme_bw()+
  theme(axis.title = element_text(family = "sans", size = 12),
        axis.text = element_text(family = "sans", size = 12),
        axis.ticks.length.y = unit(0.15, "cm"),
        legend.position = "none",
        strip.text = element_text(family = "sans", size = 12),
        panel.grid = element_blank())

#Saving figure
ggsave("Scripts_Vert/Figures/DepthTimeMonth.tiff", DepthTimeMonth, device = "tiff", dpi = 400,
       width = 12, height = 7.5)
rm(x)

# Summary Statistics Table ------------------------------------------------
#Summary statistics of temperature amd depth per individual per month
x1 <-  {Movs %>% group_by(ptt, lubridate::month(date)) %>% 
  #Calculate summary statistics for temperature and depth
  summarise(TempMean = round(mean(temperature, na.rm = T), 2),
            TempMin = min(temperature, na.rm = T),
            TempMax = max(temperature, na.rm = T),
            DepthMean = round(mean(depth, na.rm = T), 2),
            DepthMin = min(depth, na.rm = T),
            DepthMax = max(depth, na.rm = T)) %>% 
  #Creating new columns for temperature and depth merging all statistics under one column
  mutate(Temp = paste0(TempMean, " (", TempMin, " - ", TempMax, ")"),
         Depth = paste0(DepthMean, " (", DepthMin, " - ", DepthMax, ")")) %>% 
  #Renaming column for easier data manipulation
  rename("Month"=`lubridate::month(date)`) %>% 
  #Changing numbers representing months to abbreviation for months
  mutate(Month = month.abb[Month]) %>% 
  #Select only columns containing summary of calculated statistics
  select(ptt, Month, Temp, Depth) %>% 
  pivot_longer(cols = Temp:Depth, names_to = "Measurements", 
               values_to = "val") %>% 
  pivot_wider(names_from = ptt, values_from = val)
#Saving summary statistics
write.csv(x1, "Scripts_Vert/Outputs/summary.csv", row.names = F)
#Removing unused variable
rm(x1)}

# Temperature-depth profiles -----------------------------------------------
## Daily temperature-depth profiles for tags with temperature data
#Only including tags deployed in 2016
#Dividing plots - 2 weeks of data - 157561 (15 days) and 157562 (14 days)
DT_15days <- {Movs %>% filter(ptt == 157561 | ptt == 157562) %>% 
  ggplot(aes(x = dateTime, y = depth))+
  #The line colour will change based on the temperature
  geom_line(aes(color = temperature))+
  facet_wrap(ptt~., nrow = 2)+
  #Reversing the y axis so water surface (0 m) is at the top of the graph
  #Note that limits given need to be reversed, deepest point first and then the surface
  #Limits based on the maximum recorded depth for all tags
  scale_y_reverse(limits = c(plyr::round_any(max(Movs$depth, na.rm = T), 100, 
                                                  f = ceiling), 0),
                  #Breaks do not need to be specified in reverse
                  breaks = (seq(0, plyr::round_any(max(Movs$depth, na.rm = T), 100, 
                                                   f = ceiling), by = 250)))+
  #Changing color palette so blue represents cooler temperatures and red warmer temperatures
  scale_color_distiller(type = "div", palette = "Spectral", 
                        #NA values in the Temperature column to be shown as gray 
                        na.value = "#7f7f7f",
                        #Ensure colour bar covers the entire temperature range in the database
                        limits = c(floor(min(Movs$temperature, na.rm = T)-1),
                                   ceiling(max(Movs$temperature, na.rm = T)+1)))+
  #Changing label titles
  labs(x = "", y = "Depth (m)\n")+
  #Changing plot theme to black and white
  theme_bw()+
  #Moving legend to the top
  theme(legend.position = "top", 
        #Removing grid lines along the x axis
        panel.grid = element_blank(),
        #Changing the colour of grid lines along the y axis
        panel.grid.major.y = element_line(colour = "#edeeeb"),
        #Changing the font family and size for both axes, the legend and the facet labels
        axis.text = element_text(family = "sans", size = 11),
        axis.title = element_text(family = "sans", size = 12),
        legend.text = element_text(family = "sans", size = 12),
        legend.title = element_text(family = "sans", size = 12),
        strip.text = element_text(family = "sans", size = 11),
        #Changing legend margin
        legend.margin = margin(0, 0, 0, 0),
        #Decreasing the top and lower margins of the legend box 
        legend.box.margin = margin(-1, 0, 0, 0),
        plot.margin = unit(c(5.5, 15, 2, 2), "points"),
        #Decreasing distance between y axis labels and y axis title 
        axis.title.y = element_text(vjust = -1.5),
        panel.spacing = unit(7, "points"))+
  #Changing the position of the colourbar title to the top
  guides(color = guide_colorbar(title.position = "top", 
                                #Changing the title of the colourbar
                                title = "Temperature (°C)",
                                #Changing the colour of colorbar ticks and frame to black
                                ticks.colour = "black",
                                frame.colour = "black",
                                #Decreasing rge width of frame thickness
                                frame.linewidth = 0.1))}

#Dividing plots - 52 days of data - 157560 and 157566
DT_52days <- {Movs %>% filter(ptt == 157560 | ptt == 157566) %>% 
  ggplot(aes(x = dateTime, y = depth))+
  #The line colour will change based on the temperature
  geom_line(aes(color = temperature))+
  facet_wrap(ptt~., nrow = 2)+
  #Reversing the y axis so water surface is that the top of the graph
  scale_y_continuous(trans = "reverse", 
                     #Adding labels to y axis every 200 m until the maximum depth rounded up to
                     #the nearest 100
                     breaks = (seq(0, plyr::round_any(max(Movs$depth, na.rm = T), 100, 
                                                      f = ceiling), by = 250)),
                     #Adding limits to depth based on the maximum recorded depth by all tags
                     limits = c(plyr::round_any(max(Movs$depth, na.rm = T), 100, 
                                               f = ceiling), 0))+
  #Changing color palette - blue represents cooler temperatures and red warmer temperatures
  scale_color_distiller(type = "div", palette = "Spectral", na.value = "#7f7f7f", 
                        #Changing limits to min and max temperature values recorded by all tags
                        limits = c(floor(min(Movs$temperature, na.rm = T)-1),
                                   ceiling(max(Movs$temperature, na.rm = T)+1)))+
  #Removing all axis labels prior to creating a  composite figure
  labs(x = "", y = "")+
  #Changing plot theme to black and white
  theme_bw()+
  #Moving legend to the top
  theme(legend.position = "top", 
        #Removing grid lines
        panel.grid = element_blank(),
        panel.grid.major.y = element_line(colour = "#edeeeb"),
        #Changing the font family and size for both axes, the legend and the facet labels
        axis.text.x = element_text(family = "sans", size = 11),
        axis.title = element_text(family = "sans", size = 12),
        legend.text = element_text(family = "sans", size = 12),
        legend.title = element_text(family = "sans", size = 12),
        strip.text = element_text(family = "sans", size = 11),
        #Removing y axis ticks labels 
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        #Changing legend margin
        legend.margin = margin(0, 0, 0, 0),
        #Decreasing the top and lower margins of the legend box 
        legend.box.margin = margin(-1, 0, 0, 0),
        plot.margin = unit(c(5.5, 15, 2, 2), "points"),
        #Decreasing distance between y axis labels and y axis title 
        axis.title.y = element_text(vjust = -1.5),
        panel.spacing = unit(7, "points"))+
  #Changing the position of the colourbar title to the top
  guides(color = guide_colorbar(title.position = "top", 
                                #Changing the title of the colourbar
                                title = "Temperature (°C)"))}

#Dividing plots - 140 days of data - 157567
DT_140days <- {Movs %>% filter(ptt == 157567) %>% 
  ggplot(aes(x = dateTime, y = depth))+
  #The line colour will change based on the temperature
  geom_line(aes(color = temperature))+
  facet_wrap(~ptt)+
  #Reversing the y axis so water surface is that the top of the graph
  scale_y_continuous(trans = "reverse", 
                     #Adding labels to y axis every 200 m until the maximum depth rounded up to
                     #the nearest 100
                     breaks = (seq(0, plyr::round_any(max(Movs$depth, na.rm = T), 100, 
                                                      f = ceiling), by = 250)),
                     limits = c(plyr::round_any(max(Movs$depth, na.rm = T), 100, 
                                               f = ceiling) ,0))+
  #Changing color palette - blue represents cooler temperatures and red warmer temperatures
  scale_color_distiller(type = "div", palette = "Spectral", na.value = "#7f7f7f", 
                        limits = c(floor(min(Movs$temperature, na.rm = T)-1), 
                                   ceiling(max(Movs$temperature, na.rm = T)+1)))+
  #Removing axis labels prior to creating a figure composite
  labs(x = "", y = "")+
  #Changing plot theme to black and white
  theme_bw()+
  #Moving legend to the top
  theme(legend.position = "top", 
        #Removing grid lines
        panel.grid = element_blank(),
        panel.grid.major.y = element_line(colour = "#edeeeb"),
        #Changing the font family and size for both axes, the legend and the facet labels
        axis.text = element_text(family = "sans", size = 11),
        axis.title = element_text(family = "sans", size = 12),
        legend.text = element_text(family = "sans", size = 12),
        legend.title = element_text(family = "sans", size = 12),
        strip.text = element_text(family = "sans", size = 11),
        #Removing y axis ticks and labels
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        #Changing legend margin
        legend.margin = margin(0, 0, 0, 0),
        #Decreasing the top and lower margins of the legend box 
        legend.box.margin = margin(-1, 0, 0, 0),
        plot.margin = unit(c(5.5, 15, 2, 2), "points"),
        #Decreasing distance between y axis labels and y axis title 
        axis.title.y = element_text(vjust = -1.5))+
  #Changing the position of the colourbar title to the top
  guides(color = guide_colorbar(title.position = "top", 
                                #Changing the title of the colourbar
                                title = "Temperature (°C)"))}

#Dividing plots - 31 days of data - 157568
DT_31days <- {Movs %>% filter(ptt == 157568) %>% 
  ggplot(aes(x = dateTime, y = depth))+
  #The line colour will change based on the temperature
  geom_line(aes(color = temperature))+
  facet_wrap(~ptt)+
  #Reversing the y axis so water surface is that the top of the graph
  scale_y_continuous(trans = "reverse", 
                     #Adding labels to y axis every 200 m until the maximum depth rounded up to
                     #the nearest 100
                     breaks = (seq(0, plyr::round_any(max(Movs$depth, na.rm = T), 100, 
                                                      f = ceiling), by = 250)),
                     limits = c(plyr::round_any(max(Movs$depth, na.rm = T), 100,
                                                f = ceiling), 0))+
  #Changing color palette - blue represents cooler temperatures and red warmer temperatures
  scale_color_distiller(type = "div", palette = "Spectral", na.value = "#7f7f7f",
                        limits = c(floor(min(Movs$temperature, na.rm = T)-1), 
                                   ceiling(max(Movs$temperature, na.rm = T)+1)))+
  #Changing label titles
  labs(x = "", y = "Depth (m)\n")+
  #Changing plot theme to black and white
  theme_bw()+
  #Moving legend to the top
  theme(legend.position = "top", 
        #Removing grid lines
        panel.grid = element_blank(),
        panel.grid.major.y = element_line(colour = "#edeeeb"),
        #Changing the font family and size for both axes, the legend and the facet labels
        axis.text = element_text(family = "sans", size = 11),
        axis.title = element_text(family = "sans", size = 12),
        legend.text = element_text(family = "sans", size = 12),
        legend.title = element_text(family = "sans", size = 12),
        strip.text = element_text(family = "sans", size = 11),
        #Changing legend margin
        legend.margin = margin(0, 0, 0, 0),
        #Decreasing the top and lower margins of the legend box 
        legend.box.margin = margin(-1, 0, 0, 0),
        plot.margin = unit(c(5.5, 15, 2, 2), "points"),
        #Decreasing distance between y axis labels and y axis title 
        axis.title.y = element_text(vjust = -1.5))+
  #Changing the position of the colourbar title to the top
  guides(color = guide_colorbar(title.position = "top", 
                                #Changing the title of the colourbar
                                title = "Temperature (°C)"))}

#Creating composite figure using the four plots prepared above 
DT_2016 <- {ggarrange(ggarrange(DT_15days, DT_52days, ncol = 2, common.legend = T, 
                          #Decreasing the width of the second plot because it does not have labels
                          widths = c(1, 0.85),
                          #Add labels for each pair of graphs
                          labels = c("A", "B"), font.label = list(family = "sans")),
                #Create second composite with single graphs
                ggarrange(DT_31days, DT_140days, ncol = 2, legend = "none", widths = c(1, 0.85),
                    labels = c("C", "D"), font.label = list(family = "sans")),
                #Organising composite graphs on top of each other
                #Second composite is half the height than the first
                nrow = 2, heights = c(1, 0.5))}
#Saving composite figure to disk
ggsave("Scripts_Vert/Figures/DepthTempProfiles2016.tiff", DT_2016, "tiff", dpi = 400, units = "cm", 
       width = 20, height = 18.28) #width = 17.5, height = 15.99 also works, but no less than this

#Removing graphs no longer in use (ending in days)
rm(list=ls(pattern = "*days$"))

#Tag deployed in 2019 (only one - 178974)
DT_2019 <- {Movs %>% filter(ptt == 178974) %>% 
  ggplot(aes(x  = dateTime, y = depth, color = temperature))+
  geom_line()+
  facet_grid(~ptt)+
  scale_y_reverse(breaks = (seq(0, plyr::round_any(max(Movs$depth, na.rm = T), 100, 
                                                   f = ceiling), by = 250)),
                  limits = c(plyr::round_any(max(Movs$depth, na.rm = T), 100,
                                             f = ceiling), 0))+
  scale_color_distiller(type = "div", palette = "Spectral", na.value = "#7f7f7f",
                        limits = c(floor(min(Movs$temperature, na.rm = T)-1), 
                                   ceiling(max(Movs$temperature, na.rm = T)+1)))+
  #Changing label titles
  labs(x = "", y = "Depth (m)\n")+
  #Changing plot theme to black and white
  theme_bw()+
  #Moving legend to the top
  theme(legend.position = "top", 
        #Removing grid lines
        panel.grid = element_blank(),
        panel.grid.major.y = element_line(colour = "#edeeeb"),
        #Changing the font family and size for both axes, the legend and the facet labels
        axis.text = element_text(family = "sans", size = 11),
        axis.title = element_text(family = "sans", size = 12),
        legend.text = element_text(family = "sans", size = 12),
        legend.title = element_text(family = "sans", size = 12),
        strip.text = element_text(family = "sans", size = 11),
        #Changing legend margin
        legend.margin = margin(0, 0, 0, 0),
        #Decreasing the top and lower margins of the legend box 
        legend.box.margin = margin(-1, 0, 0, 0),
        plot.margin = unit(c(5.5, 15, 2, 2), "points"),
        #Decreasing distance between y axis labels and y axis title 
        axis.title.y = element_text(vjust = -1.5))+
  #Changing the position of the colourbar title to the top
  guides(color = guide_colorbar(title.position = "top", 
                                #Changing the title of the colourbar
                                title = "Temperature (°C)",
                                ticks.colour = "black", 
                                frame.colour = "black", frame.linewidth = 0.1))}
#Saving figure
ggsave("Scripts_Vert/Figures/DepthTempProfiles2019.tiff", DT_2019, device = "tiff", dpi = 400)


## Daily diving profiles for tags with no temperature data recorded
#All sharks tagged in 2018 lasting less than 90 days 
DP_90days <- {Movs %>% filter(ptt %in% c(174048, 174049)) %>% 
  ggplot(aes(x = dateTime, y = depth))+
  #The line colour will change based on the temperature
  geom_line(color = "#7f7f7f")+
  facet_wrap(~ptt, nrow = 2)+
  #Reversing the y axis so water surface is that the top of the graph
  #Adding labels to y axis every 250 m until the maximum depth rounded up to the nearest 100
  scale_y_reverse(breaks = (seq(0, plyr::round_any(max(Movs$depth, na.rm = T), 100, f = ceiling), 
                                   by = 250)),
                  limits = c(plyr::round_any(max(Movs$depth, na.rm = T), 100, f = ceiling), 0))+
  #Changing label titles
  labs(x = "", y = "Depth (m)\n")+
  #Changing plot theme to black and white
  theme_bw()+
  #Personalising theme
  #Removing grid lines
  theme(panel.grid = element_blank(),
        panel.grid.major.y = element_line(colour = "#edeeeb"),
        #Changing the font family and size for both axes, the legend and the facet labels
        axis.text = element_text(family = "sans", size = 11),
        axis.title = element_text(family = "sans", size = 12),
        strip.text = element_text(family = "sans", size = 11),
        plot.margin = unit(c(5.5, 5.5, 0.5, 2), "points"),
        #Decreasing distance between y axis labels and y axis title 
        axis.title.y = element_text(vjust = -1.5))}

#Shark tags deployed in 2018 with over 90 days of data
DP_180days <- {Movs %>% filter(ptt %in% c(174050, 174051)) %>% 
    ggplot(aes(x = dateTime, y = depth))+
    #The line colour will change based on the temperature
    geom_line(color = "#7f7f7f")+
    facet_wrap(~ptt, nrow = 2)+
    #Reversing the y axis so water surface is that the top of the graph
    #Adding labels to y axis every 250 m until the maximum depth rounded up to the nearest 100
    scale_y_reverse(breaks = (seq(0, plyr::round_any(max(Movs$depth, na.rm = T), 100, f = ceiling), 
                                  by = 250)),
                    limits = c(plyr::round_any(max(Movs$depth, na.rm = T), 100, f = ceiling), 0))+
    #Changing label titles
    labs(x = "", y = "")+
    #Changing plot theme to black and white
    theme_bw()+
    #Personalising theme
    #Removing grid lines
    theme(panel.grid = element_blank(),
          panel.grid.major.y = element_line(colour = "#edeeeb"),
          #Changing the font family and size for both axes, the legend and the facet labels
          axis.text = element_text(family = "sans", size = 11),
          axis.title = element_text(family = "sans", size = 12),
          strip.text = element_text(family = "sans", size = 11),
          plot.margin = unit(c(5.5, 5.5, 0.5, 2), "points"),
          #Decreasing distance between y axis labels and y axis title 
          axis.title.y = element_text(vjust = -1.5),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank())}

#Creating composite figure
DP_2018 <- ggarrange(DP_90days, DP_180days, ncol = 2, widths = c(1, 0.85), labels = c("A", "B"))
#Saving composite figure
ggsave("Scripts_Vert/Figures/DivingProfiles_2018.tiff", DP_2018, "tiff", dpi = 400, units = "cm",
       width = 20, height = 18.77)

rm(list = c(ls(pattern = glob2rx("DT_*days")), ls(pattern = glob2rx("DP_*days"))))

# Thermal envelope analysis -----------------------------------------------
#Create new dataframe with only the columns that will be used in the analysis
TE <- Movs %>% as_tibble() %>% 
  #Selecting only the columns to be used in the analysis
  select(ptt, dateTime, date, depth, temperature) %>% 
  #Reclassifying temperature data in new depth bins to minimise data loss prior to clustering
  mutate(DepthBin = case_when(depth < 25 ~ "<25",
                              depth >= 25 & depth <= 50 ~ "25-50",
                              depth > 50 & depth <= 100 ~ "50-100",
                              depth > 100 & depth <= 200 ~ "100-200",
                              depth > 200 & depth <= 300 ~ "200-300",
                              depth > 300 ~ ">300")) %>% 
  #Changing newly created columns into ordered factors
  mutate(DepthBin = factor(DepthBin,
                           levels = c("<25", "25-50", "50-100", "100-200", "200-300", ">300"), 
                           ordered = T)) %>% 
  #Exclude any rows with no temperature data
  drop_na(temperature) %>% 
  #Creating new ID columns
  #ID shot version includes an unique number for each date for which a tag has information
  mutate(ID = group_indices(., ptt, date), 
         #ID long version includes the tag number and date
         ID2 = paste(ptt, date, sep = "_")) %>% 
  #Grouping variables by ID prior to calculating min and max temperatures for the entire day
  group_by(ID) %>% 
  mutate(maxT = max(temperature),
         minT = min(temperature)) %>% 
  #Grouping variables by ID and depth bin to calculate min and max temperature per depth bin
  group_by(ID2, DepthBin) %>% 
  summarise(ID = mean(ID),
            minTBin = min(temperature),
            maxTBin = max(temperature),
            minT = mean(minT),
            maxT = mean(maxT)) %>% 
  #Remove >300 bin to minimise data loss
  filter(DepthBin != ">300") %>% 
  #Manipulating dataset so it is wider
  pivot_wider(names_from = DepthBin, values_from = c(minTBin, maxTBin))

#Creating copy of above dataset for further manipulating and clustering calculations
clusterTE <- TE %>% 
  #Remove any NA values and select only variables to be used in clustering
  drop_na() %>% ungroup() %>% select(-ID2) %>% 
  #Move ID column as rownames
  column_to_rownames(var = "ID") %>% 
  #Scaling all variables to remove bias towards variables with larger values
  scale()

#Determining optimal number of clusters
fviz_nbclust(clusterTE, kmeans, method = "gap_stat")
#Calculating k-means clusters based on optimal number of clusters identified above
KM_TE <- kmeans(clusterTE, centers = 2, nstart = 25)
#Visualising clusters
KM_fig <- fviz_cluster(KM_TE, data = clusterTE, palette = "jco", 
                       ggtheme = theme_bw())+
  scale_y_continuous(expand = expansion(add = 0.5))+
  theme(axis.title = element_text(family = "sans", size = 12),
        axis.text = element_text(family = "sans", size = 12),
        title = element_text(family = "sans", size = 12),
        legend.text = element_text(family = "sans", size = 12),
        legend.title  = element_text(family = "sans", size = 12))+
  labs(shape = "Cluster", fill = "Cluster", color = "Cluster")
#Saving the plot
ggsave("Scripts_Vert/Figures/KM_fig.tiff", KM_fig, "tiff", dpi = 400)

#Extracting temperature and depth data from combined data using results from k-means
KM_Res <- data.frame(cluster = KM_TE$cluster) %>% rownames_to_column(var = "ID") %>% 
  as_tibble() %>% 
  mutate_all(as.numeric) %>% 
  #Extracting ID2 column from summary TE dataset summary dataset to obtain dates and ptt info
  left_join(TE %>% select(ID, ID2), by = "ID") %>% 
  #Separate dates and ptt to match rows with information from the Movs dataset containing all info
  separate(ID2, into = c("ptt", "date"), sep = "_", remove = T) %>%
  select(-ID) %>% 
  #Changing to date format prior to match rows
  mutate(date = as.POSIXct(date, tz = "UTC")) %>%
  left_join(Movs %>% select(ptt, date, dateTime, depth, temperature, lon, lat), 
            by = c("ptt", "date")) %>%
  #Dropping any rows containing NAs
  drop_na() %>% 
  #Keeping only information for depths of up to 300 m
  filter(depth <= 300) %>% 
  #Reclassifying data to calculate means and SE per 5 m bins per cluster across all tags 
  mutate(x = cut(depth, c(seq(0, 300, 5), Inf), seq(0,300, 5), include.lowest = T)) %>% 
  group_by(x, cluster) %>% 
  summarise(meanT = mean(temperature),
            seT = plotrix::std.error(temperature)) %>% 
  ungroup() %>% 
  #Changing depth bins to numeric prior to graphing info
  mutate(x = as.numeric(as.character(x)))

#Creating graph with mean temperatures per cluster as identified by K-means
MeanTempDepth <- KM_Res %>%
  ggplot(aes(x = meanT))+
  geom_point(data = subset(KM_Res, cluster == 1), aes(y = x), color = "#EFC000FF", na.rm = T)+
  #Geom path creates a line that connects points based on the order in which they appear
  geom_path(data = subset(KM_Res, cluster == 1), aes(y = x), color = "#EFC000FF", na.rm = T)+
  #Adding error bars based on calculated SE values
  geom_errorbarh(data = subset(KM_Res, cluster == 1), 
                 aes(xmin = meanT-seT, xmax = meanT+seT, y = x), color = "#EFC000FF")+
  geom_point(data = subset(KM_Res, cluster == 2), aes(y = x), color = "#0073C2FF")+
  geom_path(data = subset(KM_Res, cluster == 2), aes(y = x), color = "#0073C2FF")+
  geom_errorbarh(data = subset(KM_Res, cluster == 2), 
                 aes(xmin = meanT-seT, xmax = meanT+seT, y = x), color = "#0073C2FF")+
  scale_y_reverse()+
  theme_bw()+
  labs(x = "Temperature (°C)", y = "Depth (m)")+
  theme(axis.text = element_text(family = "sans", size = 12),
        axis.title = element_text(family = "sans", size = 12))
#Saving graph
ggsave("Scripts_Vert/Figures/MeanTempDepth.tiff", MeanTempDepth, "tiff", dpi = 400)


#Hierarchical clusters
#Calculating gap statistics to estimate number of clusters
BinHclust <- cluster::clusGap(clusterTE, FUNcluster = hcut, K.max = 10)
fviz_gap_stat(BinHclust) #5 clusters identified as ideal number
BinHclust #5 clusters identified as ideal number
rm(BinHclust)
#Calculating hierarchical clusters - based on Euclidean disimilarity matrix
res.hc <- hclust(dist(clusterTE), method = "ward.D2")
#Visualising hierarchical clusters using the ideal number of clusters previously identified
Hclus <- fviz_dend(res.hc, cex = 0.75, k = 5, palette = "jco")
Hclus <- Hclus + theme(axis.title.y = element_text(family = "sans", size = 12),
              axis.text.y = element_text(family = "sans", size = 12),
              title = element_text(family = "sans", size = 12))
#Saving hierarchical cluster
ggsave("Scripts_Vert/Figures/HierarchicalCluster.tiff", Hclus, device = "tiff", dpi = 400, 
      width = 46.16, height = 20.36, units = "cm")
#Extract the IDs that make up every cluster
HClus_res <- data.frame(cluster = cutree(res.hc, 5)) %>% 
  rownames_to_column("ID") %>% 
  as_tibble() %>% 
  mutate_all(as.numeric) %>% 
  #Extracting ID2 column from summary TE dataset summary dataset to obtain dates and ptt info
  left_join(TE %>% select(ID, ID2), by = "ID") %>% 
  #Separate dates and ptt to match rows with information from the Movs dataset containing all info
  separate(ID2, into = c("ptt", "date"), sep = "_", remove = T) %>%
  select(-ID) %>% 
  #Changing to date format prior to match rows
  mutate(date = as.POSIXct(date, tz = "UTC")) %>%
  left_join(Movs %>% select(ptt, date, dateTime, depth, temperature, lon, lat), 
            by = c("ptt", "date")) %>%
  #Dropping any rows containing NAs
  drop_na() %>% 
  #Keeping only information for depths of up to 300 m
  filter(depth <= 300) %>% 
  #Reclassifying data to calculate means and SE per 5 m bins per cluster across all tags 
  mutate(x = cut(depth, c(seq(0, 300, 5), Inf), seq(0,300, 5), include.lowest = T)) %>% 
  group_by(x, cluster) %>% 
  summarise(meanT = mean(temperature),
            seT = plotrix::std.error(temperature)) %>% 
  ungroup() %>% 
  #Changing depth bins to numeric prior to graphing info
  mutate(x = as.numeric(as.character(x)),
         cluster = factor(cluster))
  
#Creating graph with mean temperatures per cluster as identified by hierarchical cluster
MeanTempDepth_HC <- HClus_res %>%
  ggplot(aes(x = meanT, y = x, color = cluster))+
  geom_point()+geom_path()+geom_errorbarh(aes(xmin = meanT-seT, xmax = meanT+seT))+
  scale_color_jco()+
  scale_y_reverse()+
  theme_bw()+
  labs(x = "Temperature (°C)", y = "Depth (m)")+
  theme(axis.text = element_text(family = "sans", size = 12),
        axis.title = element_text(family = "sans", size = 12),
        legend.text = element_text(family = "sans", size = 12),
        legend.title = element_text(family = "sans", size = 12),
        legend.background = element_blank(),
        legend.position = c(0.85, 0.15),
        legend.direction = "horizontal")+
  guides(color = guide_legend(title = "Cluster", title.position = "top", ncol = 2))
#Saving graph
ggsave("Scripts_Vert/Figures/MeanTempDepth_HClus.tiff", MeanTempDepth_HC, "tiff", dpi = 400)


#Determining how many tags and how many unique days were used in the cluster analysis
ClusterTags <- data.frame(cluster = cutree(res.hc, 5)) %>% 
  rownames_to_column("ID") %>% 
  as_tibble() %>% 
  mutate_all(as.numeric) %>% 
  #Extracting ID2 column from summary TE dataset summary dataset to obtain dates and ptt info
  left_join(TE %>% select(ID, ID2), by = "ID") %>%
  separate(ID2, c("ptt", "date"), sep = "_") 
#Unique tags used
ClusterTags %>% distinct(ptt)
#Number of unique dates used
ClusterTags %>% distinct(date) %>% count()


# Recorded depths comparisons by tag --------------------------------------
#Boxplot of distribution of depths recorded by each tag
RecDepths <- Movs %>% ggplot(aes(x = ptt, y = depth, colour = ptt))+
  #Outliers presented as grey dots
  geom_boxplot(outlier.colour = "grey60")+
  #Flipping coordinates
  coord_flip()+theme_bw()+
  labs(y = "Depth (m)")+
  theme(legend.position = "none",
        axis.title.y = element_blank(),
        axis.text = element_text(family = "sans", size = 12),
        axis.title.x = element_text(family = "sans", size = 12))
  
#Saving plot
ggsave("Scripts_Vert/Figures/RecDepthsTag.tiff", RecDepths, device = "tiff", dpi = 400)




