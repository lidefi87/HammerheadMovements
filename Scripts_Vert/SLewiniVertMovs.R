# Analysis of Hammerhead sharks (Sphyrna lewini) vertical movements -------
# Date of latest update: 2021-01-04
# Version: 1
# Prepared for the Sharks Ecology Project of the Charles Darwin Foundation (CDF)
# Script related to publication entitled "xxxx" in xxxx magazine.

# Uploading relevant libraries --------------------------------------------
library(tidyverse)
library(stringr)
library(suncalc)
library(ggpubr)
library(factoextra)
library(ggsci)

# Calling relevant scripts ------------------------------------------------
source("Scripts_Vert/1_Clean_Data.R")
source("Scripts_Vert/2_SummarisingMovements.R")
source("Scripts_Vert/3_StatisticalAnalyses.R")
source("Scripts_Vert/4_GeneratingFigures.R")



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




