# Analysis of Hammerhead sharks (Sphyrna lewini) vertical movements
# Step 4: Generating figures
# Date of latest update: 2021-01-08
# Version: 1
# Function: Making figures of vertical shark movements
# Prepared for the Sharks Ecology Project of the Charles Darwin Foundation (CDF)
# Script related to publication entitled "xxxx" in xxxx magazine.

# Uploading relevant libraries --------------------------------------------
library(tidyverse)
library(stringr)
library(ggpubr)
library(factoextra)
library(ggsci)

#Making composite figure with depth histograms for each tag
#Extracting legend from first figure in the list
leg <- get_legend(plot_list[[1]])
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

