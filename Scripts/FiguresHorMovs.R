#Creating graphs for HHs horizontal movements publication
#Author: Denisse Fierro Arcos
#Date: 2020-08-20

# Libraries ---------------------------------------------------------------
library(tidyverse)
library(patchwork)


# Setting up data ---------------------------------------------------------
#Creating data frame with supporting information for tags of interest
NameShark <- data.frame(Tags = c(157567, 174049, 178974, 198362, 198367),
                        #Shorten names for tags
                        ID = c("PAT5", "PAT8", "PAT11", "SPOT1", "SPOT3"),
                        #Start dates for graphs
                        SDate = as.POSIXct(c("2017-04-27", "2018-04-21", "2019-04-02", "2020-02-22",
                                  "2020-02-12")),
                        #End dates for graphs
                        EDate = as.POSIXct(c("2017-05-03", "2018-04-26", "2019-04-08", "2020-03-02", 
                                  "2020-02-18")),
                        #Likely date animal was fished/killed
                        Fished = as.POSIXct(c("2017-04-30", "2018-04-23", "2019-04-04", "2020-02-25", 
                                   "2020-02-15")))

#Accessing -Histo files for tags of interest
HistoFiles <- list.files("Data", pattern = paste(paste0(NameShark$Tags, "-Histos.csv"), 
                                         collapse = "|"), 
                         recursive = T, full.names = T)

# Uploading data ----------------------------------------------------------
for(i in 1:length(HistoFiles)){
  #Opeining file
  x <- read.csv(HistoFiles[i])
  x <- x %>% 
    #Selecting only rows with TAT (Temperature at Time) data
    filter(grepl("TAT", HistType)) %>%
    #Keeping only the number of bins where data has been recorded (based on NumBins information)
    select(Ptt, Date, paste0("Bin",1:max(x$NumBins, na.rm = T))) %>% 
    #Changing Date column to date format
    mutate(Date = as.Date(lubridate::parse_date_time(Date, orders = "HMS dbY"))) %>%
    #Keeping only rows between the start and end dates in NameShark data frame
    #Also keeping NA rows as this contains information about temperature bin
    filter(Date >= NameShark$SDate[i] & Date <= NameShark$EDate[i] | is.na(Date)) %>% 
    #Add the shorten name to the data frame
    mutate(ID = NameShark$ID[i])
  #Changing column names starting with 'Bin' so it contains the temperature associated w/that bin
  colnames(x)[grepl("Bin", colnames(x))] <- paste0("Bin", x[1, grepl("Bin", colnames(x))])
  #Remove first row, which only contains temperature for each bin
  x <- x[-1,] %>% 
    #Change dataset so there is a row for each temperature bin/date
    pivot_longer(cols = starts_with("Bin"), names_to = "TempBin", names_prefix = "Bin",
                 values_to = "TimeProp") %>% 
    #Change new TempBin column to ordered factor
    mutate(TempBin = factor(TempBin, levels = unique(.$TempBin), ordered = T)) %>% 
    #Calculating the mean and standard error for each date and temperature bin
    group_by(ID, Date, TempBin) %>% 
    summarise(Ptt = mean(Ptt), MeanTime = mean(TimeProp), 
              SETime = plotrix::std.error(TimeProp)) %>% 
    #Create new column with tags to be used to color graphs (prior to being fished, fished, and after)
    mutate(Fish = case_when(Date < NameShark$Fished[i] ~ "Before",
                            Date == NameShark$Fished[i] ~ "Fished",
                            Date > NameShark$Fished[i] ~ "After")) %>% 
    mutate(Fish = factor(Fish, levels = c("Before", "Fished", "After"), ordered = T))
  #Assign result to shorten tag name
  assign(paste(NameShark$ID[i]), x)
  rm(x)
}

# Creating plots ----------------------------------------------------------
#Creating bar graphs for each date in the data frame and colours it based on its tag (Before, Fished, After)
plotGraph2 <- function(x){
  x %>% 
    #Mean proportion of time in a day on y axis and temperature on x axis
    ggplot(aes(x = TempBin, y = MeanTime, fill = Fish))+
    #Bars showing actual number included in the temperature bin with a black outline
    geom_bar(stat = "identity", color = "black", position = position_dodge())+
    #Error bar showing standard error using SETime information
    geom_errorbar(aes(ymin = MeanTime, ymax = MeanTime+SETime), width = 0.2,
                  position = position_dodge(.9))+
    #Separate graphs in a panel for each date
    facet_wrap(Date~.)+
    #Manually setting the colours to be used for each label under Fish (Before, Fished, After)
    scale_fill_manual(values = c("#0abde3", "#e30a36", "#f4a35f"))+
    #Including a main title w/shorten tag name, and changing x and y axes titles
    labs(title = x$ID, y = "Time at temperature (%)", x = "Temperature (°C) \n")+
    #Applying minimal theme - no background
    theme_minimal()+
    #Removing legend
    theme(legend.position = "none",
          #Changing the size and font of axes text and title
          axis.text = element_text(family = "sans", size = 11),
          axis.title = element_text(family = "sans", size = 12),
          #Changing the size and font of main graph title
          title = element_text(family = "sans", size = 12), 
          #Removing grids
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          #Changing the size and font of strip text
          strip.text = element_text(family = "sans", size = 12))
  }

#Applying plot function to each tag
P1 <- plotGraph2(SPOT1)+
  #Removing axes titles
  theme(axis.title = element_blank())
P2 <- plotGraph2(SPOT3)+
  #Removing axes titles
  theme(axis.title = element_blank())
P3 <- plotGraph2(PAT5)+
  #Removing x axis title, so we only keep title for y axis
  theme(axis.title.x = element_blank())
P4 <- plotGraph2(PAT8)+
  #Removing axes titles
  theme(axis.title = element_blank())
P5 <- plotGraph2(PAT11)+
  #Removing y axis titles, so we only keep title for x axis
  theme(axis.title.y = element_blank())

#To save individual plots use the line below. Change plot name accordingly
# ggsave("../Outputs/Reports/PlotsHorMovs/SPOT1.tiff", P1, "tiff", dpi = 400, 
#        height = 16, width = 22, units = "cm")


#Creating mosaic - with graphs stacked on top of one another
Pall <- P1/P2/P3/P4/P5
#Saving mosaic
ggsave("../Outputs/Reports/PlotsHorMovs/AllTags.tiff", Pall, "tiff", dpi = 400, 
       height = 50, width = 25, units = "cm")


# Horizontal movement graph -----------------------------------------------
library(ggforce)
Movs %>% filter(ptt == 157567) %>%
  mutate(Col = case_when(date < "2017-04-25" ~ "Gray",
                         TRUE ~ "Red")) %>% 
  ggplot(aes(x = dateTime, y = temperature, color = Col))+
  geom_line()+geom_point()+
  #Zooming in on a particular date range
  facet_zoom(x = dateTime >= "2017-04-25 16:00:00" &
               dateTime <= "2017-04-30 23:50:00")

