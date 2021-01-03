
# Libraries ---------------------------------------------------------------
library(tidyverse)

#Creating data frame with supporting information for tags of interest
NameShark <- data.frame(Tags = c(157567, 174049, 178974, 198362, 198367),
                        ID = c("PAT5", "PAT8", "PAT11", "SPOT1", "SPOT3"),
                        SDate = c("2017-04-27", "2018-04-21", "2019-04-02", "2020-02-22",
                                  "2020-02-12"),
                        EDate = c("2017-05-03", "2018-04-26", "2019-04-08", "2020-03-02", 
                                  "2020-02-18"),
                        Fished = c("2017-04-30", "2018-04-23", "2019-04-04", "2020-02-25", 
                                   "2020-02-15"))


#Accessing -Histo data for the tags of interest
HistoFiles <- list.files("../../HH_Movements/Data/", 
                         pattern = paste(paste0(NameShark$Tags, "-Histos.csv"), 
                                         collapse = "|"), 
                         recursive = T, full.names = T)

# Uploading data ----------------------------------------------------------
for(i in 1:length(HistoFiles)){
  x <- read.csv(HistoFiles[i])
  x <- x %>% 
    filter(grepl("TAT", HistType)) %>%
    select(Ptt, Date, paste0("Bin",1:max(x$NumBins, na.rm = T))) %>% 
    mutate(Date = readr::parse_date(Date, format = "%H:%M:%S %d-%b-%Y")) %>% 
    filter(Date >= NameShark$SDate[i] & Date <= NameShark$EDate[i] | is.na(Date)) %>% 
    mutate(ID = NameShark$ID[i])
  colnames(x)[grepl("Bin", colnames(x))] <- paste0("Bin", x[1, grepl("Bin", colnames(x))])
  x <- x[-1,] %>% 
    pivot_longer(cols = starts_with("Bin"), names_to = "TempBin", names_prefix = "Bin",
                 values_to = "TimeProp") %>% 
    mutate(TempBin = factor(TempBin, levels = unique(.$TempBin), ordered = T)) %>% 
    group_by(ID, Date, TempBin) %>% 
    summarise(Ptt = mean(Ptt), MeanTime = mean(TimeProp), 
              SETime = plotrix::std.error(TimeProp)) %>% 
    mutate(Fish = case_when(Date < NameShark$Fished[i] ~ "Before",
                            Date == NameShark$Fished[i] ~ "Fished",
                            Date > NameShark$Fished[i] ~ "After")) %>% 
    mutate(Fish = factor(Fish, levels = c("Before", "Fished", "After"), ordered = T))
  assign(NameShark$ID[i], x)
  rm(x)
}


# Cleaning Data -----------------------------------------------------------
{SPOT1 %>% ggplot(aes(x = TempBin, y = TimeProp))+
  geom_histogram(stat = "identity")+
  facet_wrap(~Date)+
  labs(title = SPOT1$ID)

SPOT1 <- SPOT1 %>% mutate(Fished = case_when(Date == "2020-02-24" ~ "Before",
                                    TRUE ~ "After")) %>% 
  group_by(ID, TempBin, Fished) %>% 
  summarise(MeanTime = mean(TimeProp),
            SETime = plotrix::std.error(TimeProp))


SPOT3 %>% ggplot(aes(x = TempBin, y = TimeProp))+
  geom_histogram(stat = "identity")+
  facet_wrap(~Date)+
  labs(title = SPOT3$ID)

SPOT3 <- SPOT3 %>% mutate(Fished = case_when(Date <= "2020-02-14" ~ "Before",
                                             TRUE ~ "After")) %>% 
  group_by(ID, TempBin, Fished) %>% 
  summarise(MeanTime = mean(TimeProp),
            SETime = plotrix::std.error(TimeProp))


PAT5 %>% ggplot(aes(x = TempBin, y = TimeProp))+
  geom_histogram(stat = "identity")+
  facet_wrap(~Date)+
  labs(title = PAT5$ID)

PAT5 <- PAT5 %>% mutate(Fished = case_when(Date <= "2017-04-30" ~ "Before",
                                             TRUE ~ "After")) %>% 
  group_by(ID, TempBin, Fished) %>% 
  summarise(MeanTime = mean(TimeProp),
            SETime = plotrix::std.error(TimeProp))


PAT8 %>% ggplot(aes(x = TempBin, y = TimeProp))+
  geom_histogram(stat = "identity")+
  facet_wrap(~Date)+
  labs(title = PAT8$ID)

PAT8 <- PAT8 %>% mutate(Fished = case_when(Date <= "2018-04-22" ~ "Before",
                                           TRUE ~ "After")) %>% 
  group_by(ID, TempBin, Fished) %>% 
  summarise(MeanTime = mean(TimeProp),
            SETime = plotrix::std.error(TimeProp))


PAT11 %>% ggplot(aes(x = TempBin, y = TimeProp))+
  geom_histogram(stat = "identity")+
  facet_wrap(~Date)+
  labs(title = PAT11$ID)

PAT11 <- PAT11 %>% mutate(Fished = case_when(Date <= "2019-04-04" ~ "Before",
                                           TRUE ~ "After")) %>% 
  group_by(ID, TempBin, Fished) %>% 
  summarise(MeanTime = mean(TimeProp),
            SETime = plotrix::std.error(TimeProp))}

# Creating plots ----------------------------------------------------------
plotGraph <- function(x){
  x %>% 
    ggplot(aes(reorder(TempBin, desc(TempBin))))+
    geom_bar(data = subset(x, Fished == "Before"), 
                      aes(y = -MeanTime, fill = Fished), stat = "identity",
                      position = "dodge", color = "black")+
    geom_errorbar(data = subset(x, Fished == "Before"), 
                  aes(ymin = -MeanTime, ymax = -MeanTime-SETime), 
                  width = 0.2, position = position_dodge(.9))+
    geom_bar(data = subset(x, Fished == "After"), 
                                  aes(y = MeanTime, fill = Fished), stat = "identity",
                                  position = "dodge", color = "black")+
    geom_errorbar(data = subset(x, Fished == "After"), 
                  aes(ymin = MeanTime, ymax = MeanTime+SETime), 
                  width = 0.2, position = position_dodge(.9))+
    coord_flip()+
    scale_y_continuous(limits = c(-100, 100),
                       breaks = seq(-100, 100, 25), 
                       labels = c(seq(100, 25, -25), seq(0, 100, 25)))+
    labs(title = x$ID, y = "Time at temperature (%)", x = "Temperature (°C) \n")+
    theme_minimal()+
    theme(legend.position = "none",
          axis.text = element_text(family = "sans", size = 12),
          axis.title = element_text(family = "sans", size = 12),
          title = element_text(family = "sans", size = 12), 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())

p1 <- plotGraph(SPOT1)
p2 <- plotGraph(SPOT3)+theme(axis.title.y = element_blank(), 
                             axis.text.y = element_blank())
p3 <- plotGraph(PAT5)
p4 <- plotGraph(PAT8)+theme(axis.title.y = element_blank(), 
                            axis.text.y = element_blank())
p5 <- plotGraph(PAT11)+theme(axis.title.y = element_blank(), 
                             axis.text.y = element_blank())

library(patchwork)
p_all <- (p1+p2)/(p3+p4+p5)
ggsave("../../Outputs/Reports/TATs.tiff", p_all, "tiff", dpi = 400, units = "cm", 
       width = 19.65, height = 18)}


plotGraph2 <- function(x){
  x %>% 
    ggplot(aes(x = TempBin, y = MeanTime, fill = Fish))+
    geom_bar(stat = "identity", color = "black")+
    geom_errorbar(aes(ymin = MeanTime, ymax = MeanTime+SETime), width = 0.2,
                  position = position_dodge(.9))+
    facet_wrap(Date~.)+
    scale_fill_manual(values = c("#0abde3", "#e30a36", "#f4a35f"))+
    labs(title = x$ID, y = "Time at temperature (%)", x = "Temperature (°C) \n")+
    theme_minimal()+
    theme(legend.position = "none",
          axis.text = element_text(family = "sans", size = 11),
          axis.title = element_text(family = "sans", size = 12),
          title = element_text(family = "sans", size = 12), 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.text = element_text(family = "sans", size = 12))
  }
  
P1 <- plotGraph2(SPOT1)+
  theme(axis.title = element_blank())
P2 <- plotGraph2(SPOT3)+
  theme(axis.title = element_blank())
P3 <- plotGraph2(PAT5)+
  theme(axis.title.x = element_blank())
P4 <- plotGraph2(PAT8)+
  theme(axis.title = element_blank())
P5 <- plotGraph2(PAT11)+
  theme(axis.title.y = element_blank())

#Saving individual plots
ggsave("../../Outputs/Reports/PlotsHorMovs/SPOT1.tiff", P1, "tiff", dpi = 400, 
       height = 16, width = 22, units = "cm")
ggsave("../../Outputs/Reports/PlotsHorMovs/SPOT3.tiff", P2, "tiff", dpi = 400, 
       height = 16, width = 22, units = "cm")
ggsave("../../Outputs/Reports/PlotsHorMovs/PAT5.tiff", P3, "tiff", dpi = 400, 
       height = 16, width = 22, units = "cm")
ggsave("../../Outputs/Reports/PlotsHorMovs/PAT8.tiff", P4, "tiff", dpi = 400, 
       height = 16, width = 22, units = "cm")
ggsave("../../Outputs/Reports/PlotsHorMovs/PAT11.tiff", P5, "tiff", dpi = 400, 
       height = 16, width = 22, units = "cm")


#Creating mosaic
library(patchwork)
Pall <- P1/P2/P3/P4/P5
#Saving mosaic
ggsave("../../Outputs/Reports/PlotsHorMovs/AllTags.tiff", Pall, "tiff", dpi = 400, 
       height = 50, width = 25, units = "cm")

# Horizontal movement graph -----------------------------------------------
library(ggforce)
Movs %>% filter(ptt == 157567) %>%
  mutate(Col = case_when(date < "2017-04-25" ~ "Gray",
                         TRUE ~ "Red")) %>% 
  ggplot(aes(x = dateTime, y = temperature, color = Col))+
  geom_line()+geom_point()+
  facet_zoom(x = dateTime >= "2017-04-25 16:00:00" &
               dateTime <= "2017-04-30 23:50:00")


Movs %>% filter(ptt == 157567 & date >= "2017-04-23") %>%
  mutate(Col = case_when(date < "2017-04-25" ~ "Gray",
                         TRUE ~ "Red")) %>% 
  ggplot(aes(x = dateTime, y = temperature, color = Col))+
  geom_line()+geom_point()





Movs %>% filter(ptt == 157561 | ptt == 157562) %>% 
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
                                frame.linewidth = 0.1))