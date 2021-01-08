# Analysis of Hammerhead sharks (Sphyrna lewini) vertical movements
# Step 2: Summarising movements
# Date of latest update: 2021-01-08
# Version: 1
# Function: Summarising vertical movement data of S. lewini.
# Prepared for the Sharks Ecology Project of the Charles Darwin Foundation (CDF)
# Script related to publication entitled "xxxx" in xxxx magazine.

# Uploading relevant libraries --------------------------------------------
library(tidyverse)
library(rstatix)

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

# Summarising movement data -----------------------------------------------
Movs %>% group_by(ptt, diel) %>% 
  get_summary_stats(depth, temperature, type = "mean_se")

#Obtaining summary for depth and temperature per tag
Movs %>% select(ptt, diel, date, depth, temperature) %>% 
  group_by(ptt, diel, date) %>% 
  mutate(minDepth = min(depth), 
         maxDepth = max(depth)) %>% 
  group_by(ptt, diel) %>% 
  get_summary_stats(minDepth, maxDepth, temperature, type = "mean_se") %>% 
  select(-n) %>% 
  pivot_wider(names_from = variable, values_from = c(mean, se)) %>% 
  unite("MaxDepth", mean_maxDepth, se_maxDepth, sep = " ± ") %>% 
  unite("MinDepth", mean_minDepth, se_minDepth, sep = " ± ") %>% 
  unite("Temp", mean_temperature, se_temperature, sep = " ± ") %>% 
  pivot_wider(names_from = diel, values_from = c(MaxDepth, MinDepth, Temp)) %>% 
  write.csv("Scripts_Vert/Outputs/minMaxDepth.csv", row.names = F)

Movs %>% select(ptt, diel, date, depth) %>% 
  group_by(ptt, diel) %>% 
  get_summary_stats(depth, type = "mean_se") %>% 
  ggplot(aes(x = diel, y = mean, col = ptt))+geom_point()

select(-n) %>% 
  pivot_wider(names_from = variable, values_from = c(mean, se)) %>% 
  unite("MeanDepth", mean_depth, se_depth, sep = " ± ") %>% 
  pivot_wider(names_from = diel, values_from = MeanDepth)

# Summary Statistics Table ------------------------------------------------
#Summary statistics of temperature and depth per individual per month
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
