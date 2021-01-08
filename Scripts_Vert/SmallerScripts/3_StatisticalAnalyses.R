# Analysis of Hammerhead sharks (Sphyrna lewini) vertical movements
# Step 3: Statistical Analyses
# Date of latest update: 2021-01-08
# Version: 1
# Function: Applying statistical analyses to vertical movement data
# Prepared for the Sharks Ecology Project of the Charles Darwin Foundation (CDF)
# Script related to publication entitled "xxxx" in xxxx magazine.

library(tidyverse)
library(vegan)

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
adonis(minDepth ~ diel, data = SumDepTemp, permutations = 999, 
              method = "euclidean")
#Differences among animals - Non significant (p = 0.342)
adonis(minDepth ~ ptt, data = SumDepTemp, permutations = 999, 
              method = "euclidean")
#Differences across years - Non significant (p = 0.209)
adonis(minDepth ~ year, data = SumDepTemp, permutations = 999, 
              method = "euclidean")
#Differences across months - Non significant (p = 0.974)
adonis(minDepth ~ month, data = SumDepTemp, permutations = 999, 
              method = "euclidean")

#Maximum depth
#Differences between day and night - Significant (p = 0.001)
adonis(maxDepth ~ diel, data = SumDepTemp, permutations = 999, 
              method = "euclidean")
#Differences among animals - Significant (p = 0.001)
adonis(maxDepth ~ ptt, data = SumDepTemp, permutations = 999, 
              method = "euclidean")
#Differences across years - Non significant (p = 0.641)
adonis(maxDepth ~ year, data = SumDepTemp, permutations = 999, 
              method = "euclidean")
#Differences across months - Significant (p = 0.001)
adonis(maxDepth ~ month, data = SumDepTemp, permutations = 999, 
              method = "euclidean")
#Differences between day and night across individuals - Interaction non-sig (p = 0.721)
adonis(maxDepth ~ diel*ptt, data = SumDepTemp, permutations = 999, 
              method = "euclidean")
#Differences between day and night across months - Interaction non-sig (p = 0.798)
adonis(maxDepth ~ diel*month, data = SumDepTemp, permutations = 999, 
              method = "euclidean")
#Differences among individuals across months - Interaction sig (p = 0.001)
adonis(maxDepth ~ ptt*month, data = SumDepTemp, permutations = 999, 
              method = "euclidean")
#Differences between day and night across individuals and months
adonis(maxDepth ~ diel*ptt*month, data = SumDepTemp, permutations = 999, 
              method = "euclidean")

#Mean temperatures
#Differences between day and night - Non significant (p = 0.421)
adonis(meanTemp ~ diel, data = SumDepTemp %>% drop_na(meanTemp), 
              permutations = 999, method = "euclidean")
#Differences among animals - Significant (p = 0.001)
adonis(meanTemp ~ ptt, data = SumDepTemp %>% drop_na(meanTemp),  
              method = "euclidean")
#Differences across years - Non significant (p = 0.092)
adonis(meanTemp ~ year, data = SumDepTemp %>% drop_na(meanTemp), 
              method = "euclidean")
#Differences across months - Non significant (p = 0.738)
adonis(meanTemp ~ month, data = SumDepTemp %>% drop_na(meanTemp), 
              method = "euclidean")