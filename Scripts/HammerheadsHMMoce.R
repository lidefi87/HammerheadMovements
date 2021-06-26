############################################################################################
# Title: HMMoce applied to miniPAT data
# Author: Jeremy Vaudo
# Adapted by: Denisse Fierro Arcos
# Version: 1.3
# Date last updated: 2020-04-28
# Updates: Use of tidyverse, removal of parallel analysis, automatic identification of 
# values for different variables, change URL for SST data (ncdcOisst dataset) and added
# option to access data from THREDDS website. However, note that this code is not compatible
# with data from this site
############################################################################################

#Clear global environment
rm(list = ls())

################################ CALLING RELEVANT LIBRARIES ################################
library(tidyverse)
library(raster)
library(HMMoce)
library(curl)
library(ncdf4)

############################ MANUALLY SET VARIABLES FOR ANALYSIS ###########################
#Set variable to remove or keep the Atlantic Ocean
removeATL <- F
#Set spatial limits manually if needed, otherwise limits will be calculated from data in
#uploaded datasets
#Global limits
sp.lim.global <- list(lonmim = -95, lonmax = -75, latmin = -5, latmax = 11)
#Local limits
sp.lim <- list(lonmim = -95, lonmax = -75, latmin = -5, latmax = 11)

# MODELLING PARAMETERS
#Indicate the parameters to be run in the model:
# 1 = light, 2 = sst, 3 = OCH, 4 = WOA, 5 = hycom
#Can be combined with if() statements around calc functions: 
#if (any(likVec == 5) & !exists('L.5')){calc.hycom(...)}
likVec <- c(1, 2, 3, 5) #Light, SST, OCH and hycom
#Indicate which parameter combinations (from list above) you want the model to run
# 1: 1,2      #6: 3,5
# 2: 1,3      #7: 1,2,3
# 3: 1,5      #8: 1,2,5
# 4: 2,3      #9: 1,3,5
# 5: 2,5      #10: 2,3,5
run.idx <- c(1, 2, 3, 7, 8)
#Indicate if tracks produced by the model are to be saved. Note that sometimes saving the 
# track plots causes an error. It is recommended to be set to FALSE.
plottracks <- T
# vector of appropriate bounding in filter. see ?hmm.filter for more info
bndVec <- c(NA, 5, 10)
# vector of appropriate migr kernel speed. see ?makePar for more info.
#Speed < 2 have been added
parVec <- c(1.5, 1.75, 2, 2.25, 2.5, 2.75, 3)

#############################################################################################
########################## DO NOT CHANGE CODE UNDERNEATH THIS LINE ##########################
#############################################################################################

########################## DEFINING FUNCTIONS TO BE USED IN MODEL ###########################
# GET ENVIRONMENTAL DATA
get.oi.sst <- function(limits, time, filename = '', download.file = TRUE, dir = getwd()){
  options(warn = -1)
  dir.create(file.path(dir), recursive = TRUE, showWarnings = FALSE)
  setwd(dir)
  ## Set the base URL based on the start date. If the ending date exceeds the
  ## period for this experiment, then print a warning and truncate the output
  ## early.
  expts <- data.frame(
    start = c(as.Date('1981-09-01')),
    end = c(Sys.Date() + 1),
    url = c('http://coastwatch.pfeg.noaa.gov/erddap/griddap/ncdcOisst2Agg.nc?sst'))
    #THREDDS website access - DFA
    # url = c("https://www.ncei.noaa.gov/thredds/ncss/OisstBase/NetCDF/V2.0/AVHRR/"))
  
  #Activate if using THREDDS website - DFA
  #Change longitude from -180 +180 to 360
  # limits[[1]] <- limits[[1]]+360
  # limits[[2]] <- limits[[2]]+360
  
  if(time[1] < expts$start[1])
    stop('Data begins at %s and is not available at %s.',
         strftime(expts$start[1], '%d %b %Y'),
         strftime(time[1], '%d %b %Y'))
  if(time[1] > expts$end[nrow(expts)])
    stop('Data ends at %s and is not available at %s.',
         strftime(expts$end[nrow(expts)], '%d %b %Y'),
         strftime(time[1], '%d %b %Y'))
  for(i in seq(nrow(expts))) {
    if((time[1] >= expts$start[i]) & (time[1] <= expts$end[i]))
      url = expts$url[i]}
  
  ## Add the time domain.
  if(length(time) == 2){
    url <- sprintf('%s[(%s:00:00Z):1:(%s:00:00Z)]',
                  url, strftime(time[1], '%Y-%m-%dT00'),
                  strftime(time[2], '%Y-%m-%dT00'))
  } else if(length(time) == 1){
    url <- sprintf('%s[(%s:00:00Z):1:(%s:00:00Z)]',
                  url, strftime(time[1], '%Y-%m-%dT00'),
                  strftime(time[1], '%Y-%m-%dT00'))}
    #Activate if using THREDDS website - DFA
    # url <- sprintf('%s%s/avhrr-only-v2.%s.nc?', 
                   # url1, strftime(time[1], '%Y%m'), strftime(time[1], '%Y%m%d'))}
  
  ## Add the spatial domain.
  url <- sprintf('%s[(%s):1:(%s)][(%s):1:(%s)]',
                url, limits[[3]], limits[[4]], limits[[1]], limits[[2]])
  #Activate if using THREDDS website - DFA
  # url <- sprintf("%svar=sst&north=%s&west=%s&east=%s&south=%s",
  #                url, limits[[4]], limits[[1]], limits[[2]], limits[[3]])
  #Add second part of THREDDS url - DFA
  # url <- sprintf("%s&disableProjSubset=on&horizStride=1&time_start=%s&time_end=%s&timeStride=1&vertCoord=",
  #         url, strftime(time[1], '%Y-%m-%dT00%%3A00%%3A00Z'), 
  #         strftime(time[1], '%Y-%m-%dT00%%3A00%%3A00Z'))
  print(url)
  
  ## Download the data if a filename was provided.
  if(filename != ''){
    if(download.file == TRUE){
      #download.file(url, filename, method = 'auto')
      curl_download(url, filename, quiet=FALSE)
    } else if(download.file == FALSE){
      system(sprintf('curl -o "%s" "%s"', filename, url))}}
  
  options(warn = 0)
  return(url)}

get.env <- function(uniqueDates = NULL, filename = NULL, type = NULL, spatLim = NULL, 
                    resol = NULL, save.dir = getwd(), sst.type = NULL, depLevels = NULL){
  if(is.null(type)){
    stop('Type of environmental data desired not specified.')
  } else if(type == 'sst'){
    if(is.null(sst.type)){
      warning('Warning: if type=sst then sst.type should be specified. Using default GHRsst.')
      sst.type = 'ghr'}
    
    if(sst.type == 'oi'){
      for(i in 1:length(uniqueDates)){
        time <- as.Date(uniqueDates[i])
        repeat{
          get.oi.sst(spatLim, time, filename = paste(filename, '_', time, '.nc', sep = ''), 
                     download.file = TRUE, dir = save.dir) 
          # filenames based on dates from above
          tryCatch({
            err <- try(RNetCDF::open.nc(paste(save.dir, filename, '_', time, '.nc', 
                                              sep = '')), silent = T)}, 
            error = function(e){print(paste('ERROR: Download of data at ', time, 
                                            ' failed. Trying call to server again.', 
                                            sep = ''))})
          if(class(err) != 'try-error') break}}
    } else if(sst.type == 'ghr'){
      for(i in 1:length(uniqueDates)){
        time <- as.Date(uniqueDates[i])
        repeat{
          HMMoce::get.ghr.sst(spatLim, time, 
                              filename = paste(filename, '_', time, '.nc', sep = ''), 
                              download.file = TRUE, dir = save.dir) 
          # filenames based on dates from above
          tryCatch({
            err <- try(RNetCDF::open.nc(paste(save.dir, filename, '_', time, '.nc', 
                                              sep = '')), silent = T)}, 
            error = function(e){print(paste('ERROR: Download of data at ', time, 
                                            ' failed. Trying call to server again.', 
                                            sep = ''))})
          if(class(err) != 'try-error') break}}}
  } else if(type == 'hycom'){
    for(i in 1:length(uniqueDates)){
      time <- as.Date(uniqueDates[i])
      repeat{
        get.hycom2(spatLim, time, filename = paste(filename, '_', time, '.nc', sep = ''),
                   download.file = TRUE, dir = save.dir, depLevels = depLevels) 
        tryCatch({
          err <- try(RNetCDF::open.nc(paste(save.dir,'/', filename, '_', time, '.nc', 
                                            sep = '')), silent = T)}, 
          error = function(e){print(paste('ERROR: Download of data at ', time, 
                                          ' failed. Trying call to server again.', 
                                          sep = ''))})
        if(class(err) != 'try-error') break}}
  } else if (type == 'woa'){
    if(is.null(resol)){
      stop('Error: If type = woa then resol must be specified. See ?get.env for help.')}
    
    filename <- get.woa(save.dir = save.dir, resol = resol)
    print(paste('WOA data downloaded to ', filename,'...', sep=''))}}


get.hycom2 <- function(limits, time, vars = c('water_temp'), include_latlon = TRUE,
                       filename = '', download.file = TRUE, dir = getwd(), depLevels = NULL){
  dir.create(file.path(dir), recursive = TRUE, showWarnings = FALSE)
  setwd(dir)
  ## Set the base URL based on the start date. If the ending date exceeds the
  ## period for this experiment, then print a warning and truncate the output
  ## early.
  expts = data.frame(
    start <- c(as.Date('1992-10-02'), as.Date('1995-08-01'),
            as.Date('2013-01-01'), as.Date('2013-08-21'),
            as.Date('2014-04-05'), as.Date('2016-04-18'),
            as.Date('2018-11-21')),
    end <- c(as.Date('1995-07-31'), as.Date('2012-12-31'),
          as.Date('2013-08-20'), as.Date('2014-04-04'),
          as.Date('2016-04-17'), as.Date('2018-11-20'),
          Sys.Date() + 1),
    url <- c('http://ncss.hycom.org/thredds/ncss/GLBu0.08/expt_19.0/',
          'http://ncss.hycom.org/thredds/ncss/GLBu0.08/expt_19.1/',
          'http://ncss.hycom.org/thredds/ncss/GLBu0.08/expt_90.9?',
          'http://ncss.hycom.org/thredds/ncss/GLBu0.08/expt_91.0?',
          'http://ncss.hycom.org/thredds/ncss/GLBu0.08/expt_91.1?',
          'http://ncss.hycom.org/thredds/ncss/GLBu0.08/expt_91.2?',
          'http://ncss.hycom.org/thredds/ncss/GLBv0.08/expt_93.0?'))
  
  if(time[1] < expts$start[1])
    stop('Data begins at %s and is not available at %s.',
         strftime(expts$start[1], '%d %b %Y'),
         strftime(time[1], '%d %b %Y'))
  if(time[1] > expts$end[nrow(expts)])
    stop('Data ends at %s and is not available at %s.',
         strftime(expts$end[nrow(expts)], '%d %b %Y'),
         strftime(time[1], '%d %b %Y'))
  for(i in seq(nrow(expts))) {
    if((time[1] >= expts$start[i]) & (time[1] <= expts$end[i]))
      url = expts$url[i]}
  
  if(any(grep('19', url))) url = sprintf('%s%s?', url, as.numeric(format(time, '%Y')))
  
  ## Add the variables.
  for(var in vars)
    url = sprintf('%svar=%s&', url, var)
  ## Add the spatial domain.
  url = sprintf('%snorth=%f&west=%f&east=%f&south=%f&horizStride=1&',
                url, limits[[4]], limits[[1]], limits[[2]], limits[[3]])
  # north, west, east, south
  
  ## Add the time domain.
  if(length(time) == 2){
    url = sprintf('%stime_start=%s%%3A00%%3A00Z&time_end=%s%%3A00%%3A00Z&timeStride=1&',
                  url, strftime(time[1], '%Y-%m-%dT00'),
                  strftime(time[2], '%Y-%m-%dT00'))
  } else if(length(time) == 1){
    url = sprintf('%stime_start=%s%%3A00%%3A00Z&time_end=%s%%3A00%%3A00Z&timeStride=1&',
                  url, strftime(time[1], '%Y-%m-%dT00'),
                  strftime(time[1], '%Y-%m-%dT00'))}
  
  ## Add the lat-lon points if requested.
  if(include_latlon)
    url = sprintf('%saddLatLon=true&', url)
  ## Finish the URL.
  if (is.null(depLevels)){
    url = sprintf('%sdisableProjSubset=on&vertCoord=&accept=netcdf', url)
  } else{
    url = paste(url,'disableProjSubset=on&vertCoord=', depLevels, '&accept=netcdf', sep='')}
  
  print(url)
  
  ## Download the data if a filename was provided.
  if(filename != ''){
    if(download.file == TRUE){
      #download.file(url, filename, method = 'auto')
      curl_download(url, filename, quiet=FALSE)
    } else if(download.file == FALSE){
      system(sprintf('curl -o "%s" "%s"', filename, url))}}
  
  return(url)}

#Modification to fix make.L
make.L <- function (L1, L2 = NULL, L3 = NULL, known.locs = NULL, L.mle.res, 
                    dateVec = NULL, locs.grid = NULL, iniloc = NULL, bathy = NULL, 
                    pdt = NULL){
  if (class(L1) == "list") {
    if (length(L1) == 3) {
      L3 <- L1[[3]]
      L2 <- L1[[2]]
      L1 <- L1[[1]]}
    else if (length(L1) == 2) {
      L2 <- L1[[2]]
      L1 <- L1[[1]]}}
  
  if (!is.null(bathy) & is.null(pdt)) 
    stop("Error: if bathy is not NULL then a maxDep vector must be supplied.")
  if (!is.null(known.locs)) {
    print("Input known locations are being used...")
    L.locs <- L1 * 0
    known.locs$date <- as.Date(known.locs$date)
    if (is.null(dateVec)) {
      stop("Error: dateVec is null.")}
    if (is.null(locs.grid)) {
      stop("Error: locs.grid is null.")}
    
    lon <- locs.grid$lon[1, ]
    lat <- locs.grid$lat[, 1]
    kn.idx <- which(dateVec %in% known.locs$date)
    for (i in kn.idx) {
      known.locs.i <- known.locs[which(known.locs$date %in% 
                                         dateVec[i]), ]
      if (length(known.locs.i[, 1]) > 1) {
        known.locs.i <- known.locs.i[1, ]}
      x = which.min((known.locs.i$lon - lon)^2)
      y = which.min((known.locs.i$lat - rlat)^2)
      L.locs[[i]][raster::cellFromXY(L.locs[[i]], 
                                     known.locs.i[, c(3, 2)])] <- 1}}
  if (is.null(L2) & is.null(L3)) {
    print("One likelihood raster has been specified...")
    L <- L1}
  else if (!is.null(L2) & is.null(L3)) {
    print("Two likelihood rasters have been specified...")
    L <- L1 * 0
    L1[is.na(L1)] <- 0
    L2[is.na(L2)] <- 0
    L1[is.infinite(L1)] <- 0
    L2[is.infinite(L2)] <- 0
    naL1idx = raster::cellStats(L1, sum, na.rm = T) != 0
    naL2idx = raster::cellStats(L2, sum, na.rm = T) != 0
    naLidx = naL1idx + naL2idx
    idx1 = which(naLidx == 1)
    idx2 = which(naLidx == 2)
    for (ii in idx1) {
      L[[ii]] = L1[[ii]] + L2[[ii]]}
    for (ii in idx2) {
      L[[ii]] = L1[[ii]] * L2[[ii]]
      if (raster::cellStats(L[[ii]], sum, na.rm = T) == 
          0) {
        L[[ii]] = L2[[ii]]}}}
  else if (!is.null(L2) & !is.null(L3)) {
    print("Three likelihood rasters have been specified...")
    L <- L1 * 0
    L1[is.na(L1)] <- 0
    L2[is.na(L2)] <- 0
    L3[is.na(L3)] <- 0
    L1[is.infinite(L1)] <- 0
    L2[is.infinite(L2)] <- 0
    L3[is.infinite(L3)] <- 0
    naL1idx = raster::cellStats(L1, sum, na.rm = T) != 0
    naL2idx = raster::cellStats(L2, sum, na.rm = T) != 0
    naL3idx = raster::cellStats(L3, sum, na.rm = T) != 0
    naLidx = naL1idx + naL2idx + naL3idx
    idx1 = which(naLidx == 1)
    idx2 = which(naLidx == 2)
    idx3 = which(naLidx == 3)
    for (ii in idx1) {
      L[[ii]] = L1[[ii]] + L2[[ii]] + L3[[ii]]}
    for (ii in idx3) {
      L[[ii]] = L1[[ii]] * L2[[ii]] * L3[[ii]]
      if (raster::cellStats(L[[ii]], sum, na.rm = T) == 
          0 & naL1idx[ii] & naL2idx[ii]) {
        L[[ii]] <- L1[[ii]] * L2[[ii]]}
      if (raster::cellStats(L[[ii]], sum, na.rm = T) == 
          0 & naL1idx[ii] & naL3idx[ii]) {
        L[[ii]] <- L1[[ii]] * L3[[ii]]}
      if (raster::cellStats(L[[ii]], sum, na.rm = T) == 
          0 & naL2idx[ii] & naL3idx[ii]) {
        L[[ii]] <- L2[[ii]] * L3[[ii]]}
      if (raster::cellStats(L[[ii]], sum, na.rm = T) == 0) {
        L[[ii]] = L1[[ii]] + L2[[ii]] + L3[[ii]]}}
    for (ii in idx2) {
      if (raster::cellStats(L[[ii]], sum, na.rm = T) == 
          0 & naL1idx[ii] & naL2idx[ii]) {
        L[[ii]] <- L1[[ii]] * L2[[ii]]}
      if (raster::cellStats(L[[ii]], sum, na.rm = T) == 
          0 & naL1idx[ii] & naL3idx[ii]) {
        L[[ii]] <- L1[[ii]] * L3[[ii]]}
      if (raster::cellStats(L[[ii]], sum, na.rm = T) == 
          0 & naL2idx[ii] & naL3idx[ii]) {
        L[[ii]] <- L2[[ii]] * L3[[ii]]}
      if (raster::cellStats(L[[ii]], sum, na.rm = T) == 0) {
        L[[ii]] = L1[[ii]] + L2[[ii]] + L3[[ii]]}}}
  sumIdx <- which(raster::cellStats(L, sum, na.rm = T) != 0)
  for (i in sumIdx) {
    L[[i]] <- L[[i]]/(raster::cellStats(L[[i]], max, na.rm = T) * 
                        1.25)}
  
  if (!is.null(bathy)) {
    print("Starting bathymetry mask...")
    maxDep <- data.frame(dplyr::summarise(dplyr::group_by(pdt,Date), max(Depth)))
    maxDep$Date <- as.Date(maxDep$Date)
    maxDep.df <- data.frame(Date = dateVec)
    maxDep.df <- merge(maxDep.df, maxDep, by = "Date", all.x = T)
    maxDep.df[which(maxDep.df[, 2] <= 0), 2] <- 1
    maxDep.df[which(is.na(maxDep.df[, 2])), 2] <- 1
    bathy <- raster::resample(bathy, L)
    naLidx = which(raster::cellStats(L, sum, na.rm = T) != 0)
    for (i in naLidx) {
      b.i <- bathy
      b.i[b.i <= -maxDep.df[i, 2]] <- 1
      b.i[b.i != 1] <- NA
      L[[i]] <- L[[i]] * b.i}}
  
  if (!is.null(iniloc)) {
    if (!exists("L.locs")) {
      L.locs <- L1 * 0
      kn.idx <- c(1, length(dateVec))
    } else {
      kn.idx <- c(1, kn.idx, length(dateVec))}
    lon <- locs.grid$lon[1, ]
    lat <- locs.grid$lat[, 1]
    x = which.min((iniloc$lon[1] - lon)^2)
    y = which.min((iniloc$lat[1] - lat)^2)
    L.locs[[1]][raster::cellFromXY(L.locs[[1]], iniloc[1, 
                                                       c(5, 4)])] <- 1
    x = which.min((iniloc$lon[2] - lon)^2)
    y = which.min((iniloc$lat[2] - lat)^2)
    L.locs[[length(dateVec)]][raster::cellFromXY(L.locs[[length(dateVec)]], 
                                                 iniloc[2, c(5, 4)])] <- 1
    for (bb in kn.idx) {
      L[[bb]] <- L.locs[[bb]]}}
  
  print("Starting raster resample...")
  L.mle <- raster::resample(L, L.mle.res)
  L[L <= 1e-15] <- 1e-15
  L[is.na(L)] <- 1e-15
  L.mle[L.mle <= 1e-15] <- 1e-15
  L.mle[is.na(L.mle)] <- 1e-15
  #L[is.na(L)] <- 0
  #L.mle[is.na(L.mle)] <- 0
  try(assign(paste0(tagnumbers[i],'_idx',tt,'-Likelihoods'),L))
  try(save(list = paste0(tagnumbers[i],'_idx',tt,'-Likelihoods'), 
           file = paste0(tagnumbers[i],'_idx',tt,'-Likelihoods.rda')))
  L <- aperm(raster::as.array(raster::flip(L, direction = "y")), 
             c(3, 2, 1))
  L.mle <- aperm(raster::as.array(raster::flip(L.mle, direction = "y")), 
                 c(3, 2, 1))
  print("Finishing make.L...", sep = "")
  return(list(L = L, L.mle = L.mle))}

############ GETTING DATA FOR FISH TAGS, TAGGING AND POP UP DATES AND LOCATIONS #############  
#Getting directories and necessary datasets for all tagged animals
tagnumbers <- dir("Data/MiniPAT_Tags/", pattern = "[0-9]{6}", recursive = F)
#Uploading file containing tagging and pop up dates and locations
TagPopLocs <- read.csv("Data/MiniPAT_Tags/TagPopLocs.csv")

###################### SETTING SPATIAL LIMITS AND GETTING BATHYMETRY ########################
# SPATIAL LIMITS CALCULATION - IF NOT MANUALLY STATED ABOVE
if(!exists("sp.lim.global")){
  setwd("Data/MiniPAT_Tags/")
  locs.global <- NULL
  #Getting a list of file directories of datasets containing location data
  fishlocs <- list.files(pattern = "-Locations.csv", recursive = T)
  for(i in seq_along(fishlocs)){
    dat <- read.csv(fishlocs[i])
    locs.global <- plyr::rbind.fill(locs.global, dat)
    rm(dat)}
  sp.lim.global <- list(lonmin = min(locs.global$Longitude, na.rm = T), 
                        lonmax = max(locs.global$Longitude, na.rm = T),
                        latmin = min(locs.global$Latitude, na.rm = T), 
                        latmax = max(locs.global$Latitude, na.rm = T))
  setwd("../..")
  rm(locs.global, fishlocs, i)}

# BATHYMETRY
#If bathymetry is not downloaded, do it now, otherwise upload to environment
if(!dir.exists("Data/bathy/")){
  dir.create("Data/bathy", recursive = T)
  bathy <- try(get.bath.data(sp.lim.global$lonmin, sp.lim.global$lonmax, 
                             sp.lim.global$latmin, sp.lim.global$latmax, 
                             folder = "Data/bathy/"))
  }else{bathy <- raster::raster("Data/bathy/bathy.nc")} #request.nc may need to be changed to bathy.nc
rm(sp.lim.global)

# REMOVING ATLANTIC OCEAN FROM RASTER DATA - ONLY IF removeATL SET TO T
if(removeATL == T){
  bathy <- raster::raster(paste(bathy.dir, 'bathy.nc', sep='/'))
  bathy <- raster::reclassify(bathy,c(0,10000,4000))
  xmin <- bathy@extent@xmin*-1
  xmax <- bathy@extent@xmax*-1
  ymin <- bathy@extent@ymin
  ymax <- bathy@extent@ymax
  resolution <- (xmin-xmax)/bathy@ncols
  r1 <- ((ymax-9.1)/resolution)
  bathy[1:r1,] <- 4000
  r2 <- ((xmin-80)/resolution)
  r3 <- ((ymax-8.5)/resolution)
  bathy[r1:r3,1:r2] <- 4000
  r4 <- ((ymax-7.5)/resolution)
  r5 <- ((xmin-78)/resolution)
  bathy[r1:r4,r5:bathy@ncols] <- 4000
  rm("xmin","xmax","ymin","ymax","resolution","r1","r2","r3","r4","r5")}

############################ APPLYING HMMoce MODEL TO TAG DATA ############################
# SETTING UP TAG DATA
# Changing directory to where data is located
setwd("Data/MiniPAT_Tags/")
#Creating loop to read data from all tags
for(i in seq_along(tagnumbers)){
  #Set directory to tag folder
  setwd(tagnumbers[i])
  #Check if HMMoce folder exists, if not, create one
  if(!dir.exists("HMMoce")){
    dir.create("HMMoce")}
  #Extract dates when sharks were tagged and when tags popped off using the TagPopLocs dataset
  TagPop <- TagPopLocs %>% filter(Tag == tagnumbers[i]) %>% 
    dplyr::select(Status, Date) %>% mutate(Date = lubridate::parse_date_time(Date, "dmy"))
  #VECTOR OF DATES FROM DATA. THIS WILL BE THE TIME STEPS, T, IN THE LIKELIHOODS
  dateVec <- as.Date(seq(TagPop$Date[TagPop$Status == "Tag"], 
                 TagPop$Date[TagPop$Status == "Pop"], by = "days"))
  #Getting Sea Surface Temperature (SST) data
  tag.sst <- read.wc(ptt = tagnumbers[i], 
                     filename = list.files(pattern = "SST.csv$"), 
                     type = 'sst', 
                     tag = TagPop$Date[TagPop$Status == "Tag"], 
                     pop = TagPop$Date[TagPop$Status == "Pop"], verbose = T)
  tag.sst <- tag.sst$data %>% drop_na(Date)
  #Getting Depth-Temperature profile (PDT) data
  pdt <- read.wc(ptt = tagnumbers[i], 
                 filename = list.files(pattern = "PDTs.csv$"), 
                 type = 'pdt', 
                 tag = TagPop$Date[TagPop$Status == "Tag"], 
                 pop = TagPop$Date[TagPop$Status == "Pop"], verbose = T)
  pdt.udates <- pdt$udates[!is.na(pdt$udates)]
  pdt <- pdt$data %>% drop_na(Date)
  #Getting light-based positions from GPE2 (INSTEAD OF RAW LIGHTLOCS FROM PREVIOUS)
  locs <- read.csv(list.files(pattern = "Locations.csv$"), na.strings = "") %>% 
    drop_na(Date)
  locDates <- readr::parse_date(as.character(locs$Date), format = "%m/%d/%Y %H:%M")
  #SET SPATIAL LIMITS
  #If local limits have not been specified then calculate from data
  if(!exists("sp.lim")){
  locs.grid <- setup.locs.grid(locs)
  sp.lim <- list(lonmin = min(locs.grid$lon[1,]), lonmax = max(locs.grid$lon[1,]),
                 latmin = min(locs.grid$lat[,1]), latmax = max(locs.grid$lat[,1]))}
  rm(TagPop)
  # DOWNLOADING ENVIRONMENTAL DATA
  #Find if folder named HMMoce exist, if not create one. Set directory to this folder
  setwd("HMMoce/")
  if(!dir.exists("sst")){
    dir.create("sst", recursive = TRUE)
    try(get.env(uniqueDates = dateVec, type = 'sst', sst.type = 'oi', 
                filename = tagnumbers[i], spatLim = sp.lim, 
                save.dir = paste0(getwd(),"/sst/")))
    setwd("..")}
  #remove Atlantic from north of Panama
  if(removeATL == T){
    sstfiles <- list.files(path = "sst", recursive = T, full.names = T)
    for(x in sstfiles){
      nc <- nc_open(x, write = T)
      lon <- ncvar_get(nc, 'longitude')
      lat <- ncvar_get(nc, 'latitude')
      sst <- ncvar_get(nc, 'sst')
      xmin = lon[1]
      xmax = lon[length(lon)]
      ymin = lat[1]
      ymax = lat[length(lat)]
      resolution = (xmax-xmin)/length(lon)
      r1 = (abs(ymin-15)/resolution)
      sst[,r1:length(lat)] = NaN
      r2 = (abs(ymin-9.25)/resolution)
      r3 = (abs(xmin+84)/resolution)
      sst[r3:length(lon),r2:r1] = NaN
      r4 = (abs(ymin-8)/resolution)
      r5 = (abs(xmin+77.75)/resolution)
      sst[r5:length(lon),r4:r2] = NaN
      r6 = (abs(ymin-8.75)/resolution)
      r7 = (abs(xmin+82.5)/resolution)
      r8 = (abs(xmin+80)/resolution)
      sst[r7:r8,r6:r2] = NaN
      ncvar_put(nc, "sst", sst)
      nc_sync(nc)
      nc_close(nc)
      rm(lat, lon, sst, xmin, xmax, ymin, ymax, resolution, r1, r2, r3, r4, r5, r6, r7, r8,
         nc, x)}
    rm(sstfiles)}
  #YOU NEED SOME REPRESENTATION OF ENVIRONMENTAL DEPTH-TEMPERATURE
  # HYCOM DATA
  if(!dir.exists("hycom")){
    dir.create("hycom", recursive = T)
    try(get.env(pdt.udates, filename = 'hycom', type = 'hycom', spatLim = sp.lim, 
              save.dir = paste0(getwd(), "/hycom/")))
    setwd("..")}
  rm(pdt.udates)
  #remove Atlantic from north of Panama
  if(removeATL == T){
    hycomfiles <- list.files(path = "hycom", recursive = T, full.names = T)
    for(x in hycomfiles){
      nc <- nc_open(x, write = T)
      lon <- ncvar_get(nc, 'lon')
      lat <- ncvar_get(nc, 'lat')
      water_temp <- ncvar_get(nc, 'water_temp')
      xmin = lon[1]
      xmax = lon[length(lon)]
      ymin = lat[1]
      ymax = lat[length(lat)]
      resolution = (xmax-xmin)/length(lon)
      r1 = (abs(ymin-15)/resolution)
      water_temp[,r1:length(lat),] = NA
      r2 = (abs(ymin-9.25)/resolution)
      r3 = (abs(xmin-276)/resolution)
      water_temp[r3:length(lon),r2:r1,] = NA
      r4 = (abs(ymin-8)/resolution)
      r5 = (abs(xmin-282.25)/resolution)
      water_temp[r5:length(lon),r4:r2,] = NA
      r6 = (abs(ymin-8.75)/resolution)
      r7 = (abs(xmin-277.5)/resolution)
      r8 = (abs(xmin-280)/resolution)
      water_temp[r7:r8,r6:r2,] = NA
      water_temp <- (water_temp-20)/0.00100000004749745
      ncvar_put(nc, "water_temp", water_temp)
      nc_sync(nc)
      nc_close(nc)
      rm(lat, lon, water_temp, xmin, xmax, ymin, ymax, resolution, r1, r2, r3, r4, r5, r6,
         r7, r8, nc, x)}
    rm(hycomfiles)}
  # CALCULATE LIKELIHOODS
  #.par functions are the same calculations as those lacking .par, 
  #except they have been parallelized to leverage multiple CPUs
  locs.grid <- setup.locs.grid(sp.lim)
  
  # LIGHT-BASED LIKELIHOODS
  #L.1 <- calc.srss(light, locs.grid = locs.grid, dateVec = dateVec, res=0.25) 
  # if trying to use raw light levels, not currently recommended (v0.2)
  L.1 <- calc.gpe2(locs, locDates, locs.grid = locs.grid, dateVec = dateVec, errEll = F, 
                   gpeOnly = T)
  rm(locs,locDates)
  #library(fields);library(raster)
  #plot(L.1[[12]]); world(add=T)
  
  # SST LIKELIHOODS
  #L.2 <- calc.sst(tag.sst, filename='oisst', sst.dir = sst.dir, dateVec = dateVec, 
  # sens.err = 1)
  L.2 <- try(calc.sst.par(tag.sst, filename = tagnumbers[i], sst.dir = "sst/", 
                          dateVec = dateVec, sens.err = 1))
  if(!exists("L.2")) L.2 <- L.1*0
  rm(tag.sst)
  #good idea to save after these larger calculations in case the next one causes problems
  save.image(paste0(tagnumbers[i], "-HMMoce.RData")) 
  # also good to do garbage collection and kill any straggling processes that are running
  gc()
  closeAllConnections()
  
  # PDT LIKELIHOODS
  # OCEAN HEAT CONTENT (INTEGRATED PDTS)
  L.3 <- try(calc.ohc.par(pdt, filename = 'hycom', ohc.dir = "hycom/", dateVec = dateVec, 
                          isotherm = '', use.se = F))
  if(!exists("L.3")) L.3 <- L.1*0
  # good idea to save after these larger calculations in case the next one causes problems
  save.image(paste0(tagnumbers[i], "-HMMoce.RData")) 
  #also good to do garbage collection and kill any straggling processes that are running
  gc()
  closeAllConnections() 
  
  # WORLD OCEAN ATLAS-BASED LIKELIHOODS
  #L.4 <- calc.woa.par(pdt, ptt=ptt, woa.data = woa.quarter, sp.lim=sp.lim, focalDim = 9, 
  # dateVec = dateVec, use.se = T)
  if(!exists("L.4")) L.4 <- L.1*0
  # good idea to save after these larger calculations in case the next one causes problems
  # save.image()
  # also good to do garbage collection and kill any straggling processes that are running
  # gc()
  # closeAllConnections() 
  
  # HYCOM PROFILE BASED LIKELIHOODS
  L.5 <- try(calc.hycom.par(pdt, filename = 'hycom', "hycom/", focalDim = 9, 
                            dateVec = dateVec, use.se = T))
  if(!exists("L.5")) L.5 <- L.1*0
  # good idea to save after these larger calculations in case the next one causes problems
  save.image(paste0(tagnumbers[i], "-HMMoce.RData"))
  # also good to do garbage collection and kill any straggling processes that are running
  gc()
  closeAllConnections()
  #save.image('~/ebs/example.rda')
  
  # PREPARE TO RUN THE MODEL
  #Only likelihood outputs calculated above included
  L.rasters <- try(mget(ls(pattern = "L\\.[0-9]{1}"))) 
  resamp.idx <- try(which.max(lapply(L.rasters, FUN = function(x) raster::res(x)[1])))
  L.res <- try(resample.grid(L.rasters, L.rasters[[resamp.idx]]))
  
  # Figure out appropriate L combinations
  # use this if you have a vector (likVec) indicating which likelihoods you are calculating
  # for example, likVec <- c(1,2,5) for light, sst, and hycom likelihoods
  if (length(likVec) > 2){
    L.idx <- c(utils::combn(likVec, 2, simplify = F), 
               utils::combn(likVec, 3, simplify = F))
  }else{
    L.idx <- utils::combn(likVec, 2, simplify = F)}
  
  # GOOD IDEA TO AND SAVE THINGS UP
  save.image(paste0(tagnumbers[i], "-HMMoce.RData"))
  
  # RUN THE MODEL
  if(!dir.exists("modeloutput")){
    dir.create("modeloutput", recursive = T)}
  setwd("modeloutput/")
  
  iniloc <- TagPopLocs %>% filter(Tag == tagnumbers[i]) %>% dplyr::select(Date, Lat, Lon) %>% 
    mutate(Date = lubridate::parse_date_time(Date, "dmy")) %>% 
    separate(Date, c("year", "month", "day"), sep = "-") %>% 
    dplyr::select(day, month, year, everything())
  resultsHMMoce <- data.frame()
  ptt <- tagnumbers[i]
  
  # RUN THE MODEL
  try(for (tt in run.idx){
    for (bnd in bndVec){
      for (j in parVec){
        runName <- paste0(ptt, '_idx', tt, '_bnd', bnd, '_par', j)
        
        library(raster)
        library(HMMoce)
        # COMBINE LIKELIHOOD MATRICES
        # L.idx combination indicates likelihood surfaces to consider
        L <- make.L(L1 = L.res[[1]][L.idx[[tt]]],
                    L.mle.res = L.res$L.mle.res, dateVec = dateVec,
                    locs.grid = locs.grid, iniloc = iniloc, bathy = bathy,
                    pdt = pdt)
        L.mle <- L$L.mle
        L <- L$L
        g <- L.res$g
        g.mle <- L.res$g.mle
        lon <- g$lon[1,]
        lat <- g$lat[,1]
        
        # GET MOVEMENT KERNELS AND SWITCH PROB FOR COARSE GRID
        par0 <- makePar(migr.spd = j, grid = g.mle, L.arr = L.mle, p.guess = c(.7, .8), 
                        calcP = T)
        P.final <- par0$P.final
        
        # GET MOVEMENT KERNELS AND SWITCH PROB FOR FINER GRID
        par0 <- makePar(migr.spd = j, grid = g, L.arr = L, p.guess = c(.7, .8), calcP = F)
        K1 <- par0$K1
        K2 <- par0$K2
        
        # RUN THE FILTER STEP
        if(!is.na(bnd)){
          f <- hmm.filter(g, L, K1, K2, maskL = T, P.final, minBounds = bnd)
          maskL.logical <- TRUE
        }else{
          f <- hmm.filter(g, L, K1, K2, P.final, maskL = F)
          maskL.logical <- FALSE}
        nllf <- -sum(log(f$psi[f$psi>0])) # negative log-likelihood
        
        # RUN THE SMOOTHING STEP
        s <- hmm.smoother(f, K1, K2, L, P.final)
        
        # GET THE MOST PROBABLE TRACK
        tr <- calc.track(s, g, dateVec, iniloc)
        
        if(plottracks == TRUE){
          plotHMM(s, tr, dateVec, ptt = runName, behav.pts = T, save.plot = T)}
        
        # WRITE OUT RESULTS
        outVec <- matrix(c(ptt = ptt, minBounds = bnd, migr.spd = j,
                           Lidx = paste(L.idx[[tt]], collapse=''), P1 = P.final[1,1], 
                           P2 = P.final[2,2], spLims = sp.lim[1:4], 
                           resol = raster::res(L.rasters[[resamp.idx]]),
                           maskL = maskL.logical, NLL = nllf, name = runName), ncol = 15)
        names(outVec) <- c('ptt','bnd','migr.spd','Lidx','P1','P2','spLims-LongLow',
                           'spLims-LongHigh','spLims-LatLow','spLims-LatHigh','resol-1',
                           'resol-2','maskL','nll','name')
        res <- list(outVec = outVec, s = s, g = g, tr = tr, dateVec = dateVec, 
                    iniloc = iniloc, grid = raster::res(L.res[[1]]$L.5)[1])
        save(res, file = paste0(runName, '-HMMoce.rda'))
        #hmm.diagnose(res, L.idx, L.res, dateVec, locs.grid, iniloc, bathy, pdt, plot=T)
        outVec <- t(unlist(outVec)) %>% as.data.frame()
        resultsHMMoce <- rbind(resultsHMMoce, outVec)
        print(runName)} # parVec loop
      gc()} # bndVec loop
    gc()}) # L.idx loop
  # #Arranging results based on nll values (from lowest to largest)
  # resultsHMMoce <- resultsHMMoce %>% 
  #   mutate(nll = as.numeric(as.character(nll))) %>% arrange(nll)
  #Saving results in csv format
  write.csv(resultsHMMoce, file = 'HMMoce_results.csv', row.names = F)
  rm(list = c('L.rasters','L.res','L.idx', paste0("L.", rep(1:5))))
  setwd("../../..")}
