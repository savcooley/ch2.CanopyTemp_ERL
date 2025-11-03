#######   Overview & Install packages ######

#' Goal: Quantify Land Surface Temperature summary statistics (5th percentile, Q25, Q50 (median), Q75, 95th percentile) 
#'        for boosted gradient trees land cover classification of Mato Grosso, Brazil
#'        
#'        I downloaded all Landsat 8 thermal (band10) and QA images available for June, July and Aug 2019 in my study area.  
#'        "path" is the variable storing the file path location of all of these images. 
#'        For any given location in the study area, there is a Landsat overpass every 16 days.
#'        (With the launch of Landsat 9 in 2021 this repeat time for TIR data is now every 8 days).
#'        
#'        This R script creates a set of summary statistic rasters computed based on 
#'        all available spatially overlapping Landsat scenes. Spatially overlapping scenes can be 
#'        determined based on the WRS path and row information stored in each image file name (see note below on naming convention).   
#'        Summary statistics of interest include the following: percentiles 5, 25, 50 (median), 75, 95; mean, sd, n (number of pixels used in calculation for each pixel) 

library(raster)
library(dplyr)
library(purrr)

library(data.table); library(sf); library(raster); library(terra); 
library(FNN); library(ggplot2); library(dplyr); 
library(tidyverse); library(ncdf4); 
library(colorspace); library(ggplot2);library(RColorBrewer)
path="/Users/savannah/Documents/Ch2/MT_classification_data/Landsat_LST" 

#'        ####################
#'        Detailed description: 
#'        
#'        Step 1. Prep data with QA filter & convert temperature to Celsius
#'        
#'        For each Landsat 8 thermal band10 image in directory:
#'        Mask out clouds and noise by applying a QA filter with corresponding the QA_PIXEL image 
#'        (pixel quality of 1 = no data, pixel quality >~50000 is a cloud / noise. So good data seem to have values >1 and <50000)
#'        Convert from DN to Kelvin: tempK = 0.00341802 * DN + 149.0 
#'        Convert from Kelvin to Celsius: tempC = tempK − 273.15
#'        Write processed raster to file in new folder
#'        
#'        Step 2. Create a single raster for each summary statistic of interest. 
#'        
#'        For each processed band10 image: 
#'        Find all spatially overlapping band10 images for the current Landsat scene - these are all captured 16 days apart. 
#'        Compute summary statistics for each unique path-row combination. 
#'        Summary statistics will include the following: percentiles 5, 25, 50 (median), 75, 95; mean, sd, n (number of pixels used in calculation for each pixel) 
#'        
#'        ####################

#' Notes:
#' 
#' The naming convention for Landsat Collections Level-1 scenes:

# The Landsat Collection 1 Level-1 product identifier includes the Collection processing levels, processing date, collection number, and collection tier category:
#   
#   LXSS_LLLL_PPPRRR_YYYYMMDD_yyyymmdd_CC_TX
# 
# Where:
#   
# L = Landsat
# X = Sensor (“C”=OLI/TIRS combined, “O”=OLI-only, “T”=TIRS-only, “E”=ETM+, “T”=“TM, “M”=MSS)
# SS = Satellite (”07”=Landsat 7, “08”=Landsat 8)
# LLL = Processing correction level (L1TP/L1GT/L1GS)
# PPP = WRS path
# RRR = WRS row
# YYYYMMDD = Acquisition year, month, day
# yyyymmdd - Processing year, month, day
# CC = Collection number (01, 02, …)
# TX = Collection category (“RT”=Real-Time, “T1”=Tier 1, “T2”=Tier 2)
# 
# Example:  LC08_L1GT_029030_20151209_20160131_01_RT
# 
# Means: Landsat 8; OLI/TIRS combined; processing correction level L1GT; path 029; row 030; acquired December 9, 2015; processed January 31, 2016; Collection 1; Real-Time

#####     Step 1. Process data: QA and conversion to Celsius     #####

# Function to process Landsat 8 thermal band 10 image
process_landsat8 <- function(lst_image_path, QA_path) {
  # Read band 10 and QA images
  band10 <- raster(lst_image_path) # Testing: raster(band10_files[86]) #
  qa <- raster(QA_path) # Testing: raster(qa_files[86]) #
  
  # Apply QA filter: Mask out clouds and noise with QA filter using the QA_PIXEL image 
  # (pixel quality of 1 = no data, pixel quality >~50000 is a cloud / noise. So good data seem to have values >1 and <50000)
  # Testing with: [85] "/Users/savannah/Documents/Ch2/MT_classification_data/Landsat_LST/B10_LST/LC08_L2SP_225068_20190608_20200828_02_T1_ST_B10.TIF"
  
  good_qa <- qa > 1 & qa < 50000
  band10_filtered <- band10*good_qa # mask(band10, good_qa) didn't work as expected
  
  # Convert DN to Kelvin and then to Celsius
  tempK <- 0.00341802 * band10_filtered + 149.0
  tempC <- tempK - 273.15
  
  # need additional QA bc good_qa misses some clouds
  # all values > -Inf and <= 0.1 become NA, etc. 
  rclmat <- matrix(c(-Inf, 20, NA), ncol=3, byrow=TRUE) # Check MT min and max temp for June-Aug
  tempC <- reclassify(tempC, rclmat)
  
  return(tempC)
}

# Step 1: Prep data with QA filter & convert temperature to Celsius
processed_folder <- "/Users/savannah/Documents/Ch2/MT_classification_data/Landsat_LST/_processed_LST" 
band10_files <- list.files(paste(path,"/B10_LST", sep = ""), full.names = TRUE) # path, pattern = "_B10.TIF$"
qa_files <- list.files(paste(path,"/QA_PIXEL", sep = ""), full.names = TRUE) # path, pattern = "_QA_PIXEL.TIF$"

band10_files <- band10_files[85:96] # then 1:85, then 97:327 
qa_files <- qa_files[85:96]
processed_files <- mapply(function(band10, qa) {
  processed <- process_landsat8(band10, qa)
  file_name <- gsub("_B10.TIF", "_processed.tif", basename(band10))
  writeRaster(processed, file.path(processed_folder, file_name), overwrite = TRUE)
  file.path(processed_folder, file_name)
}, band10_files, qa_files, SIMPLIFY = FALSE)

#####     Step 2. Compute statistics: Create a single raster for each summary statistic     #####

get_path_row <- function(file) {
  basename(file) %>%
    gsub("^.*_(\\d{3})(\\d{3})_.*$", "\\1_\\2", .)
}
processed_folder_stats<-"/Users/savannah/Documents/Ch2/MT_classification_data/Landsat_LST/_processed_LST/_Sum_stats_v2"  
processed_folder <- "/Users/savannah/Documents/Ch2/MT_classification_data/Landsat_LST/_processed_LST"  
#summary_stats <- c("p5", "p25", "p50", "p75", "p95") #, "mean", "sd", "n")

lst_files <- list.files(path=processed_folder, pattern="*.tif", full.names = TRUE, recursive = F) 
path_row = get_path_row(lst_files)
path_row<- gsub("_", "", path_row); path_row<- as.factor(path_row)

for(i in 16:16){ #3:14, then 17:length(levels(path_row)) 
  print(levels(path_row)[i])
  scene_i_files <- list.files(path=processed_folder, pattern=paste("*",levels(path_row)[i], sep=""), full.names = TRUE, recursive = F) 
  # raster_stack <- stack(scene_i_files) # Need to resample rasters to address: raster stack() Error in compareRaster(rasters) : different extent
  
  r1<- raster(scene_i_files[1])
  raster_stack<-stack(r1)
  for(j in 2:length(scene_i_files)){
    r2<- raster(scene_i_files[j])
    r2<- resample(r2, r1) # takes ~1 min
    raster_stack <- stack(raster_stack,r2) 
  }
  
  # Create a single raster for each summary statistic
  # Computing one stat at a time and saving files individually could be faster.
  #Sys.time() # [1] "2024-04-03 17:02:48 EDT"
  p5 = calc(raster_stack, fun = function(x) quantile(x, 0.05, na.rm = TRUE))
  writeRaster(p5, file.path(processed_folder_stats, paste0(levels(path_row)[i], "_", "p5", ".tif")), overwrite = TRUE)
  #Sys.time() "2024-04-03 18:40:55 EDT"
  
    #p25 = calc(raster_stack, fun = function(x) quantile(x, 0.25, na.rm = TRUE))
  p50 = calc(raster_stack, fun = function(x) quantile(x, 0.50, na.rm = TRUE))
  writeRaster(p50, file.path(processed_folder_stats, paste0(levels(path_row)[i], "_", "p50", ".tif")), overwrite = TRUE)
  
    #p75 = calc(raster_stack, fun = function(x) quantile(x, 0.75, na.rm = TRUE))
  p95 = calc(raster_stack, fun = function(x) quantile(x, 0.95, na.rm = TRUE))
  writeRaster(p95, file.path(processed_folder_stats, paste0(levels(path_row)[i], "_", "p95", ".tif")), overwrite = TRUE)
  
    #mean = calc(raster_stack, fun = mean, na.rm = TRUE),
    #sd = calc(raster_stack, fun = sd, na.rm = TRUE),
    #n = calc(raster_stack, fun = function(x) sum(!is.na(x)))
  
}

###############
# Identify which Landsat scene (path_row) corresponds to the Feliz Natal site
#' Feliz Natal Municipality center lat, lon: -12.050804046818172, -54.29033607522494
#' WRS Path: 225, WRS Row: 069. Process i = 16: [16] "225069". Also i =15: "225068"
#' 
#' WRS Path: 226, WRS Row: 068. Process i=23. Also i = 24 "226069" 



##############

# Create a single raster for each summary statistic

# This steps took a very long time (start time for i=1 of ~3:30pm, end time: ). 
# For i=2, I could test computing one stat at a time and saving files individually 

summary_list <- list(
  p5 = calc(raster_stack, fun = function(x) quantile(x, 0.05, na.rm = TRUE)),
  #p25 = calc(raster_stack, fun = function(x) quantile(x, 0.25, na.rm = TRUE)),
  p50 = calc(raster_stack, fun = function(x) quantile(x, 0.50, na.rm = TRUE)),
  #p75 = calc(raster_stack, fun = function(x) quantile(x, 0.75, na.rm = TRUE)),
  p95 = calc(raster_stack, fun = function(x) quantile(x, 0.95, na.rm = TRUE))
  #mean = calc(raster_stack, fun = mean, na.rm = TRUE),
  #sd = calc(raster_stack, fun = sd, na.rm = TRUE),
  #n = calc(raster_stack, fun = function(x) sum(!is.na(x)))
) 

# Write summary rasters to file
for (j in 1:length(summary_list)) {
  file_name <- paste0(levels(path_row)[i], "_", summary_stats[j], ".tif")
  writeRaster(summary_list[j], file.path(processed_folder, file_name), overwrite = TRUE)
}

###


# Step 2: Create a single raster for each summary statistic
summary_stats <- c("p5", "p25", "p50", "p75", "p95", "mean", "sd", "n")

processed_files_df <- data.frame(file = unlist(processed_files), stringsAsFactors = FALSE) %>%
  mutate(path_row = get_path_row(file))

summary_rasters <- processed_files_df %>%
  split(.$path_row) %>%
  map(function(df) {
    raster_stack <- stack(df$file)
    raster_list <- as.list(raster_stack)
    
    summary_list <- list(
      p5 = calc(raster_stack, fun = function(x) quantile(x, 0.05, na.rm = TRUE)),
      p25 = calc(raster_stack, fun = function(x) quantile(x, 0.25, na.rm = TRUE)),
      p50 = calc(raster_stack, fun = function(x) quantile(x, 0.50, na.rm = TRUE)),
      p75 = calc(raster_stack, fun = function(x) quantile(x, 0.75, na.rm = TRUE)),
      p95 = calc(raster_stack, fun = function(x) quantile(x, 0.95, na.rm = TRUE)),
      mean = calc(raster_stack, fun = mean, na.rm = TRUE),
      sd = calc(raster_stack, fun = sd, na.rm = TRUE),
      n = calc(raster_stack, fun = function(x) sum(!is.na(x)))
    )
    
    return(summary_list)
  })

# Write summary rasters to file
for (path_row in names(summary_rasters)) {
  for (stat in names(summary_rasters[[path_row]])) {
    file_name <- paste0(path_row, "_", stat, ".tif")
    writeRaster(summary_rasters[[path_row]][[stat]], file.path(processed_folder, file_name), overwrite = TRUE)
  }
}