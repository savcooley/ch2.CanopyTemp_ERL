#######   Notes     #######   

#' GOAL: This script determines temperature differences among forest degradation types. 
#' 1. Pre-process ECOSTRESS LST with QC product
#' 2. Compute p90 and p95 percentile LST stats (dry season months 6/1 - 8/31 for 2018, 2019, 2020)
#' 3. Hard classification of degradation type probability data from Pinage et al 20203
#' 4. Box plots of LST across forest degradation classes (burned, logged, intact)
#' 5. Add covariates, including distance to forest edge
#' 6. Statistically test for differences in LST among classes  

# ECOSTRESS Level2 LST documentation: https://lpdaac.usgs.gov/documents/423/ECO2_User_Guide_V1.pdf

# Bits 1&0 Mandatory QA flags 
#' 00 = Pixel produced, best quality
#' 01 = Pixel produced, nominal quality. Either one or more of the following conditions are met:
#' 1. Emissivity in both bands 4 and 5 < 0.95, i.e. possible cloud contamination
#' 2. Low transmissivity due to high water vapor loading (<0.4), check PWV values and error estimates
#' 3. Pixel falls on missing scan line in bands 1&5, and filled using spatial neural net. Check error estimates.
#' Recommend more detailed analysis of other QC information
#' 10 = Pixel produced, but cloud detected
#' 11 = Pixel not produced due to missing/bad data, or TES divergence, user should check data quality flag bits.

# Bits 15 & 14 LST accuracy 
#' 00 = >2 K (Poor performance)
#' 01 = 1.5 - 2 K (Marginal performance)
#' 10 = 1 - 1.5 K (Good performance) 
#' 11 = <1 K (Excellent performance)

# BIO1 = Annual Mean Temperature 
# BIO2 = Mean Diurnal Range (Mean of monthly (max temp - min temp))
# BIO3 = Isothermality (BIO2/BIO7) (×100)
# BIO4 = Temperature Seasonality (standard deviation ×100)
# BIO5 = Max Temperature of Warmest Month
# BIO6 = Min Temperature of Coldest Month
# BIO7 = Temperature Annual Range (BIO5-BIO6)
# BIO8 = Mean Temperature of Wettest Quarter
# BIO9 = Mean Temperature of Driest Quarter
# BIO10 = Mean Temperature of Warmest Quarter
# BIO11 = Mean Temperature of Coldest Quarter
# BIO12 = Annual Precipitation
# BIO13 = Precipitation of Wettest Month
# BIO14 = Precipitation of Driest Month
# BIO15 = Precipitation Seasonality (Coefficient of Variation)
# BIO16 = Precipitation of Wettest Quarter
# BIO17 = Precipitation of Driest Quarter
# BIO18 = Precipitation of Warmest Quarter
# BIO19 = Precipitation of Coldest Quarter

#' 3. Create raster of hard classification categories with probability data from Pinage et al 20203
#FNA_Probs<- raster(paste(root,"Pinage_class_data/FNA_Probabilities.tif", sep=""))
#plot(FNA_Probs) 
# plot each band 

#' => Bands of the probability raster (each band represents the probability for one forest class):
#' 1) Probability of the grid cell being intact forest
#' 2) Probability of the grid cell being logged forest
#' 3) Probability of the grid cell being burned forest
#' 
#' Hard classification categories are: 
#' grid cell = 1) intact forest if Band 1 > 75
#' grid cell = 2) logged forest if Band 2 > 75
#' grid cell = 3) burned forest if Band 3 > 75
#' NOTE: if a grid cell has probability >75 for two or more classes, then assign it to the class with the highest probability 
#' such that the result is a single band raster with pixel values of 1, 2 and 3 based on degradation classes of burned (3), logged (2), intact (1).


# Longo et al 2019: Differences along biomass gradients exceed 9°C in degraded forests 
# during the dry season at all regions except GYF (Figure 6)
# 
# LST ~ AGB [kg C m-2] 
# LST ~ RH95 [m]
# LST ~ RH50 [m] # or other metric reflecting distribution of biomass along the vertical 3D forest structure profile
# 
# Fig 6. Feliz Natal curve saturates at ~39 C at 3 kg C m-2 or less for dry season months
# Plot LST vs structure colored by forest degradation type

#######   Install packages     #######   
library(raster)
library(dplyr)
library(tidyr)
library(ggplot2)
library(spatstat)  # for distance calculations
library(sp)
library(spdep)
library(nlme); library(lme4)
library('sjPlot') 
library("reshape2")

root="/Users/savannah/Documents/Ch2/FNA_pilot/"

set_theme(base = theme_classic(), #theme_minimal(), #axis.tickslen = 0, # hides tick marks
          axis.title.size = 2.1, axis.textsize = 1.2, legend.size = 1.7,  legend.title.size = 1.8,
          geom.label.size = 2.5, title.size = 2.1)

#######   1. ECOSTRESS Data Pre-processing     #######   

# qc_binary = 0 IF Bits 1 & 0 = 00 (Pixel produced, best quality, may or may not have clouds)
# AND if Bits 15 & 14 = 10 or 11 (LST accuracy of Good and Excellent performance)
# qc_binary=1 otherwise. 

# ECOSTRESS L2 User Guide Table 7 details the data sets included in the L2_CLOUD output. Users can interpret the data sets as follows:
# 1. Cloud_confidence contains the results of the brightness temperature LUT test with: 
#     0 = confident clear, 1 =  probably clear, 2 = probably cloudy, and 3 = confident cloudy.

# 2. Cloud_final contains a final cloud mask (1=cloud, 0=clear) based on the following criteria:
#     a. For elevations < 2km, cloud = probably cloud + confident cloudy pixels
#     b. For elevations >2km, cloud = confident cloudy pixels

# Note: Bit 0 is the least significant bit, i.e. the right-most digit.

LST_files <- list.files(paste(root,"ECOSTRESS_LST_vDec2024/LST", sep=""), full.names = T)
LST_files <- c("/Users/savannah/Documents/Ch2/FNA_pilot/ECOSTRESS_LST_vDec2024/LST/ECO_L2_LSTE.002_LST_doy2018243173803_aid0001.tif",
               "/Users/savannah/Documents/Ch2/FNA_pilot/ECOSTRESS_LST_vDec2024/LST/ECO_L2_LSTE.002_LST_doy2018212164454_aid0001.tif")

for(i in 1:length(LST_files)){ #LST_files ){ ## i = 14 error: Error in if (x == "" | x == ".") { : argument is of length zero
  
  r_LST <- raster(LST_files[i]) #;plot(r_LST)
  
  cloud_file <- list.files(paste(root,"ECOSTRESS_LST_vDec2024/cloud_mask", sep=""), pattern = substr(LST_files[i], 89, 103), full.names = T)
  r_cloudMask <- raster(cloud_file)
  
  # Apply cloud mask
  lst_masked <- mask(r_LST, r_cloudMask, maskvalue = 1)
  
  QC_file <- list.files(paste(root,"ECOSTRESS_LST_vDec2024/QC", sep=""), pattern = substr(LST_files[i],  89, 103), full.names = T)
  r_QC <- raster(QC_file)
  
  # Create binary mask based on QC criteria
  qc_binary <- overlay(r_QC, fun = function(x) {
    # Check Bits 1 & 0 (00 for best quality)
    quality_bits <- bitwAnd(x, 3)
    
    # Check Bits 15 & 14 (01, 10 or 11 for Marginal, Good or Excellent LST accuracy)
    accuracy_bits <- bitwShiftR(bitwAnd(x, 49152), 14)
    
    # Combine conditions
    ifelse(quality_bits == 0 & (accuracy_bits == 1 | accuracy_bits == 2 | accuracy_bits == 3), 0, 1)
  })
  #plot(qc_binary)
  
  # Apply QC mask to LST raster
  lst_masked <- mask(lst_masked, qc_binary, maskvalue = 1)
  
  # Apply scale factor 0.02 and subtract 273.15 for Kelvin to Celsius conversion
  lst_masked_degC <- (lst_masked * 0.02) - 273.15 # vDec2024 data do not fall within expected temp range - why??
  plot(lst_masked_degC)
  writeRaster(lst_masked_degC, paste(root,"ECOSTRESS_LST_vDec2024/LST_processed_v3/LST_", substr(LST_files[i],89, 103), ".tif", sep="") )
}

#######   3. Hard classification with FNA data     #######   

# Process FNA_Probs. Hard classification categories are: 
#' grid cell = 1) intact forest if Band 1 > 75
#' grid cell = 2) logged forest if Band 2 > 75
#' grid cell = 3) burned forest if Band 3 > 75
#' NOTE: if a grid cell has prob>75 for two or more classes, then assign it to the class with the highest probability 
 
FNA_Probs <- stack(paste0(root, "Pinage_class_data/FNA_Probabilities.tif"))

hard_classify <- function(x) {
  if (all(is.na(x))) return(NA)
  if (any(x > 75, na.rm = TRUE)) {
    return(which.max(ifelse(x > 75, x, -Inf)))
  } else {
    return(which.max(x))
  }
}
hard_class <- calc(FNA_Probs, hard_classify)
plot(hard_class, legend = FALSE, col = c("green4","orange","red3")) # rev(topo.colors(3))
legend("topright", legend = c("Intact", "Logged", "Burned"), fill = c("green4","orange","red3"))

writeRaster(hard_class, paste0(root, "Pinage_class_data/FNA_HardClass.tif"))

#######   5. Plot Day of Year vs. Canopy temperature     #######   

# Process LST files
LST_files <- list.files(paste(root, "ECOSTRESS_LST/LST_processed_v2", sep = ""), full.names = TRUE, pattern = "*.tif")
LST_files_short <- list.files(paste(root, "ECOSTRESS_LST/LST_processed_v2", sep = ""), full.names = FALSE, pattern = "*.tif")
df_all<- data.frame()

#' For each LST file
for(i in 1:length(LST_files)){
  
  r_LST <- raster(LST_files[i])
  
  # Ensure rasters have the same crs (extent and resolution ok to differ)
  r_LST <- projectRaster(r_LST, hard_class)
  
  # Spatial join LST stats with classification data
  lst_class <- stack(r_LST, hard_class)
  lst_df <- as.data.frame(lst_class, xy = TRUE)
  colnames(lst_df) <- c("lon", "lat", "r_LST", "forestClass")
  lst_df_sub<- lst_df %>%
    mutate(forestClassCh = case_when((forestClass == 1) ~ as.character("Intact"),
                                     (forestClass == 2) ~ as.character("Logged"),
                                     (forestClass == 3) ~ as.character("Burned"),
                                     TRUE ~ "Other" )) # (forestClass == NA) ~ "Other", # TRUE ~ NA )) 
  
  #' Compute the mean, median, p75, p90, p95 and p99 percentiles and n observations for each forestClass factor level
  #lst_df_sub$forestClassCh <- as.factor(lst_df_sub$forestClassCh)
  
  lst_df_sub<- lst_df_sub[lst_df_sub$forestClassCh!="Other",]
  result <- lst_df_sub %>% #group_by(forestClassCh) %>% 
    summarise(n=n(), meanLST=mean(r_LST, na.rm=T), medLST=median(r_LST, na.rm=T), p75_LST=as.numeric(quantile(r_LST, na.rm=T)[4]),
              maxLST=max(r_LST, na.rm=T), minLST=min(r_LST, na.rm=T), sdLST=sd(r_LST, na.rm=T))
  
  #' Add in entries for acquisition date (substring of the file name) 
  result$aquDate<- substr(LST_files_short[i], 8, 20) # YYYYdddHHMMSS - last 6 characters are Time of Acquisition (HHMMSS) (in UTC)
    
  #result <-dcast(melt(result, id.vars=c( "forestClassCh")), aquDate~variable+forestClassCh)
  df_all<- rbind(df_all, result)
   
  }
df_all$DoY<- substr(df_all$aquDate, 5, 7)
write.csv(df_all, paste(root, "ECOSTRESS_LST/LST_processed/LST_df_all_DoY_stats_v2.csv", sep = ""))

df_all<- na.omit(df_all)
df_all$DoY<- as.numeric(df_all$DoY)
ggplot(df_all, aes(x = DoY, y = medLST)) +
  geom_point() +
  labs(x = "Day of year", y = "Median Canopy Temp (degrees C)") #, title = "Maximum Canopy Temperature") 

ggplot(df_all, aes(x = DoY, y = maxLST)) +
  geom_point() +
  labs(x = "Day of year", y = "Max. Canopy Temp (degrees C)") #, title = "Maximum Canopy Temperature") 

#######   6. Doughty Fig 1d analog - Canopy temperature distribution     #######   

# Process LST files
LST_files <- list.files(paste0(root, "ECOSTRESS_LST/LST_processed_v2"), full.names = TRUE, pattern = "*.tif")
LST_files_short <- list.files0(paste(root, "ECOSTRESS_LST/LST_processed_v2"), full.names = FALSE, pattern = "*.tif")
df_all<- data.frame()

i=2 # 2018212164454 - looks fine 

#' For each LST file
for(i in 1:length(LST_files)){
  
  r_LST <- raster(LST_files[i])
  
  # Ensure rasters have the same crs (extent and resolution ok to differ)
  r_LST <- projectRaster(r_LST, hard_class)
  
  # Spatial join LST stats with classification data
  lst_class <- stack(r_LST, hard_class)
  lst_df <- as.data.frame(lst_class, xy = TRUE)
  colnames(lst_df) <- c("lon", "lat", "r_LST", "forestClass")
  lst_df_sub<- lst_df %>%
    mutate(forestClassCh = case_when((forestClass == 1) ~ as.character("Intact"),
                                     (forestClass == 2) ~ as.character("Logged"),
                                     (forestClass == 3) ~ as.character("Burned"),
                                     TRUE ~ "Other" )) # (forestClass == NA) ~ "Other", # TRUE ~ NA )) 
  
  #' Compute the mean, median, p75, p90, p95 and p99 percentiles and n observations for each forestClass factor level
  #lst_df_sub$forestClassCh <- as.factor(lst_df_sub$forestClassCh)
  
  lst_df_sub<- lst_df_sub[lst_df_sub$forestClassCh!="Other",]
  result <- lst_df_sub %>% #group_by(forestClassCh) %>% 
    summarise(n=n(), meanLST=mean(r_LST, na.rm=T), medLST=median(r_LST, na.rm=T), p75_LST=as.numeric(quantile(r_LST, na.rm=T)[4]),
              maxLST=max(r_LST, na.rm=T), minLST=min(r_LST, na.rm=T), sdLST=sd(r_LST, na.rm=T))
  
  #' Add in entries for acquisition date (substring of the file name) 
  result$aquDate<- substr(LST_files_short[i], 8, 20) # YYYYdddHHMMSS - last 6 characters are Time of Acquisition (HHMMSS) (in UTC)
  
  #result <-dcast(melt(result, id.vars=c( "forestClassCh")), aquDate~variable+forestClassCh)
  df_all<- rbind(df_all, result)
  
}

write.csv(lst_df_sub, paste0(root, "ECOSTRESS_LST/LST_processed_v2/",
                                    "_forest_LST_UTM21S_doy",substr(LST_files_short[i], 8, 20), ".csv"))

ggplot(lst_df_sub, aes(x = factor(forestClassCh), y = r_LST, fill= forestClassCh)) +
  geom_boxplot() + scale_fill_manual(values=c("red","green","blue"))+
  theme(legend.position="none") +
  labs(x = "Forest Class", y = "Canopy Temperature (degrees C)" ) #, 

#r_lst_df <- rasterFromXYZ(lst_df)  #Convert first two columns as lon-lat and third as value                

# Apply forest mask to LST raster (include only forest pixels)
m <- c(0, 4, 1,  -Inf, 0, NA,  4, Inf, NA)
rclmat <- matrix(m, ncol=3, byrow=TRUE)
class_binary <- reclassify(hard_class, rclmat)
lst_masked <- class_binary*r_LST
writeRaster(lst_masked, paste(root, "ECOSTRESS_LST/lst_forestOnly_doy2018212164454.tif", sep = ""))

#   Doughty Fig 1d analog - Canopy temperature distribution
ggplot(lst_df_sub, aes(x=r_LST, fill=factor(forestClassCh))) +
  geom_histogram(alpha = 0.7, position="dodge", colour=NA) + scale_fill_manual(values=c("red","green","blue"))+
  #theme(legend.position="none") +
  labs(x = "Canopy Temperature (degrees C)", y = "Count (n pixels)" ) 



###   Notes on methods & peak canopy temp of specific ECOSTRESS scenes     #######   

'''

Data pre-processing steps for all_vars.csv

1. FNA_probabilities.tif: Raster to vector then Dissolve 
2. Prepare for addressing edge effects: Buffer -0.0002250002 degrees = -25 m  / 111111 degrees per m

3. Convert GEDI L2A to points: convert tifs to df, filter out NAs, merge & save as RH_df_all.csv [Orig approach] Raster pixels to points for just RH_98_FNA_pilot.tif 
4. Remove edge effects from GEDI L2A: Intersect RH_98_FNA_pilot.tif with Buffer_FNA_probabilities.geojson  (takes ~20 min to run)
5.Compile final data set with all covariates: Run Sample Raster Values (annoying bc no batch functionality; but SAGA Add Raster Values to Points didnt work, kept getting error:  Loading resulting layers The following layers were not correctly generated. )


5. [Orig approach] Filter out NAs in GEDI L2A: Select by attribute all points !=0 & save as all_vars.csv
6. [Orig approach] Compile final data set with all covariates: Run SAGA Add Raster Values to Points. 
#  (use all_vars.csv and add in all covariate rasters of interest, including ECOSTRESS LST scene w max mean temp )

Statistical explorations
1. Comute spatial variogram, noting the range distance value to use in the next step
2. To minimize spatial autocorrelation and avoid pseudo replicates, Sub-sample points to with spatialEco::subsample.distance. 
3. Run RF variable importance on LST 
 
 
 
 Notes on ECOSTRESS high LST images
 
2020259101050 - strange artifacts influencing temp
2018243173711 - all NA 
2018243173803 - looks fine
2020333151658 - mostly NA (cloud masked) 
2018 212 164454 - looks fine. Picked this one for first test.
2020090152422 - mostly NA (cloud masked)
2020262192004 - looks fine 

'''


#######   1 & 2. Spatial subset of points: minimize sp autocorr     #######   
library(gstat); library(sp)
df=read.csv("/Users/savannah/Documents/Ch2/FNA_pilot/all_vars_FNA.csv")

# convert simple data frame into a spatial data frame object
coordinates(df)= ~ x+y
proj4string(df) <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")

TheVariogram=variogram(rh_98~1, data=df) # takes a few min
plot(TheVariogram) 

TheVariogramModel <- vgm(psill=38, model="Gau", nugget=25, range=0.08)
FittedModel <- fit.variogram(TheVariogram, model=TheVariogramModel)   
# Warning message: No convergence after 200 iterations: try different initial values? 
# I tried a bunch of possibilities. I'm going to approximate the range via visually inspecting plot(TheVariogram) 
# So: psill= ~38, nugget= ~25, range= ~0.08

#plot(TheVariogram, model=FittedModel)

devtools::install_github("jeffreyevans/spatialEco")
library(sp); library(spatialEco)

numSamp = 5000 #1672 # if d = 0.005 degrees (FNA_propoabilities.tif resolution)
#sub.df <- sample.distance(df, n = numSamp, d = 0.005093893051948002702, trace = TRUE)  
df<- st_as_sf(df) 
sub.df <- spatialEco::subsample.distance(df, size = numSamp, d = 0.08) # not be possible to reach nSamp. See note below  
# start time 3:32, took less than an hour to run for 167k pts in df, 
# resulted in nrow(sub.df) of 1672. Weird bc I was expecting fewer. See NOTE below.
# just because you specify a condition for your data does not mean that it can actually be met. Here is an example using the meuse data. The data cannot meet the condition of a 500m minimum sampling distance for greater than ~15 points (n for 50% sample is 78). This is obviously dictated by the configuration of the randomization but n should not vary that much. I added error checking for non-convergence and the function will return the subsample on however many samples can be identified using the given conditions. 

# NOTE: Needed to convert back to df before writing to csv
write.csv(sub.df, "/Users/savannah/Documents/Ch2/FNA_pilot/all_vars_FNA_sampDist.csv")


#######   3. RF var importance plot    ####### 
library(randomForest); library(corrplot)

df <- read.csv("/Users/savannah/Documents/Ch2/FNA_pilot/all_vars_FNA.csv")
sub.df <- df
#sub.df<- read.csv("/Users/savannah/Documents/Ch2/FNA_pilot/all_vars_FNA_sampDist_numSamp5k.csv")

# Remove unnecessary variables
#sub.df$LST_2018212_16 <- NULL
sub.df$LST_2018_243_17 <- NULL
sub.df$LST_2020_262_19 <- NULL
sub.df$MeanMaxLST<- NULL

# Remove rows with missing values & outliers
sub.df <- na.omit(sub.df)
sub.df <- subset(sub.df, agbd_Mg_ha > 0 & agbd_Mg_ha < 600)
hist(sub.df$agbd_Mg_ha)

# Add burned and logged dummy variables
sub.df$Forest_cat <- factor(sub.df$Forest_cat, levels = c("Intact", "Logged", "Burned"))
sub.df$burned <- as.numeric(sub.df$Forest_cat == "Burned")
sub.df$logged <- as.numeric(sub.df$Forest_cat == "Logged")

#sub_dt<- as.data.frame(sapply( sub.df, as.numeric) )
#sub_dt<-as.data.frame(scale(sub_dt, center = T, scale = T)) # To avoid Warning message: In checkConv(attr(opt, "derivs"), opt$par, ctrl = control$checkConv,  : Model is nearly unidentifiable: very large eigenvalue - Rescale variables?

rf <-randomForest(LST_2018212_16 ~.,data=sub.df, ntree=500, na.action = na.exclude)
#print(rf) # If a dependent variable is a factor, classification is assumed, otherwise regression is assumed. If omitted, randomForest will run in unsupervised mode.

# In example code, the number of variables tried at each split is based on the following formula. -1 is used as dataset contains dependent variable as well.
floor(sqrt(ncol(sub.df) - 1))

# The number of variables selected at each split is denoted by mtry in randomforest function. 
# Find the optimal mtry value by Selecting mtry value with minimum out of bag(OOB) error.
#sub.df<-sub.df[!is.na(sub.df$CH4.C.conv_pos),]
m<-na.omit(sub.df)
mtry <- tuneRF(m[1:(length(m)-1)],
               m$LST_2018212_16, ntreeTry=500,
               stepFactor=1.5,improve=0.01, trace=TRUE, plot=FALSE)

best.m <- mtry[mtry[, 2] == min(mtry[, 2]), 1]
#print(mtry); print(best.m) # best mtry=2 

# Parameters in tuneRF function:  
#The stepFactor specifies at each iteration, mtry is inflated (or deflated) by this value
#The improve specifies the (relative) improvement in OOB error must be by this much for the search to continue
#The trace specifies whether to print the progress of the search
#The plot specifies whether to plot the OOB error as function of mtry

#Build model again using best mtry value.
#set.seed(71)
rf <-randomForest(LST_2018212_16 ~.,data=sub.df, mtry=best.m, importance=TRUE,ntree=500, na.action = na.exclude) #; print(rf)
#Evaluate variable importance
m<-importance(rf); m<-data.frame(m)

# Sort the data frame by the "X.IncMSE" column in descending order
sorted_m <- m[order(-m$X.IncMSE), ]; print(sorted_m)

varImpPlot(rf, main = NULL) #"LST Variable Importance")

# Higher the value of mean decrease accuracy or mean decrease gini score , higher the importance of the variable in the model. In the plot shown above, Account Balance is most important variable.
# %IncMSE: Mean Decrease Accuracy - How much the model accuracy decreases if we drop that variable.
# IncNodePurity: Mean Decrease Gini - Measure of variable importance based on the Gini impurity index used for the calculation of splits in trees.

##    Corr plot based on RF results
sub.df<-subset(sub.df, select = c(rownames(sorted_m)[1:13])) # Retrieve the row names of the top 13 rows
sub.df$Forest_cat<- NULL # can't compute corr w factor, not going to turn to dummy var bc I already know high corr exists w Age variable
corrplot.mixed(cor(sub.df, use="complete.obs"), main = "Variable correlations", order = 'AOE') #order = 'AOE' for the angular order of the eigenvectors. # order="alphabet"



#######   4. Gamma regression models    ####### 
sub.df<- read.csv("/Users/savannah/Documents/Ch2/FNA_pilot/all_vars_FNA_sampDist_numSamp5k.csv")
#sub.df$LST_2020212_16<- NULL
sub.df$LST_2020262_19<- NULL
sub.df$LST_2018243_17<- NULL
sub.df$MeanMaxLST<- NULL
#sub.df$MeanMaxLST<- as.numeric(sub.df$MeanMaxLST)
sub_dt<- na.omit(sub.df) # from 5k obs to 4.9k (MeanMaxLST) and 4.7k (with LST_2020262_19)  
sub_dt$Forest_cat <- factor(sub_dt$Forest_cat, levels = c("Intact", "Logged", "Burned" ))
b_sub_dt<- sub_dt
sub_dt<- subset(b_sub_dt, select=c("rh95", "air_temp", "SM_surf_anom", "elevation", "aspect"))
sub_dt<-as.data.frame(scale(sub_dt, center = T, scale = T)) # To avoid Warning message: In checkConv(attr(opt, "derivs"), opt$par, ctrl = control$checkConv,  : Model is nearly unidentifiable: very large eigenvalue - Rescale variables?

sub_dt$Forest_cat<- b_sub_dt$Forest_cat
sub_dt$LST_2020212_16<- b_sub_dt$LST_2020212_16

m1 <- glm(LST_2020212_16 ~ Forest_cat, data= sub_dt, family = Gamma(link = "log"), na.action = na.exclude) 
#, nAGQ=0, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=1000)) ) # only for glmer

m2 <- glm(LST_2020212_16 ~ rh95, data= sub_dt, family = Gamma(link = "log"), na.action = na.exclude)
m3 <- glm(LST_2020212_16 ~ rh95+air_temp, data= sub_dt, family = Gamma(link = "log"), na.action = na.exclude)
m4 <- glm(LST_2020212_16 ~ rh95+air_temp+SM_surf_anom, data= sub_dt, family = Gamma(link = "log"), na.action = na.exclude)
m5 <- glm(LST_2020212_16 ~ rh95+air_temp+SM_surf_anom+elevation, data= sub_dt, family = Gamma(link = "log"), na.action = na.exclude)
m6 <- glm(LST_2020212_16 ~ rh95+air_temp+SM_surf_anom+elevation+aspect, data= sub_dt, family = Gamma(link = "log"), na.action = na.exclude)

print(summary(m1)); print(summary(m2)); print(summary(m3));
print(summary(m4)); print(summary(m5)); print(summary(m6));

kruskal.test(rh95 ~ Forest_cat, data = sub_dt)
kruskal.test(rh50 ~ Forest_cat, data = sub.df)
# Kruskal-Wallis chi-squared = 101.73 and 178.36, df = 2, p-value < 2.2e-16
# strong correlations between canopy height, rh50 and forest degredation type
# --> can't use as covariates in the same glm

m6 <- glm(LST_2020212_16 ~ rh95+air_temp+SM_surf_anom+elevation+aspect, data= sub_dt, family = Gamma(link = "log"), na.action = na.exclude)

summary <- summary(m6); print(summary(m6))

# Creating a data frame with the desired information 
result <- data.frame(Intercept = summary$coeff[1, 1], SE_Intercept = summary$coeff[1, 2],
                     CanopyHeight_coeff = summary$coeff[2, 1], SE_CanopyHeight_coeff = summary$coeff[2, 2], p_CanopyHeight = summary$coeff[2, 4],
                     air_temp_coeff = summary$coeff[3, 1], SE_air_temp_coeff = summary$coeff[3, 2], p_air_temp = summary$coeff[3, 4],
                     SM_surf_anom_coeff = summary$coeff[4, 1], SE_SM_surf_anom_coeff = summary$coeff[4, 2], p_SM_surf_anom = summary$coeff[4, 4],
                     elevation_coeff = summary$coeff[5, 1], SE_elevation_coeff = summary$coeff[5, 2], p_elevation = summary$coeff[5, 4], #)
                     aspect_coeff = summary$coeff[6, 1], SE_aspect_coeff = summary$coeff[6, 2], p_aspect = summary$coeff[6, 4])


                     #Logged_coeff = summary$coeff[2, 1], SE_Logged_coeff = summary$coeff[2, 2], p_Logged = summary$coeff[2, 4],
                     #Burned_coeff = summary$coeff[3, 1], SE_Burned_coeff = summary$coeff[3, 2], p_Burned = summary$coeff[3, 4],
#                      CanopyHeight_coeff = summary$coeff[4, 1], SE_CanopyHeight_coeff = summary$coeff[4, 2], p_CanopyHeight = summary$coeff[4, 4],
#                      air_temp_coeff = summary$coeff[5, 1], SE_air_temp_coeff = summary$coeff[5, 2], p_air_temp = summary$coeff[5, 4],
#                      SM_surf_anom_coeff = summary$coeff[6, 1], SE_SM_surf_anom_coeff = summary$coeff[6, 2], p_SM_surf_anom = summary$coeff[6, 4],
#                      elevation_coeff = summary$coeff[7, 1], SE_elevation_coeff = summary$coeff[7, 2], p_elevation = summary$coeff[7, 4], #)
#                      aspect_coeff = summary$coeff[8, 1], SE_aspect_coeff = summary$coeff[8, 2], p_aspect = summary$coeff[8, 4])
write.csv(result,"/Users/savannah/Documents/Ch2/FNA_pilot/m6glm_SCALEDresult_FNA_sampDist_numSamp5k.csv")

# Compute vcov() matrix to use in uncertainty quantification
#print(vcov(m6))

# df$CH4.CO2e.BackTr<- exp(df$CH4.CO2e.pred) # Adjust for log link in gamma regression 

#######   5. Simulate m5 glm residuals    ####### 

# Investigate residuals
library(DHARMa); library(lme4)

simout  <-  simulateResiduals(m5,  n=1000)
plot(simout)  # plotSimulatedResiduals()  is deprecated, please switch your code to simply using the plot()
plotQQunif(simout)

'''
The QQ plot here is not against the normal distribution, 
it is against the simulation-based expected distribution of the residuals. 
That is the idea with the simulation: we do not know, theoretically, what is the null distribution of the residuals, 
so we approximate it by simulating from the fitted model. 

This have similarities with bayesian simulation-based inference. A paper giving the theoretical background is https://calhoun.nps.edu/bitstream/handle/10945/24467/NPS-OR-95-009.pdf;sequence=3 
A bayesian paper is http://www.stat.columbia.edu/~gelman/research/published/dogs.pdf The other plot, residuals against expected (fitted) value, is augmented with three red lines which should (if the model is correct) be horizontal. They are based on quantile regression of the residuals, for the three quantiles 0.25, 0.5, 0.75. In our example this plot then could indicate some problems, but then it could only be an effect of the very low sample size.
'''



#################       LST ~ structure scientific model & plot data        ###########
sub.df<- read.csv("/Users/savannah/Documents/Ch2/FNA_pilot/all_vars_FNA_sampDist_numSamp5k.csv")
sub.df$rh50_95<- sub.df$rh50/sub.df$rh98
#sub.df$LST_2018212_16<- NULL
sub.df$LST_2020262_19<- NULL
sub.df$LST_2018243_17<- NULL
sub.df$MeanMaxLST<- NULL
#sub.df$MeanMaxLST<- as.numeric(sub.df$MeanMaxLST)
sub.df<- na.omit(sub.df) # from 5k obs to 4.2k (with LST_2018212_16), 4.9k (MeanMaxLST) and 4.7k (with LST_2020262_19)  
sub.df$Forest_cat <- factor(sub.df$Forest_cat, levels = c("Intact", "Logged", "Burned" ))

sub.df %>%
  group_by(Forest_cat) %>%
  summarise(meanCanopyTemp = mean(LST_2018212_16), sd = sd(LST_2018212_16),
            medCanopyTemp = median(LST_2018212_16),
            p75CanopyTemp = quantile(LST_2018212_16, probs = 0.75, na.rm = TRUE),
            p90CanopyTemp = quantile(LST_2018212_16, probs = 0.9, na.rm = TRUE))

hist(sub.df$LST_2018212_16)

ggplot(sub.df, aes(x = rh98, y = LST_2018212_16, color=factor(Forest_cat))) +
  geom_point(size=0.5) + 
  scale_color_manual(name="Degredation type",values=c("red", "green", "blue")) + 
  labs(x = "Canopy Height [m]", y = "Canopy Temperature [ºC]") 
  #labs(x = "Height of median energy, RH50 [m]", y = "Canopy Temperature [ºC]" ) 
  #labs(x = "Ratio RH50/RH98 [m]", y = "Canopy Temperature [ºC]" ) 

sub.df %>% group_by(Forest_cat) %>% 
  summarise(n=n(), meanrh50_95=mean(rh50_95, na.rm=T), medrh50_95=median(rh50_95, na.rm=T),
            mean_rh50=mean(rh50, na.rm=T), med_rh50=median(rh50, na.rm=T), #p75_LST=as.numeric(quantile(r_LST, na.rm=T)[4]),
            max_rh50=max(rh50, na.rm=T), min_rh50=min(rh50, na.rm=T), sd_rh50=sd(rh50, na.rm=T))


#######   MLE: Estimate params for the scientific model   ##########

library(bbmle); library(rjags); library(coda) 

# Function for the scientific model based on Longo et al 2019.
yhat_m1 = function(c, z, rh98){ 
  mu = exp(c) * rh98^z
  return(mu)
} 

# Log-likelihood functions
LL1 <- function(scale, C, Z){
  mu = yhat_m1(c=C, z=Z, rh98=sub.df$rh98)
  -sum(dgamma(x=sub.df$LST_2020212_16, shape = mu/scale, scale= scale, log = TRUE), na.rm = TRUE)
}

# Calculate initial values
mean_y <- mean(sub.df$LST_2020212_16, na.rm=TRUE)
sigma_sq_y <- var(sub.df$LST_2020212_16, na.rm=TRUE)
scale_init <- sigma_sq_y / mean_y

# Fit models
sol1 = mle2(minuslogl = LL1, 
            start = list(scale=scale_init, C=3.8, Z= -0.15), 
            method="L-BFGS-B",
            lower = c(0.001, -Inf, -Inf),
            upper = c(Inf, Inf, Inf))

# Look at the results
summary(sol1)
#summary(sol2)
#summary(sol3)

# Plot the observations vs. predictions

# fitted model 1
#AGB.est<- coef(sol1)[1]
y = yhat_m1(rh98=sub.df$rh98,c=as.numeric(coef(sol1)[2]), # [1] 3.583817 
            z=as.numeric(coef(sol1)[3])) # [1] -0.03257676

df<-data.frame( cbind(sub.df$LST_2020212_16, y) )
names(df)<- c("LST.obs", "LST.pred")
#df<-sample_n(df, size = 10000)
# Plot the fitted model against the data --> need to use the function coef(model.object)
plot(df$LST.obs, df$LST.pred, cex.lab=1.5, cex=1.5, col="black", pch = 6,main="Model 1",
     xlab= "Observed Canopy Temp. [ºC]", ylab= "Modeled  Canopy Temp. [ºC]")  

#######   MLE: Plot scientific model w data   #########

# Function for the scientific model based on Longo et al 2019
yhat_m1_rh98 = function(c, z, rh98){ 
  mu = exp(c) * rh98^z 
  return(mu)
} 
yhat_m1_agbd = function(c, z, agbd){ 
  mu = c * agbd^z # NOTE: for agbd JAGS model I just used c not exp(c)
  return(mu)
} 

# Create a data frame for the yhat_m1 line
x_seq = seq(min(sub.df$agbd_Mg_ha), max(sub.df$agbd_Mg_ha), length.out = 1000)
yhat_df = data.frame(
  #rh98 = x_seq,
  agbd = x_seq,
  yhat = yhat_m1_agbd(c= 33.5, z= -0.007, agbd = x_seq),
  yhat_intact = yhat_m1_agbd(c= 33.7, z= 0.0119697, agbd = x_seq), #(c= 35.1, z= -0.008, agbd = x_seq),
  yhat_logged = yhat_m1_agbd(c= 36.037, z= -0.011447429, agbd = x_seq), #(c= 35.2, z= -0.006, agbd = x_seq),
  yhat_burned = yhat_m1_agbd(c= 43.23994, z= -0.012345989, agbd = x_seq)) # yhat_m1_rh98(c=  3.583817, z= -0.03257676, rh98 = x_seq)

# Preprocess the data
prepare_interleaved_data <- function(df, chunk_size = 20) {
  df %>%
    # Ensure Forest_cat is a factor
    mutate(Forest_cat = factor(Forest_cat, levels = c("Burned", "Intact", "Logged"))) %>%
    # Group and assign IDs
    group_by(Forest_cat) %>%
    mutate(group_id = ceiling(row_number() / chunk_size)) %>%
    ungroup() %>%
    # Create a new ID for interleaving
    mutate(interleave_id = row_number() + (as.numeric(Forest_cat) - 1) * chunk_size) %>%
    # Arrange by the new interleave ID
    arrange(interleave_id) %>%
    # Drop the helper columns
    select(-group_id, -interleave_id)
}

# Apply the function to create interleaved data
sub.df_interleaved <- prepare_interleaved_data(sub.df, chunk_size = 20)

# Create the plot with interleaved points
ggplot(sub.df, aes(x = agbd_Mg_ha, y = LST_2018212_16, color = Forest_cat)) + # x= rh98,
  geom_point(size = 0.8, alpha = 0.6) +
  #geom_line(data = yhat_df, aes(x = agbd, y = yhat_intact), color = "darkgreen", linewidth = 1.3) +
  #geom_line(data = yhat_df, aes(x = agbd, y = yhat_logged), color = "darkblue", linewidth = 1.3) +
  #geom_line(data = yhat_df, aes(x = agbd, y = yhat_burned), color = "darkred", linewidth = 1.3) +
  scale_color_manual(name = "Degradation type",
                     values = c("red", "green", "blue"),
                     labels = c("Burned", "Intact", "Logged")) +
  labs(x = "Aboveground biomass [Mg/ha]" ,#"Canopy Height [m]", 
       y = "Canopy Temperature [ºC]") +
  theme_minimal() +
  theme(legend.position = "right") +
  coord_cartesian(ylim = c(31, max(sub.df$LST_2018212_16, na.rm = TRUE)))

#################       Overview notes & prep: model Canopy temp ~ C,Z ~ structure + covs    ####### 

'''
We developed a hierarchical Bayesian model predicting canopy temperature as a function of canopy height (rh95 from spaceborne lidar data). Parameters depend on hyperparameters that are influenced by forest degradation type among other covariates. 

v3.2 - same as v3 except using LST_2018_243_17 instead of LST_2018212_16
v3.3 - Going to try this prompt also using CLaude v3 code (no priors, trace included) 

Improve the R script pasted here to address the questions and issues I raise. 

To save on credits: 
1. Only provide the relevant code snippets that changed since the last version;
2. Summarize all the changes made in written form with bullet points in the most clear yet concise way possible.


'''


#######   Packages & sub.df data prep: subset & plot AGB data vs temp    #########
library(rjags)
library(coda)
library(ggplot2)
library(gridExtra)
library(dplyr)
library(sf)

df_all <- read.csv("/Users/savannah/Documents/Ch2/FNA_pilot/all_vars_FNA_v2.csv")
sub.df <- df_all 

# Remove unnecessary variables
sub.df$LST_2018212_16 <- NULL 
#sub.df$LST_2018_243_17 <- NULL # This is the warmest, least cloud contaminated day in the 2-yr dry season ECOSTRESS record 
sub.df$LST_2020_262_19 <- NULL
sub.df$MeanMaxLST<- NULL

# Remove rows with probability of bare ground > 4% (based on Dynamic World 2018 land cover classification using Sentinel2 imagery) 
sub.df <- subset(sub.df, probBareGr<0.04) # Additional 2.1% reduction in n: from 6912 to 6767 obs
#hist(sub.df$agbd_Mg_ha_2018_2023); hist(sub.df$rh95); hist(sub.df$year)

# Create new column Forest_cat_v2 using reference data set (column called "class": CVL, BRN )
sub.df <- sub.df %>% #group_by(Biome) %>%
  mutate(Forest_cat_v2 =
           case_when(
             class=="BRN" ~ "Burned",
             class=="CVL" ~ "Logged",
             is.na(class) ~ Forest_cat, # Assumed to be classification prediction
             TRUE ~ Forest_cat # If none of the prev conditions apply, assign Forest_cat classification (Mostly intact forest pixels)
           )) 

# First, ensure time since disturbance is properly prepared for all observations
# For intact forests with NA in year, set a reasonable value (e.g., 50 years)
sub.df$year_adj <- sub.df$year
sub.df$year_adj[is.na(sub.df$year_adj) & sub.df$Forest_cat_v2 == "Intact"] <- 1968 # 50 years before 2018
sub.df$year_adj[is.na(sub.df$year_adj) & sub.df$Forest_cat_v2 != "Intact"] <- 2012 # Mixed forest types outside reference data set - assume classification can only detect degradation <5 years from classification date
sub.df$timeSinceDist <- 2018 - sub.df$year_adj

# Remove all rows with missing values in LST_2018_243_17 column (rows w NAs from other columns must be kept)
sub.df <- sub.df[!is.na(sub.df$LST_2018_243_17), ] ## From ~23k to ~20k obs. 
# Columns that have some NAs: [class, year, Area_ha] & agbd_Mg_ha_2018_2020, agbd_Mg_ha_2018_2022
sub.df <- sub.df[!is.na(sub.df$rh95), ] ## Does not change n obs 
sub.df <- sub.df[!is.na(sub.df$agbd_Mg_ha_2018_2023), ] ## From ~20k to ~18k obs. 

# Define the thresholds based on Forest_cat_v2 
sub.df <- sub.df %>% # From 8k to 6912 obs
  filter(  # Based on Pinage et al 2023
    # Lower bounds - remove observations with potential land use change since 17-June-2018 classification date
    (Forest_cat_v2 == "Intact" & agbd_Mg_ha_2018_2023 > 60) | 
    (Forest_cat_v2 == "Logged" & agbd_Mg_ha_2018_2023 > 40) | 
    (Forest_cat_v2 == "Burned" & agbd_Mg_ha_2018_2023 > 5) )

sub.df %>%
  group_by(Forest_cat) %>%
  summarise(
    heigh_mean = mean(rh95, na.rm = T),
    heigh_median = median(rh95, na.rm = T),
    height_sd = sd(rh95, na.rm = T),
    rh50_median = median(rh50, na.rm = T),
    rh50_mean = mean(rh50, na.rm = T),
    rh50_sd = sd(rh50, na.rm = T),
    LST_median = median(LST_2018_243_17, na.rm = T),
    LST_mean = mean(LST_2018_243_17, na.rm = T),
    LST_sd = sd(LST_2018_243_17, na.rm = T),
    LST_p75 = quantile(LST_2018_243_17, probs = c(0.75), na.rm = T),
    LST_p90 = quantile(LST_2018_243_17, probs = c(0.9), na.rm = T),
    n = n()) 
ggplot(sub.df, aes(x = rh95, y = LST_2018_243_17, color=factor(Forest_cat))) +
  geom_point(size=0.1) + 
  scale_color_manual(name="Forest type",values=c("red", "green", "blue")) + 
  labs(x = "Canopy Height [m]", y = "Canopy Temperature [ºC]") 
#labs(x = "Height of median energy, RH50 [m]", y = "Canopy Temperature [ºC]" ) 
#labs(x = "Ratio RH50/RH98 [m]", y = "Canopy Temperature [ºC]" ) 


#######      Prep data      #######

# Add burned and logged dummy variables
sub.df$Forest_cat_v2 <- factor(sub.df$Forest_cat_v2, levels = c("Intact", "Logged", "Burned"))
sub.df$burned <- as.numeric(sub.df$Forest_cat_v2 == "Burned")
sub.df$logged <- as.numeric(sub.df$Forest_cat_v2 == "Logged")

#   Doughty Fig 1d analog - Canopy temperature distribution
ggplot(sub.df, aes(x=LST_2018_243_17, fill=factor(Forest_cat_v2))) +
  geom_histogram(alpha = 0.7, position="dodge", colour=NA) + 
  scale_fill_manual(values=c("green","blue","red"))+
  theme(legend.title =NA) + theme_minimal() +
  labs(x = "Canopy Temperature (°C)", y = "Count (n pixels)" )

# Standardize continuous variables
sub.df$air_temp_std <- as.numeric(scale(sub.df$air_temp))
sub.df$SM_surf_anom_std <- as.numeric(scale(sub.df$SM_surf_an))
sub.df$elevation_std <- as.numeric(scale(sub.df$elevation))
sub.df$aspect_std <- as.numeric(scale(sub.df$aspect))
sub.df$probBG_std <- as.numeric(scale(sub.df$probBareGr))
sub.df$timeSinceDist_std <- as.numeric(scale(sub.df$timeSinceDist))

# Create spatial quadrats
n_quadrats <- 10
sub.df$quadrat_x <- cut(sub.df$x, breaks = n_quadrats, labels = FALSE)
sub.df$quadrat_y <- cut(sub.df$y, breaks = n_quadrats, labels = FALSE)
sub.df$quadrat <- paste(sub.df$quadrat_x, sub.df$quadrat_y, sep = "_")

# Compute quadrat-level mean covariates
quad_levels <- levels(factor(sub.df$quadrat))
quad_means <- sub.df %>%
  group_by(quadrat) %>%
  summarise(
    air_temp_mean = mean(air_temp_std),
    SM_surf_anom_mean = mean(SM_surf_anom_std),
    elevation_mean = mean(elevation_std),
    aspect_mean = mean(aspect_std),
    probBG_mean = mean(probBG_std),
    #timeSinceDist_mean = mean(timeSinceDist_std),
    burned_mean = mean(burned),
    logged_mean = mean(logged)
  ) %>%
  ungroup() %>%
  mutate(quadrat = as.character(quadrat)) %>%
  arrange(factor(quadrat, levels = quad_levels))

# Prepare data for JAGS
jags_data <- list(
  N = nrow(sub.df),
  y = sub.df$LST_2018_243_17,
  rh95 = sub.df$rh95,
  quadrat = as.numeric(factor(sub.df$quadrat)),
  N_quadrats = length(unique(sub.df$quadrat)),
  air_temp_q = quad_means$air_temp_mean,
  SM_surf_anom_q = quad_means$SM_surf_anom_mean,
  elevation_q = quad_means$elevation_mean,
  aspect_q = quad_means$aspect_mean,
  probBG_q = quad_means$probBG_mean,
  # ADD in next model run: timeSinceDist_q = quad_means$timeSinceDist_mean,
  burned_q = quad_means$burned_mean,
  logged_q = quad_means$logged_mean
)

#######      Define JAGS Model      #######
# JAGS Model
jags_model <- "model {
  # Likelihood
  for (i in 1:N) {
    y[i] ~ dnorm(mu[i], tau)
    mu[i] <- c_adj[quadrat[i]] * pow(rh95[i], z_adj[quadrat[i]])
  }
  
  # Priors for c and z
  for (q in 1:N_quadrats) {
    c[q] ~ dnorm(c_mean, c_tau)
    z[q] ~ dnorm(z_mean, z_tau)
    
    # Adjusted c and z based on quadrat-level covariates
    c_adj[q] <- c[q] + beta_c_burned * burned_q[q] + beta_c_logged * logged_q[q] +
                beta_c_air_temp * air_temp_q[q] + beta_c_SM_surf_anom * SM_surf_anom_q[q] +
                beta_c_elevation * elevation_q[q] + beta_c_aspect * aspect_q[q] + 
                beta_c_ground * probBG_q[q] #+ beta_c_timeSinceDist*timeSinceDist_q[q]
    
    z_adj[q] <- z[q] + beta_z_burned * burned_q[q] + beta_z_logged * logged_q[q] +
                beta_z_air_temp * air_temp_q[q] + beta_z_SM_surf_anom * SM_surf_anom_q[q] +
                beta_z_elevation * elevation_q[q] + beta_z_aspect * aspect_q[q] +
                beta_z_ground * probBG_q[q] #+ beta_z_timeSinceDist*timeSinceDist_q[q]
  }
  
  # Hyperpriors for c and z
  c_mean ~ dnorm(0, 0.001)
  c_tau ~ dgamma(0.001, 0.001)
  z_mean ~ dnorm(0, 0.001)
  z_tau ~ dgamma(0.001, 0.001)
  
  # Priors for beta coefficients
  beta_c_air_temp ~ dnorm(0, 0.001)
  beta_c_SM_surf_anom ~ dnorm(0, 0.001)
  beta_c_elevation ~ dnorm(0, 0.001)
  beta_c_aspect ~ dnorm(0, 0.001)
  beta_c_burned ~ dnorm(0, 0.001)
  beta_c_logged ~ dnorm(0, 0.001)
  beta_c_ground ~ dnorm(0, 0.001)
  
  beta_z_air_temp ~ dnorm(0, 0.001)
  beta_z_SM_surf_anom ~ dnorm(0, 0.001)
  beta_z_elevation ~ dnorm(0, 0.001)
  beta_z_aspect ~ dnorm(0, 0.001)
  beta_z_burned ~ dnorm(0, 0.001)
  beta_z_logged ~ dnorm(0, 0.001)
  beta_z_ground ~ dnorm(0, 0.001)
  
  # Prior for precision
  tau ~ dgamma(0.001, 0.001)
  
  # Derived quantities
  sigma <- 1 / sqrt(tau)
  
  # Calculate mean c and z for each forest type
  c_intact <- c_mean
  z_intact <- z_mean
  c_logged <- c_mean + beta_c_logged
  z_logged <- z_mean + beta_z_logged
  c_burned <- c_mean + beta_c_burned
  z_burned <- z_mean + beta_z_burned
  
  # Predict temperatures for a range of canopy height values
  for (i in 1:101) {
    height_pred[i] <- 1 + (i - 1) * 0.5
    pred_temp_intact[i] <- c_intact * pow(height_pred[i], z_intact)
    pred_temp_logged[i] <- c_logged * pow(height_pred[i], z_logged)
    pred_temp_burned[i] <- c_burned * pow(height_pred[i], z_burned)
  }
}"
#######      Run JAGS Model      #######

# Initialize and run the model
jags_init <- function() {
  list(
    c = rnorm(jags_data$N_quadrats, 0, 0.1),
    z = rnorm(jags_data$N_quadrats, 0, 0.1),
    tau = rgamma(1, 1, 1)
  )
}

# Run JAGS model with updated iterations
jags_fit <- jags.model(textConnection(jags_model), 
                       data = jags_data, 
                       inits = jags_init, 
                       n.chains = 3)

update(jags_fit, 15000)  # Updated burn-in

mcmc_samples <- coda.samples(jags_fit,
                             variable.names = c("pred_temp_intact", "pred_temp_logged", "pred_temp_burned",
                                                "height_pred", "c_intact", "z_intact", "c_logged", "z_logged",
                                                "c_burned", "z_burned", "beta_c_ground", "beta_z_ground"),
                             n.iter = 40000)  # Updated iterations

# Process MCMC samples and create plots (rest of the code remains the same)
calc_mean_ci <- function(x) {
  mean_val <- mean(x)
  ci <- quantile(x, c(0.05, 0.95))
  return(c(mean = mean_val, lower = ci[1], upper = ci[2]))
}

mcmc_combined <- do.call(rbind, mcmc_samples)
height_pred <- mcmc_combined[1, grep("height_pred", colnames(mcmc_combined))]

# Process predictions
process_predictions <- function(mcmc_data, prefix) {
  t(apply(mcmc_data[, grep(prefix, colnames(mcmc_data))], 2, calc_mean_ci))
}

pred_intact <- process_predictions(mcmc_combined, "pred_temp_intact")
pred_logged <- process_predictions(mcmc_combined, "pred_temp_logged")
pred_burned <- process_predictions(mcmc_combined, "pred_temp_burned")

# Prepare plot data
plot_data <- data.frame(
  Height = rep(height_pred, 3),
  Temperature = c(pred_intact[,"mean"], pred_logged[,"mean"], pred_burned[,"mean"]),
  Lower = c(pred_intact[,"lower.5%"], pred_logged[,"lower.5%"], pred_burned[,"lower.5%"]),
  Upper = c(pred_intact[,"upper.95%"], pred_logged[,"upper.95%"], pred_burned[,"upper.95%"]),
  Forest_type = factor(rep(c("Intact", "Logged", "Burned"), each = length(height_pred)))
)

# Create and print plots
fig3 <- ggplot(plot_data, aes(x = Height, y = Temperature, color = Forest_type, fill = Forest_type)) +
  geom_line(linewidth=1.1) +
  geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha = 0.2, linetype = 0) +
  scale_color_manual(values = c("Intact" = "green", "Logged" = "blue", "Burned" = "red")) +
  scale_fill_manual(values = c("Intact" = "green", "Logged" = "blue", "Burned" = "red")) +
  labs(x = "Canopy Height (m)", y = "Predicted Canopy Temperature (°C)") +
  theme_minimal() +
  theme(legend.position = "bottom")

print(fig3)

# Print effect estimates for all covariates
# (currently the code below only does this for bare ground effect)
bg_effects <- data.frame(
  Parameter = c("beta_c_ground", "beta_z_ground"),
  Mean = c(mean(mcmc_combined[, "beta_c_ground"]), 
           mean(mcmc_combined[, "beta_z_ground"])),
  Lower = c(quantile(mcmc_combined[, "beta_c_ground"], 0.05),
            quantile(mcmc_combined[, "beta_z_ground"], 0.05)),
  Upper = c(quantile(mcmc_combined[, "beta_c_ground"], 0.95),
            quantile(mcmc_combined[, "beta_z_ground"], 0.95))
)

print("Bare Ground Effects:")
print(bg_effects)
#######      DIC & diagnostic plots      #######

# Compute DIC and create diagnostic plots
full_run_dic <- dic.samples(jags_fit, n.iter = 1000)
dic_results <- data.frame(
  DIC = sum(full_run_dic$deviance) + sum(full_run_dic$penalty),
  pD = sum(full_run_dic$penalty)
)
print(dic_results)

# Calculate fitted values and residuals
fitted_values <- numeric(jags_data$N)
mcmc_matrix <- as.matrix(mcmc_combined)

# Get mean parameter values
c_intact_mean <- mean(mcmc_matrix[, "c_intact"])
z_intact_mean <- mean(mcmc_matrix[, "z_intact"])
c_logged_mean <- mean(mcmc_matrix[, "c_logged"])
z_logged_mean <- mean(mcmc_matrix[, "z_logged"])
c_burned_mean <- mean(mcmc_matrix[, "c_burned"])
z_burned_mean <- mean(mcmc_matrix[, "z_burned"])

# Calculate fitted values using mean parameter estimates
for (i in 1:jags_data$N) {
  forest_type <- sub.df$Forest_cat[i]
  if (forest_type == "Intact") {
    fitted_values[i] <- c_intact_mean * (jags_data$rh95[i]^z_intact_mean)
  } else if (forest_type == "Logged") {
    fitted_values[i] <- c_logged_mean * (jags_data$rh95[i]^z_logged_mean)
  } else {  # Burned
    fitted_values[i] <- c_burned_mean * (jags_data$rh95[i]^z_burned_mean)
  }
}

residuals <- jags_data$y - fitted_values

# Create diagnostic plots (residuals, spatial residuals, posterior distributions)
residual_plot <- ggplot(data.frame(fitted = fitted_values,
                                   residuals = residuals,
                                   Forest_type = sub.df$Forest_cat),
                        aes(x = fitted, y = residuals, color = Forest_type)) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_color_manual(values = c("Intact" = "green", "Logged" = "blue", "Burned" = "red")) +
  labs(x = "Fitted Values (°C)", y = "Residuals (°C)") +
  theme_minimal()

print(residual_plot)

# Print summary statistics for model evaluation
cat("\nSummary of Residuals by Forest Type:\n")
print(tapply(residuals, sub.df$Forest_cat, function(x) {
  c(mean = mean(x), sd = sd(x),
    q25 = quantile(x, 0.25),
    q75 = quantile(x, 0.75))
}))

# Spatial residual plot
spatial_residual_plot <- ggplot(data.frame(x = sub.df$x, 
                                           y = sub.df$y, 
                                           residuals = residuals,
                                           Forest_type = sub.df$Forest_cat), 
                                aes(x = x, y = y, color = residuals)) +
  geom_point(aes(shape = Forest_type)) +
  scale_color_gradient2(low = "blue", mid = "white", high = "red", 
                        midpoint = 0, name = "Residual (°C)") +
  labs(x = "Longitude", y = "Latitude", 
       title = "")+ #"Spatial Distribution of Residuals") +
  theme_minimal()

# Posterior distribution plots
# Convert MCMC samples to long format for plotting
params_to_plot <- c("c_intact", "z_intact", "c_logged", "z_logged", 
                    "c_burned", "z_burned")

posterior_data <- data.frame(
  Parameter = rep(params_to_plot, each = nrow(mcmc_matrix)),
  Value = c(
    mcmc_matrix[, "c_intact"],
    mcmc_matrix[, "z_intact"],
    mcmc_matrix[, "c_logged"],
    mcmc_matrix[, "z_logged"],
    mcmc_matrix[, "c_burned"],
    mcmc_matrix[, "z_burned"]
  )
)

# Create posterior distribution plot
posterior_plot <- ggplot(posterior_data, aes(x = Value, color = Parameter)) + # , fill = Parameter
  scale_color_manual(values = c("c_intact" = "green", "c_logged" = "blue", "c_burned" = "red1",
                                "z_intact" = "green3", "z_logged" = "blue4", "z_burned" = "red4")) +
  geom_density(alpha = 0.5, fill="grey88", linewidth=1.2) +
  facet_wrap(~Parameter, scales = "free") +
  theme_minimal() +
  labs(title = "", #"Posterior Distributions of Key Parameters",
       x = "Parameter Value",
       y = "Density")
posterior_plot

# Arrange plots in a grid
grid.arrange(residual_plot, spatial_residual_plot, #posterior_plot, 
             ncol = 2, nrow = 1)







