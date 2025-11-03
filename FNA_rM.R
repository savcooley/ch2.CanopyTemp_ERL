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
sub.df$LST_2018212_16 <- NULL;sub.df$LST_2020_262_19 <- NULL; sub.df$MeanMaxLST<- NULL
#sub.df$LST_2018_243_17 <- NULL # This is the warmest, least cloud contaminated day in the 2-yr dry season ECOSTRESS record 
#hist(sub.df$agbd_Mg_ha_2018_2023); hist(sub.df$rh95); hist(sub.df$year)

# Remove all rows with missing values in LST_2018_243_17 column (rows w NAs from other columns must be kept)
sub.df <- sub.df[!is.na(sub.df$LST_2018_243_17), ] ## From ~23k to 17k #~20k obs. 

# Columns that have some NAs: [class, year, Area_ha] & agbd_Mg_ha_2018_2020, agbd_Mg_ha_2018_2022
sub.df <- sub.df[!is.na(sub.df$rh95), ] ## Does not change n obs 
sub.df <- sub.df[!is.na(sub.df$agbd_Mg_ha_2018_2023), ] #~20k to ~18k obs. ## From 17k to 15k 
sub.df <- sub.df[!is.na(sub.df$agbd_Mg_ha_2018_2020), ] #From ~18k to ~10906 obs. ## From 15k to 8665 

# Create new column Forest_cat_v2 using reference data set (column called "class": CVL, BRN )
sub.df <- sub.df %>% 
  mutate(Forest_cat_v2 = # For some reason this leaves no more intact forest obs
           case_when(
             class=="BRN" ~ "Burned",
             class=="CVL" ~ "Logged",
             Forest_cat=="Intact" ~ "Intact",
             is.na(class) ~ Forest_cat, # Assumed to be classification prediction
             TRUE ~ Forest_cat # Catch all clause - If none of the prev conditions apply (which will not happen), assign Forest_cat classification (Mostly intact forest pixels)
           )) 


# Define the thresholds based on Forest_cat_v2 
sub.df <- sub.df %>% # From 10906 to 9177 # From 8665 to 7526  # From 18,705 to 15,745 obs
  filter(  # Based on Pinage et al 2023
    # Lower bounds - remove observations with potential land use change since 17-June-2018 classification date
    (Forest_cat_v2 == "Intact" & agbd_Mg_ha_2018_2020 > 60) | 
      (Forest_cat_v2 == "Logged" & agbd_Mg_ha_2018_2020 > 40) | 
      (Forest_cat_v2 == "Burned" & agbd_Mg_ha_2018_2020 > 5) )

# Remove rows with probability of bare ground > 4% (based on Dynamic World 2018 land cover classification using Sentinel2 imagery) 
sub.df <- subset(sub.df, probBareGr<0.04) # From 9177 to 8980 #7526 to 7331 #  Additional 1.4% reduction in n: from 15745 to 15519 obs

# Ensure time since disturbance is properly prepared for all observations
sub.df$timeSinceDist <- 2018 - sub.df$year # intact forests would be assigned 2018
sub.df$timeSinceDist[sub.df$Forest_cat_v2 == "Intact"] <- 100 # Placeholder for Intact forests - not actually used in analysis
sub.df$timeSinceDist[is.na(sub.df$timeSinceDist)] <- 4 # Burned and logged Forest types outside reference data set - assume classification can only detect degradation <5 years from classification date
sub.df <- subset(sub.df,timeSinceDist > -1) # Final n = 8916. Removed ~20 obs with disturbance after classification date

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

ggplot(sub.df, aes(x = rh95, y = LST_2018_243_17, color=factor(Forest_cat))) +
  geom_point(size=0.2) + theme_minimal() +
  scale_color_manual(name="Forest type",values=c("red", "green", "blue")) + 
  labs(x = "Canopy Height [m]", y = "Canopy Temperature [ºC]") 
#labs(x = "Height of median energy, RH50 [m]", y = "Canopy Temperature [ºC]" ) 
#labs(x = "Ratio RH50/RH98 [m]", y = "Canopy Temperature [ºC]" ) 

# Standardize continuous variables
sub.df$air_temp_std <- as.numeric(scale(sub.df$air_temp))
sub.df$SM_surf_anom_std <- as.numeric(scale(sub.df$SM_surf_an))
sub.df$elevation_std <- as.numeric(scale(sub.df$elevation))
sub.df$aspect_std <- as.numeric(scale(sub.df$aspect))
sub.df$probBG_std <- as.numeric(scale(sub.df$probBareGr))
sub.df$timeSinceDist_std <- as.numeric(scale(sub.df$timeSinceDist))

#######      Models 1-5 specifications in JAGS      #######
# Modified code to compare different model specifications

# Create different model specifications
# Model 1: Original model + time since disturbance
jags_model1 <- "model {
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
              beta_c_ground * probBG_q[q] + beta_c_timeSinceDist * timeSinceDist_q[q]
    
    z_adj[q] <- z[q] + beta_z_burned * burned_q[q] + beta_z_logged * logged_q[q] +
              beta_z_air_temp * air_temp_q[q] + beta_z_SM_surf_anom * SM_surf_anom_q[q] +
              beta_z_elevation * elevation_q[q] + beta_z_aspect * aspect_q[q] +
              beta_z_ground * probBG_q[q] + beta_z_timeSinceDist * timeSinceDist_q[q]
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
  beta_c_timeSinceDist ~ dnorm(0, 0.001)
  
  beta_z_air_temp ~ dnorm(0, 0.001)
  beta_z_SM_surf_anom ~ dnorm(0, 0.001)
  beta_z_elevation ~ dnorm(0, 0.001)
  beta_z_aspect ~ dnorm(0, 0.001)
  beta_z_burned ~ dnorm(0, 0.001)
  beta_z_logged ~ dnorm(0, 0.001)
  beta_z_ground ~ dnorm(0, 0.001)
  beta_z_timeSinceDist ~ dnorm(0, 0.001)
  
  # Prior for precision
  tau ~ dgamma(0.001, 0.001)
  
  # Derived quantities
  sigma <- 1 / sqrt(tau)
  
  # Calculate mean c and z for each forest type at median time since disturbance
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

# Model 2: Simplified model - remove aspect as it may be less important
jags_model2 <- "model {
  # Likelihood
  for (i in 1:N) {
    y[i] ~ dnorm(mu[i], tau)
    mu[i] <- c_adj[quadrat[i]] * pow(rh95[i], z_adj[quadrat[i]])
  }
  
  # Priors for c and z
  for (q in 1:N_quadrats) {
    c[q] ~ dnorm(c_mean, c_tau)
    z[q] ~ dnorm(z_mean, z_tau)
    
    # Adjusted c and z based on quadrat-level covariates - without aspect
    c_adj[q] <- c[q] + beta_c_burned * burned_q[q] + beta_c_logged * logged_q[q] +
              beta_c_air_temp * air_temp_q[q] + beta_c_SM_surf_anom * SM_surf_anom_q[q] +
              beta_c_elevation * elevation_q[q] + 
              beta_c_ground * probBG_q[q] + beta_c_timeSinceDist * timeSinceDist_q[q]
    
    z_adj[q] <- z[q] + beta_z_burned * burned_q[q] + beta_z_logged * logged_q[q] +
              beta_z_air_temp * air_temp_q[q] + beta_z_SM_surf_anom * SM_surf_anom_q[q] +
              beta_z_elevation * elevation_q[q] + 
              beta_z_ground * probBG_q[q] + beta_z_timeSinceDist * timeSinceDist_q[q]
  }
  
  # Same hyperpriors and priors as model 1, without aspect parameters
  c_mean ~ dnorm(0, 0.001)
  c_tau ~ dgamma(0.001, 0.001)
  z_mean ~ dnorm(0, 0.001)
  z_tau ~ dgamma(0.001, 0.001)
  
  beta_c_air_temp ~ dnorm(0, 0.001)
  beta_c_SM_surf_anom ~ dnorm(0, 0.001)
  beta_c_elevation ~ dnorm(0, 0.001)
  beta_c_burned ~ dnorm(0, 0.001)
  beta_c_logged ~ dnorm(0, 0.001)
  beta_c_ground ~ dnorm(0, 0.001)
  beta_c_timeSinceDist ~ dnorm(0, 0.001)
  
  beta_z_air_temp ~ dnorm(0, 0.001)
  beta_z_SM_surf_anom ~ dnorm(0, 0.001)
  beta_z_elevation ~ dnorm(0, 0.001)
  beta_z_burned ~ dnorm(0, 0.001)
  beta_z_logged ~ dnorm(0, 0.001)
  beta_z_ground ~ dnorm(0, 0.001)
  beta_z_timeSinceDist ~ dnorm(0, 0.001)
  
  tau ~ dgamma(0.001, 0.001)
  sigma <- 1 / sqrt(tau)
  
  # Same derived quantities as model 1
  c_intact <- c_mean
  z_intact <- z_mean
  c_logged <- c_mean + beta_c_logged
  z_logged <- z_mean + beta_z_logged
  c_burned <- c_mean + beta_c_burned
  z_burned <- z_mean + beta_z_burned
  
  for (i in 1:101) {
    height_pred[i] <- 1 + (i - 1) * 0.5
    pred_temp_intact[i] <- c_intact * pow(height_pred[i], z_intact)
    pred_temp_logged[i] <- c_logged * pow(height_pred[i], z_logged)
    pred_temp_burned[i] <- c_burned * pow(height_pred[i], z_burned)
  }
}"

# Model 3: Focused model with disturbance and key environmental variables
jags_model3 <- "model {
  # Likelihood
  for (i in 1:N) {
    y[i] ~ dnorm(mu[i], tau)
    mu[i] <- c_adj[quadrat[i]] * pow(rh95[i], z_adj[quadrat[i]])
  }
  
  # Priors for c and z
  for (q in 1:N_quadrats) {
    c[q] ~ dnorm(c_mean, c_tau)
    z[q] ~ dnorm(z_mean, z_tau)
    
    # Focused set of covariates: disturbance, time since disturbance, air temp, and elevation
    c_adj[q] <- c[q] + beta_c_burned * burned_q[q] + beta_c_logged * logged_q[q] +
              beta_c_air_temp * air_temp_q[q] + beta_c_elevation * elevation_q[q] +
              beta_c_timeSinceDist * timeSinceDist_q[q]
    
    z_adj[q] <- z[q] + beta_z_burned * burned_q[q] + beta_z_logged * logged_q[q] +
              beta_z_air_temp * air_temp_q[q] + beta_z_elevation * elevation_q[q] +
              beta_z_timeSinceDist * timeSinceDist_q[q]
  }
  
  # Hyperpriors
  c_mean ~ dnorm(0, 0.001)
  c_tau ~ dgamma(0.001, 0.001)
  z_mean ~ dnorm(0, 0.001)
  z_tau ~ dgamma(0.001, 0.001)
  
  # Reduced set of beta priors
  beta_c_air_temp ~ dnorm(0, 0.001)
  beta_c_elevation ~ dnorm(0, 0.001)
  beta_c_burned ~ dnorm(0, 0.001)
  beta_c_logged ~ dnorm(0, 0.001)
  beta_c_timeSinceDist ~ dnorm(0, 0.001)
  
  beta_z_air_temp ~ dnorm(0, 0.001)
  beta_z_elevation ~ dnorm(0, 0.001)
  beta_z_burned ~ dnorm(0, 0.001)
  beta_z_logged ~ dnorm(0, 0.001)
  beta_z_timeSinceDist ~ dnorm(0, 0.001)
  
  tau ~ dgamma(0.001, 0.001)
  sigma <- 1 / sqrt(tau)
  
  # Derived quantities
  c_intact <- c_mean
  z_intact <- z_mean
  c_logged <- c_mean + beta_c_logged
  z_logged <- z_mean + beta_z_logged
  c_burned <- c_mean + beta_c_burned
  z_burned <- z_mean + beta_z_burned
  
  for (i in 1:101) {
    height_pred[i] <- 1 + (i - 1) * 0.5
    pred_temp_intact[i] <- c_intact * pow(height_pred[i], z_intact)
    pred_temp_logged[i] <- c_logged * pow(height_pred[i], z_logged)
    pred_temp_burned[i] <- c_burned * pow(height_pred[i], z_burned)
  }
}"

# Model 4: Interaction model - with interactions between disturbance and time since disturbance
jags_model4 <- "model {
  # Likelihood
  for (i in 1:N) {
    y[i] ~ dnorm(mu[i], tau)
    mu[i] <- c_adj[quadrat[i]] * pow(rh95[i], z_adj[quadrat[i]])
  }
  
  # Priors for c and z
  for (q in 1:N_quadrats) {
    c[q] ~ dnorm(c_mean, c_tau)
    z[q] ~ dnorm(z_mean, z_tau)
    
    # Add interactions between disturbance type and time since disturbance
    c_adj[q] <- c[q] + beta_c_burned * burned_q[q] + beta_c_logged * logged_q[q] +
              beta_c_air_temp * air_temp_q[q] + beta_c_SM_surf_anom * SM_surf_anom_q[q] +
              beta_c_elevation * elevation_q[q] + beta_c_aspect * aspect_q[q] +
              beta_c_ground * probBG_q[q] + beta_c_timeSinceDist * timeSinceDist_q[q] +
              beta_c_burned_time * (burned_q[q] * timeSinceDist_q[q]) +
              beta_c_logged_time * (logged_q[q] * timeSinceDist_q[q])
    
    z_adj[q] <- z[q] + beta_z_burned * burned_q[q] + beta_z_logged * logged_q[q] +
              beta_z_air_temp * air_temp_q[q] + beta_z_SM_surf_anom * SM_surf_anom_q[q] +
              beta_z_elevation * elevation_q[q] + beta_z_aspect * aspect_q[q] +
              beta_z_ground * probBG_q[q] + beta_z_timeSinceDist * timeSinceDist_q[q] +
              beta_z_burned_time * (burned_q[q] * timeSinceDist_q[q]) +
              beta_z_logged_time * (logged_q[q] * timeSinceDist_q[q])
  }
  
  # Hyperpriors
  c_mean ~ dnorm(0, 0.001)
  c_tau ~ dgamma(0.001, 0.001)
  z_mean ~ dnorm(0, 0.001)
  z_tau ~ dgamma(0.001, 0.001)
  
  # Priors for all beta coefficients, including interaction terms
  beta_c_air_temp ~ dnorm(0, 0.001)
  beta_c_SM_surf_anom ~ dnorm(0, 0.001)
  beta_c_elevation ~ dnorm(0, 0.001)
  beta_c_aspect ~ dnorm(0, 0.001)
  beta_c_burned ~ dnorm(0, 0.001)
  beta_c_logged ~ dnorm(0, 0.001)
  beta_c_ground ~ dnorm(0, 0.001)
  beta_c_timeSinceDist ~ dnorm(0, 0.001)
  beta_c_burned_time ~ dnorm(0, 0.001)
  beta_c_logged_time ~ dnorm(0, 0.001)
  
  beta_z_air_temp ~ dnorm(0, 0.001)
  beta_z_SM_surf_anom ~ dnorm(0, 0.001)
  beta_z_elevation ~ dnorm(0, 0.001)
  beta_z_aspect ~ dnorm(0, 0.001)
  beta_z_burned ~ dnorm(0, 0.001)
  beta_z_logged ~ dnorm(0, 0.001)
  beta_z_ground ~ dnorm(0, 0.001)
  beta_z_timeSinceDist ~ dnorm(0, 0.001)
  beta_z_burned_time ~ dnorm(0, 0.001)
  beta_z_logged_time ~ dnorm(0, 0.001)
  
  tau ~ dgamma(0.001, 0.001)
  sigma <- 1 / sqrt(tau)
  
  # Derived quantities
  c_intact <- c_mean
  z_intact <- z_mean
  c_logged <- c_mean + beta_c_logged
  z_logged <- z_mean + beta_z_logged
  c_burned <- c_mean + beta_c_burned
  z_burned <- z_mean + beta_z_burned
  
  # Include predictions at different times since disturbance
  for (i in 1:101) {
    height_pred[i] <- 1 + (i - 1) * 0.5
    pred_temp_intact[i] <- c_intact * pow(height_pred[i], z_intact)
    pred_temp_logged[i] <- c_logged * pow(height_pred[i], z_logged)
    pred_temp_burned[i] <- c_burned * pow(height_pred[i], z_burned)
  }
}"


# Model 5a: Simplified interactions - interactions only in z (exponent)
jags_model5a <- "model {
  # Likelihood
  for (i in 1:N) {
    y[i] ~ dnorm(mu[i], tau)
    mu[i] <- c_adj[quadrat[i]] * pow(rh95[i], z_adj[quadrat[i]])
  }
  
  # Priors for c and z
  for (q in 1:N_quadrats) {
    c[q] ~ dnorm(c_mean, c_tau)
    z[q] ~ dnorm(z_mean, z_tau)
    
    # c without interaction terms
    c_adj[q] <- c[q] + beta_c_burned * burned_q[q] + beta_c_logged * logged_q[q] +
              beta_c_air_temp * air_temp_q[q] + beta_c_elevation * elevation_q[q] 
    
    # z with interaction terms
    z_adj[q] <- z[q] + beta_z_burned * burned_q[q] + beta_z_logged * logged_q[q] +
              beta_z_air_temp * air_temp_q[q] + beta_z_elevation * elevation_q[q] +
              beta_z_burned_time * (burned_q[q] * timeSinceDist_q[q]) +
              beta_z_logged_time * (logged_q[q] * timeSinceDist_q[q])
  }
  
  # Hyperpriors
  c_mean ~ dnorm(0, 0.001)
  c_tau ~ dgamma(0.001, 0.001)
  z_mean ~ dnorm(0, 0.001)
  z_tau ~ dgamma(0.001, 0.001)
  
  # Reduced set of beta priors
  beta_c_air_temp ~ dnorm(0, 0.001)
  beta_c_elevation ~ dnorm(0, 0.001)
  beta_c_burned ~ dnorm(0, 0.001)
  beta_c_logged ~ dnorm(0, 0.001)
  beta_c_timeSinceDist ~ dnorm(0, 0.001)
  
  beta_z_air_temp ~ dnorm(0, 0.001)
  beta_z_elevation ~ dnorm(0, 0.001)
  beta_z_burned ~ dnorm(0, 0.001)
  beta_z_logged ~ dnorm(0, 0.001)
  beta_z_timeSinceDist ~ dnorm(0, 0.001)
  beta_z_burned_time ~ dnorm(0, 0.001)
  beta_z_logged_time ~ dnorm(0, 0.001)
  
  tau ~ dgamma(0.001, 0.001)
  sigma <- 1 / sqrt(tau)
  
  # Derived quantities
  c_intact <- c_mean
  z_intact <- z_mean
  c_logged <- c_mean + beta_c_logged
  z_logged <- z_mean + beta_z_logged
  c_burned <- c_mean + beta_c_burned
  z_burned <- z_mean + beta_z_burned
  
  for (i in 1:101) {
    height_pred[i] <- 1 + (i - 1) * 0.5
    pred_temp_intact[i] <- c_intact * pow(height_pred[i], z_intact)
    pred_temp_logged[i] <- c_logged * pow(height_pred[i], z_logged)
    pred_temp_burned[i] <- c_burned * pow(height_pred[i], z_burned)
  }
}"

# Model 5b: Interactions only in c (coefficient)
jags_model5b <- "model {
  # Likelihood
  for (i in 1:N) {
    y[i] ~ dnorm(mu[i], tau)
    mu[i] <- c_adj[quadrat[i]] * pow(rh95[i], z_adj[quadrat[i]])
  }
  
  # Priors for c and z
  for (q in 1:N_quadrats) {
    c[q] ~ dnorm(c_mean, c_tau)
    z[q] ~ dnorm(z_mean, z_tau)
    
    # c with interaction terms
    c_adj[q] <- c[q] + beta_c_burned * burned_q[q] + beta_c_logged * logged_q[q] +
              beta_c_air_temp * air_temp_q[q] + beta_c_elevation * elevation_q[q] +
              beta_c_burned_time * (burned_q[q] * timeSinceDist_q[q]) +
              beta_c_logged_time * (logged_q[q] * timeSinceDist_q[q])
    
    # z without interaction terms
    z_adj[q] <- z[q] + beta_z_burned * burned_q[q] + beta_z_logged * logged_q[q] +
              beta_z_air_temp * air_temp_q[q] + beta_z_elevation * elevation_q[q] 
  }
  
  # Hyperpriors
  c_mean ~ dnorm(0, 0.001)
  c_tau ~ dgamma(0.001, 0.001)
  z_mean ~ dnorm(0, 0.001)
  z_tau ~ dgamma(0.001, 0.001)
  
  # Reduced set of beta priors
  beta_c_air_temp ~ dnorm(0, 0.001)
  beta_c_elevation ~ dnorm(0, 0.001)
  beta_c_burned ~ dnorm(0, 0.001)
  beta_c_logged ~ dnorm(0, 0.001)
  beta_c_burned_time ~ dnorm(0, 0.001)
  beta_c_logged_time ~ dnorm(0, 0.001)
  
  beta_z_air_temp ~ dnorm(0, 0.001)
  beta_z_elevation ~ dnorm(0, 0.001)
  beta_z_burned ~ dnorm(0, 0.001)
  beta_z_logged ~ dnorm(0, 0.001)
  
  tau ~ dgamma(0.001, 0.001)
  sigma <- 1 / sqrt(tau)
  
  # Derived quantities
  c_intact <- c_mean
  z_intact <- z_mean
  c_logged <- c_mean + beta_c_logged
  z_logged <- z_mean + beta_z_logged
  c_burned <- c_mean + beta_c_burned
  z_burned <- z_mean + beta_z_burned
  
  for (i in 1:101) {
    height_pred[i] <- 1 + (i - 1) * 0.5
    pred_temp_intact[i] <- c_intact * pow(height_pred[i], z_intact)
    pred_temp_logged[i] <- c_logged * pow(height_pred[i], z_logged)
    pred_temp_burned[i] <- c_burned * pow(height_pred[i], z_burned)
  }
}"

# Model 5c Revised: Remove time since disturbance effect for intact forests
jags_model5c <- "model {
  # Likelihood
  for (i in 1:N) {
    y[i] ~ dnorm(mu[i], tau)
    mu[i] <- c_adj[quadrat[i]] * pow(rh95[i], z_adj[quadrat[i]])
  }
  
  # Priors for c and z
  for (q in 1:N_quadrats) {
    c[q] ~ dnorm(c_mean, c_tau)
    z[q] ~ dnorm(z_mean, z_tau)
    
    # Only air temperature for environmental variables in c
    # Remove timeSinceDist main effect - only include interactions with degradation
    c_adj[q] <- c[q] + beta_c_burned * burned_q[q] + beta_c_logged * logged_q[q] +
              beta_c_air_temp * air_temp_q[q] +
              beta_c_burned_time * (burned_q[q] * timeSinceDist_q[q]) +
              beta_c_logged_time * (logged_q[q] * timeSinceDist_q[q])
    
    # Only elevation for environmental variables in z
    # Remove timeSinceDist main effect - only include interactions with degradation
    z_adj[q] <- z[q] + beta_z_burned * burned_q[q] + beta_z_logged * logged_q[q] +
              beta_z_elevation * elevation_q[q] +
              beta_z_burned_time * (burned_q[q] * timeSinceDist_q[q]) +
              beta_z_logged_time * (logged_q[q] * timeSinceDist_q[q])
  }
  
  # Hyperpriors
  c_mean ~ dnorm(0, 0.001)
  c_tau ~ dgamma(0.001, 0.001)
  z_mean ~ dnorm(0, 0.001)
  z_tau ~ dgamma(0.001, 0.001)
  
  # Reduced set of beta priors
  # Remove timeSinceDist parameters
  beta_c_air_temp ~ dnorm(0, 0.001)
  beta_c_burned ~ dnorm(0, 0.001)
  beta_c_logged ~ dnorm(0, 0.001)
  beta_c_burned_time ~ dnorm(0, 0.001)
  beta_c_logged_time ~ dnorm(0, 0.001)
  
  beta_z_elevation ~ dnorm(0, 0.001)
  beta_z_burned ~ dnorm(0, 0.001)
  beta_z_logged ~ dnorm(0, 0.001)
  beta_z_burned_time ~ dnorm(0, 0.001)
  beta_z_logged_time ~ dnorm(0, 0.001)
  
  tau ~ dgamma(0.001, 0.001)
  sigma <- 1 / sqrt(tau)
  
  # Derived quantities
  c_intact <- c_mean
  z_intact <- z_mean
  c_logged <- c_mean + beta_c_logged
  z_logged <- z_mean + beta_z_logged
  c_burned <- c_mean + beta_c_burned
  z_burned <- z_mean + beta_z_burned
  
  for (i in 1:101) {
    height_pred[i] <- 1 + (i - 1) * 0.5
    pred_temp_intact[i] <- c_intact * pow(height_pred[i], z_intact)
    pred_temp_logged[i] <- c_logged * pow(height_pred[i], z_logged)
    pred_temp_burned[i] <- c_burned * pow(height_pred[i], z_burned)
  }
}"


#######      Functions to prepare, run and test models      #######
# Function to prepare quadrat-level data with time since disturbance
prepare_quadrat_data <- function(sub.df, n_quadrats) {
  # Create spatial quadrats
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
      timeSinceDist_mean = mean(timeSinceDist_std),
      burned_mean = mean(burned),
      logged_mean = mean(logged)
    ) %>%
    ungroup() %>%
    mutate(quadrat = as.character(quadrat)) %>%
    arrange(factor(quadrat, levels = quad_levels))
  
  # Prepare JAGS data
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
    timeSinceDist_q = quad_means$timeSinceDist_mean,
    burned_q = quad_means$burned_mean,
    logged_q = quad_means$logged_mean
  )
  
  return(jags_data)
}

# Function to run a model and calculate DIC
run_model_and_get_dic <- function(jags_model_text, jags_data, n_burn = 15000, n_iter = 45000) {
  # Initialize model
  jags_init <- function() {
    list(
      c = rnorm(jags_data$N_quadrats, 0, 0.1),
      z = rnorm(jags_data$N_quadrats, 0, 0.1),
      tau = rgamma(1, 1, 1)
    )
  }
  
  # Run JAGS model
  jags_fit <- jags.model(textConnection(jags_model_text),
                         data = jags_data,
                         inits = jags_init,
                         n.chains = 3)
  
  update(jags_fit, n_burn)  # Burn-in
  
  # Calculate DIC
  dic_result <- dic.samples(jags_fit, n.iter = 1000)
  dic_value <- sum(dic_result$deviance) + sum(dic_result$penalty)
  pd_value <- sum(dic_result$penalty)
  
  # Calculate WAIC (if available in your version of JAGS/rjags)
  # This requires additional packages like 'loo'
  # waic_result <- try(waic(jags_fit), silent = TRUE)
  # waic_value <- if(class(waic_result) != "try-error") waic_result$waic else NA
  
  # Return DIC and pD
  return(list(
    DIC = dic_value,
    pD = pd_value
    # WAIC = waic_value  # Uncomment if you have WAIC calculation capability
  ))
}

# Function to test different numbers of quadrats
test_quadrat_numbers <- function(sub.df, quad_numbers = c(5, 10, 15, 20)) {
  results <- data.frame(
    n_quadrats = integer(),
    DIC = numeric(),
    pD = numeric()
  )
  
  # Use a simpler model for quadrat testing to save computation time
  simple_model <- jags_model3
  
  for (n_quad in quad_numbers) {
    cat("Testing with", n_quad, "quadrats...\n")
    
    # Prepare data with current number of quadrats
    test_data <- prepare_quadrat_data(sub.df, n_quad)
    
    # Run model with shorter iterations for testing
    model_result <- run_model_and_get_dic(simple_model, test_data, n_burn = 5000, n_iter = 10000)
    
    # Store results
    results <- rbind(results, data.frame(
      n_quadrats = n_quad,
      DIC = model_result$DIC,
      pD = model_result$pD
    ))
  }
  
  return(results)
}

# Function to compare all model specifications
compare_models <- function(sub.df, n_quadrats) {
  # Prepare data
  jags_data <- prepare_quadrat_data(sub.df, n_quadrats)
  
  # List of models to compare
  model_list <- list(
    #"Model 1: Full + Time" = jags_model1,
    #"Model 2: No Aspect" = jags_model2,
    "Model 3: Focused" = jags_model3,
    #"Model 4: Interactions" = jags_model4,
    #"Model 5: Focused+Inter" = jags_model5,
    "Model 5a: z-Inter Only" = jags_model5a,
    "Model 5b: c-Inter Only" = jags_model5b,
    "Model 5c: Simplified" = jags_model5c
  )
  
  # Run each model and collect results
  results <- data.frame(
    Model = character(),
    DIC = numeric(),
    pD = numeric(),
    Parameters = integer()
  )
  
  for (model_name in names(model_list)) {
    cat("Running", model_name, "...\n")
    
    # Count number of parameters (approximate)
    param_count <- length(gregexpr("beta_[cz]_", model_list[[model_name]])[[1]])
    if (param_count < 0) param_count <- 0  # Handle case of no matches
    
    # Run model
    model_result <- run_model_and_get_dic(model_list[[model_name]], jags_data)
    
    # Store results
    results <- rbind(results, data.frame(
      Model = model_name,
      DIC = model_result$DIC,
      pD = model_result$pD,
      Parameters = param_count
    ))
  }
  
  # Calculate AIC-like measure (DIC + 2*k)
  results$AIC_like <- results$DIC + 2 * results$Parameters
  
  # Sort by DIC
  results <- results[order(results$DIC), ]
  
  return(results)
}

# Function to extract and process the interaction effects
analyze_interactions <- function(mcmc_samples, model_name) {
  # Convert to matrix for easier manipulation
  mcmc_matrix <- as.matrix(do.call(rbind, mcmc_samples))
  
  # Extract parameters of interest
  params <- c(
    # Add c_mean and z_mean to the list of parameters
    "c_mean", "z_mean", 
    "beta_c_burned", "beta_c_logged", "beta_c_timeSinceDist", 
    "beta_z_burned", "beta_z_logged", "beta_z_timeSinceDist"
  )
  
  # Check if interaction terms exist
  if("beta_c_burned_time" %in% colnames(mcmc_matrix)) {
    params <- c(params, "beta_c_burned_time", "beta_c_logged_time")
  }
  
  if("beta_z_burned_time" %in% colnames(mcmc_matrix)) {
    params <- c(params, "beta_z_burned_time", "beta_z_logged_time")
  }
  
  # Filter parameters that exist in the model
  params <- params[params %in% colnames(mcmc_matrix)]
  
  # Calculate mean and credible intervals
  param_summary <- data.frame(
    Parameter = params,
    Mean = sapply(params, function(p) mean(mcmc_matrix[, p])),
    Q2.5 = sapply(params, function(p) quantile(mcmc_matrix[, p], 0.025)),
    Q97.5 = sapply(params, function(p) quantile(mcmc_matrix[, p], 0.975))
  )
  
  param_summary$Significant <- (param_summary$Q2.5 * param_summary$Q97.5) > 0
  
  # Add model name for comparison
  param_summary$Model <- model_name
  
  return(param_summary)
}

# Faster approach to generate Fig 3 compared to predict_recovery_trajectory
calculate_temperature_differences <- function(mcmc_samples, height_range = seq(5, 40, by=5), recovery_time = 5) {
  # Convert to matrix for easier manipulation
  mcmc_matrix <- as.matrix(do.call(rbind, mcmc_samples))
  
  # Number of posterior samples to use (can reduce for speed if needed)
  n_samples <- min(5000, nrow(mcmc_matrix))
  
  # Get mean and sd for time standardization
  time_mean <- mean(sub.df$timeSinceDist, na.rm=TRUE)
  time_sd <- sd(sub.df$timeSinceDist, na.rm=TRUE)
  
  # Standardize recovery time
  time_std <- (recovery_time - time_mean) / time_sd
  
  # Initialize matrices to store results
  temp_intact <- matrix(NA, nrow = n_samples, ncol = length(height_range))
  temp_logged <- matrix(NA, nrow = n_samples, ncol = length(height_range))
  temp_burned <- matrix(NA, nrow = n_samples, ncol = length(height_range))
  
  # Calculate temperatures for each forest type across all heights and samples
  for (i in 1:n_samples) {
    # Extract parameters for this sample
    c_mean <- mcmc_matrix[i, "c_mean"]
    z_mean <- mcmc_matrix[i, "z_mean"]
    beta_c_burned <- mcmc_matrix[i, "beta_c_burned"]
    beta_c_logged <- mcmc_matrix[i, "beta_c_logged"]
    beta_z_burned <- mcmc_matrix[i, "beta_z_burned"]
    beta_z_logged <- mcmc_matrix[i, "beta_z_logged"]
    
    # Extract interaction effects if they exist
    beta_c_burned_time <- if("beta_c_burned_time" %in% colnames(mcmc_matrix))
      mcmc_matrix[i, "beta_c_burned_time"] else 0
    beta_c_logged_time <- if("beta_c_logged_time" %in% colnames(mcmc_matrix))
      mcmc_matrix[i, "beta_c_logged_time"] else 0
    beta_z_burned_time <- if("beta_z_burned_time" %in% colnames(mcmc_matrix))
      mcmc_matrix[i, "beta_z_burned_time"] else 0
    beta_z_logged_time <- if("beta_z_logged_time" %in% colnames(mcmc_matrix))
      mcmc_matrix[i, "beta_z_logged_time"] else 0
    
    # Calculate c and z values with time effects
    c_intact <- c_mean
    z_intact <- z_mean
    
    c_logged <- c_mean + beta_c_logged + beta_c_logged_time * time_std
    z_logged <- z_mean + beta_z_logged + beta_z_logged_time * time_std
    
    c_burned <- c_mean + beta_c_burned + beta_c_burned_time * time_std
    z_burned <- z_mean + beta_z_burned + beta_z_burned_time * time_std
    
    for (j in 1:length(height_range)) {
      h <- height_range[j]
      
      # Calculate temperatures for each forest type
      temp_intact[i, j] <- c_intact * (h^z_intact)
      temp_logged[i, j] <- c_logged * (h^z_logged)
      temp_burned[i, j] <- c_burned * (h^z_burned)
    }
  }
  
  # Calculate differences between forest types
  diff_burned_intact <- temp_burned - temp_intact
  diff_logged_intact <- temp_logged - temp_intact
  diff_burned_logged <- temp_burned - temp_logged
  
  # Create summary dataframes for each forest type
  summarize_temps <- function(temp_matrix, type) {
    data.frame(
      Forest_Type = type,
      Height = rep(height_range, each = 1),
      Mean_Temp = colMeans(temp_matrix),
      Lower_CI = apply(temp_matrix, 2, function(x) quantile(x, 0.025)),
      Upper_CI = apply(temp_matrix, 2, function(x) quantile(x, 0.975))
    )
  }
  
  # Create summary dataframes for differences
  summarize_diffs <- function(diff_matrix, type) {
    data.frame(
      Comparison = type,
      Height = rep(height_range, each = 1),
      Mean_Diff = colMeans(diff_matrix),
      Lower_CI = apply(diff_matrix, 2, function(x) quantile(x, 0.025)),
      Upper_CI = apply(diff_matrix, 2, function(x) quantile(x, 0.975))
    )
  }
  
  # Generate all summary tables
  temps_intact <- summarize_temps(temp_intact, "Intact")
  temps_logged <- summarize_temps(temp_logged, "Logged")
  temps_burned <- summarize_temps(temp_burned, "Burned")
  
  diff_bi <- summarize_diffs(diff_burned_intact, "Burned vs Intact")
  diff_li <- summarize_diffs(diff_logged_intact, "Logged vs Intact")
  diff_bl <- summarize_diffs(diff_burned_logged, "Burned vs Logged")
  
  # Combine temperature tables
  all_temps <- rbind(temps_intact, temps_logged, temps_burned)
  
  # Combine difference tables
  all_diffs <- rbind(diff_bi, diff_li, diff_bl)
  
  # Calculate overall summaries across all heights
  overall_diffs <- data.frame(
    Comparison = c("Burned vs Intact", "Logged vs Intact", "Burned vs Logged"),
    Avg_Diff = c(mean(colMeans(diff_burned_intact)), mean(colMeans(diff_logged_intact)), mean(colMeans(diff_burned_logged))),
    Min_Diff = c(min(colMeans(diff_burned_intact)), min(colMeans(diff_logged_intact)), min(colMeans(diff_burned_logged))),
    Max_Diff = c(max(colMeans(diff_burned_intact)), max(colMeans(diff_logged_intact)), max(colMeans(diff_burned_logged))),
    Min_CI = c(min(apply(diff_burned_intact, 2, function(x) quantile(x, 0.025))), 
               min(apply(diff_logged_intact, 2, function(x) quantile(x, 0.025))),
               min(apply(diff_burned_logged, 2, function(x) quantile(x, 0.025)))),
    Max_CI = c(max(apply(diff_burned_intact, 2, function(x) quantile(x, 0.975))),
               max(apply(diff_logged_intact, 2, function(x) quantile(x, 0.975))),
               max(apply(diff_burned_logged, 2, function(x) quantile(x, 0.975))))
  )
  
  # Return all results
  return(list(
    Temperatures = all_temps,
    Differences = all_diffs,
    Overall_Summary = overall_diffs,
    RecoveryTime = recovery_time
  ))
}

# Function to predict temperature at different canopy heights and times since disturbance
# Corrected version of predict_recovery_trajectory with proper error handling for intact forests
predict_recovery_trajectory <- function(mcmc_samples, model_name, max_samples = 3000) {
  # Convert to matrix for easier manipulation - use only a small subsample
  mcmc_matrix <- as.matrix(do.call(rbind, mcmc_samples))
  
  # Determine subsample size - use at most max_samples iterations
  n_iterations <- nrow(mcmc_matrix)
  subsample_indices <- seq(1, n_iterations, length.out = min(max_samples, n_iterations))
  subsample_indices <- round(subsample_indices)
  
  # Extract the subsample
  mcmc_matrix <- mcmc_matrix[subsample_indices, ]
  
  # Set up recovery times to evaluate (years since disturbance)
  recovery_times <- c(1, 5, 10, 20, 30)
  
  # Set up heights to evaluate
  heights <- c(5, 10,  15, 20, 25, 30, 35, 40) #7.5, 12.5,17.5,22.5, 27.5, 32.5, 37.5,
  
  # Pre-allocate results matrix for speed
  n_samples <- nrow(mcmc_matrix)
  n_combinations <- length(recovery_times) * length(heights) * 3  # 3 forest types
  results <- data.frame(
    Sample = integer(n_samples * n_combinations),
    Forest = character(n_samples * n_combinations),
    RecoveryTime = numeric(n_samples * n_combinations),
    Height = numeric(n_samples * n_combinations),
    Temperature = numeric(n_samples * n_combinations),
    stringsAsFactors = FALSE
  )
  
  # Use vectorization where possible
  row_index <- 1
  
  # Get mean and sd for time standardization (do once outside the loop)
  time_mean <- mean(sub.df$timeSinceDist, na.rm=TRUE)
  time_sd <- sd(sub.df$timeSinceDist, na.rm=TRUE)
  
  # Process each sample
  for(i in 1:n_samples) {
    # Extract parameters for this iteration
    c_mean <- mcmc_matrix[i, "c_mean"]
    z_mean <- mcmc_matrix[i, "z_mean"]
    beta_c_burned <- mcmc_matrix[i, "beta_c_burned"]
    beta_c_logged <- mcmc_matrix[i, "beta_c_logged"]
    beta_z_burned <- mcmc_matrix[i, "beta_z_burned"]
    beta_z_logged <- mcmc_matrix[i, "beta_z_logged"]
    
    # Extract time effects if they exist (we'll ignore these for intact forests)
    beta_c_time <- if("beta_c_timeSinceDist" %in% colnames(mcmc_matrix)) 
      mcmc_matrix[i, "beta_c_timeSinceDist"] else 0
    beta_z_time <- if("beta_z_timeSinceDist" %in% colnames(mcmc_matrix)) 
      mcmc_matrix[i, "beta_z_timeSinceDist"] else 0
    
    # Extract interaction effects if they exist
    beta_c_burned_time <- if("beta_c_burned_time" %in% colnames(mcmc_matrix)) 
      mcmc_matrix[i, "beta_c_burned_time"] else 0
    beta_c_logged_time <- if("beta_c_logged_time" %in% colnames(mcmc_matrix)) 
      mcmc_matrix[i, "beta_c_logged_time"] else 0
    beta_z_burned_time <- if("beta_z_burned_time" %in% colnames(mcmc_matrix)) 
      mcmc_matrix[i, "beta_z_burned_time"] else 0
    beta_z_logged_time <- if("beta_z_logged_time" %in% colnames(mcmc_matrix)) 
      mcmc_matrix[i, "beta_z_logged_time"] else 0
    
    # Process combinations more efficiently
    for(forest in c("Intact", "Burned", "Logged")) {
      for(time in recovery_times) {
        # Standardize time (once per time value)
        time_std <- (time - time_mean) / time_sd
        
        # Calculate c and z values based on forest type and time
        if(forest == "Intact") {
          # For intact forests, use only the base parameters
          # Ignore time since disturbance effects completely
          c_val <- c_mean
          z_val <- z_mean
        } else if(forest == "Burned") {
          # For burned forests, include time effects and interactions
          c_val <- c_mean + beta_c_burned + beta_c_time * time_std + 
            beta_c_burned_time * time_std
          z_val <- z_mean + beta_z_burned + beta_z_time * time_std + 
            beta_z_burned_time * time_std
        } else { # Logged
          # For logged forests, include time effects and interactions
          c_val <- c_mean + beta_c_logged + beta_c_time * time_std + 
            beta_c_logged_time * time_std
          z_val <- z_mean + beta_z_logged + beta_z_time * time_std + 
            beta_z_logged_time * time_std
        }
        
        # Process all heights at once
        for(height in heights) {
          # Calculate predicted temperature
          temp <- c_val * (height^z_val)
          
          # Add to results dataframe
          results[row_index, ] <- list(i, forest, time, height, temp)
          row_index <- row_index + 1
        }
      }
    }
  }
  
  # Remove any unused rows if we pre-allocated too many
  if(row_index <= nrow(results)) {
    results <- results[1:(row_index-1), ]
  }
  
  # Summarize results - use data.table for faster aggregation if available
  if(requireNamespace("data.table", quietly = TRUE)) {
    library(data.table)
    dt <- as.data.table(results)
    summary_results <- dt[, .(
      MeanTemp = mean(Temperature),
      Q2.5 = quantile(Temperature, 0.025),
      Q97.5 = quantile(Temperature, 0.975)
    ), by = .(Forest, RecoveryTime, Height)]
    summary_results <- as.data.frame(summary_results)
  } else {
    # Fallback to dplyr if data.table not available
    summary_results <- results %>%
      group_by(Forest, RecoveryTime, Height) %>%
      summarise(
        MeanTemp = mean(Temperature),
        Q2.5 = quantile(Temperature, 0.025),
        Q97.5 = quantile(Temperature, 0.975),
        .groups = 'drop'
      )
  }
  
  summary_results$Model <- model_name
  
  return(summary_results)
}

# The plot function remains the same, it will now work with the corrected data
plot_height_temp_relationship <- function(trajectory_data, recovery_time = 10) { # recovery_times could be 1,5,10,20,30 
  # Filter data for the specified recovery time, keeping all forest types
  plot_data <- trajectory_data %>% 
    filter(RecoveryTime == recovery_time)
  
  # Create plot
  ggplot(plot_data, aes(x = Height, y = MeanTemp, color = Forest, fill = Forest)) +
    geom_line(linewidth = 1) +
    geom_ribbon(aes(ymin = Q2.5, ymax = Q97.5), alpha = 0.2, linetype = 0) +
    scale_color_manual(values = c("Intact" = "green", "Logged" = "blue", "Burned" = "red")) +
    scale_fill_manual(values = c("Intact" = "green", "Logged" = "blue", "Burned" = "red")) +
    labs(
      x = "Canopy Height (m)",
      y = "Predicted Temperature (°C)",
      title = paste0("Temperature vs Height Relationship at ", recovery_time, " Years Since Disturbance")
    ) +
    theme_minimal()
}
# Function to plot recovery trajectories
plot_recovery_trajectories <- function(trajectory_data, height_value = 15) {
  # Filter data for the specified height
  plot_data <- trajectory_data %>% 
    filter(Height == height_value)
  
  # Create plot
  ggplot(plot_data, aes(x = RecoveryTime, y = MeanTemp, color = Forest, fill = Forest)) +
    geom_line(linewidth = 1) +
    geom_ribbon(aes(ymin = Q2.5, ymax = Q97.5), alpha = 0.2, linetype = 0) +
    scale_color_manual(values = c("Intact" = "green", "Logged" = "blue", "Burned" = "red")) +
    scale_fill_manual(values = c("Intact" = "green", "Logged" = "blue", "Burned" = "red")) +
    labs(
      x = "Time Since Disturbance (years)",
      y = paste0("Predicted Temperature at ", height_value, "m Height (°C)"),
      title = paste0("Forest Recovery Trajectory at ", height_value, "m Canopy Height")
    ) +
    theme_minimal()
}

#######      Implementation: evaluate models      #######
#samp.df<-sample_n(sub.df, 3000)
# 1. Test different numbers of quadrats
#quadrat_results <- test_quadrat_numbers(samp.df, c(5, 10, 15));print(quadrat_results)
#optimal_quadrats <- quadrat_results$n_quadrats[which.min(quadrat_results$DIC)]

# 2. Compare all model specifications with optimal number of quadrats
# model_comparison <- compare_models(samp.df, 10) #(sub.df, optimal_quadrats)
# print(model_comparison)
# setwd("/Users/savannah/Documents/Ch2/FNA_pilot/ch2_Figs_and_tables/")
# write.csv(model_comparison, "Table_S2_models_1to5_comparison_3ksamp.csv")

# 3. Select best model and run with full iterations
best_model_name <- "Model 5c: Simplified" #"Model 5b: c-Inter Only" ##model_comparison$Model[1]
best_model <- get(paste0("jags_model", substr(best_model_name, 7, 8)))

#######      Implementation: Run best model      #######

# Initialize and run the selected model
jags_data <- prepare_quadrat_data(sub.df, 20) #optimal_quadrats)
jags_init <- function() {
  list(
    c = rnorm(jags_data$N_quadrats, 0, 0.1),
    z = rnorm(jags_data$N_quadrats, 0, 0.1),
    tau = rgamma(1, 1, 1)
  )
}

# Run JAGS model with full iterations
jags_fit <- jags.model(textConnection(best_model),
                       data = jags_data, inits = jags_init,n.chains = 3)

update(jags_fit, 15000)  # Burn-in

# Define parameters based on the specific model
if(best_model_name == "Model 5c: Simplified") {
  variable_names <- c(
    # Basic parameters
    "c_mean", "z_mean", "c_intact", "z_intact", "c_logged", "z_logged",
    "c_burned", "z_burned", "sigma",
    
    # Model 5c specific parameters
    "beta_c_air_temp", "beta_z_elevation", 
    "beta_c_burned", "beta_c_logged", "beta_z_burned", "beta_z_logged",
    "beta_c_burned_time", "beta_c_logged_time", 
    "beta_z_burned_time", "beta_z_logged_time"
  )
} else if(best_model_name =="Model 5b: c-Inter Only") {
  variable_names <- c(
    # Basic parameters
    "c_mean", "z_mean", "c_intact", "z_intact", "c_logged", "z_logged",
    "c_burned", "z_burned", "sigma",
    
    # Model 5b specific parameters
    "beta_c_air_temp", "beta_c_elevation", "beta_z_air_temp","beta_z_elevation", 
    "beta_c_burned", "beta_c_logged", "beta_z_burned", "beta_z_logged",
    "beta_c_burned_time", "beta_c_logged_time"
  )
}else {
  # For other models, start with basic parameters
  variable_names <- c(
    "c_mean", "z_mean", "c_intact", "z_intact", "c_logged", "z_logged",
    "c_burned", "z_burned", "sigma"
  )
  
  # Add all beta parameters from a test run
  test_samples <- coda.samples(jags_fit, variable.names = "beta_c_burned", n.iter = 10)
  beta_pattern <- paste0("beta_[cz]_", collapse="|")
  mcmc_matrix <- as.matrix(do.call(rbind, test_samples))
  beta_cols <- grep(beta_pattern, colnames(mcmc_matrix), value=TRUE)
  variable_names <- c(variable_names, beta_cols)
}

# Sample from posterior with specific parameters
mcmc_samples <- coda.samples(jags_fit, variable.names = variable_names, n.iter = 15000)

# 4. Analyze model outputs
param_summary <- analyze_interactions(mcmc_samples, best_model_name)
print(param_summary)
setwd("/Users/savannah/Documents/Ch2/FNA_pilot/ch2_Figs_and_tables/")
write.csv(param_summary, "Table_S3_m5c_param_summary.csv")

#######      Implementation: Plot results of best model      #######

# 5. Predict recovery trajectories - takes a long time bc I added more height vals to predict
recovery_data <- predict_recovery_trajectory(mcmc_samples, best_model_name)

# 6. Plot recovery trajectories
recovery_plot <- plot_recovery_trajectories(recovery_data, height_value = 25)
print(recovery_plot)

# 7. Create height vs temperature plots at different recovery times
plot_height_temp_relationship <- function(trajectory_data, recovery_time = 10) {
  # Filter data for the specified recovery time
  plot_data <- trajectory_data %>% 
    filter(RecoveryTime == recovery_time)
  
  # Create plot
  ggplot(plot_data, aes(x = Height, y = MeanTemp, color = Forest, fill = Forest)) +
    geom_line(linewidth = 1) +
    geom_ribbon(aes(ymin = Q2.5, ymax = Q97.5), alpha = 0.2, linetype = 0) +
    scale_color_manual(values = c("Intact" = "green", "Logged" = "blue", "Burned" = "red")) +
    scale_fill_manual(values = c("Intact" = "green", "Logged" = "blue", "Burned" = "red")) +
    labs(
      x = "Canopy Height (m)",
      y = "Predicted Temperature (°C)",
      title = paste0("Temperature vs Height Relationship at ", recovery_time, " Years Since Disturbance")
    ) +
    theme_minimal()
}
# recovery_times could be 1,5,10,20,30 
height_temp_plot <- plot_height_temp_relationship(recovery_data, recovery_time = 5)
print(height_temp_plot)

#######      Implementation: Explore & evaluate best model      #######

# 8. Plot parameter posterior distributions - focus only on relevant parameters
plot_posterior_distributions <- function(mcmc_samples) {
  # Convert to matrix for easier manipulation
  mcmc_matrix <- as.matrix(do.call(rbind, mcmc_samples))
  
  # Select parameters including c_mean and z_mean
  interaction_params <- c(
    "c_mean", "z_mean",
    "beta_c_burned", "beta_c_logged", "beta_z_burned", "beta_z_logged",
    "beta_c_burned_time", "beta_c_logged_time", "beta_z_burned_time", "beta_z_logged_time"
  )
  
  # Keep only parameters that exist in the model
  interaction_params <- interaction_params[interaction_params %in% colnames(mcmc_matrix)]
  
  # Create a dataframe for ggplot
  posterior_data <- data.frame()
  
  for(param in interaction_params) {
    temp_data <- data.frame(
      Parameter = param,
      Value = mcmc_matrix[, param]
    )
    posterior_data <- rbind(posterior_data, temp_data)
  }
  
  # Create plot
  ggplot(posterior_data, aes(x = Value, fill = Parameter)) +
    geom_density(alpha = 0.7) +
    facet_wrap(~ Parameter, scales = "free") +
    geom_vline(xintercept = 0, linetype = "dashed") +
    labs(
      x = "Parameter Value",
      y = "Density",
      title = "Posterior Distributions of Base Parameters and Disturbance Effects"
    ) +
    theme_minimal() +
    theme(legend.position = "none")
}

posterior_plot <- plot_posterior_distributions(mcmc_samples)
print(posterior_plot)

# 9. Calculate recovery time to reach intact forest conditions - ignoring main time effect
estimate_recovery_time <- function(mcmc_samples) {
  # Convert to matrix for easier manipulation
  mcmc_matrix <- as.matrix(do.call(rbind, mcmc_samples))
  
  # Placeholder for results
  results <- data.frame(
    ForestType = c("Burned", "Logged"),
    RecoveryYears_c = NA,
    RecoveryYears_z = NA
  )
  
  # Number of samples to use
  n_samples <- min(1000, nrow(mcmc_matrix))
  
  # For burned forests - using only the interaction effect
  # beta_c_burned + beta_c_burned_time * time_std = 0
  recovery_times_c <- -mcmc_matrix[1:n_samples, "beta_c_burned"] / 
    mcmc_matrix[1:n_samples, "beta_c_burned_time"]
  
  # Convert from standardized to actual years
  recovery_times_c <- recovery_times_c * sd(sub.df$timeSinceDist, na.rm=TRUE) + 
    mean(sub.df$timeSinceDist, na.rm=TRUE)
  
  # Use median of estimated recovery times
  results$RecoveryYears_c[1] <- median(recovery_times_c[recovery_times_c > 0 & 
                                                          recovery_times_c < 100],
                                       na.rm = TRUE)
  
  # Estimate recovery time for z parameter - using only the interaction effect
  recovery_times_z <- -mcmc_matrix[1:n_samples, "beta_z_burned"] / 
    mcmc_matrix[1:n_samples, "beta_z_burned_time"]
  
  recovery_times_z <- recovery_times_z * sd(sub.df$timeSinceDist, na.rm=TRUE) + 
    mean(sub.df$timeSinceDist, na.rm=TRUE)
  
  results$RecoveryYears_z[1] <- median(recovery_times_z[recovery_times_z > 0 & 
                                                          recovery_times_z < 100],
                                       na.rm = TRUE)
  
  # Repeat for logged forests - using only the interaction effect
  recovery_times_c <- -mcmc_matrix[1:n_samples, "beta_c_logged"] / 
    mcmc_matrix[1:n_samples, "beta_c_logged_time"]
  
  recovery_times_c <- recovery_times_c * sd(sub.df$timeSinceDist, na.rm=TRUE) + 
    mean(sub.df$timeSinceDist, na.rm=TRUE)
  
  results$RecoveryYears_c[2] <- median(recovery_times_c[recovery_times_c > 0 & 
                                                          recovery_times_c < 100],
                                       na.rm = TRUE)
  
  recovery_times_z <- -mcmc_matrix[1:n_samples, "beta_z_logged"] / 
    mcmc_matrix[1:n_samples, "beta_z_logged_time"]
  
  recovery_times_z <- recovery_times_z * sd(sub.df$timeSinceDist, na.rm=TRUE) + 
    mean(sub.df$timeSinceDist, na.rm=TRUE)
  
  results$RecoveryYears_z[2] <- median(recovery_times_z[recovery_times_z > 0 & 
                                                          recovery_times_z < 100],
                                       na.rm = TRUE)
  
  return(results)
}

recovery_time_estimates <- estimate_recovery_time(mcmc_samples)
print(recovery_time_estimates)

# 10. Calculate model fit statistics and create diagnostic plots
analyze_model_fit <- function(jags_fit, jags_data, sub.df) {
  # Sample fitted values
  fitted_samples <- coda.samples(jags_fit,
                                 variable.names = c("mu"),
                                 n.iter = 1000)
  
  # Calculate mean fitted values
  fitted_matrix <- as.matrix(do.call(rbind, fitted_samples))
  fitted_values <- colMeans(fitted_matrix)
  
  # Calculate residuals
  residuals <- jags_data$y - fitted_values
  
  # Create diagnostic dataframe
  diag_data <- data.frame(
    Fitted = fitted_values,
    Residuals = residuals,
    Forest_type = sub.df$Forest_cat_v2,
    x = sub.df$x,
    y = sub.df$y,
    Height = sub.df$rh95,
    TimeSinceDist = sub.df$timeSinceDist
  )
  
  # Create residual plot - separate by forest type
  residual_plot <- ggplot(diag_data, aes(x = Fitted, y = Residuals, color = Forest_type)) +
    geom_point(alpha = 0.5) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    scale_color_manual(values = c("Intact" = "green", "Logged" = "blue", "Burned" = "red")) +
    labs(
      x = "Fitted Values",
      y = "Residuals",
      title = "Residuals vs Fitted Values"
    ) +
    theme_minimal()
  
  # Calculate summary statistics for model evaluation
  rmse <- sqrt(mean(residuals^2))
  mae <- mean(abs(residuals))
  r_squared <- 1 - sum(residuals^2) / sum((jags_data$y - mean(jags_data$y))^2)
  
  # Summary by forest type
  forest_stats <- diag_data %>%
    group_by(Forest_type) %>%
    summarise(
      RMSE = sqrt(mean(Residuals^2)),
      MAE = mean(abs(Residuals)),
      Mean_Residual = mean(Residuals),
      .groups = 'drop'
    )
  
  # Create spatial residual plot
  spatial_plot <- ggplot(diag_data, aes(x = x, y = y, color = Residuals)) +
    geom_point() +
    scale_color_gradient2(low = "blue", mid = "yellow", high = "red", midpoint = 0) +
    labs(
      x = "Longitude",
      y = "Latitude",
      title = "Spatial Distribution of Residuals"
    ) +
    theme_minimal()
  
  # Create residuals vs height plot
  height_plot <- ggplot(diag_data, aes(x = Height, y = Residuals, color = Forest_type)) +
    geom_point(alpha = 0.5) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_smooth(method = "loess", se = FALSE) +
    scale_color_manual(values = c("Intact" = "green", "Logged" = "blue", "Burned" = "red")) +
    labs(
      x = "Canopy Height (m)",
      y = "Residuals",
      title = "Residuals vs Canopy Height"
    ) +
    theme_minimal()
  
  # Create residuals vs time since disturbance plot - only for disturbed forests
  time_plot <- ggplot(diag_data[diag_data$Forest_type != "Intact", ], 
                      aes(x = TimeSinceDist, y = Residuals, color = Forest_type)) +
    geom_point(alpha = 0.5) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_smooth(method = "loess", se = FALSE) +
    scale_color_manual(values = c("Logged" = "blue", "Burned" = "red")) +
    labs(
      x = "Time Since Disturbance (years)",
      y = "Residuals",
      title = "Residuals vs Time Since Disturbance (Disturbed Forests Only)"
    ) +
    theme_minimal()
  
  # Return results
  return(list(
    fit_stats = data.frame(RMSE = rmse, MAE = mae, R_squared = r_squared),
    forest_stats = forest_stats,
    residual_plot = residual_plot,
    spatial_plot = spatial_plot,
    height_plot = height_plot,
    time_plot = time_plot
  ))
}

# Create complete model summary report with revised interpretation
create_model_report <- function(param_summary, recovery_data, recovery_time_estimates, fit_stats) {
  
  cat("\n==== BEST MODEL PARAMETERS ====\n")
  # Filter out the main time since disturbance effects
  param_summary_filtered <- param_summary[!grepl("^beta_[cz]_timeSinceDist$", param_summary$Parameter), ]
  print(param_summary_filtered)
  
  cat("\n==== RECOVERY TIME ESTIMATES ====\n")
  print(recovery_time_estimates)
  
  cat("\n==== MODEL FIT STATISTICS ====\n")
  print(fit_stats)
  
  cat("\n==== ECOLOGICAL INTERPRETATION ====\n")
  
  # Extract significant parameters, excluding main timeSinceDist effects
  sig_params <- param_summary_filtered[param_summary_filtered$Significant, ]
  
  cat("Significant effects on canopy temperature:\n")
  for(i in 1:nrow(sig_params)) {
    param <- sig_params$Parameter[i]
    effect <- sig_params$Mean[i]
    
    # Interpret parameter
    if(grepl("^beta_c_burned$", param)) {
      cat("- ", param, ": Immediate effect of burning on baseline temperature = ", round(effect, 4), "\n")
    } else if(grepl("^beta_c_logged$", param)) {
      cat("- ", param, ": Immediate effect of logging on baseline temperature = ", round(effect, 4), "\n")
    } else if(grepl("^beta_z_burned$", param)) {
      cat("- ", param, ": Immediate effect of burning on temperature-height scaling = ", round(effect, 4), "\n")
    } else if(grepl("^beta_z_logged$", param)) {
      cat("- ", param, ": Immediate effect of logging on temperature-height scaling = ", round(effect, 4), "\n")
    } else if(grepl("^beta_c_burned_time$", param)) {
      cat("- ", param, ": Recovery rate of baseline temperature in burned forests = ", round(effect, 4), "\n")
    } else if(grepl("^beta_c_logged_time$", param)) {
      cat("- ", param, ": Recovery rate of baseline temperature in logged forests = ", round(effect, 4), "\n")
    } else if(grepl("^beta_z_burned_time$", param)) {
      cat("- ", param, ": Recovery rate of temperature-height scaling in burned forests = ", round(effect, 4), "\n")
    } else if(grepl("^beta_z_logged_time$", param)) {
      cat("- ", param, ": Recovery rate of temperature-height scaling in logged forests = ", round(effect, 4), "\n")
    } else if(grepl("^beta_c_air_temp$", param)) {
      cat("- ", param, ": Effect of air temperature on baseline canopy temperature = ", round(effect, 4), "\n")
    } else if(grepl("^beta_z_elevation$", param)) {
      cat("- ", param, ": Effect of elevation on temperature-height scaling = ", round(effect, 4), "\n")
    }
  }
  
  # Interpret recovery patterns
  cat("\nBurned forest recovery:\n")
  if("beta_c_burned_time" %in% sig_params$Parameter) {
    effect <- sig_params$Mean[sig_params$Parameter == "beta_c_burned_time"]
    if(effect > 0) {
      cat("- Baseline temperature increases with recovery time\n")
    } else {
      cat("- Baseline temperature decreases with recovery time\n")
    }
  } else {
    cat("- No significant change in baseline temperature with recovery time\n")
  }
  
  if("beta_z_burned_time" %in% sig_params$Parameter) {
    effect <- sig_params$Mean[sig_params$Parameter == "beta_z_burned_time"]
    if(effect > 0) {
      cat("- Temperature-height relationship becomes steeper with recovery time\n")
    } else {
      cat("- Temperature-height relationship becomes flatter with recovery time\n")
    }
  } else {
    cat("- No significant change in temperature-height relationship with recovery time\n")
  }
  
  cat("\nLogged forest recovery:\n")
  if("beta_c_logged_time" %in% sig_params$Parameter) {
    effect <- sig_params$Mean[sig_params$Parameter == "beta_c_logged_time"]
    if(effect > 0) {
      cat("- Baseline temperature increases with recovery time\n")
    } else {
      cat("- Baseline temperature decreases with recovery time\n")
    }
  } else {
    cat("- No significant change in baseline temperature with recovery time\n")
  }
  
  if("beta_z_logged_time" %in% sig_params$Parameter) {
    effect <- sig_params$Mean[sig_params$Parameter == "beta_z_logged_time"]
    if(effect > 0) {
      cat("- Temperature-height relationship becomes steeper with recovery time\n")
    } else {
      cat("- Temperature-height relationship becomes flatter with recovery time\n")
    }
  } else {
    cat("- No significant change in temperature-height relationship with recovery time\n")
  }
  
  # Compare recovery rates between forest types
  if(!is.na(recovery_time_estimates$RecoveryYears_c[1]) && !is.na(recovery_time_estimates$RecoveryYears_c[2])) {
    cat("\nEstimated years to recover baseline temperature (c):\n")
    cat("- Burned forests:", round(recovery_time_estimates$RecoveryYears_c[1], 1), "years\n")
    cat("- Logged forests:", round(recovery_time_estimates$RecoveryYears_c[2], 1), "years\n")
    
    if(recovery_time_estimates$RecoveryYears_c[1] < recovery_time_estimates$RecoveryYears_c[2]) {
      cat("- Burned forests recover baseline temperature faster than logged forests\n")
    } else {
      cat("- Logged forests recover baseline temperature faster than burned forests\n")
    }
  }
  
  if(!is.na(recovery_time_estimates$RecoveryYears_z[1]) && !is.na(recovery_time_estimates$RecoveryYears_z[2])) {
    cat("\nEstimated years to recover temperature-height relationship (z):\n")
    cat("- Burned forests:", round(recovery_time_estimates$RecoveryYears_z[1], 1), "years\n")
    cat("- Logged forests:", round(recovery_time_estimates$RecoveryYears_z[2], 1), "years\n")
    
    if(recovery_time_estimates$RecoveryYears_z[1] < recovery_time_estimates$RecoveryYears_z[2]) {
      cat("- Burned forests recover temperature-height relationship faster than logged forests\n")
    } else {
      cat("- Logged forests recover temperature-height relationship faster than burned forests\n")
    }
  }
}

# Run diagnostic analysis - takes a while 
fit_diagnostics <- analyze_model_fit(jags_fit, jags_data, sub.df)
print(fit_diagnostics$fit_stats)
print(fit_diagnostics$residual_plot)
print(fit_diagnostics$spatial_plot)
print(fit_diagnostics$height_plot)
print(fit_diagnostics$time_plot)

# Generate model report
model_report <- create_model_report(
  param_summary,
  recovery_data,
  recovery_time_estimates,
  fit_diagnostics$fit_stats
)

#write(model_report, "Table_S5_m5c_model_report.txt")



#######      Fig 3 v2 (faster): Temperatures by forest type and height    ########
# Calculate for a specific recovery time (e.g., 5 years)
temp_differences_5yr <- calculate_temperature_differences(mcmc_samples, recovery_time = 5)

# Print results
print(paste("Results for", temp_differences_5yr$RecoveryTime, "years since disturbance:"))
print("Temperatures by forest type and height:");print(temp_differences_5yr$Temperatures)
print("Temperature differences between forest types:");print(temp_differences_5yr$Differences)
print("Overall temperature difference summary across heights:");print(temp_differences_5yr$Overall_Summary)

# Write results to CSV files
write.csv(temp_differences_5yr$Temperatures, 
          paste0("Table_S4_prep_temperatures_by_height_", temp_differences_5yr$RecoveryTime, "yr.csv"), 
          row.names = FALSE)
write.csv(temp_differences_5yr$Differences, 
          paste0("temperature_differences_by_height_", temp_differences_5yr$RecoveryTime, "yr.csv"), 
          row.names = FALSE)
write.csv(temp_differences_5yr$Overall_Summary, 
          paste0("Table_S4_overall_temperature_differences_", temp_differences_5yr$RecoveryTime, "yr.csv"), 
          row.names = FALSE)

# Plot temperatures by forest type
ggplot(temp_differences_5yr$Temperatures, aes(x = Height, y = Mean_Temp, color = Forest_Type, fill = Forest_Type)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = Lower_CI, ymax = Upper_CI), alpha = 0.2, linetype = 0) +
  scale_color_manual(values = c("Intact" = "green", "Logged" = "blue", "Burned" = "red")) +
  scale_fill_manual(values = c("Intact" = "green", "Logged" = "blue", "Burned" = "red")) +
  labs(
    x = "Canopy Height (m)",
    y = "Predicted Temperature (°C)",
    title = paste("Predicted Temperatures at", temp_differences_5yr$RecoveryTime, "Years Since Disturbance")
  ) +
  theme_minimal()

# Plot temperature differences
ggplot(temp_differences_5yr$Differences, aes(x = Height, y = Mean_Diff, color = Comparison, fill = Comparison)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = Lower_CI, ymax = Upper_CI), alpha = 0.2, linetype = 0) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_color_manual(values = c("Burned vs Intact" = "darkred", "Logged vs Intact" = "darkblue", "Burned vs Logged" = "purple")) +
  scale_fill_manual(values = c("Burned vs Intact" = "darkred", "Logged vs Intact" = "darkblue", "Burned vs Logged" = "purple")) +
  labs(
    x = "Canopy Height (m)",
    y = "Temperature Difference (°C)",
    title = paste("Temperature Differences at", temp_differences_5yr$RecoveryTime, "Years Since Disturbance")
  ) +
  theme_minimal()
