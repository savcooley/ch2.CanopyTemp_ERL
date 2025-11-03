#################       Install packages & Notes on forest canopy temp distributions        ###########
library(fitdistrplus)
library(rjags)
library(ggplot2)
library('sjPlot') 
library(dplyr)
root="/Users/savannah/Documents/Ch2/FNA_pilot/"

set_theme(base = theme_classic(), #theme_minimal(), #axis.tickslen = 0, # hides tick marks
          axis.title.size = 1.8, axis.textsize = 1.3, #legend.size = .7,  legend.title.size = .8,
          geom.label.size = 2.5, title.size = 2.1)

# ch2 Upper Canopy Temperature project. 
# This code models upper canopy leaf temperature distributions of warmest forest patches among intact, logged and burned forests.
# "Warmest" forest patches are defined as top 25th percentile pixels from an ECOSTRESS land surface temperature image acquired 
# on one of the warmest, least cloudy days observed in the 2018-2020 dry season time series in the study area:

# 31-August-2018 (2018_243 Julian date) at 17:38:03 UTC / 13:38:03 local time. 
# Alternative image: 31-July-2018 (2018212 Julian date) captured by ECOSTRESS at 16:44:54 UTC

# The JAGS code parameterizes three gamma distributions of leaf temperatures for each forest type with mean temperatires based on 
# p90 ECOSTRESS LST values from 31-August-2018 aquistion over Feliz Natal study area. 
# The spread of each gamma distribution is estimated based on leaf thermocouple measurements from numerous studies aggregated by Doughty et al (2023) gamma curves.  

df_all <- read.csv("/Users/savannah/Documents/Ch2/FNA_pilot/all_vars_FNA_v2.csv")
sub.df <- df_all 

# Remove unnecessary variables
sub.df$LST_2018212_16 <- NULL 
#sub.df$LST_2018_243_17 <- NULL # This is the warmest, least cloud contaminated day in the 2-yr dry season ECOSTRESS record 
sub.df$LST_2020_262_19 <- NULL
sub.df$MeanMaxLST<- NULL

# Remove rows with missing values & outliers
sub.df <- na.omit(sub.df) # From 18k to 8k obs
sub.df %>%
  group_by(Forest_cat) %>%
  summarise(
    AGB_mean = mean(agbd_Mg_ha, na.rm = T),
    AGB_sd = sd(agbd_Mg_ha, na.rm = T),
    LST_mean = mean(LST_2018_243_17, na.rm = T),
    LST_sd = sd(LST_2018_243_17, na.rm = T),
    LST_p90 = quantile(LST_2018_243_17, probs = c(0.9), na.rm = T),
    n = n()) 

# Define the thresholds based on Forest_cat
sub.df <- sub.df %>% # From 8k to 6912 obs
  filter(  # Based on Pinage et al 2023
    # Lower bounds - remove observations with potential land use change since 2020 classification
      (Forest_cat == "Intact" & agbd_Mg_ha > 60) | 
      (Forest_cat == "Logged" & agbd_Mg_ha > 40) | 
      (Forest_cat == "Burned" & agbd_Mg_ha > 5) )

sub.df %>%
  group_by(Forest_cat) %>%
  summarise(
    AGB_mean = mean(agbd_Mg_ha, na.rm = T),
    AGB_sd = sd(agbd_Mg_ha, na.rm = T),
    rh95_mean = mean(rh95, na.rm = T),
    rh95_sd = sd(rh95, na.rm = T),
    # LST_median = median(LST_2018_243_17, na.rm = T),
    # LST_mean = mean(LST_2018_243_17, na.rm = T),
    # LST_sd = sd(LST_2018_243_17, na.rm = T),
    # LST_p75 = quantile(LST_2018_243_17, probs = c(0.75), na.rm = T),
    # LST_p90 = quantile(LST_2018_243_17, probs = c(0.9), na.rm = T),
    n = n()) 


# Remove rows with probability of bare ground > 4% (based on Dynamic World land cover classification using Sentinel2 imagery) 
sub.df <- subset(sub.df, probBareGround.2018<0.04) # Additional 2.1% reduction in n: from 6912 to 6767 obs

# sub.df data - for LST_2018_243_17 

# Forest_cat AGB_mean AGB_sd LST_median LST_sd LST_p75 LST_p90     n
# 1 Burned         102.   79.3       36.2  1.81     37.3    38.1   911
# 2 Intact         122.   72.1       34.8  0.877    35.2    35.9  1859
# 3 Logged         120.   72.0       35.0  0.915    35.4    35.9  3997

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# df_all data - for LST_2018212_16
# Forest_cat AGB_mean AGB_sd LST_mean LST_sd LST_p90 LST_p95
# 1 Burned         88.4   79.9     34.1   2.02    36.3    37.2
# 2 Intact         99.0   74.6     32.8   1.32    33.8    34.3
# 3 Logged        107.    73.7     32.5   1.62    33.5    34.3

#'
###########       v8 - JAGS only method        ###########
# Leverage the leaf thermocouple data from Doughty et al 

# Started with using Doughty et al Fig 1c WPD extracted data 
csv_data <- read.csv("/Users/savannah/Documents/Ch2/Doughty_et_al_data/leaftempall_hist_WPD.csv")
combined_data <- rep(csv_data$LeafTemp, round(csv_data$Count)) 
combined_mean <- mean(combined_data) # 25.43 degrees C
combined_sd <- sd(combined_data) # 2.85 degrees C

# Using Fig 1b data (single site) instead of 1c (multiple thermocouple studies)
# b, Histogram of individual canopy top leaf thermocouples from 11 individual leaves from the same site as in a over 54 sunny periods lasting 20 min
LeafTemp<- c(29, 32, 35, 38, 41, 44, 47)
Count<- c(75, 350, 200, 38, 72, 10, 2) 
csv_data <- data.frame(LeafTemp, Count) 
combined_data <- rep(csv_data$LeafTemp, round(csv_data$Count)) 

# Calculate mean and sd of combined data
combined_mean <- mean(combined_data) # 33.87 degrees C
combined_sd <- sd(combined_data) # 3.46 degrees C

# This provides evidence that higher mean leaf temp --> higher SD in leaf temp. 
# I will assume a linear relationship and use this to modify the SD of predicted leaf temp among forest degradation types

# Create linear model from known points
temp_means <- c(25.43, 33.87)
temp_sds <- c(2.85, 3.46)

# Fit linear model
lm_temp <- lm(temp_sds ~ temp_means)

# Extract slope and intercept
slope <- coef(lm_temp)[2]  # 0.07236
intercept <- coef(lm_temp)[1]  # 1.01168

# Create updated forest data with predicted SDs
forest_data <- data.frame(
  Forest_cat = c("Intact", "Logged", "Burned"),
  meanCanopyTemp = c(35.2, 35.4, 37.3), # p75  # c(35.9, 35.9, 38.1) is p90 for AGB-filtered LST_2018_243_17.
  sd = c(35.2 * slope + intercept,
         35.4 * slope + intercept,
         37.3 * slope + intercept)
)

# Function to generate sample data based on combined data and forest-specific stats
generate_sample <- function(combined_data, target_mean, target_sd, n = 10000) {
  shift <- target_mean - mean(combined_data)
  scale <- target_sd / sd(combined_data)
  sample <- (combined_data - mean(combined_data)) * scale + target_mean
  sample[sample <= 0] <- target_mean  # Replace negative values with target mean
  return(sample[1:n])  # Ensure we return exactly n samples
}

# Estimate parameters for each forest type using JAGS
forest_params <- lapply(1:nrow(forest_data), function(i) {
  sample_data <- generate_sample(combined_data, forest_data$meanCanopyTemp[i], forest_data$sd[i])
  
  # JAGS method
  jags_model <- "
  model {
    for (i in 1:n) {
      y[i] ~ dgamma(shape, rate)
    }
    shape ~ dunif(0, 100)
    rate ~ dunif(0, 10)
  }
  "
  jags_data <- list(y = sample_data, n = length(sample_data))
  jags_fit <- jags.model(textConnection(jags_model), data = jags_data, n.chains = 3)
  update(jags_fit, 1000)  # Burn-in
  mcmc_samples <- coda.samples(jags_fit, variable.names = c("shape", "rate"), n.iter = 5000)
  jags_results <- summary(mcmc_samples)$statistics[, "Mean"]
  
  return(list(shape = jags_results["shape"], rate = jags_results["rate"]))
})

# Create plot data
temp_seq <- seq(25, 54, by = 0.1)
plot_data <- expand.grid(Temperature = temp_seq,
                         Forest_cat = forest_data$Forest_cat)
plot_data$Density <- NA

# Add estimated distributions for each forest type
for (i in 1:nrow(forest_data)) {
  shape <- forest_params[[i]]$shape
  rate <- forest_params[[i]]$rate
  idx <- plot_data$Forest_cat == forest_data$Forest_cat[i]
  plot_data$Density[idx] <- dgamma(temp_seq, shape = shape, rate = rate)
}

# Add actual data for combined forests
actual_data <- data.frame(
  Temperature = csv_data$LeafTemp,
  Density = csv_data$Count / sum(csv_data$Count),
  Forest_cat = "Doughty et al."
)
#plot_data <- rbind(plot_data, actual_data)

# Add temperature thresholds
plot_data$Tar <- 33.0
plot_data$Tcrit <- 46.7
plot_data$T50 <- 49.9

# Create the plot
ggplot(plot_data, aes(x = Temperature, y = Density, color = Forest_cat)) +
  geom_line(size = 1) +
  geom_point(data = subset(plot_data, Forest_cat == "Doughty et al."), size = 2) +
  geom_vline(aes(xintercept = Tar), color = "gold", size = 1.5,linetype="dotted") +
  geom_vline(aes(xintercept = Tcrit), color = "orange2", size = 1.5,linetype="dotted") +
  geom_vline(aes(xintercept = T50), color = "red4", size = 1.5,linetype="dotted") +
  scale_color_manual(values = c("Intact" = "green", "Logged" = "blue", "Burned" = "red", "Doughty et al." = "purple")) +
  labs(title = "Leaf Temperature Distributions by Forest Type",
       #subtitle = "JAGS estimates of max. ECOSTRESS LST with Thermocouple data from Doughty et al. (2023)",
       x = "Temperature (°C)",
       y = "Density",
       color = "Forest Category") +
  theme_minimal() +
  theme(legend.position = "bottom")

# Print JAGS results
results <- data.frame(
  Forest_cat = forest_data$Forest_cat,
  JAGS_shape = sapply(forest_params, function(x) x$shape),
  JAGS_rate = sapply(forest_params, function(x) x$rate)
)

print(results)


# Calculate and print the mean and variance for each fitted distribution
calculate_mean_var <- function(shape, rate) {
  mean <- shape / rate
  variance <- shape / (rate^2)
  return(c(mean = mean, variance = variance))
}

mean_var <- t(mapply(calculate_mean_var, 
                     results$JAGS_shape, 
                     results$JAGS_rate))
results$JAGS_mean <- mean_var[, "mean"]
results$JAGS_variance <- mean_var[, "variance"]

print(results[, c("Forest_cat", "JAGS_mean", "JAGS_variance")])
print(forest_data[, c("Forest_cat", "meanCanopyTemp", "sd")])

# Print combined data statistics
cat("Combined data: mean =", combined_mean, ", sd =", combined_sd, "\n")

#################       Compute Tcrit probabilities        ###########

# compute p(leafTemp > 46.7 ºC 0.0092%) across forest types
Tcrit= 49.9 # 46.7 #33.0 # #;method = "jags" #"moments" # "jags" "mle"
 
for(i in 1:3){ # 1 is intact, 2 is logged, 3 is burned
  shape <- forest_params[[i]]$shape; rate <- forest_params[[i]]$rate # forest_params[[i]][[method]]$rate
  
  # probability that a Gamma(shape = shape, rate = rate) random variable is less than Tcrit:
  cat(1-pgamma(Tcrit, shape = shape, rate = rate)) # subtract from 1 to get greater than
  
}

# probability that a Gamma(1,1) random variable lies between 1 and 3:
pgamma(3,1)-pgamma(1,1) # or more compactly diff(pgamma(c(1,3),1))
# This is the integral of the Gamma between these points.

#' Possible next steps: 
#' determine function to modify sd as a funciton of mean canopy temp
#' map of p(leafTemp > 46.7 ºC 0.0092%) across forest types

