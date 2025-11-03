#########       Install packages, prep & notes        ###########
library(ggplot2)
library(lubridate)
library(dplyr)
library('sjPlot') 
root="/Users/savannah/Documents/Ch2/FNA_pilot/INMET_station_data/"
setwd(root)
set_theme(base = theme_classic(), #theme_minimal(), #axis.tickslen = 0, # hides tick marks
          axis.title.size = 2.4, axis.textsize = 2.8, #legend.size = .7,  legend.title.size = .8,
          geom.label.size = 2.5, title.size = 2.1)

# Goal: Plot daily maximum and minimum air temperature averaged from two meteorological stations located near study area 
# (Sinop, lat: -11.98, lon: -55.57; and Sorriso, lat: -12.55, lon: -55.72) 
# from 07/01/2018 through 07/01/2024. 

#########       Prep INMET station data        ###########
df_all<- data.frame()

# Aggregate data into df_all
for(year in 2018:2023){
  print(year)
  df1<- read.csv(paste0("Sorriso_",year,".csv"), header = T, sep=";") 
  df1<- df1[,1:5] # Only focused on temperature for now
  df1<- data.frame(df1)
  names(df1)<- c("Date", "Hour_UTC", "Temp_inst","Temp_max", "Temp_min")
  
  # Need to convert numbers formatted as characters such as "26,7" to numeric values (eg., 26.7)
  # df1 columns 3-5 have numbers with commas as decimal separators
  df1[,3:5] <- apply(df1[,3:5], 2, function(x) as.numeric(gsub(",", ".", x)))
  df1$station<- "Sorriso"
  
  df2<- read.csv(paste0("Sinop_",year,".csv"), header = T, sep=";") 
  df2<- df2[,1:5] # Only focused on temperature for now
  df2<- data.frame(df2)
  names(df2)<- c("Date", "Hour_UTC", "Temp_inst","Temp_max", "Temp_min")
  df2[,3:5] <- apply(df2[,3:5], 2, function(x) as.numeric(gsub(",", ".", x)))
  
  df2$station<- "Sinop"
  
  
  df<- rbind(df1, df2)
  df_all<- rbind(df_all, df)
}

write.csv(df_all, "INMET_dfall_temp_2018_2024.csv")
#########       Plot INMET station data        ###########

df_all<- read.csv("INMET_dfall_temp_2018_2024.csv")

# Data preparation - now keeping station-specific maximum temperatures
daily_temps <- df_all %>%
  # Convert Date to Date type
  mutate(Date = dmy(Date)) %>%
  # Group by date and station and get daily max only
  group_by(Date, station) %>%
  summarize(
    Daily_max = if(all(is.na(Temp_max))) NA_real_ else max(Temp_max, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  # Filter out the NA period
  filter(!(Date >= as.Date("2020-11-01") & Date <= as.Date("2022-11-01")))

highlight_date <- as.Date("2018-08-31")

# Create the plot
ggplot(daily_temps, aes(x = Date, color = station)) +
  # Add lines for station-specific max temperatures
  geom_line(aes(y = Daily_max), size = 0.8) +
  # Add vertical line for highlighted date
  geom_vline(xintercept = as.numeric(highlight_date),
             color = "darkred", linetype = "dashed", size = 0.8) +
  # Add annotation for highlighted date
  annotate("text",
           x = highlight_date,
           y = max(daily_temps$Daily_max, na.rm = TRUE),
           label = "August 31, 2018",
           hjust = -0.1,
           vjust = 1,
           color = "darkred") +
  # Customize colors
  scale_color_manual(values = c("Sinop" = "purple3", "Sorriso" = "red")) +
  # Labels and title
  labs(
    title = "Daily Maximum Temperature (Sorriso and Sinop)",
    x = "Date",
    y = "Temperature (°C)",
    color = "Station"
  ) +
  # Theme customization
  theme_minimal() +
  theme(
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12),
    panel.grid.minor = element_line(color = "gray90")
  ) +
  # Format date axis
  scale_x_date(
    date_breaks = "3 months",
    date_labels = "%b %Y",
    expand = expansion(mult = c(0.02, 0.1))
  )
