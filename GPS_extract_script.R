
##Load Packages
library(dplyr)
library(tidyr)
library(lubridate)
library(ggplot2)
library(ggthemes)
library(raster)
library(sp)
library(rgdal)
library(sf)

#write csv so we can save the steps above in future
trajectory <- read.csv("C:/Users/sean.konkolics/Downloads/trajectory.csv")

# Prepare trajectory df
# Convert FormattedDateTime to POSIXct
trajectory$FormattedDateTime <- as.POSIXct(trajectory$FormattedDateTime, format = "%d.%m.%Y %H:%M")

##asign lat/lon CRS and then convert to UTM CRS to match the TIFF Files
coordinates(trajectory) <- c("Longitude", "Latitude")
proj4string(trajectory) <- CRS("+proj=longlat +datum=WGS84")
trajectory <- spTransform(trajectory, CRS("+proj=utm +zone=11 ellps=WGS84"))
st_crs(trajectory)


#check for coords with decimal nubmers
unique_bears_with_decimal_coords <- unique(final[grepl("\\.\\d+", final$coords.x1), "Bear"])
print(unique_bears_with_decimal_coords)


##I think from here you would be able to switch back to your script for the extractions. I did something slightly different below
##that made more sense to me. 

##import aspect TIFF
asp <- raster("Z:/_Shared/Research/Volunteers/Volunteer Agreements/ACTIVE Volunteers/Alyssa Bohart/DenningTransfer/TIFFS_RSFLayers/Aspect1.tif")
st_crs(asp)
plot(asp)

##import shrub TIFF
shrub <- raster("Z:/_Shared/Research/Volunteers/Volunteer Agreements/ACTIVE Volunteers/Alyssa Bohart/DenningTransfer/TIFFS_RSFLayers/Shrub1.tif")
st_crs(shrub)
plot(shrub)

##import all of the TIFFs you want to use
next_raster <- raster("Z:/_Shared/Research/Volunteers/Volunteer Agreements/ACTIVE Volunteers/Alyssa Bohart/DenningTransfer/TIFFS_RSFLayers/Aspect1.tif")
st_crs(next_raster)
plot(next_raster)

##extract the TIFF values to the GPS points. Do each extraction seperatly and then combine into a single dataframe.
asp_extract <- extract(asp, trajectory)
shrub_extract <- extract(shrub,trajectory)
combine_df <- cbind(trajectory,asp_extract, shrub_extract)
final <- as.data.frame(combine_df)

##Rename the columns to match that TIFF file names
final <- final %>% 
  rename.........

