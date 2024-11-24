
# Load the necessary packages
# library(rgdal)
# library(amt)
library(plyr)
library(dplyr)
library(ggplot2)
library(sf)
library(lubridate)
library(tibble)
library(sp)
library(raster)
library(gridExtra)
library(lme4)
library(caret)
library(tidyverse)
library(Hmisc)
library(AICcmodavg)
library(dotwhisker)
library(odbc)
library(anytime)
library(traj)
library(glmmTMB)
library(adehabitatHR) #important
library(adehabitatLT)
library(readxl)
library(geosphere)
library(stringr)
library(readr)
library(survival)
library(MuMIn)
library(HDInterval)
library(bayesplot) # For posterior predictive checks and hdi calculations, adjust as necessary if using another package

# LOAD BRMS
# install.packages("brms")
library(brms)
# install_cmdstan()
# install.packages('raster')
# install all the packages needed

# install.packages(c("rgdal", "amt", "plyr", "dplyr", "ggplot2", "sf", "lubridate", "tibble", "sp", "raster", "gridExtra", "lme4", "caret", "tidyverse", "Hmisc", "AICcmodavg", "dotwhisker", "odbc", "anytime", "traj", "glmmTMB", "adehabitatHR", "adehabitatLT", "readxl", "geosphere", "stringr", "readr", "survival", "MuMIn"))

#working directory
setwd("C:/Users/nanda/Desktop/Grizzly project/R/Data")




#write.csv(final2, file = "final2.csv", row.names = FALSE)   
final2 <- read.csv("final2.csv")       


#refinements
# Ensure categories and Bear IDs are factors
final2$Sex <- as.factor(final2$Sex)
final2$Day_periods <- as.factor(final2$Day_periods)
final2$Bear <- as.factor(final2$Bear)
final2$aspen <- as.factor(final2$aspen)
final2$barren <- as.factor(final2$barren)
final2$conifer <- as.factor(final2$conifer)
final2$shrub <- as.factor(final2$shrub)
final2$grass <- as.factor(final2$grass)
final2$Season <- as.factor(final2$Season)
final2$Day_periods <- as.factor(final2$Day_periods)

### recalculate things     

###recalculate time_diff (in h)
final2 <- final2 %>%
  arrange(Bear, FormattedDateTime)


# transform FormattedDateTime
final2$FormattedDateTime <- as.POSIXct(final2$FormattedDateTime, format = "%Y-%m-%d %H:%M:%S")

# Now, you can proceed with the time difference calculation
final2_time_diff <- final2 %>%
  arrange(Bear, FormattedDateTime) %>%
  group_by(Bear) %>%
  mutate(
    time_diff_hours = c(0, diff(na.omit(as.numeric(FormattedDateTime)))/3600)
  ) %>%
  ungroup()


#round hours
final2$FormattedDateTime <- as.POSIXct(final2$FormattedDateTime, format = "%Y-%m-%d %H:%M:%S")

# Calculate time differences in hours
final2_time_diff <- final2 %>%
  arrange(Bear, FormattedDateTime) %>%
  group_by(Bear) %>%
  mutate(
    time_diff = c(0, diff(na.omit(as.numeric(FormattedDateTime)))/3600),
    rounded_time_diff = round((time_diff * 60) / 15) * 15 / 60
  ) %>%
  ungroup()


#add time diff to final2
final2$time_diff <- final2_time_diff$rounded_time_diff



###recalculate step_length
#final2$step_length <- NULL
final2 <- final2 %>%
  arrange(Bear, FormattedDateTime)

# Calculate step length for each bear
final2_step_length <- final2 %>%
  group_by(Bear) %>%
  mutate(
    # Calculate differences in Easting and Northing
    delta_easting = c(NA, diff(easting)),
    delta_northing = c(NA, diff(northing)),
    
    # Calculate Euclidean distance (step length)
    step_length = sqrt(delta_easting^2 + delta_northing^2)
  ) %>%
  ungroup()

#add step_length to final2
final2$step_length <- final2_step_length$step_length

#replace each NA with 0         #don't do it, you will need the NAs to sort out entries with 0 time_diff!
# final2_no_na <- final2 %>%
#  mutate(step_length = ifelse(is.na(step_length), 0, step_length))

#add it to final2 
# final2$step_length <- final2_no_na$step_length


###recalculate turn angles 
final2$FormattedDateTime <- as.POSIXct(final2$FormattedDateTime, format = "%Y-%m-%d %H:%M:%S")

# Calculate turn angles based on easting and northing
final2_turn_angles <- final2 %>%
  group_by(Bear) %>%
  mutate(
    # Calculate the differences in easting and northing
    delta_easting = c(0, diff(easting)),
    delta_northing = c(0, diff(northing)),
    
    # Calculate the turn angles using atan2
    turn_angles = atan2(delta_northing, delta_easting) * 180 / pi
  ) %>%
  ungroup()

# add turn angles to final2
final2$turn_angles <- final2_turn_angles$turn_angles    

###calculate kmh
final2 <- final2 %>%
  group_by(Bear) %>%
  mutate(
    time_diff_seconds = c(NA, diff(as.numeric(FormattedDateTime))),
    speed_kmh = ifelse(is.na(time_diff_seconds), 0, (step_length / 1000) / (time_diff_seconds / 3600))
  ) %>%
  ungroup()

#remove time_diff in seconds
final2$time_diff_seconds <- NULL

#remove wrong movementspeed column
final2$movement_speed_kmh <- NULL


# Remove NAs
final2 <- final2[complete.cases(final2), ]

# dist ohv not needed
final2$distance_ohv <- NULL

#dist cover seems to have many errors
final2$distance_cover <- NULL


# Ensure categories and Bear IDs are factors
final2$Sex <- as.factor(final2$Sex)
final2$Day_periods <- as.factor(final2$Day_periods)
final2$Bear <- as.factor(final2$Bear)
final2$aspen <- as.factor(final2$aspen)
final2$barren <- as.factor(final2$barren)
final2$conifer <- as.factor(final2$conifer)
final2$shrub <- as.factor(final2$shrub)
final2$grass <- as.factor(final2$grass)
final2$Season <- as.factor(final2$Season)
final2$Day_periods <- as.factor(final2$Day_periods)



# create trajectory and home ranges
if (nrow(final2) > 0) {
  
  # Create separate plots for each bear
  unique_bears <- unique(final2$Bear)
  
  for (bear_id in unique_bears) {
    bear_data <- final2[final2$Bear == bear_id, ]
    
    # Create SpatialPoints
    bear_coords <- SpatialPoints(bear_data[, c("easting", "northing")], proj4string = CRS(as.character(NA)))
    
    # Set up smaller margins
    par(mar = c(3, 3, 1, 1))  # Adjust the values as needed
    
    # Plot the trajectory
    plot(bear_coords, pch = 16, col = "blue", main = paste("Bear Trajectory - ID", bear_id))
    
    # Calculate KDE for home range
    kde_result <- kernelUD(bear_coords, h = "href")
    
    # Plot home range
    plot(kde_result, col = "red", main = paste("Home Range - Bear ID", bear_id))
  }
  
} else {
  print("Data frame is empty.")
}




head(final2)
dim(final2)

names(final2)

# How many unique bears are in the dataset?
length(unique(final2$Bear)) # 23

# what is the earliest and latest date in the dataset?
min(final2$FormattedDateTime) # 2000-09-16 19:30:00
max(final2$FormattedDateTime) # 2022-11-12 22:01:00

# calculate the time span of the dataset
max(final2$FormattedDateTime) - min(final2$FormattedDateTime) # Time difference of 8033.52 days




head(final2)
dim(final2)
levels(final2$Season)

final2$FormattedDateTime[1:100]

subset(final2, Season == "Fall")



# Create a multivariate model to model a resource selection function (RSF) for aspen, barren, conifer, grass, and shrub
# use Season, Greeness as fixed effects and Bear as a random effect
# l1 <- bf(X1 | mi() ~ 0, family = gaussian())

# make a PCA 
pcaStru <- prcomp(final2[, c("distance_facility", "distance_road", "distance_trail")], scale = TRUE)

# Extract the loadings
pcaStru$rotation


# Extract explained variance
pcaStru$sdev^2 / sum(pcaStru$sdev^2)

# get principal components for habitat quality
# make a PCA 
pcaHQ <- prcomp(final2[, c("distance_stream", "distance_river", "greenness")], scale = TRUE)

# Extract the loadings
pcaHQ$rotation
# Extract explained variance
pcaHQ$sdev^2 / sum(pcaHQ$sdev^2)


# add the principal components to the data frame
final2$Str1 <- pcaStru$x[, 1]
final2$Str2 <- pcaStru$x[, 2]

final2$HQ1 <- pcaHQ$x[, 1]
final2$HQ2 <- pcaHQ$x[, 2]

# try with one individual
bear1 <- final2[final2$Bear == unique(final2$Bear)[1], ]

final2$df2 <- final2$distance_facility^2

head(final2)

# aggregate the data to get the mean of the variables for each seven days interval

# Create a new column for the week number
final2$week_number <- as.numeric(format(final2$FormattedDateTime, "%V"))
# get year from date
final2$year <- as.numeric(format(final2$FormattedDateTime, "%Y"))


# convert aspen, barren, conifer, grass, and shrub to integer
final2$aspen <- as.integer(final2$aspen)
final2$barren <- as.integer(final2$barren)
final2$conifer <- as.integer(final2$conifer)
final2$grass <- as.integer(final2$grass)
final2$shrub <- as.integer(final2$shrub)

mean(final2$greenness)
sd(final2$greenness)



# Aggregate the data by year,  week number, bear, sex, and season; and sum aspen, barren, conifer, grass, and shrub
final2_agg <- final2 %>%
  group_by(year, week_number, Bear, Sex) %>%
  summarise(
    aspen = sum(aspen),
    barren = sum(barren),
    conifer = sum(conifer),
    grass = sum(grass),
    shrub = sum(shrub),
    facility = mean(log10(distance_facility)),
    road = mean(log10(distance_road)),
    trail = mean(log10(distance_trail)),
    stream = mean(log10(distance_stream)),
    river = mean(log10(distance_river)),
    greenness = mean(log10(greenness+1)),
    
    
    # now sd for each habitat
    facility_sd = sd(log10(distance_facility)),
    road_sd = sd(log10(distance_road)),
    trail_sd = sd(log10(distance_trail)),
    stream_sd = sd(log10(distance_stream)),
    river_sd = sd(log10(distance_river)),
    greenness_sd = sd(log10(greenness+1), na.rm = TRUE),
  )

head(final2_agg) 
dim(final2_agg)
unique(final2_agg$Bear)

final2_agg$sex <- factor(final2_agg$Sex)

levels(final2_agg$Sex) <- c(0,1)

# remove NAs from the data frame with zero
final2_agg[is.na(final2_agg)] <- 0

# save as csv
write.csv(final2_agg, file = "final2_agg.csv", row.names = FALSE)

N <- 15

# dat <- data.frame(
#     y1 = rbinom(N, 10, 0.3), y2 = rbinom(N, 10, 0.5), 
#     y3 = rbinom(N, 10, 0.7), x = rnorm(N)
# )

# dat$size <- with(dat, y1 + y2 + y3)
# dat$y <- with(dat, cbind(y1, y2, y3))



dat <- final2
dat$size <- with(dat, aspen + barren + conifer + grass + shrub)
dat$y <- with(dat, cbind(aspen, barren, conifer, grass, shrub))

# prior <- prior(normal(0, 10), "b", dpar = muy2) +
#     prior(cauchy(0, 1), "Intercept") +
#     prior(normal(0, 2), "Intercept", dpar = muy3)

names(dat)

# NORMALIZE THE PREDICTORS
names(dat)
dat$facility <- (dat$distance_facility - mean(dat$distance_facility)) / sd(dat$distance_facility)
dat$road <- (dat$distance_road - mean(dat$distance_road)) / sd(dat$distance_road)
dat$trail <- (dat$distance_trail - mean(dat$distance_trail)) / sd(dat$distance_trail)
dat$stream <- (dat$distance_stream - mean(dat$distance_stream)) / sd(dat$distance_stream)
dat$river <- (dat$distance_river - mean(dat$distance_river)) / sd(dat$distance_river)
dat$greenness <- (dat$greenness - mean(dat$greenness)) / sd(dat$greenness)


# check ranges

range(dat$distance_facility)
range(dat$facility)

mean(dat$distance_facility)
sd(dat$distance_facility)

mean(dat$distance_facility) + sd(dat$distance_facility) * 6.269179
mean(dat$distance_facility) + sd(dat$distance_facility) * 0
mean(dat$distance_facility) + sd(dat$distance_facility) * 1


range(dat$greenness) #-72.7711, 89.9671
mean(dat$greenness) #27.57094
sd(dat$greenness) #13.19787
range(dat$greenness) # -7.602895  4.727744







# set plot dimensions to 1 x 1
par(mfrow = c(1, 1))

plot(dm, sm, type = "l", xlab = "distance_facility", ylab = "standardized distance_facility")

# fit <- brm(bf(y | trials(size) ~ facility + road + trail + stream + river + greenness + Sex + Season + (1|Bear)), data = dat, 
#                 family = multinomial())

# summary(fit)

# saveRDS(fit, "fit_season.rds")

fit = readRDS("fit_season.rds")

summary(fit)

# so, the significant predictors are facility, road, trail, stream, river indicate whether the probability of the bear selecting a habitat is different from the reference level (0)
# the values are in a logit scale, so we need to transform them to get the probabilities



# # plot predictions

conditional_effects(fit, conditions = "Season:facility", categorical = TRUE, ndraws = 100)







#change 0 to sm and vice versa for each covariate plot (facility, road, etc.) seperately
#Facility: Male x Season
        
        # check ranges
        range(dat$distance_facility)
        range(dat$facility)
        
        mean(dat$distance_facility)
        sd(dat$distance_facility)
        
        mean(dat$distance_facility) + sd(dat$distance_facility) * 6.269179
        mean(dat$distance_facility) + sd(dat$distance_facility) * 0
        mean(dat$distance_facility) + sd(dat$distance_facility) * -1.338  
        
        #determine xlim
        sm = seq(from=-2, to=7, by=0.1) 
        
      
        
        # Get mean and standard deviation of greenness from your dataset
        mean_greenness <- mean(dat$greenness, na.rm = TRUE)
        sd_greenness <- sd(dat$greenness, na.rm = TRUE)
        
        # Define greenness levels as mean, -2 SD, and +2 SD
        greenness_levels <- c(mean_greenness - 2 * sd_greenness, mean_greenness, mean_greenness + 2 * sd_greenness)
        
        # Function to create data frames with adjusted greenness levels
        create_new_data <- function(season, sex) {
          expand.grid(
            facility = sm,
            road = 0,
            trail = 0,
            stream = 0,
            river = 0,
            greenness = greenness_levels,
            Sex = sex,
            size = 1,
            Season = season
          )
        }
        
        # Create data frames for each season
        new_MF <- create_new_data("Fall", "M")
        new_MSP <- create_new_data("Spring", "M")
        new_MSU <- create_new_data("Summer", "M")
        
        # Generate posterior predictions for each dataset
        pdf_MF <- posterior_epred(fit, newdata = new_MF, re_formula = NA, ndraws = 1000)
        pdf_MSP <- posterior_epred(fit, newdata = new_MSP, re_formula = NA, ndraws = 1000)
        pdf_MSU <- posterior_epred(fit, newdata = new_MSU, re_formula = NA, ndraws = 1000)
        
        # Function to calculate means and confidence intervals
        calc_means_cis <- function(pdf) {
          list(
            ciAspen = t(hdi(pdf[, , 1])),
            ciBarren = t(hdi(pdf[, , 2])),
            ciConifer = t(hdi(pdf[, , 3])),
            ciGrass = t(hdi(pdf[, , 4])),
            ciShrub = t(hdi(pdf[, , 5])),
            meanAspen = apply(pdf[, , 1], 2, mean),
            meanBarren = apply(pdf[, , 2], 2, mean),
            meanConifer = apply(pdf[, , 3], 2, mean),
            meanGrass = apply(pdf[, , 4], 2, mean),
            meanShrub = apply(pdf[, , 5], 2, mean)
          )
        }
        
        means_cis_MF <- calc_means_cis(pdf_MF)
        means_cis_MSP <- calc_means_cis(pdf_MSP)
        means_cis_MSU <- calc_means_cis(pdf_MSU)
        
        # Create data frames for plotting
        create_df <- function(sm, means_cis, season) {
          data.frame(
            x = sm,
            Aspen = means_cis$meanAspen,
            lAspen = means_cis$ciAspen[, 1],
            upAspen = means_cis$ciAspen[, 2],
            Barren = means_cis$meanBarren,
            lBarren = means_cis$ciBarren[, 1],
            upBarren = means_cis$ciBarren[, 2],
            Conifer = means_cis$meanConifer,
            lConifer = means_cis$ciConifer[, 1],
            upConifer = means_cis$ciConifer[, 2],
            Grass = means_cis$meanGrass,
            lGrass = means_cis$ciGrass[, 1],
            upGrass = means_cis$ciGrass[, 2],
            Shrub = means_cis$meanShrub,
            lShrub = means_cis$ciShrub[, 1],
            upShrub = means_cis$ciShrub[, 2],
            Sex = "M",
            Season = season,
            greenness = rep(greenness_levels, each = length(sm))
          )
        }
        
        df_MF <- create_df(sm, means_cis_MF, "Fall")
        df_MSP <- create_df(sm, means_cis_MSP, "Spring")
        df_MSU <- create_df(sm, means_cis_MSU, "Summer")
        
        # Combine all data frames for plotting
        df <- bind_rows(df_MF, df_MSP, df_MSU)
        
        # Set greenness factor levels for correct ordering
        df$greenness <- factor(df$greenness, levels = greenness_levels, labels = c("Low", "Mean", "High"))
        
        # Define colors for the plots
        colors <- c("Aspen" = "limegreen",
                    "Barren" = "brown",
                    "Conifer" = "yellow2",
                    "Grass" = "magenta",
                    "Shrub" = "coral")
        
        # Plot for Aspen + Conifer with Title
        ggplot(data = df, aes(x = x)) +
          geom_ribbon(aes(ymin = lAspen, ymax = upAspen, fill = "Aspen"), alpha = 0.3) +
          geom_ribbon(aes(ymin = lConifer, ymax = upConifer, fill = "Conifer"), alpha = 0.3) +
          geom_line(aes(y = Aspen, color = "Aspen"), size = 1.5) +
          geom_line(aes(y = Conifer, color = "Conifer"), size = 1.5) +
          facet_grid(greenness ~ Season, labeller = labeller(Sex = c(M = "Male"))) +
          xlab("Distance To Facilities (m)") +
          ylab("Probability Of Occurrence (%)") +
          scale_color_manual(values = c("Aspen" = "limegreen", "Conifer" = "yellow2"), name = "Habitat") +
          scale_fill_manual(values = c("Aspen" = "limegreen", "Conifer" = "yellow2"), name = "Habitat") +
          scale_x_continuous(breaks = c(-1.338257, 0.661743, 2.661743, 4.661743, 6.61743),
                             labels = c(0, 2687, 5373, 8060, 10747)) +
          scale_y_continuous(breaks = c(0.15, 0.20, 0.25),
                             labels = c(15, 20, 25)) +
          ggtitle("Male - Facilities with Greenness Variations") +
          theme_bw() +
          coord_cartesian(xlim = c(-1.338257, 6.61743), expand = FALSE) +
          theme(panel.spacing = unit(2.5, "lines"),
                axis.text.x = element_text(size = 10, margin = margin(b = 5)),
                axis.text.y = element_text(size = 10, margin = margin(r = 5)),
                axis.title = element_text(size = 12),
                strip.text = element_text(size = 12),
                legend.position = "right",
                legend.title = element_text(size = 12),
                legend.text = element_text(size = 10))
        
        # Plot for Barren + Grass + Shrub with Title
        ggplot(data = df, aes(x = x)) +
          geom_ribbon(aes(ymin = lBarren, ymax = upBarren, fill = "Barren"), alpha = 0.3) +
          geom_ribbon(aes(ymin = lGrass, ymax = upGrass, fill = "Grass"), alpha = 0.3) +
          geom_ribbon(aes(ymin = lShrub, ymax = upShrub, fill = "Shrub"), alpha = 0.3) +
          geom_line(aes(y = Barren, color = "Barren"), size = 1.5) +
          geom_line(aes(y = Grass, color = "Grass"), size = 1.5) +
          geom_line(aes(y = Shrub, color = "Shrub"), size = 1.5) +
          facet_grid(greenness ~ Season, labeller = labeller(Sex = c(M = "Male"))) +
          xlab("Distance To Facilities (m)") +
          ylab("Probability Of Occurrence (%)") +
          scale_color_manual(values = c("Barren" = "brown", "Grass" = "magenta", "Shrub" = "coral"), name = "Habitat") +
          scale_fill_manual(values = c("Barren" = "brown", "Grass" = "magenta", "Shrub" = "coral"), name = "Habitat") +
          scale_x_continuous(breaks = c(-1.338257, 0.661743, 2.661743, 4.661743, 6.61743),
                             labels = c(0, 2687, 5373, 8060, 10747)) +
          scale_y_continuous(breaks = c(0.15, 0.20, 0.25),
                             labels = c(15, 20, 25)) +
          ggtitle("Male - Facilities with Greenness Variations") +
          theme_bw() +
          coord_cartesian(xlim = c(-1.338257, 6.61743), expand = FALSE) +
          theme(panel.spacing = unit(2.5, "lines"),
                axis.text.x = element_text(size = 10, margin = margin(b = 5)),
                axis.text.y = element_text(size = 10, margin = margin(r = 5)),
                axis.title = element_text(size = 12),
                strip.text = element_text(size = 12),
                legend.position = "right",
                legend.title = element_text(size = 12),
                legend.text = element_text(size = 10))
        
  

  
  #overall probabilities between curves          
  pdf_MF[, , 2] # barren
  
  dAB = pdf_MSP[, , 4] - pdf_MSU[, , 4]
  
  
  LOS <- function(x) {
    return(100 * length(which(x > 0))/length(x))
  }
  
  
  PP <- apply(dAB, 2, LOS)
  LOS(dAB) #overall probability of different curves between habitat types, if the number's close to 50, they're similar
  #closer to 0 -> second curve is larger and vice versa
  
  
  
  
  
  
    
  
  
  
  
#Road: Male x Season
        # check ranges
        range(dat$distance_road)
        range(dat$road)
        
        mean(dat$distance_road)
        sd(dat$distance_road)
        
        mean(dat$distance_road) + sd(dat$distance_road) * 6.5549247
        mean(dat$distance_road) + sd(dat$distance_road) * 0
        mean(dat$distance_road) + sd(dat$distance_road) * -0.867
        
        
        # Define sequence for 'sm' based on the road data
        sm <- seq(from = -2, to = 7.5, by = 0.1) # Adjust based on your specific data ranges
        
        # Get mean and standard deviation of greenness from your dataset
        mean_greenness <- mean(dat$greenness, na.rm = TRUE)
        sd_greenness <- sd(dat$greenness, na.rm = TRUE)
        
        # Define greenness levels as mean, -2 SD, and +2 SD
        greenness_levels <- c(mean_greenness - 2 * sd_greenness, mean_greenness, mean_greenness + 2 * sd_greenness)
        
        # Function to create data frames with adjusted greenness levels for road
        create_new_data <- function(season, greenness) {
          expand.grid(
            facility = 0,
            road = sm, # Road is now varied
            trail = 0,
            stream = 0,
            river = 0,
            greenness = greenness, 
            Sex = "M",
            size = 1, 
            Season = season
          )
        }
        
        # Create data frames for each season and greenness level
        new_MRF <- bind_rows(
          create_new_data("Fall", greenness_levels[1]),
          create_new_data("Fall", greenness_levels[2]),
          create_new_data("Fall", greenness_levels[3])
        )
        
        new_MRSP <- bind_rows(
          create_new_data("Spring", greenness_levels[1]),
          create_new_data("Spring", greenness_levels[2]),
          create_new_data("Spring", greenness_levels[3])
        )
        
        new_MRSU <- bind_rows(
          create_new_data("Summer", greenness_levels[1]),
          create_new_data("Summer", greenness_levels[2]),
          create_new_data("Summer", greenness_levels[3])
        )
        
        # Generate posterior predictions for each dataset
        pdf_MRF <- posterior_epred(fit, newdata = new_MRF, re_formula = NA, ndraws = 1000)
        pdf_MRSP <- posterior_epred(fit, newdata = new_MRSP, re_formula = NA, ndraws = 1000)
        pdf_MRSU <- posterior_epred(fit, newdata = new_MRSU, re_formula = NA, ndraws = 1000)
        
        # Calculate means and confidence intervals for each greenness level
        means_cis_MRF <- calc_means_cis(pdf_MRF)
        means_cis_MRSP <- calc_means_cis(pdf_MRSP)
        means_cis_MRSU <- calc_means_cis(pdf_MRSU)
        
        # Function to create data frame for plotting
        create_df <- function(sm, means_cis, season, greenness) {
          data.frame(
            x = sm,
            Aspen = means_cis$meanAspen,
            lAspen = means_cis$ciAspen[, 1],
            upAspen = means_cis$ciAspen[, 2],
            Barren = means_cis$meanBarren,
            lBarren = means_cis$ciBarren[, 1],
            upBarren = means_cis$ciBarren[, 2],
            Conifer = means_cis$meanConifer,
            lConifer = means_cis$ciConifer[, 1],
            upConifer = means_cis$ciConifer[, 2],
            Grass = means_cis$meanGrass,
            lGrass = means_cis$ciGrass[, 1],
            upGrass = means_cis$ciGrass[, 2],
            Shrub = means_cis$meanShrub,
            lShrub = means_cis$ciShrub[, 1],
            upShrub = means_cis$ciShrub[, 2],
            Sex = "M",
            Season = season,
            Greenness = greenness
          )
        }
        
        # Create data frames for plotting with the specified greenness levels
        df_MRF <- create_df(sm, means_cis_MRF, "Fall", rep(greenness_levels, each = length(sm)))
        df_MRSP <- create_df(sm, means_cis_MRSP, "Spring", rep(greenness_levels, each = length(sm)))
        df_MRSU <- create_df(sm, means_cis_MRSU, "Summer", rep(greenness_levels, each = length(sm)))
        
        # Combine all data frames for plotting
        df <- bind_rows(df_MRF, df_MRSP, df_MRSU)
        df$Season <- factor(df$Season, levels = c("Spring", "Summer", "Fall"))
        df$Greenness <- factor(df$Greenness, levels = greenness_levels, labels = c("Low", "Mean", "High"))
        
        # Define colors for the plots
        colors <- c("Aspen" = "limegreen",
                    "Barren" = "brown",
                    "Conifer" = "yellow2",
                    "Grass" = "magenta",
                    "Shrub" = "coral")
        
        # Plot for Aspen + Conifer
        ggplot(data = df, aes(x = x)) +
          geom_ribbon(aes(ymin = lAspen, ymax = upAspen, fill = "Aspen"), alpha = 0.3) +
          geom_ribbon(aes(ymin = lConifer, ymax = upConifer, fill = "Conifer"), alpha = 0.3) +
          geom_line(aes(y = Aspen, color = "Aspen"), size = 1.5) +
          geom_line(aes(y = Conifer, color = "Conifer"), size = 1.5) +
          facet_grid(Greenness ~ Season, labeller = labeller(Sex = c(M = "Male"))) +
          xlab("Distance To Roads (m)") +
          ylab("Probability Of Occurrence (%)") +
          scale_color_manual(values = c("Aspen" = "limegreen", "Conifer" = "yellow2"), name = "Habitat") +
          scale_fill_manual(values = c("Aspen" = "limegreen", "Conifer" = "yellow2"), name = "Habitat") +
          scale_x_continuous(breaks = c(-0.867, 1.133, 3.133, 5.133, 7.133),
                             labels = c(0, 2497, 4992, 7487, 9983)) +
          scale_y_continuous(breaks = c(0.15, 0.20, 0.25, 0.30, 0.35),
                             labels = c(15, 20, 25, 30, 35)) +
          ggtitle("Male - Roads with Greenness Variations") +
          theme_bw() +
          coord_cartesian(xlim = c(-0.867, 7.133), expand = FALSE) + 
          theme(panel.spacing.x = unit(1.5, "lines"),
                axis.text.x = element_text(margin = margin(b = 2)),
                axis.title.y = element_text(margin = margin(t = 0, r = 5, b = 0, l = 0)))
        
        # Plot for Barren + Grass + Shrub
        ggplot(data = df, aes(x = x)) +
          geom_ribbon(aes(ymin = lBarren, ymax = upBarren, fill = "Barren"), alpha = 0.3) +
          geom_ribbon(aes(ymin = lGrass, ymax = upGrass, fill = "Grass"), alpha = 0.3) +
          geom_ribbon(aes(ymin = lShrub, ymax = upShrub, fill = "Shrub"), alpha = 0.3) +
          geom_line(aes(y = Barren, color = "Barren"), size = 1.5) +
          geom_line(aes(y = Grass, color = "Grass"), size = 1.5) +
          geom_line(aes(y = Shrub, color = "Shrub"), size = 1.5) +
          facet_grid(Greenness ~ Season, labeller = labeller(Sex = c(M = "Male"))) +
          xlab("Distance To Roads (m)") +
          ylab("Probability Of Occurrence (%)") +
          scale_color_manual(values = c("Barren" = "brown", "Grass" = "magenta", "Shrub" = "coral"), name = "Habitat") +
          scale_fill_manual(values = c("Barren" = "brown", "Grass" = "magenta", "Shrub" = "coral"), name = "Habitat") +
          scale_x_continuous(breaks = c(-0.867, 1.133, 3.133, 5.133, 7.133),
                             labels = c(0, 2497, 4992, 7487, 9983)) +
          scale_y_continuous(breaks = c(0.15, 0.20, 0.25, 0.30, 0.35),
                             labels = c(15, 20, 25, 30, 35)) +
          ggtitle("Male - Roads with Greenness Variations") +
          theme_bw() +
          coord_cartesian(xlim = c(-0.867, 7.133), expand = FALSE) + 
          theme(panel.spacing.x = unit(1.5, "lines"),
                axis.text.x = element_text(margin = margin(b = 2)),
                axis.title.y = element_text(margin = margin(t = 0, r = 5, b = 0, l = 0)))
        
  
    
  
  #overall probabilities between curves          
    pdf_FF[, , 2] # barren
    
    dAB = pdf_MRF[, , 3] - pdf_MRSP[, , 3]
    
    LOS <- function(x) {
      return(100 * length(which(x > 0))/length(x))
    }
    
    
    PP <- apply(dAB, 2, LOS)
    LOS(dAB) #overall probability of different curves between habitat types, if the number's close to 50, they're similar
    #closer to 0 -> second curve is larger and vice versa
    
    
    
    
    
  
            
#Trail: Male x Season
    
        # Get mean and standard deviation of greenness from your dataset
        mean_greenness <- mean(dat$greenness, na.rm = TRUE)
        sd_greenness <- sd(dat$greenness, na.rm = TRUE)
        
        # Define greenness levels as mean, -2 SD, and +2 SD
        greenness_levels <- c(mean_greenness - 2 * sd_greenness, mean_greenness, mean_greenness + 2 * sd_greenness)
        
        # Function to create data frames with adjusted greenness levels for trail
        create_new_data <- function(season, greenness) {
          expand.grid(
            facility = 0,
            road = 0,
            trail = sm, # Trail is now varied
            stream = 0,
            river = 0,
            greenness = greenness, 
            Sex = "M",
            size = 1, 
            Season = season
          )
        }
        
        # Create data frames for each season and greenness level
        new_MTF <- bind_rows(
          create_new_data("Fall", greenness_levels[1]),
          create_new_data("Fall", greenness_levels[2]),
          create_new_data("Fall", greenness_levels[3])
        )
        
        new_MTSP <- bind_rows(
          create_new_data("Spring", greenness_levels[1]),
          create_new_data("Spring", greenness_levels[2]),
          create_new_data("Spring", greenness_levels[3])
        )
        
        new_MTSU <- bind_rows(
          create_new_data("Summer", greenness_levels[1]),
          create_new_data("Summer", greenness_levels[2]),
          create_new_data("Summer", greenness_levels[3])
        )
        
        # Generate posterior predictions for each dataset
        pdf_MTF <- posterior_epred(fit, newdata = new_MTF, re_formula = NA, ndraws = 1000)
        pdf_MTSP <- posterior_epred(fit, newdata = new_MTSP, re_formula = NA, ndraws = 1000)
        pdf_MTSU <- posterior_epred(fit, newdata = new_MTSU, re_formula = NA, ndraws = 1000)
        
        # Function to calculate means and confidence intervals
        calc_means_cis <- function(pdf) {
          list(
            ciAspen = t(hdi(pdf[, , 1])),
            ciBarren = t(hdi(pdf[, , 2])),
            ciConifer = t(hdi(pdf[, , 3])),
            ciGrass = t(hdi(pdf[, , 4])),
            ciShrub = t(hdi(pdf[, , 5])),
            meanAspen = apply(pdf[, , 1], 2, mean),
            meanBarren = apply(pdf[, , 2], 2, mean),
            meanConifer = apply(pdf[, , 3], 2, mean),
            meanGrass = apply(pdf[, , 4], 2, mean),
            meanShrub = apply(pdf[, , 5], 2, mean)
          )
        }
        
        means_cis_MTF <- calc_means_cis(pdf_MTF)
        means_cis_MTSP <- calc_means_cis(pdf_MTSP)
        means_cis_MTSU <- calc_means_cis(pdf_MTSU)
        
        # Function to create data frame for plotting
        create_df <- function(sm, means_cis, season, greenness) {
          data.frame(
            x = sm,
            Aspen = means_cis$meanAspen,
            lAspen = means_cis$ciAspen[, 1],
            upAspen = means_cis$ciAspen[, 2],
            Barren = means_cis$meanBarren,
            lBarren = means_cis$ciBarren[, 1],
            upBarren = means_cis$ciBarren[, 2],
            Conifer = means_cis$meanConifer,
            lConifer = means_cis$ciConifer[, 1],
            upConifer = means_cis$ciConifer[, 2],
            Grass = means_cis$meanGrass,
            lGrass = means_cis$ciGrass[, 1],
            upGrass = means_cis$ciGrass[, 2],
            Shrub = means_cis$meanShrub,
            lShrub = means_cis$ciShrub[, 1],
            upShrub = means_cis$ciShrub[, 2],
            Sex = "M",
            Season = season,
            Greenness = greenness
          )
        }
        
        # Create data frames for plotting with the specified greenness levels
        df_MTF <- create_df(sm, means_cis_MTF, "Fall", rep(greenness_levels, each = length(sm)))
        df_MTSP <- create_df(sm, means_cis_MTSP, "Spring", rep(greenness_levels, each = length(sm)))
        df_MTSU <- create_df(sm, means_cis_MTSU, "Summer", rep(greenness_levels, each = length(sm)))
        
        # Combine all data frames for plotting
        df <- bind_rows(df_MTF, df_MTSP, df_MTSU)
        df$Season <- factor(df$Season, levels = c("Spring", "Summer", "Fall"))
        df$Greenness <- factor(df$Greenness, levels = greenness_levels, labels = c("Low", "Mean", "High"))
        
        # Define colors for the plots
        colors <- c("Aspen" = "limegreen",
                    "Barren" = "brown",
                    "Conifer" = "yellow2",
                    "Grass" = "magenta",
                    "Shrub" = "coral")
        
        # Plot for Aspen + Conifer
        ggplot(data = df, aes(x = x)) +
          geom_ribbon(aes(ymin = lAspen, ymax = upAspen, fill = "Aspen"), alpha = 0.3) +
          geom_ribbon(aes(ymin = lConifer, ymax = upConifer, fill = "Conifer"), alpha = 0.3) +
          geom_line(aes(y = Aspen, color = "Aspen"), size = 1.5) +
          geom_line(aes(y = Conifer, color = "Conifer"), size = 1.5) +
          facet_grid(Greenness ~ Season, labeller = labeller(Sex = c(M = "Male"))) +
          xlab("Distance To Trails (m)") +
          ylab("Probability Of Occurrence (%)") +
          scale_color_manual(values = c("Aspen" = "limegreen", "Conifer" = "yellow2"), name = "Habitat") +
          scale_fill_manual(values = c("Aspen" = "limegreen", "Conifer" = "yellow2"), name = "Habitat") +
          scale_x_continuous(breaks = c(-0.9726994, 1.027301, 3.027301, 5.027301, 7.027301),
                             labels = c(0, 1422, 2843, 4265, 5686)) +
          scale_y_continuous(breaks = c(0.15, 0.20, 0.25, 0.30),
                             labels = c(15, 20, 25, 30)) +
          ggtitle("Male - Trails with Greenness Variations") +
          theme_bw() +
          coord_cartesian(xlim = c(-0.9726994, 7.027301), expand = FALSE) + 
          theme(panel.spacing.x = unit(1.5, "lines"),
                axis.text.x = element_text(margin = margin(b = 2)),
                axis.title.y = element_text(margin = margin(t = 0, r = 5, b = 0, l = 0)))
        
        # Plot for Barren + Grass + Shrub
        ggplot(data = df, aes(x = x)) +
          geom_ribbon(aes(ymin = lBarren, ymax = upBarren, fill = "Barren"), alpha = 0.3) +
          geom_ribbon(aes(ymin = lGrass, ymax = upGrass, fill = "Grass"), alpha = 0.3) +
          geom_ribbon(aes(ymin = lShrub, ymax = upShrub, fill = "Shrub"), alpha = 0.3) +
          geom_line(aes(y = Barren, color = "Barren"), size = 1.5) +
          geom_line(aes(y = Grass, color = "Grass"), size = 1.5) +
          geom_line(aes(y = Shrub, color = "Shrub"), size = 1.5) +
          facet_grid(Greenness ~ Season, labeller = labeller(Sex = c(M = "Male"))) +
          xlab("Distance To Trails (m)") +
          ylab("Probability Of Occurrence (%)") +
          scale_color_manual(values = c("Barren" = "brown", "Grass" = "magenta", "Shrub" = "coral"), name = "Habitat") +
          scale_fill_manual(values = c("Barren" = "brown", "Grass" = "magenta", "Shrub" = "coral"), name = "Habitat") +
          scale_x_continuous(breaks = c(-0.9726994, 1.027301, 3.027301, 5.027301, 7.027301),
                             labels = c(0, 1422, 2843, 4265, 5686)) +
          scale_y_continuous(breaks = c(0.15, 0.20, 0.25, 0.30),
                             labels = c(15, 20, 25, 30)) +
          ggtitle("Male - Trails with Greenness Variations") +
          theme_bw() +
          coord_cartesian(xlim = c(-0.9726994, 7.027301), expand = FALSE) + 
          theme(panel.spacing.x = unit(1.5, "lines"),
                axis.text.x = element_text(margin = margin(b = 2)),
                axis.title.y = element_text(margin = margin(t = 0, r = 5, b = 0, l = 0)))
                
            
  #         
  #
  #
  #
  #
  #
  #
  #
  #
  #
  #
  #
  #
  #
            
#Facility: Female x Season
        
        # Define sequence for 'sm' based on the facility data
        sm <- seq(from = -2, to = 7, by = 0.1)
        
        # Get mean and standard deviation of greenness from your dataset
        mean_greenness <- mean(dat$greenness, na.rm = TRUE)
        sd_greenness <- sd(dat$greenness, na.rm = TRUE)
        
        # Define greenness levels as mean, -2 SD, and +2 SD
        greenness_levels <- c(mean_greenness - 2 * sd_greenness, mean_greenness, mean_greenness + 2 * sd_greenness)
        
        # Function to create data frames with adjusted greenness levels for facilities
        create_new_data <- function(season, greenness) {
          expand.grid(
            facility = sm,
            road = 0,
            trail = 0,
            stream = 0,
            river = 0,
            greenness = greenness, 
            Sex = "F",
            size = 1, 
            Season = season
          )
        }
        
        # Create data frames for each season and greenness level
        new_FF <- bind_rows(
          create_new_data("Fall", greenness_levels[1]),
          create_new_data("Fall", greenness_levels[2]),
          create_new_data("Fall", greenness_levels[3])
        )
        
        new_FSP <- bind_rows(
          create_new_data("Spring", greenness_levels[1]),
          create_new_data("Spring", greenness_levels[2]),
          create_new_data("Spring", greenness_levels[3])
        )
        
        new_FSU <- bind_rows(
          create_new_data("Summer", greenness_levels[1]),
          create_new_data("Summer", greenness_levels[2]),
          create_new_data("Summer", greenness_levels[3])
        )
        
        # Generate posterior predictions for each dataset
        pdf_FF <- posterior_epred(fit, newdata = new_FF, re_formula = NA, ndraws = 1000)
        pdf_FSP <- posterior_epred(fit, newdata = new_FSP, re_formula = NA, ndraws = 1000)
        pdf_FSU <- posterior_epred(fit, newdata = new_FSU, re_formula = NA, ndraws = 1000)
        
        # Function to calculate means and confidence intervals
        calc_means_cis <- function(pdf) {
          list(
            ciAspen = t(hdi(pdf[, , 1])),
            ciBarren = t(hdi(pdf[, , 2])),
            ciConifer = t(hdi(pdf[, , 3])),
            ciGrass = t(hdi(pdf[, , 4])),
            ciShrub = t(hdi(pdf[, , 5])),
            meanAspen = apply(pdf[, , 1], 2, mean),
            meanBarren = apply(pdf[, , 2], 2, mean),
            meanConifer = apply(pdf[, , 3], 2, mean),
            meanGrass = apply(pdf[, , 4], 2, mean),
            meanShrub = apply(pdf[, , 5], 2, mean)
          )
        }
        
        means_cis_FF <- calc_means_cis(pdf_FF)
        means_cis_FSP <- calc_means_cis(pdf_FSP)
        means_cis_FSU <- calc_means_cis(pdf_FSU)
        
        # Function to create data frame for plotting
        create_df <- function(sm, means_cis, season, greenness) {
          data.frame(
            x = sm,
            Aspen = means_cis$meanAspen,
            lAspen = means_cis$ciAspen[, 1],
            upAspen = means_cis$ciAspen[, 2],
            Barren = means_cis$meanBarren,
            lBarren = means_cis$ciBarren[, 1],
            upBarren = means_cis$ciBarren[, 2],
            Conifer = means_cis$meanConifer,
            lConifer = means_cis$ciConifer[, 1],
            upConifer = means_cis$ciConifer[, 2],
            Grass = means_cis$meanGrass,
            lGrass = means_cis$ciGrass[, 1],
            upGrass = means_cis$ciGrass[, 2],
            Shrub = means_cis$meanShrub,
            lShrub = means_cis$ciShrub[, 1],
            upShrub = means_cis$ciShrub[, 2],
            Sex = "F",
            Season = season,
            Greenness = greenness
          )
        }
        
        # Create data frames for plotting with the specified greenness levels
        df_FF <- create_df(sm, means_cis_FF, "Fall", rep(greenness_levels, each = length(sm)))
        df_FSP <- create_df(sm, means_cis_FSP, "Spring", rep(greenness_levels, each = length(sm)))
        df_FSU <- create_df(sm, means_cis_FSU, "Summer", rep(greenness_levels, each = length(sm)))
        
        # Combine all data frames for plotting
        df <- bind_rows(df_FF, df_FSP, df_FSU)
        df$Season <- factor(df$Season, levels = c("Spring", "Summer", "Fall"))
        df$Greenness <- factor(df$Greenness, levels = greenness_levels, labels = c("Low", "Mean", "High"))
        
        # Define colors for the plots
        colors <- c("Aspen" = "limegreen",
                    "Barren" = "brown",
                    "Conifer" = "yellow2",
                    "Grass" = "magenta",
                    "Shrub" = "coral")
        
        # Plot for Aspen + Conifer
        ggplot(data = df, aes(x = x)) +
          geom_ribbon(aes(ymin = lAspen, ymax = upAspen, fill = "Aspen"), alpha = 0.3) +
          geom_ribbon(aes(ymin = lConifer, ymax = upConifer, fill = "Conifer"), alpha = 0.3) +
          geom_line(aes(y = Aspen, color = "Aspen"), size = 1.5) +
          geom_line(aes(y = Conifer, color = "Conifer"), size = 1.5) +
          facet_grid(Greenness ~ Season, labeller = labeller(Sex = c(F = "Female"))) +
          xlab("Distance To Facilities (m)") +
          ylab("Probability Of Occurrence (%)") +
          scale_color_manual(values = c("Aspen" = "limegreen", "Conifer" = "yellow2"), name = "Habitat") +
          scale_fill_manual(values = c("Aspen" = "limegreen", "Conifer" = "yellow2"), name = "Habitat") +
          scale_x_continuous(breaks = c(-1.338257, 0.661743, 2.661743, 4.661743, 6.61743),
                             labels = c(0, 2687, 5373, 8060, 10747)) +
          scale_y_continuous(breaks = c(0.15, 0.20, 0.25),
                             labels = c(15, 20, 25)) +
          ggtitle("Female - Facilities with Greenness Variations") +
          theme_bw() +
          coord_cartesian(xlim = c(-1.338257, 6.61743), expand = FALSE) + 
          theme(panel.spacing.x = unit(1.5, "lines"),
                axis.text.x = element_text(margin = margin(b = 2)),
                axis.title.y = element_text(margin = margin(t = 0, r = 5, b = 0, l = 0)))
        
        # Plot for Barren + Grass + Shrub
        ggplot(data = df, aes(x = x)) +
          geom_ribbon(aes(ymin = lBarren, ymax = upBarren, fill = "Barren"), alpha = 0.3) +
          geom_ribbon(aes(ymin = lGrass, ymax = upGrass, fill = "Grass"), alpha = 0.3) +
          geom_ribbon(aes(ymin = lShrub, ymax = upShrub, fill = "Shrub"), alpha = 0.3) +
          geom_line(aes(y = Barren, color = "Barren"), size = 1.5) +
          geom_line(aes(y = Grass, color = "Grass"), size = 1.5) +
          geom_line(aes(y = Shrub, color = "Shrub"), size = 1.5) +
          facet_grid(Greenness ~ Season, labeller = labeller(Sex = c(F = "Female"))) +
          xlab("Distance To Facilities (m)") +
          ylab("Probability Of Occurrence (%)") +
          scale_color_manual(values = c("Barren" = "brown", "Grass" = "magenta", "Shrub" = "coral"), name = "Habitat") +
          scale_fill_manual(values = c("Barren" = "brown", "Grass" = "magenta", "Shrub" = "coral"), name = "Habitat") +
          scale_x_continuous(breaks = c(-1.338257, 0.661743, 2.661743, 4.661743, 6.61743),
                             labels = c(0, 2687, 5373, 8060, 10747)) +
          scale_y_continuous(breaks = c(0.15, 0.20, 0.25),
                             labels = c(15, 20, 25)) +
          ggtitle("Female - Facilities with Greenness Variations") +
          theme_bw() +
          coord_cartesian(xlim = c(-1.338257, 6.61743), expand = FALSE) + 
          theme(panel.spacing.x = unit(1.5, "lines"),
                axis.text.x = element_text(margin = margin(b = 2)),
                axis.title.y = element_text(margin = margin(t = 0, r = 5, b = 0, l = 0)))
  
  
  
  
  
  
  #Road: Female x Season
            # check ranges
            range(dat$distance_road)
            range(dat$road)
            
            mean(dat$distance_road)
            sd(dat$distance_road)
            
            mean(dat$distance_road) + sd(dat$distance_road) * 6.5549247
            mean(dat$distance_road) + sd(dat$distance_road) * 0
            mean(dat$distance_road) + sd(dat$distance_road) * -0.867
            
            #determine xlim
            sm = seq(from=-2, to=7.5, by=0.1) #check the ranges of the covariates above; e.g. facility: -1.338 to 6.269
            
            
            # Get mean and standard deviation of greenness from your dataset
            mean_greenness <- mean(dat$greenness, na.rm = TRUE)
            sd_greenness <- sd(dat$greenness, na.rm = TRUE)
            
            # Define greenness levels as mean, -2 SD, and +2 SD
            greenness_levels <- c(mean_greenness - 2 * sd_greenness, mean_greenness, mean_greenness + 2 * sd_greenness)
            
            # Function to create data frames with adjusted greenness levels for road
            create_new_data <- function(season, greenness) {
              expand.grid(
                facility = 0,
                road = sm, # Road is now varied
                trail = 0,
                stream = 0,
                river = 0,
                greenness = greenness, 
                Sex = "F",
                size = 1, 
                Season = season
              )
            }
            
            # Create data frames for each season and greenness level
            new_FRF <- bind_rows(
              create_new_data("Fall", greenness_levels[1]),
              create_new_data("Fall", greenness_levels[2]),
              create_new_data("Fall", greenness_levels[3])
            )
            
            new_FRSP <- bind_rows(
              create_new_data("Spring", greenness_levels[1]),
              create_new_data("Spring", greenness_levels[2]),
              create_new_data("Spring", greenness_levels[3])
            )
            
            new_FRSU <- bind_rows(
              create_new_data("Summer", greenness_levels[1]),
              create_new_data("Summer", greenness_levels[2]),
              create_new_data("Summer", greenness_levels[3])
            )
            
            # Generate posterior predictions for each dataset
            pdf_FRF <- posterior_epred(fit, newdata = new_FRF, re_formula = NA, ndraws = 1000)
            pdf_FRSP <- posterior_epred(fit, newdata = new_FRSP, re_formula = NA, ndraws = 1000)
            pdf_FRSU <- posterior_epred(fit, newdata = new_FRSU, re_formula = NA, ndraws = 1000)
            
            # Function to calculate means and confidence intervals
            calc_means_cis <- function(pdf) {
              list(
                ciAspen = t(hdi(pdf[, , 1])),
                ciBarren = t(hdi(pdf[, , 2])),
                ciConifer = t(hdi(pdf[, , 3])),
                ciGrass = t(hdi(pdf[, , 4])),
                ciShrub = t(hdi(pdf[, , 5])),
                meanAspen = apply(pdf[, , 1], 2, mean),
                meanBarren = apply(pdf[, , 2], 2, mean),
                meanConifer = apply(pdf[, , 3], 2, mean),
                meanGrass = apply(pdf[, , 4], 2, mean),
                meanShrub = apply(pdf[, , 5], 2, mean)
              )
            }
            
            means_cis_FRF <- calc_means_cis(pdf_FRF)
            means_cis_FRSP <- calc_means_cis(pdf_FRSP)
            means_cis_FRSU <- calc_means_cis(pdf_FRSU)
            
            # Function to create data frame for plotting
            create_df <- function(sm, means_cis, season, greenness) {
              data.frame(
                x = sm,
                Aspen = means_cis$meanAspen,
                lAspen = means_cis$ciAspen[, 1],
                upAspen = means_cis$ciAspen[, 2],
                Barren = means_cis$meanBarren,
                lBarren = means_cis$ciBarren[, 1],
                upBarren = means_cis$ciBarren[, 2],
                Conifer = means_cis$meanConifer,
                lConifer = means_cis$ciConifer[, 1],
                upConifer = means_cis$ciConifer[, 2],
                Grass = means_cis$meanGrass,
                lGrass = means_cis$ciGrass[, 1],
                upGrass = means_cis$ciGrass[, 2],
                Shrub = means_cis$meanShrub,
                lShrub = means_cis$ciShrub[, 1],
                upShrub = means_cis$ciShrub[, 2],
                Sex = "F",
                Season = season,
                Greenness = greenness
              )
            }
            
            # Create data frames for plotting with the specified greenness levels
            df_FRF <- create_df(sm, means_cis_FRF, "Fall", rep(greenness_levels, each = length(sm)))
            df_FRSP <- create_df(sm, means_cis_FRSP, "Spring", rep(greenness_levels, each = length(sm)))
            df_FRSU <- create_df(sm, means_cis_FRSU, "Summer", rep(greenness_levels, each = length(sm)))
            
            # Combine all data frames for plotting
            df <- bind_rows(df_FRF, df_FRSP, df_FRSU)
            df$Season <- factor(df$Season, levels = c("Spring", "Summer", "Fall"))
            df$Greenness <- factor(df$Greenness, levels = greenness_levels, labels = c("Low", "Mean", "High"))
            
            # Define colors for the plots
            colors <- c("Aspen" = "limegreen",
                        "Barren" = "brown",
                        "Conifer" = "yellow2",
                        "Grass" = "magenta",
                        "Shrub" = "coral")
            
            # Plot for Aspen + Conifer
            ggplot(data = df, aes(x = x)) +
              geom_ribbon(aes(ymin = lAspen, ymax = upAspen, fill = "Aspen"), alpha = 0.3) +
              geom_ribbon(aes(ymin = lConifer, ymax = upConifer, fill = "Conifer"), alpha = 0.3) +
              geom_line(aes(y = Aspen, color = "Aspen"), size = 1.5) +
              geom_line(aes(y = Conifer, color = "Conifer"), size = 1.5) +
              facet_grid(Greenness ~ Season, labeller = labeller(Sex = c(F = "Female"))) +
              xlab("Distance To Roads (m)") +
              ylab("Probability Of Occurrence (%)") +
              scale_color_manual(values = c("Aspen" = "limegreen", "Conifer" = "yellow2"), name = "Habitat") +
              scale_fill_manual(values = c("Aspen" = "limegreen", "Conifer" = "yellow2"), name = "Habitat") +
              scale_x_continuous(breaks = c(-0.8670253, 1.132975, 3.132975, 5.132975, 7.132975),
                                 labels = c(0, 2496, 4991, 7487, 9982)) +
              scale_y_continuous(breaks = c(0.15, 0.20, 0.25, 0.30, 0.35),
                                 labels = c(15, 20, 25, 30, 35)) +
              ggtitle("Female - Roads with Greenness Variations") +
              theme_bw() +
              coord_cartesian(xlim = c(-0.8670253, 7.132975), expand = FALSE) + 
              theme(panel.spacing.x = unit(1.5, "lines"),
                    axis.text.x = element_text(margin = margin(b = 2)),
                    axis.title.y = element_text(margin = margin(t = 0, r = 5, b = 0, l = 0)))
            
            # Plot for Barren + Grass + Shrub
            ggplot(data = df, aes(x = x)) +
              geom_ribbon(aes(ymin = lBarren, ymax = upBarren, fill = "Barren"), alpha = 0.3) +
              geom_ribbon(aes(ymin = lGrass, ymax = upGrass, fill = "Grass"), alpha = 0.3) +
              geom_ribbon(aes(ymin = lShrub, ymax = upShrub, fill = "Shrub"), alpha = 0.3) +
              geom_line(aes(y = Barren, color = "Barren"), size = 1.5) +
              geom_line(aes(y = Grass, color = "Grass"), size = 1.5) +
              geom_line(aes(y = Shrub, color = "Shrub"), size = 1.5) +
              facet_grid(Greenness ~ Season, labeller = labeller(Sex = c(F = "Female"))) +
              xlab("Distance To Roads (m)") +
              ylab("Probability Of Occurrence (%)") +
              scale_color_manual(values = c("Barren" = "brown", "Grass" = "magenta", "Shrub" = "coral"), name = "Habitat") +
              scale_fill_manual(values = c("Barren" = "brown", "Grass" = "magenta", "Shrub" = "coral"), name = "Habitat") +
              scale_x_continuous(breaks = c(-0.8670253, 1.132975, 3.132975, 5.132975, 7.132975),
                                 labels = c(0, 2496, 4991, 7487, 9982)) +
              scale_y_continuous(breaks = c(0.15, 0.20, 0.25, 0.30, 0.35),
                                 labels = c(15, 20, 25, 30, 35)) +
              ggtitle("Female - Roads with Greenness Variations") +
              theme_bw() +
              coord_cartesian(xlim = c(-0.8670253, 7.132975), expand = FALSE) + 
              theme(panel.spacing.x = unit(1.5, "lines"),
                    axis.text.x = element_text(margin = margin(b = 2)),
                    axis.title.y = element_text(margin = margin(t = 0, r = 5, b = 0, l = 0)))
            
  
  
  
  
  #Trail: Female x Season
        # check ranges
        range(dat$distance_trail)
        range(dat$trail)
        
        mean(dat$distance_trail)
        sd(dat$distance_trail)
        
        mean(dat$distance_trail) + sd(dat$distance_trail) * 8.026
        mean(dat$distance_trail) + sd(dat$distance_trail) * 0
        mean(dat$distance_trail) + sd(dat$distance_trail) * -0.9726994
        
        #determine xlim
        sm = seq(from=-1, to=8.5, by=0.1) #check the ranges of the covariates above; e.g. facility: -1.338 to 6.269
        
        
        # Get mean and standard deviation of greenness from your dataset
        mean_greenness <- mean(dat$greenness, na.rm = TRUE)
        sd_greenness <- sd(dat$greenness, na.rm = TRUE)
        
        # Define greenness levels as mean, -2 SD, and +2 SD
        greenness_levels <- c(mean_greenness - 2 * sd_greenness, mean_greenness, mean_greenness + 2 * sd_greenness)
        
        # Function to create data frames with adjusted greenness levels for trail
        create_new_data <- function(season, greenness) {
          expand.grid(
            facility = 0,
            road = 0,
            trail = sm, # Trail is now varied
            stream = 0,
            river = 0,
            greenness = greenness, 
            Sex = "F",
            size = 1, 
            Season = season
          )
        }
        
        # Create data frames for each season and greenness level
        new_FTF <- bind_rows(
          create_new_data("Fall", greenness_levels[1]),
          create_new_data("Fall", greenness_levels[2]),
          create_new_data("Fall", greenness_levels[3])
        )
        
        new_FTSP <- bind_rows(
          create_new_data("Spring", greenness_levels[1]),
          create_new_data("Spring", greenness_levels[2]),
          create_new_data("Spring", greenness_levels[3])
        )
        
        new_FTSU <- bind_rows(
          create_new_data("Summer", greenness_levels[1]),
          create_new_data("Summer", greenness_levels[2]),
          create_new_data("Summer", greenness_levels[3])
        )
        
        # Generate posterior predictions for each dataset
        pdf_FTF <- posterior_epred(fit, newdata = new_FTF, re_formula = NA, ndraws = 1000)
        pdf_FTSP <- posterior_epred(fit, newdata = new_FTSP, re_formula = NA, ndraws = 1000)
        pdf_FTSU <- posterior_epred(fit, newdata = new_FTSU, re_formula = NA, ndraws = 1000)
        
        # Function to calculate means and confidence intervals
        calc_means_cis <- function(pdf) {
          list(
            ciAspen = t(hdi(pdf[, , 1])),
            ciBarren = t(hdi(pdf[, , 2])),
            ciConifer = t(hdi(pdf[, , 3])),
            ciGrass = t(hdi(pdf[, , 4])),
            ciShrub = t(hdi(pdf[, , 5])),
            meanAspen = apply(pdf[, , 1], 2, mean),
            meanBarren = apply(pdf[, , 2], 2, mean),
            meanConifer = apply(pdf[, , 3], 2, mean),
            meanGrass = apply(pdf[, , 4], 2, mean),
            meanShrub = apply(pdf[, , 5], 2, mean)
          )
        }
        
        means_cis_FTF <- calc_means_cis(pdf_FTF)
        means_cis_FTSP <- calc_means_cis(pdf_FTSP)
        means_cis_FTSU <- calc_means_cis(pdf_FTSU)
        
        # Function to create data frame for plotting
        create_df <- function(sm, means_cis, season, greenness) {
          data.frame(
            x = sm,
            Aspen = means_cis$meanAspen,
            lAspen = means_cis$ciAspen[, 1],
            upAspen = means_cis$ciAspen[, 2],
            Barren = means_cis$meanBarren,
            lBarren = means_cis$ciBarren[, 1],
            upBarren = means_cis$ciBarren[, 2],
            Conifer = means_cis$meanConifer,
            lConifer = means_cis$ciConifer[, 1],
            upConifer = means_cis$ciConifer[, 2],
            Grass = means_cis$meanGrass,
            lGrass = means_cis$ciGrass[, 1],
            upGrass = means_cis$ciGrass[, 2],
            Shrub = means_cis$meanShrub,
            lShrub = means_cis$ciShrub[, 1],
            upShrub = means_cis$ciShrub[, 2],
            Sex = "F",
            Season = season,
            Greenness = greenness
          )
        }
        
        # Create data frames for plotting with the specified greenness levels
        df_FTF <- create_df(sm, means_cis_FTF, "Fall", rep(greenness_levels, each = length(sm)))
        df_FTSP <- create_df(sm, means_cis_FTSP, "Spring", rep(greenness_levels, each = length(sm)))
        df_FTSU <- create_df(sm, means_cis_FTSU, "Summer", rep(greenness_levels, each = length(sm)))
        
        # Combine all data frames for plotting
        df <- bind_rows(df_FTF, df_FTSP, df_FTSU)
        df$Season <- factor(df$Season, levels = c("Spring", "Summer", "Fall"))
        df$Greenness <- factor(df$Greenness, levels = greenness_levels, labels = c("Low", "Mean", "High"))
        
        # Define colors for the plots
        colors <- c("Aspen" = "limegreen",
                    "Barren" = "brown",
                    "Conifer" = "yellow2",
                    "Grass" = "magenta",
                    "Shrub" = "coral")
        
        # Plot for Aspen + Conifer
        ggplot(data = df, aes(x = x)) +
          geom_ribbon(aes(ymin = lAspen, ymax = upAspen, fill = "Aspen"), alpha = 0.3) +
          geom_ribbon(aes(ymin = lConifer, ymax = upConifer, fill = "Conifer"), alpha = 0.3) +
          geom_line(aes(y = Aspen, color = "Aspen"), size = 1.5) +
          geom_line(aes(y = Conifer, color = "Conifer"), size = 1.5) +
          facet_grid(Greenness ~ Season, labeller = labeller(Sex = c(F = "Female"))) +
          xlab("Distance To Trails (m)") +
          ylab("Probability Of Occurrence (%)") +
          scale_color_manual(values = c("Aspen" = "limegreen", "Conifer" = "yellow2"), name = "Habitat") +
          scale_fill_manual(values = c("Aspen" = "limegreen", "Conifer" = "yellow2"), name = "Habitat") +
          scale_x_continuous(breaks = c(-0.9726994, 1.027301, 3.027301, 5.027301, 7.027301),
                             labels = c(0, 1422, 2843, 4265, 5686)) +
          scale_y_continuous(breaks = c(0.15, 0.20, 0.25, 0.30),
                             labels = c(15, 20, 25, 30)) +
          ggtitle("Female - Trails with Greenness Variations") +
          theme_bw() +
          coord_cartesian(xlim = c(-0.9726994, 7.027301), expand = FALSE) + 
          theme(panel.spacing.x = unit(1.5, "lines"),
                axis.text.x = element_text(margin = margin(b = 2)),
                axis.title.y = element_text(margin = margin(t = 0, r = 5, b = 0, l = 0)))
        
        # Plot for Barren + Grass + Shrub
        ggplot(data = df, aes(x = x)) +
          geom_ribbon(aes(ymin = lBarren, ymax = upBarren, fill = "Barren"), alpha = 0.3) +
          geom_ribbon(aes(ymin = lGrass, ymax = upGrass, fill = "Grass"), alpha = 0.3) +
          geom_ribbon(aes(ymin = lShrub, ymax = upShrub, fill = "Shrub"), alpha = 0.3) +
          geom_line(aes(y = Barren, color = "Barren"), size = 1.5) +
          geom_line(aes(y = Grass, color = "Grass"), size = 1.5) +
          geom_line(aes(y = Shrub, color = "Shrub"), size = 1.5) +
          facet_grid(Greenness ~ Season, labeller = labeller(Sex = c(F = "Female"))) +
          xlab("Distance To Trails (m)") +
          ylab("Probability Of Occurrence (%)") +
          scale_color_manual(values = c("Barren" = "brown", "Grass" = "magenta", "Shrub" = "coral"), name = "Habitat") +
          scale_fill_manual(values = c("Barren" = "brown", "Grass" = "magenta", "Shrub" = "coral"), name = "Habitat") +
          scale_x_continuous(breaks = c(-0.9726994, 1.027301, 3.027301, 5.027301, 7.027301),
                             labels = c(0, 1422, 2843, 4265, 5686)) +
          scale_y_continuous(breaks = c(0.15, 0.20, 0.25, 0.30),
                             labels = c(15, 20, 25, 30)) +
          ggtitle("Female - Trails with Greenness Variations") +
          theme_bw() +
          coord_cartesian(xlim = c(-0.9726994, 7.027301), expand = FALSE) + 
          theme(panel.spacing.x = unit(1.5, "lines"),
                axis.text.x = element_text(margin = margin(b = 2)),
                axis.title.y = element_text(margin = margin(t = 0, r = 5, b = 0, l = 0)))
        
        
        
        
        
        
        
        
        
        
        
        
        
        
# male x female -facilities
        
        # Define the sequence for 'sm' based on the provided range
        sm <- seq(from = -2, to = 7, by = 0.1)
        
        # Get mean and standard deviation of greenness from your dataset
        mean_greenness <- mean(dat$greenness, na.rm = TRUE)
        sd_greenness <- sd(dat$greenness, na.rm = TRUE)
        
        # Define greenness levels as mean, -2 SD, and +2 SD
        greenness_levels <- c(mean_greenness - 2 * sd_greenness, mean_greenness, mean_greenness + 2 * sd_greenness)
        
        # Function to create data frames with adjusted greenness levels
        create_new_data <- function(season, sex) {
          expand.grid(
            facility = sm,
            road = 0,
            trail = 0,
            stream = 0,
            river = 0,
            greenness = greenness_levels,
            Sex = sex,
            size = 1,
            Season = season
          )
        }
        
        # Create data frames for each season and sex
        new_MF <- create_new_data("Fall", "M")
        new_MSP <- create_new_data("Spring", "M")
        new_MSU <- create_new_data("Summer", "M")
        new_FF <- create_new_data("Fall", "F")
        new_FSP <- create_new_data("Spring", "F")
        new_FSU <- create_new_data("Summer", "F")
        
        # Generate posterior predictions for each dataset
        pdf_MF <- posterior_epred(fit, newdata = new_MF, re_formula = NA, ndraws = 1000)
        pdf_MSP <- posterior_epred(fit, newdata = new_MSP, re_formula = NA, ndraws = 1000)
        pdf_MSU <- posterior_epred(fit, newdata = new_MSU, re_formula = NA, ndraws = 1000)
        pdf_FF <- posterior_epred(fit, newdata = new_FF, re_formula = NA, ndraws = 1000)
        pdf_FSP <- posterior_epred(fit, newdata = new_FSP, re_formula = NA, ndraws = 1000)
        pdf_FSU <- posterior_epred(fit, newdata = new_FSU, re_formula = NA, ndraws = 1000)
        
        # Function to calculate means and confidence intervals
        calc_means_cis <- function(pdf) {
          list(
            ciAspen = t(hdi(pdf[, , 1])),
            ciBarren = t(hdi(pdf[, , 2])),
            ciConifer = t(hdi(pdf[, , 3])),
            ciGrass = t(hdi(pdf[, , 4])),
            ciShrub = t(hdi(pdf[, , 5])),
            meanAspen = apply(pdf[, , 1], 2, mean),
            meanBarren = apply(pdf[, , 2], 2, mean),
            meanConifer = apply(pdf[, , 3], 2, mean),
            meanGrass = apply(pdf[, , 4], 2, mean),
            meanShrub = apply(pdf[, , 5], 2, mean)
          )
        }
        
        # Calculate means and confidence intervals for each dataset
        means_cis_MF <- calc_means_cis(pdf_MF)
        means_cis_MSP <- calc_means_cis(pdf_MSP)
        means_cis_MSU <- calc_means_cis(pdf_MSU)
        means_cis_FF <- calc_means_cis(pdf_FF)
        means_cis_FSP <- calc_means_cis(pdf_FSP)
        means_cis_FSU <- calc_means_cis(pdf_FSU)
        
        # Function to create data frames for plotting
        create_df <- function(sm, means_cis, season, sex) {
          data.frame(
            x = sm,
            Aspen = means_cis$meanAspen,
            lAspen = means_cis$ciAspen[, 1],
            upAspen = means_cis$ciAspen[, 2],
            Barren = means_cis$meanBarren,
            lBarren = means_cis$ciBarren[, 1],
            upBarren = means_cis$ciBarren[, 2],
            Conifer = means_cis$meanConifer,
            lConifer = means_cis$ciConifer[, 1],
            upConifer = means_cis$ciConifer[, 2],
            Grass = means_cis$meanGrass,
            lGrass = means_cis$ciGrass[, 1],
            upGrass = means_cis$ciGrass[, 2],
            Shrub = means_cis$meanShrub,
            lShrub = means_cis$ciShrub[, 1],
            upShrub = means_cis$ciShrub[, 2],
            Sex = sex,
            Season = season,
            greenness = rep(greenness_levels, each = length(sm))
          )
        }
        
        # Create data frames for plotting
        df_MF <- create_df(sm, means_cis_MF, "Fall", "M")
        df_MSP <- create_df(sm, means_cis_MSP, "Spring", "M")
        df_MSU <- create_df(sm, means_cis_MSU, "Summer", "M")
        df_FF <- create_df(sm, means_cis_FF, "Fall", "F")
        df_FSP <- create_df(sm, means_cis_FSP, "Spring", "F")
        df_FSU <- create_df(sm, means_cis_FSU, "Summer", "F")
        
        # Combine all data frames for plotting
        df <- bind_rows(df_MF, df_MSP, df_MSU, df_FF, df_FSP, df_FSU)
        
        # Set greenness factor levels for correct ordering
        df$greenness <- factor(df$greenness, levels = greenness_levels, labels = c("Low Greenness", "Mean Greenness", "High Greenness"))
        
        # Define colors for the plots
        colors <- c("Aspen" = "limegreen",
                    "Barren" = "brown",
                    "Conifer" = "yellow2",
                    "Grass" = "magenta",
                    "Shrub" = "coral")
        
        # Function to create individual season plots with Female and Male side-by-side
        create_season_plot <- function(season) {
          ggplot(data = df %>% filter(Season == season), aes(x = x)) +
            geom_ribbon(aes(ymin = lAspen, ymax = upAspen, fill = "Aspen"), alpha = 0.3) +
            geom_ribbon(aes(ymin = lBarren, ymax = upBarren, fill = "Barren"), alpha = 0.3) +
            geom_ribbon(aes(ymin = lConifer, ymax = upConifer, fill = "Conifer"), alpha = 0.3) +
            geom_ribbon(aes(ymin = lShrub, ymax = upShrub, fill = "Shrub"), alpha = 0.3) +
            geom_ribbon(aes(ymin = lGrass, ymax = upGrass, fill = "Grass"), alpha = 0.3) +
            geom_line(aes(y = Aspen, color = "Aspen"), size = 1) +
            geom_line(aes(y = Barren, color = "Barren"), size = 1) + 
            geom_line(aes(y = Conifer, color = "Conifer"), size = 1) +  
            geom_line(aes(y = Shrub, color = "Shrub"), size = 1) +  
            geom_line(aes(y = Grass, color = "Grass"), size = 1) +  
            facet_grid(greenness ~ Sex) +
            xlab("Distance To Facilities (m)") +
            ylab("Probability Of Occurrence (%)") +
            scale_color_manual(values = colors, labels = names(colors), name = "Habitat") +
            scale_fill_manual(values = colors, labels = names(colors), name = "Habitat") +
            scale_x_continuous(breaks = c(-1.338257, 0.661743, 2.661743, 4.661743, 6.661743),
                               labels = c(0, 2687, 5373, 8060, 10747)) +
            scale_y_continuous(breaks = c(0.15, 0.20, 0.25),
                               labels = c(15, 20, 25)) +
            theme_bw() +
            coord_cartesian(xlim = c(-1.338257, 6.661743), expand = FALSE) +
            theme(panel.spacing.x = unit(1.5, "lines"),
                  axis.text.x = element_text(margin = margin(b = 2)),
                  axis.title.y = element_text(margin = margin(t = 0, r = 5, b = 0, l = 0))) +
            ggtitle(paste(season, "- Female and Male"))
        }
        
        # Generate plots for each season separately
        plot_spring <- create_season_plot("Spring")
        plot_summer <- create_season_plot("Summer")
        plot_fall <- create_season_plot("Fall")
        
        # Display plots individually
        print(plot_spring)
        print(plot_summer)
        print(plot_fall)
        
        
#Male x Female - Road
        # Load necessary libraries
        library(ggplot2)
        library(dplyr)
        
        # Define the sequence for 'sm' based on the road data
        sm <- seq(from = -2, to = 7.5, by = 0.1) # Adjust based on your specific data ranges
        
        # Get mean and standard deviation of greenness from your dataset
        mean_greenness <- mean(dat$greenness, na.rm = TRUE)
        sd_greenness <- sd(dat$greenness, na.rm = TRUE)
        
        # Define greenness levels as mean, -2 SD, and +2 SD
        greenness_levels <- c(mean_greenness - 2 * sd_greenness, mean_greenness, mean_greenness + 2 * sd_greenness)
        
        # Function to create data frames with adjusted greenness levels for road
        create_new_data <- function(season, sex) {
          expand.grid(
            facility = 0,
            road = sm, # Road is now varied
            trail = 0,
            stream = 0,
            river = 0,
            greenness = greenness_levels, 
            Sex = sex,
            size = 1, 
            Season = season
          )
        }
        
        # Create data frames for each season and greenness level for both males and females
        new_Fall_M <- create_new_data("Fall", "M")
        new_Fall_F <- create_new_data("Fall", "F")
        new_Spring_M <- create_new_data("Spring", "M")
        new_Spring_F <- create_new_data("Spring", "F")
        new_Summer_M <- create_new_data("Summer", "M")
        new_Summer_F <- create_new_data("Summer", "F")
        
        # Generate posterior predictions for each dataset
        pdf_Fall_M <- posterior_epred(fit, newdata = new_Fall_M, re_formula = NA, ndraws = 1000)
        pdf_Fall_F <- posterior_epred(fit, newdata = new_Fall_F, re_formula = NA, ndraws = 1000)
        pdf_Spring_M <- posterior_epred(fit, newdata = new_Spring_M, re_formula = NA, ndraws = 1000)
        pdf_Spring_F <- posterior_epred(fit, newdata = new_Spring_F, re_formula = NA, ndraws = 1000)
        pdf_Summer_M <- posterior_epred(fit, newdata = new_Summer_M, re_formula = NA, ndraws = 1000)
        pdf_Summer_F <- posterior_epred(fit, newdata = new_Summer_F, re_formula = NA, ndraws = 1000)
        
        # Function to calculate means and confidence intervals
        calc_means_cis <- function(pdf) {
          list(
            ciAspen = t(hdi(pdf[, , 1])),
            ciBarren = t(hdi(pdf[, , 2])),
            ciConifer = t(hdi(pdf[, , 3])),
            ciGrass = t(hdi(pdf[, , 4])),
            ciShrub = t(hdi(pdf[, , 5])),
            meanAspen = apply(pdf[, , 1], 2, mean),
            meanBarren = apply(pdf[, , 2], 2, mean),
            meanConifer = apply(pdf[, , 3], 2, mean),
            meanGrass = apply(pdf[, , 4], 2, mean),
            meanShrub = apply(pdf[, , 5], 2, mean)
          )
        }
        
        # Calculate means and confidence intervals for each dataset
        means_cis_Fall_M <- calc_means_cis(pdf_Fall_M)
        means_cis_Fall_F <- calc_means_cis(pdf_Fall_F)
        means_cis_Spring_M <- calc_means_cis(pdf_Spring_M)
        means_cis_Spring_F <- calc_means_cis(pdf_Spring_F)
        means_cis_Summer_M <- calc_means_cis(pdf_Summer_M)
        means_cis_Summer_F <- calc_means_cis(pdf_Summer_F)
        
        # Function to create data frame for plotting
        create_df <- function(sm, means_cis, season, sex) {
          data.frame(
            x = sm,
            Aspen = means_cis$meanAspen,
            lAspen = means_cis$ciAspen[, 1],
            upAspen = means_cis$ciAspen[, 2],
            Barren = means_cis$meanBarren,
            lBarren = means_cis$ciBarren[, 1],
            upBarren = means_cis$ciBarren[, 2],
            Conifer = means_cis$meanConifer,
            lConifer = means_cis$ciConifer[, 1],
            upConifer = means_cis$ciConifer[, 2],
            Grass = means_cis$meanGrass,
            lGrass = means_cis$ciGrass[, 1],
            upGrass = means_cis$ciGrass[, 2],
            Shrub = means_cis$meanShrub,
            lShrub = means_cis$ciShrub[, 1],
            upShrub = means_cis$ciShrub[, 2],
            Sex = sex,
            Season = season,
            greenness = rep(greenness_levels, each = length(sm))
          )
        }
        
        # Create data frames for plotting
        df_Fall_M <- create_df(sm, means_cis_Fall_M, "Fall", "M")
        df_Fall_F <- create_df(sm, means_cis_Fall_F, "Fall", "F")
        df_Spring_M <- create_df(sm, means_cis_Spring_M, "Spring", "M")
        df_Spring_F <- create_df(sm, means_cis_Spring_F, "Spring", "F")
        df_Summer_M <- create_df(sm, means_cis_Summer_M, "Summer", "M")
        df_Summer_F <- create_df(sm, means_cis_Summer_F, "Summer", "F")
        
        # Combine all data frames for plotting
        df <- bind_rows(df_Fall_M, df_Fall_F, df_Spring_M, df_Spring_F, df_Summer_M, df_Summer_F)
        df$Season <- factor(df$Season, levels = c("Spring", "Summer", "Fall"))
        df$greenness <- factor(df$greenness, levels = greenness_levels, labels = c("Low Greenness", "Mean Greenness", "High Greenness"))
        
        # Define colors for the plots
        colors <- c("Aspen" = "limegreen", "Barren" = "brown", "Conifer" = "yellow2", "Grass" = "magenta", "Shrub" = "coral")
        
        # Function to create individual season plots with Female and Male side-by-side
        create_season_plot <- function(season) {
          ggplot(data = df %>% filter(Season == season), aes(x = x)) +
            geom_ribbon(aes(ymin = lAspen, ymax = upAspen, fill = "Aspen"), alpha = 0.3) +
            geom_ribbon(aes(ymin = lBarren, ymax = upBarren, fill = "Barren"), alpha = 0.3) +
            geom_ribbon(aes(ymin = lConifer, ymax = upConifer, fill = "Conifer"), alpha = 0.3) +
            geom_ribbon(aes(ymin = lShrub, ymax = upShrub, fill = "Shrub"), alpha = 0.3) +
            geom_ribbon(aes(ymin = lGrass, ymax = upGrass, fill = "Grass"), alpha = 0.3) +
            geom_line(aes(y = Aspen, color = "Aspen"), size = 1) +
            geom_line(aes(y = Barren, color = "Barren"), size = 1) + 
            geom_line(aes(y = Conifer, color = "Conifer"), size = 1) +  
            geom_line(aes(y = Shrub, color = "Shrub"), size = 1) +  
            geom_line(aes(y = Grass, color = "Grass"), size = 1) +  
            facet_grid(greenness ~ Sex) +
            xlab("Distance To Roads (m)") +
            ylab("Probability Of Occurrence (%)") +
            scale_color_manual(values = colors, labels = names(colors), name = "Habitat") +
            scale_fill_manual(values = colors, labels = names(colors), name = "Habitat") +
            scale_x_continuous(breaks = c(-0.867, 1.133, 3.133, 5.133, 7.133),
                               labels = c(0, 2497, 4992, 7487, 9983)) +
            scale_y_continuous(breaks = c(0.15, 0.20, 0.25, 0.30, 0.35),
                               labels = c(15, 20, 25, 30, 35)) +
            theme_bw() +
            coord_cartesian(xlim = c(-0.867, 7.133), expand = FALSE) +
            theme(panel.spacing.x = unit(1.5, "lines"), 
                  axis.text.x = element_text(margin = margin(b = 2)),
                  axis.title.y = element_text(margin = margin(t = 0, r = 5, b = 0, l = 0))) +
            ggtitle(paste(season, "- Female and Male"))
        }
        
        # Generate plots for each season separately
        plot_spring <- create_season_plot("Spring")
        plot_summer <- create_season_plot("Summer")
        plot_fall <- create_season_plot("Fall")
        
        # Display plots individually
        print(plot_spring)
        print(plot_summer)
        print(plot_fall)
        
        
#Male x Female - Trail        
        # Load necessary libraries
        library(ggplot2)
        library(dplyr)
        
        # Define the sequence for 'sm' based on the trail data
        sm <- seq(from = -2, to = 7.5, by = 0.1) # Adjust based on your specific data ranges
        
        # Get mean and standard deviation of greenness from your dataset
        mean_greenness <- mean(dat$greenness, na.rm = TRUE)
        sd_greenness <- sd(dat$greenness, na.rm = TRUE)
        
        # Define greenness levels as mean, -2 SD, and +2 SD
        greenness_levels <- c(mean_greenness - 2 * sd_greenness, mean_greenness, mean_greenness + 2 * sd_greenness)
        
        # Function to create data frames with adjusted greenness levels for trail
        create_new_data <- function(season, sex) {
          expand.grid(
            facility = 0,
            road = 0,
            trail = sm, # Trail is now varied
            stream = 0,
            river = 0,
            greenness = greenness_levels, 
            Sex = sex,
            size = 1, 
            Season = season
          )
        }
        
        # Create data frames for each season and greenness level for both males and females
        new_Fall_M <- create_new_data("Fall", "M")
        new_Fall_F <- create_new_data("Fall", "F")
        new_Spring_M <- create_new_data("Spring", "M")
        new_Spring_F <- create_new_data("Spring", "F")
        new_Summer_M <- create_new_data("Summer", "M")
        new_Summer_F <- create_new_data("Summer", "F")
        
        # Generate posterior predictions for each dataset
        pdf_Fall_M <- posterior_epred(fit, newdata = new_Fall_M, re_formula = NA, ndraws = 1000)
        pdf_Fall_F <- posterior_epred(fit, newdata = new_Fall_F, re_formula = NA, ndraws = 1000)
        pdf_Spring_M <- posterior_epred(fit, newdata = new_Spring_M, re_formula = NA, ndraws = 1000)
        pdf_Spring_F <- posterior_epred(fit, newdata = new_Spring_F, re_formula = NA, ndraws = 1000)
        pdf_Summer_M <- posterior_epred(fit, newdata = new_Summer_M, re_formula = NA, ndraws = 1000)
        pdf_Summer_F <- posterior_epred(fit, newdata = new_Summer_F, re_formula = NA, ndraws = 1000)
        
        # Function to calculate means and confidence intervals
        calc_means_cis <- function(pdf) {
          list(
            ciAspen = t(hdi(pdf[, , 1])),
            ciBarren = t(hdi(pdf[, , 2])),
            ciConifer = t(hdi(pdf[, , 3])),
            ciGrass = t(hdi(pdf[, , 4])),
            ciShrub = t(hdi(pdf[, , 5])),
            meanAspen = apply(pdf[, , 1], 2, mean),
            meanBarren = apply(pdf[, , 2], 2, mean),
            meanConifer = apply(pdf[, , 3], 2, mean),
            meanGrass = apply(pdf[, , 4], 2, mean),
            meanShrub = apply(pdf[, , 5], 2, mean)
          )
        }
        
        # Calculate means and confidence intervals for each dataset
        means_cis_Fall_M <- calc_means_cis(pdf_Fall_M)
        means_cis_Fall_F <- calc_means_cis(pdf_Fall_F)
        means_cis_Spring_M <- calc_means_cis(pdf_Spring_M)
        means_cis_Spring_F <- calc_means_cis(pdf_Spring_F)
        means_cis_Summer_M <- calc_means_cis(pdf_Summer_M)
        means_cis_Summer_F <- calc_means_cis(pdf_Summer_F)
        
        # Function to create data frame for plotting
        create_df <- function(sm, means_cis, season, sex) {
          data.frame(
            x = sm,
            Aspen = means_cis$meanAspen,
            lAspen = means_cis$ciAspen[, 1],
            upAspen = means_cis$ciAspen[, 2],
            Barren = means_cis$meanBarren,
            lBarren = means_cis$ciBarren[, 1],
            upBarren = means_cis$ciBarren[, 2],
            Conifer = means_cis$meanConifer,
            lConifer = means_cis$ciConifer[, 1],
            upConifer = means_cis$ciConifer[, 2],
            Grass = means_cis$meanGrass,
            lGrass = means_cis$ciGrass[, 1],
            upGrass = means_cis$ciGrass[, 2],
            Shrub = means_cis$meanShrub,
            lShrub = means_cis$ciShrub[, 1],
            upShrub = means_cis$ciShrub[, 2],
            Sex = sex,
            Season = season,
            greenness = rep(greenness_levels, each = length(sm))
          )
        }
        
        # Create data frames for plotting
        df_Fall_M <- create_df(sm, means_cis_Fall_M, "Fall", "M")
        df_Fall_F <- create_df(sm, means_cis_Fall_F, "Fall", "F")
        df_Spring_M <- create_df(sm, means_cis_Spring_M, "Spring", "M")
        df_Spring_F <- create_df(sm, means_cis_Spring_F, "Spring", "F")
        df_Summer_M <- create_df(sm, means_cis_Summer_M, "Summer", "M")
        df_Summer_F <- create_df(sm, means_cis_Summer_F, "Summer", "F")
        
        # Combine all data frames for plotting
        df <- bind_rows(df_Fall_M, df_Fall_F, df_Spring_M, df_Spring_F, df_Summer_M, df_Summer_F)
        df$Season <- factor(df$Season, levels = c("Spring", "Summer", "Fall"))
        df$greenness <- factor(df$greenness, levels = greenness_levels, labels = c("Low Greenness", "Mean Greenness", "High Greenness"))
        
        # Define colors for the plots
        colors <- c("Aspen" = "limegreen", "Barren" = "brown", "Conifer" = "yellow2", "Grass" = "magenta", "Shrub" = "coral")
        
        # Function to create individual season plots with Male and Female side-by-side
        create_season_plot <- function(season) {
          ggplot(data = df %>% filter(Season == season), aes(x = x)) +
            geom_ribbon(aes(ymin = lAspen, ymax = upAspen, fill = "Aspen"), alpha = 0.3) +
            geom_ribbon(aes(ymin = lBarren, ymax = upBarren, fill = "Barren"), alpha = 0.3) +
            geom_ribbon(aes(ymin = lConifer, ymax = upConifer, fill = "Conifer"), alpha = 0.3) +
            geom_ribbon(aes(ymin = lShrub, ymax = upShrub, fill = "Shrub"), alpha = 0.3) +
            geom_ribbon(aes(ymin = lGrass, ymax = upGrass, fill = "Grass"), alpha = 0.3) +
            geom_line(aes(y = Aspen, color = "Aspen"), size = 1) +
            geom_line(aes(y = Barren, color = "Barren"), size = 1) + 
            geom_line(aes(y = Conifer, color = "Conifer"), size = 1) +  
            geom_line(aes(y = Shrub, color = "Shrub"), size = 1) +  
            geom_line(aes(y = Grass, color = "Grass"), size = 1) +  
            facet_grid(greenness ~ Sex) +
            xlab("Distance To Trails (m)") +
            ylab("Probability Of Occurrence (%)") +
            scale_color_manual(values = colors, labels = names(colors), name = "Habitat") +
            scale_fill_manual(values = colors, labels = names(colors), name = "Habitat") +
            scale_x_continuous(breaks = c(-0.9726994, 1.027301, 3.027301, 5.027301, 7.027301),
                               labels = c(0, 1422, 2843, 4265, 5686)) +
            scale_y_continuous(breaks = c(0.15, 0.20, 0.25, 0.30),
                               labels = c(15, 20, 25, 30)) +
            theme_bw() +
            coord_cartesian(xlim = c(-0.9726994, 7.027301), expand = FALSE) +
            theme(panel.spacing.x = unit(1.5, "lines"), 
                  axis.text.x = element_text(margin = margin(b = 2)),
                  axis.title.y = element_text(margin = margin(t = 0, r = 5, b = 0, l = 0))) +
            ggtitle(paste(season, "- Male and Female"))
        }
        
        # Generate plots for each season separately
        plot_spring <- create_season_plot("Spring")
        plot_summer <- create_season_plot("Summer")
        plot_fall <- create_season_plot("Fall")
        
        # Display plots individually
        print(plot_spring)
        print(plot_summer)
        print(plot_fall)
        
        
        