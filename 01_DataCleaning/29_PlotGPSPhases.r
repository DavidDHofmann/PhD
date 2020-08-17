############################################################
#### Plotting GPS Phases
############################################################
# Description: Preparation of a plot that depicts the time for which we collected
# GPS data for each individual

# Clear R's brain
rm(list = ls())

# Set the working directory
wd <- "/home/david/ownCloud/University/15. PhD/00_WildDogs"
setwd(wd)

# Load required packages
library(data.table)
library(viridis)
library(adehabitatHR)
library(raster)
library(tidyverse)
library(ggplot2)
library(lubridate)
library(rasterVis)
library(rgdal)
library(rgeos)

# Load custom functions
source("Functions.r")

############################################################
#### Plotting GPS Observation Phases
############################################################
# Before we can calculate a social landscape we need to identify the fixes that we can use for this task. Let's load the data containing all GPS fixes
dat <- read_csv("03_Data/02_CleanData/00_General_Dispersers_Popecol(Regular).csv")

# Replace NAs in the current packs after dispersal with character NA (i.e. "NA"). This will avoid some errors.
dat$CurrentPack[is.na(dat$CurrentPack) & dat$State == "Resident"] <- "NA"

# We want to plot a graph that depicts the period for which each individual was collared and the time during which the individual was dispersing. The colours should indicate to which pack each individual belongs. As a first step we get the earliest and latest gps recording for each dog and collar, and current pack
collars <- dat %>%

  # Group by CollarID, DogName and Currentpack
  group_by(., CollarID, DogName, CurrentPack) %>%

  # Now we can calculate the first and last observation for each group.
  summarize(.
    , FirstDate = min(Timestamp)
    , LastDate  = max(Timestamp)
  ) %>%

  # We get rid of the "non-character" NAs since they refer to the either Abrahms' individuals, or the dispersal part of our individuals. We don't want these in our plot anyways so let's remove them
  subset(., !is.na(CurrentPack)) %>%

  # Order the data. We put the focus on the pack identities since we assume that members of the same pack are at the same location. The individual is therefore not really interesting for us.
  arrange(., CurrentPack, FirstDate)

# Make the dog names factorial. Note that the factors are now ordered as we specified them in the ordering above. This is helpful for plotting
collars$DogName <- factor(collars$DogName, levels = unique(collars$DogName))

# We could now make the character "NAs" true NAs again
collars$CurrentPack[collars$CurrentPack == "NA"] <- NA

# Since we also want to depict the dispersal phase, we need to get the dispersal dates. Since these dates are stored in the cutoff dates table we can load it again
cut <- read_csv("03_Data/02_CleanData/00_General_Dispersers_Popecol_CutoffDates.csv")

# Prepare the plot
p <- ggplot(collars, aes(x = rbind(FirstDate, LastDate), y = DogName)) +

  # Add segments for the gps observation phase
  geom_segment(data = collars, aes(
      x       = FirstDate
    , xend    = LastDate
    , y       = DogName
    , yend    = DogName
    , colour  = CurrentPack
  ), size = 3.5) +

  scale_fill_viridis() +

  # Add segments for the dispersal phase
  geom_segment(data = cut, aes(
      x       = StartDate
    , xend    = EndDate
    , y       = DogName
    , yend    = DogName
  ), size = 2, colour = "black") +

  # Revert the y scale (top to bottom)
  scale_y_discrete(limits = rev(levels(collars$DogName))) +

  # Put a useful title and axis labels
  ggtitle("GPS Observations") +
  xlab("Date") +
  ylab("Name") +

  # Use the viridis colour scheme for the entire plot
  scale_color_viridis(
      discrete  = TRUE
    , begin     = 0.3
    , na.value  = "darkgrey"
  )

# Store the plot
# ggsave("GPSObservations.png"
#   , plot = last_plot()
#   , bg = "transparent"
#   , device = "png"
# )

############################################################
#### Same Plot but With Transparency
############################################################
# Prepare the plot
ggplot(collars, aes(x = rbind(FirstDate, LastDate), y = DogName)) +

  # Add segments for the gps observation phase
  geom_segment(data = collars, aes(
      x       = FirstDate
    , xend    = LastDate
    , y       = DogName
    , yend    = DogName
    , colour  = CurrentPack
  ), size = 3.5) +

  scale_fill_viridis() +

  # Add segments for the dispersal phase
  geom_segment(data = cut, aes(
      x       = StartDate
    , xend    = EndDate
    , y       = DogName
    , yend    = DogName
  ), size = 2, colour = "white") +

  # Revert the y scale (top to bottom)
  scale_y_discrete(limits = rev(levels(collars$DogName))) +

  # Put a useful title
  ggtitle("GPS Observations") +

  # Use the viridis colour scheme for the entire plot
  scale_color_viridis(
      discrete  = TRUE
    , begin     = 0.3
    , na.value  = "darkgrey"
  ) +

  theme(

    # Make all rectangles white
    rect = element_rect(fill = "transparent")

    # Remove the axis names
    , axis.title.x = element_blank()
    , axis.title.y = element_blank()

    # Make the axis text white
    , axis.text.x = element_text(colour = "white")
    , axis.text.y = element_text(colour = "white")

    # Make the legend text white and the legend background transparent
    , legend.text = element_text(colour = "white")
    , legend.title = element_text(colour = "white")
    , legend.background = element_rect(colour = "transparent")
    , legend.key = element_rect(fill = "transparent", colour = NA)

    # Make the plot title white and centered
    , plot.title = element_text(hjust = 0.5, colour = "white")

    # Remove the grid
    , panel.grid.major = element_line(colour = "white", size = 0.1)
    , panel.grid.minor = element_line(colour = "white", size = 0.1)
    , panel.background = element_rect(fill = "transparent")

    # Make the axis ticks white
    , axis.ticks = element_line(colour = "white")
  )

# Store the plot
# ggsave("GPSObservations(Transparent).png"
#   , plot = last_plot()
#   , bg = "transparent"
#   , device = "png"
# )

############################################################
#### Same Plot but only with dispersers
############################################################
disp <- "03_Data/02_CleanData/00_General_Dispersers_Popecol(Regular).csv" %>%
  read_csv() %>%
  subset(State == "Disperser") %>%
  select(DogName) %>%
  unique() %>%
  .[[1]]

# Subset to these individuals
collars_sub <- subset(collars, DogName %in% disp) %>% ungroup()
cut_sub <- subset(cut, DogName %in% disp) %>% ungroup()
collars_sub$DogName <- as.character(collars_sub$DogName)
cut_sub$DogName <- as.character(cut_sub$DogName)

# Prepare the plot
p <- ggplot(collars_sub, aes(x = rbind(FirstDate, LastDate), y = DogName)) +

  # Add segments for the gps observation phase
  geom_segment(data = collars_sub, aes(
      x       = FirstDate
    , xend    = LastDate
    , y       = DogName
    , yend    = DogName
    , colour  = CurrentPack
  ), size = 3.5) +

  scale_fill_viridis() +

  # Add segments for the dispersal phase
  geom_segment(data = cut_sub, aes(
      x       = StartDate
    , xend    = EndDate
    , y       = DogName
    , yend    = DogName
  ), size = 2, colour = "black") +

  # Revert the y scale (top to bottom)
  scale_y_discrete(limits = rev(levels(collars_sub$DogName))) +

  # Put a useful title and axis labels
  ggtitle("GPS Observations") +
  xlab("Date") +
  ylab("Coalition") +

  # Use the viridis colour scheme for the entire plot
  scale_color_viridis(
      discrete  = TRUE
    , begin     = 0.3
    , na.value  = "darkgrey"
    , name = "Pack Identity"
  )

# Store the plot
ggsave("99_GPSObservations.png"
  , plot = last_plot()
  , bg = "transparent"
  , device = "png"
  , width = 7
  , height = 5
)
