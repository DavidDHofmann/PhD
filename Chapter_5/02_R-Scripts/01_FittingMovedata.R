################################################################################
#### Fitting Distributions to Movement Data
################################################################################
# Description: Fitting parametric distributions to the turning angles and step
# lengths of wolf data obtained through movebank.
# Authors: David, Eva, & Dilsad
# Date: November 2020

# Clear R's brain
rm(list = ls())

# Load required packages
library(raster)       # To handle spatial data
library(rgeos)        # To manipulate spatial data
library(tidyverse)    # For data wrangling
library(amt)          # To generate steps from gps data
library(fitdistrplus) # To fit step length and turning angle distributions
library(cowplot)      # To add image to ggplot
library(rphylopic)    # To add image to ggplot

# Set the working directory
setwd("/home/david/ownCloud/University/15. PhD/Chapter_5")

# Load custom functions
source("00_Functions.R")

# Set seed for reproducability
set.seed(1234)

################################################################################
#### Fit Movement Parameters to Wolf Data
################################################################################
# Load movement data
dat <- dir(pattern = "csv$", path = "Data/Input", full.names = T)
dat <- lapply(dat, function(x){
  data <- read.csv(x)
  data <- dplyr::select(data, c(
      Timestamp = timestamp
    , x         = location.long
    , y         = location.lat
    , Species   = individual.taxon.canonical.name
    , ID        = individual.local.identifier
    , Source    = study.name
  ))
  return(data)
})
dat <- do.call(rbind, dat)

# Let's make sure that each individual gets a unique ID
dat$ID <- dat %>% group_indices(Species, ID)

# Make proper timestamps
dat$Timestamp <- as.POSIXct(as.character(dat$Timestamp), tz = "UTC")

# Remove NA coordinates
dat <- subset(dat, !is.na(x) & !is.na(y))

# Nest data by the species
dat <- dat %>% group_by(Species) %>% nest()

# We want to create steps that are regularly spaced in time. The data of the two
# species were also collected in different utm zones so we need to make sure to
# project them correctly.
dat$Sampling <- c(4, 4)
dat$CRS <- c(
    CRS("+proj=utm +zone=10 +datum=WGS84") # UTM zone where wolf data comes from
  , CRS("+proj=utm +zone=32 +datum=WGS84") # UTM zone where deer data comes from
)

# Let's also add an image from phylopic to each species
dat$Picture <- list(
    image_data("8cad2b22-30d3-4cbd-86a3-a6d2d004b201", size = "512")[[1]]
  , image_data("b36a215a-adb3-445d-b364-1e63dddd6950", size = "512")[[1]]
)

# Switch from a point to a step representation.
dat$Steps <- lapply(1:nrow(dat), function(x){
  dat$data[[x]] %>%
    make_track(
        .x      = x
      , .y      = y
      , .t      = Timestamp
      , id      = ID
      , crs     = sp::CRS("+init=epsg:4326")
    ) %>%
    transform_coords(dat$CRS[[x]]) %>%
    nest(data = -"id") %>%
    mutate(data = map(data, function(y){
      y %>%
        track_resample(rate = hours(dat$Sampling[x]), tolerance = minutes(15)) %>%
        steps_by_burst(keep_cols = "start")
    })) %>%
    unnest(cols = data) %>%
    mutate(dt_ = as.numeric(dt_, units = "hours")) %>%
    subset(dt_ > dat$Sampling[x] - 0.25 & dt_ < dat$Sampling[x] + 0.25)
})

# Plot step lengths of the two species
p <- dat %>%
  dplyr::select(Species, Steps) %>%
  unnest(cols = Steps) %>%
  ggplot(aes(x = sl_), col = Species) +
    geom_histogram(aes(y = ..density..), fill = "cornflowerblue", color = "black") +
    geom_rug(color = "cornflowerblue") +
    facet_wrap("Species", scales = "free") +
    xlab("Step Length (m)") +
    ylab("Density") +
    theme_cowplot() +
    theme(strip.background = element_rect(fill = "gray95"))

# Add silhouettes
ggdraw() +
  draw_image(dat$Picture[[1]] / 8, x = 0.3, y = 0.4, scale = 0.2, hjust = 0.5, vjust = 0.5) +
  draw_image(dat$Picture[[2]] / 8, x = 0.8, y = 0.4, scale = 0.2, hjust = 0.5, vjust = 0.5) +
  draw_plot(p)

# Plot turning angles of the two species
p <- dat %>%
  dplyr::select(Species, Steps) %>%
  unnest(cols = Steps) %>%
  ggplot(aes(x = ta_)) +
    geom_histogram(aes(y = ..density..), fill = "cornflowerblue", color = "black") +
    geom_rug(color = "cornflowerblue") +
    facet_wrap("Species", scales = "free") +
    xlab("Turning Angle") +
    ylab("Density") +
    scale_x_continuous(breaks  = c(seq(-pi, +pi, pi)), labels = c("-\u03c0", "0", "+\u03c0")) +
    theme_cowplot() +
    theme(strip.background = element_rect(fill = "gray95"))

# Add silhouettes
ggdraw() +
  draw_image(dat$Picture[[1]] / 8, x = 0.2, y = 0.8, scale = 0.1, hjust = 0.5, vjust = 0.5) +
  draw_image(dat$Picture[[2]] / 8, x = 0.7, y = 0.8, scale = 0.1, hjust = 0.5, vjust = 0.5) +
  draw_plot(p)

# Let's create seperate columns where we store the step lengths and turning
# angles. Note that steps need to cover a positive distance and that we need to
# remove any NAs.
dat$sl <- lapply(dat$Steps, function(x){
  as.numeric(na.omit(pmax(x$sl_, 0.1)))
})
dat$ta <- lapply(dat$Steps, function(x){
  as.numeric(na.omit(x$ta_))
})

# Now we can fit a gamma distribution to the step lengths and a mixed von mises
# distribution to the turning angles
dat$dist_sl <- lapply(dat$sl, function(x){
  dist_sl <- fitdist(x, "gamma", method = "mle", lower = 0)
  dist_sl <- list(
      scale = 1 / dist_sl$estimate[["rate"]]
    , shape = dist_sl$estimate[["shape"]]
  )
  return(dist_sl)
})
dat$dist_ta <- lapply(dat$ta, function(x){
  dist_ta <- fitdist(x, "vonMises", start = list(k1 = 1, k2 = 2))
  dist_ta <- list(
      k1 = dist_ta$estimate[["k1"]]
    , k2 = dist_ta$estimate[["k2"]]
  )
  return(dist_ta)
})

################################################################################
#### Visualizations
################################################################################
# Plot the fit of the step lengths
par(mfrow = c(2, 2))
for (n in 1:2){
  x <- seq(0, max(dat$sl[[n]]), by = 1)
  y <- dgamma(x
    , scale = dat$dist_sl[[n]]$scale
    , shape = dat$dist_sl[[n]]$shape
  )
  hist(dat$sl[[n]]
    , freq   = F
    , breaks = 100
    , main   = dat$Species[[n]]
    , xlab   = "Step Length"
  )
  lines(y ~ x, col = "red")
}

# Plot the fit of the turning angles
for (n in 1:2){
  x <- seq(-pi, +pi, by = 0.01)
  y <- dvonMises(x
    , k1 = dat$dist_ta[[n]]$k1
    , k2 = dat$dist_ta[[n]]$k2
  )
  hist(dat$ta[[n]]
    , freq   = F
    , breaks = 50
    , xlim   = c(-pi, +pi)
    , main   = dat$Species[[n]]
    , xlab   = "Turning Angle"
  )
  lines(y ~ x, col = "red")
}

# Store the parameters of the distributions to file
dat %>%
  dplyr::select(Species, dist_sl, dist_ta) %>%
  write_rds("Data/Output/FittedDistributions.rds")

# Store the parameters of the distributions to file, I'll also store a copy to
# the shiny app folder
dat %>%
  dplyr::select(Species, dist_sl, dist_ta) %>%
  write_rds("Shiny/data/FittedDistributions.rds")

# I'll also create a copy of the "00_Functions.R" file, so that the functions
# can be used by the shiny app
file.copy("00_Functions.R", "Shiny")
