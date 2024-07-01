################################################################################
#### Covariates
################################################################################
# Description: Plot of the covariates

# Clear R's brain
rm(list = ls())

# Load required packages
library(raster)        # To handle spatial data
library(tidyverse)     # For data wrangling
library(sf)            # For plotting
library(pbmcapply)     # For multicore abilities
library(ggh4x)         # For nested facet_wraps

# Set working directory
setwd("/home/david/ownCloud/University/15. PhD/Chapter_4")

# Load observed movement data and covariates
# obs <- read_csv("03_Data/01_RawData/ObservedMovements.csv")
dat <- "03_Data/Simulation.rds" %>%
  read_rds() %>%
  subset(Replicate == 2) %>%
  dplyr::select(AutocorrRange, Covariates, Movement)

# Extent in which animals were released
ext <- extent(c(50, 250, 50, 250))
ext <- as(ext, "SpatialPolygons")

################################################################################
#### Plot of Covariates and Simulated Trajectories
################################################################################
# Prepare the data
dat_cov <- dat %>%
  subset(AutocorrRange == 10) %>%
  pull(Covariates) %>%
  stack() %>%
  as.data.frame(xy = T) %>%
  gather(key = covariate, value = value, 3:5) %>%
  mutate(covariate = gsub(covariate, pattern = "dist", replacement = "Dist")) %>%
  mutate(covariate = gsub(covariate, pattern = "elev", replacement = "Elev")) %>%
  mutate(covariate = gsub(covariate, pattern = "forest", replacement = "Forest"))
dat_mov <- dat %>%
  subset(AutocorrRange == 10) %>%
  dplyr::select(Movement) %>%
  unnest(Movement)

# Plot
p1 <- ggplot(dat_cov, aes(x = x, y = y, fill = value)) +
    geom_raster() +
    geom_sf(data = st_as_sf(ext), inherit.aes = F, fill = NA, col = "white", lty = 2, linewidth = 0.5) +
    geom_path(data = dat_mov, inherit.aes = F, aes(x = x, y = y), linewidth = 0.2) +
    scale_fill_viridis_c(option = "viridis", name = "Value") +
    coord_sf() +
    theme_minimal() +
    facet_wrap(~ covariate) +
    scale_x_continuous(breaks = seq(0, 300, by = 100), limits = c(0, 300)) +
    scale_y_continuous(breaks = seq(0, 300, by = 100), limits = c(0, 300)) +
    theme(axis.title.y = element_text(angle = 0, vjust = 0.5)) +
    theme(strip.background = element_rect(fill = "gray95", color = "white"))

# Prepare the data
dat_cov <- dat %>%
  dplyr::select(AutocorrRange, Covariates) %>%
  mutate(Covariates = map(Covariates, function(x) {
    as.data.frame(x, xy = T)
  })) %>%
  unnest(Covariates) %>%
  gather(key = covariate, value = value, 4:6) %>%
  mutate(covariate = gsub(covariate, pattern = "dist", replacement = "Dist")) %>%
  mutate(covariate = gsub(covariate, pattern = "elev", replacement = "Elev")) %>%
  mutate(covariate = gsub(covariate, pattern = "forest", replacement = "Forest"))
dat_mov <- dat %>%
  dplyr::select(AutocorrRange, Movement) %>%
  unnest(Movement)

# Plot
p2 <- ggplot(dat_cov, aes(x = x, y = y, fill = value)) +
    geom_raster() +
    geom_sf(data = st_as_sf(ext), inherit.aes = F, fill = NA, col = "white", lty = 2) +
    geom_path(data = dat_mov, inherit.aes = F, aes(x = x, y = y), linewidth = 0.2) +
    scale_fill_viridis_c(option = "viridis", name = "Value") +
    coord_sf() +
    theme_minimal() +
    facet_grid(AutocorrRange ~ covariate) +
    scale_x_continuous(breaks = seq(0, 300, by = 100), limits = c(0, 300)) +
    scale_y_continuous(breaks = seq(0, 300, by = 100), limits = c(0, 300)) +
    theme(axis.title.y = element_text(angle = 0, vjust = 0.5)) +
    theme(
        strip.background = element_rect(fill = "gray95", color = "white")
      , strip.text.y = element_text(angle = 0)
    )

# Store to file
ggsave("04_Manuscript/99_Covariates.png"
  , plot   = p1
  , width  = 6
  , height = 2.5
  , bg     = "white"
  , scale  = 1.2
  , device = png
)
ggsave("04_Manuscript/99_CovariatesAppendix.png"
  , plot   = p2
  , width  = 6
  , height = 5
  , bg     = "white"
  , scale  = 1.2
  , device = png
)
