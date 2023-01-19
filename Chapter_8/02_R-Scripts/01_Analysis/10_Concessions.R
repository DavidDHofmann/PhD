################################################################################
#### Consession Differences
################################################################################
# Description: Computing differences in connectivity across concessions

# Clear R's brain
rm(list = ls())

# Change the working directory
setwd("/home/david/ownCloud/University/15. PhD/Chapter_8")

# Load required packages
library(terra)     # To handle spatial data
library(raster)    # To handle spatial data
library(tidyverse) # To wrangle data
library(rgdal)

# Load custom functions
source("02_R-Scripts/00_Functions.R")

# Load concession borders
con <- ogrListLayers("/home/david/Schreibtisch/Concession Boundaries.kml")
con <- con[-1]
con <- lapply(con, function(x) {
  bound <- readOGR("/home/david/Schreibtisch/Concession Boundaries.kml", layer = x)
  bound <- vect(bound)
  bound <- as.polygons(bound)
  return(bound)
})
con <- do.call(c, con)

# Load heatmaps
maps1 <- "03_Data/03_Results/HeatmapsBetweennessGlobal.rds" %>%
  read_rds() %>%
  mutate(SourceArea = NA, Level = "Global") %>%
  select(Steps, SourceArea, FloodLevel, Level, Data = Heatmap) %>%
  arrange(Steps, SourceArea, FloodLevel, Level)
maps2 <- "03_Data/03_Results/HeatmapsBetweennessLocal.rds" %>%
  read_rds() %>%
  mutate(Level = "Local") %>%
  select(Steps, SourceArea, FloodLevel, Level, Data = Heatmap) %>%
  arrange(Steps, SourceArea, FloodLevel, Level)
heat <- rbind(maps1, maps2)

# I also want to plot a difference map
diff <- subset(heat, Level == "Global" & Steps == 2000)
diff <- diff$Data[[1]] - diff$Data[[2]]
diff <- rast(diff)

extr <- terra::extract(diff, con, fun = "mean")
values(con) <- data.frame(Difference = extr)
library(sf)
library(scales)
names(con) <- c("ID", "Difference")
ext <- c(21, 25.1, -20.7, -17.9)
consf <- st_as_sf(con)
ggplot(consf, aes(fill = Difference)) +
  geom_sf(col = "gray95") +
  scale_fill_gradientn(
      colors  = colorRampPalette(c("orange", "white", "cornflowerblue"))(100)
    , limits  = c(-300, 300)
    , oob     = squish
    , guide   = guide_colorbar(
      , title          = "Test"
      , show.limits    = T
      , title.position = "bottom"
      , title.hjust    = 0.5
      , ticks          = F
      , barheight      = unit(0.2, "cm")
      , barwidth       = unit(12, "cm")
    )
  ) +
  theme_minimal() +
  theme(legend.position = "bottom") +
    coord_sf(
        crs    = 4326
      , xlim   = ext[1:2]
      , ylim   = ext[3:4]
      , expand = F
    )

diff <- tibble(
    Steps      = 2000
  , SourceArea = NA
  , FloodLevel = "Difference"
  , Level      = "Difference"
  , Data       = list(diff)
)
heat <- rbind(heat, diff)

# Convert the maps to dataframes
heat <- heat %>%
  mutate(Data = map(Data, function(x) {
    result <- as.data.frame(x, xy = T)
    names(result) <- c("x", "y", "Value")
    return(result)
  })) %>% unnest(Data)

# Make sure the levels are correctly ordered
heat$FloodLevel <- factor(heat$FloodLevel, levels = c("Min", "Max", "Difference"))
heat$Metric     <- "Heatmap"
