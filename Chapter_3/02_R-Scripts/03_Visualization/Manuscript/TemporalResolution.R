################################################################################
#### Temporal Resolution of Covariates
################################################################################
# Clear R's brain
rm(list = ls())

# Set working directory
wd <- "/home/david/ownCloud/University/15. PhD/Chapter_3"
setwd(wd)

# Load required packages
library(tidyverse)    # For plotting and data wrangling
library(terra)        # To handle spatial data
library(lubridate)    # To handle dates
library(hms)          # To handle times
library(ggpubr)       # To arrange multiple plots
library(sf)           # For spatial data
library(stars)        # For spatial data
library(ggnewscale)   # For multiple colorscales
library(colorspace)   # For colors

# Load custom functions
source("02_R-Scripts/00_Functions.R")

# Load all spatial covariates
covs1 <- "03_Data/02_CleanData/Covariates.rds" %>%
  read_rds() %>%
  subset(Type == "Dynamic") %>%
  mutate(Dates = map(Dates, function(x) {
    x$Layerdate
  })) %>%
  dplyr::select(-c(Filename, Type))

# Also load the nightly statistics
covs2 <- "03_Data/02_CleanData/Moonlight.rds" %>%
  read_rds() %>%
  pull(Timestamp) %>%
  expand_grid(Date = ., Covariate = c("MoonIllumination", "TimeOfDay")) %>%
  nest(Dates = -Covariate) %>%
  mutate(Dates = map(Dates, function(x) {
    x$Date
  }))

# Put them together
covs <- rbind(covs1, covs2)
print(covs)

################################################################################
#### Covariate Date-Ranges
################################################################################
# Function that provides a "stylized" set of timestamps given our covariate
# timestamps. This is just to get a nice figure and doesn't need to be 100%
# perfect
styledDates <- function(dates, spacing = 0.1, max_days = 365, dt_min = 1) {

  # If there is only a single date, span the entire year
  if (length(dates) == 1) {
    newdays = tibble(Start = 1, End = 366)
    return(newdays)
  }

  # Drop dates coming from leap years
  remove <- yday(dates) > 365
  dates  <- dates[!remove]

  # Convert dates to days
  dates_days <- yday(dates) + hour(dates) / 24

  # Find median time difference between dates
  dt <- lead(dates_days) - dates_days
  dt <- median(dt, na.rm = T)

  # Prepare dataframe of stylized dates
  newdays <- tibble(
      Start = seq(min(dates_days), max(dates_days), by = dt)
    , End   = Start + (1 - spacing) * dt
  )

  # Return stylized dates
  return(newdays)
}

# Define spacings
covs$Spacing = c(0.1, 0.1, 0.02, 0.02, 0.02, 0.02, 0.02, 0.01, rep(0.1, 4))
covs$Styleddates <- lapply(1:nrow(covs), function(i) {
  styledDates(covs$Dates[[i]], spacing = covs$Spacing[i], dt_min = 1)
})

# Prepare figure
p1 <- covs %>%
  dplyr::select(Covariate, Styleddates) %>%
  unnest(Styleddates) %>%
  mutate(Covariate = factor(Covariate
    , levels = c("Humans", "Forest", "Trees", "Shrubs", "Water", "DistanceToWater", "DistanceToPans", "NDVI", "Temperature", "Precipitation", "MoonIllumination", "TimeOfDay")
  )) %>%
  ggplot(aes(x = Start, xend = End, y = Covariate, ystart = Covariate, yend = Covariate)) +
    geom_segment(size = 3, col = "gray75") +
    # geom_point(aes(x = (Day + ToDay) / 2), pch = "|") +
    theme_minimal() +
    xlab("Day of the Year") +
    ylab("Covariate") +
    scale_x_continuous(breaks = seq(0, 350, by = 5), limits = c(0, 366), position = "top") +
    scale_y_discrete(limits = rev) +
    coord_cartesian(xlim = c(10, 50)) +
    theme(
        panel.grid.minor = element_blank()
      , legend.position  = "none"
      , strip.background = element_rect(fill = "gray95", color = "white")
      , axis.line.x      = element_line()
      , axis.ticks.x     = element_line()
      , axis.text.y      = element_text(size = 6)
    )

################################################################################
#### Covariate Stack
################################################################################
# Simualte some covariates to visualize
set.seed(123)
dim  <- 10
n_layers <- length(covs$Covariate)
covars <- lapply(1:n_layers, function(x) {
  r   <- rast(xmin = 0, xmax = dim, ymin = 0, ymax = dim, res = 1)
  r[] <- rnorm(ncell(r))
  r   <- focal(r, fun = mean, expand = T)
  r   <- normalizeRaster(r)
  return(r)
}) %>% rast() %>% setNames(c(covs$Covariate))

# Also derive a frame so we can smooth out the edges of the layers in the plot
frame <- as.polygons(ext(covars), crs = crs(covars))

# Function to rotate data
rotateSF <- function(data, x_add = 0, y_add = 0, value = 0) {

  # Convert object to type sf
  if (inherits(data, "SpatVector")) {
      data <- st_as_sf(data)
    } else if (inherits(data, "SpatRaster")) {
      data <- st_as_stars(data)
      data <- st_as_sf(data)
  }

  # Create necessary matrices
  shear_matrix <- function() {
    rbind(
        c(2, 0)
      , c(1, 0.5)
    )
  }
  rotate_matrix <- function(value) {
    rbind(
        c(cos(value), -sin(value))
      , c(sin(value), cos(value))
    )
  }

  # Run transformation
  data <- mutate(data, geometry = geometry *
    shear_matrix() *
    rotate_matrix(value) +
    c(x_add, y_add)
  )

  # Return
  return(data)
}

# Create a random covariate layer for the redistribution kernel
set.seed(123)
r <- rast(xmin = 0, xmax = dim, ymin = 0, ymax = dim)
r[] <- rnorm(ncell(r))
r <- focal(r, fun = mean, expand = T, w = focalMat(r, d = 0.1, type = "Gauss"))

# Creat a grid of points around the center of the raster layers
set.seed(123)
steps <- expand_grid(
    x_to = seq(0, dim, length.out = 100)
  , y_to = seq(0, dim, length.out = 100)
  ) %>%
  mutate(
      x_from = dim / 2
    , y_from = dim / 2
    , sl     = sqrt((x_to - x_from) ** 2 + (y_to - y_from) ** 2)
    , absta  = (atan2(y_to - y_from, x_to - x_from) - pi / 2) * (-1)
    , absta  = ifelse(absta < 0, 2 * pi + absta, absta)
    , ta     = absta - runif(1, min = 0, max = 2 * pi)
    , ta     = ifelse(ta > +pi, ta - 2 * pi, ta)
    , ta     = ifelse(ta < -pi, 2 * pi + ta, ta)
    , cos_ta = cos(ta)
    , log_sl = log(sl)
  )

# Extract covariate
steps$covariate <- extract(r, steps[, c("x_to", "y_to")])$focal_mean

# Compute a step-selection score for each step
formula     <- ~ sl + log_sl + cos_ta + covariate
betas       <- c(-0.5, 1, 0.5, 0.5)
steps$probs <- as.vector(exp(model.matrix(formula, steps)[, -1] %*% betas))
steps       <- vect(steps, geom = c("x_to", "y_to")) %>% st_as_sf()

# Plot values
spacing    <- 1.5
framewidth <- 0.5
framecolor <- "gray30"
gridcolor  <- "gray80"

# Prepare the plot
set.seed(123)
p2 <- ggplot() +
  geom_sf(
      data    = rotateSF(steps, y_add = -15)
    , mapping = aes(col = probs)
    , size    = 1
  ) +
  geom_sf(
      data        = rotateSF(frame, y_add = -15)
    , fill        = "transparent"
    , col         = framecolor
    , linewidth   = 2
    , show.legend = F
  ) +
  theme_void() +
  scale_color_viridis_c(name = "Step-Probability") +
  guides(
    color = guide_colorbar(
      , title.position = "top"
      , title.hjust    = 0.5
      , ticks          = F
      , barheight      = unit(0.2, "cm")
      , barwidth       = unit(6, "cm")
      , label          = F
    )
  ) +
  theme(legend.position = "bottom")

# Loop through the covariate layers and add them ontop
for (i in 1:nlyr(covars)) {
  col <- hcl.colors(n_layers)[i]
  flip <- sample(c(-1, 0, 1), size = 1)
  if (flip == -1) {
      layer <- disagg(covars[[i]], fact = 2, method = "bilinear")
    } else if (flip == 1) {
      layer <- aggregate(covars[[i]], fact = 2, fun = "mean")
    } else {
      layer <- covars[[i]]
  }
  names(layer) <- "Value"
  p2 <- p2 +
    geom_sf(
        data        = rotateSF(layer, y_add = spacing * (i - 1))
      , mapping     = aes(fill = Value)
      , col         = gridcolor
      , show.legend = F
    ) +
    geom_sf(
        data        = rotateSF(frame, y_add = spacing * (i - 1))
      , fill        = "transparent"
      , col         = framecolor
      , linewidth   = framewidth
      , show.legend = F
    ) +
    scale_fill_gradient(low = col, high = lighten(col, 0.9)) +
    new_scale_fill()
}

################################################################################
#### Store the Plots
################################################################################
# Store the plots
ggsave("04_Manuscript/Figures/Schematic1.png"
  , plot   = p1
  , width  = 10
  , height = 1.75
  , scale  = 1
  , bg     = "white"
  , device = png
)
ggsave("04_Manuscript/Figures/Schematic2.png"
  , plot   = p2
  , width  = 3
  , height = 4
  , scale  = 1.5
  , bg     = "white"
  , device = png
)
