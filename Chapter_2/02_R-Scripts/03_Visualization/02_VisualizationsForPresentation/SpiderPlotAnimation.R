################################################################################
#### Spider Plot Animation
################################################################################
# Clear R's brain
rm(list = ls())

# Change the working directory
wd <- "/home/david/ownCloud/University/15. PhD/Chapter_2"
setwd(wd)

# Load required packages
library(tidyverse)   # To wrangle data
library(scales)      # To rescale data
library(gganimate)   # To animate plots
library(ggdark)      # For dark themes

################################################################################
#### Still Image
################################################################################
# Dataframe containing seasonal habitat preferences
df_wet <- tibble(
    Covariate   = c("Water", "DistanceToWater", "Grassland", "Woodland", "Humans", "cos_ta", "sl", "log_sl")
#   , Beta        = c(-0.5, 0, -0.2, 0.5, -0.2, -0.5)
  , Beta        = rnorm(length(Covariate), mean = 0, sd = 2)
  , Season      = "Wet"
)
df_dry <- tibble(
    Covariate   = c("Water", "DistanceToWater", "Grassland", "Woodland", "Humans", "cos_ta", "sl", "log_sl")
#   , Beta        = c(-0.2, -0.8, 0.2, +0.2, -0.8)
  , Beta        = rnorm(length(Covariate), mean = 0, sd = 2)
  , Season      = "Dry"
)

# Put the two dataframes together
dat      <- rbind(df_wet, df_dry)
dat$Beta <- rescale(dat$Beta)

# Function to create circular coordinates
spiderCoords <- function(value, distance, angle) {
  x <- (value + distance) * sin(angle)
  y <- (value + distance) * cos(angle)
  xy <- data.frame(x = x, y = y)
  return(xy)
}

# Function to create circles
spiderCircles <- function(ymin = 1, ymax = 2, breaks = 3) {
  r <- seq(ymin, ymax, length.out = breaks)
  angles <- seq(0, 2 * pi, by = 0.05)
  crds <- lapply(r, function(i) {
    xy <- spiderCoords(i, distance = 0, angle = angles)
    xy$r <- i
    return(xy)
  })
  crds <- do.call(rbind, crds)
  return(crds)
}

# Function to create lines
spiderLines <- function(ymin = 1, ymax = 2, breaks = NULL) {
  angles <- seq(0, 2 * pi, length.out = breaks + 1)
  angles <- angles[-length(angles)]
  crds <- lapply(angles, function(i) {
    from <- spiderCoords(ymin, distance = 0, angle = i)
    to   <- spiderCoords(ymax, distance = 0, angle = i)
    xy   <- data.frame(x1 = from[, 1], y1 = from[, 2], x2 = to[, 1], y2 = to[, 2])
    xy$r <- i
    return(xy)
  })
  crds <- do.call(rbind, crds)
  return(crds)
}


# Function to prepare data for spider plot
spiderData <- function(data, value, group, distance = 0) {

  # Create unique group ids
  data <- data %>%
    group_by_(group) %>%
    mutate(GroupID = cur_group_id()) %>%
    ungroup()

  # Calculate angle for each group
  angles <- seq(0, 2 * pi, length.out = length(unique(data$GroupID)) + 1)
  angles <- angles[-length(angles)]
  angles <- data.frame(GroupID = 1:length(angles), Angle = angles)

  # Put together
  data <- left_join(data, angles, by = "GroupID")

  # Compute coordinates for each entry
  xy <- spiderCoords(value = data$Beta, distance = distance, angle = data$Angle)
  data <- cbind(data, xy)
  data <- arrange(data, GroupID)

  # Return it
  return(data)
}

# Function to generate all data for the spider plot
spiderPlotData <- function(data
    , group        = NULL
    , value        = NULL
    , distance     = 0
    , label_offset = 0
  ) {

  # Testing
#   data <- dat
#   value <- "Beta"
#   group <- "Covariate"
#   distance <- 1
#   label_offset <- 1

  # Prepare data
  spider_dat <- spiderData(data, value, group, distance)

  # Get necessary metrics
  groups <- unique(pull(spider_dat, group))
  ngroups <- length(groups)
  maxval <- max(pull(data, value))

  # Prepare circles, lines, and labels
  spider_cir <- spiderCircles(
      ymin = distance
    , ymax = maxval + distance
  )
  spider_lin <- spiderLines(
      ymin   = distance
    , ymax   = maxval + distance
    , breaks = ngroups
  )
  spider_lab <- spiderLines(
      ymin   = distance
    , ymax   = maxval + distance + label_offset
    , breaks = ngroups
  )
  spider_lab$Label <- groups

  # Return all
  return(list(
        Data    = spider_dat
      , Circles = spider_cir
      , Lines   = spider_lin
      , Labels  = spider_lab
    )
  )
}

# Prepare the plot
dat_plot <- spiderPlotData(dat
  , group        = "Covariate"
  , value        = "Beta"
  , distance     = 1
  , label_offset = 0.2
)

# Plot it
ggplot() + 
  geom_polygon(data = dat_plot$Circles
    , aes(x = x, y = y, group = r)
    , fill  = NA
    , col   = "gray80"
    , lty   = 2
  ) +
  geom_segment(data = dat_plot$Lines
    , aes(x = x1, y = y1, xend = x2, yend = y2, group = r)
    , col  = "gray80"
    , lty  = 2
  ) +
  geom_point(data = dat_plot$Data
    , aes(x = x, y = y, col = Season)
    , size = 5
  ) +
  geom_polygon(data = dat_plot$Data
    , aes(x = x, y = y, fill = Season, col = Season)
    , linewidth = 1
    , alpha     = 0.2
  ) +
  geom_text(data = dat_plot$Labels
    , aes(x = x2, y = y2, label = Label)
    , size = 5
  ) +
  coord_equal() +
  dark_theme_void() +
  theme(
      legend.position  = "none"
    , panel.background = element_blank()
    , plot.background  = element_blank()
  ) +
  scale_fill_manual(values = c("cornflowerblue", "orange")) +
  scale_color_manual(values = c("cornflowerblue", "orange")) +
  theme(strip.text = element_blank())

# Store the pot
ggsave("05_Presentation/SpiderPlot.png", bg = "transparent")

################################################################################
#### Animation
################################################################################
# Animate it
p <- ggplot() +
    geom_polygon(data = dat_plot$Circles
      , aes(x = x, y = y, group = r)
      , fill  = NA
      , col   = "gray80"
      , lty   = 2
    ) +
    geom_segment(data = dat_plot$Lines
      , aes(x = x1, y = y1, xend = x2, yend = y2, group = r)
      , col  = "gray80"
      , lty  = 2
    ) +
    geom_polygon(data = dat_plot$Data
      , aes(x = x, y = y)
      , fill      = "orange"
      , col       = "orange"
      , linewidth = 1
      , alpha     = 0.5
    ) +
    geom_point(data = dat_plot$Data
      , aes(x = x, y = y)
      , col  = "orange"
      , size = 5
    ) +
    geom_text(data = dat_plot$Labels
      , aes(x = x2, y = y2, label = Label)
      , size = 5
      , col  = "white"
    ) +
    coord_equal() +
    transition_states(Season, transition_length = 5, state_length = 0) +
    theme_void() +
    theme(
        strip.text = element_blank()
      , plot.background = element_rect(fill = "gray30")
      , panel.background = element_rect(fill = "gray30")
    )
anim_save(filename = "05_Presentation/SpiderAnimation.gif", animation = p)
