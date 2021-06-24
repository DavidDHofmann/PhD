################################################################################
#### Distributing Cameras & Quantify Success
################################################################################
# Description: In this script we distribute the camera traps in the study area
# Authors: David, Eva, & Dilsad
# Date: November 2020

# Clear R's brain
rm(list = ls())

# Load required packages
library(raster)       # To handle spatial data
library(rgeos)        # To manipulate spatial data
library(tidyverse)    # For data wrangling
library(rgdal)        # To read and write spatial data
library(pbmcapply)    # For multicore abilities with progress bar
library(cowplot)      # For nice ggplots

# Set the working directory
setwd("/home/david/ownCloud/Dokumente/Bibliothek/Wissen/R-Scripts/ProgrammingProject")

# Load custom functions
source("00_Functions.R")

# Set seed for reproducability
set.seed(1234)

################################################################################
#### Distribute Camera Traps
################################################################################
# Let's reload the study area, the boundaries of switzerland, and our paths
area_cams <- readOGR("Data/Output/CameraArea.shp")
area_move <- readOGR("Data/Output/MovementArea.shp")
che <- readOGR("Data/Output/Switzerland.shp")
paths <- readOGR("Data/Output/Simulations.shp")

# Distribute camera traps and vary the number of cameras
cams <- lapply(1:12, function(x){
  distributeCams(area_cams, n = x, range = 100)
})

# Store the geometries in a nice tibble
cams <- tibble(
    NCams  = (1:length(cams)) ** 2
  , Grid   = lapply(cams, function(x){x$Grid})
  , Points = lapply(cams, function(x){x$Points})
  , Range  = lapply(cams, function(x){x$Range})
)

# Visualize some of the grids
n <- 6
plot(cams$Grid[[n]]
  , border = "black"
  , lwd    = 1
  , main = paste0("#Cameras: ", n ** 2)
)
plot(cams$Range[[n]]
  , add     = T
  , col     = adjustcolor("red"
  , alpha.f = 0.2)
  , border  = "red"
)
plot(cams$Points[[n]]
  , add = T
  , pch = 20
  , cex = 0.1
)

# Store the grids to file
write_rds(cams, "Data/Output/CameraGrids.rds")

################################################################################
#### Quantify Success
################################################################################
# Let's check how many individuals we detected
detections <- suppressWarnings(
  pbmclapply(1:nrow(cams), ignore.interactive = T,
    function(x){
      camranges <- cams$Range[[x]]
      camranges <- aggregate(camranges)
      as.vector(gIntersects(camranges, paths, byid = T))
  })
)

# Do some data cleaning
detections <- detections %>%
  do.call(cbind, .) %>%
  as.data.frame() %>%
  setNames(as.character((1:ncol(.)) ** 2)) %>%
  cbind(paths@data, .) %>%
  gather(key = "NCams", value = "Detection", 3:ncol(.)) %>%
  mutate(NCams = as.numeric(NCams)) %>%
  group_by(Species, NCams) %>%
  summarize(
      Success = sum(Detection)
    , Failure = sum(!Detection)
    , SuccessRate = Hmisc::binconf(sum(Detection), n(), alpha = 0.05, method = "exact")[, 1]
    , LowerConf = Hmisc::binconf(sum(Detection), n(), alpha = 0.05, method = "exact")[, 2]
    , UpperConf = Hmisc::binconf(sum(Detection), n(), alpha = 0.05, method = "exact")[, 3]
  )

# Visualize them
ggplot(detections, aes(x = NCams, y = SuccessRate, fill = Species, col = Species)) +
  geom_ribbon(aes(ymin = LowerConf, ymax = UpperConf), alpha = 0.2, col = NA) +
  geom_point() +
  geom_line() +
  theme_cowplot() +
  scale_color_manual(values = c("chartreuse4", "cornflowerblue")) +
  scale_fill_manual(values = c("chartreuse4", "cornflowerblue"))

################################################################################
#### Determine Costs
################################################################################
# function for the total cost of the cameras
camCosts <- function(cameras, uni, camera_price){

  # Specify speeds (in km/h)
  speed_travel <- 50
  speed_setup <- 10
  time_installation <- 0.5
  salary <- 50
  gas_money <- 1

  # Get the number of cameras
  ncams <- length(cameras)

  # Calculate distance from the uni to study area
  distance2area <- gDistance(gCentroid(cameras), gCentroid(uni)) / 1000

  # Calculate total distance to drive between cams
  distance2cams <- getDistance(cameras) / 1000

  # Travel gas costs (to get to study area)
  gas_travel <- 2 * distance2area * gas_money

  # Setup gas costs (to set up the cameras)
  gas_setup <- 2 * distance2cams * gas_money

  # Travel time (to get to the study area)
  travel_time <- 2 * distance2area / speed_travel

  # Setup time (to set up the cameras)
  setup_time <- 2 * distance2cams / speed_setup + ncams * time_installation

  # Salary costs (for travel and for setting up cameras)
  salary_travel <- salary * travel_time
  salary_setup <- salary * setup_time

  # Calculate discount for the cameras (the more cameras, the higher the
  # discount, yet limited to 40%)
  discount <- (-exp(-0.05 * ncams) + 1) * 40 / 100

  # Calculate costs for buying the cameras
  costs_cameras <- ncams * camera_price * (1 - discount)

  # Calculate total costs
  costs_fix <- salary_travel + gas_travel
  costs_var <- salary_setup + costs_cameras + gas_setup
  costs_tot <- costs_fix + costs_var

  # Prepare a dataframe with the different costs
  costs <- data.frame(
      Description = c("Salary Travel", "Salary Setup", "Gas Travel", "Gas Setup", "Camera Costs", "Total")
    , TypeI       = c("Fix", "Variable", "Fix", "Variable", "Variable", NA)
    , TypeII      = c("Salary", "Salary", "Gas", "Gas", "Camera", "Total")
    , Value       = c(salary_travel, salary_setup, gas_travel, gas_setup, costs_cameras, costs_tot)
  )

  # Return the dataframe
  return(costs)
}

# Specify location of uni
uni <- SpatialPoints(data.frame(x = 465942, y = 5249470))

# Apply the function to the different camera grids
cams$Costs <- lapply(cams$Points, function(x){
  camCosts(x, uni = uni, camera_price = 200)
})

# Create a dataframe for the costs
costs <- cams %>%
  dplyr::select(NCams, Costs) %>%
  unnest(cols = Costs) %>%
  group_by(NCams, Type = TypeII) %>%
  summarize(Cost = sum(Value)) %>%
  mutate(LowerConf = Cost - 0.1 * Cost) %>%
  mutate(UpperConf = Cost + 0.1 * Cost) %>%
  mutate(Type = factor(Type, levels = c("Total", "Camera", "Salary", "Gas")))

# Visualize them
ggplot(costs, aes(x = NCams, y = Cost, col = Type, fill = Type)) +
  geom_ribbon(aes(ymin = LowerConf, ymax = UpperConf), alpha = 0.2, col = NA) +
  geom_line() +
  geom_point() +
  theme_cowplot() +
  scale_y_continuous(labels = function(x) format(x, big.mark = "'", scientific = FALSE)) +
  scale_color_manual(values = c("cornflowerblue", "chartreuse4", "darkorchid", "darkgoldenrod1")) +
  scale_fill_manual(values = c("cornflowerblue", "chartreuse4", "darkorchid", "darkgoldenrod1"))

################################################################################
#### Compare Benefits and Costs
################################################################################
# Compile dataframe with costs vs benefits
comp <- detections %>%
  group_by(NCams) %>%
  summarize(SuccessRate = mean(SuccessRate)) %>%
  left_join(., subset(costs, Type == "Total") %>% dplyr::select(c(NCams, Cost)))

# Let's see how we could translate the success rate into "money"
x <- seq(0, 1, 0.01)
y <- x ** 0.4 * 50000
plot(y ~ x, type = "l")

# Calculate monetary benefit
comp$Benefit <- comp$SuccessRate ** 0.4 * 50000

# Calculate net benefit
comp$NetBenefit <- comp$Benefit - comp$Cost

# Identify at how many cameras the maximum is reached
opt_ben <- comp[comp$NetBenefit == max(comp$NetBenefit), ]
opt_ben$Type <- "Total"

# Let's also check the detection rate that we'd get
opt_det1 <- detections %>%
  subset(NCams == opt_ben$NCams) %>%
  group_by(NCams, Species)

# Plot net benefit
ggplot(comp, aes(x = NCams, y = NetBenefit)) +
  geom_line() +
  geom_point() +
  theme_cowplot()

################################################################################
#### Identify Number of Cameras We Can Buy
################################################################################
# Specify your budget
budget <- 10000

# Let's see how many cams we could buy with this
opt_cost <- costs %>%
  subset(Type == "Total" & Cost <= budget) %>%
  subset(Cost == max(Cost))

# Let's also check the detection rate that we'd get
opt_det2 <- detections %>%
  subset(NCams == opt_cost$NCams) %>%
  group_by(NCams, Species)

################################################################################
#### Visualizations
################################################################################
# Visualize detections
p1 <- ggplot(detections, aes(x = NCams, y = SuccessRate, fill = Species, col = Species)) +
  geom_ribbon(aes(ymin = LowerConf, ymax = UpperConf), alpha = 0.2, col = NA) +
  geom_point() +
  geom_line() +
  theme_cowplot() +
  scale_color_manual(values = c("chartreuse4", "cornflowerblue")) +
  scale_fill_manual(values = c("chartreuse4", "cornflowerblue")) +
  geom_segment(data = opt_det1, aes(xend = -Inf, yend = SuccessRate), col = "gray20", lty = 2) +
  geom_segment(data = opt_det1, aes(xend = NCams, yend = -Inf), col = "gray20", lty = 2) +
  geom_segment(data = opt_det2, aes(xend = -Inf, yend = SuccessRate), col = "gray50", lty = 2) +
  geom_segment(data = opt_det2, aes(xend = NCams, yend = -Inf), col = "gray50", lty = 2)

# Visualize the costs
p2 <- ggplot(costs, aes(x = NCams, y = Cost, col = Type, fill = Type)) +
  geom_ribbon(aes(ymin = LowerConf, ymax = UpperConf), alpha = 0.2, col = NA) +
  geom_line() +
  geom_point() +
  theme_cowplot() +
  scale_y_continuous(labels = function(x) format(x, big.mark = "'", scientific = FALSE)) +
  scale_color_manual(values = c("cornflowerblue", "chartreuse4", "darkorchid", "darkgoldenrod1")) +
  scale_fill_manual(values = c("cornflowerblue", "chartreuse4", "darkorchid", "darkgoldenrod1")) +
  geom_segment(data = opt_ben, aes(xend = -Inf, yend = Cost), col = "gray20", lty = 2) +
  geom_segment(data = opt_ben, aes(xend = NCams, yend = -Inf), col = "gray20", lty = 2) +
  geom_segment(data = opt_cost, aes(xend = -Inf, yend = Cost), col = "gray50", lty = 2) +
  geom_segment(data = opt_cost, aes(xend = NCams, yend = -Inf), col = "gray50", lty = 2)
