################################################################################
#### Visualizations
################################################################################
# Description: Visualizations for the manuscript

# Clear R's brain
rm(list = ls())

# Load required packages
library(raster)        # To handle spatial data
library(tidyverse)     # For data wrangling
library(sf)            # For plotting
library(ggpubr)        # To arrange plots
library(pbmcapply)     # For multicore abilities

# Set working directory
setwd("/home/david/ownCloud/University/15. PhD/Chapter_4")

# Load observed movement data and covariates
obs <- read_csv("03_Data/01_RawData/ObservedMovements.csv")
covars <- stack("03_Data/01_RawData/CovariateLayers.grd")

# Extent in which animals were released
ext <- extent(c(50, 250, 50, 250))
ext <- as(ext, "SpatialPolygons")

################################################################################
#### Plot of Covariates and Simulated Trajectories
################################################################################
# Plot the covariates
p1 <- as.data.frame(covars, xy = T) %>%
  gather(key = covariate, value = value, 3:5) %>%
  ggplot(aes(x = x, y = y, fill = value)) +
    geom_raster() +
    geom_sf(data = st_as_sf(ext), inherit.aes = F, fill = NA, col = "white", lty = 2) +
    scale_fill_viridis_c(option = "viridis") +
    coord_sf() +
    theme_minimal() +
    facet_wrap("covariate") +
    theme(axis.title.y = element_text(angle = 0, vjust = 0.5)) +
    scale_x_continuous(breaks = seq(0, 300, by = 100)) +
    scale_y_continuous(breaks = seq(0, 300, by = 100))

# Plot of trajectories
p2 <- ggplot(obs, aes(x = x, y = y, col = as.factor(ID))) +
  geom_sf(data = st_as_sf(ext), col = "cornflowerblue", fill = NA, inherit.aes = F, lwd = 1.20) +
  geom_path(size = 0.1) +
  geom_point(size = 0.1) +
  coord_sf() +
  theme_minimal() +
  theme(legend.position = "none") +
  scale_color_viridis_d()

# Store to file
ggsave(plot = p1, "04_Manuscript/99_Covariates.png", width = 6, height = 2.5, bg = "white")
ggsave(plot = p2, "04_Manuscript/99_Simulations.png", scale = 0.75, bg = "white")

################################################################################
#### Lost Fixes
################################################################################
# Simulate a dataset of coordinates
set.seed(12345)
dat <- data.frame(
    x           = runif(1000)
  , y           = runif(1000)
  , step_number = 1:1000
)

# Function to create a dataset with rarified observations (the same function is
# also used to subsample x individuals from all individuals)
rarifyData <- function(data, missingness) {
  rarified <- data[sort(sample(1:nrow(data), size = nrow(data) * (1 - missingness))), ]
  return(rarified)
}

# Function to determine the number of steps resulting from a dataframe. Also
# indicate the step durations.
numberSteps <- function(data, forgiveness) {

  # Determine bursts
  data$duration <- lead(data$step_number) - data$step_number
  data$irregular <- data$duration > forgiveness
  data$burst <- NA
  data$burst[1] <- 1
  data$burst[2:nrow(data)] <- lag(cumsum(data$irregular) + 1)[-1]

  # Count total number of fixes and number of valid fixes (that will give a
  # relative turning angle) per burst
  n_steps <- data %>%
    count(burst) %>%
    subset(n > 2) %>%
    mutate(nvalid = n - 2) %>%
    pull(nvalid) %>%
    sum()

  # Return it
  return(n_stepts)
}

# Design matrix to loop through
design <- expand_grid(
    Missingness    = seq(0, 0.9, by = 0.01)  # Fraction of the fixes that is removed
  , Forgiveness    = 1:5                     # Allowed lag of steps (in steps)
  , Replicate      = 1:100                   # Number of replicates for each combination
)

# Run through the design
design$NumberSteps <- pbmclapply(1:nrow(design), ignore.interactive = T, mc.cores = detectCores() - 1, function(x) {
  sub <- rarifyData(dat, missingness = design$Missingness[x])
  res <- numberSteps(sub, forgiveness = design$Forgiveness[x])
  return(res)
}) %>% do.call(c, .)

# Compute summary stats
summaries <- design %>%
  group_by(Missingness, Forgiveness) %>%
  summarize(
      LCI         = quantile(NumberSteps, 0.025)
    , UCI         = quantile(NumberSteps, 0.975)
    , NumberSteps = mean(NumberSteps)
    , .groups     = "drop"
  )

# Visualize
p <- ggplot(summaries, aes(x = Missingness, y = NumberSteps, col = as.factor(Forgiveness), fill = as.factor(Forgiveness), ymin = LCI, ymax = UCI)) +
  geom_ribbon(alpha = 0.2, lwd = 0.4) +
  geom_line(size = 0.4) +
  scale_fill_viridis_d(name = "Forgiveness") +
  scale_color_viridis_d(name = "Forgiveness") +
  theme_minimal() +
  ylab("Number of Steps") +
  theme(legend.position = "bottom")
  # geom_point() +
  # geom_line() +
  # geom_errorbar(width = 0.05)

# Store the plot
ggsave("04_Manuscript/99_NumberOfSteps.png", plot = p, width = 5, height = 3)

################################################################################
#### Step Length Distribution Parameters
################################################################################
# Load fitted distributions
dists_dynamic <- read_rds("03_Data/03_Results/StepLengthDynamic.rds")
dists_means <- read_rds("03_Data/03_Results/StepLengthMeans.rds")

# Visualize everything
labelfun <- function(l) {
  rl <- round(l, 3)
  l <- ifelse(rl == 0 & l != 0, format(l, scientific = T), l)
}
p <- dists_dynamic %>%
  pivot_longer(shape:mu, names_to = "Parameter", values_to = "Value") %>%
  ggplot(aes(x = duration, y = Value)) +
    geom_jitter(width = 0.1, alpha = 0.2, size = 0.5) +
    geom_point(data = dists_means, aes(x = duration, y = mean), col = "orange", size = 5) +
    geom_line(data = dists_means, aes(x = duration, y = mean), col = "orange") +
    scale_y_continuous(labels = labelfun) +
    facet_wrap(~ Parameter, nrow = 2, scales = "free") +
    theme_minimal() +
    xlab("Step Duration") +
    ylab("Parameter Estimate") +
    theme(
        strip.background = element_rect(fill = "gray95", color = NA)
    )

# Store the plot
ggsave("04_Manuscript/99_DistributionParameters.png", plot = p, width = 5, height = 5)
