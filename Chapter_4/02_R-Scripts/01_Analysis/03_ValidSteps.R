################################################################################
#### Valid Steps
################################################################################
# Description: Analyse how number of valid steps changes as the missingness and
# forgiveness increase.

# Clear R's brain
rm(list = ls())

# Load required packages
library(tidyverse)     # For data wrangling
library(pbmcapply)     # For multicore abilities

# Set working directory
setwd("/home/david/ownCloud/University/15. PhD/Chapter_4")

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
  return(n_steps)
}

# Design matrix to loop through
design <- expand_grid(
    Missingness    = seq(0, 0.9, by = 0.01)  # Fraction of the fixes that is removed
  , Forgiveness    = 1:5                     # Allowed lag of steps (in steps)
  , Replicate      = 1:1000                  # Number of replicates for each combination
)

# Simulate rarified datasets for the above design matrix and compute number of
# usable steps for different levels of missingness
design$NumberSteps <- pbmclapply(1:nrow(design), ignore.interactive = T, mc.cores = detectCores() - 1, function(x) {
  sub <- rarifyData(dat, missingness = design$Missingness[x])
  res <- numberSteps(sub, forgiveness = design$Forgiveness[x])
  return(res)
}) %>% do.call(c, .)

# Write results to file
write_csv(design, "03_Data/ValidSteps.csv")
