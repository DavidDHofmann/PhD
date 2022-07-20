################################################################################
#### Lost Fixes
################################################################################
# Description: Plot of fixes lost due to irregularity

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
ggsave("04_Manuscript/99_NumberOfSteps.png"
  , plot   = p
  , width  = 5
  , height = 3
  , bg     = "white"
  , scale  = 1.4
)
